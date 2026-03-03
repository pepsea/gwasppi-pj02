"""
GWAS Fetcher Module
===================
GWAS Catalog REST API から SNP → 遺伝子マッピングデータを取得
"""

import requests
import pandas as pd
import time

import config


def search_studies_by_disease(disease_trait: str, max_results: int = 100) -> list:
    """
    疾患名で GWAS Catalog のスタディを検索

    Parameters
    ----------
    disease_trait : str
        疾患名 (例: "Type 2 Diabetes")
    max_results : int
        最大取得件数

    Returns
    -------
    list : スタディ情報のリスト
    """
    url = f"{config.GWAS_CATALOG_BASE_URL}/studies/search/findByDiseaseTrait"
    params = {
        "diseaseTrait": disease_trait,
        "page": 0,
        "size": max_results,
    }

    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        studies = []
        embedded = data.get("_embedded", {})
        for study in embedded.get("studies", []):
            studies.append({
                "study_id": study.get("accessionId", ""),
                "title": study.get("publicationInfo", {}).get("title", ""),
                "pubmed_id": study.get("publicationInfo", {}).get("pubmedId", ""),
            })

        print(f"[GWAS] '{disease_trait}' で {len(studies)} スタディを取得")
        return studies

    except Exception as e:
        print(f"[GWAS] スタディ検索エラー: {e}")
        return []


def fetch_associations_by_study(study_id: str, p_value_threshold: float = None) -> pd.DataFrame:
    """
    スタディ ID から関連アソシエーションを取得

    Parameters
    ----------
    study_id : str
        GWAS Catalog スタディ ID (例: "GCST001234")
    p_value_threshold : float
        P値の閾値 (None の場合は config のデフォルトを使用)

    Returns
    -------
    pd.DataFrame : アソシエーション情報
    """
    if p_value_threshold is None:
        p_value_threshold = config.GWAS_P_VALUE_THRESHOLD

    url = f"{config.GWAS_CATALOG_BASE_URL}/studies/{study_id}/associations"
    params = {"size": 1000}

    try:
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        associations = []
        embedded = data.get("_embedded", {})
        for assoc in embedded.get("associations", []):
            # P値の取得 (pvalue or pvalueMantissa * 10^pvalueExponent)
            p_value = assoc.get("pvalue", None)
            if p_value is None:
                mantissa = assoc.get("pvalueMantissa", None)
                exponent = assoc.get("pvalueExponent", None)
                if mantissa is not None and exponent is not None:
                    try:
                        p_value = float(mantissa) * (10 ** int(exponent))
                    except (ValueError, TypeError):
                        p_value = 1.0
                else:
                    p_value = 1.0
            else:
                try:
                    p_value = float(p_value)
                except (ValueError, TypeError):
                    p_value = 1.0

            if p_value > p_value_threshold:
                continue

            # Loci から rsID と遺伝子を取得
            loci = assoc.get("loci", [])
            for locus in loci:
                # rsID は strongestRiskAlleles の riskAlleleName から抽出
                risk_alleles = locus.get("strongestRiskAlleles", [])
                rsids = []
                for ra in risk_alleles:
                    risk_name = ra.get("riskAlleleName", "")
                    if risk_name:
                        # "rs10906115-A" → "rs10906115"
                        rsid = risk_name.split("-")[0].strip()
                        if rsid.startswith("rs"):
                            rsids.append(rsid)

                # 遺伝子は authorReportedGenes から取得
                reported_genes = locus.get("authorReportedGenes", [])
                gene_names = []
                for gene in reported_genes:
                    gn = gene.get("geneName", "")
                    if gn and gn.lower() not in ("intergenic", "nr", "n/a", ""):
                        gene_names.append(gn)

                # rsID × gene の組み合わせを作成
                if not rsids:
                    rsids = ["unknown"]
                if not gene_names:
                    gene_names = ["intergenic"]

                for rsid in rsids:
                    for gene_name in gene_names:
                        associations.append({
                            "rsid": rsid,
                            "gene_symbol": gene_name.upper(),
                            "p_value": p_value,
                            "study_id": study_id,
                        })

        df = pd.DataFrame(associations)
        if not df.empty:
            df = df.drop_duplicates(subset=["rsid", "gene_symbol"])
            # "intergenic" や "unknown" を除外
            df = df[df["gene_symbol"] != "INTERGENIC"]
            df = df[df["rsid"] != "unknown"]
        return df

    except Exception as e:
        print(f"[GWAS] アソシエーション取得エラー (Study {study_id}): {e}")
        return pd.DataFrame()


def fetch_snps_for_disease(disease_trait: str, max_studies: int = None,
                           p_value_threshold: float = None) -> pd.DataFrame:
    """
    疾患名から SNPs と関連遺伝子を一括取得 (メイン関数)

    Parameters
    ----------
    disease_trait : str
        疾患名
    max_studies : int
        検索するスタディの最大数 (None=configの設定値)
    p_value_threshold : float
        P値の閾値

    Returns
    -------
    pd.DataFrame : columns=[rsid, gene_symbol, p_value, study_id]
    """
    if max_studies is None:
        max_studies = getattr(config, "DEFAULT_MAX_GWAS_STUDIES", 50)

    print(f"\n{'='*60}")
    print(f"GWAS Catalog: '{disease_trait}' の SNPs を検索中 (最大 {max_studies} スタディ)...")
    print(f"{'='*60}")

    studies = search_studies_by_disease(disease_trait, max_results=max_studies)
    if not studies:
        print("[GWAS] スタディが見つかりませんでした")
        return pd.DataFrame(columns=["rsid", "gene_symbol", "p_value", "study_id"])

    all_associations = []
    for i, study in enumerate(studies[:max_studies]):
        study_id = study["study_id"]
        print(f"  [{i+1}/{min(len(studies), max_studies)}] Study {study_id} を処理中...")
        df = fetch_associations_by_study(study_id, p_value_threshold)
        if not df.empty:
            all_associations.append(df)
            print(f"    → {len(df)} アソシエーション")
        else:
            print(f"    → 0 アソシエーション (閾値未満)")
        time.sleep(0.5)  # rate limiting

    if not all_associations:
        print("[GWAS] アソシエーションが見つかりませんでした")
        return pd.DataFrame(columns=["rsid", "gene_symbol", "p_value", "study_id"])

    result = pd.concat(all_associations, ignore_index=True)
    result = result.drop_duplicates(subset=["rsid", "gene_symbol"])
    result = result.sort_values("p_value").reset_index(drop=True)

    # intergenic / LOC 遺伝子を除外
    result = result[~result["gene_symbol"].str.contains(
        r"^$|INTERGENIC|NR_|LOC\d+", case=False, na=False
    )].reset_index(drop=True)

    print(f"\n[GWAS] 結果サマリー:")
    print(f"  - SNPs 数: {result['rsid'].nunique()}")
    print(f"  - 関連遺伝子数: {result['gene_symbol'].nunique()}")
    print(f"  - アソシエーション数: {len(result)}")

    return result


def fetch_snps_by_efo_trait(efo_trait_id: str, p_value_threshold: float = None) -> pd.DataFrame:
    """
    EFO Trait ID を使って関連アソシエーションを検索する代替メソッド

    Parameters
    ----------
    efo_trait_id : str
        EFO Trait ID (例: "EFO_0001360" for Type 2 Diabetes)
    p_value_threshold : float
        P値の閾値

    Returns
    -------
    pd.DataFrame : columns=[rsid, gene_symbol, p_value, efo_trait_id]
    """
    if p_value_threshold is None:
        p_value_threshold = config.GWAS_P_VALUE_THRESHOLD

    print(f"\n[GWAS/EFO] {efo_trait_id} のアソシエーションを取得中...")

    url = f"{config.GWAS_CATALOG_BASE_URL}/efoTraits/{efo_trait_id}/associations"
    params = {"size": 1000}

    try:
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        associations = []
        embedded = data.get("_embedded", {})
        for assoc in embedded.get("associations", []):
            p_value = assoc.get("pvalue", None)
            if p_value is None:
                mantissa = assoc.get("pvalueMantissa", None)
                exponent = assoc.get("pvalueExponent", None)
                if mantissa is not None and exponent is not None:
                    try:
                        p_value = float(mantissa) * (10 ** int(exponent))
                    except (ValueError, TypeError):
                        p_value = 1.0
                else:
                    p_value = 1.0
            else:
                try:
                    p_value = float(p_value)
                except (ValueError, TypeError):
                    p_value = 1.0

            if p_value > p_value_threshold:
                continue

            loci = assoc.get("loci", [])
            for locus in loci:
                risk_alleles = locus.get("strongestRiskAlleles", [])
                rsids = []
                for ra in risk_alleles:
                    risk_name = ra.get("riskAlleleName", "")
                    if risk_name:
                        rsid = risk_name.split("-")[0].strip()
                        if rsid.startswith("rs"):
                            rsids.append(rsid)

                reported_genes = locus.get("authorReportedGenes", [])
                gene_names = []
                for gene in reported_genes:
                    gn = gene.get("geneName", "")
                    if gn and gn.lower() not in ("intergenic", "nr", "n/a", ""):
                        gene_names.append(gn.upper())

                if not rsids:
                    continue
                if not gene_names:
                    continue

                for rsid in rsids:
                    for gene_name in gene_names:
                        associations.append({
                            "rsid": rsid,
                            "gene_symbol": gene_name,
                            "p_value": p_value,
                            "efo_trait_id": efo_trait_id,
                        })

        df = pd.DataFrame(associations)
        if not df.empty:
            df = df.drop_duplicates(subset=["rsid", "gene_symbol"])
            df = df.sort_values("p_value").reset_index(drop=True)

        print(f"[GWAS/EFO] {efo_trait_id}: {len(df)} アソシエーション取得")
        return df

    except Exception as e:
        print(f"[GWAS/EFO] エラー: {e}")
        return pd.DataFrame(columns=["rsid", "gene_symbol", "p_value", "efo_trait_id"])
