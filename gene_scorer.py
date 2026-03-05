"""
Gene Scorer Module
==================
GWAS P値と遺伝子アノテーション情報に基づくスコアリング

スコアリング要素:
1. GWAS P値: -log10(p_value) を 0-10 に正規化して使用
2. 遺伝子機能アノテーション:
   - LoF (Loss of Function): stop_gained, frameshift 等
   - GoF (Gain of Function): missense + damaging
   - Missense: ミスセンスバリアント
   - Coding: コーディング領域バリアント
   - Regulatory: 調節領域バリアント
"""

# P値スコアの最大値 (これ以上は飽和)
PVALUE_SCORE_CAP = 10.0

import time
import requests
import pandas as pd
import numpy as np

import config

# LoF consequence terms
LOF_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
    "transcript_ablation",
}

# GoF / damaging missense 判定用
GOF_CONSEQUENCES = {
    "missense_variant",  # SIFT+PolyPhen で damaging の場合
}

MISSENSE_CONSEQUENCES = {
    "missense_variant",
}

CODING_CONSEQUENCES = {
    "synonymous_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "incomplete_terminal_codon_variant",
    "coding_sequence_variant",
}

REGULATORY_CONSEQUENCES = {
    "regulatory_region_variant", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "TFBS_ablation", "TFBS_amplification",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
}


def calculate_gwas_pvalue_score(p_value: float) -> float:
    """
    GWAS P値からスコアを計算
    -log10 変換後そのままスコアとして返す (上限なし)

    Parameters
    ----------
    p_value : float
        GWAS P値

    Returns
    -------
    float : スコア (-log10(P値))
    """
    if p_value is None or p_value <= 0:
        return 0.0
    return -np.log10(max(p_value, 1e-300))


def fetch_variant_consequences(rsids: list) -> pd.DataFrame:
    """
    Ensembl VEP REST API でバリアントの機能的影響を取得

    Parameters
    ----------
    rsids : list
        rs ID のリスト

    Returns
    -------
    pd.DataFrame : columns=[rsid, gene_symbol, consequence, impact,
                            sift_prediction, polyphen_prediction]
    """
    if not rsids:
        return pd.DataFrame()

    print(f"[VEP] {len(rsids)} バリアントの機能アノテーション取得中...")
    records = []

    # バッチ処理 (VEP POST は最大200件)
    batch_size = 200
    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i + batch_size]

        url = f"{config.ENSEMBL_REST_URL}/vep/human/id"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        payload = {"ids": batch}

        max_retries = 3
        for attempt in range(max_retries):
            try:
                resp = requests.post(url, json=payload, headers=headers, timeout=180)
                
                if resp.status_code == 429:
                    retry_after = int(resp.headers.get("Retry-After", 5))
                    print(f"  → Rate limit, {retry_after}秒待機...")
                    time.sleep(retry_after)
                    continue
                    
                resp.raise_for_status()
                data = resp.json()

                for variant in data:
                    rsid = variant.get("id", "")
                    for tc in variant.get("transcript_consequences", []):
                        gene = tc.get("gene_symbol", "")
                        consequences = tc.get("consequence_terms", [])
                        impact = tc.get("impact", "")
                        sift = tc.get("sift_prediction", "")
                        polyphen = tc.get("polyphen_prediction", "")

                        if gene:
                            for csq in consequences:
                                records.append({
                                    "rsid": rsid,
                                    "gene_symbol": gene.upper(),
                                    "consequence": csq,
                                    "impact": impact,
                                    "sift_prediction": sift,
                                    "polyphen_prediction": polyphen,
                                })

                time.sleep(0.5)
                break  # 成功したらループを抜ける

            except requests.exceptions.RequestException as e:
                print(f"[VEP] バッチ {i // batch_size + 1} ({len(batch)}件) リトライ {attempt+1}/{max_retries} エラー: {e}")
                if attempt < max_retries - 1:
                    time.sleep(5 * (attempt + 1))  # Exponential backoff
                else:
                    print(f"[VEP] バッチ {i // batch_size + 1} 取得失敗 (最大リトライ到達)")

    result = pd.DataFrame(records)
    if not result.empty:
        print(f"[VEP] {len(result)} アノテーション取得 ({result['rsid'].nunique()} バリアント)")
    return result


def classify_variant_effects(consequences_df: pd.DataFrame) -> pd.DataFrame:
    """
    バリアント結果を機能カテゴリに分類

    Parameters
    ----------
    consequences_df : pd.DataFrame
        fetch_variant_consequences の出力

    Returns
    -------
    pd.DataFrame : columns=[rsid, gene_symbol, is_missense, is_lof, is_gof,
                             is_coding, is_regulatory, impact_severity]
    """
    if consequences_df.empty:
        return pd.DataFrame()

    # 遺伝子ごとに集約
    grouped = consequences_df.groupby(["rsid", "gene_symbol"])

    records = []
    for (rsid, gene), group in grouped:
        csq_set = set(group["consequence"])

        is_lof = bool(csq_set & LOF_CONSEQUENCES)

        # GoF: missense + damaging prediction
        is_missense = bool(csq_set & MISSENSE_CONSEQUENCES)
        sift_vals = set(group["sift_prediction"].dropna())
        polyphen_vals = set(group["polyphen_prediction"].dropna())
        is_damaging = ("deleterious" in sift_vals or
                       "probably_damaging" in polyphen_vals or
                       "possibly_damaging" in polyphen_vals)
        is_gof = is_missense and is_damaging

        is_coding = bool(csq_set & CODING_CONSEQUENCES) or is_missense or is_lof
        is_regulatory = bool(csq_set & REGULATORY_CONSEQUENCES)

        # 影響度スコア
        impact_vals = set(group["impact"])
        if "HIGH" in impact_vals:
            impact_severity = 4
        elif "MODERATE" in impact_vals:
            impact_severity = 3
        elif "LOW" in impact_vals:
            impact_severity = 2
        elif "MODIFIER" in impact_vals:
            impact_severity = 1
        else:
            impact_severity = 0

        records.append({
            "rsid": rsid,
            "gene_symbol": gene,
            "is_missense": is_missense,
            "is_lof": is_lof,
            "is_gof": is_gof,
            "is_coding": is_coding,
            "is_regulatory": is_regulatory,
            "impact_severity": impact_severity,
        })

    return pd.DataFrame(records)


def score_genes(gwas_df: pd.DataFrame,
                variant_classifications: pd.DataFrame = None,
                ppi_df: pd.DataFrame = None) -> pd.DataFrame:
    """
    遺伝子の総合スコアを計算 (GWAS P値 + VEPアノテーション)

    Parameters
    ----------
    gwas_df : pd.DataFrame
        GWAS データ (columns: rsid, gene_symbol, p_value)
    variant_classifications : pd.DataFrame
        classify_variant_effects の出力
    ppi_df : pd.DataFrame
        PPI データ (互換性のため保持、スコアには使用しない)

    Returns
    -------
    pd.DataFrame : 遺伝子ごとの総合スコア (is_gwas 列でGWAS/PPI区別)
    """
    weights = config.SCORING_WEIGHTS

    # 遺伝子ごとに最良 P値を取得
    gene_pvalues = gwas_df.groupby("gene_symbol")["p_value"].min().to_dict()
    gwas_gene_set = set(gene_pvalues.keys())

    # PPI遺伝子も含めた全遺伝子セット
    all_genes = set(gene_pvalues.keys())
    if ppi_df is not None and not ppi_df.empty:
        all_genes |= set(ppi_df["gene_a"]) | set(ppi_df["gene_b"])

    # スコア計算
    records = []

    for gene in sorted(all_genes):
        is_gwas = gene in gwas_gene_set

        # GWAS P値スコア
        p_val = gene_pvalues.get(gene, 1.0)
        pval_score = calculate_gwas_pvalue_score(p_val) * weights["gwas_pvalue"] if is_gwas else 0.0

        # バリアント影響度スコア (VEP)
        lof_score = 0.0
        gof_score = 0.0
        missense_score = 0.0
        coding_score = 0.0
        regulatory_score = 0.0

        if variant_classifications is not None and not variant_classifications.empty:
            gene_vars = variant_classifications[
                variant_classifications["gene_symbol"] == gene
            ]
            if not gene_vars.empty:
                if gene_vars["is_lof"].any():
                    lof_score = weights["lof"]
                if gene_vars["is_gof"].any():
                    gof_score = weights["gof"]
                if gene_vars["is_missense"].any():
                    missense_score = weights["missense"]
                if gene_vars["is_coding"].any():
                    coding_score = weights["coding_variant"]
                if gene_vars["is_regulatory"].any():
                    regulatory_score = weights["regulatory"]

        # PPI マルチソーススコア
        ppi_score = 0.0
        if ppi_df is not None and not ppi_df.empty:
            gene_ppi = ppi_df[
                (ppi_df["gene_a"] == gene) | (ppi_df["gene_b"] == gene)
            ]
            n_sources = gene_ppi["source"].nunique()
            ppi_score = n_sources * weights["ppi_multi_source"]

        total = pval_score + lof_score + gof_score + missense_score + coding_score + regulatory_score + ppi_score

        records.append({
            "gene_symbol": gene,
            "is_gwas": is_gwas,
            "p_value": p_val if is_gwas else None,
            "pval_score": round(pval_score, 3),
            "lof_score": round(lof_score, 3),
            "gof_score": round(gof_score, 3),
            "missense_score": round(missense_score, 3),
            "coding_score": round(coding_score, 3),
            "regulatory_score": round(regulatory_score, 3),
            "ppi_score": round(ppi_score, 3),
            "total_score": round(total, 3),
        })

    result = pd.DataFrame(records).sort_values("total_score", ascending=False).reset_index(drop=True)
    n_gwas = result["is_gwas"].sum()
    n_ppi = len(result) - n_gwas
    print(f"[Scorer] {len(result)} 遺伝子をスコアリング (GWAS: {n_gwas}, PPI: {n_ppi})")
    return result
