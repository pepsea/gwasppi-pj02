"""
PPI Fetcher Module
==================
4つのPPIデータベース (SIGNOR, BioGRID, STRING, Reactome) から
タンパク質間相互作用データを取得・統合
1～3階層の繋がりに対応
"""

import io
import requests
import pandas as pd
import time

import config

# ============================================================
# SIGNOR
# ============================================================
_signor_cache = None


def _load_signor_human_data():
    """SIGNOR の全ヒトインタラクションデータを TSV でダウンロード (キャッシュ付き)"""
    global _signor_cache
    if _signor_cache is not None:
        return _signor_cache

    print("  → SIGNOR 全データダウンロード中 (初回のみ)...")
    url = "https://signor.uniroma2.it/getData.php"
    params = {"organism": "9606"}

    try:
        resp = requests.get(url, params=params, timeout=120)
        if resp.status_code != 200:
            print(f"  → SIGNOR ダウンロード失敗: HTTP {resp.status_code}")
            return pd.DataFrame()

        # 手動でTSVをパース (pd.read_csvでは失敗するため)
        lines = resp.text.strip().split("\n")
        rows = []
        for line in lines:
            cols = line.split("\t")
            if len(cols) < 28:
                continue

            entity_a = cols[0].strip()
            type_a = cols[1].strip()
            entity_b = cols[4].strip()
            type_b = cols[5].strip()
            effect = cols[8].strip()
            mechanism = cols[9].strip()
            organism = cols[12].strip()

            try:
                score = float(cols[27]) if cols[27].strip() else 0.5
            except (ValueError, IndexError):
                score = 0.5

            rows.append({
                "entity_a": entity_a,
                "type_a": type_a,
                "entity_b": entity_b,
                "type_b": type_b,
                "effect": effect,
                "mechanism": mechanism,
                "organism": organism,
                "score": score,
            })

        df = pd.DataFrame(rows)
        # ヒト protein-protein のみフィルタ
        df = df[df["organism"] == "9606"]
        protein_types = {"protein", "complex"}
        df = df[df["type_a"].isin(protein_types) & df["type_b"].isin(protein_types)]

        print(f"  → SIGNOR: {len(df)} protein-protein インタラクション (human)")
        _signor_cache = df
        return df

    except Exception as e:
        print(f"  → SIGNOR ダウンロードエラー: {e}")
        return pd.DataFrame()


def fetch_signor_interactions(gene_symbols: list) -> pd.DataFrame:
    """
    SIGNOR データベースから指定遺伝子の相互作用を取得

    Parameters
    ----------
    gene_symbols : list
        遺伝子シンボルのリスト

    Returns
    -------
    pd.DataFrame : columns=[gene_a, gene_b, source, interaction_type, score]
    """
    print("\n[PPI/SIGNOR] データ取得中...")

    signor_all = _load_signor_human_data()
    if signor_all.empty:
        print("[PPI/SIGNOR] 0 インタラクションを取得")
        return pd.DataFrame(columns=["gene_a", "gene_b", "source", "interaction_type", "score"])

    gene_set = set(g.upper() for g in gene_symbols)

    # 指定遺伝子が entity_a または entity_b に含まれるものを抽出
    mask = (signor_all["entity_a"].str.upper().isin(gene_set) |
            signor_all["entity_b"].str.upper().isin(gene_set))
    filtered = signor_all[mask].copy()

    interactions = []
    for _, row in filtered.iterrows():
        effect = row["effect"]
        mechanism = row["mechanism"]
        itype = f"{mechanism} ({effect})" if mechanism else effect if effect else "interaction"

        interactions.append({
            "gene_a": row["entity_a"].upper(),
            "gene_b": row["entity_b"].upper(),
            "source": "SIGNOR",
            "interaction_type": itype,
            "score": row["score"],
        })

    df = pd.DataFrame(interactions)
    if not df.empty:
        df = df.drop_duplicates(subset=["gene_a", "gene_b"])

    print(f"[PPI/SIGNOR] {len(df)} インタラクションを取得")
    return df


# ============================================================
# BioGRID
# ============================================================
def fetch_biogrid_interactions(gene_symbols: list, access_key: str = None) -> pd.DataFrame:
    """BioGRID REST API から相互作用を取得"""
    if access_key is None:
        access_key = config.BIOGRID_ACCESS_KEY

    if not access_key:
        print("[PPI/BioGRID] アクセスキーなし、スキップ")
        return pd.DataFrame(columns=["gene_a", "gene_b", "source", "interaction_type", "score"])

    records = []
    # バッチ処理 (API制限回避)
    batch_size = 50
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i + batch_size]
        params = {
            "accesskey": access_key,
            "format": "json",
            "searchNames": "true",
            "geneList": "|".join(batch),
            "organism": config.DEFAULT_TAXID,
            "includeInteractors": "true",
            "includeInteractorInteractions": "false",
        }

        try:
            resp = requests.get(config.BIOGRID_BASE_URL, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()

            for _, interaction in data.items():
                if isinstance(interaction, dict):
                    gene_a = interaction.get("OFFICIAL_SYMBOL_A", "").upper()
                    gene_b = interaction.get("OFFICIAL_SYMBOL_B", "").upper()
                    if gene_a and gene_b and gene_a != gene_b:
                        records.append({
                            "gene_a": gene_a,
                            "gene_b": gene_b,
                            "source": "BioGRID",
                            "interaction_type": interaction.get("EXPERIMENTAL_SYSTEM", "physical"),
                            "score": 1.0,
                        })

            time.sleep(0.3)
        except Exception as e:
            print(f"[PPI/BioGRID] バッチ {i // batch_size + 1} エラー: {e}")

    result = pd.DataFrame(records)
    if not result.empty:
        result = result.drop_duplicates(subset=["gene_a", "gene_b"]).reset_index(drop=True)
    print(f"[PPI/BioGRID] {len(result)} インタラクション")
    return result


# ============================================================
# STRING
# ============================================================
def fetch_string_interactions(gene_symbols: list,
                              species: int = None,
                              min_score: int = None) -> pd.DataFrame:
    """STRING API からタンパク質間相互作用を取得"""
    if species is None:
        species = config.STRING_SPECIES
    if min_score is None:
        min_score = config.STRING_MIN_SCORE

    records = []
    batch_size = 100
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i + batch_size]

        url = f"{config.STRING_BASE_URL}/json/network"
        params = {
            "identifiers": "%0d".join(batch),
            "species": species,
            "required_score": min_score,
            "caller_identity": "gwasppi_pj02",
        }

        try:
            resp = requests.get(url, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()

            for item in data:
                gene_a = item.get("preferredName_A", "").upper()
                gene_b = item.get("preferredName_B", "").upper()
                score = item.get("score", 0)

                if gene_a and gene_b and gene_a != gene_b:
                    records.append({
                        "gene_a": gene_a,
                        "gene_b": gene_b,
                        "source": "STRING",
                        "interaction_type": "functional",
                        "score": score,
                    })

            time.sleep(0.5)
        except Exception as e:
            print(f"[PPI/STRING] バッチ {i // batch_size + 1} エラー: {e}")

    result = pd.DataFrame(records)
    if not result.empty:
        result = result.drop_duplicates(subset=["gene_a", "gene_b"]).reset_index(drop=True)
    print(f"[PPI/STRING] {len(result)} インタラクション (score >= {min_score / 1000:.1f})")
    return result


# ============================================================
# Reactome
# ============================================================
def fetch_reactome_interactions(gene_symbols: list) -> pd.DataFrame:
    """Reactome Content Service から相互作用を取得"""
    records = []

    for gene in gene_symbols:
        try:
            # まず遺伝子を Reactome エンティティに変換
            url = f"{config.REACTOME_CONTENT_URL}/search/query"
            params = {"query": gene, "species": "Homo sapiens", "types": "Protein"}
            resp = requests.get(url, params=params, timeout=30)

            if resp.status_code != 200:
                continue

            data = resp.json()
            entries = data.get("results", [])

            for entry in entries[:1]:  # 最初のカテゴリのみ
                for item in entry.get("entries", [])[:1]:
                    entity_id = item.get("stId", "")
                    if entity_id:
                        _fetch_reactome_interactors(entity_id, gene, records)

            time.sleep(0.3)
        except Exception as e:
            continue

    result = pd.DataFrame(records) if records else pd.DataFrame(
        columns=["gene_a", "gene_b", "source", "interaction_type", "score"]
    )

    if not result.empty:
        result = result.drop_duplicates(subset=["gene_a", "gene_b"]).reset_index(drop=True)

    print(f"[PPI/Reactome] {len(result)} インタラクション")
    return result


def _fetch_reactome_interactors(entity_id: str, source_gene: str, interactions: list):
    """Reactome のエンティティからインタラクターを取得する内部関数"""
    try:
        url = f"{config.REACTOME_CONTENT_URL}/interactors/static/molecule/{entity_id}/details"
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            return

        data = resp.json()
        for entity in data.get("entities", []):
            interactors = entity.get("interactors", [])
            for interactor in interactors:
                alias = interactor.get("alias", "")
                if alias and alias.upper() != source_gene.upper():
                    interactions.append({
                        "gene_a": source_gene.upper(),
                        "gene_b": alias.upper(),
                        "source": "Reactome",
                        "interaction_type": "physical",
                        "score": interactor.get("score", 1.0),
                    })
    except Exception:
        pass


# ============================================================
# 統合・階層展開
# ============================================================
def fetch_all_ppi(gene_symbols: list,
                  sources: list = None,
                  string_min_score: int = None,
                  biogrid_key: str = None) -> pd.DataFrame:
    """
    選択された PPI データベースから統合的にインタラクションを取得

    Parameters
    ----------
    gene_symbols : list
        遺伝子シンボルのリスト
    sources : list
        使用するPPIソース (デフォルト: config.DEFAULT_PPI_SOURCES)
        ["signor", "reactome", "string", "biogrid"] から選択
    string_min_score : int
        STRING API の最小スコア (デフォルト: config.STRING_MIN_SCORE)
    biogrid_key : str
        BioGRID アクセスキー

    Returns
    -------
    pd.DataFrame : columns=[gene_a, gene_b, source, interaction_type, score]
    """
    if sources is None:
        sources = config.DEFAULT_PPI_SOURCES
    sources = [s.lower() for s in sources]

    print(f"\n[PPI] 統合PPI取得: {len(gene_symbols)} 遺伝子, ソース: {sources}")

    all_dfs = []

    if "signor" in sources:
        df = fetch_signor_interactions(gene_symbols)
        if not df.empty:
            all_dfs.append(df)

    if "biogrid" in sources:
        df = fetch_biogrid_interactions(gene_symbols, access_key=biogrid_key)
        if not df.empty:
            all_dfs.append(df)

    if "string" in sources:
        df = fetch_string_interactions(gene_symbols, min_score=string_min_score)
        if not df.empty:
            all_dfs.append(df)

    if "reactome" in sources:
        df = fetch_reactome_interactions(gene_symbols)
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        print("[PPI] インタラクションなし")
        return pd.DataFrame(columns=["gene_a", "gene_b", "source", "interaction_type", "score"])

    result = pd.concat(all_dfs, ignore_index=True)
    print(f"[PPI] 統合結果: {len(result)} インタラクション, "
          f"{len(set(result['gene_a']) | set(result['gene_b']))} 遺伝子")
    return result


def fetch_multi_layer_ppi(seed_genes: list,
                          layers: int = None,
                          sources: list = None,
                          string_min_score: int = None,
                          biogrid_key: str = None) -> pd.DataFrame:
    """
    多階層 PPI 取得

    Parameters
    ----------
    seed_genes : list
        シード遺伝子 (SNP関連遺伝子)
    layers : int
        階層数 (1～3, デフォルト: config.DEFAULT_PPI_LAYERS)
    sources : list
        PPI ソース
    string_min_score : int
        STRING 最小スコア
    biogrid_key : str
        BioGRID キー

    Returns
    -------
    pd.DataFrame : 全階層統合の PPI インタラクション
    """
    if layers is None:
        layers = config.DEFAULT_PPI_LAYERS
    layers = max(1, min(3, layers))  # 1～3 にクランプ

    print(f"\n{'='*60}")
    print(f"多階層 PPI 取得: {layers} 階層")
    print(f"シード遺伝子: {len(seed_genes)}")
    print(f"{'='*60}")

    all_interactions = []
    current_genes = list(set(g.upper() for g in seed_genes))
    seen_genes = set(current_genes)

    for layer in range(1, layers + 1):
        print(f"\n--- 階層 {layer}/{layers}: クエリ遺伝子 {len(current_genes)} ---")

        layer_ppi = fetch_all_ppi(
            current_genes,
            sources=sources,
            string_min_score=string_min_score,
            biogrid_key=biogrid_key,
        )

        if layer_ppi.empty:
            print(f"[PPI] 階層 {layer}: インタラクションなし、終了")
            break

        # 階層情報を追加
        layer_ppi["layer"] = layer
        all_interactions.append(layer_ppi)

        # 次の階層用に新しい遺伝子を収集
        new_genes = set()
        for _, row in layer_ppi.iterrows():
            new_genes.add(row["gene_a"])
            new_genes.add(row["gene_b"])

        new_genes = new_genes - seen_genes
        if not new_genes:
            print(f"[PPI] 階層 {layer}: 新規遺伝子なし、終了")
            break

        seen_genes.update(new_genes)
        current_genes = list(new_genes)
        print(f"[PPI] 階層 {layer}: {len(new_genes)} 新規遺伝子を次の階層へ")

    if not all_interactions:
        return pd.DataFrame(columns=["gene_a", "gene_b", "source", "interaction_type", "score", "layer"])

    result = pd.concat(all_interactions, ignore_index=True)
    result = result.drop_duplicates(subset=["gene_a", "gene_b", "source"]).reset_index(drop=True)

    all_genes = set(result["gene_a"]) | set(result["gene_b"])
    print(f"\n[PPI] 多階層結果: {len(result)} インタラクション, {len(all_genes)} 遺伝子")
    return result


def identify_ppi_neighbors(gwas_genes: list, ppi_df: pd.DataFrame) -> pd.DataFrame:
    """
    PPI データからGWAS遺伝子の相互作用パートナーを特定

    Parameters
    ----------
    gwas_genes : list
        GWAS 遺伝子リスト
    ppi_df : pd.DataFrame
        PPI インタラクション

    Returns
    -------
    pd.DataFrame : columns=[gwas_gene, ppi_gene, source, layer]
    """
    gwas_set = set(g.upper() for g in gwas_genes)
    records = []

    for _, row in ppi_df.iterrows():
        gene_a = row["gene_a"]
        gene_b = row["gene_b"]
        source = row["source"]
        layer = row.get("layer", 1)

        if gene_a in gwas_set and gene_b not in gwas_set:
            records.append({"gwas_gene": gene_a, "ppi_gene": gene_b, "source": source, "layer": layer})
        elif gene_b in gwas_set and gene_a not in gwas_set:
            records.append({"gwas_gene": gene_b, "ppi_gene": gene_a, "source": source, "layer": layer})

    result = pd.DataFrame(records)
    if not result.empty:
        result = result.drop_duplicates(subset=["gwas_gene", "ppi_gene"]).reset_index(drop=True)
    return result
