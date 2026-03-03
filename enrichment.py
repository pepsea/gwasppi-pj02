"""
Enrichment Module (Local)
=========================
ローカルで Over-Representation Analysis を実行
対応データベース: Reactome, Gene Ontology (GO), HPO

Fisher's exact test + Benjamini-Hochberg FDR 補正
"""

import os
import io
import zipfile
import requests
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict

import config

# ============================================================
# GMT キャッシュ
# ============================================================
_reactome_cache = None
_go_cache = None
_hpo_cache = None


# ============================================================
# Reactome GMT
# ============================================================
def _download_reactome_gmt() -> dict:
    """Reactome GMT をダウンロードし、パスウェイ→遺伝子セット辞書を返す"""
    global _reactome_cache
    if _reactome_cache is not None:
        return _reactome_cache

    local_path = os.path.join(config.DATA_DIR, "ReactomePathways.gmt")

    gmt_text = None
    if os.path.exists(local_path):
        print(f"  → ローカル GMT 読み込み: {local_path}")
        with open(local_path, "r") as f:
            gmt_text = f.read()
    else:
        print(f"  → Reactome GMT ダウンロード中...")
        try:
            resp = requests.get(config.REACTOME_GMT_URL, timeout=60)
            resp.raise_for_status()
            z = zipfile.ZipFile(io.BytesIO(resp.content))
            gmt_text = z.read(z.namelist()[0]).decode("utf-8")
            with open(local_path, "w") as f:
                f.write(gmt_text)
            print(f"  → 保存: {local_path}")
        except Exception as e:
            print(f"  → Reactome GMT ダウンロードエラー: {e}")
            return {}

    pathways = {}
    for line in gmt_text.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        pathway_name = parts[0]
        pathway_id = parts[1]
        genes = set(g.upper() for g in parts[2:] if g.strip())
        if pathway_id.startswith("R-HSA-"):
            pathways[pathway_id] = {"name": pathway_name, "genes": genes}

    print(f"  → Reactome パスウェイ: {len(pathways)} (ヒト)")
    _reactome_cache = pathways
    return pathways


# ============================================================
# Gene Ontology GMT
# ============================================================
def _download_go_gmt() -> dict:
    """Gene Ontology GMT (MSigDB c5) をダウンロード"""
    global _go_cache
    if _go_cache is not None:
        return _go_cache

    local_path = os.path.join(config.DATA_DIR, "go_c5.gmt")

    gmt_text = None
    if os.path.exists(local_path):
        print(f"  → ローカル GO GMT 読み込み: {local_path}")
        with open(local_path, "r") as f:
            gmt_text = f.read()
    else:
        print(f"  → GO GMT ダウンロード中...")
        try:
            resp = requests.get(config.GO_GMT_URL, timeout=120)
            resp.raise_for_status()
            gmt_text = resp.text
            with open(local_path, "w") as f:
                f.write(gmt_text)
            print(f"  → 保存: {local_path}")
        except Exception as e:
            print(f"  → GO GMT ダウンロードエラー: {e}")
            return {}

    pathways = {}
    for line in gmt_text.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        go_name = parts[0]
        go_url = parts[1]
        genes = set(g.upper() for g in parts[2:] if g.strip())
        pathways[go_name] = {"name": go_name, "genes": genes, "url": go_url}

    print(f"  → GO ターム: {len(pathways)}")
    _go_cache = pathways
    return pathways


# ============================================================
# HPO (Human Phenotype Ontology)
# ============================================================
def _download_hpo_data() -> dict:
    """HPO genes_to_phenotype.txt をダウンロードし、HPO term → 遺伝子セットに変換"""
    global _hpo_cache
    if _hpo_cache is not None:
        return _hpo_cache

    local_path = os.path.join(config.DATA_DIR, "genes_to_phenotype.txt")

    lines = None
    if os.path.exists(local_path):
        print(f"  → ローカル HPO 読み込み: {local_path}")
        with open(local_path, "r") as f:
            lines = f.readlines()
    else:
        print(f"  → HPO データダウンロード中...")
        try:
            resp = requests.get(config.HPO_ANNOTATION_URL, timeout=60)
            resp.raise_for_status()
            lines = resp.text.strip().split("\n")
            with open(local_path, "w") as f:
                f.write(resp.text)
            print(f"  → 保存: {local_path}")
        except Exception as e:
            print(f"  → HPO ダウンロードエラー: {e}")
            return {}

    # HPO term → 遺伝子セット構築
    hpo_terms = defaultdict(lambda: {"name": "", "genes": set()})
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            gene_symbol = parts[1].upper().strip() if len(parts) > 1 else ""
            hpo_id = parts[2].strip() if len(parts) > 2 else ""
            hpo_name = parts[3].strip() if len(parts) > 3 else ""
            if gene_symbol and hpo_id:
                hpo_terms[hpo_id]["name"] = hpo_name
                hpo_terms[hpo_id]["genes"].add(gene_symbol)

    result = dict(hpo_terms)
    print(f"  → HPO ターム: {len(result)}")
    _hpo_cache = result
    return result


# ============================================================
# Fisher's Exact Test + BH FDR 補正
# ============================================================
def _benjamini_hochberg(p_values: list) -> list:
    """Benjamini-Hochberg FDR 補正"""
    n = len(p_values)
    if n == 0:
        return []

    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    fdr = [0.0] * n
    min_so_far = 1.0

    for rank_idx in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1
        adjusted = pval * n / rank
        min_so_far = min(min_so_far, adjusted)
        fdr[orig_idx] = min(min_so_far, 1.0)

    return fdr


def _run_enrichment(gene_list: list, pathways: dict, db_name: str,
                    fdr_threshold: float = None,
                    min_size: int = 5, max_size: int = 500,
                    background_size: int = None) -> pd.DataFrame:
    """
    Fisher's exact test によるエンリッチメント解析 (共通関数)
    """
    if fdr_threshold is None:
        fdr_threshold = config.ENRICHMENT_FDR_THRESHOLD
    if background_size is None:
        background_size = config.BACKGROUND_GENE_COUNT

    gene_set = set(g.upper().strip() for g in gene_list if g and g.strip())
    n_input = len(gene_set)

    if not pathways or n_input == 0:
        return pd.DataFrame()

    results = []
    for pw_id, pw_data in pathways.items():
        pw_genes = pw_data["genes"]
        pw_size = len(pw_genes)

        if pw_size < min_size or pw_size > max_size:
            continue

        overlap = gene_set & pw_genes
        k = len(overlap)
        if k == 0:
            continue

        a = k
        b = n_input - k
        c = pw_size - k
        d = max(0, background_size - pw_size - b)

        _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        results.append({
            "term_id": pw_id,
            "term_name": pw_data["name"],
            "database": db_name,
            "p_value": p_value,
            "gene_count": k,
            "total_genes": pw_size,
            "gene_ratio": k / pw_size if pw_size > 0 else 0,
            "fold_enrichment": (k / n_input) / (pw_size / background_size) if pw_size > 0 else 0,
            "genes": ", ".join(sorted(overlap)),
        })

    if not results:
        return pd.DataFrame()

    # BH FDR 補正
    p_values = [r["p_value"] for r in results]
    fdr_values = _benjamini_hochberg(p_values)
    for r, fdr_val in zip(results, fdr_values):
        r["fdr"] = fdr_val

    df = pd.DataFrame(results)
    df = df.sort_values("fdr").reset_index(drop=True)
    sig_df = df[df["fdr"] <= fdr_threshold].reset_index(drop=True)

    return sig_df


# ============================================================
# 公開 API
# ============================================================
def reactome_enrichment(gene_list: list, **kwargs) -> pd.DataFrame:
    """Reactome Pathway Enrichment"""
    print(f"\n{'='*60}")
    print(f"Reactome Enrichment: {len(gene_list)} 遺伝子")
    print(f"{'='*60}")
    pathways = _download_reactome_gmt()
    result = _run_enrichment(gene_list, pathways, "Reactome", **kwargs)
    print(f"[Reactome] 有意パスウェイ: {len(result)}")
    return result


def go_enrichment(gene_list: list, **kwargs) -> pd.DataFrame:
    """Gene Ontology Enrichment"""
    print(f"\n{'='*60}")
    print(f"GO Enrichment: {len(gene_list)} 遺伝子")
    print(f"{'='*60}")
    pathways = _download_go_gmt()
    result = _run_enrichment(gene_list, pathways, "GO", **kwargs)
    print(f"[GO] 有意ターム: {len(result)}")
    return result


def hpo_enrichment(gene_list: list, min_size: int = 3, **kwargs) -> pd.DataFrame:
    """HPO (Human Phenotype Ontology) Enrichment"""
    print(f"\n{'='*60}")
    print(f"HPO Enrichment: {len(gene_list)} 遺伝子")
    print(f"{'='*60}")
    pathways = _download_hpo_data()
    result = _run_enrichment(gene_list, pathways, "HPO", min_size=min_size, **kwargs)
    print(f"[HPO] 有意ターム: {len(result)}")
    return result


def run_all_enrichment(gene_list: list,
                       databases: list = None,
                       **kwargs) -> pd.DataFrame:
    """
    全データベースで Enrichment 解析を実行

    Parameters
    ----------
    gene_list : list
        遺伝子シンボルのリスト
    databases : list
        使用するデータベース (デフォルト: ["reactome", "go", "hpo"])

    Returns
    -------
    pd.DataFrame : 統合 Enrichment 結果
    """
    if databases is None:
        databases = ["reactome", "go", "hpo"]

    all_results = []

    if "reactome" in databases:
        df = reactome_enrichment(gene_list, **kwargs)
        if not df.empty:
            all_results.append(df)

    if "go" in databases:
        df = go_enrichment(gene_list, **kwargs)
        if not df.empty:
            all_results.append(df)

    if "hpo" in databases:
        df = hpo_enrichment(gene_list, **kwargs)
        if not df.empty:
            all_results.append(df)

    if not all_results:
        return pd.DataFrame()

    result = pd.concat(all_results, ignore_index=True)
    result = result.sort_values("fdr").reset_index(drop=True)

    print(f"\n[Enrichment] 全体結果: {len(result)} 有意ターム")
    for db in result["database"].unique():
        n = len(result[result["database"] == db])
        print(f"  - {db}: {n} ターム")

    return result
