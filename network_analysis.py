"""
Network Analysis Module
=======================
PPI ネットワークのフローベース中心性解析

readme.md の要件:
- SNP関連遺伝子からのフローと考えて、最初の流量を定義
  1. GWAS P値 (-log10)
  2. LoF/GoF/Missense の variant 影響度
- フローをベースにネットワークの流れを計算し、重要遺伝子をスコア化
"""

import networkx as nx
import numpy as np
import pandas as pd


def build_ppi_network(ppi_df: pd.DataFrame, gwas_genes: list,
                      gene_scores: dict = None) -> nx.Graph:
    """
    PPI DataFrame からネットワークを構築し、GWAS 遺伝子スコアを付与

    Parameters
    ----------
    ppi_df : pd.DataFrame
        PPI インタラクション (columns: gene_a, gene_b, source, score)
    gwas_genes : list
        GWAS 関連遺伝子リスト
    gene_scores : dict
        遺伝子スコア {gene: total_score}

    Returns
    -------
    nx.Graph : ノード属性付き PPI グラフ
    """
    G = nx.Graph()
    gwas_set = set(g.upper() for g in gwas_genes)

    if gene_scores is None:
        gene_scores = {}

    # エッジ追加
    for _, row in ppi_df.iterrows():
        gene_a = row["gene_a"]
        gene_b = row["gene_b"]
        source = row.get("source", "unknown")
        score = row.get("score", 1.0)
        layer = row.get("layer", 1)

        if G.has_edge(gene_a, gene_b):
            # 既存エッジにソース追加
            existing = G[gene_a][gene_b]
            existing_sources = existing.get("sources", set())
            existing_sources.add(source)
            existing["sources"] = existing_sources
            existing["weight"] = max(existing.get("weight", 0), score)
        else:
            G.add_edge(gene_a, gene_b,
                       weight=score,
                       sources={source},
                       layer=layer)

    # ノード属性
    for node in G.nodes():
        G.nodes[node]["is_gwas"] = node in gwas_set
        G.nodes[node]["gene_score"] = gene_scores.get(node, 0.0)
        G.nodes[node]["type"] = "gwas" if node in gwas_set else "ppi"

    n_gwas_in_net = sum(1 for n in G.nodes() if G.nodes[n]["is_gwas"])
    print(f"[Network] グラフ: {G.number_of_nodes()} ノード, "
          f"{G.number_of_edges()} エッジ "
          f"(GWAS遺伝子: {n_gwas_in_net})")

    return G


def compute_centrality_metrics(G: nx.Graph) -> pd.DataFrame:
    """
    ネットワーク中心性指標を計算

    Parameters
    ----------
    G : nx.Graph
        PPI ネットワーク

    Returns
    -------
    pd.DataFrame : 各ノードの中心性指標
    """
    if G.number_of_nodes() == 0:
        return pd.DataFrame()

    print("[Network] 中心性指標を計算中...")

    # Degree centrality
    degree = nx.degree_centrality(G)

    # Betweenness centrality
    betweenness = nx.betweenness_centrality(G, weight="weight")

    # Closeness centrality
    closeness = nx.closeness_centrality(G)

    # PageRank
    try:
        pagerank = nx.pagerank(G, weight="weight", max_iter=200)
    except Exception:
        pagerank = {n: 1.0 / G.number_of_nodes() for n in G.nodes()}

    # Eigenvector centrality
    try:
        eigenvector = nx.eigenvector_centrality(G, weight="weight", max_iter=200)
    except Exception:
        eigenvector = {n: 0.0 for n in G.nodes()}

    records = []
    for node in G.nodes():
        records.append({
            "gene_symbol": node,
            "is_gwas": G.nodes[node].get("is_gwas", False),
            "degree": G.degree(node),
            "degree_centrality": degree.get(node, 0),
            "betweenness_centrality": betweenness.get(node, 0),
            "closeness_centrality": closeness.get(node, 0),
            "pagerank": pagerank.get(node, 0),
            "eigenvector_centrality": eigenvector.get(node, 0),
        })

    df = pd.DataFrame(records).sort_values("pagerank", ascending=False).reset_index(drop=True)
    print(f"[Network] {len(df)} ノードの中心性指標を計算完了")
    return df


def propagate_flow_scores(G: nx.Graph, gene_scores_df: pd.DataFrame,
                          damping: float = 0.5, iterations: int = 5) -> dict:
    """
    GWAS 遺伝子スコアを PPI ネットワーク上でフロー伝播させる

    readme.md のコンセプト:
    - SNP関連遺伝子からの「フロー」として、GWAS P値 (-log10) + variant影響度を
      初期流量として定義
    - ネットワーク上のエッジを通じて、隣接遺伝子にスコアを伝播
    - damping factor で減衰しながら拡散

    Parameters
    ----------
    G : nx.Graph
        PPI ネットワーク
    gene_scores_df : pd.DataFrame
        遺伝子スコア (columns: gene_symbol, total_score)
    damping : float
        減衰係数 (0-1, 高いほど遠い遺伝子にも伝播)
    iterations : int
        伝播の反復回数

    Returns
    -------
    dict : {gene_symbol: flow_score}
    """
    if G.number_of_nodes() == 0:
        return {}

    print(f"[Flow] フロー伝播: damping={damping}, iterations={iterations}")

    # 初期スコア設定
    flow_scores = {}
    score_dict = {}
    if gene_scores_df is not None and not gene_scores_df.empty:
        score_dict = dict(zip(gene_scores_df["gene_symbol"], gene_scores_df["total_score"]))

    for node in G.nodes():
        if G.nodes[node].get("is_gwas", False):
            flow_scores[node] = score_dict.get(node, 1.0)
        else:
            flow_scores[node] = 0.0

    # 反復的フロー伝播
    for it in range(iterations):
        new_scores = {}
        for node in G.nodes():
            neighbors = list(G.neighbors(node))
            if not neighbors:
                new_scores[node] = flow_scores[node]
                continue

            # 隣接ノードからの流入
            incoming_flow = 0.0
            for neighbor in neighbors:
                n_neighbors = max(G.degree(neighbor), 1)
                edge_weight = G[node][neighbor].get("weight", 1.0)
                if isinstance(edge_weight, (int, float)):
                    w = edge_weight
                else:
                    w = 1.0
                incoming_flow += flow_scores[neighbor] * w / n_neighbors

            # 自分の初期スコア + 減衰した流入スコア
            initial = score_dict.get(node, 0.0) if G.nodes[node].get("is_gwas", False) else 0.0
            new_scores[node] = initial + damping * incoming_flow

        flow_scores = new_scores

    print(f"[Flow] 非ゼロスコア: {sum(1 for v in flow_scores.values() if v > 0)} / {len(flow_scores)}")
    return flow_scores


def compute_rwr_scores(G: nx.Graph, gene_scores_df: pd.DataFrame,
                       restart_prob: float = 0.3, max_iter: int = 100,
                       tol: float = 1e-6) -> dict:
    """
    Random Walk with Restart (RWR) を計算し、シグナルをネットワーク上に伝播させる。
    
    Parameters
    ----------
    G : nx.Graph
        PPI ネットワーク
    gene_scores_df : pd.DataFrame
        シードとする GWAS 遺伝子スコア
    restart_prob : float
        シードに戻る確率 (0〜1)
    max_iter : int
        最大反復回数
    tol : float
        収束判定の許容誤差
        
    Returns
    -------
    dict : {gene_symbol: rwr_score}
    """
    if G.number_of_nodes() == 0:
        return {}
        
    print(f"[Network] RWR計算中: restart_prob={restart_prob}, max_iter={max_iter}")
    
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)
    
    # 初期確率分布 p0 (シードスコア) の構築
    p0 = np.zeros(n_nodes)
    if gene_scores_df is not None and not gene_scores_df.empty:
        score_dict = dict(zip(gene_scores_df["gene_symbol"], gene_scores_df["total_score"]))
        for n in nodes:
            if G.nodes[n].get("is_gwas", False):
                p0[node_idx[n]] = score_dict.get(n, 1.0)
                
    sum_p0 = p0.sum()
    if sum_p0 > 0:
        p0 = p0 / sum_p0  # 確率分布として正規化
    else:
        p0 = np.ones(n_nodes) / n_nodes
        
    # 推移確率行列 M の構築
    A = nx.adjacency_matrix(G, weight="weight")
    deg = np.array(A.sum(axis=1)).flatten()
    deg[deg == 0] = 1.0  # ゼロ割回避
    
    from scipy import sparse
    D_inv = sparse.diags(1.0 / deg)
    M = A @ D_inv  # M_ij = A_ij / deg_j
    
    # 反復計算
    p_last = p0.copy()
    for i in range(max_iter):
        p_next = (1 - restart_prob) * M.dot(p_last) + restart_prob * p0
        err = np.linalg.norm(p_next - p_last, ord=1)
        if err < tol:
            print(f"  -> {i+1} 回の反復で収束 (err: {err:.2e})")
            p_last = p_next
            break
        p_last = p_next
    else:
        print(f"  -> 最大反復回数 {max_iter} に達しました (err: {err:.2e})")
        
    rwr_scores = {nodes[i]: p_last[i] for i in range(n_nodes)}
    return rwr_scores


def compute_integrated_scores(G: nx.Graph,
                              centrality_df: pd.DataFrame,
                              flow_scores: dict,
                              rwr_scores: dict = None,
                              weights: dict = None) -> pd.DataFrame:
    """
    中心性指標 + フロー + RWR スコアを統合した総合スコアを計算

    Parameters
    ----------
    G : nx.Graph
        PPI ネットワーク
    centrality_df : pd.DataFrame
        中心性指標テーブル
    flow_scores : dict
        単純フロー伝播スコア
    rwr_scores : dict
        Random Walk with Restart スコア
    weights : dict
        各指標の重み
        
    Returns
    -------
    pd.DataFrame : 統合スコア
    """
    if weights is None:
        weights = {
            "rwr_score": 0.40,
            "flow_score": 0.30,
            "pagerank": 0.10,
            "betweenness_centrality": 0.10,
            "degree_centrality": 0.05,
            "closeness_centrality": 0.05,
            "eigenvector_centrality": 0.00,
        }

    if centrality_df.empty:
        return pd.DataFrame()

    df = centrality_df.copy()
    df["flow_score"] = df["gene_symbol"].map(flow_scores).fillna(0)
    if rwr_scores:
        df["rwr_score"] = df["gene_symbol"].map(rwr_scores).fillna(0)
    else:
        df["rwr_score"] = 0.0

    # 正規化 (min-max)
    score_columns = ["rwr_score", "flow_score", "pagerank", "betweenness_centrality",
                     "degree_centrality", "closeness_centrality", "eigenvector_centrality"]

    for col in score_columns:
        if col in df.columns:
            col_min = df[col].min()
            col_max = df[col].max()
            if col_max > col_min:
                df[f"{col}_norm"] = (df[col] - col_min) / (col_max - col_min)
            else:
                df[f"{col}_norm"] = 0.0

    # 統合スコア計算
    df["integrated_score"] = 0.0
    for col, w in weights.items():
        norm_col = f"{col}_norm"
        if norm_col in df.columns:
            df["integrated_score"] += df[norm_col] * w

    df = df.sort_values("integrated_score", ascending=False).reset_index(drop=True)
    print(f"[Network] 統合スコア計算完了: {len(df)} 遺伝子")
    return df


def select_key_genes(integrated_scores_df: pd.DataFrame,
                     top_n: int = 150,
                     min_score_percentile: float = 60,
                     include_all_gwas: bool = True) -> pd.DataFrame:
    """
    統合スコアに基づいて重要遺伝子を選択

    Parameters
    ----------
    integrated_scores_df : pd.DataFrame
        compute_integrated_scores の出力
    top_n : int
        最大選択遺伝子数
    min_score_percentile : float
        最小スコアパーセンタイル
    include_all_gwas : bool
        全GWAS遺伝子を含めるか

    Returns
    -------
    pd.DataFrame : 重要遺伝子リスト
    """
    if integrated_scores_df.empty:
        return pd.DataFrame()

    df = integrated_scores_df.copy()

    # パーセンタイル閾値
    threshold = np.percentile(df["integrated_score"], min_score_percentile)

    # GWAS 遺伝子は常に含める
    if include_all_gwas:
        gwas_mask = df["is_gwas"] == True
        score_mask = df["integrated_score"] >= threshold
        selected = df[gwas_mask | score_mask].copy()
    else:
        selected = df[df["integrated_score"] >= threshold].copy()

    # top_n で制限
    selected = selected.sort_values("integrated_score", ascending=False).head(top_n)
    selected = selected.reset_index(drop=True)

    n_gwas = selected["is_gwas"].sum()
    n_ppi = len(selected) - n_gwas
    print(f"[Network] 重要遺伝子: {len(selected)} (GWAS: {n_gwas}, PPI: {n_ppi})")

    return selected
