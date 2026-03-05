"""
Visualization Module
====================
Plotly による可視化:
1. サンキーダイアグラム (SNP → SNP関連遺伝子 → PPI重要遺伝子)
2. ネットワーク図 (Plotly)
3. エンリッチメント解析結果の図示
4. パスウェイ-遺伝子マトリックス
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx

import config


# ============================================================
# 1. サンキーダイアグラム
# ============================================================
def create_sankey_diagram(gwas_df: pd.DataFrame,
                          ppi_df: pd.DataFrame = None,
                          key_genes_df: pd.DataFrame = None,
                          enrichment_df: pd.DataFrame = None,
                          top_n_pathways: int = 15,
                          save_path: str = None) -> go.Figure:
    """
    SNP → GWAS関連遺伝子 → PPI重要遺伝子 のサンキーダイアグラム

    Parameters
    ----------
    gwas_df : pd.DataFrame
        GWAS SNP データ
    ppi_df : pd.DataFrame
        PPI インタラクション
    key_genes_df : pd.DataFrame
        重要遺伝子リスト
    enrichment_df : pd.DataFrame
        (未使用、互換性のため保持)
    top_n_pathways : int
        (未使用、互換性のため保持)
    save_path : str
        保存パス

    Returns
    -------
    go.Figure : Plotly Figure
    """
    labels = []
    sources = []
    targets = []
    values = []
    colors = []

    node_index = {}

    def get_node_idx(name):
        if name not in node_index:
            node_index[name] = len(labels)
            labels.append(name)
        return node_index[name]

    # Layer 1: SNP → GWAS遺伝子
    gwas_genes = set()
    if not gwas_df.empty:
        top_snps = gwas_df.sort_values("p_value").head(50)
        for _, row in top_snps.iterrows():
            snp = row["rsid"]
            gene = row["gene_symbol"]
            gwas_genes.add(gene)

            src = get_node_idx(f"SNP:{snp}")
            tgt = get_node_idx(gene)
            sources.append(src)
            targets.append(tgt)
            values.append(max(1, -np.log10(max(row["p_value"], 1e-300))))
            colors.append("rgba(31, 119, 180, 0.4)")

    # Layer 2: GWAS遺伝子 → PPI重要遺伝子
    key_gene_set = set()
    if key_genes_df is not None and not key_genes_df.empty:
        key_gene_set = set(key_genes_df["gene_symbol"])

    if ppi_df is not None and not ppi_df.empty:
        for _, row in ppi_df.iterrows():
            gene_a = row["gene_a"]
            gene_b = row["gene_b"]

            # GWAS遺伝子 → 重要PPI遺伝子のフロー
            if gene_a in gwas_genes and gene_b in key_gene_set and gene_b not in gwas_genes:
                src = get_node_idx(gene_a)
                tgt = get_node_idx(gene_b)
                sources.append(src)
                targets.append(tgt)
                values.append(1)
                colors.append("rgba(255, 127, 14, 0.4)")
            elif gene_b in gwas_genes and gene_a in key_gene_set and gene_a not in gwas_genes:
                src = get_node_idx(gene_b)
                tgt = get_node_idx(gene_a)
                sources.append(src)
                targets.append(tgt)
                values.append(1)
                colors.append("rgba(255, 127, 14, 0.4)")

    # ノードの色分け
    node_colors = []
    for label in labels:
        if label.startswith("SNP:"):
            node_colors.append("#1f77b4")   # 青: SNP
        elif label in gwas_genes:
            node_colors.append("#ff7f0e")   # オレンジ: GWAS遺伝子
        else:
            node_colors.append("#d62728")   # 赤: PPI重要遺伝子

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=node_colors,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=colors,
        ),
    )])

    fig.update_layout(
        title="GWAS Flow: SNP → 関連遺伝子 → PPI重要遺伝子",
        font_size=11,
        height=max(500, len(labels) * 12),
        width=1000,
        margin=dict(l=20, r=20, t=40, b=20),
    )

    if save_path:
        fig.write_html(save_path)
        print(f"[Viz] サンキーダイアグラム保存: {save_path}")

    return fig


# ============================================================
# 2. ネットワーク図 (Plotly) — フローベースレイアウト
# ============================================================
def create_network_plot(G: nx.Graph,
                        key_genes: list = None,
                        gene_scores: dict = None,
                        known_targets: dict = None,
                        save_path: str = None) -> go.Figure:
    """
    Plotly による PPI ネットワーク可視化 (ターゲット強調)
    """
    if key_genes is None: key_genes = []
    if gene_scores is None: gene_scores = {}
    if known_targets is None: known_targets = {}

    key_gene_set = set(g.upper() for g in key_genes)

    n_nodes = G.number_of_nodes()
    if n_nodes > 200:
        k_val, iters = 2.5 / (n_nodes ** 0.5), 50
    else:
        k_val, iters = 4.0 / max(n_nodes ** 0.5, 1), 100

    pos = nx.spring_layout(G, k=k_val, iterations=iters, seed=42, weight='weight')
    gwas_set = set(n for n in G.nodes() if G.nodes[n].get("is_gwas", False))

    gwas_edge_x, gwas_edge_y, ppi_edge_x, ppi_edge_y = [], [], [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        if u in gwas_set or v in gwas_set:
            gwas_edge_x.extend([x0, x1, None])
            gwas_edge_y.extend([y0, y1, None])
        else:
            ppi_edge_x.extend([x0, x1, None])
            ppi_edge_y.extend([y0, y1, None])

    edge_traces = []
    if gwas_edge_x:
        edge_traces.append(go.Scatter(x=gwas_edge_x, y=gwas_edge_y, line=dict(width=1.0, color="rgba(29,53,87,0.3)"), hoverinfo="none", mode="lines", name="GWAS ↔ PPI"))
    if ppi_edge_x:
        edge_traces.append(go.Scatter(x=ppi_edge_x, y=ppi_edge_y, line=dict(width=0.8, color="rgba(231,111,81,0.3)"), hoverinfo="none", mode="lines", name="PPI ↔ PPI"))

    gwas_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": [], "line_color": [], "line_width": []}
    ppi_key_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": [], "line_color": [], "line_width": []}
    ppi_other_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": [], "line_color": [], "line_width": []}

    for node in G.nodes():
        x, y = pos[node]
        score = gene_scores.get(node, 0)
        is_gwas, is_key, is_target = node in gwas_set, node in key_gene_set, node in known_targets
        target_phase = known_targets.get(node, "")

        hover = f"{node}<br>Score: {score:.2f}<br>Degree: {G.degree(node)}"
        if is_target: hover += f"<br><b>💊 Drug Target (Phase {target_phase})</b>"
        if is_gwas: hover += "<br>Type: GWAS (SNP関連)"
        elif is_key: hover += "<br>Type: Key PPI遺伝子"

        label = f"💊{node}" if is_target else (node if (is_key or is_gwas) else "")
        base_size = 18 if is_key else 16 if is_gwas else 8
        if score > 0: base_size = max(base_size, min(30, 8 + score * 2))
        if is_target: base_size = max(base_size, 22)

        l_color = "#FF00FF" if is_target else ("white" if (is_gwas or is_key) else "grey")
        l_width = 3 if is_target else (1.5 if is_gwas else (1 if is_key else 0.5))

        if is_gwas:
            target_dict = gwas_nodes
            target_dict["color"].append("#E63946" if is_key else "#1D3557")
        elif is_key:
            target_dict = ppi_key_nodes
            target_dict["color"].append("#E76F51")
        else:
            target_dict = ppi_other_nodes
            target_dict["color"].append("#A8DADC")

        target_dict["x"].append(x); target_dict["y"].append(y); target_dict["text"].append(label)
        target_dict["hover"].append(hover); target_dict["size"].append(base_size)
        target_dict["line_color"].append(l_color); target_dict["line_width"].append(l_width)

    traces = [*edge_traces]
    if ppi_other_nodes["x"]:
        traces.append(go.Scatter(x=ppi_other_nodes["x"], y=ppi_other_nodes["y"], mode="markers+text", hoverinfo="text", text=ppi_other_nodes["text"], textposition="top center", textfont=dict(size=10, color="#FF00FF"), hovertext=ppi_other_nodes["hover"], name="PPI遺伝子", marker=dict(symbol="circle", size=ppi_other_nodes["size"], color=ppi_other_nodes["color"], line=dict(width=ppi_other_nodes["line_width"], color=ppi_other_nodes["line_color"]))))
    if ppi_key_nodes["x"]:
        traces.append(go.Scatter(x=ppi_key_nodes["x"], y=ppi_key_nodes["y"], mode="markers+text", hoverinfo="text", text=ppi_key_nodes["text"], textposition="top center", textfont=dict(size=10, color=ppi_key_nodes["line_color"]), hovertext=ppi_key_nodes["hover"], name="Key PPI遺伝子 ★", marker=dict(symbol="star", size=[s*1.2 for s in ppi_key_nodes["size"]], color="#E76F51", line=dict(width=ppi_key_nodes["line_width"], color=ppi_key_nodes["line_color"]))))
    if gwas_nodes["x"]:
        traces.append(go.Scatter(x=gwas_nodes["x"], y=gwas_nodes["y"], mode="markers+text", hoverinfo="text", text=gwas_nodes["text"], textposition="top center", textfont=dict(size=10, color=gwas_nodes["line_color"]), hovertext=gwas_nodes["hover"], name="GWAS遺伝子 ◆", marker=dict(symbol="diamond", size=gwas_nodes["size"], color="#1D3557", line=dict(width=gwas_nodes["line_width"], color=gwas_nodes["line_color"]))))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title="PPI Network with Targets", showlegend=True, hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
        width=1100, height=1000, plot_bgcolor="white",
        legend=dict(x=1.02, y=1, bordercolor="black", borderwidth=1, font=dict(size=12)),
        annotations=[dict(x=0.01, y=1.02, xref="paper", yref="paper", text="💊 Drug Targets (Magenta Outlines)", showarrow=False, font=dict(size=12, color="#FF00FF"))]
    )
    if save_path: fig.write_html(save_path)
    return fig

def create_pathway_network_plot(G: nx.Graph,
                                pathway_name: str,
                                pathway_genes: list,
                                gene_scores: dict = None,
                                known_targets: dict = None,
                                save_path: str = None) -> go.Figure:
    """
    指定されたパスウェイに含まれる遺伝子だけを抽出してPPIサブグラフを描画
    (スコアで色分け、ターゲットを強調)
    """
    if gene_scores is None: gene_scores = {}
    if known_targets is None: known_targets = {}
    
    pathway_nodes = set(g.upper() for g in pathway_genes) & set(G.nodes())
    if not pathway_nodes: return go.Figure()
        
    subG = G.subgraph(pathway_nodes)
    pos = nx.spring_layout(subG, k=2.0/max(len(pathway_nodes)**0.5, 1), iterations=50, seed=42)
    
    edge_x, edge_y = [], []
    for u, v in subG.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        
    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=1.0, color="#888"), hoverinfo="none", mode="lines", name="PPI")
    
    node_x, node_y, node_text, node_hover, node_color, node_size, line_color, line_width = [], [], [], [], [], [], [], []
    for node in subG.nodes():
        node_x.append(pos[node][0]); node_y.append(pos[node][1])
        score = gene_scores.get(node, 0)
        is_target = node in known_targets
        
        node_text.append(f"💊{node}" if is_target else node)
        
        hover = f"{node}<br>Score: {score:.2f}"
        if is_target: hover += f"<br><b>💊 Drug Target (Phase {known_targets.get(node)})</b>"
        node_hover.append(hover)
        
        node_color.append(score)
        node_size.append(max(10, min(30, 8 + score * 2)))
        line_color.append("#FF00FF" if is_target else ("#E63946" if subG.nodes[node].get("is_gwas", False) else "white"))
        line_width.append(3 if is_target else (2 if subG.nodes[node].get("is_gwas", False) else 0.5))
        
    node_trace = go.Scatter(
        x=node_x, y=node_y, mode="markers+text", text=node_text, textposition="top center", hoverinfo="text", hovertext=node_hover,
        marker=dict(showscale=True, colorscale="Viridis", color=node_color, size=node_size, colorbar=dict(thickness=15, title="Gene Score"), line=dict(width=line_width, color=line_color))
    )
    
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title=f"Pathway Network: {pathway_name[:50]}...", showlegend=False, 
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        width=800, height=800, plot_bgcolor="white"
    )
    if save_path: fig.write_html(save_path)
    return fig




# ============================================================
# 3. エンリッチメント可視化
# ============================================================
def create_enrichment_bubble_chart(enrichment_df: pd.DataFrame,
                                   top_n: int = 20,
                                   save_path: str = None) -> go.Figure:
    """エンリッチメント結果のバブルチャート"""
    if enrichment_df.empty:
        return go.Figure()

    df = enrichment_df.head(top_n).copy()
    df["neg_log_fdr"] = -np.log10(df["fdr"].clip(lower=1e-50))

    # データベース別に色分け
    color_map = {"Reactome": "#1f77b4", "GO": "#2ca02c", "HPO": "#d62728"}
    df["color"] = df["database"].map(color_map).fillna("#7f7f7f")

    fig = go.Figure()

    for db_name in df["database"].unique():
        db_df = df[df["database"] == db_name]
        fig.add_trace(go.Scatter(
            x=db_df["fold_enrichment"],
            y=db_df["term_name"].str[:50],
            mode="markers",
            name=db_name,
            marker=dict(
                size=db_df["gene_count"] * 3 + 5,
                color=db_df["neg_log_fdr"],
                colorscale="Viridis",
                showscale=True if db_name == df["database"].unique()[0] else False,
                colorbar=dict(title="-log10(FDR)"),
                line=dict(width=1, color="white"),
            ),
            text=[f"{row['term_name']}<br>FDR: {row['fdr']:.2e}<br>"
                  f"Genes: {row['gene_count']}/{row['total_genes']}<br>"
                  f"Fold: {row['fold_enrichment']:.1f}x"
                  for _, row in db_df.iterrows()],
            hoverinfo="text",
        ))

    fig.update_layout(
        title="Enrichment Analysis Results",
        xaxis_title="Fold Enrichment",
        yaxis_title="",
        height=max(500, top_n * 30),
        width=900,
        yaxis=dict(autorange="reversed"),
    )

    if save_path:
        fig.write_html(save_path)
        print(f"[Viz] バブルチャート保存: {save_path}")

    return fig


def create_enrichment_barplot(enrichment_df: pd.DataFrame,
                               top_n: int = 20,
                               save_path: str = None) -> go.Figure:
    """エンリッチメント結果の棒グラフ"""
    if enrichment_df.empty:
        return go.Figure()

    df = enrichment_df.head(top_n).copy()
    df["neg_log_fdr"] = -np.log10(df["fdr"].clip(lower=1e-50))

    color_map = {"Reactome": "#1f77b4", "GO": "#2ca02c", "HPO": "#d62728"}
    df["color"] = df["database"].map(color_map).fillna("#7f7f7f")

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df["neg_log_fdr"],
        y=df["term_name"].str[:50],
        orientation="h",
        marker_color=df["color"],
        text=[f"{row['database']} | Genes: {row['gene_count']}" for _, row in df.iterrows()],
        textposition="auto",
    ))

    fig.update_layout(
        title="Enrichment Analysis (-log10 FDR)",
        xaxis_title="-log10(FDR)",
        yaxis_title="",
        height=max(500, top_n * 30),
        width=900,
        yaxis=dict(autorange="reversed"),
    )

    if save_path:
        fig.write_html(save_path)
        print(f"[Viz] 棒グラフ保存: {save_path}")

    return fig


# ============================================================
# 4. パスウェイ-遺伝子マトリックス
# ============================================================
def create_pathway_gene_matrix(enrichment_df: pd.DataFrame,
                               key_genes: list = None,
                               gene_scores: dict = None,
                               known_targets: dict = None,
                               top_n_pathways: int = 15,
                               save_path: str = None) -> go.Figure:
    """
    パスウェイごとに重要遺伝子をリストアップしたヒートマップ (ターゲット強調付き)
    """
    if enrichment_df.empty: return go.Figure()
    if key_genes is None: key_genes = []
    if gene_scores is None: gene_scores = {}
    if known_targets is None: known_targets = {}
        
    key_set = set(g.upper() for g in key_genes)
    top_pw = enrichment_df.head(top_n_pathways)

    all_genes = set()
    pw_gene_dict = {}
    for _, row in top_pw.iterrows():
        pw_name = row["term_name"][:50]
        genes = set(g.strip().upper() for g in row["genes"].split(","))
        relevant = genes & key_set if key_set else genes
        pw_gene_dict[pw_name] = relevant
        all_genes.update(relevant)

    if not all_genes: return go.Figure()

    gene_list = sorted(all_genes)
    pw_list = list(pw_gene_dict.keys())

    matrix, text_matrix = [], []
    for pw in pw_list:
        row, text_row = [], []
        for gene in gene_list:
            if gene in pw_gene_dict[pw]:
                score = gene_scores.get(gene, 1.0)
                row.append(score)
                text_row.append(f"{score:.2f}")
            else:
                row.append(0)
                text_row.append("")
        matrix.append(row)
        text_matrix.append(text_row)

    x_labels = [f"💊{g} (Ph{known_targets[g]})" if g in known_targets else g for g in gene_list]

    fig = go.Figure(data=go.Heatmap(
        z=matrix, x=x_labels, y=pw_list, text=text_matrix, texttemplate="%{text}",
        colorscale="Blues", showscale=True, colorbar=dict(title="Gene Score"),
        hovertemplate="Pathway: %{y}<br>Gene: %{x}<br>Score: %{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="Pathway-Gene Score Matrix (💊 = Known Target)",
        xaxis_title="Gene", yaxis_title="Pathway",
        height=max(400, top_n_pathways * 30), width=max(800, len(gene_list) * 30),
        xaxis=dict(tickangle=45),
    )

    if save_path:
        fig.write_html(save_path)
    return fig

# ============================================================
# 5. 遺伝子スコア可視化
# ============================================================

# ============================================================
def create_gene_score_plot(gene_scores_df: pd.DataFrame,
                           top_n: int = 30,
                           save_path: str = None) -> go.Figure:
    """遺伝子スコアの棒グラフ (スコア内訳表示)"""
    if gene_scores_df.empty:
        return go.Figure()

    df = gene_scores_df.head(top_n).copy()

    score_components = ["pval_score", "lof_score", "gof_score", "missense_score",
                        "coding_score", "regulatory_score", "ppi_score"]
    component_colors = {
        "pval_score": "#1f77b4",
        "lof_score": "#d62728",
        "gof_score": "#ff7f0e",
        "missense_score": "#2ca02c",
        "coding_score": "#9467bd",
        "regulatory_score": "#8c564b",
        "ppi_score": "#e377c2",
    }
    component_labels = {
        "pval_score": "GWAS P値",
        "lof_score": "LoF",
        "gof_score": "GoF",
        "missense_score": "Missense",
        "coding_score": "Coding",
        "regulatory_score": "Regulatory",
        "ppi_score": "PPI Multi-Source",
    }

    fig = go.Figure()

    for comp in score_components:
        if comp in df.columns:
            fig.add_trace(go.Bar(
                name=component_labels.get(comp, comp),
                x=df["gene_symbol"],
                y=df[comp],
                marker_color=component_colors.get(comp, "#7f7f7f"),
            ))

    fig.update_layout(
        barmode="stack",
        title="Gene Scores (Top Genes)",
        xaxis_title="Gene",
        yaxis_title="Score",
        height=500,
        width=max(800, top_n * 25),
        xaxis=dict(tickangle=45),
    )

    if save_path:
        fig.write_html(save_path)
        print(f"[Viz] 遺伝子スコア図保存: {save_path}")

    return fig
