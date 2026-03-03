import re
import os
import json

filepath = '/Users/yoshinorisatomi/Documents/antigravity/gwasppi/pj02/visualization.py'
with open(filepath, 'r') as f:
    content = f.read()

# Add the missing create_pathway_gene_matrix function
matrix_func = '''
# ============================================================
# 4. パスウェイ-遺伝子マトリックス
# ============================================================
def create_pathway_gene_matrix(enrichment_df: pd.DataFrame,
                               key_genes: list = None,
                               gene_scores: dict = None,
                               top_n_pathways: int = 15,
                               save_path: str = None) -> go.Figure:
    """
    パスウェイごとに重要遺伝子をリストアップしたヒートマップ (スコア表示対応)

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        エンリッチメント結果
    key_genes : list
        重要遺伝子リスト
    gene_scores : dict
        遺伝子のスコア辞書 {gene: score}
    top_n_pathways : int
        表示パスウェイ数
    save_path : str
        保存パス
    """
    if enrichment_df.empty:
        return go.Figure()

    if key_genes is None:
        key_genes = []
    if gene_scores is None:
        gene_scores = {}
        
    key_set = set(g.upper() for g in key_genes)

    # 上位パスウェイ
    top_pw = enrichment_df.head(top_n_pathways)

    # マトリックス構築
    all_genes = set()
    pw_gene_dict = {}
    for _, row in top_pw.iterrows():
        pw_name = row["term_name"][:50]
        genes = set(g.strip().upper() for g in row["genes"].split(","))
        relevant = genes & key_set if key_set else genes
        pw_gene_dict[pw_name] = relevant
        all_genes.update(relevant)

    if not all_genes:
        return go.Figure()

    gene_list = sorted(all_genes)
    pw_list = list(pw_gene_dict.keys())

    # マトリックス (1/0 ではなくスコアを入れる)
    matrix = []
    text_matrix = []
    for pw in pw_list:
        row = []
        text_row = []
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

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=gene_list,
        y=pw_list,
        text=text_matrix,
        texttemplate="%{text}",
        colorscale="Blues",
        showscale=True,
        colorbar=dict(title="Gene Score"),
        hovertemplate="Pathway: %{y}<br>Gene: %{x}<br>Score: %{z:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="Pathway-Gene Score Matrix (Key Genes × Top Pathways)",
        xaxis_title="Gene",
        yaxis_title="Pathway",
        height=max(400, top_n_pathways * 30),
        width=max(800, len(gene_list) * 25),
        xaxis=dict(tickangle=45),
    )

    if save_path:
        fig.write_html(save_path)
        print(f"[Viz] パスウェイ-遺伝子マトリックス保存: {save_path}")

    return fig

# ============================================================
# 5. 遺伝子スコア可視化
# ============================================================
'''

# The previous patch accidentally removed the pathway matrix function entirely!
# Let's add it back right before the gene score plot # 5. 遺伝子スコア可視化
new_content = content.replace('# ============================================================\n# 5. 遺伝子スコア可視化', matrix_func)

with open(filepath, 'w') as f:
    f.write(new_content)

print("Matrix function restored")
