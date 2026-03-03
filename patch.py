import re
import os
import json

filepath = '/Users/yoshinorisatomi/Documents/antigravity/gwasppi/pj02/visualization.py'
with open(filepath, 'r') as f:
    content = f.read()

# Replace create_network_plot
new_func = '''def create_network_plot(G: nx.Graph,
                        key_genes: list = None,
                        gene_scores: dict = None,
                        save_path: str = None) -> go.Figure:
    """
    Plotly による PPI ネットワーク可視化 (スプリングレイアウト・中心配置)
    """
    if key_genes is None:
        key_genes = []
    if gene_scores is None:
        gene_scores = {}

    key_gene_set = set(g.upper() for g in key_genes)

    # ---- 散開レイアウト (Spring Layout) ----
    n_nodes = G.number_of_nodes()
    # 広めに配置するためkを大きくする
    if n_nodes > 200:
        k_val = 2.5 / (n_nodes ** 0.5)
        iters = 50
    else:
        k_val = 4.0 / max(n_nodes ** 0.5, 1)
        iters = 100

    pos = nx.spring_layout(G, k=k_val, iterations=iters, seed=42, weight='weight')

    # ---- エッジを GWAS関連 と PPI-PPI に分離 ----
    gwas_set = set(n for n in G.nodes() if G.nodes[n].get("is_gwas", False))

    gwas_edge_x, gwas_edge_y = [], []
    ppi_edge_x, ppi_edge_y = [], []

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
        edge_traces.append(go.Scatter(
            x=gwas_edge_x, y=gwas_edge_y,
            line=dict(width=1.0, color="rgba(29,53,87,0.3)"),
            hoverinfo="none", mode="lines", name="GWAS ↔ PPI"
        ))
    if ppi_edge_x:
        edge_traces.append(go.Scatter(
            x=ppi_edge_x, y=ppi_edge_y,
            line=dict(width=0.8, color="rgba(231,111,81,0.3)"),
            hoverinfo="none", mode="lines", name="PPI ↔ PPI"
        ))

    # ---- ノード ----
    gwas_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": []}
    ppi_key_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": []}
    ppi_other_nodes = {"x": [], "y": [], "text": [], "hover": [], "size": [], "color": []}

    for node in G.nodes():
        x, y = pos[node]
        score = gene_scores.get(node, 0)
        is_gwas = node in gwas_set
        is_key = node in key_gene_set

        hover = f"{node}<br>Score: {score:.2f}<br>Degree: {G.degree(node)}"
        if is_gwas: hover += "<br>Type: GWAS (SNP関連)"
        elif is_key: hover += "<br>Type: Key PPI遺伝子"
        else: hover += "<br>Type: PPI遺伝子"

        label = node if (is_key or is_gwas) else ""
        base_size = 18 if is_key else 16 if is_gwas else 8
        if score > 0:
            base_size = max(base_size, min(30, 8 + score * 2))

        if is_gwas:
            gwas_nodes["x"].append(x); gwas_nodes["y"].append(y); gwas_nodes["text"].append(label)
            gwas_nodes["hover"].append(hover); gwas_nodes["size"].append(base_size)
            gwas_nodes["color"].append("#E63946" if is_key else "#1D3557")
        elif is_key:
            ppi_key_nodes["x"].append(x); ppi_key_nodes["y"].append(y); ppi_key_nodes["text"].append(label)
            ppi_key_nodes["hover"].append(hover); ppi_key_nodes["size"].append(base_size)
            ppi_key_nodes["color"].append("#E76F51")
        else:
            ppi_other_nodes["x"].append(x); ppi_other_nodes["y"].append(y); ppi_other_nodes["text"].append(label)
            ppi_other_nodes["hover"].append(hover); ppi_other_nodes["size"].append(base_size)
            ppi_other_nodes["color"].append("#A8DADC")

    traces = [*edge_traces]

    if ppi_other_nodes["x"]:
        traces.append(go.Scatter(x=ppi_other_nodes["x"], y=ppi_other_nodes["y"], mode="markers",
            hoverinfo="text", hovertext=ppi_other_nodes["hover"], name="PPI遺伝子",
            marker=dict(symbol="circle", size=ppi_other_nodes["size"], color=ppi_other_nodes["color"], line=dict(width=0.5, color="grey"))))
    
    if ppi_key_nodes["x"]:
        traces.append(go.Scatter(x=ppi_key_nodes["x"], y=ppi_key_nodes["y"], mode="markers+text",
            hoverinfo="text", text=ppi_key_nodes["text"], textposition="top center", textfont=dict(size=9, color="#E76F51"),
            hovertext=ppi_key_nodes["hover"], name="Key PPI遺伝子 ★",
            marker=dict(symbol="star", size=[s*1.2 for s in ppi_key_nodes["size"]], color="#E76F51", line=dict(width=1, color="white"))))
            
    if gwas_nodes["x"]:
        traces.append(go.Scatter(x=gwas_nodes["x"], y=gwas_nodes["y"], mode="markers+text",
            hoverinfo="text", text=gwas_nodes["text"], textposition="top center", textfont=dict(size=9, color="#1D3557"),
            hovertext=gwas_nodes["hover"], name="GWAS遺伝子 ◆",
            marker=dict(symbol="diamond", size=gwas_nodes["size"], color="#1D3557", line=dict(width=1.5, color="white"))))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title="PPI Network (Spring Layout)",
        showlegend=True, hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
        width=1100, height=1000, plot_bgcolor="white",
        legend=dict(x=1.02, y=1, bordercolor="black", borderwidth=1, font=dict(size=12))
    )

    if save_path:
        fig.write_html(save_path)
    return fig
'''

pattern = r'def create_network_plot\(.*?return fig'
new_content = re.sub(pattern, new_func, content, flags=re.DOTALL)

with open(filepath, 'w') as f:
    f.write(new_content)

# Update Notebook to pass flow_scores to pathway matrix
nb_path = '/Users/yoshinorisatomi/Documents/antigravity/gwasppi/pj02/notebooks/gwas_ppi_analysis.ipynb'
with open(nb_path, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    src = ''.join(cell.get('source', []))
    if 'fig_matrix = visualization.create_pathway_gene_matrix' in src:
        new_src = src.replace(
            'key_genes=key_gene_list,',
            'key_genes=key_gene_list,\n    gene_scores=flow_scores,'
        )
        cell['source'] = [line + '\n' for line in new_src.split('\n')][:-1]

with open(nb_path, 'w') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print("Patch applied")
