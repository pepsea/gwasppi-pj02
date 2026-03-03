import json
import os

NOTEBOOK_PATH = '/Users/yoshinorisatomi/Documents/antigravity/gwasppi/pj02/notebooks/gwas_ppi_analysis.ipynb'

with open(NOTEBOOK_PATH, 'r') as f:
    nb = json.load(f)

# 1. ターゲット取得セルを作成
target_md_cell = {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 4. 創薬ターゲットの取得 (機能拡張)\n",
    "OpenTargets/ChEMBL情報を利用し、指定疾患に対する既存の創薬ターゲット（承認済・開発中）を取得します。\n",
    "※現在API接続制約のため `target_fetcher.py` 内のフォールバック・データセットを使用します。"
   ]
}
target_code_cell = {
   "cell_type": "code",
   "execution_count": None,
   "metadata": {},
   "outputs": [],
   "source": [
    "import target_fetcher\n",
    "\n",
    "# 疾患名に関連するターゲットを取得\n",
    "known_targets = target_fetcher.get_chembl_targets(DISEASE_TRAIT)\n",
    "target_fetcher.save_targets(known_targets)\n",
    "\n",
    "print(f\"取得したターゲット数: {len(known_targets)}\")\n",
    "if known_targets:\n",
    "    display(pd.DataFrame(list(known_targets.items()), columns=['Gene', 'Max Phase']))"
   ]
}

# 2. パスウェイ別ネットワーク可視化セルを作成
pathway_net_md_cell = {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 10. パスウェイ別ネットワーク可視化\n",
    "各エンリッチメントパスウェイに含まれる遺伝子だけを抽出し、PPIサブネットワークを描画します。"
   ]
}
pathway_net_code_cell = {
   "cell_type": "code",
   "execution_count": None,
   "metadata": {},
   "outputs": [],
   "source": [
    "top_pw = enrichment_df.head(5)  # 上位5つのパスウェイを描画\n",
    "\n",
    "for i, row in top_pw.iterrows():\n",
    "    pw_name = row['term_name']\n",
    "    pw_genes = [g.strip().upper() for g in row['genes'].split(',')]\n",
    "    \n",
    "    safe_name = \"\".join([c if c.isalnum() else \"_\" for c in pw_name[:30]])\n",
    "    save_path = os.path.join(config.OUTPUT_DIR, f'pathway_net_{i+1}_{safe_name}.html')\n",
    "    \n",
    "    fig_pw = visualization.create_pathway_network_plot(\n",
    "        G=ppi_graph,\n",
    "        pathway_name=pw_name,\n",
    "        pathway_genes=pw_genes,\n",
    "        gene_scores=flow_scores,\n",
    "        known_targets=known_targets,\n",
    "        save_path=save_path\n",
    "    )\n",
    "    if len(fig_pw.data) > 0:\n",
    "        fig_pw.show()"
   ]
}

# 変更適用
new_cells = []
for cell in nb['cells']:
    # gwas_fetcher 呼び出しの修正 (max_studies 追加)
    if cell['cell_type'] == 'code' and 'gwas_fetcher.fetch_snps_for_disease' in "".join(cell.get('source', [])):
        src = "".join(cell['source'])
        if 'max_studies=config.DEFAULT_MAX_GWAS_STUDIES' not in src:
            src = src.replace("p_value_threshold=P_VALUE_THRESHOLD", 
                              "max_studies=config.DEFAULT_MAX_GWAS_STUDIES, p_value_threshold=P_VALUE_THRESHOLD")
            cell['source'] = [src]
            
    # 設定セルにインポート追加
    if cell['cell_type'] == 'code' and 'import config' in "".join(cell.get('source', [])):
        src = "".join(cell['source'])
        if 'import target_fetcher' not in src:
            cell['source'] = [src + "\nimport target_fetcher\n"]

    new_cells.append(cell)
    
    # "3. 遺伝子スコアリング" の前に ターゲット取得を差し込む (元々の3の前)
    if cell['cell_type'] == 'markdown' and '## 3. 遺伝子スコアリング' in "".join(cell.get('source', [])):
        # すでに挿入されていなければ挿入
        new_cells.insert(-1, target_md_cell)
        new_cells.insert(-1, target_code_cell)
        
    # create_network_plot の修正
    if cell['cell_type'] == 'code' and 'create_network_plot' in "".join(cell.get('source', [])):
        src = "".join(cell['source'])
        if 'known_targets=' not in src:
            src = src.replace("gene_scores=flow_scores,", "gene_scores=flow_scores,\n        known_targets=known_targets,")
            cell['source'] = [src]

    # create_pathway_gene_matrix の修正
    if cell['cell_type'] == 'code' and 'create_pathway_gene_matrix' in "".join(cell.get('source', [])):
        src = "".join(cell['source'])
        if 'known_targets=' not in src:
            src = src.replace("gene_scores=flow_scores,", "gene_scores=flow_scores,\n        known_targets=known_targets,")
            cell['source'] = [src]

# 最後のセルの後に パスウェイ別ネットワーク可視化 を追加
# 重複追加を防ぐ
has_pathway_net = any('10. パスウェイ別ネットワーク可視化' in "".join(c.get('source', [])) for c in new_cells)
if not has_pathway_net:
    new_cells.append(pathway_net_md_cell)
    new_cells.append(pathway_net_code_cell)

nb['cells'] = new_cells

with open(NOTEBOOK_PATH, 'w') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print("Notebook patched successfully!")
