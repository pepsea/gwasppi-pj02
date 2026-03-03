import json

NOTEBOOK_PATH = '/Users/yoshinorisatomi/Documents/antigravity/gwasppi/pj02/notebooks/gwas_ppi_analysis.ipynb'

with open(NOTEBOOK_PATH, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code' and 'target_fetcher.get_chembl_targets' in "".join(cell.get('source', [])):
        src = "".join(cell['source'])
        new_src = src.replace(
            "known_targets = target_fetcher.get_chembl_targets(DISEASE_TRAIT)",
            "DISEASE_EFO = 'EFO_0001360'\nDISEASE_MESH = 'D003924'\nknown_targets = target_fetcher.get_combined_targets(DISEASE_TRAIT, mesh_id=DISEASE_MESH, efo_id=DISEASE_EFO)"
        )
        cell['source'] = [new_src]

with open(NOTEBOOK_PATH, 'w') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print("Notebook patched to use combined targets!")
