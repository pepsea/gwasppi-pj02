import target_fetcher

disease = "type 2 diabetes mellitus"
efo = "EFO_0001360"
mesh = "D003924"

targets = target_fetcher.get_combined_targets(disease, mesh_id=mesh, efo_id=efo)
print(f"Total targets: {len(targets)}")
print(list(targets.items())[:10])
