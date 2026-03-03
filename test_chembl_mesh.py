import requests

def search_chembl_mesh(mesh_id="D003924"):
    # ChEMBL drug_indication API using MESH_ID
    url = f"https://www.ebi.ac.uk/chembl/api/data/drug_indication?mesh_id={mesh_id}&format=json"
    
    print(f"Querying ChEMBL for MeSH ID: {mesh_id}")
    try:
        resp = requests.get(url)
        resp.raise_for_status()
        data = resp.json()
        
        indications = data.get("drug_indications", [])
        print(f"Found {len(indications)} indications using MeSH ID.")
            
        drugs = {}
        for ind in indications:
            mol_id = ind.get("molecule_chembl_id")
            phase = ind.get("max_phase_for_ind", 0)
            if phase is None:
                phase = 0
            
            if mol_id:
                drugs[mol_id] = max(drugs.get(mol_id, 0), int(phase))
                
        print(f"Found {len(drugs)} unique drugs.")
        
        if not drugs:
            return
            
        # Get mechanisms for top 50 drugs to test
        mol_filter = ",".join(list(drugs.keys())[:50])
        mech_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id__in={mol_filter}&format=json&limit=1000"
        m_resp = requests.get(mech_url)
        m_data = m_resp.json()
        mechanisms = m_data.get("mechanisms", [])
        
        targets = {}
        for mech in mechanisms:
            t_id = mech.get("target_chembl_id")
            m_id = mech.get("molecule_chembl_id")
            if t_id and m_id:
                phase = drugs.get(m_id, 0)
                targets[t_id] = max(targets.get(t_id, 0), phase)
                
        print(f"Found {len(targets)} targets for the first 50 drugs.")
        
        # Resolve target ChEMBL IDs to Gene Symbols
        target_symbols = {}
        t_ids = list(targets.keys())
        
        chunk_size = 50
        for i in range(0, len(t_ids), chunk_size):
            chunk = t_ids[i:i+chunk_size]
            t_filter = ",".join(chunk)
            t_url = f"https://www.ebi.ac.uk/chembl/api/data/target?target_chembl_id__in={t_filter}&format=json&limit=1000"
            t_resp = requests.get(t_url)
            if t_resp.status_code == 200:
                t_data = t_resp.json()
                for t_info in t_data.get("targets", []):
                    t_id = t_info.get("target_chembl_id")
                    
                    # 遺伝子シンボルは target_components 内にある
                    components = t_info.get("target_components", [])
                    for comp in components:
                        synonyms = comp.get("target_component_synonyms", [])
                        for syn in synonyms:
                            if syn.get("syn_type") == "GENE_SYMBOL":
                                sym = syn.get("component_synonym")
                                if sym:
                                    target_symbols[sym.upper()] = targets[t_id]
                                    break
        
        print(f"Final resolved gene symbols: {len(target_symbols)}")
        print(list(target_symbols.items())[:20])
        
    except Exception as e:
        print(f"Error: {e}")

search_chembl_mesh()
