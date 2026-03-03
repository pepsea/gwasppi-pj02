import requests
import json

def search_chembl_targets(efo_id):
    # EFO ID format: EFO_0001360 for Type 2 Diabetes
    print(f"Querying ChEMBL for EFO ID: {efo_id}")
    
    # 1. 疾患 (indication) から薬剤を引く
    # https://www.ebi.ac.uk/chembl/api/data/drug_indication?efo_id=EFO_0001360&format=json
    url = f"https://www.ebi.ac.uk/chembl/api/data/drug_indication?efo_id={efo_id}&format=json&limit=1000"
    
    try:
        resp = requests.get(url)
        resp.raise_for_status()
        data = resp.json()
        
        indications = data.get("drug_indications", [])
        print(f"Found {len(indications)} indications.")
        
        # 薬剤のchembl_idとmax_phaseを取得
        drugs = {}
        for ind in indications:
            molecule_id = ind.get("molecule_chembl_id")
            phase = ind.get("max_phase_for_ind", 0)
            if molecule_id:
                drugs[molecule_id] = max(drugs.get(molecule_id, 0), phase)
                
        print(f"Found {len(drugs)} unique drugs.")
        
        # 2. 薬剤の作用機序 (mechanism) からターゲットを取得
        # Note: bulk querying mechanism is better but let's do a grouped filter
        targets = {}
        
        # ChEMBL mechanism API allows filtering by molecule_chembl_id
        # To avoid 100+ requests, we can chunk them
        chembl_ids = list(drugs.keys())
        chunk_size = 50
        
        for i in range(0, len(chembl_ids), chunk_size):
            chunk = chembl_ids[i:i+chunk_size]
            mol_filter = ",".join(chunk)
            
            mech_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id__in={mol_filter}&format=json&limit=1000"
            m_resp = requests.get(mech_url)
            if m_resp.status_code == 200:
                m_data = m_resp.json()
                mechanisms = m_data.get("mechanisms", [])
                
                for mech in mechanisms:
                    target_chembl_id = mech.get("target_chembl_id")
                    mol_id = mech.get("molecule_chembl_id")
                    
                    if target_chembl_id and mol_id:
                        phase = drugs.get(mol_id, 0)
                        targets[target_chembl_id] = max(targets.get(target_chembl_id, 0), phase)
                        
        print(f"Found {len(targets)} unique target ChEMBL IDs.")
        
        # 3. ターゲット ChEMBL ID から Gene Symbol を取得
        target_symbols = {}
        t_ids = list(targets.keys())
        
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

search_chembl_targets("EFO_0001360")
