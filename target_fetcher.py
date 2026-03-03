"""
Target Fetcher Module
=====================
疾患に対する既知の創薬ターゲット（承認済、または臨床開発中）を取得する機能。
OpenTargets APIが利用できない環境に対応するため、Pharos API と ChEMBL API の
両方を組み合わせてターゲット遺伝子を抽出し、カバー率を高めます。

要件:
- 指定された疾患名 (および EFO ID, MeSH ID 等) に関連するターゲットを抽出
- Pharos GraphQL API を使用して Tclin/Tchem クラスのターゲットを取得
- ChEMBL API を使用して MeSH/EFO または疾患文字列に合致する適応症から薬剤を辿り、ターゲットを取得
"""

import requests
import os
import json
import time

def get_pharos_targets(disease_name: str) -> dict:
    """
    Pharos GraphQL APIを使用して、疾患名からTclin(承認/臨床段階)ターゲットを取得
    """
    url = "https://pharos-api.ncats.io/graphql"
    query = """
    query searchDiseases($name: String!) {
      diseases(filter: {term: $name}) {
        diseases {
          name
          targets {
            sym
            tdl
          }
        }
      }
    }
    """
    variables = {"name": disease_name}
    targets = {}

    print(f"[Target] Pharos APIを検索中: {disease_name}")
    try:
        response = requests.post(url, json={'query': query, 'variables': variables}, timeout=30)
        if response.status_code == 200:
            data = response.json()
            diseases = data.get("data", {}).get("diseases", {}).get("diseases", [])
            
            if diseases:
                # 最初の一致を使用
                found_targets = diseases[0].get("targets", [])
                for t in found_targets:
                    sym = t.get("sym")
                    tdl = t.get("tdl")
                    if sym and tdl == "Tclin":
                        targets[sym.upper()] = 4 # Tclin は原則承認薬扱い(Phase 4と同等)としてマッピング
                    elif sym and tdl == "Tchem":
                        targets[sym.upper()] = max(targets.get(sym.upper(), 0), 1) # Tchem は初期フェーズ扱い

                print(f"[Target] Pharos APIから {len(targets)} 件のターゲットを取得しました。")
        else:
            print(f"[Target] Pharos APIエラー: HTTP {response.status_code}")
    except Exception as e:
        print(f"[Target] Pharos API接続エラー: {e}")

    return targets

def get_chembl_targets(query_str: str, mesh_id: str = None, efo_id: str = None) -> dict:
    """
    ChEMBL APIを使用して、適応症(Indication)→薬剤→メカニズム→ターゲットの流れで抽出
    """
    print(f"[Target] ChEMBL APIを検索中: {query_str} (MeSH: {mesh_id}, EFO: {efo_id})")
    targets = {}
    
    # 1. 適応症(drug_indication)から化合物を引く
    indications_url = "https://www.ebi.ac.uk/chembl/api/data/drug_indication?format=json&limit=1000"
    if mesh_id:
        indications_url += f"&mesh_id={mesh_id}"
    elif efo_id:
        indications_url += f"&efo_id={efo_id}"
    else:
        # String fallback
        indications_url += f"&mesh_heading__icontains={query_str.split()[0]}" 
        
    try:
        resp = requests.get(indications_url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        indications = data.get("drug_indications", [])
        
        drugs = {}
        for ind in indications:
            mol_id = ind.get("molecule_chembl_id")
            phase_val = ind.get("max_phase_for_ind", 0)
            phase = float(phase_val) if phase_val is not None else 0.0
            if mol_id:
                drugs[mol_id] = max(drugs.get(mol_id, 0), int(phase))
                
        if not drugs:
            print("[Target] ChEMBL APIから該当疾患の承認薬が見つかりませんでした。")
            return targets
            
        # 上限を設けてAPIリクエスト過多を防ぐ
        mol_ids = list(drugs.keys())[:200]
        chunk_size = 50
        target_chembl_ids = {}
        
        # 2. 化合物からメカニズムを引く
        for i in range(0, len(mol_ids), chunk_size):
            chunk = mol_ids[i:i+chunk_size]
            mol_filter = ",".join(chunk)
            mech_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id__in={mol_filter}&format=json&limit=1000"
            
            m_resp = requests.get(mech_url, timeout=30)
            if m_resp.status_code == 200:
                m_data = m_resp.json()
                for mech in m_data.get("mechanisms", []):
                    t_id = mech.get("target_chembl_id")
                    m_id = mech.get("molecule_chembl_id")
                    if t_id and m_id:
                        phase = drugs.get(m_id, 0)
                        target_chembl_ids[t_id] = max(target_chembl_ids.get(t_id, 0), phase)
            time.sleep(0.5) # レートリミット対策

        # 3. ターゲット ChEMBL ID から 遺伝子シンボルを引く
        t_ids = list(target_chembl_ids.keys())
        for i in range(0, len(t_ids), chunk_size):
            chunk = t_ids[i:i+chunk_size]
            t_filter = ",".join(chunk)
            t_url = f"https://www.ebi.ac.uk/chembl/api/data/target?target_chembl_id__in={t_filter}&format=json&limit=1000"
            
            t_resp = requests.get(t_url, timeout=30)
            if t_resp.status_code == 200:
                t_data = t_resp.json()
                for t_info in t_data.get("targets", []):
                    t_id = t_info.get("target_chembl_id")
                    phase = target_chembl_ids.get(t_id, 0)
                    
                    components = t_info.get("target_components", [])
                    for comp in components:
                        synonyms = comp.get("target_component_synonyms", [])
                        for syn in synonyms:
                            if syn.get("syn_type") == "GENE_SYMBOL":
                                sym = syn.get("component_synonym")
                                if sym:
                                    # 最高フェーズを更新
                                    targets[sym.upper()] = max(targets.get(sym.upper(), 0), phase)
                                    break
            time.sleep(0.5)

        print(f"[Target] ChEMBL APIから {len(targets)} 件のターゲットを取得しました。")

    except Exception as e:
        print(f"[Target] ChEMBL API接続エラー: {e}")
        
    return targets

def get_combined_targets(disease_name: str, mesh_id: str = None, efo_id: str = None) -> dict:
    """
    PharosとChEMBLの両方からターゲットを取得し統合する
    """
    combined_targets = {}
    
    # 1. Pharos から取得
    pharos_targets = get_pharos_targets(disease_name)
    for sym, phase in pharos_targets.items():
        combined_targets[sym] = phase
        
    # 2. ChEMBL から取得
    chembl_targets = get_chembl_targets(disease_name, mesh_id=mesh_id, efo_id=efo_id)
    for sym, phase in chembl_targets.items():
        # 両方に存在する場合はフェーズが高い方（より開発が進んでいる方）を採用
        combined_targets[sym] = max(combined_targets.get(sym, 0), phase)
        
    print(f"[Target] 【統合結果】合計 {len(combined_targets)} 件の固有ターゲットを抽出しました。")
    return combined_targets

def save_targets(targets_dict: dict, save_dir: str = "data"):
    """
    ターゲット辞書をCSVに保存
    """
    import pandas as pd
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, "known_targets.csv")
    
    if not targets_dict:
        print("[Target] 保存するターゲットがありません。")
        return
        
    df = pd.DataFrame(list(targets_dict.items()), columns=["gene_symbol", "max_phase"])
    df.to_csv(save_path, index=False)
    print(f"[Target] 既知ターゲット情報を保存しました: {save_path}")
