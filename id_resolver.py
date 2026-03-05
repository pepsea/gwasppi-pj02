import requests
import json
import urllib.parse
import time

def auto_resolve_disease(disease_id: str) -> dict:
    """
    EFOまたはMONDO IDから疾患名とMeSH IDを自動解決する
    
    Parameters
    ----------
    disease_id : str
        EFO ID (e.g., EFO_0001360) または MONDO ID (e.g., MONDO_0005148)
        
    Returns
    -------
    dict
        {"trait": "Type 2 diabetes", "mesh_id": "D003924"}
    """
    result = {"trait": None, "mesh_id": None}
    print(f"[Resolver] IDを解決中: {disease_id} ...")
    
    # GWAS Catalog から直接 Trait名 を引く (Primary source of truth for GWAS)
    gwas_url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{disease_id}"
    try:
        resp = requests.get(gwas_url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            result["trait"] = data.get("trait")
            print(f"  [GWAS Catalog] Trait名を取得: {result['trait']}")
            
            # EFO/MONDO -> MeSH is sometimes returned in cross-references or we can use string search
    except Exception as e:
        print(f"  [GWAS Catalog] APIエラー: {e}")
        
    # OLS4 API を使用して MeSH ID などのクロスリファレンスとラベルを取得
    normalized_id = disease_id.replace("_", ":")
    ontology = "efo" if disease_id.startswith("EFO") else "mondo"
    ols_url = f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}/terms?obo_id={normalized_id}"
    
    try:
        resp = requests.get(ols_url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            terms = data.get("_embedded", {}).get("terms", [])
            
            if terms:
                term = terms[0]
                
                # Trait名がGWASから取れなかった場合のフォールバック
                if not result["trait"]:
                    result["trait"] = term.get("label")
                    print(f"  [OLS] ラベルを取得: {result['trait']}")
                
                # MeSH ID の抽出
                mesh_id = None
                
                # パターン1: annotation -> database_cross_reference (MONDO)
                xrefs = term.get("annotation", {}).get("database_cross_reference", [])
                for ref in xrefs:
                    if isinstance(ref, str) and (ref.startswith("MESH:") or ref.startswith("MeSH:")):
                        mesh_id = ref.split(":")[1].strip()
                        break
                        
                # パターン2: xrefs / obo_xref
                if not mesh_id:
                    obo_xrefs = term.get("obo_xref", [])
                    for xref in obo_xrefs:
                        if isinstance(xref, dict) and xref.get("database", "").upper() == "MESH":
                            mesh_id = xref.get("id", "").replace("MESH:", "").replace("MeSH:", "").strip()
                            break
                        elif isinstance(xref, str) and (xref.startswith("MESH:") or xref.startswith("MeSH:")):
                            mesh_id = xref.split(":")[1].strip()
                            break
                            
                result["mesh_id"] = mesh_id
                if mesh_id:
                    print(f"  [OLS] MeSH IDを特定: {mesh_id}")
                else:
                    print("  [OLS] MeSH クロスリファレンスが見つかりません。")
            else:
                print(f"  [OLS] '{normalized_id}' に一致するタームが見つかりません。")
    except Exception as e:
        print(f"  [OLS] APIエラー: {e}")
        
    return result

if __name__ == "__main__":
    print("-" * 40)
    print("Testing obsolete EFO ID:")
    auto_resolve_disease("EFO_0001360")
    print("-" * 40)
    print("Testing active MONDO ID for Type 2 Diabetes:")
    auto_resolve_disease("MONDO_0005148")

def search_disease_id(query: str, max_results: int = 5):
    """
    キーワードから GWAS Catalog 推奨の EFO/MONDO ID を検索して表示するユーティリティ。
    Jupyter Notebook上でユーザーが適切なIDを見つけるために使用する。
    
    Parameters
    ----------
    query : str
        検索キーワード (例: "diabetes", "alzheimer")
    max_results : int
        表示する最大件数
    """
    print(f"'{query}' に関連するGWAS Catalog上の疾患IDを検索中...\n")
    
    # EBI OLS4の検索APIを使用して疾患(EFO / MONDO)を検索する
    encoded_query = urllib.parse.quote(query)
    url = f"https://www.ebi.ac.uk/ols4/api/search?q={encoded_query}&ontology=efo,mondo&exact=false"
    
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            docs = data.get("response", {}).get("docs", [])
            
            if not docs:
                print("該当する疾患は見つかりませんでした。別のキーワード(英語)で試してください。")
                return
                
            print(f"=== 検索結果 (上位 {min(len(docs), max_results)} 件) ===")
            for i, d in enumerate(docs[:max_results]):
                trait_name = d.get("label", "Unknown")
                short_form = d.get("short_form", "Unknown ID")
                # 疾患名とIDを表示 (GWAS Catalog 検索用のアノテーションとして)
                print(f"{i+1}. ID: {short_form} | 疾患名: {trait_name}")
                
            print("\n→ 使用するIDを下のセルの TARGET_DISEASE_ID にコピーして設定してください。")
        else:
            print(f"検索APIエラー: ステータスコード {resp.status_code}")
    except Exception as e:
        print(f"検索中にエラーが発生しました: {e}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "search":
        search_disease_id(sys.argv[2])
