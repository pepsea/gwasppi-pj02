import requests
import json

def search_pharos_targets(disease_name):
    # Pharos GraphQL API: https://pharos.nih.gov/api
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

    variables = {
        "name": disease_name
    }

    print(f"Querying Pharos for: {disease_name}")
    try:
        response = requests.post(url, json={'query': query, 'variables': variables})
        response.raise_for_status()
        data = response.json()
        
        diseases = data.get("data", {}).get("diseases", {}).get("diseases", [])
        if not diseases:
            print("No diseases found in Pharos.")
            return

        d = diseases[0]
        print(f"Found disease: {d['name']}")
        
        targets = d.get("targets", [])
        num_targets = len(targets)
        tclin = [t['sym'] for t in targets if t.get('tdl') == 'Tclin']
        
        print(f"Total targets: {num_targets}")
        print(f"Tclin (Approved/Clinical) targets: {len(tclin)}")
        print(tclin[:10])
        
    except Exception as e:
        print(f"Error connecting to Pharos: {e}")

search_pharos_targets("Type 2 diabetes")
