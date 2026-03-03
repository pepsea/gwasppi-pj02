import requests

def get_targets(disease_id='EFO_0001360'):
    query = """
    query getTargets($efoId: String!) {
      disease(efoId: $efoId) {
        knownDrugs {
          rows {
            drug { name }
            target { approvedSymbol }
            phase
            status
          }
        }
      }
    }
    """
    url = 'https://api.opentargets.io/api/v4/graphql'
    resp = requests.post(url, json={'query': query, 'variables': {'efoId': disease_id}})
    if resp.status_code != 200:
        print("API Error:", resp.status_code)
        return
        
    data = resp.json()
    rows = data.get('data', {}).get('disease', {}).get('knownDrugs', {}).get('rows', [])
    targets = {}
    for r in rows:
        target = r.get('target', {}).get('approvedSymbol')
        phase = r.get('phase', 0)
        if target:
            targets[target] = max(targets.get(target, 0), phase)
                
    print(f'Found {len(targets)} targets')
    print(list(targets.items())[:10])
    
get_targets()
