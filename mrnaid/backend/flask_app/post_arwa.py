from requests import post, get

# add ../common/arw_mrna/src to sys.path
import sys
sys.path.append("../common/arw_mrna/src")

import protein
from awalk import adaptive_random_walk, WalkConfig
def get_arwa_task(id):
    url = f"http://localhost:5000/api/v1/status/{id}"
    return get(url)

def post_arwa_async():
    url = "http://localhost:5000/api/v1/arwa"
    headers = {
        'Content-Type': 'application/json'
    }
    data = {
        "aa_seq": "MVSKGEELFTGVVPILVELDGDVNGH",
        "steps": 10,
        "verbose": True,
        "stability": "efe",
        "cai_exp_scale": 1.0,
        "cai_threshold": 0.9,
        "freq_table_path": "homosapiens.txt",
    }
    return post(url, headers=headers, json=data).json()

def post_arwa_sync():
    url = "http://localhost:5000/api/v1/arwa_sync"
    headers = {
        'Content-Type': 'application/json'
    }
    data = {
        "aa_seq": "MVSKGEELFTGVVPILVELDGDVNGH",
        "steps": 100,
        "verbose": True,
        "stability": "efe",
        "cai_exp_scale": 1.0,
        "cai_threshold": 0.9,
        "freq_table_path": "homosapiens.txt",
    }
    response = post(url, headers=headers, json=data)
    print(response.text)

def get_file(path):
    url = "http://localhost:5000/api/v1/file"
    headers = {
        'Content-Type': 'application/json'
    }
    data = {
        "path": path
    }
    response = post(url, headers=headers, json=data)
    print(response.text)

if __name__ == "__main__":
    post_arwa_async()