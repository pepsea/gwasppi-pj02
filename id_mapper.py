"""
ID Mapper Module
================
遺伝子 ID 変換ユーティリティ
Gene Symbol ↔ UniProt ↔ Entrez Gene ID
"""

import time
import requests
import pandas as pd

import config


def gene_symbols_to_uniprot(gene_symbols: list, species: int = 9606) -> dict:
    """
    Gene Symbol → UniProt ID 変換 (UniProt ID Mapping API)

    Parameters
    ----------
    gene_symbols : list
        遺伝子シンボルのリスト
    species : int
        NCBI Taxonomy ID (default: 9606 = Homo sapiens)

    Returns
    -------
    dict : {gene_symbol: uniprot_id}
    """
    if not gene_symbols:
        return {}

    url = f"{config.UNIPROT_ID_MAPPING_URL}/run"
    params = {
        "from": "Gene_Name",
        "to": "UniProtKB",
        "ids": ",".join(gene_symbols),
        "taxId": species,
    }

    try:
        resp = requests.post(url, data=params, timeout=30)
        resp.raise_for_status()
        job_id = resp.json().get("jobId")

        if not job_id:
            print("[ID Mapper] ジョブID取得失敗")
            return {}

        # ポーリングで結果待ち
        result_url = f"{config.UNIPROT_ID_MAPPING_URL}/results/{job_id}"
        for _ in range(30):
            time.sleep(2)
            status_resp = requests.get(
                f"{config.UNIPROT_ID_MAPPING_URL}/status/{job_id}", timeout=10
            )
            status_data = status_resp.json()
            if "results" in status_data or "jobStatus" not in status_data:
                break
            if status_data.get("jobStatus") == "FINISHED":
                break

        result_resp = requests.get(
            result_url, params={"format": "json", "size": 500}, timeout=30
        )
        result_resp.raise_for_status()
        results = result_resp.json().get("results", [])

        mapping = {}
        for r in results:
            from_id = r.get("from", "")
            to_entry = r.get("to", {})
            if isinstance(to_entry, dict):
                uniprot_id = to_entry.get("primaryAccession", "")
            else:
                uniprot_id = str(to_entry)
            if from_id and uniprot_id:
                if from_id not in mapping:
                    mapping[from_id] = uniprot_id

        print(f"[ID Mapper] {len(mapping)}/{len(gene_symbols)} Gene Symbol → UniProt ID")
        return mapping

    except Exception as e:
        print(f"[ID Mapper] UniProt ID Mapping エラー: {e}")
        return {}


def uniprot_to_gene_symbols(uniprot_ids: list) -> dict:
    """
    UniProt ID → Gene Symbol 変換

    Parameters
    ----------
    uniprot_ids : list
        UniProt ID のリスト

    Returns
    -------
    dict : {uniprot_id: gene_symbol}
    """
    if not uniprot_ids:
        return {}

    url = f"{config.UNIPROT_ID_MAPPING_URL}/run"
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "Gene_Name",
        "ids": ",".join(uniprot_ids),
    }

    try:
        resp = requests.post(url, data=params, timeout=30)
        resp.raise_for_status()
        job_id = resp.json().get("jobId")

        if not job_id:
            return {}

        for _ in range(30):
            time.sleep(2)
            status_resp = requests.get(
                f"{config.UNIPROT_ID_MAPPING_URL}/status/{job_id}", timeout=10
            )
            status_data = status_resp.json()
            if "results" in status_data or "jobStatus" not in status_data:
                break
            if status_data.get("jobStatus") == "FINISHED":
                break

        result_resp = requests.get(
            f"{config.UNIPROT_ID_MAPPING_URL}/results/{job_id}",
            params={"format": "json", "size": 500},
            timeout=30,
        )
        result_resp.raise_for_status()
        results = result_resp.json().get("results", [])

        mapping = {}
        for r in results:
            from_id = r.get("from", "")
            to_id = r.get("to", "")
            if from_id and to_id:
                if from_id not in mapping:
                    mapping[from_id] = to_id

        print(f"[ID Mapper] {len(mapping)}/{len(uniprot_ids)} UniProt ID → Gene Symbol")
        return mapping

    except Exception as e:
        print(f"[ID Mapper] UniProt → Gene Symbol エラー: {e}")
        return {}


def normalize_gene_list(genes: list) -> list:
    """遺伝子リストを正規化 (大文字、重複除去、ソート)"""
    return sorted(set(g.upper().strip() for g in genes if g and g.strip()))
