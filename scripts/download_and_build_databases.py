#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

UNIPROT_SPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
PFAM_A_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def download_file(url, output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    print(f"Pobieranie: {url}")
    print(f"Do pliku: {output_file}")
    urllib.request.urlretrieve(url, output_file)


def gunzip_file(gz_file, output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    print(f"Rozpakowywanie: {gz_file}")
    with gzip.open(gz_file, "rb") as source:
        with open(output_file, "wb") as target:
            shutil.copyfileobj(source, target)


def run_command(command):
    print("Komenda:")
    print(" ".join(str(item) for item in command))
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def build_diamond_database(output_db):
    output_db = Path(output_db)
    fasta_gz = PROJECT_ROOT / "data/databases/downloads/uniprot_sprot.fasta.gz"
    fasta = PROJECT_ROOT / "data/databases/downloads/uniprot_sprot.fasta"

    download_file(UNIPROT_SPROT_URL, fasta_gz)
    gunzip_file(fasta_gz, fasta)

    run_command([
        "diamond",
        "makedb",
        "--in", str(fasta),
        "--db", str(output_db)
    ])

    print(f"Gotowa baza DIAMOND: {output_db}")


def build_hmm_database(output_hmm):
    output_hmm = Path(output_hmm)
    hmm_gz = PROJECT_ROOT / "data/databases/downloads/Pfam-A.hmm.gz"

    download_file(PFAM_A_HMM_URL, hmm_gz)
    gunzip_file(hmm_gz, output_hmm)

    run_command([
        "hmmpress",
        str(output_hmm)
    ])

    print(f"Gotowa baza HMM: {output_hmm}")


def main():
    parser = argparse.ArgumentParser(description="Pobieranie i budowanie baz danych dla pipeline'u.")
    parser.add_argument("--kind", choices=["diamond", "hmm"], required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.kind == "diamond":
        build_diamond_database(args.output)
    elif args.kind == "hmm":
        build_hmm_database(args.output)


if __name__ == "__main__":
    main()
