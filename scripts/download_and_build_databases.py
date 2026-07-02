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

HYDROLASE_KEYWORDS = [
    "hydrolase",
    "glycosidase",
    "glucosidase",
    "cellulase",
    "xylanase",
    "amylase",
    "lipase",
    "esterase",
    "protease",
    "peptidase",
    "chitinase",
    "mannanase",
    "pectinase",
    "polygalacturonase",
    "cutinase",
    "laccase",
    "cazyme",
    "carbohydrate-active",
]


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


def ensure_swissprot_fasta():
    fasta_gz = PROJECT_ROOT / "data/databases/downloads/uniprot_sprot.fasta.gz"
    fasta = PROJECT_ROOT / "data/databases/downloads/uniprot_sprot.fasta"

    if not fasta.exists():
        download_file(UNIPROT_SPROT_URL, fasta_gz)
        gunzip_file(fasta_gz, fasta)

    return fasta


def write_record_if_hydrolase(header, sequence_lines, output_handle):
    lower_header = header.lower()
    if any(keyword in lower_header for keyword in HYDROLASE_KEYWORDS):
        output_handle.write(header)
        output_handle.writelines(sequence_lines)
        return True
    return False


def filter_hydrolase_fasta(input_fasta, output_fasta):
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    total = 0
    header = None
    sequence_lines = []

    with open(input_fasta, "r", encoding="utf-8", errors="replace") as source, open(output_fasta, "w", encoding="utf-8") as target:
        for line in source:
            if line.startswith(">"):
                if header is not None:
                    total += 1
                    if write_record_if_hydrolase(header, sequence_lines, target):
                        kept += 1
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header is not None:
            total += 1
            if write_record_if_hydrolase(header, sequence_lines, target):
                kept += 1

    if kept == 0:
        raise RuntimeError("Filtrowanie nie znalazło żadnych rekordów hydrolaz w Swiss-Prot.")

    print(f"Przefiltrowano rekordy hydrolaz: {kept} z {total}")
    print(f"FASTA hydrolaz: {output_fasta}")

    return output_fasta


def build_diamond_database(output_db):
    output_db = Path(output_db)
    fasta = ensure_swissprot_fasta()

    run_command([
        "diamond",
        "makedb",
        "--in", str(fasta),
        "--db", str(output_db)
    ])

    print(f"Gotowa baza DIAMOND: {output_db}")


def build_hydrolase_diamond_database(output_db):
    output_db = Path(output_db)
    swissprot_fasta = ensure_swissprot_fasta()
    hydrolase_fasta = PROJECT_ROOT / "data/databases/downloads/uniprot_sprot_hydrolases.fasta"

    filter_hydrolase_fasta(swissprot_fasta, hydrolase_fasta)

    run_command([
        "diamond",
        "makedb",
        "--in", str(hydrolase_fasta),
        "--db", str(output_db)
    ])

    print(f"Gotowa baza DIAMOND hydrolaz: {output_db}")


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
    parser.add_argument("--kind", choices=["diamond", "diamond_hydrolase", "hmm"], required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.kind == "diamond":
        build_diamond_database(args.output)
    elif args.kind == "diamond_hydrolase":
        build_hydrolase_diamond_database(args.output)
    elif args.kind == "hmm":
        build_hmm_database(args.output)


if __name__ == "__main__":
    main()
