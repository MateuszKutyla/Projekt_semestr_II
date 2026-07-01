#!/usr/bin/env python3
import argparse
import subprocess
import urllib.request
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROTEINS = PROJECT_ROOT / "data/predicted_genes/predicted_proteins.faa"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/functional_annotation"
PUBLISHED_DIR = PROJECT_ROOT / "data/functional_annotation"
DIAMOND_DB_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
DEFAULT_DIAMOND_DB_PREFIX = PUBLISHED_DIR / "diamond_database" / "uniprot_sprot"


def run_command(name, command, log_file):
    with open(log_file, "a", encoding="utf-8") as log:
        log.write("\n" + "=" * 80 + "\n")
        log.write(name + "\n")
        log.write(" ".join(str(item) for item in command) + "\n")
        log.write("=" * 80 + "\n")

        process = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )

    if process.returncode != 0:
        raise RuntimeError(f"Etap zakonczyl sie bledem: {name}. Szczegoly: {log_file}")


def download_file(url, output_file, log_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "a", encoding="utf-8") as log:
        log.write("\n" + "=" * 80 + "\n")
        log.write("Pobieranie aktualnej bazy Swiss-Prot z UniProt\n")
        log.write(f"URL: {url}\n")
        log.write(f"Plik: {output_file}\n")
        log.write("=" * 80 + "\n")

    urllib.request.urlretrieve(url, output_file)


def prepare_diamond_database(output_dir):
    tool_dir = output_dir / "diamond"
    tool_dir.mkdir(parents=True, exist_ok=True)
    log_file = tool_dir / "diamond_database.log"

    fasta_file = DEFAULT_DIAMOND_DB_PREFIX.with_suffix(".fasta.gz")
    db_prefix = DEFAULT_DIAMOND_DB_PREFIX
    db_file = DEFAULT_DIAMOND_DB_PREFIX.with_suffix(".dmnd")

    if not fasta_file.exists():
        download_file(DIAMOND_DB_URL, fasta_file, log_file)

    if not db_file.exists():
        run_command(
            "Budowanie bazy DIAMOND z aktualnego Swiss-Prot",
            ["diamond", "makedb", "--in", fasta_file, "--db", db_prefix],
            log_file
        )

    return db_file, log_file


def run_diamond(proteins, output_dir, database, threads):
    tool_dir = output_dir / "diamond"
    tool_dir.mkdir(parents=True, exist_ok=True)

    output_file = tool_dir / "diamond_matches.tsv"
    log_file = tool_dir / "diamond.log"

    command = [
        "diamond",
        "blastp",
        "--query", proteins,
        "--db", database,
        "--out", output_file,
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
        "--threads", str(threads),
        "--max-target-seqs", "1",
        "--evalue", "1e-5"
    ]

    run_command("Annotacja funkcjonalna - DIAMOND blastp", command, log_file)
    return output_file, log_file


def run_eggnog(proteins, output_dir, threads):
    tool_dir = output_dir / "eggnog_mapper"
    tool_dir.mkdir(parents=True, exist_ok=True)

    log_file = tool_dir / "eggnog_mapper.log"

    command = [
        "emapper.py",
        "-i", proteins,
        "-o", "eggnog_annotation",
        "--output_dir", tool_dir,
        "--cpu", str(threads)
    ]

    run_command("Annotacja funkcjonalna - eggNOG-mapper", command, log_file)
    return tool_dir / "eggnog_annotation.emapper.annotations", log_file


def run_interproscan(proteins, output_dir, threads):
    tool_dir = output_dir / "interproscan"
    tool_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = tool_dir / "interproscan_annotation"
    log_file = tool_dir / "interproscan.log"

    command = [
        "interproscan.sh",
        "-i", proteins,
        "-f", "TSV,GFF3",
        "-o", str(output_prefix) + ".tsv",
        "-cpu", str(threads),
        "-goterms",
        "-pa"
    ]

    run_command("Annotacja funkcjonalna - InterProScan", command, log_file)
    return output_prefix.with_suffix(".tsv"), log_file


def write_report(tool, proteins, result_file, log_file, output_dir):
    report_file = output_dir / "functional_annotation_report.txt"
    with open(report_file, "w", encoding="utf-8") as report:
        report.write("Raport annotacji funkcjonalnej\n")
        report.write("==============================\n\n")
        report.write(f"Data analizy: {datetime.now().isoformat(timespec='seconds')}\n")
        report.write(f"Narzędzie: {tool}\n")
        report.write(f"Plik bialek: {proteins}\n")
        report.write(f"Plik wynikowy: {result_file}\n")
        report.write(f"Log: {log_file}\n")
    return report_file


def main():
    parser = argparse.ArgumentParser(description="Annotacja funkcjonalna przewidywanych bialek.")
    parser.add_argument("--tool", choices=["diamond", "eggnog", "interproscan", "all"], required=True)
    parser.add_argument("--proteins", default=str(DEFAULT_PROTEINS), help="Plik FASTA z bialkami.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Katalog wynikowy.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--diamond-db", default=None, help="Baza DIAMOND, wymagana dla --tool diamond.")
    parser.add_argument("--download-diamond-db", action="store_true", help="Pobiera aktualny Swiss-Prot z UniProt i buduje baze DIAMOND.")
    args = parser.parse_args()

    proteins = Path(args.proteins)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    PUBLISHED_DIR.mkdir(parents=True, exist_ok=True)

    if not proteins.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku bialek: {proteins}")

    results = []

    if args.tool in ["diamond", "all"]:
        if args.download_diamond_db:
            diamond_db, db_log_file = prepare_diamond_database(output_dir)
            results.append(("DIAMOND database", diamond_db, db_log_file))
        elif args.diamond_db:
            diamond_db = args.diamond_db
        else:
            raise ValueError("Dla DIAMOND trzeba podac --diamond-db albo uzyc --download-diamond-db.")

        result_file, log_file = run_diamond(proteins, output_dir, diamond_db, args.threads)
        results.append(("DIAMOND", result_file, log_file))

    if args.tool in ["eggnog", "all"]:
        result_file, log_file = run_eggnog(proteins, output_dir, args.threads)
        results.append(("eggNOG-mapper", result_file, log_file))

    if args.tool in ["interproscan", "all"]:
        result_file, log_file = run_interproscan(proteins, output_dir, args.threads)
        results.append(("InterProScan", result_file, log_file))

    report_file = output_dir / "functional_annotation_report.txt"
    with open(report_file, "w", encoding="utf-8") as report:
        report.write("Raport annotacji funkcjonalnej\n")
        report.write("==============================\n\n")
        report.write(f"Data analizy: {datetime.now().isoformat(timespec='seconds')}\n")
        report.write(f"Plik bialek: {proteins}\n")
        report.write(f"Tryb: {args.tool}\n\n")
        for tool_name, result_file, log_file in results:
            report.write(f"Narzędzie: {tool_name}\n")
            report.write(f"Plik wynikowy: {result_file}\n")
            report.write(f"Log: {log_file}\n\n")

    print("Annotacja funkcjonalna zakonczona.")
    print(f"Tryb: {args.tool}")
    for tool_name, result_file, log_file in results:
        print(f"{tool_name}: {result_file}")
    print(f"Raport: {report_file}")


if __name__ == "__main__":
    main()


