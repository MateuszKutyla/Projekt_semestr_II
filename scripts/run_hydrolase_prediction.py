#!/usr/bin/env python3
import argparse
import subprocess
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROTEINS = PROJECT_ROOT / "data/predicted_genes/predicted_proteins.faa"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/hydrolase_prediction"
PUBLISHED_DIR = PROJECT_ROOT / "data/hydrolases"


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


def run_dbcan(proteins, output_dir, threads):
    tool_dir = output_dir / "dbcan"
    tool_dir.mkdir(parents=True, exist_ok=True)
    log_file = tool_dir / "dbcan.log"

    command = [
        "run_dbcan",
        str(proteins),
        "protein",
        "--out_dir", str(tool_dir),
        "--db_dir", str(PROJECT_ROOT / "data/databases/dbcan"),
        "--dia_cpu", str(threads),
        "--hmm_cpu", str(threads),
        "--tf_cpu", str(threads)
    ]

    run_command("Predykcja hydrolaz / CAZymes - dbCAN", command, log_file)
    return tool_dir, log_file


def run_hmmer(proteins, output_dir, hmm_database, threads):
    tool_dir = output_dir / "hmmer"
    tool_dir.mkdir(parents=True, exist_ok=True)

    table_file = tool_dir / "hydrolase_hmmer.tbl"
    domain_file = tool_dir / "hydrolase_hmmer.domtbl"
    log_file = tool_dir / "hmmer.log"

    command = [
        "hmmscan",
        "--cpu", str(threads),
        "--tblout", str(table_file),
        "--domtblout", str(domain_file),
        str(hmm_database),
        str(proteins)
    ]

    run_command("Predykcja hydrolaz - HMMER", command, log_file)
    return domain_file, log_file


def run_diamond_hydrolases(proteins, output_dir, diamond_database, threads):
    tool_dir = output_dir / "diamond_hydrolases"
    tool_dir.mkdir(parents=True, exist_ok=True)

    output_file = tool_dir / "hydrolase_diamond_matches.tsv"
    log_file = tool_dir / "diamond_hydrolases.log"

    command = [
        "diamond",
        "blastp",
        "--query", str(proteins),
        "--db", str(diamond_database),
        "--out", str(output_file),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
        "--threads", str(threads),
        "--max-target-seqs", "5",
        "--evalue", "1e-5"
    ]

    run_command("Predykcja hydrolaz - DIAMOND", command, log_file)
    return output_file, log_file


def run_signalp(proteins, output_dir):
    tool_dir = output_dir / "signalp"
    tool_dir.mkdir(parents=True, exist_ok=True)
    log_file = tool_dir / "signalp.log"

    command = [
        "signalp6",
        "--fastafile", str(proteins),
        "--organism", "eukarya",
        "--output_dir", str(tool_dir),
        "--format", "txt"
    ]

    run_command("Predykcja peptydow sygnalowych - SignalP", command, log_file)
    return tool_dir, log_file


def run_deeptmhmm(proteins, output_dir):
    tool_dir = output_dir / "deeptmhmm"
    tool_dir.mkdir(parents=True, exist_ok=True)
    log_file = tool_dir / "deeptmhmm.log"

    command = [
        "biolib",
        "run",
        "DTU/DeepTMHMM",
        "--fasta", str(proteins)
    ]

    run_command("Predykcja domen transblonowych - DeepTMHMM", command, log_file)
    return tool_dir, log_file


def add_prediction(predictions, protein_id, source, hit_id, description, evalue="", score=""):
    predictions.append({
        "protein_id": protein_id,
        "source": source,
        "hit_id": hit_id,
        "description": description,
        "evalue": evalue,
        "score": score
    })


def parse_diamond_results(result_file, predictions):
    if not Path(result_file).exists():
        return

    with open(result_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 7:
                continue
            add_prediction(
                predictions,
                protein_id=columns[0],
                source="DIAMOND",
                hit_id=columns[1],
                description=columns[6],
                evalue=columns[4],
                score=columns[5]
            )


def parse_hmmer_results(result_file, predictions):
    if not Path(result_file).exists():
        return

    with open(result_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) < 7:
                continue
            add_prediction(
                predictions,
                protein_id=columns[3],
                source="HMMER",
                hit_id=columns[0],
                description=" ".join(columns[22:]) if len(columns) > 22 else columns[0],
                evalue=columns[6],
                score=columns[7] if len(columns) > 7 else ""
            )


def parse_dbcan_results(result_dir, predictions):
    result_dir = Path(result_dir)
    candidates = [
        result_dir / "overview.txt",
        result_dir / "hmmer.out",
        result_dir / "diamond.out",
        result_dir / "Hotpep.out"
    ]

    for candidate in candidates:
        if not candidate.exists():
            continue

        with open(candidate, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.lower().startswith("gene") or line.startswith("#"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 2:
                    columns = line.split()
                if len(columns) < 2:
                    continue

                add_prediction(
                    predictions,
                    protein_id=columns[0],
                    source="dbCAN",
                    hit_id=columns[1],
                    description=";".join(columns[2:]) if len(columns) > 2 else "CAZyme/hydrolase candidate"
                )


def write_combined_hydrolase_table(predictions, output_file):
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    with open(output_file, "w", encoding="utf-8") as output:
        output.write("protein_id\tsource\thit_id\tdescription\tevalue\tscore\n")
        for item in predictions:
            key = (
                item["protein_id"],
                item["source"],
                item["hit_id"],
                item["description"],
                item["evalue"],
                item["score"]
            )
            if key in seen:
                continue
            seen.add(key)
            output.write(
                f"{item['protein_id']}\t{item['source']}\t{item['hit_id']}\t"
                f"{item['description']}\t{item['evalue']}\t{item['score']}\n"
            )

    return output_file


def write_report(results, proteins, output_dir):
    report_file = output_dir / "hydrolase_prediction_report.txt"

    with open(report_file, "w", encoding="utf-8") as report:
        report.write("Raport predykcji hydrolaz\n")
        report.write("=========================\n\n")
        report.write(f"Data analizy: {datetime.now().isoformat(timespec='seconds')}\n")
        report.write(f"Plik bialek: {proteins}\n\n")

        for tool_name, result_path, log_file in results:
            report.write(f"Narzedzie: {tool_name}\n")
            report.write(f"Wynik: {result_path}\n")
            report.write(f"Log: {log_file}\n\n")

    return report_file


def main():
    parser = argparse.ArgumentParser(description="Predykcja hydrolaz na podstawie przewidywanych bialek.")
    parser.add_argument("--proteins", default=str(DEFAULT_PROTEINS), help="Plik FASTA z bialkami.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Katalog wynikowy.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--tool", choices=["dbcan", "hmmer", "diamond", "signalp", "deeptmhmm", "all"], default="all")
    parser.add_argument("--hmm-db", default=None, help="Baza HMM dla HMMER, np. Pfam-A.hmm albo baza hydrolaz.")
    parser.add_argument("--diamond-db", default=None, help="Baza DIAMOND z sekwencjami hydrolaz.")
    args = parser.parse_args()

    proteins = Path(args.proteins)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    PUBLISHED_DIR.mkdir(parents=True, exist_ok=True)

    if not proteins.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku bialek: {proteins}")

    results = []
    predictions = []

    if args.tool in ["dbcan", "all"]:
        result_path, log_file = run_dbcan(proteins, output_dir, args.threads)
        results.append(("dbCAN", result_path, log_file))
        parse_dbcan_results(result_path, predictions)

    if args.tool in ["hmmer", "all"]:
        if not args.hmm_db:
            raise ValueError("Dla HMMER trzeba podac --hmm-db.")
        result_path, log_file = run_hmmer(proteins, output_dir, Path(args.hmm_db), args.threads)
        results.append(("HMMER", result_path, log_file))
        parse_hmmer_results(result_path, predictions)

    if args.tool in ["diamond", "all"]:
        if not args.diamond_db:
            raise ValueError("Dla DIAMOND trzeba podac --diamond-db.")
        result_path, log_file = run_diamond_hydrolases(proteins, output_dir, Path(args.diamond_db), args.threads)
        results.append(("DIAMOND hydrolases", result_path, log_file))
        parse_diamond_results(result_path, predictions)

    if args.tool in ["signalp", "all"]:
        result_path, log_file = run_signalp(proteins, output_dir)
        results.append(("SignalP", result_path, log_file))

    if args.tool in ["deeptmhmm", "all"]:
        result_path, log_file = run_deeptmhmm(proteins, output_dir)
        results.append(("DeepTMHMM", result_path, log_file))

    combined_file = write_combined_hydrolase_table(predictions, PUBLISHED_DIR / "predicted_hydrolases.tsv")
    report_file = write_report(results, proteins, output_dir)

    print("Predykcja hydrolaz zakonczona.")
    print(f"Raport: {report_file}")
    print(f"Zbiorczy plik predykcji hydrolaz: {combined_file}")


if __name__ == "__main__":
    main()

