#!/usr/bin/env python3
import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from pipeline_config import get_value, load_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run_step(name, command, log_handle):
    log_handle.write("\n" + "=" * 80 + "\n")
    log_handle.write(name + "\n")
    log_handle.write(" ".join(command) + "\n")
    log_handle.write("=" * 80 + "\n")
    log_handle.flush()

    process = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True
    )

    if process.returncode != 0:
        raise RuntimeError(f"Etap zakonczyl sie bledem: {name}")


def main():
    parser = argparse.ArgumentParser(description="Pelny pipeline analizy genomu grzyba.")
    parser.add_argument("--assembly-mode", choices=["illumina", "ont", "hybrid"], default=None)
    parser.add_argument("--species", default=None)
    parser.add_argument("--no-species", action="store_true")
    parser.add_argument("--functional-diamond-db", default=None)
    parser.add_argument("--hydrolase-hmm-db", default=None)
    parser.add_argument("--hydrolase-diamond-db", default=None)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--config", default=None, help="Plik konfiguracyjny YAML.")
    args = parser.parse_args()

    config = load_config(args.config)

    args.assembly_mode = args.assembly_mode or get_value(config, "assembly.mode", "hybrid")
    args.species = args.species or get_value(config, "gene_prediction.species", None)
    args.functional_diamond_db = args.functional_diamond_db or get_value(config, "functional_annotation.databases.diamond_db", None)
    args.hydrolase_hmm_db = args.hydrolase_hmm_db or get_value(config, "hydrolase_prediction.databases.hmm_db", None)
    args.hydrolase_diamond_db = args.hydrolase_diamond_db or get_value(config, "hydrolase_prediction.databases.diamond_db", None)
    args.threads = args.threads or int(get_value(config, "project.threads", 8))

    if not args.no_species and not args.species:
        raise ValueError("Podaj --species albo ustaw gene_prediction.species w config/config.yaml.")

    if not args.functional_diamond_db:
        raise ValueError("Podaj --functional-diamond-db albo ustaw functional_annotation.databases.diamond_db w config/config.yaml.")
    if not args.hydrolase_hmm_db:
        raise ValueError("Podaj --hydrolase-hmm-db albo ustaw hydrolase_prediction.databases.hmm_db w config/config.yaml.")
    if not args.hydrolase_diamond_db:
        raise ValueError("Podaj --hydrolase-diamond-db albo ustaw hydrolase_prediction.databases.diamond_db w config/config.yaml.")

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "full_pipeline.log"

    with open(log_file, "w", encoding="utf-8") as log:
        log.write("Pelny pipeline analizy genomu grzyba\n")
        log.write(f"Start: {datetime.now().isoformat(timespec='seconds')}\n")

        run_step(
            "1. Assemblacja genomu",
            [
                sys.executable,
                "scripts/run_denovo_assembly.py",
                "--mode",
                args.assembly_mode,
                "--threads",
                str(args.threads)
            ],
            log
        )

        gene_command = [
            sys.executable,
            "scripts/run_gene_prediction.py",
            "--genome",
            "data/assemble_genome/latest_assembly.fasta"
        ]

        if args.no_species:
            gene_command.append("--no-species")
        else:
            gene_command.extend(["--species", args.species])

        run_step("2. Predykcja genow", gene_command, log)

        run_step(
            "3. Annotacja funkcjonalna",
            [
                sys.executable,
                "scripts/run_functional_annotation.py",
                "--tool",
                "all",
                "--proteins",
                "data/predicted_genes/predicted_proteins.faa",
                "--diamond-db",
                args.functional_diamond_db,
                "--threads",
                str(args.threads)
            ],
            log
        )

        run_step(
            "4. Predykcja hydrolaz",
            [
                sys.executable,
                "scripts/run_hydrolase_prediction.py",
                "--tool",
                "all",
                "--proteins",
                "data/predicted_genes/predicted_proteins.faa",
                "--hmm-db",
                args.hydrolase_hmm_db,
                "--diamond-db",
                args.hydrolase_diamond_db,
                "--threads",
                str(args.threads)
            ],
            log
        )

        run_step(
            "5. Raport koncowy",
            [
                sys.executable,
                "scripts/build_final_report.py"
            ],
            log
        )

        log.write(f"\nKoniec: {datetime.now().isoformat(timespec='seconds')}\n")
        log.write("Pelny pipeline zakonczony powodzeniem.\n")

    print(f"Pelny pipeline zakonczony. Log: {log_file}")


if __name__ == "__main__":
    main()




