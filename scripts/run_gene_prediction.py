#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOME = PROJECT_ROOT / "data/assemble_genome/latest_assembly.fasta"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/gene_prediction/augustus"


def run_augustus(genome, species, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    gff_file = output_dir / f"augustus_{species}.gff"
    log_file = output_dir / f"augustus_{species}.log"

    command = [
        "augustus",
        f"--species={species}",
        "--gff3=on",
        str(genome)
    ]

    with open(gff_file, "w", encoding="utf-8") as gff, open(log_file, "w", encoding="utf-8") as log:
        log.write("Komenda:\n")
        log.write(" ".join(command) + "\n\n")
        process = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            stdout=gff,
            stderr=log,
            text=True
        )

    if process.returncode != 0:
        raise RuntimeError(f"Augustus zakonczyl sie bledem. Szczegoly: {log_file}")

    print(f"Predykcja genow zakonczona.")
    print(f"Wynik GFF3: {gff_file}")
    print(f"Log: {log_file}")


def main():
    parser = argparse.ArgumentParser(description="Predykcja genow narzedziem Augustus.")
    parser.add_argument("--genome", default=str(DEFAULT_GENOME), help="Plik FASTA ze zlozonym genomem.")
    parser.add_argument("--species", required=True, help="Model gatunkowy Augustusa.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Katalog wynikowy.")
    args = parser.parse_args()

    genome = Path(args.genome)
    if not genome.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku genomu: {genome}")

    run_augustus(genome, args.species, Path(args.output_dir))


if __name__ == "__main__":
    main()
