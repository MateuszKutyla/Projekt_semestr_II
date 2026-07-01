#!/usr/bin/env python3
import argparse
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOME = PROJECT_ROOT / "data/assemble_genome/latest_assembly.fasta"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results/gene_prediction/augustus"
PUBLISHED_DIR = PROJECT_ROOT / "data/predicted_genes"


def safe_label(value):
    if not value:
        return "bez_modelu"
    return value.replace("/", "_").replace("\\", "_").replace(" ", "_")


def count_predicted_genes(gff_file):
    count = 0
    with open(gff_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) >= 3 and columns[2] == "gene":
                count += 1
    return count


def wrap_sequence(sequence, width=60):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def extract_proteins_from_augustus_gff(gff_file, output_faa):
    proteins = []
    current = []
    in_protein = False

    with open(gff_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()

            if stripped.startswith("# protein sequence = ["):
                in_protein = True
                first_part = stripped.split("[", 1)[1]
                if first_part.endswith("]"):
                    proteins.append(first_part.rstrip("]"))
                    current = []
                    in_protein = False
                else:
                    current = [first_part]
                continue

            if in_protein:
                part = stripped.lstrip("# ").strip()
                if part.endswith("]"):
                    current.append(part.rstrip("]"))
                    proteins.append("".join(current))
                    current = []
                    in_protein = False
                else:
                    current.append(part)

    with open(output_faa, "w", encoding="utf-8") as output:
        for index, protein in enumerate(proteins, start=1):
            protein = protein.replace(" ", "").replace("\t", "")
            output.write(f">augustus_protein_{index}\n")
            output.write(wrap_sequence(protein) + "\n")

    return len(proteins)


def run_augustus(genome, species, no_species, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    PUBLISHED_DIR.mkdir(parents=True, exist_ok=True)

    model_label = "bez_modelu" if no_species else safe_label(species)
    gff_file = output_dir / f"augustus_{model_label}.gff3"
    log_file = output_dir / f"augustus_{model_label}.log"
    report_file = output_dir / "gene_prediction_report.txt"
    published_gff = PUBLISHED_DIR / "predicted_genes.gff3"
    output_proteins = output_dir / f"augustus_{model_label}.faa"
    published_proteins = PUBLISHED_DIR / "predicted_proteins.faa"

    command = ["augustus", "--gff3=on"]
    if not no_species:
        command.append(f"--species={species}")
    command.append(str(genome))

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

    shutil.copy2(gff_file, published_gff)
    protein_count = extract_proteins_from_augustus_gff(gff_file, output_proteins)
    shutil.copy2(output_proteins, published_proteins)
    gene_count = count_predicted_genes(gff_file)

    with open(report_file, "w", encoding="utf-8") as report:
        report.write("Raport predykcji genow\n")
        report.write("======================\n\n")
        report.write(f"Data analizy: {datetime.now().isoformat(timespec='seconds')}\n")
        report.write(f"Narzedzie: Augustus\n")
        report.write(f"Plik genomu: {genome}\n")
        report.write(f"Model referencyjny: {'bez modelu' if no_species else species}\n")
        report.write(f"Liczba przewidzianych genow: {gene_count}\n")
        report.write(f"Liczba wyekstrahowanych bialek: {protein_count}\n")
        report.write(f"Wynik GFF3: {gff_file}\n")
        report.write(f"Plik GFF3 dla kolejnych etapow: {published_gff}\n")
        report.write(f"Plik bialek dla annotacji funkcjonalnej: {published_proteins}\n")
        report.write(f"Log: {log_file}\n")

    print("Predykcja genow zakonczona.")
    print(f"Liczba przewidzianych genow: {gene_count}")
    print(f"Wynik GFF3: {gff_file}")
    print(f"Plik GFF3 dla kolejnych etapow: {published_gff}")
    print(f"Plik bialek dla annotacji funkcjonalnej: {published_proteins}")
    print(f"Raport: {report_file}")
    print(f"Log: {log_file}")


def main():
    parser = argparse.ArgumentParser(description="Predykcja genow narzedziem Augustus.")
    parser.add_argument("--genome", default=str(DEFAULT_GENOME), help="Plik FASTA ze zlozonym genomem.")
    parser.add_argument("--species", default=None, help="Model gatunkowy Augustusa.")
    parser.add_argument("--no-species", action="store_true", help="Uruchamia Augustusa bez wskazywania modelu gatunkowego.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Katalog wynikowy.")
    args = parser.parse_args()

    genome = Path(args.genome)
    if not genome.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku genomu: {genome}")

    if not args.no_species and not args.species:
        raise ValueError("Podaj --species albo uzyj opcji --no-species.")

    run_augustus(genome, args.species, args.no_species, Path(args.output_dir))


if __name__ == "__main__":
    main()

