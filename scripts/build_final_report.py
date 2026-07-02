#!/usr/bin/env python3
import argparse
import json
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def read_fasta_stats(path):
    path = Path(path)
    if not path.exists():
        return None

    lengths = []
    current = 0

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)

    if current > 0:
        lengths.append(current)

    if not lengths:
        return {
            "file": str(path),
            "contigs": 0,
            "total_length": 0,
            "n50": 0,
            "longest_contig": 0,
        }

    sorted_lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    running = 0
    n50 = 0

    for length in sorted_lengths:
        running += length
        if running >= half:
            n50 = length
            break

    return {
        "file": str(path),
        "contigs": len(lengths),
        "total_length": sum(lengths),
        "n50": n50,
        "longest_contig": max(lengths),
    }


def count_gff_genes(path):
    path = Path(path)
    if not path.exists():
        return None

    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) >= 3 and columns[2] == "gene":
                count += 1
    return count


def count_fasta_records(path):
    path = Path(path)
    if not path.exists():
        return None

    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def count_table_rows(path):
    path = Path(path)
    if not path.exists():
        return None

    rows = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                rows += 1

    if rows > 0:
        rows -= 1

    return max(rows, 0)


def find_first_existing(patterns):
    for pattern in patterns:
        found = sorted(PROJECT_ROOT.glob(pattern))
        if found:
            return found[0]
    return None


def read_short_text(path, max_lines=12):
    if not path or not Path(path).exists():
        return None

    lines = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            if index >= max_lines:
                break
            lines.append(line.rstrip())
    return lines


def build_summary():
    assembly_file = PROJECT_ROOT / "data" / "assemble_genome" / "latest_assembly.fasta"
    genes_gff = PROJECT_ROOT / "data" / "predicted_genes" / "predicted_genes.gff3"
    proteins_faa = PROJECT_ROOT / "data" / "predicted_genes" / "predicted_proteins.faa"
    hydrolases_tsv = PROJECT_ROOT / "data" / "hydrolases" / "predicted_hydrolases.tsv"

    functional_file = find_first_existing([
        "results/functional_annotation/*combined*.tsv",
        "results/functional_annotation/*.tsv",
        "data/functional_annotation/*.tsv",
    ])

    quast_file = find_first_existing([
        "results/assembly_qc/*/quast/report.tsv",
        "results/assembly_qc/*/quast/transposed_report.tsv",
    ])

    busco_file = find_first_existing([
        "results/assembly_qc/*/busco_*/short_summary*.txt",
        "results/assembly_qc/*/short_summary*.txt",
    ])

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "assembly": read_fasta_stats(assembly_file),
        "gene_prediction": {
            "gff": str(genes_gff),
            "predicted_genes": count_gff_genes(genes_gff),
            "proteins": str(proteins_faa),
            "predicted_proteins": count_fasta_records(proteins_faa),
        },
        "functional_annotation": {
            "file": str(functional_file) if functional_file else "",
            "annotated_rows": count_table_rows(functional_file) if functional_file else None,
        },
        "hydrolase_prediction": {
            "file": str(hydrolases_tsv),
            "predicted_hydrolase_rows": count_table_rows(hydrolases_tsv),
        },
        "quality_reports": {
            "quast_report": str(quast_file) if quast_file else "",
            "busco_summary": str(busco_file) if busco_file else "",
            "busco_preview": read_short_text(busco_file),
        },
    }

    return summary


def value_or_missing(value):
    if value is None or value == "":
        return "brak danych"
    return str(value)


def write_text_report(summary, output_file):
    lines = []
    lines.append("Raport końcowy analizy genomu grzyba")
    lines.append("====================================")
    lines.append("")
    lines.append(f"Data utworzenia raportu: {summary['generated_at']}")
    lines.append("")

    assembly = summary["assembly"]
    lines.append("Assemblacja genomu")
    lines.append("------------------")
    if assembly:
        lines.append(f"Plik genomu: {assembly['file']}")
        lines.append(f"Liczba kontigów: {assembly['contigs']}")
        lines.append(f"Łączna długość złożenia: {assembly['total_length']} bp")
        lines.append(f"N50: {assembly['n50']} bp")
        lines.append(f"Najdłuższy kontig: {assembly['longest_contig']} bp")
    else:
        lines.append("Brak pliku data/assemble_genome/latest_assembly.fasta")
    lines.append("")

    gene = summary["gene_prediction"]
    lines.append("Predykcja genów")
    lines.append("---------------")
    lines.append(f"Plik GFF3: {gene['gff']}")
    lines.append(f"Liczba przewidzianych genów: {value_or_missing(gene['predicted_genes'])}")
    lines.append(f"Plik białek: {gene['proteins']}")
    lines.append(f"Liczba przewidzianych białek: {value_or_missing(gene['predicted_proteins'])}")
    lines.append("")

    functional = summary["functional_annotation"]
    lines.append("Annotacja funkcjonalna")
    lines.append("----------------------")
    lines.append(f"Plik wynikowy: {value_or_missing(functional['file'])}")
    lines.append(f"Liczba wpisów annotacji: {value_or_missing(functional['annotated_rows'])}")
    lines.append("")

    hydro = summary["hydrolase_prediction"]
    lines.append("Predykcja hydrolaz")
    lines.append("------------------")
    lines.append(f"Plik wynikowy: {hydro['file']}")
    lines.append(f"Liczba kandydatów hydrolaz / wpisów: {value_or_missing(hydro['predicted_hydrolase_rows'])}")
    lines.append("")

    quality = summary["quality_reports"]
    lines.append("Kontrola jakości assemblacji")
    lines.append("----------------------------")
    lines.append(f"Raport QUAST: {value_or_missing(quality['quast_report'])}")
    lines.append(f"Podsumowanie BUSCO: {value_or_missing(quality['busco_summary'])}")

    if quality["busco_preview"]:
        lines.append("")
        lines.append("Podgląd BUSCO:")
        lines.extend(quality["busco_preview"])

    lines.append("")
    lines.append("Interpretacja")
    lines.append("-------------")
    lines.append("Ten raport zbiera najważniejsze wyniki z poszczególnych etapów. Szczegółowe raporty pozostają w katalogu results/, a pliki używane przez kolejne etapy w katalogu data/.")

    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Tworzy końcowy raport zbiorczy pipeline'u.")
    parser.add_argument("--output-dir", default="results/final_report")
    args = parser.parse_args()

    output_dir = PROJECT_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = build_summary()

    text_file = output_dir / "final_report.txt"
    json_file = output_dir / "final_report.json"

    write_text_report(summary, text_file)
    json_file.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Raport końcowy zapisano w: {text_file}")
    print(f"Dane raportu JSON zapisano w: {json_file}")


if __name__ == "__main__":
    main()
