#!/usr/bin/env python3
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def write_file(path, content):
    path = PROJECT_ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"Utworzono: {path}")


def main():
    now = datetime.now().isoformat(timespec="seconds")

    write_file(
        "data/assemble_genome/latest_assembly.fasta",
        """>demo_contig_1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGATCGATCGATCGA
>demo_contig_2
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAACGATCGATCGATCGATCGATCGATCGATCGA
"""
    )

    write_file(
        "data/predicted_genes/predicted_genes.gff3",
        """##gff-version 3
demo_contig_1\tAugustus\tgene\t1\t72\t.\t+\t.\tID=gene1
demo_contig_1\tAugustus\tmRNA\t1\t72\t.\t+\t.\tID=transcript1;Parent=gene1
demo_contig_1\tAugustus\tCDS\t1\t72\t.\t+\t0\tID=cds1;Parent=transcript1
demo_contig_2\tAugustus\tgene\t1\t72\t.\t-\t.\tID=gene2
demo_contig_2\tAugustus\tmRNA\t1\t72\t.\t-\t.\tID=transcript2;Parent=gene2
demo_contig_2\tAugustus\tCDS\t1\t72\t.\t-\t0\tID=cds2;Parent=transcript2
"""
    )

    write_file(
        "data/predicted_genes/predicted_proteins.faa",
        """>augustus_protein_1
MRTVASLADRSIDRSIDRS
>augustus_protein_2
MALWMRLLPLLALLALWGPG
"""
    )

    write_file(
        "results/functional_annotation/demo_functional_annotation.tsv",
        """protein_id\ttool\tbest_hit\tfunction\tevalue
augustus_protein_1\tDIAMOND\tsp|P001|HYDROLASE_DEMO\tputative hydrolase\t1e-40
augustus_protein_2\tInterProScan\tIPR000001\tmembrane protein domain\t2e-12
"""
    )

    write_file(
        "data/hydrolases/predicted_hydrolases.tsv",
        """protein_id\tsource\tprediction\tevidence
augustus_protein_1\tDIAMOND\tputative hydrolase\tSwiss-Prot hydrolase-like hit
augustus_protein_1\tHMMER\tglycoside hydrolase domain\tPfam-like demo domain
"""
    )

    write_file(
        "results/assembly_qc/demo/quast/report.tsv",
        """Assembly\t# contigs\tTotal length\tN50
demo\t2\t144\t72
"""
    )

    write_file(
        "results/assembly_qc/demo/busco_demo/short_summary_demo.txt",
        """# BUSCO demo summary
C:95.0%[S:94.0%,D:1.0%],F:3.0%,M:2.0%,n:100
95 Complete BUSCOs
3 Fragmented BUSCOs
2 Missing BUSCOs
"""
    )

    write_file(
        "results/demo_mode/demo_report.txt",
        f"""Tryb demonstracyjny
===================

Data utworzenia: {now}

Utworzono małe przykładowe wyniki imitujące działanie pipeline'u.
Pliki te służą do sprawdzenia działania GUI, raportu końcowego i struktury projektu.
"""
    )

    print("Tryb demonstracyjny zakończony.")
    print("Możesz teraz uruchomić scripts/build_final_report.py albo kafelek Raport końcowy w GUI.")


if __name__ == "__main__":
    main()
