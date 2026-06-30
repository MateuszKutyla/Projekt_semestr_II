#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_R1 = "data/raw/Unknown_CK982-002R0001_1.fq.gz"
RAW_R2 = "data/raw/Unknown_CK982-002R0001_2.fq.gz"
RAW_ONT = "data/raw/ONT_CK982-001N0001_raw.fq.gz"

CLEAN_R1 = "data/clean/illumina_R1.clean.fq.gz"
CLEAN_R2 = "data/clean/illumina_R2.clean.fq.gz"
CLEAN_ONT = "data/clean/ont.adapters_removed.fq.gz"
FILTERED_ONT = "data/clean/ont.clean.fq"


def path(relative_path):
    return PROJECT_ROOT / relative_path


def ensure_dirs():
    for directory in [
        "data/clean",
        "logs",
        "results/qc/illumina_before",
        "results/qc/illumina_after",
        "results/qc/ont_before",
        "results/qc/ont_after",
        "results/assembly",
        "results/assembly_qc"
    ]:
        path(directory).mkdir(parents=True, exist_ok=True)


def run_step(name, command, log_handle, stdout_file=None):
    log_handle.write("\n" + "=" * 80 + "\n")
    log_handle.write(name + "\n")
    log_handle.write(" ".join(str(item) for item in command) + "\n")
    log_handle.write("=" * 80 + "\n")
    log_handle.flush()

    if stdout_file:
        with open(stdout_file, "w", encoding="utf-8") as output:
            process = subprocess.run(
                command,
                cwd=PROJECT_ROOT,
                stdout=output,
                stderr=log_handle,
                text=True
            )
    else:
        process = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True
        )

    if process.returncode != 0:
        raise RuntimeError(f"Etap zakonczyl sie bledem: {name}")


def run_fastqc_before_after_illumina(log_handle):
    run_step(
        "QC Illumina przed usunieciem adapterow - FastQC",
        ["fastqc", path(RAW_R1), path(RAW_R2), "-o", path("results/qc/illumina_before")],
        log_handle
    )
    run_step(
        "Raport zbiorczy Illumina przed czyszczeniem - MultiQC",
        ["multiqc", path("results/qc/illumina_before"), "-o", path("results/qc/illumina_before")],
        log_handle
    )

    run_step(
        "Usuniecie adapterow i czyszczenie Illumina - fastp",
        [
            "fastp",
            "-i", path(RAW_R1),
            "-I", path(RAW_R2),
            "-o", path(CLEAN_R1),
            "-O", path(CLEAN_R2),
            "--html", path("results/qc/illumina_after/fastp.html"),
            "--json", path("results/qc/illumina_after/fastp.json")
        ],
        log_handle
    )

    run_step(
        "QC Illumina po usunieciu adapterow - FastQC",
        ["fastqc", path(CLEAN_R1), path(CLEAN_R2), "-o", path("results/qc/illumina_after")],
        log_handle
    )
    run_step(
        "Raport zbiorczy Illumina po czyszczeniu - MultiQC",
        ["multiqc", path("results/qc/illumina_after"), "-o", path("results/qc/illumina_after")],
        log_handle
    )


def run_qc_before_after_ont(log_handle):
    run_step(
        "QC ONT przed usunieciem adapterow - NanoPlot",
        ["NanoPlot", "--fastq", path(RAW_ONT), "--outdir", path("results/qc/ont_before")],
        log_handle
    )

    run_step(
        "Usuniecie adapterow ONT - Porechop",
        ["porechop", "-i", path(RAW_ONT), "-o", path(CLEAN_ONT)],
        log_handle
    )

    run_step(
        "Filtrowanie ONT - Filtlong",
        ["filtlong", "--min_length", "1000", path(CLEAN_ONT)],
        log_handle,
        stdout_file=path(FILTERED_ONT)
    )

    run_step(
        "QC ONT po czyszczeniu - NanoPlot",
        ["NanoPlot", "--fastq", path(FILTERED_ONT), "--outdir", path("results/qc/ont_after")],
        log_handle
    )


def run_quast_busco(assembly_file, label, lineage, threads, log_handle):
    qc_dir = path(f"results/assembly_qc/{label}")
    qc_dir.mkdir(parents=True, exist_ok=True)

    run_step(
        f"Kontrola assemblacji {label} - QUAST",
        ["quast", assembly_file, "-o", qc_dir / "quast"],
        log_handle
    )

    run_step(
        f"Kontrola kompletnosci {label} - BUSCO",
        [
            "busco",
            "-i", assembly_file,
            "-o", f"busco_{label}",
            "--out_path", qc_dir,
            "-m", "genome",
            "-l", lineage,
            "-c", str(threads)
        ],
        log_handle
    )


def run_illumina_pipeline(args, log_handle):
    run_fastqc_before_after_illumina(log_handle)

    out_dir = path("results/assembly/spades_illumina")
    run_step(
        "Assemblacja de novo Illumina - SPAdes",
        ["spades.py", "-1", path(CLEAN_R1), "-2", path(CLEAN_R2), "-o", out_dir, "-t", str(args.threads)],
        log_handle
    )

    run_quast_busco(out_dir / "contigs.fasta", "spades_illumina", args.busco_lineage, args.threads, log_handle)


def run_ont_pipeline(args, log_handle):
    run_qc_before_after_ont(log_handle)

    out_dir = path("results/assembly/flye_ont")
    run_step(
        "Assemblacja de novo ONT - Flye",
        ["flye", "--nano-raw", path(FILTERED_ONT), "--out-dir", out_dir, "--threads", str(args.threads)],
        log_handle
    )

    run_quast_busco(out_dir / "assembly.fasta", "flye_ont", args.busco_lineage, args.threads, log_handle)


def run_hybrid_pipeline(args, log_handle):
    run_fastqc_before_after_illumina(log_handle)
    run_qc_before_after_ont(log_handle)

    out_dir = path("results/assembly/spades_hybrid")
    run_step(
        "Assemblacja hybrydowa Illumina + ONT - SPAdes",
        [
            "spades.py",
            "-1", path(CLEAN_R1),
            "-2", path(CLEAN_R2),
            "--nanopore", path(FILTERED_ONT),
            "-o", out_dir,
            "-t", str(args.threads)
        ],
        log_handle
    )

    run_quast_busco(out_dir / "contigs.fasta", "spades_hybrid", args.busco_lineage, args.threads, log_handle)


def main():
    parser = argparse.ArgumentParser(description="Pipeline assemblacji de novo z QC, czyszczeniem, QUAST i BUSCO.")
    parser.add_argument("--mode", choices=["illumina", "ont", "hybrid"], required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--busco-lineage", default="fungi_odb10")
    args = parser.parse_args()

    ensure_dirs()
    log_file = path(f"logs/assembly_{args.mode}.log")

    with open(log_file, "w", encoding="utf-8") as log_handle:
        log_handle.write(f"Tryb assemblacji: {args.mode}\n")

        if args.mode == "illumina":
            run_illumina_pipeline(args, log_handle)
        elif args.mode == "ont":
            run_ont_pipeline(args, log_handle)
        elif args.mode == "hybrid":
            run_hybrid_pipeline(args, log_handle)

        log_handle.write("\nPipeline zakonczony powodzeniem.\n")

    print(f"Pipeline zakonczony. Log zapisano w: {log_file}")


if __name__ == "__main__":
    main()
