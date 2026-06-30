#!/usr/bin/env python3
import argparse
import gzip
import json
from datetime import datetime
from pathlib import Path

VALID_BASES = set("ACGTNacgtn")
VALID_QUAL_MIN = 33
VALID_QUAL_MAX = 126


def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_fastq(path):
    with open_fastq(path) as handle:
        number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            number += 1
            yield number, header.rstrip(), seq.rstrip(), plus.rstrip(), qual.rstrip()


def validate_fastq(path, label, max_errors=20):
    result = {
        "label": label,
        "path": str(path),
        "exists": Path(path).exists(),
        "records": 0,
        "min_length": None,
        "max_length": None,
        "total_bases": 0,
        "has_quality": True,
        "errors": [],
        "warnings": []
    }

    if not result["exists"]:
        result["errors"].append("Plik nie istnieje.")
        return result

    try:
        for number, header, seq, plus, qual in read_fastq(path):
            if len(result["errors"]) >= max_errors:
                result["warnings"].append(f"Przerwano po {max_errors} bledach.")
                break

            result["records"] += 1
            result["total_bases"] += len(seq)
            result["min_length"] = len(seq) if result["min_length"] is None else min(result["min_length"], len(seq))
            result["max_length"] = len(seq) if result["max_length"] is None else max(result["max_length"], len(seq))

            if not header.startswith("@"):
                result["errors"].append(f"Rekord {number}: naglowek nie zaczyna sie od @.")
            if not plus.startswith("+"):
                result["errors"].append(f"Rekord {number}: separator nie zaczyna sie od +.")
            if not seq:
                result["errors"].append(f"Rekord {number}: pusta sekwencja.")
            if not qual:
                result["has_quality"] = False
                result["errors"].append(f"Rekord {number}: brak informacji o jakosci.")
            if len(seq) != len(qual):
                result["errors"].append(f"Rekord {number}: dlugosc sekwencji i jakosci jest rozna.")

            bad_bases = sorted(set(seq) - VALID_BASES)
            if bad_bases:
                result["errors"].append(f"Rekord {number}: niedozwolone znaki w sekwencji: {bad_bases}.")

            if any(ord(ch) < VALID_QUAL_MIN or ord(ch) > VALID_QUAL_MAX for ch in qual):
                result["errors"].append(f"Rekord {number}: niedozwolone znaki jakosci.")

    except gzip.BadGzipFile:
        result["errors"].append("Plik ma rozszerzenie .gz, ale nie jest poprawnym gzip.")
    except UnicodeDecodeError:
        result["errors"].append("Nie udalo sie odczytac pliku jako tekst FASTQ.")

    if result["records"] == 0 and not result["errors"]:
        result["errors"].append("Plik nie zawiera rekordow FASTQ.")

    result["average_length"] = round(result["total_bases"] / result["records"], 2) if result["records"] else None
    result["valid"] = len(result["errors"]) == 0
    return result


def validate_pairs(r1, r2, max_errors=20):
    result = {
        "checked_pairs": 0,
        "read_count_match": True,
        "read_length_mismatches": 0,
        "errors": [],
        "warnings": []
    }

    it1 = read_fastq(r1)
    it2 = read_fastq(r2)

    while True:
        rec1 = next(it1, None)
        rec2 = next(it2, None)

        if rec1 is None and rec2 is None:
            break
        if rec1 is None or rec2 is None:
            result["read_count_match"] = False
            result["errors"].append("Pliki Illumina R1 i R2 maja rozna liczbe odczytow.")
            break

        result["checked_pairs"] += 1
        if len(rec1[2]) != len(rec2[2]):
            result["read_length_mismatches"] += 1
            if len(result["warnings"]) < max_errors:
                result["warnings"].append(
                    f"Para {result['checked_pairs']}: R1 ma {len(rec1[2])} nt, R2 ma {len(rec2[2])} nt."
                )

    if result["read_length_mismatches"] > 0:
        result["errors"].append("Wykryto pary Illumina o roznej dlugosci odczytow.")

    result["valid"] = len(result["errors"]) == 0
    return result


def make_markdown(report):
    lines = [
        "# Raport kontroli poprawnosci plikow FASTQ",
        "",
        f"Data analizy: {report['created_at']}",
        f"Status koncowy: {'DANE POPRAWNE' if report['all_valid'] else 'WYKRYTO PROBLEMY'}",
        "",
        "## Co sprawdzono",
        "",
        "- format rekordow FASTQ",
        "- obecnosc sekwencji i informacji o jakosci",
        "- zgodnosc dlugosci sekwencji i jakosci",
        "- dozwolone znaki w sekwencji DNA",
        "- dozwolone znaki jakosci",
        "- zgodnosc liczby odczytow Illumina R1/R2",
        "- zgodnosc dlugosci sparowanych odczytow Illumina",
        ""
    ]

    for item in report["files"]:
        lines += [
            f"## {item['label']}",
            "",
            f"- plik: {item['path']}",
            f"- liczba rekordow: {item['records']}",
            f"- minimalna dlugosc: {item['min_length']}",
            f"- maksymalna dlugosc: {item['max_length']}",
            f"- srednia dlugosc: {item['average_length']}",
            f"- informacja o jakosci: {'tak' if item['has_quality'] else 'nie'}",
            f"- status: {'OK' if item['valid'] else 'BLAD'}",
            ""
        ]
        if item["errors"]:
            lines.append("Bledy:")
            lines += [f"- {error}" for error in item["errors"]]
            lines.append("")

    pair = report["illumina_pairs"]
    lines += [
        "## Pary Illumina",
        "",
        f"- sprawdzone pary: {pair['checked_pairs']}",
        f"- taka sama liczba odczytow R1/R2: {'tak' if pair['read_count_match'] else 'nie'}",
        f"- pary o roznej dlugosci: {pair['read_length_mismatches']}",
        f"- status: {'OK' if pair['valid'] else 'BLAD'}",
        ""
    ]

    if pair["errors"]:
        lines.append("Bledy:")
        lines += [f"- {error}" for error in pair["errors"]]
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Walidacja plikow FASTQ Illumina i ONT.")
    parser.add_argument("--illumina-r1", required=True)
    parser.add_argument("--illumina-r2", required=True)
    parser.add_argument("--ont", required=True)
    parser.add_argument("--report-dir", default="results/read_validation")
    parser.add_argument("--clean-dir", default="data/clean")
    args = parser.parse_args()

    report_dir = Path(args.report_dir)
    clean_dir = Path(args.clean_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)

    files = [
        validate_fastq(args.illumina_r1, "Illumina R1"),
        validate_fastq(args.illumina_r2, "Illumina R2"),
        validate_fastq(args.ont, "ONT")
    ]
    pairs = validate_pairs(args.illumina_r1, args.illumina_r2)

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "files": files,
        "illumina_pairs": pairs,
        "clean_directory": str(clean_dir),
        "all_valid": all(item["valid"] for item in files) and pairs["valid"]
    }

    (report_dir / "read_validation_report.md").write_text(make_markdown(report), encoding="utf-8")
    (report_dir / "read_validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8"
    )

    if report["all_valid"]:
        print(f"Dane sa poprawne. Raport zapisano w: {report_dir}")
        return 0

    print(f"Wykryto problemy. Raport zapisano w: {report_dir}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
