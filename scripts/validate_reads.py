#!/usr/bin/env python3
import argparse
import gzip
import json
import shutil
from datetime import datetime
from pathlib import Path

VALID_BASES = set("ACGTNacgtn")
VALID_QUAL_MIN = 33
VALID_QUAL_MAX = 126


def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def open_output_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "wt", encoding="utf-8")


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


def validate_fastq(path, label, max_errors=20, max_records=None):
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
            if max_records is not None and result["records"] >= max_records:
                result["warnings"].append(f"Tryb testowy: sprawdzono pierwsze {max_records} rekordow.")
                break
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


def validate_pairs(r1, r2, max_errors=20, max_records=None):
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
        if max_records is not None and result["checked_pairs"] >= max_records:
            result["warnings"].append(f"Tryb testowy: sprawdzono pierwsze {max_records} par odczytow.")
            break
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


def clean_output_path(input_path, clean_dir, suffix="clean"):
    path = Path(input_path)
    name = path.name
    if name.endswith(".fastq.gz"):
        out_name = name.replace(".fastq.gz", f".{suffix}.fastq.gz")
    elif name.endswith(".fq.gz"):
        out_name = name.replace(".fq.gz", f".{suffix}.fq.gz")
    elif name.endswith(".fastq"):
        out_name = name.replace(".fastq", f".{suffix}.fastq")
    elif name.endswith(".fq"):
        out_name = name.replace(".fq", f".{suffix}.fq")
    else:
        out_name = f"{name}.{suffix}.fq"
    return Path(clean_dir) / out_name


def repair_record(header, seq, plus, qual):
    changes = []

    if not header.startswith("@"):
        header = "@" + header.lstrip("@")
        changes.append("naprawiono naglowek")

    seq = "".join(base.upper() if base in VALID_BASES else "N" for base in seq)
    if "N" in seq:
        changes.append("zamieniono niedozwolone zasady na N")

    plus = "+"

    qual = "".join(ch if VALID_QUAL_MIN <= ord(ch) <= VALID_QUAL_MAX else "!" for ch in qual)
    if len(seq) != len(qual):
        common_length = min(len(seq), len(qual))
        seq = seq[:common_length]
        qual = qual[:common_length]
        changes.append("przycieto sekwencje i jakosc do wspolnej dlugosci")

    if not seq or not qual:
        return None, changes + ["pominieto pusty rekord"]

    return (header, seq, plus, qual), changes


def write_record(handle, record):
    header, seq, plus, qual = record
    handle.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


def repair_single_fastq(input_path, output_path, max_records=None):
    stats = {"input": str(input_path), "output": str(output_path), "written": 0, "skipped": 0, "changes": []}

    with open_output_fastq(output_path) as output:
        for number, header, seq, plus, qual in read_fastq(input_path):
            if max_records is not None and stats["written"] >= max_records:
                stats["changes"].append(f"tryb testowy: zapisano pierwsze {max_records} rekordow")
                break

            repaired, changes = repair_record(header, seq, plus, qual)
            stats["changes"].extend([f"rekord {number}: {change}" for change in changes])

            if repaired is None:
                stats["skipped"] += 1
                continue

            write_record(output, repaired)
            stats["written"] += 1

    return stats


def repair_illumina_pairs(r1_path, r2_path, clean_dir, max_records=None):
    r1_output = clean_output_path(r1_path, clean_dir)
    r2_output = clean_output_path(r2_path, clean_dir)
    stats = {
        "input_r1": str(r1_path),
        "input_r2": str(r2_path),
        "output_r1": str(r1_output),
        "output_r2": str(r2_output),
        "written_pairs": 0,
        "skipped_pairs": 0,
        "changes": []
    }

    with open_output_fastq(r1_output) as out1, open_output_fastq(r2_output) as out2:
        it1 = read_fastq(r1_path)
        it2 = read_fastq(r2_path)

        while True:
            if max_records is not None and stats["written_pairs"] >= max_records:
                stats["changes"].append(f"tryb testowy: zapisano pierwsze {max_records} par")
                break

            rec1 = next(it1, None)
            rec2 = next(it2, None)

            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                stats["changes"].append("pominieto niesparowany odczyt na koncu pliku")
                break

            number = rec1[0]
            repaired1, changes1 = repair_record(rec1[1], rec1[2], rec1[3], rec1[4])
            repaired2, changes2 = repair_record(rec2[1], rec2[2], rec2[3], rec2[4])

            if repaired1 is None or repaired2 is None:
                stats["skipped_pairs"] += 1
                stats["changes"].append(f"para {number}: pominieto pare z pustym rekordem")
                continue

            common_length = min(len(repaired1[1]), len(repaired2[1]))
            repaired1 = (repaired1[0], repaired1[1][:common_length], repaired1[2], repaired1[3][:common_length])
            repaired2 = (repaired2[0], repaired2[1][:common_length], repaired2[2], repaired2[3][:common_length])

            if changes1 or changes2 or len(rec1[2]) != len(rec2[2]):
                stats["changes"].append(f"para {number}: naprawiono i zsynchronizowano odczyty")

            write_record(out1, repaired1)
            write_record(out2, repaired2)
            stats["written_pairs"] += 1

    return stats


def copy_valid_files_to_clean(files, clean_dir):
    actions = []
    for input_path in files:
        input_path = Path(input_path)
        output_path = Path(clean_dir) / input_path.name
        shutil.copy2(input_path, output_path)
        actions.append({"action": "copy", "input": str(input_path), "output": str(output_path)})
    return actions


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
        if item["warnings"]:
            lines.append("Ostrzezenia:")
            lines += [f"- {warning}" for warning in item["warnings"]]
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
    if pair["warnings"]:
        lines.append("Ostrzezenia:")
        lines += [f"- {warning}" for warning in pair["warnings"]]
        lines.append("")

    if report["clean_actions"]:
        lines += ["## Dzialania na danych czystych", ""]
        for action in report["clean_actions"]:
            lines.append(f"- {action}")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Walidacja plikow FASTQ Illumina i ONT.")
    parser.add_argument("--illumina-r1", required=True)
    parser.add_argument("--illumina-r2", required=True)
    parser.add_argument("--ont", required=True)
    parser.add_argument("--report-dir", default="results/read_validation")
    parser.add_argument("--clean-dir", default="data/clean")
    parser.add_argument("--max-records", type=int, default=None, help="Tryb testowy: sprawdza tylko podana liczbe rekordow z kazdego pliku.")
    parser.add_argument("--copy-to-clean", action="store_true", help="Kopiuje poprawne pliki do data/clean po udanej walidacji.")
    parser.add_argument("--repair-to-clean", action="store_true", help="Tworzy naprawione pliki FASTQ w data/clean do kolejnych etapow.")
    args = parser.parse_args()

    report_dir = Path(args.report_dir)
    clean_dir = Path(args.clean_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)

    files = [
        validate_fastq(args.illumina_r1, "Illumina R1", max_records=args.max_records),
        validate_fastq(args.illumina_r2, "Illumina R2", max_records=args.max_records),
        validate_fastq(args.ont, "ONT", max_records=args.max_records)
    ]
    pairs = validate_pairs(args.illumina_r1, args.illumina_r2, max_records=args.max_records)

    report = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "files": files,
        "illumina_pairs": pairs,
        "clean_directory": str(clean_dir),
        "all_valid": all(item["valid"] for item in files) and pairs["valid"],
        "clean_actions": []
    }

    input_files = [args.illumina_r1, args.illumina_r2, args.ont]

    if args.copy_to_clean:
        if report["all_valid"]:
            actions = copy_valid_files_to_clean(input_files, clean_dir)
            report["clean_actions"].extend([f"skopiowano {item['input']} do {item['output']}" for item in actions])
        else:
            report["clean_actions"].append("nie skopiowano plikow, poniewaz walidacja wykryla problemy")

    if args.repair_to_clean:
        pair_repair = repair_illumina_pairs(args.illumina_r1, args.illumina_r2, clean_dir, args.max_records)
        ont_output = clean_output_path(args.ont, clean_dir)
        ont_repair = repair_single_fastq(args.ont, ont_output, args.max_records)

        report["repair"] = {
            "illumina_pairs": pair_repair,
            "ont": ont_repair
        }
        report["clean_actions"].append(f"zapisano naprawione R1: {pair_repair['output_r1']}")
        report["clean_actions"].append(f"zapisano naprawione R2: {pair_repair['output_r2']}")
        report["clean_actions"].append(f"zapisano naprawione ONT: {ont_repair['output']}")

    (report_dir / "read_validation_report.txt").write_text(make_markdown(report), encoding="utf-8")
    (report_dir / "read_validation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8"
    )

    if report["all_valid"]:
        print(f"Dane sa poprawne. Raport zapisano w: {report_dir}")
        return 0

    if args.repair_to_clean:
        print(f"Wykryto problemy, ale zapisano naprawione pliki w: {clean_dir}")
        print(f"Raport zapisano w: {report_dir}")
        return 0

    print(f"Wykryto problemy. Raport zapisano w: {report_dir}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
