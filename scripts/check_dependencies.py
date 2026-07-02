#!/usr/bin/env python3
import argparse
import json
import shutil
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

TOOLS = {
    "Kontrola jakości i czyszczenie odczytów": [
        ("FastQC", "fastqc", True),
        ("MultiQC", "multiqc", True),
        ("NanoPlot", "NanoPlot", True),
        ("fastp", "fastp", True),
        ("Filtlong", "filtlong", True),
        ("Porechop", "porechop", True),
    ],
    "Assemblacja genomu": [
        ("SPAdes", "spades.py", True),
        ("Flye", "flye", True),
    ],
    "Ocena jakości assemblacji": [
        ("QUAST", "quast", True),
        ("BUSCO", "busco", True),
    ],
    "Predykcja genów": [
        ("Augustus", "augustus", True),
    ],
    "Annotacja funkcjonalna": [
        ("DIAMOND", "diamond", True),
        ("eggNOG-mapper", "emapper.py", False),
        ("InterProScan", "interproscan.sh", False),
    ],
    "Predykcja hydrolaz": [
        ("HMMER hmmscan", "hmmscan", True),
        ("HMMER hmmpress", "hmmpress", True),
        ("dbCAN", "run_dbcan.py", False),
        ("SignalP", "signalp", False),
        ("BioLib", "biolib", False),
    ],
}


def check_tool(command):
    found_path = shutil.which(command)
    return {
        "command": command,
        "available": found_path is not None,
        "path": found_path or "",
    }


def build_report(results):
    lines = []
    lines.append("Raport sprawdzenia narzędzi")
    lines.append("===========================")
    lines.append("")
    lines.append(f"Data sprawdzenia: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")

    missing_required = 0
    missing_optional = 0

    for group, tools in results.items():
        lines.append(group)
        lines.append("-" * len(group))

        for item in tools:
            status = "OK" if item["available"] else "BRAK"
            requirement = "wymagane" if item["required"] else "opcjonalne"
            location = f" ({item['path']})" if item["path"] else ""
            lines.append(f"- {item['name']}: {status} [{requirement}] - komenda: {item['command']}{location}")

            if not item["available"] and item["required"]:
                missing_required += 1
            if not item["available"] and not item["required"]:
                missing_optional += 1

        lines.append("")

    lines.append("Podsumowanie")
    lines.append("------------")
    lines.append(f"Brakujące narzędzia wymagane: {missing_required}")
    lines.append(f"Brakujące narzędzia opcjonalne: {missing_optional}")

    if missing_required == 0:
        lines.append("Status: środowisko ma wszystkie wymagane narzędzia do głównych etapów analizy.")
    else:
        lines.append("Status: środowisko wymaga uzupełnienia przed pełną analizą.")

    return "\n".join(lines) + "\n", missing_required, missing_optional


def main():
    parser = argparse.ArgumentParser(description="Sprawdza dostępność narzędzi wymaganych w pipeline.")
    parser.add_argument("--output-dir", default="results/dependency_check")
    parser.add_argument("--strict", action="store_true", help="Zwraca kod błędu, jeśli brakuje narzędzi wymaganych.")
    args = parser.parse_args()

    output_dir = PROJECT_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for group, tools in TOOLS.items():
        results[group] = []
        for name, command, required in tools:
            checked = check_tool(command)
            checked["name"] = name
            checked["required"] = required
            results[group].append(checked)

    report_text, missing_required, missing_optional = build_report(results)

    report_file = output_dir / "dependency_check_report.txt"
    json_file = output_dir / "dependency_check_report.json"

    report_file.write_text(report_text, encoding="utf-8")
    json_file.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")

    print(report_text)
    print(f"Raport zapisano w: {report_file}")

    if args.strict and missing_required > 0:
        raise SystemExit(1)

    raise SystemExit(0)


if __name__ == "__main__":
    main()
