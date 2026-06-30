# Projekt semestralny II

## Temat

Złożenie genomu grzyba na podstawie danych sekwencjonowania Illumina oraz Oxford Nanopore Technologies (ONT).

## Cel projektu

Celem projektu jest przygotowanie powtarzalnego pipeline'u bioinformatycznego do kontroli jakości odczytów, oczyszczania danych, składania genomu, polishingu oraz oceny jakości otrzymanego złożenia.

## Dane wejściowe

Projekt zakłada wykorzystanie:

- krótkich odczytów Illumina paired-end,
- długich odczytów ONT.

Surowych danych nie dodajemy do repozytorium Git. Należy je trzymać lokalnie w `data/raw/`.

## Struktura repozytorium

- `config/` - konfiguracja projektu
- `data/raw/` - surowe dane sekwencjonowania
- `data/processed/` - oczyszczone dane
- `docs/` - dokumentacja i opis pipeline'u
- `envs/` - środowisko Conda/Mamba
- `results/` - wyniki analiz
- `scripts/` - skrypty pomocnicze
- `workflow/` - workflow analityczny

## Plan analizy

1. Kontrola jakości Illumina i ONT.
2. Oczyszczenie i filtrowanie odczytów.
3. Złożenie genomu z odczytów ONT.
4. Korekcja złożenia.
5. Polishing z użyciem Illumina.
6. Ocena jakości złożenia przy użyciu QUAST i BUSCO.
