# Walidacja plików FASTQ

Pierwszy etap projektu sprawdza poprawność plików ONT i Illumina.

Program kontroluje:

- format FASTQ,
- dozwolone znaki w sekwencjach,
- obecność informacji o jakości,
- zgodność długości sekwencji i jakości,
- zgodność liczby sparowanych odczytów Illumina R1/R2,
- zgodność długości sparowanych odczytów Illumina.

## Uruchomienie w Linuxie

Polecenie uruchomienia:

python3 scripts/validate_reads.py --illumina-r1 data/raw/illumina_R1.fastq.gz --illumina-r2 data/raw/illumina_R2.fastq.gz --ont data/raw/ont_reads.fastq.gz

Raport tekstowy zostanie zapisany jako results/read_validation/read_validation_report.txt.

Folder data/clean/ jest przeznaczony na późniejsze czyste dane.

