# Projekt semestralny II

## Temat

Złożenie genomu grzyba na podstawie danych sekwencjonowania Illumina oraz Oxford Nanopore Technologies (ONT).

## Cel projektu

Celem projektu jest przygotowanie powtarzalnego pipeline'u bioinformatycznego do kontroli jakości odczytów, oczyszczania danych, składania genomu, polishingu oraz oceny jakości otrzymanego złożenia.

## System operacyjny

Wszystkie etapy analizy są przeznaczone do uruchamiania w systemie Linux. Komendy w projekcie będą zapisywane w składni bash, a ścieżki do plików będą podawane w formacie linuxowym, np. data/raw/illumina_R1.fastq.gz.

Projekt można uruchamiać na przykład w:

- natywnym systemie Linux,
- WSL na Windowsie,
- maszynie wirtualnej z Linuxem,
- serwerze obliczeniowym z Linuxem.

Zalecane jest użycie środowiska Conda/Mamba zdefiniowanego w pliku envs/environment.yaml.

Przykładowe przygotowanie środowiska:

    mamba env create -f envs/environment.yaml
    mamba activate fungus-genome-pipeline


Przed uruchomieniem właściwego programu do składania genomu trzeba aktywować środowisko z narzędziami bioinformatycznymi:

    mamba activate fungus-genome-pipeline

Dopiero po aktywacji środowiska należy uruchamiać GUI albo skrypty pipeline'u, ponieważ etapy assemblacji korzystają z narzędzi takich jak FastQC, MultiQC, fastp, Porechop, Filtlong, SPAdes, Flye, QUAST i BUSCO.

## Dane wejściowe

Projekt zakłada wykorzystanie:

- krótkich odczytów Illumina paired-end,
- długich odczytów ONT.

Surowych danych nie dodajemy do repozytorium Git. Należy je trzymać lokalnie w data/raw/.

## Struktura repozytorium

- config/ - konfiguracja projektu
- data/raw/ - surowe dane sekwencjonowania
- data/processed/ - oczyszczone dane
- docs/ - dokumentacja i opis pipeline'u
- envs/ - środowisko Conda/Mamba
- results/ - wyniki analiz
- scripts/ - skrypty pomocnicze
- workflow/ - workflow analityczny

## Plan analizy

1. Filtracja programem uruchamianym z konsoli: scripts/validate_reads.py.
2. Assemblacja genomu w programie GUI: scripts/genome_pipeline_gui.py.
3. Predykcja genów w programie GUI: scripts/genome_pipeline_gui.py.
4. Annotacja funkcjonalna w programie GUI: scripts/genome_pipeline_gui.py.
5. Predykcja hydrolaz w programie GUI: scripts/genome_pipeline_gui.py.



