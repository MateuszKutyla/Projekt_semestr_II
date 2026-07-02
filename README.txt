Pipeline do składania i analizy funkcjonalnej genomu grzyba
======================

Temat
-----

Złożenie genomu grzyba na podstawie danych sekwencjonowania Illumina oraz Oxford Nanopore Technologies (ONT), a następnie predykcja genów, annotacja funkcjonalna i predykcja hydrolaz.

Cel projektu
------------

Celem projektu jest przygotowanie powtarzalnego pipeline'u bioinformatycznego do:

1. sprawdzenia poprawności plików FASTQ,
2. oczyszczenia danych sekwencjonowania,
3. assemblacji de novo genomu,
4. predykcji genów,
5. annotacji funkcjonalnej białek,
6. predykcji potencjalnych hydrolaz.

System operacyjny
-----------------

Właściwe etapy analizy są przeznaczone do uruchamiania w systemie Linux.

GUI można testować na Windowsie, ale narzędzia bioinformatyczne używane w pipeline'ie, takie jak FastQC, MultiQC, fastp, Porechop, Filtlong, SPAdes, Flye, QUAST, BUSCO, Augustus, DIAMOND, eggNOG-mapper, InterProScan, HMMER, dbCAN, SignalP i DeepTMHMM, powinny być uruchamiane w Linuxie albo WSL.

Środowisko
----------

Plik środowiska znajduje się tutaj:

    envs/environment.yaml

Utworzenie środowiska:

    mamba env create -f envs/environment.yaml

Aktywacja środowiska przed uruchomieniem właściwego pipeline'u:

    mamba activate fungus-genome-pipeline

Aktualizacja istniejącego środowiska:

    mamba env update -f envs/environment.yaml --prune

Uwaga: SignalP nie jest instalowany przez envs/environment.yaml. Jeśli analiza peptydów sygnałowych ma być używana, SignalP trzeba zainstalować osobno zgodnie z instrukcją producenta i udostępnić komendę signalp6 w PATH.

Dane wejściowe
--------------

Surowe dane należy umieścić w katalogu:

    data/raw/

Aktualnie projekt zakłada pliki:

    data/raw/Unknown_CK982-002R0001_1.fq.gz
    data/raw/Unknown_CK982-002R0001_2.fq.gz
    data/raw/ONT_CK982-001N0001_raw.fq.gz

Surowych danych nie dodajemy do repozytorium Git.

Program konsolowy do walidacji danych
-------------------------------------

Pierwszy etap wykonuje program konsolowy:

    scripts/validate_reads.py

Ten etap wykonuje się przed uruchomieniem GUI. Pełny pipeline dostępny w GUI zaczyna się od assemblacji genomu i zakłada, że dane zostały już sprawdzone oraz zapisane w data/clean/.

Przykład uruchomienia:

    python3 scripts/validate_reads.py --illumina-r1 data/raw/Unknown_CK982-002R0001_1.fq.gz --illumina-r2 data/raw/Unknown_CK982-002R0001_2.fq.gz --ont data/raw/ONT_CK982-001N0001_raw.fq.gz --copy-to-clean

Program sprawdza:

- format FASTQ,
- obecność informacji o jakości,
- poprawność znaków sekwencji,
- zgodność długości sekwencji i jakości,
- zgodność par Illumina R1/R2.

Czyste dane lub dane po naprawie trafiają do:

    data/clean/

Główny program GUI
------------------

Główny program uruchamia się komendą:

    python3 scripts/genome_pipeline_gui.py

W Windowsie do testowania okien można użyć:

    python scripts\genome_pipeline_gui.py

Główne moduły GUI:

- Assemblacja de novo,
- Predykcja genów,
- Annotacja funkcjonalna,
- Predykcja hydrolaz,
- Pełny pipeline.

Assemblacja de novo
-------------------

Moduł uruchamia:

    scripts/run_denovo_assembly.py

Dostępne tryby:

1. Na bazie odczytów Illumina:
   - QC przed czyszczeniem: FastQC, MultiQC,
   - usuwanie adapterów i czyszczenie: fastp,
   - QC po czyszczeniu: FastQC, MultiQC,
   - assemblacja: SPAdes,
   - kontrola assemblacji: QUAST, BUSCO.

2. Na bazie odczytów ONT:
   - QC przed czyszczeniem: NanoPlot,
   - usuwanie adapterów: Porechop,
   - filtrowanie: Filtlong,
   - QC po czyszczeniu: NanoPlot,
   - assemblacja: Flye,
   - kontrola assemblacji: QUAST, BUSCO.

3. Assemblacja hybrydowa:
   - wykorzystuje odczyty Illumina i ONT,
   - Illumina jest czyszczona programem fastp,
   - ONT jest czyszczone Porechop i Filtlong,
   - assemblacja jest wykonywana przez SPAdes z opcją nanopore,
   - jakość złożenia jest oceniana przez QUAST i BUSCO.

Finalny genom jest kopiowany do:

    data/assemble_genome/latest_assembly.fasta

Predykcja genów
---------------

Moduł uruchamia:

    scripts/run_gene_prediction.py

Narzędzie:

    Augustus

Wejście domyślne:

    data/assemble_genome/latest_assembly.fasta

GUI pozwala wybrać:

- plik FASTA ze złożonym genomem,
- model referencyjny Augustusa,
- własny model Augustusa.

Wyniki:

    data/predicted_genes/predicted_genes.gff3
    data/predicted_genes/predicted_proteins.faa
    results/gene_prediction/augustus/gene_prediction_report.txt

Annotacja funkcjonalna
----------------------

Moduł uruchamia:

    scripts/run_functional_annotation.py

Wejście domyślne:

    data/predicted_genes/predicted_proteins.faa

Dostępne narzędzia:

- DIAMOND,
- eggNOG-mapper,
- InterProScan.

GUI pozwala uruchomić wybrane narzędzie albo wszystkie analizy funkcjonalne po kolei.
Z poziomu GUI można też pobrać i zbudować bazę DIAMOND do annotacji funkcjonalnej.

Wyniki:

    results/functional_annotation/

Predykcja hydrolaz
------------------

Moduł uruchamia:

    scripts/run_hydrolase_prediction.py

Wejście domyślne:

    data/predicted_genes/predicted_proteins.faa

Dostępne analizy:

- dbCAN,
- HMMER,
- DIAMOND względem bazy hydrolaz,
- SignalP,
- DeepTMHMM.

GUI pozwala zaznaczyć dowolny zestaw analiz albo wszystkie narzędzia.
Z poziomu GUI można przygotować bazę HMM dla HMMER oraz zbudować bazę DIAMOND używaną w etapie predykcji hydrolaz.
SignalP jest narzędziem opcjonalnym i wymaga osobnej instalacji poza środowiskiem Conda/Mamba.

Główny zbiorczy plik wynikowy:

    data/hydrolases/predicted_hydrolases.tsv

Raport:

    results/hydrolase_prediction/hydrolase_prediction_report.txt

Pełny pipeline
--------------

Główny program GUI zawiera przycisk:

    Pełny pipeline

Uruchamia on kolejno:

1. assemblację genomu,
2. predykcję genów,
3. annotację funkcjonalną,
4. predykcję hydrolaz.

Walidacja plików FASTQ nie jest uruchamiana przez ten przycisk. Należy wykonać ją wcześniej z konsoli programem scripts/validate_reads.py.

Skrypt konsolowy pełnego pipeline'u:

    scripts/run_full_pipeline.py

Log pełnego pipeline'u:

    logs/full_pipeline.log

Konfiguracja projektu
---------------------

Główne założenia i domyślne ścieżki są opisane w pliku:

    config/config.yaml

Plik konfiguracyjny porządkuje ustawienia projektu, ale część skryptów nadal korzysta z domyślnych ścieżek zapisanych bezpośrednio w kodzie. W praktyce najważniejsze pliki wejściowe powinny znajdować się w data/raw/ i data/clean/.

Dokument workflow
-----------------

Szczegółowy opis przebiegu analizy znajduje się w pliku:

    docs/workflow.txt

Struktura najważniejszych katalogów
-----------------------------------

    data/raw/                  surowe dane wejściowe
    data/clean/                dane oczyszczone
    data/assemble_genome/      finalny złożony genom
    data/predicted_genes/      geny i białka po predykcji
    data/functional_annotation/ dane pomocnicze annotacji funkcjonalnej
    data/hydrolases/           zbiorcze wyniki predykcji hydrolaz
    data/databases/            lokalne bazy DIAMOND, HMM i inne bazy pomocnicze
    data/example_results/      przykładowe pliki wynikowe pokazujące format rezultatów
    results/                   raporty i wyniki poszczególnych analiz
    logs/                      logi wykonania pipeline'u
    scripts/                   skrypty programu
    envs/                      definicja środowiska Conda/Mamba

Plan analizy
------------

1. Filtracja programem uruchamianym z konsoli: scripts/validate_reads.py.
2. Assemblacja genomu w programie GUI: scripts/genome_pipeline_gui.py.
3. Predykcja genów w programie GUI: scripts/genome_pipeline_gui.py.
4. Annotacja funkcjonalna w programie GUI: scripts/genome_pipeline_gui.py.
5. Predykcja hydrolaz w programie GUI: scripts/genome_pipeline_gui.py.
6. Opcjonalnie: uruchomienie całej analizy jednym przyciskiem Pełny pipeline.



