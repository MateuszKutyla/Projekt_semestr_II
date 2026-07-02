#!/usr/bin/env python3
import subprocess
import sys
import threading
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

PROJECT_ROOT = Path(__file__).resolve().parents[1]


class GenomePipelineApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Projekt semestralny II - analiza genomu grzyba")
        self.geometry("1040x700")
        self.minsize(900, 600)

        self.project_dir = tk.StringVar(value=str(PROJECT_ROOT))
        self.status_text = tk.StringVar(value="Gotowe do pracy")
        self.current_process = None

        self.configure(bg="#f4f6f8")
        self.build_ui()

    def build_ui(self):
        header = tk.Frame(self, bg="#1f2937", padx=22, pady=18)
        header.pack(fill="x")

        tk.Label(header, text="Analiza genomu grzyba", bg="#1f2937", fg="white", font=("Segoe UI", 20, "bold")).pack(anchor="center")

        main = ttk.Frame(self, padding=20)
        main.pack(fill="both", expand=True)

        ttk.Label(main, text="Folder projektu").pack(anchor="w")
        path_row = ttk.Frame(main)
        path_row.pack(fill="x", pady=(6, 16))

        ttk.Entry(path_row, textvariable=self.project_dir).pack(side="left", fill="x", expand=True)
        ttk.Button(path_row, text="Wybierz", command=self.choose_project_dir).pack(side="left", padx=(8, 0))

        button_grid = ttk.Frame(main)
        button_grid.pack(fill="x")
        button_grid.columnconfigure(0, weight=1)
        button_grid.columnconfigure(1, weight=1)
        button_grid.columnconfigure(2, weight=1)

        self.add_module(button_grid, 0, 0, "Assemblacja de novo", "Wybór typu assemblacji: Illumina, ONT albo hybrydowa.", self.open_denovo_window)
        self.add_module(button_grid, 0, 1, "Predykcja genów", "Moduł do wykrywania genów w złożonym genomie.", self.run_gene_prediction)
        self.add_module(button_grid, 1, 0, "Annotacja funkcjonalna", "Moduł do przypisywania funkcji przewidywanym genom i białkom.", self.run_annotation)
        self.add_module(button_grid, 1, 1, "Predykcja hydrolaz", "Moduł do wyszukiwania potencjalnych enzymów hydrolitycznych.", self.run_hydrolases)
        self.add_module(button_grid, 2, 0, "Pełny pipeline", "Uruchamia assemblację, predykcję genów, annotację funkcjonalną i predykcję hydrolaz.", self.open_full_pipeline_window)

        ttk.Label(main, text="Log programu").pack(anchor="w", pady=(18, 6))
        self.log = tk.Text(main, height=12, bg="#111827", fg="#e5e7eb", insertbackground="white", wrap="word")
        self.log.pack(fill="both", expand=True)

        bottom = ttk.Frame(self, padding=(12, 8))
        bottom.pack(fill="x")
        ttk.Label(bottom, textvariable=self.status_text).pack(side="left")
        ttk.Button(bottom, text="Wyczyść log", command=self.clear_log).pack(side="right")

        self.write_log("Uruchomiono graficzny rdzeń projektu.")

    def add_module(self, parent, row, column, title, description, command):
        frame = ttk.Frame(parent, padding=18, relief="ridge", cursor="hand2")
        frame.grid(row=row, column=column, sticky="nsew", padx=8, pady=8)

        title_label = ttk.Label(frame, text=title, font=("Segoe UI", 15, "bold"), cursor="hand2")
        title_label.pack(anchor="w")

        description_label = ttk.Label(frame, text=description, wraplength=390, cursor="hand2")
        description_label.pack(anchor="w", pady=(8, 0))

        frame.bind("<Button-1>", lambda event: command())
        title_label.bind("<Button-1>", lambda event: command())
        description_label.bind("<Button-1>", lambda event: command())

    def choose_project_dir(self):
        selected = filedialog.askdirectory(initialdir=self.project_dir.get(), title="Wybierz folder projektu")
        if selected:
            self.project_dir.set(selected)
            self.write_log(f"Ustawiono folder projektu: {selected}")

    def project_path(self, relative_path):
        return Path(self.project_dir.get()) / relative_path

    def first_existing(self, paths):
        for path in paths:
            full_path = self.project_path(path)
            if full_path.exists():
                return str(full_path)
        return str(self.project_path(paths[0]))

    def open_full_pipeline_window(self):
        window = tk.Toplevel(self)
        window.title("Pełny pipeline")
        window.geometry("860x720")
        window.minsize(780, 640)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Pełny pipeline analizy genomu", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(
            frame,
            text="Ten tryb wykona kolejno: assemblację genomu, predykcję genów, annotację funkcjonalną i predykcję hydrolaz.",
            wraplength=800
        ).pack(anchor="w", pady=(6, 16))

        assembly_mode = tk.StringVar(value="hybrid")
        species = tk.StringVar(value="aspergillus_nidulans")
        custom_species = tk.StringVar(value="")
        functional_diamond_db = tk.StringVar(value="")
        hydrolase_hmm_db = tk.StringVar(value="")
        hydrolase_diamond_db = tk.StringVar(value="")
        threads = tk.StringVar(value="8")

        ttk.Label(frame, text="Tryb assemblacji:").pack(anchor="w")
        ttk.Combobox(
            frame,
            textvariable=assembly_mode,
            values=["illumina", "ont", "hybrid"],
            state="readonly"
        ).pack(fill="x", pady=(6, 10))

        ttk.Label(frame, text="Model Augustusa:").pack(anchor="w")
        ttk.Combobox(
            frame,
            textvariable=species,
            values=[
                "aspergillus_nidulans",
                "botrytis_cinerea",
                "candida_albicans",
                "fusarium_graminearum",
                "neurospora_crassa",
                "saccharomyces_cerevisiae_S288C",
                "ustilago_maydis",
                "yarrowia_lipolytica",
                "własny model"
            ],
            state="readonly"
        ).pack(fill="x", pady=(6, 8))

        ttk.Label(frame, text="Własny model Augustusa:").pack(anchor="w")
        custom_row = ttk.Frame(frame)
        custom_row.pack(fill="x", pady=(6, 10))
        ttk.Entry(custom_row, textvariable=custom_species).pack(side="left", fill="x", expand=True)
        ttk.Button(
            custom_row,
            text="Wybierz",
            command=lambda: self.choose_augustus_model_file(custom_species)
        ).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Baza DIAMOND do annotacji funkcjonalnej:").pack(anchor="w")
        row1 = ttk.Frame(frame)
        row1.pack(fill="x", pady=(6, 8))
        ttk.Entry(row1, textvariable=functional_diamond_db).pack(side="left", fill="x", expand=True)
        ttk.Button(row1, text="Wybierz", command=lambda: self.choose_diamond_database(functional_diamond_db)).pack(side="left", padx=(8, 0))
        ttk.Button(row1, text="Zbuduj", command=lambda: self.build_diamond_database(functional_diamond_db, "Zbuduj bazę DIAMOND do annotacji funkcjonalnej")).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Baza HMM do predykcji hydrolaz:").pack(anchor="w")
        row2 = ttk.Frame(frame)
        row2.pack(fill="x", pady=(6, 8))
        ttk.Entry(row2, textvariable=hydrolase_hmm_db).pack(side="left", fill="x", expand=True)
        ttk.Button(row2, text="Wybierz", command=lambda: self.choose_hmm_database(hydrolase_hmm_db)).pack(side="left", padx=(8, 0))
        ttk.Button(row2, text="Przygotuj HMM", command=lambda: self.prepare_hmm_database(hydrolase_hmm_db)).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Baza DIAMOND hydrolaz:").pack(anchor="w")
        row3 = ttk.Frame(frame)
        row3.pack(fill="x", pady=(6, 8))
        ttk.Entry(row3, textvariable=hydrolase_diamond_db).pack(side="left", fill="x", expand=True)
        ttk.Button(row3, text="Wybierz", command=lambda: self.choose_diamond_database(hydrolase_diamond_db)).pack(side="left", padx=(8, 0))
        ttk.Button(row3, text="Zbuduj", command=lambda: self.build_diamond_database(hydrolase_diamond_db, "Zbuduj bazę DIAMOND hydrolaz")).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Liczba wątków:").pack(anchor="w")
        ttk.Entry(frame, textvariable=threads).pack(fill="x", pady=(6, 14))

        ttk.Button(
            frame,
            text="Uruchom pełny pipeline",
            command=lambda: self.run_full_pipeline(
                assembly_mode.get(),
                species.get(),
                custom_species.get(),
                functional_diamond_db.get(),
                hydrolase_hmm_db.get(),
                hydrolase_diamond_db.get(),
                threads.get()
            )
        ).pack(anchor="e")

    def choose_augustus_model_file(self, custom_species):
        selected = filedialog.askopenfilename(
            initialdir=str(self.project_path("data")),
            title="Wybierz plik własnego modelu Augustusa",
            filetypes=[
                ("Augustus/model files", "*.*"),
                ("Wszystkie pliki", "*.*")
            ]
        )
        if selected:
            custom_species.set(selected)
            self.write_log(f"Wybrano własny model Augustusa: {selected}")

    def build_diamond_database(self, target_variable, title):
        diamond_dir = self.project_path("data/databases/diamond")
        diamond_dir.mkdir(parents=True, exist_ok=True)

        if "hydrolaz" in title.lower():
            output_db = diamond_dir / "hydrolases.dmnd"
        else:
            output_db = diamond_dir / "uniprot_sprot.dmnd"

        target_variable.set(str(output_db))

        command = [
            sys.executable,
            "scripts/download_and_build_databases.py",
            "--kind",
            "diamond",
            "--output",
            str(output_db)
        ]
        self.run_command(title + " - pobieranie i budowanie", command)
    def prepare_hmm_database(self, target_variable):
        hmm_dir = self.project_path("data/databases/hmmer")
        hmm_dir.mkdir(parents=True, exist_ok=True)

        output_hmm = hmm_dir / "Pfam-A.hmm"
        target_variable.set(str(output_hmm))

        command = [
            sys.executable,
            "scripts/download_and_build_databases.py",
            "--kind",
            "hmm",
            "--output",
            str(output_hmm)
        ]
        self.run_command("Pobieranie i przygotowanie bazy HMM", command)
    def run_full_pipeline(
        self,
        assembly_mode,
        species,
        custom_species,
        functional_diamond_db,
        hydrolase_hmm_db,
        hydrolase_diamond_db,
        threads
    ):
        if species == "własny model":
            selected_species = custom_species.strip()
            if not selected_species:
                messagebox.showerror("Pełny pipeline", "Dla opcji własny model wskaż plik albo wpisz nazwę modelu Augustusa.")
                return

            possible_path = Path(selected_species)
            if possible_path.exists():
                selected_species = possible_path.stem
                self.write_log(f"Użyto nazwy modelu Augustusa na podstawie pliku: {selected_species}")
        else:
            selected_species = species

        if not functional_diamond_db.strip():
            messagebox.showerror("Pełny pipeline", "Wskaż albo zbuduj bazę DIAMOND do annotacji funkcjonalnej.")
            return
        if not hydrolase_hmm_db.strip():
            messagebox.showerror("Pełny pipeline", "Wskaż albo przygotuj bazę HMM do predykcji hydrolaz.")
            return
        if not hydrolase_diamond_db.strip():
            messagebox.showerror("Pełny pipeline", "Wskaż albo zbuduj bazę DIAMOND hydrolaz.")
            return

        command = [
            sys.executable,
            "scripts/run_full_pipeline.py",
            "--assembly-mode",
            assembly_mode,
            "--species",
            selected_species,
            "--functional-diamond-db",
            functional_diamond_db.strip(),
            "--hydrolase-hmm-db",
            hydrolase_hmm_db.strip(),
            "--hydrolase-diamond-db",
            hydrolase_diamond_db.strip(),
            "--threads",
            threads.strip() or "8"
        ]

        self.run_command("Pełny pipeline analizy genomu", command)
    def open_denovo_window(self):
        window = tk.Toplevel(self)
        window.title("Assemblacja de novo")
        window.geometry("860x560")
        window.minsize(760, 500)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Assemblacja de novo", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(frame, text="Wybierz typ assemblacji. Program uruchomi odpowiednie narzędzie w systemie Linux.", wraplength=760).pack(anchor="w", pady=(6, 18))

        self.add_assembly_option(
            frame,
            "Na bazie odczytów Illumina",
            "Narzędzie: SPAdes. Wejście: sparowane odczyty Illumina R1/R2 z folderu data/clean.",
            self.run_illumina_assembly
        )
        self.add_assembly_option(
            frame,
            "Na bazie odczytów ONT",
            "Narzędzie: Flye. Wejście: długie odczyty ONT z folderu data/clean.",
            self.run_ont_assembly
        )
        self.add_assembly_option(
            frame,
            "Assemblacja hybrydowa",
            "Narzędzie: SPAdes w trybie hybrydowym. Wejście: odczyty Illumina R1/R2 oraz ONT.",
            self.run_hybrid_assembly
        )

    def add_assembly_option(self, parent, title, description, command):
        box = ttk.Frame(parent, padding=12, relief="ridge")
        box.pack(fill="x", pady=7)
        ttk.Label(box, text=title, font=("Segoe UI", 12, "bold")).pack(anchor="w")
        ttk.Label(box, text=description, wraplength=740).pack(anchor="w", pady=(4, 8))
        ttk.Button(box, text="Uruchom", command=command).pack(anchor="e")

    def default_illumina_r1(self):
        return self.first_existing([
            "data/clean/Unknown_CK982-002R0001_1.clean.fq.gz",
            "data/clean/Unknown_CK982-002R0001_1.fq.gz",
            "data/raw/Unknown_CK982-002R0001_1.fq.gz"
        ])

    def default_illumina_r2(self):
        return self.first_existing([
            "data/clean/Unknown_CK982-002R0001_2.clean.fq.gz",
            "data/clean/Unknown_CK982-002R0001_2.fq.gz",
            "data/raw/Unknown_CK982-002R0001_2.fq.gz"
        ])

    def default_ont(self):
        return self.first_existing([
            "data/clean/ONT_CK982-001N0001_raw.clean.fq.gz",
            "data/clean/ONT_CK982-001N0001_raw.fq.gz",
            "data/raw/ONT_CK982-001N0001_raw.fq.gz"
        ])

    def run_illumina_assembly(self):
        command = [sys.executable, "scripts/run_denovo_assembly.py", "--mode", "illumina"]
        self.run_command("Pipeline assemblacji de novo Illumina", command)

    def run_ont_assembly(self):
        command = [sys.executable, "scripts/run_denovo_assembly.py", "--mode", "ont"]
        self.run_command("Pipeline assemblacji de novo ONT", command)

    def run_hybrid_assembly(self):
        command = [sys.executable, "scripts/run_denovo_assembly.py", "--mode", "hybrid"]
        self.run_command("Pipeline assemblacji hybrydowej Illumina + ONT", command)

    def run_command(self, title, command):
        if self.current_process is not None:
            messagebox.showwarning("Proces w toku", "Inna analiza jest już uruchomiona.")
            return

        self.write_log("")
        self.write_log(f"Start: {title}")
        self.write_log("Komenda:")
        self.write_log(" ".join(command))

        thread = threading.Thread(target=self._run_command_worker, args=(title, command), daemon=True)
        thread.start()

    def _run_command_worker(self, title, command):
        try:
            self.current_process = subprocess.Popen(
                command,
                cwd=self.project_dir.get(),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            for line in self.current_process.stdout:
                self.after(0, self.write_log, line.rstrip())

            return_code = self.current_process.wait()

            if return_code == 0:
                self.after(0, self.write_log, f"Zakończono powodzeniem: {title}")
                self.after(0, messagebox.showinfo, title, "Assemblacja zakończona powodzeniem.")
            else:
                self.after(0, self.write_log, f"Proces zakończył się błędem. Kod: {return_code}")
                self.after(0, messagebox.showerror, title, f"Proces zakończył się błędem. Kod: {return_code}")

        except FileNotFoundError:
            tool = command[0]
            self.after(0, self.write_log, f"Nie znaleziono narzędzia: {tool}")
            self.after(0, messagebox.showerror, title, f"Nie znaleziono narzędzia: {tool}. Sprawdź, czy jest zainstalowane w Linuxie i dostępne w PATH.")
        except Exception as error:
            self.after(0, self.write_log, f"Błąd uruchomienia: {error}")
            self.after(0, messagebox.showerror, title, f"Błąd uruchomienia: {error}")
        finally:
            self.current_process = None

    def write_log(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log.insert("end", f"[{timestamp}] {message}\n")
        self.log.see("end")
        self.status_text.set(message if message else "Praca programu")

    def clear_log(self):
        self.log.delete("1.0", "end")
        self.status_text.set("Log wyczyszczony")

    def placeholder(self, module_name):
        self.write_log(f"Wybrano moduł: {module_name}. Funkcja zostanie dodana później.")
        messagebox.showinfo(module_name, "Funkcja tego modułu zostanie dodana w kolejnym etapie.")

    def open_gene_prediction_window(self):
        window = tk.Toplevel(self)
        window.title("Predykcja genów")
        window.geometry("720x500")
        window.minsize(640, 440)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Predykcja genów", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(
            frame,
            text="Wybierz model grzyba dla narzędzia Augustus albo wskaż własny model.",
            wraplength=660
        ).pack(anchor="w", pady=(6, 16))

        species_list = [
            "aspergillus_nidulans",
            "botrytis_cinerea",
            "candida_albicans",
            "candida_guilliermondii",
            "candida_tropicalis",
            "cryptococcus_neoformans_gattii",
            "fusarium_graminearum",
            "laccaria_bicolor",
            "neurospora_crassa",
            "saccharomyces_cerevisiae_S288C",
            "ustilago_maydis",
            "yarrowia_lipolytica",
            "własny model"
        ]

        selected_species = tk.StringVar(value="aspergillus_nidulans")
        custom_species = tk.StringVar(value="")
        genome_path = tk.StringVar(value=str(self.project_path("data/assemble_genome/latest_assembly.fasta")))

        ttk.Label(frame, text="Referencyjny gatunek/model:").pack(anchor="w")
        combo = ttk.Combobox(frame, textvariable=selected_species, values=species_list, state="readonly")
        combo.pack(fill="x", pady=(6, 12))

        ttk.Label(frame, text="Własny model Augustusa:").pack(anchor="w")
        custom_row = ttk.Frame(frame)
        custom_row.pack(fill="x", pady=(6, 12))
        ttk.Entry(custom_row, textvariable=custom_species).pack(side="left", fill="x", expand=True)
        ttk.Button(
            custom_row,
            text="Wybierz",
            command=lambda: self.choose_augustus_model_file(custom_species)
        ).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Plik wejściowy FASTA ze złożonym genomem:").pack(anchor="w")
        genome_row = ttk.Frame(frame)
        genome_row.pack(fill="x", pady=(6, 14))

        ttk.Entry(genome_row, textvariable=genome_path).pack(side="left", fill="x", expand=True)
        ttk.Button(
            genome_row,
            text="Wybierz",
            command=lambda: self.choose_genome_file(genome_path)
        ).pack(side="left", padx=(8, 0))

        ttk.Button(
            frame,
            text="Uruchom predykcję genów",
            command=lambda: self.run_gene_prediction_with_species(
                selected_species.get(),
                custom_species.get(),
                genome_path.get()
            )
        ).pack(anchor="e")

    def choose_genome_file(self, genome_path):
        selected = filedialog.askopenfilename(
            initialdir=str(self.project_path("data/assemble_genome")),
            title="Wybierz plik FASTA ze złożonym genomem",
            filetypes=[
                ("FASTA", "*.fasta *.fa *.fna"),
                ("Wszystkie pliki", "*.*")
            ]
        )
        if selected:
            genome_path.set(selected)
            self.write_log(f"Wybrano plik genomu do predykcji genów: {selected}")

    def run_gene_prediction_with_species(self, selected_species, custom_species, genome_path):
        if selected_species == "własny model":
            species = custom_species.strip()
            if not species:
                messagebox.showerror("Predykcja genów", "Dla opcji własny model wskaż plik albo wpisz nazwę modelu Augustusa.")
                return

            possible_path = Path(species)
            if possible_path.exists():
                species = possible_path.stem
                self.write_log(f"Użyto nazwy modelu Augustusa na podstawie pliku: {species}")
        else:
            species = selected_species

        command = [
            sys.executable,
            "scripts/run_gene_prediction.py",
            "--genome",
            genome_path,
            "--species",
            species
        ]

        self.run_command(f"Predykcja genów - Augustus ({species})", command)
    def run_gene_prediction(self):
        self.open_gene_prediction_window()

    def open_functional_annotation_window(self):
        window = tk.Toplevel(self)
        window.title("Annotacja funkcjonalna")
        window.geometry("760x500")
        window.minsize(680, 440)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Annotacja funkcjonalna", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(
            frame,
            text="Wybierz narzędzie do annotacji funkcjonalnej przewidywanych białek.",
            wraplength=700
        ).pack(anchor="w", pady=(6, 16))

        proteins_path = tk.StringVar(value=str(self.project_path("data/predicted_genes/predicted_proteins.faa")))
        selected_tool = tk.StringVar(value="diamond")
        diamond_db = tk.StringVar(value="")

        ttk.Label(frame, text="Plik wejściowy FASTA z białkami:").pack(anchor="w")
        proteins_row = ttk.Frame(frame)
        proteins_row.pack(fill="x", pady=(6, 14))

        ttk.Entry(proteins_row, textvariable=proteins_path).pack(side="left", fill="x", expand=True)
        ttk.Button(
            proteins_row,
            text="Wybierz",
            command=lambda: self.choose_proteins_file(proteins_path)
        ).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Narzędzie:").pack(anchor="w")
        tool_combo = ttk.Combobox(
            frame,
            textvariable=selected_tool,
            values=["diamond", "eggnog", "interproscan"],
            state="readonly"
        )
        tool_combo.pack(fill="x", pady=(6, 14))

        ttk.Label(frame, text="Baza DIAMOND, wymagana tylko dla DIAMOND:").pack(anchor="w")
        db_row = ttk.Frame(frame)
        db_row.pack(fill="x", pady=(6, 14))

        ttk.Entry(db_row, textvariable=diamond_db).pack(side="left", fill="x", expand=True)
        ttk.Button(
            db_row,
            text="Wybierz",
            command=lambda: self.choose_diamond_database(diamond_db)
        ).pack(side="left", padx=(8, 0))

        button_row = ttk.Frame(frame)
        button_row.pack(fill="x", pady=(6, 0))

        ttk.Button(
            button_row,
            text="Uruchom wybrane narzędzie",
            command=lambda: self.run_functional_annotation_tool(
                selected_tool.get(),
                proteins_path.get(),
                diamond_db.get()
            )
        ).pack(side="right")

        ttk.Button(
            button_row,
            text="Uruchom wszystkie analizy",
            command=lambda: self.run_functional_annotation_tool(
                "all",
                proteins_path.get(),
                diamond_db.get(),
                False
            )
        ).pack(side="right", padx=(0, 8))

        ttk.Button(
            button_row,
            text="Pobierz bazę DIAMOND i uruchom wszystkie",
            command=lambda: self.run_functional_annotation_tool(
                "all",
                proteins_path.get(),
                diamond_db.get(),
                True
            )
        ).pack(side="right", padx=(0, 8))

    def choose_proteins_file(self, proteins_path):
        selected = filedialog.askopenfilename(
            initialdir=str(self.project_path("data/predicted_genes")),
            title="Wybierz plik FASTA z białkami",
            filetypes=[
                ("FASTA protein", "*.faa *.fasta *.fa"),
                ("Wszystkie pliki", "*.*")
            ]
        )
        if selected:
            proteins_path.set(selected)
            self.write_log(f"Wybrano plik białek do annotacji: {selected}")

    def choose_diamond_database(self, diamond_db):
        selected = filedialog.askopenfilename(
            initialdir=str(self.project_path("data")),
            title="Wybierz bazę DIAMOND",
            filetypes=[
                ("DIAMOND database", "*.dmnd"),
                ("Wszystkie pliki", "*.*")
            ]
        )
        if selected:
            diamond_db.set(selected)
            self.write_log(f"Wybrano bazę DIAMOND: {selected}")

    def run_functional_annotation_tool(self, tool, proteins_path, diamond_db, download_diamond_db=False):
        command = [
            sys.executable,
            "scripts/run_functional_annotation.py",
            "--tool",
            tool,
            "--proteins",
            proteins_path
        ]

        if download_diamond_db:
            command.append("--download-diamond-db")
        elif tool in ["diamond", "all"]:
            if not diamond_db.strip():
                messagebox.showerror("DIAMOND", "Dla DIAMOND trzeba wskazać bazę .dmnd.")
                return
            command.extend(["--diamond-db", diamond_db.strip()])

        self.run_command(f"Annotacja funkcjonalna - {tool}", command)

    def run_annotation(self):
        self.open_functional_annotation_window()

    def open_hydrolase_prediction_window(self):
        window = tk.Toplevel(self)
        window.title("Predykcja hydrolaz")
        window.geometry("820x640")
        window.minsize(740, 560)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Predykcja hydrolaz", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(
            frame,
            text="Wybierz analizy, które mają zostać uruchomione dla przewidywanych białek.",
            wraplength=760
        ).pack(anchor="w", pady=(6, 16))

        proteins_path = tk.StringVar(value=str(self.project_path("data/predicted_genes/predicted_proteins.faa")))
        hmm_db = tk.StringVar(value="")
        diamond_db = tk.StringVar(value="")

        use_dbcan = tk.BooleanVar(value=True)
        use_hmmer = tk.BooleanVar(value=True)
        use_diamond = tk.BooleanVar(value=True)
        use_signalp = tk.BooleanVar(value=True)
        use_deeptmhmm = tk.BooleanVar(value=True)

        ttk.Label(frame, text="Plik wejściowy FASTA z białkami:").pack(anchor="w")
        proteins_row = ttk.Frame(frame)
        proteins_row.pack(fill="x", pady=(6, 12))

        ttk.Entry(proteins_row, textvariable=proteins_path).pack(side="left", fill="x", expand=True)
        ttk.Button(
            proteins_row,
            text="Wybierz",
            command=lambda: self.choose_proteins_file(proteins_path)
        ).pack(side="left", padx=(8, 0))

        options_box = ttk.Frame(frame, padding=12, relief="ridge")
        options_box.pack(fill="x", pady=(6, 14))

        ttk.Label(options_box, text="Analizy do uruchomienia:", font=("Segoe UI", 11, "bold")).pack(anchor="w")
        ttk.Checkbutton(options_box, text="dbCAN - predykcja CAZymes/hydrolaz węglowodanowych", variable=use_dbcan).pack(anchor="w", pady=(6, 0))
        ttk.Checkbutton(options_box, text="HMMER - wyszukiwanie domen w bazie HMM", variable=use_hmmer).pack(anchor="w")
        ttk.Checkbutton(options_box, text="DIAMOND - porównanie z bazą sekwencji hydrolaz", variable=use_diamond).pack(anchor="w")
        ttk.Checkbutton(options_box, text="SignalP - przewidywanie peptydów sygnałowych", variable=use_signalp).pack(anchor="w")
        ttk.Checkbutton(options_box, text="DeepTMHMM - przewidywanie domen transbłonowych", variable=use_deeptmhmm).pack(anchor="w")

        select_row = ttk.Frame(options_box)
        select_row.pack(fill="x", pady=(8, 0))
        ttk.Button(
            select_row,
            text="Zaznacz wszystkie",
            command=lambda: self.set_hydrolase_options(
                [use_dbcan, use_hmmer, use_diamond, use_signalp, use_deeptmhmm],
                True
            )
        ).pack(side="left")
        ttk.Button(
            select_row,
            text="Odznacz wszystkie",
            command=lambda: self.set_hydrolase_options(
                [use_dbcan, use_hmmer, use_diamond, use_signalp, use_deeptmhmm],
                False
            )
        ).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Baza HMM dla HMMER:").pack(anchor="w")
        hmm_row = ttk.Frame(frame)
        hmm_row.pack(fill="x", pady=(6, 12))

        ttk.Entry(hmm_row, textvariable=hmm_db).pack(side="left", fill="x", expand=True)
        ttk.Button(
            hmm_row,
            text="Wybierz",
            command=lambda: self.choose_hmm_database(hmm_db)
        ).pack(side="left", padx=(8, 0))
        ttk.Button(
            hmm_row,
            text="Przygotuj HMM",
            command=lambda: self.prepare_hmm_database(hmm_db)
        ).pack(side="left", padx=(8, 0))

        ttk.Label(frame, text="Baza DIAMOND hydrolaz:").pack(anchor="w")
        diamond_row = ttk.Frame(frame)
        diamond_row.pack(fill="x", pady=(6, 12))

        ttk.Entry(diamond_row, textvariable=diamond_db).pack(side="left", fill="x", expand=True)
        ttk.Button(
            diamond_row,
            text="Wybierz",
            command=lambda: self.choose_diamond_database(diamond_db)
        ).pack(side="left", padx=(8, 0))
        ttk.Button(
            diamond_row,
            text="Zbuduj",
            command=lambda: self.build_diamond_database(diamond_db, "Zbuduj bazę DIAMOND hydrolaz")
        ).pack(side="left", padx=(8, 0))

        ttk.Button(
            frame,
            text="Uruchom wybrane analizy hydrolaz",
            command=lambda: self.run_selected_hydrolase_predictions(
                proteins_path.get(),
                hmm_db.get(),
                diamond_db.get(),
                use_dbcan.get(),
                use_hmmer.get(),
                use_diamond.get(),
                use_signalp.get(),
                use_deeptmhmm.get()
            )
        ).pack(anchor="e", pady=(8, 0))

    def set_hydrolase_options(self, variables, value):
        for variable in variables:
            variable.set(value)

    def choose_hmm_database(self, hmm_db):
        selected = filedialog.askopenfilename(
            initialdir=str(self.project_path("data")),
            title="Wybierz bazę HMM",
            filetypes=[
                ("HMM database", "*.hmm"),
                ("Wszystkie pliki", "*.*")
            ]
        )
        if selected:
            hmm_db.set(selected)
            self.write_log(f"Wybrano bazę HMM: {selected}")

    def run_selected_hydrolase_predictions(
        self,
        proteins_path,
        hmm_db,
        diamond_db,
        use_dbcan,
        use_hmmer,
        use_diamond,
        use_signalp,
        use_deeptmhmm
    ):
        selected_tools = []
        if use_dbcan:
            selected_tools.append("dbcan")
        if use_hmmer:
            selected_tools.append("hmmer")
        if use_diamond:
            selected_tools.append("diamond")
        if use_signalp:
            selected_tools.append("signalp")
        if use_deeptmhmm:
            selected_tools.append("deeptmhmm")

        if not selected_tools:
            messagebox.showerror("Predykcja hydrolaz", "Wybierz przynajmniej jedną analizę.")
            return

        if "hmmer" in selected_tools and not hmm_db.strip():
            messagebox.showerror("HMMER", "Dla HMMER trzeba wskazać bazę HMM.")
            return

        if "diamond" in selected_tools and not diamond_db.strip():
            messagebox.showerror("DIAMOND", "Dla DIAMOND trzeba wskazać bazę .dmnd.")
            return

        for tool in selected_tools:
            command = [
                sys.executable,
                "scripts/run_hydrolase_prediction.py",
                "--tool",
                tool,
                "--proteins",
                proteins_path
            ]

            if tool == "hmmer":
                command.extend(["--hmm-db", hmm_db.strip()])
            if tool == "diamond":
                command.extend(["--diamond-db", diamond_db.strip()])

            self.run_command(f"Predykcja hydrolaz - {tool}", command)
    def run_hydrolases(self):
        self.open_hydrolase_prediction_window()


def main():
    app = GenomePipelineApp()
    app.mainloop()


if __name__ == "__main__":
    main()


















