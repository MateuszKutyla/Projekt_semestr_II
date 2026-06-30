#!/usr/bin/env python3
import subprocess
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

        tk.Label(header, text="Analiza genomu grzyba", bg="#1f2937", fg="white", font=("Segoe UI", 20, "bold")).pack(anchor="w")
        tk.Label(header, text="Illumina + ONT | rdzeń programu do kolejnych etapów analizy", bg="#1f2937", fg="#d1d5db", font=("Segoe UI", 10)).pack(anchor="w", pady=(4, 0))

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

        self.add_module(button_grid, 0, 0, "Assemblacja de novo", "Wybór typu assemblacji: Illumina, ONT albo hybrydowa.", self.open_denovo_window)
        self.add_module(button_grid, 0, 1, "Predykcja genów", "Moduł do wykrywania genów w złożonym genomie.", self.run_gene_prediction)
        self.add_module(button_grid, 1, 0, "Annotacja funkcjonalna", "Moduł do przypisywania funkcji przewidywanym genom i białkom.", self.run_annotation)
        self.add_module(button_grid, 1, 1, "Predykcja hydrolaz", "Moduł do wyszukiwania potencjalnych enzymów hydrolitycznych.", self.run_hydrolases)

        ttk.Label(main, text="Log programu").pack(anchor="w", pady=(18, 6))
        self.log = tk.Text(main, height=12, bg="#111827", fg="#e5e7eb", insertbackground="white", wrap="word")
        self.log.pack(fill="both", expand=True)

        bottom = ttk.Frame(self, padding=(12, 8))
        bottom.pack(fill="x")
        ttk.Label(bottom, textvariable=self.status_text).pack(side="left")
        ttk.Button(bottom, text="Wyczyść log", command=self.clear_log).pack(side="right")

        self.write_log("Uruchomiono graficzny rdzeń projektu.")

    def add_module(self, parent, row, column, title, description, command):
        frame = ttk.Frame(parent, padding=14, relief="ridge")
        frame.grid(row=row, column=column, sticky="nsew", padx=8, pady=8)
        ttk.Label(frame, text=title, font=("Segoe UI", 15, "bold")).pack(anchor="w")
        ttk.Label(frame, text=description, wraplength=390).pack(anchor="w", pady=(6, 14))
        ttk.Button(frame, text=title, command=command).pack(fill="x")

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
        output_dir = self.project_path("results/assembly/spades_illumina")
        command = [
            "spades.py",
            "-1", self.default_illumina_r1(),
            "-2", self.default_illumina_r2(),
            "-o", str(output_dir),
            "-t", "8"
        ]
        self.run_command("Assemblacja de novo Illumina - SPAdes", command)

    def run_ont_assembly(self):
        output_dir = self.project_path("results/assembly/flye_ont")
        command = [
            "flye",
            "--nano-raw", self.default_ont(),
            "--out-dir", str(output_dir),
            "--threads", "8"
        ]
        self.run_command("Assemblacja de novo ONT - Flye", command)

    def run_hybrid_assembly(self):
        output_dir = self.project_path("results/assembly/spades_hybrid")
        command = [
            "spades.py",
            "-1", self.default_illumina_r1(),
            "-2", self.default_illumina_r2(),
            "--nanopore", self.default_ont(),
            "-o", str(output_dir),
            "-t", "8"
        ]
        self.run_command("Assemblacja hybrydowa Illumina + ONT - SPAdes", command)

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

    def run_gene_prediction(self):
        self.placeholder("Predykcja genów")

    def run_annotation(self):
        self.placeholder("Annotacja funkcjonalna")

    def run_hydrolases(self):
        self.placeholder("Predykcja hydrolaz")


def main():
    app = GenomePipelineApp()
    app.mainloop()


if __name__ == "__main__":
    main()
