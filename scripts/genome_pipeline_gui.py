#!/usr/bin/env python3
import tkinter as tk
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

PROJECT_ROOT = Path(__file__).resolve().parents[1]


class GenomePipelineApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Projekt semestralny II - analiza genomu grzyba")
        self.geometry("980x640")
        self.minsize(860, 560)

        self.project_dir = tk.StringVar(value=str(PROJECT_ROOT))
        self.status_text = tk.StringVar(value="Gotowe do pracy")

        self.configure(bg="#f4f6f8")
        self.build_ui()

    def build_ui(self):
        header = tk.Frame(self, bg="#1f2937", padx=22, pady=18)
        header.pack(fill="x")

        tk.Label(
            header,
            text="Analiza genomu grzyba",
            bg="#1f2937",
            fg="white",
            font=("Segoe UI", 20, "bold")
        ).pack(anchor="w")

        tk.Label(
            header,
            text="Illumina + ONT | rdzeń programu do kolejnych etapów analizy",
            bg="#1f2937",
            fg="#d1d5db",
            font=("Segoe UI", 10)
        ).pack(anchor="w", pady=(4, 0))

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

        self.log = tk.Text(main, height=10, bg="#111827", fg="#e5e7eb", insertbackground="white", wrap="word")
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

    def open_denovo_window(self):
        window = tk.Toplevel(self)
        window.title("Assemblacja de novo")
        window.geometry("720x430")
        window.minsize(640, 380)
        window.transient(self)

        frame = ttk.Frame(window, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="Assemblacja de novo", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(
            frame,
            text="Wybierz typ assemblacji. Na tym etapie przyciski tworzą rdzeń logiki; właściwe komendy zostaną podłączone później.",
            wraplength=650
        ).pack(anchor="w", pady=(6, 18))

        self.add_assembly_option(
            frame,
            "Na bazie odczytów Illumina",
            "Narzędzie: SPAdes. Ten tryb będzie używał sparowanych odczytów Illumina R1/R2.",
            self.run_illumina_assembly
        )

        self.add_assembly_option(
            frame,
            "Na bazie odczytów ONT",
            "Narzędzie: Flye. Ten tryb będzie używał długich odczytów Oxford Nanopore.",
            self.run_ont_assembly
        )

        self.add_assembly_option(
            frame,
            "Assemblacja hybrydowa",
            "Tryb wykorzystujący dwa rodzaje odczytów: Illumina oraz ONT.",
            self.run_hybrid_assembly
        )

    def add_assembly_option(self, parent, title, description, command):
        box = ttk.Frame(parent, padding=12, relief="ridge")
        box.pack(fill="x", pady=6)

        ttk.Label(box, text=title, font=("Segoe UI", 12, "bold")).pack(anchor="w")
        ttk.Label(box, text=description, wraplength=620).pack(anchor="w", pady=(4, 8))
        ttk.Button(box, text="Wybierz", command=command).pack(anchor="e")

    def choose_project_dir(self):
        selected = filedialog.askdirectory(initialdir=self.project_dir.get(), title="Wybierz folder projektu")
        if selected:
            self.project_dir.set(selected)
            self.write_log(f"Ustawiono folder projektu: {selected}")

    def write_log(self, message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log.insert("end", f"[{timestamp}] {message}\n")
        self.log.see("end")
        self.status_text.set(message)

    def clear_log(self):
        self.log.delete("1.0", "end")
        self.status_text.set("Log wyczyszczony")

    def placeholder(self, module_name, tool_name=None):
        if tool_name:
            self.write_log(f"Wybrano: {module_name}. Planowane narzędzie: {tool_name}.")
            messagebox.showinfo(module_name, f"Wybrano moduł: {module_name}\nPlanowane narzędzie: {tool_name}\nFunkcja zostanie podłączona w kolejnym etapie.")
        else:
            self.write_log(f"Wybrano moduł: {module_name}. Funkcja zostanie dodana później.")
            messagebox.showinfo(module_name, "Rdzeń interfejsu działa. Funkcja tego modułu zostanie dodana w kolejnym etapie.")

    def run_illumina_assembly(self):
        self.placeholder("Assemblacja de novo na bazie odczytów Illumina", "SPAdes")

    def run_ont_assembly(self):
        self.placeholder("Assemblacja de novo na bazie odczytów ONT", "Flye")

    def run_hybrid_assembly(self):
        self.placeholder("Assemblacja hybrydowa Illumina + ONT", "SPAdes/Flye oraz integracja odczytów Illumina i ONT")

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
