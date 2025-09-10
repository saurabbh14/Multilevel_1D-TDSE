"""
Simple Tkinter-based UI for the Multi-level 1D TDSE project.

Features:
- Open / edit / save input.ini (raw text editor + quick fields for common options)
- Build project (nix build) and run simulation (./result/bin/ML-TDSE input.ini)
- Capture and display stdout/stderr in a log pane
- Simple plotting of common output file (e.g. output_data/density_1d.out)
- Start / Stop running simulation (uses subprocess.Popen)

Usage:
$ python3 tools/gui.py

Dependencies:
- Python 3.8+
- tkinter (standard on most Linux distros)
- matplotlib, numpy

Install Python deps:
$ python3 -m pip install matplotlib numpy
"""
import os
import re
import threading
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

PROJECT_ROOT = Path(__file__).resolve().parents[1]  # .../TDSE_1D
DEFAULT_INPUT = PROJECT_ROOT / "input.ini"
EXEC_PATH = PROJECT_ROOT / "result" / "bin" / "ML-TDSE"

# --- helpers for namelist editing (basic, robust enough for simple edits) ---


def find_namelist_block(text, name):
    """Return (start_idx, end_idx) of &name ... / block or (None, None)."""
    pat = re.compile(r"(^\s*&" + re.escape(name) + r"\b.*?^\s*/\s*$)", re.M | re.S)
    m = pat.search(text)
    return (m.start(1), m.end(1)) if m else (None, None)


def set_or_add_key(text, section, key, value):
    """
    Set key=value in &section ... / block. If block or key missing, add it.
    Returns modified text.
    """
    s_start, s_end = find_namelist_block(text, section)
    key_line = f"{key}={value}"
    if s_start is None:
        # add a new namelist at end
        block = f"\n&{section}\n  {key_line}\n/\n"
        return text + block
    block = text[s_start:s_end]
    # if key present replace whole line
    key_pat = re.compile(r"(^\s*" + re.escape(key) + r"\s*=.*?$)", re.M)
    if key_pat.search(block):
        block2 = key_pat.sub(key_line, block)
    else:
        # insert before trailing / (i.e., before end)
        block2 = block.rstrip()
        block2 = block2[:-1] + f"  {key_line}\n/\n"
    return text[:s_start] + block2 + text[s_end:]


def read_text_file(path):
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def write_text_file(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


# --- GUI App ---


class TDSEApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ML-1D-TDSE — UI")
        self.geometry("1100x700")
        self._proc = None
        self._proc_lock = threading.Lock()
        self._build_thread = None

        self._create_widgets()
        self.load_input(DEFAULT_INPUT if DEFAULT_INPUT.exists() else None)

    def _create_widgets(self):
        # Top frame: buttons
        frm_top = ttk.Frame(self)
        frm_top.pack(fill="x", padx=6, pady=4)

        ttk.Button(frm_top, text="Open input.ini", command=self.open_input).pack(side="left")
        ttk.Button(frm_top, text="Save input.ini", command=self.save_input).pack(side="left")
        ttk.Button(frm_top, text="Build (nix build)", command=self.build_project).pack(side="left")
        ttk.Button(frm_top, text="Run", command=self.run_simulation).pack(side="left")
        ttk.Button(frm_top, text="Stop", command=self.stop_simulation).pack(side="left")
        ttk.Button(frm_top, text="Plot density_1d", command=self.plot_density).pack(side="left")
        ttk.Button(frm_top, text="Refresh quick fields", command=self.refresh_quick_fields).pack(side="left")

        # Left: text editor for input.ini and quick fields
        pan = ttk.PanedWindow(self, orient="horizontal")
        pan.pack(fill="both", expand=True, padx=6, pady=4)

        left = ttk.Frame(pan)
        pan.add(left, weight=3)
        right = ttk.Frame(pan)
        pan.add(right, weight=1)

        # Text editor
        lbl = ttk.Label(left, text="input.ini (raw editor)")
        lbl.pack(anchor="w")
        self.txt = tk.Text(left, wrap="none", undo=True)
        self.txt.pack(fill="both", expand=True)
        txt_scroll_y = ttk.Scrollbar(left, orient="vertical", command=self.txt.yview)
        txt_scroll_y.pack(side="right", fill="y")
        self.txt.configure(yscrollcommand=txt_scroll_y.set)

        # Quick fields
        frm_q = ttk.Labelframe(right, text="Quick fields")
        frm_q.pack(fill="both", expand=False, padx=4, pady=4)

        # A few commonly changed fields
        self.q_vars = {}
        for label, section, key in [
            ("NR", "grid", "NR"),
            ("dt (fs)", "time_grid", "dt"),
            ("Nt", "time_grid", "Nt"),
            ("Nstates", "elec_states", "Nstates"),
            ("initial_distribution", "ini_state", "initial_distribution"),
            ("v_ini", "ini_state", "v_ini"),
            ("absorber", "absorber_choice", "absorber"),
            ("output_data_dir", "output_files", "output_data_dir"),
        ]:
            row = ttk.Frame(frm_q)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=label, width=18).pack(side="left")
            v = tk.StringVar()
            ent = ttk.Entry(row, textvariable=v, width=28)
            ent.pack(side="left", padx=2)
            self.q_vars[(section, key)] = v

        ttk.Button(frm_q, text="Apply quick fields → editor", command=self.apply_quick_fields).pack(pady=6)

        # Log pane
        frm_log = ttk.Labelframe(right, text="Run log")
        frm_log.pack(fill="both", expand=True, padx=4, pady=4)
        self.log = tk.Text(frm_log, height=12, wrap="none", bg="#111", fg="#ddd")
        self.log.pack(fill="both", expand=True)

        # Plot area (bottom)
        frm_plot = ttk.Labelframe(self, text="Plot")
        frm_plot.pack(fill="both", expand=True, padx=6, pady=4)
        self.fig = Figure(figsize=(6, 3))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm_plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # ---- file operations ----
    def open_input(self):
        fn = filedialog.askopenfilename(title="Open input.ini", initialdir=PROJECT_ROOT, filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if fn:
            self.load_input(fn)

    def load_input(self, path):
        if path is None:
            self.txt.delete("1.0", tk.END)
            return
        text = read_text_file(path)
        self.current_input_path = Path(path)
        self.txt.delete("1.0", tk.END)
        self.txt.insert("1.0", text)
        self.refresh_quick_fields()
        self.append_log(f"Loaded: {path}\n")

    def save_input(self):
        path = getattr(self, "current_input_path", DEFAULT_INPUT)
        write_text_file(path, self.txt.get("1.0", tk.END))
        self.append_log(f"Saved input: {path}\n")
        messagebox.showinfo("Saved", f"Saved: {path}")

    # ---- quick fields ----
    def refresh_quick_fields(self):
        text = self.txt.get("1.0", tk.END)
        for (section, key), var in self.q_vars.items():
            # regex: find key=VALUE inside &section ... /
            s_start, s_end = find_namelist_block(text, section)
            val = ""
            if s_start is not None:
                block = text[s_start:s_end]
                m = re.search(r"\b" + re.escape(key) + r"\s*=\s*([^!\n\r]+)", block)
                if m:
                    val = m.group(1).strip().strip('"').strip("'")
            var.set(val)

    def apply_quick_fields(self):
        text = self.txt.get("1.0", tk.END)
        for (section, key), var in self.q_vars.items():
            v = var.get().strip()
            if v == "":
                continue
            # ensure strings are quoted in ini
            if not re.match(r"^[0-9eE\.\+\-]*$", v):
                v2 = '"' + v + '"'
            else:
                v2 = v
            text = set_or_add_key(text, section, key, v2)
        self.txt.delete("1.0", tk.END)
        self.txt.insert("1.0", text)
        self.append_log("Applied quick fields to editor.\n")

    # ---- build / run ----
    def build_project(self):
        def run_build():
            self.append_log("Starting build: nix build\n")
            try:
                p = subprocess.Popen(["nix", "build"], cwd=PROJECT_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            except FileNotFoundError:
                self.append_log("Error: nix not found on PATH\n")
                return
            for line in p.stdout:
                self.append_log(line)
            p.wait()
            self.append_log(f"Build finished. returncode={p.returncode}\n")
        t = threading.Thread(target=run_build, daemon=True)
        t.start()

    def run_simulation(self):
        if not EXEC_PATH.exists():
            self.append_log(f"Executable not found: {EXEC_PATH}\nRun build first.\n")
            messagebox.showwarning("Not found", f"Executable not found: {EXEC_PATH}\nRun nix build first.")
            return
        # ensure input saved
        self.save_input()

        cmd = [str(EXEC_PATH), str(self.current_input_path)]
        self.append_log(f"Starting simulation: {' '.join(cmd)}\n")
        def target():
            with self._proc_lock:
                self._proc = subprocess.Popen(cmd, cwd=PROJECT_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            try:
                for line in self._proc.stdout:
                    self.append_log(line)
            finally:
                with self._proc_lock:
                    rc = self._proc.poll() if self._proc else None
                    self.append_log(f"Simulation finished (returncode={rc}).\n")
                    self._proc = None
        t = threading.Thread(target=target, daemon=True)
        t.start()

    def stop_simulation(self):
        with self._proc_lock:
            if self._proc:
                self._proc.terminate()
                self.append_log("Terminated running simulation.\n")
            else:
                self.append_log("No running simulation to stop.\n")

    # ---- plotting ----
    def plot_density(self):
        # locate output_data_dir in the input text
        text = self.txt.get("1.0", tk.END)
        m = re.search(r"\boutput_data_dir\s*=\s*['\"]?([^'\"]+?)['\"\s\r\n/]*\n", text)
        outdir = PROJECT_ROOT
        if m:
            outdir = (PROJECT_ROOT / m.group(1).strip()).resolve()
        # common file name
        f = outdir / "density_1d.out"
        if not f.exists():
            self.append_log(f"Plot file not found: {f}\n")
            messagebox.showinfo("Not found", f"{f} not found. Ensure simulation produced outputs.")
            return
        try:
            data = np.loadtxt(f)
        except Exception as e:
            self.append_log(f"Error reading {f}: {e}\n")
            return
        if data.ndim == 1:
            x = np.arange(data.size)
            y = data
        else:
            x = data[:, 0]
            y = data[:, 1] if data.shape[1] > 1 else data[:, 0]
        self.ax.clear()
        self.ax.plot(x, y, lw=1.0)
        self.ax.set_title(str(f.name))
        self.ax.set_xlabel("R (units in file)")
        self.ax.set_ylabel("density")
        self.canvas.draw()
        self.append_log(f"Plotted: {f}\n")

    # ---- logging ----
    def append_log(self, txt):
        self.log.insert(tk.END, txt)
        self.log.see(tk.END)


if __name__ == "__main__":
    app = TDSEApp()
    app.mainloop()