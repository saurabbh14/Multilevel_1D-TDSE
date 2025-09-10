"""
PySide6-based GUI for ML-1D-TDSE

Features:
- Edit input.ini (raw editor)
- Quick fields (common params)
- Build (nix build) and Run/Stop simulation (captures stdout/stderr)
- Run log viewer
- Simple plot of output_data/density_1d.out using matplotlib

Install:
  python3 -m pip install -r tools/requirements_pyqt.txt

Run:
  python3 tools/pyqt_ui.py
"""
from pathlib import Path
import sys
import re
import subprocess
import threading

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QFileDialog, QMessageBox, QPlainTextEdit, QLineEdit, QLabel, QFormLayout,
    QSplitter, QTextEdit, QSizePolicy
)
from PySide6.QtCore import Qt, Slot, Signal, QObject

import numpy as np
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "input.ini"
EXEC_PATH = PROJECT_ROOT / "result" / "bin" / "ML-TDSE"


# --- simple namelist helpers (same basic heuristics as the tkinter version) ---
def find_namelist_block(text, name):
    pat = re.compile(r"(^\s*&" + re.escape(name) + r"\b.*?^\s*/\s*$)", re.M | re.S)
    m = pat.search(text)
    return (m.start(1), m.end(1)) if m else (None, None)


def set_or_add_key(text, section, key, value):
    s_start, s_end = find_namelist_block(text, section)
    key_line = f"{key}={value}"
    if s_start is None:
        block = f"\n&{section}\n  {key_line}\n/\n"
        return text + block
    block = text[s_start:s_end]
    key_pat = re.compile(r"(^\s*" + re.escape(key) + r"\s*=.*?$)", re.M)
    if key_pat.search(block):
        block2 = key_pat.sub(key_line, block)
    else:
        block2 = block.rstrip()
        block2 = block2[:-1] + f"  {key_line}\n/\n"
    return text[:s_start] + block2 + text[s_end:]


# --- Worker to run subprocess and emit lines ---
class ProcRunner(QObject):
    line = Signal(str)
    finished = Signal(int)

    def __init__(self, cmd, cwd):
        super().__init__()
        self.cmd = cmd
        self.cwd = cwd
        self._proc = None

    def run(self):
        try:
            self._proc = subprocess.Popen(self.cmd, cwd=self.cwd, stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT, text=True)
        except FileNotFoundError as e:
            self.line.emit(f"Error: {e}\n")
            self.finished.emit(-1)
            return
        for ln in self._proc.stdout:
            self.line.emit(ln)
        rc = self._proc.wait()
        self.finished.emit(rc)

    def terminate(self):
        if self._proc:
            self._proc.terminate()


# --- Main Window ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ML-1D-TDSE — PySide6 UI")
        self.resize(1200, 800)
        self.current_input = DEFAULT_INPUT if DEFAULT_INPUT.exists() else None
        self.proc_runner = None
        self.proc_thread = None

        self._build_ui()
        if self.current_input:
            self.load_input(self.current_input)

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # Top buttons
        btn_layout = QHBoxLayout()
        for txt, cb in [
            ("Open input.ini", self.open_input),
            ("Save input.ini", self.save_input),
            ("Build (nix build)", self.build_project),
            ("Run", self.run_simulation),
            ("Stop", self.stop_simulation),
            ("Plot density_1d", self.plot_density),
        ]:
            b = QPushButton(txt)
            b.clicked.connect(cb)
            btn_layout.addWidget(b)
        btn_layout.addStretch()
        main_layout.addLayout(btn_layout)

        # Split editor / quick fields + log
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter, stretch=3)

        # Left: editor
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        self.editor = QPlainTextEdit()
        left_layout.addWidget(QLabel("input.ini (raw editor)"))
        left_layout.addWidget(self.editor)
        splitter.addWidget(left_widget)

        # Right: quick fields and log
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        # quick fields
        qform = QFormLayout()
        self.qfields = {}
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
            le = QLineEdit()
            qform.addRow(QLabel(label), le)
            self.qfields[(section, key)] = le
        right_layout.addLayout(qform)
        apply_btn = QPushButton("Apply quick fields → editor")
        apply_btn.clicked.connect(self.apply_quick_fields)
        right_layout.addWidget(apply_btn)

        # log
        right_layout.addWidget(QLabel("Run log"))
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setStyleSheet("background:#111;color:#ddd;")
        right_layout.addWidget(self.log, stretch=1)

        splitter.addWidget(right_widget)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Plot area
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        main_layout.addWidget(plot_widget, stretch=2)
        self.fig = Figure(figsize=(6, 3))
        self.canvas = FigureCanvas(self.fig)
        plot_layout.addWidget(self.canvas)

    # --- file ops ---
    def open_input(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Open input.ini", str(PROJECT_ROOT), "INI Files (*.ini);;All Files (*)")
        if fn:
            self.load_input(Path(fn))

    def load_input(self, path: Path):
        try:
            txt = path.read_text(encoding="utf-8")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load: {e}")
            return
        self.current_input = path
        self.editor.setPlainText(txt)
        self.refresh_quick_fields()
        self.append_log(f"Loaded: {path}\n")

    def save_input(self):
        if not self.current_input:
            fn, _ = QFileDialog.getSaveFileName(self, "Save input.ini", str(PROJECT_ROOT / "input.ini"))
            if not fn:
                return
            self.current_input = Path(fn)
        try:
            self.current_input.write_text(self.editor.toPlainText(), encoding="utf-8")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to save: {e}")
            return
        self.append_log(f"Saved: {self.current_input}\n")
        QMessageBox.information(self, "Saved", str(self.current_input))

    # --- quick fields ---
    def refresh_quick_fields(self):
        text = self.editor.toPlainText()
        for (section, key), widget in self.qfields.items():
            s_start, s_end = find_namelist_block(text, section)
            val = ""
            if s_start is not None:
                block = text[s_start:s_end]
                m = re.search(r"\b" + re.escape(key) + r"\s*=\s*([^!\n\r]+)", block)
                if m:
                    val = m.group(1).strip().strip('"').strip("'")
            widget.setText(val)

    def apply_quick_fields(self):
        text = self.editor.toPlainText()
        for (section, key), widget in self.qfields.items():
            v = widget.text().strip()
            if not v:
                continue
            if not re.match(r"^[0-9eE\.\+\-]*$", v):
                v2 = '"' + v + '"'
            else:
                v2 = v
            text = set_or_add_key(text, section, key, v2)
        self.editor.setPlainText(text)
        self.append_log("Applied quick fields to editor.\n")

    # --- build / run ---
    def build_project(self):
        self.append_log("Starting build: nix build\n")
        runner = ProcRunner(["nix", "build"], str(PROJECT_ROOT))
        runner.line.connect(self.append_log)
        runner.finished.connect(lambda rc: self.append_log(f"Build finished rc={rc}\n"))
        self._start_runner(runner)

    def run_simulation(self):
        if not EXEC_PATH.exists():
            QMessageBox.warning(self, "Not found", f"Executable not found: {EXEC_PATH}\nRun build first.")
            self.append_log(f"Executable not found: {EXEC_PATH}\n")
            return
        self.save_input()
        cmd = [str(EXEC_PATH), str(self.current_input)]
        self.append_log(f"Starting simulation: {' '.join(cmd)}\n")
        runner = ProcRunner(cmd, str(PROJECT_ROOT))
        runner.line.connect(self.append_log)
        runner.finished.connect(lambda rc: self.append_log(f"Simulation finished (rc={rc})\n"))
        self._start_runner(runner)

    def stop_simulation(self):
        if self.proc_runner:
            self.proc_runner.terminate()
            self.append_log("Requested termination of running process.\n")
        else:
            self.append_log("No running process.\n")

    def _start_runner(self, runner: ProcRunner):
        # stop existing
        if self.proc_runner:
            self.append_log("A process is already running. Stop it first.\n")
            return
        self.proc_runner = runner
        th = threading.Thread(target=runner.run, daemon=True)
        self.proc_thread = th
        th.start()
        # poll thread finish to clear proc_runner
        def poll():
            th.join()
            self.proc_runner = None
        threading.Thread(target=poll, daemon=True).start()

    # --- plotting ---
    def plot_density(self):
        text = self.editor.toPlainText()
        m = re.search(r"\boutput_data_dir\s*=\s*['\"]?([^'\"]+?)['\"\s\r\n/]*\n", text)
        outdir = PROJECT_ROOT
        if m:
            outdir = (PROJECT_ROOT / m.group(1).strip()).resolve()
        f = outdir / "density_1d.out"
        if not f.exists():
            QMessageBox.information(self, "Not found", f"{f} not found. Ensure simulation produced outputs.")
            self.append_log(f"Plot file not found: {f}\n")
            return
        try:
            data = np.loadtxt(f)
        except Exception as e:
            self.append_log(f"Error reading {f}: {e}\n")
            return
        ax = self.fig.subplots()
        ax.clear()
        if data.ndim == 1:
            ax.plot(data)
        else:
            x = data[:, 0]
            y = data[:, 1] if data.shape[1] > 1 else data[:, 0]
            ax.plot(x, y)
        ax.set_title(f.name)
        self.canvas.draw()
        self.append_log(f"Plotted: {f}\n")

    # --- logging ---
    @Slot(str)
    def append_log(self, txt):
        self.log.moveCursor(self.log.textCursor().End)
        self.log.insertPlainText(txt)
        self.log.moveCursor(self.log.textCursor().End)


def main():
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()