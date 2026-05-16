"""
Chandler Wobble Decomposer — IERS polar motion analysis.

Decomposes Earth Orientation Parameter time-series (x_pole, y_pole)
into:
    - Chandler oscillation (~433-day free Eulerian wobble)
    - Annual oscillation (~365.25-day forced wobble)
    - Linear and quadratic secular drifts

Methods: least-squares sinusoid fit plus discrete wavelet transform
(Morlet) for time-frequency analysis. Reproduces classical Chandler /
annual amplitudes (~150-200 mas / ~90-100 mas) seen in IERS EOP-C04.

Author: Dr. Mosab Hawarey (@DrHawarey)
License: MIT
"""

from __future__ import annotations

import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


CHANDLER_PERIOD_DAYS = 433.0
ANNUAL_PERIOD_DAYS = 365.25


def synthesize_polar_motion(n_years: float = 25.0, seed: int = 11) -> np.ndarray:
    """
    Generate a realistic polar-motion time-series mimicking IERS EOP-C04.

    Returns columns: [MJD, x_pole_arcsec, y_pole_arcsec].

    The IERS x/y_pole are reported in arcseconds with typical Chandler
    amplitude ~0.18", annual amplitude ~0.10", plus secular drift of
    roughly +3-4 mas/yr toward Greenland following GRACE-era ice loss.
    """
    rng = np.random.default_rng(seed)
    n = int(n_years * 365)
    mjd = 55000.0 + np.arange(n, dtype=float)
    t = (mjd - mjd[0]) / 365.25

    # Chandler (with slow phase wander to reflect excitation)
    phi_c = 2 * math.pi * (mjd - mjd[0]) / CHANDLER_PERIOD_DAYS
    a_c = 0.18 + 0.03 * np.cos(2 * math.pi * t / 6.4)   # 6.4-yr beat
    x_chan = a_c * np.cos(phi_c)
    y_chan = a_c * np.sin(phi_c)

    # Annual
    phi_a = 2 * math.pi * (mjd - mjd[0]) / ANNUAL_PERIOD_DAYS
    a_a = 0.10
    x_ann = a_a * np.cos(phi_a + 0.6)
    y_ann = a_a * np.sin(phi_a + 0.6)

    # Secular drift
    x_sec = +0.0033 * t
    y_sec = -0.0040 * t

    # Long-period (Markowitz) wobble + noise
    x_lp = 0.012 * np.cos(2 * math.pi * t / 25.0)
    y_lp = 0.010 * np.sin(2 * math.pi * t / 25.0)
    noise_amp = 0.0015
    x = x_chan + x_ann + x_sec + x_lp + rng.normal(0, noise_amp, n)
    y = y_chan + y_ann + y_sec + y_lp + rng.normal(0, noise_amp, n)
    return np.column_stack([mjd, x, y])


def least_squares_decompose(mjd: np.ndarray, signal: np.ndarray) -> dict:
    """
    Solve for amplitudes and phases of Chandler + annual + secular components.
    """
    t = (mjd - mjd[0])
    omega_c = 2 * math.pi / CHANDLER_PERIOD_DAYS
    omega_a = 2 * math.pi / ANNUAL_PERIOD_DAYS

    A = np.column_stack([
        np.cos(omega_c * t), np.sin(omega_c * t),
        np.cos(omega_a * t), np.sin(omega_a * t),
        np.ones_like(t), t / 365.25, (t / 365.25) ** 2,
    ])
    coef, *_ = np.linalg.lstsq(A, signal, rcond=None)
    chan_amp = math.hypot(coef[0], coef[1])
    ann_amp = math.hypot(coef[2], coef[3])
    chan_phase = math.degrees(math.atan2(-coef[1], coef[0]))
    ann_phase = math.degrees(math.atan2(-coef[3], coef[2]))
    return {
        "chandler_amp": chan_amp,
        "chandler_phase_deg": chan_phase,
        "annual_amp": ann_amp,
        "annual_phase_deg": ann_phase,
        "secular_mean": coef[4],
        "secular_slope_per_yr": coef[5],
        "secular_quadratic": coef[6],
        "model": A @ coef,
        "residual": signal - A @ coef,
    }


def morlet_wavelet_scalogram(signal: np.ndarray, dt_days: float,
                              periods_days: np.ndarray, w0: float = 6.0) -> np.ndarray:
    """
    Continuous wavelet transform with the Morlet basis.

    Returns |W|^2 with shape (n_periods, n_time).
    """
    n = len(signal)
    n_freqs = len(periods_days)
    sig_fft = np.fft.fft(signal - signal.mean())
    freqs = np.fft.fftfreq(n, d=dt_days)
    scalogram = np.zeros((n_freqs, n))

    for i, P in enumerate(periods_days):
        s = P / (1.0327 * 2)  # scale that matches Morlet center freq to period
        morlet_fft = (np.pi ** -0.25) * np.sqrt(2 * s) * \
                     np.exp(-(s * 2 * math.pi * freqs - w0) ** 2 / 2.0)
        morlet_fft[freqs <= 0] = 0
        wave = np.fft.ifft(sig_fft * morlet_fft)
        scalogram[i] = np.abs(wave) ** 2
    return scalogram


class WobbleApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Chandler Wobble Decomposer")
        self.root.geometry("1200x800")
        self.data: np.ndarray | None = None
        self._build_ui()
        self._load_synthetic()

    def _build_ui(self) -> None:
        left = ttk.Frame(self.root, padding=10)
        left.pack(side=tk.LEFT, fill=tk.Y)
        ttk.Label(left, text="Chandler Wobble", font=("Segoe UI", 14, "bold")).pack(anchor="w")
        ttk.Label(left, text="Decomposer", font=("Segoe UI", 14, "bold")).pack(anchor="w")
        ttk.Button(left, text="Load synthetic EOP series",
                   command=self._load_synthetic).pack(fill=tk.X, pady=8)
        ttk.Button(left, text="Load EOP-C04 CSV…",
                   command=self._load_csv).pack(fill=tk.X)
        ttk.Button(left, text="Decompose & plot",
                   command=self._analyze).pack(fill=tk.X, pady=(12, 4))

        self.out = tk.Text(left, width=38, height=28, wrap="word",
                           font=("Consolas", 9))
        self.out.pack(fill=tk.BOTH, pady=8)

        ttk.Label(left, text="Dr. Mosab Hawarey • MIT License",
                  font=("Segoe UI", 8), foreground="#777").pack(anchor="w")

        right = ttk.Frame(self.root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.fig = plt.Figure(figsize=(9, 8), dpi=100)
        gs = self.fig.add_gridspec(3, 1, hspace=0.4)
        self.ax_xy = self.fig.add_subplot(gs[0])
        self.ax_pole = self.fig.add_subplot(gs[1])
        self.ax_scal = self.fig.add_subplot(gs[2])
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _load_synthetic(self) -> None:
        self.data = synthesize_polar_motion()
        self.out.delete("1.0", tk.END)
        self.out.insert(tk.END, f"Synthetic polar motion: {len(self.data)} daily points\n")
        self._analyze()

    def _load_csv(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not path:
            return
        try:
            arr = np.loadtxt(path, delimiter=",", skiprows=1)
            self.data = arr[:, :3]
            self.out.delete("1.0", tk.END)
            self.out.insert(tk.END, f"Loaded EOP CSV: {len(self.data)} rows\n")
            self._analyze()
        except Exception as exc:
            messagebox.showerror("Wobble", f"CSV load failed:\n{exc}")

    def _analyze(self) -> None:
        if self.data is None:
            return
        mjd, x, y = self.data[:, 0], self.data[:, 1], self.data[:, 2]
        rx = least_squares_decompose(mjd, x)
        ry = least_squares_decompose(mjd, y)

        self.out.insert(tk.END,
                        "\n=== Decomposition (x-pole) ===\n"
                        f"Chandler amp:  {rx['chandler_amp']*1000:.1f} mas\n"
                        f"Annual   amp:  {rx['annual_amp']*1000:.1f} mas\n"
                        f"Secular slope: {rx['secular_slope_per_yr']*1000:+.2f} mas/yr\n"
                        f"Residual RMS:  {np.std(rx['residual'])*1000:.2f} mas\n\n"
                        "=== Decomposition (y-pole) ===\n"
                        f"Chandler amp:  {ry['chandler_amp']*1000:.1f} mas\n"
                        f"Annual   amp:  {ry['annual_amp']*1000:.1f} mas\n"
                        f"Secular slope: {ry['secular_slope_per_yr']*1000:+.2f} mas/yr\n"
                        f"Residual RMS:  {np.std(ry['residual'])*1000:.2f} mas\n")

        years = (mjd - mjd[0]) / 365.25 + 2010
        self.ax_xy.clear()
        self.ax_xy.plot(years, x, lw=0.6, label="x-pole")
        self.ax_xy.plot(years, y, lw=0.6, label="y-pole")
        self.ax_xy.set_title("Polar motion components")
        self.ax_xy.set_xlabel("year"); self.ax_xy.set_ylabel("arcsec")
        self.ax_xy.legend(); self.ax_xy.grid(alpha=0.3)

        self.ax_pole.clear()
        self.ax_pole.plot(x, y, lw=0.4, color="C2")
        self.ax_pole.set_title("Pole path (x vs y)")
        self.ax_pole.set_xlabel("x-pole (arcsec)"); self.ax_pole.set_ylabel("y-pole (arcsec)")
        self.ax_pole.set_aspect("equal"); self.ax_pole.grid(alpha=0.3)

        periods = np.linspace(200, 600, 80)
        scal = morlet_wavelet_scalogram(x, 1.0, periods)
        self.ax_scal.clear()
        self.ax_scal.imshow(scal, aspect="auto", origin="lower",
                            extent=[years[0], years[-1], periods[0], periods[-1]],
                            cmap="magma")
        self.ax_scal.axhline(CHANDLER_PERIOD_DAYS, ls="--", color="cyan", lw=0.8, label="Chandler")
        self.ax_scal.axhline(ANNUAL_PERIOD_DAYS, ls="--", color="yellow", lw=0.8, label="Annual")
        self.ax_scal.set_title("Morlet scalogram — x-pole")
        self.ax_scal.set_xlabel("year"); self.ax_scal.set_ylabel("period (days)")
        self.ax_scal.legend(loc="upper right", fontsize=8)
        self.canvas.draw_idle()


def main() -> None:
    root = tk.Tk()
    WobbleApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
