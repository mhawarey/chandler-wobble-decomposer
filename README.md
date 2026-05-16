# chandler-wobble-decomposer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A desktop Tkinter analyzer for **Earth Orientation Parameter (EOP)** polar-motion time-series. Decomposes the (x_pole, y_pole) signal into the **Chandler wobble** (≈ 433-day free Eulerian nutation), the **annual** forced wobble, and **secular drift**, then visualizes the pole path and a Morlet wavelet scalogram showing time-frequency content.

## Features

- Synthetic IERS EOP-C04-like generator (daily, multi-decadal, realistic amplitudes).
- Joint least-squares fit: Chandler + annual + linear + quadratic secular terms.
- Continuous Morlet wavelet transform for time-frequency band tracking.
- Pole-path visualizer with the characteristic spirograph pattern.
- CSV loader for real IERS EOP-C04 series.

## Why this matters

The Chandler wobble is one of geodesy's enduring puzzles: a free oscillation of the rotating Earth that should damp out in ~70 years yet persists, requiring continuous excitation from ocean-bottom pressure and atmospheric processes. Together with the annual wobble it sets the limits of high-precision GNSS, VLBI, and SLR reductions. Tooling to dissect these signals from raw EOP data is foundational for any geodetic researcher.

## Quick start

```bash
pip install -r requirements.txt
python chandler_wobble_decomposer.py
```

CSV format for real data: `MJD, x_pole_arcsec, y_pole_arcsec` (header row required).

## References

- Gross, R. S. (2007). *Earth Rotation Variations — Long Period.* In Treatise on Geophysics, vol. 3 (Geodesy), Elsevier.
- Bizouard, C. & Gambis, D. (2009). *The combined solution C04 for Earth Orientation Parameters.* IERS Tech. Note.
- Chandler, S. C. (1891). *On the variation of latitude.* Astronomical Journal, 11, 65.
- Wilson, C. R. & Vicente, R. O. (1990). *Maximum likelihood estimates of polar motion parameters.* In Variations in Earth Rotation, AGU.

## Author

**Dr. Mosab Hawarey**
>
PhD, Geodetic & Photogrammetric Engineering (ITU) | MSc, Geomatics (Purdue) | MBA (Wales) | BSc, MSc (METU)

- GitHub: https://github.com/mhawarey
- Personal: https://hawarey.org/mosab
- ORCID: https://orcid.org/0000-0001-7846-951X

## License

MIT License
