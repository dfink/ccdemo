# SDSS F Star Spectra Visualization

A Python toolkit for downloading and visualizing stellar spectra from the Sloan Digital Sky Survey (SDSS) Data Release 18, focused on F-type stars with effective temperatures around 6500 K.

## Live Demo

- [Interactive Sky Map](https://dfink.github.io/ccdemo/sky_map.html) - Explore star positions, hover to see spectra
- [Spectrum Gallery](https://dfink.github.io/ccdemo/spectrum_plots/) - Browse all spectrum plots

## Features

### Spectrum Plots
- Publication-quality plots with labeled spectral absorption lines
- Major features marked: Balmer series (Hα, Hβ, Hγ, Hδ), Ca II H & K, Na D, Mg I b
- Color-coded line markers by element type
- Noise panel showing flux uncertainty (σ = 1/√ivar) vs wavelength
- Stellar parameters displayed: Teff, log g, [Fe/H], coordinates, S/N

### Interactive Sky Map
- Plotly-based HTML visualization of star positions (RA vs Dec)
- Points colored by metallicity [Fe/H] using Viridis colorscale
- Hover over any star to see its full spectrum as a popup
- Proper astronomical convention (RA increases right-to-left)

## Requirements

- Python 3.9+
- numpy
- matplotlib
- astropy
- plotly

## Usage

Generate spectrum plots for all FITS files:
```bash
/opt/miniconda3/envs/py39/bin/python plot_spectra.py
```

Generate the interactive sky map:
```bash
/opt/miniconda3/envs/py39/bin/python sky_map.py
```

## Output

- `spectrum_plots/` - PNG plots for each spectrum
- `spectrum_plots/index.html` - Gallery view of all spectra
- `sky_map.html` - Interactive sky map with spectrum popups

## Data

The dataset contains 100 F-type stars from SDSS DR18 with:
- Effective temperature: Teff ≈ 6500 K
- Surface gravity: log g ranging from ~2.2 to ~4.4
- Metallicity: [Fe/H] from -3.1 to 0.0
- Signal-to-noise: 10 to 111

Stellar parameters are stored in `spectra_metadata.json` and were obtained from the SDSS Stellar Parameter Pipeline (SSPP).

## File Structure

```
├── plot_spectra.py        # Spectrum plotting script
├── sky_map.py             # Interactive sky map generator
├── spectra_metadata.json  # Stellar parameters (Teff, log g, [Fe/H], etc.)
├── spec-*.fits            # SDSS spectrum FITS files
├── spectrum_plots/        # Generated PNG plots
│   └── index.html         # Gallery view
└── sky_map.html           # Interactive sky map
```

## License

MIT
