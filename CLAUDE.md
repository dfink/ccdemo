# SDSS F Star Spectra Project

## Overview
This project downloads and visualizes SDSS DR18 spectra for F-type stars (Teff ~ 6500K).

## Files
- `plot_spectra.py` - Main plotting script
- `sky_map.py` - Interactive sky map generator
- `spectra_metadata.json` - Stellar parameters from SDSS
- `spectrum_plots/` - Generated PNG plots
- `spectrum_plots/index.html` - Gallery view of all spectra
- `sky_map.html` - Interactive sky map with spectrum popups

## Running
```bash
# Generate spectrum plots
/opt/miniconda3/envs/py39/bin/python plot_spectra.py

# Generate interactive sky map
/opt/miniconda3/envs/py39/bin/python sky_map.py
```

## Recent Changes

### Interactive Sky Map (2026-02-05)
- Added `sky_map.py` to generate interactive HTML visualization
- Plotly scatter plot showing star positions (RA vs Dec)
- RA axis reversed (East is left, astronomical convention)
- Points colored by metallicity [Fe/H] with Viridis colorscale
- Hover over any star to see spectrum popup with:
  - Full spectrum PNG image
  - Stellar parameters (Teff, log g, [Fe/H], S/N)
  - Plate-MJD-Fiber identifier
- Popup auto-positions to stay on screen

### Spectrum Plot Improvements (2026-02-05)

**1. Fixed Label Collisions**
- Reduced spectral lines to 8 major features only:
  - Balmer series: Hα, Hβ, Hγ, Hδ
  - Calcium: Ca II H & K
  - Sodium: Na D
  - Magnesium: Mg I b
- Moved info box from upper-left to upper-right (avoids Ca H&K region)
- Labels placed at top of plot with staggered heights for close lines
- Increased font size from 7pt to 11pt

**2. Added Color**
- Spectrum line: blue (#2E86AB) instead of black
- Subtle light blue fill under the curve (alpha=0.15)
- Color-coded spectral line markers:
  - Blue: Balmer lines
  - Orange: Calcium lines
  - Gray: Other metals (Na, Mg)

**3. Added Noise Panel**
- Two-panel layout with 4:1 height ratio
- Bottom panel shows noise (σ = 1/√ivar) vs wavelength
- Coral/light red line (#E07A5F) for noise
- Dashed horizontal line at median noise level
- Noise typically higher at blue/red ends, lower in middle
