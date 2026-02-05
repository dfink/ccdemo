#!/usr/bin/env python3
"""
Plot SDSS spectra with stellar parameters and reference spectral lines.
"""

import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving files
import matplotlib.pyplot as plt
from astropy.io import fits

# Major spectral lines for clean labeling (reduced set to avoid collisions)
SPECTRAL_LINES = {
    # Balmer series (Hydrogen) - major lines only
    r'H$\alpha$': {'wavelength': 6562.8, 'color': '#2E86AB'},
    r'H$\beta$': {'wavelength': 4861.3, 'color': '#2E86AB'},
    r'H$\gamma$': {'wavelength': 4340.5, 'color': '#2E86AB'},
    r'H$\delta$': {'wavelength': 4101.7, 'color': '#2E86AB'},

    # Calcium H & K
    'Ca II K': {'wavelength': 3933.7, 'color': '#E07A00'},
    'Ca II H': {'wavelength': 3968.5, 'color': '#E07A00'},

    # Sodium D doublet
    'Na D': {'wavelength': 5892.9, 'color': '#666666'},

    # Magnesium triplet
    'Mg I b': {'wavelength': 5175.0, 'color': '#666666'},
}


def load_spectrum(fits_path):
    """Load SDSS spectrum from FITS file."""
    with fits.open(fits_path) as hdul:
        # SDSS spectra have data in extension 1 (COADD)
        data = hdul[1].data
        header = hdul[0].header

        # Get wavelength - SDSS uses log-linear wavelength
        # loglam is log10(wavelength in Angstroms)
        if 'loglam' in data.names:
            wavelength = 10**data['loglam']
        elif 'CRVAL1' in header:
            # Reconstruct from header
            coeff0 = header.get('CRVAL1', header.get('COEFF0'))
            coeff1 = header.get('CD1_1', header.get('COEFF1', 0.0001))
            naxis1 = len(data['flux'])
            loglam = coeff0 + coeff1 * np.arange(naxis1)
            wavelength = 10**loglam
        else:
            wavelength = data['loglam'] if 'loglam' in data.names else None

        flux = data['flux']

        # Get inverse variance for error estimation
        ivar = data['ivar'] if 'ivar' in data.names else None

        # Get additional info from header/extensions
        info = {
            'ra': header.get('PLUG_RA', header.get('RA')),
            'dec': header.get('PLUG_DEC', header.get('DEC')),
            'plate': header.get('PLATEID', header.get('PLATE')),
            'mjd': header.get('MJD'),
            'fiberid': header.get('FIBERID'),
            'snr': header.get('SN_MEDIAN', header.get('SPEC1_G')),
        }

        # Try to get reddening from header
        info['ebv'] = header.get('EBV', header.get('SFD_EBV'))

        return wavelength, flux, ivar, info


def plot_spectrum(fits_path, metadata=None, output_dir=None, output_format='png'):
    """
    Create a publication-quality plot of an SDSS spectrum.

    Parameters
    ----------
    fits_path : str or Path
        Path to the FITS file
    metadata : dict, optional
        Stellar parameters from the metadata JSON
    output_dir : str or Path, optional
        Directory to save the plot
    output_format : str
        Output format ('png' or 'svg')

    Returns
    -------
    str : Path to the saved plot
    """
    fits_path = Path(fits_path)

    # Load the spectrum
    wavelength, flux, ivar, info = load_spectrum(fits_path)

    if wavelength is None or len(wavelength) == 0:
        print(f"Could not load wavelength from {fits_path}")
        return None

    # Merge metadata if provided
    if metadata:
        info.update(metadata)

    # Create figure with two panels: main spectrum and noise
    fig, (ax, ax_noise) = plt.subplots(2, 1, figsize=(14, 7), dpi=100,
                                        gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05})

    # Plot the spectrum with color
    ax.plot(wavelength, flux, color='#2E86AB', linewidth=0.6, alpha=0.9)
    # Add subtle fill under the curve
    ax.fill_between(wavelength, flux, alpha=0.15, color='#2E86AB')

    # Calculate y-axis limits (clip outliers)
    valid = np.isfinite(flux) & (flux > 0)
    if valid.sum() > 0:
        flux_valid = flux[valid]
        y_low = np.percentile(flux_valid, 1)
        y_high = np.percentile(flux_valid, 99)
        y_range = y_high - y_low
        ax.set_ylim(y_low - 0.1 * y_range, y_high + 0.35 * y_range)

    # Get y limits for label positioning
    ylim = ax.get_ylim()
    wl_min, wl_max = wavelength.min(), wavelength.max()

    # Sort spectral lines by wavelength for staggering
    sorted_lines = sorted(SPECTRAL_LINES.items(), key=lambda x: x[1]['wavelength'])

    # Add spectral line markers with staggered labels at top
    label_heights = [0.92, 0.84]  # Two stagger levels (as fraction of axis)
    prev_wl = 0
    height_idx = 0
    for name, line_info in sorted_lines:
        wl = line_info['wavelength']
        color = line_info['color']
        if wl_min < wl < wl_max:
            # Draw vertical line
            ax.axvline(wl, color=color, linestyle='-', linewidth=0.8, alpha=0.5)

            # Stagger labels if close together (within 100 Angstroms)
            if wl - prev_wl < 100:
                height_idx = (height_idx + 1) % 2
            else:
                height_idx = 0

            # Place label at top of plot
            ax.text(wl, ylim[0] + label_heights[height_idx] * (ylim[1] - ylim[0]),
                    name, fontsize=11, rotation=90,
                    ha='center', va='bottom', color=color, alpha=0.9,
                    fontweight='medium')
            prev_wl = wl

    # Labels and title (x-label on noise panel since they share x-axis)
    ax.set_ylabel('Flux (10⁻¹⁷ erg/s/cm²/Å)', fontsize=12)
    ax.tick_params(axis='x', labelbottom=False)  # Hide x labels on main panel

    # Build title with stellar parameters
    title_parts = []
    if info.get('plate') and info.get('mjd') and info.get('fiberid'):
        title_parts.append(f"SDSS spec-{info['plate']:04d}-{info['mjd']}-{info['fiberid']:04d}")

    ax.set_title(' | '.join(title_parts) if title_parts else fits_path.stem, fontsize=12)

    # Build parameter text box
    param_lines = []

    # Coordinates
    ra = info.get('ra')
    dec = info.get('dec')
    if ra is not None and dec is not None:
        param_lines.append(f"RA, Dec: {ra:.5f}°, {dec:.5f}°")

    # Stellar parameters
    teff = info.get('teff')
    if teff is not None:
        param_lines.append(f"Teff: {teff:.0f} K")

    logg = info.get('logg')
    if logg is not None:
        param_lines.append(f"log g: {logg:.2f}")

    feh = info.get('feh')
    if feh is not None:
        param_lines.append(f"[Fe/H]: {feh:+.2f}")

    # Alpha enhancement (check various possible keys)
    alpha = info.get('alphafe', info.get('alpha', info.get('afe')))
    if alpha is not None:
        param_lines.append(f"[α/Fe]: {alpha:+.2f}")

    # Reddening
    ebv = info.get('ebv')
    if ebv is not None:
        param_lines.append(f"E(B-V): {ebv:.3f}")

    # SNR
    snr = info.get('snr')
    if snr is not None:
        param_lines.append(f"S/N: {snr:.1f}")

    # Add text box with parameters - positioned upper-right to avoid Ca H&K region
    if param_lines:
        param_text = '\n'.join(param_lines)
        props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='lightgray')
        ax.text(0.98, 0.98, param_text, transform=ax.transAxes, fontsize=9,
               verticalalignment='top', horizontalalignment='right',
               fontfamily='monospace', bbox=props)

    # Style improvements for main panel
    ax.set_xlim(wl_min, wl_max)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_facecolor('#fafafa')

    # === Noise Panel ===
    if ivar is not None:
        # Calculate noise (sigma) from inverse variance
        # Avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            noise = 1.0 / np.sqrt(ivar)
            noise[~np.isfinite(noise)] = np.nan

        # Plot noise as coral/light red line
        ax_noise.plot(wavelength, noise, color='#E07A5F', linewidth=0.6, alpha=0.8)

        # Add median noise reference line
        valid_noise = noise[np.isfinite(noise)]
        if len(valid_noise) > 0:
            median_noise = np.median(valid_noise)
            ax_noise.axhline(median_noise, color='#666666', linestyle='--',
                            linewidth=0.8, alpha=0.6, label=f'median = {median_noise:.2f}')

        # Set y-limits for noise panel (clip extreme values)
        if len(valid_noise) > 0:
            noise_low = np.percentile(valid_noise, 1)
            noise_high = np.percentile(valid_noise, 98)
            noise_range = noise_high - noise_low
            ax_noise.set_ylim(0, noise_high + 0.2 * noise_range)
    else:
        # No ivar available - show empty panel with message
        ax_noise.text(0.5, 0.5, 'No noise data available',
                     transform=ax_noise.transAxes, ha='center', va='center',
                     fontsize=10, color='gray')

    # Noise panel styling
    ax_noise.set_xlim(wl_min, wl_max)  # Match main panel x-axis
    ax_noise.set_xlabel('Wavelength (Å)', fontsize=12)
    ax_noise.set_ylabel('σ (noise)', fontsize=10)
    ax_noise.tick_params(axis='both', which='major', labelsize=10)
    ax_noise.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax_noise.set_facecolor('#fafafa')

    plt.tight_layout()

    # Save the plot
    if output_dir is None:
        output_dir = fits_path.parent / 'plots'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{fits_path.stem}.{output_format}"
    plt.savefig(output_path, format=output_format, dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return str(output_path)


def plot_all_spectra(spectra_dir, metadata_file=None, output_dir=None, output_format='png'):
    """
    Plot all SDSS spectra in a directory.

    Parameters
    ----------
    spectra_dir : str or Path
        Directory containing FITS files
    metadata_file : str or Path, optional
        Path to JSON file with stellar parameters
    output_dir : str or Path, optional
        Directory to save plots
    output_format : str
        Output format ('png' or 'svg')
    """
    spectra_dir = Path(spectra_dir)

    # Load metadata if available
    metadata_lookup = {}
    if metadata_file and Path(metadata_file).exists():
        with open(metadata_file) as f:
            data = json.load(f)
            # Handle SDSS JSON format
            if isinstance(data, list) and len(data) > 0 and 'Rows' in data[0]:
                rows = data[0]['Rows']
            else:
                rows = data

            for row in rows:
                # Create lookup key from plate-mjd-fiber
                key = f"spec-{row['plate']:04d}-{row['mjd']}-{row['fiberid']:04d}"
                metadata_lookup[key] = row

    # Set output directory
    if output_dir is None:
        output_dir = spectra_dir / 'spectrum_plots'
    output_dir = Path(output_dir)

    # Find all FITS files
    fits_files = sorted(spectra_dir.glob('spec-*.fits'))
    print(f"Found {len(fits_files)} spectra to plot")
    print(f"Output directory: {output_dir}")

    for i, fits_path in enumerate(fits_files):
        # Get metadata for this spectrum
        key = fits_path.stem
        metadata = metadata_lookup.get(key, {})

        print(f"[{i+1}/{len(fits_files)}] Plotting {fits_path.name}...", end=' ')

        try:
            output_path = plot_spectrum(fits_path, metadata, output_dir, output_format)
            if output_path:
                print("done")
            else:
                print("failed")
        except Exception as e:
            print(f"error: {e}")

    print(f"\nPlots saved to: {output_dir}")

    # Create a simple HTML index
    create_html_index(output_dir, output_format)


def create_html_index(output_dir, output_format='png'):
    """Create an HTML index page for browsing the plots."""
    output_dir = Path(output_dir)
    plots = sorted(output_dir.glob(f'*.{output_format}'))

    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>SDSS F Star Spectra</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        h1 {{
            color: #333;
            border-bottom: 2px solid #4a90d9;
            padding-bottom: 10px;
        }}
        .info {{
            background: #e8f4f8;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
        }}
        .gallery {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(600px, 1fr));
            gap: 20px;
        }}
        .spectrum {{
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .spectrum img {{
            width: 100%;
            height: auto;
            display: block;
        }}
        .spectrum-name {{
            padding: 10px;
            font-size: 12px;
            color: #666;
            text-align: center;
            border-top: 1px solid #eee;
        }}
        .nav {{
            margin-bottom: 20px;
        }}
        .nav a {{
            margin-right: 15px;
            color: #4a90d9;
        }}
    </style>
</head>
<body>
    <h1>SDSS DR18 F Star Spectra (Teff ~ 6500 K)</h1>
    <div class="info">
        <strong>{len(plots)} spectra</strong> |
        Gray vertical lines mark common stellar absorption features (Balmer series, Ca H&K, Mg, Na D, Fe lines)
    </div>
    <div class="gallery">
'''

    for plot in plots:
        name = plot.stem.replace('spec-', '').replace('-', ' / ')
        html += f'''        <div class="spectrum">
            <a href="{plot.name}" target="_blank">
                <img src="{plot.name}" alt="{plot.stem}">
            </a>
            <div class="spectrum-name">{plot.stem}</div>
        </div>
'''

    html += '''    </div>
</body>
</html>
'''

    index_path = output_dir / 'index.html'
    with open(index_path, 'w') as f:
        f.write(html)

    print(f"HTML index created: {index_path}")


if __name__ == '__main__':
    import sys

    spectra_dir = Path('/Users/dfink/ccdemo')
    metadata_file = spectra_dir / 'spectra_metadata.json'
    output_dir = spectra_dir / 'spectrum_plots'

    # Use PNG for web viewing (fast loading, good quality)
    plot_all_spectra(spectra_dir, metadata_file, output_dir, output_format='png')
