#!/usr/bin/env python
"""
Generate an interactive sky map with spectrum popups.

Creates an HTML page with a Plotly scatter plot showing F star positions
on the sky (RA, Dec). Hovering over a star displays its spectrum PNG.
"""

import json
import plotly.graph_objects as go


def load_metadata(filename='spectra_metadata.json'):
    """Load stellar metadata from JSON file."""
    with open(filename, 'r') as f:
        return json.load(f)


def get_spectrum_filename(star):
    """Generate spectrum PNG filename from star metadata."""
    return f"spec-{star['plate']:04d}-{star['mjd']}-{star['fiberid']:04d}"


def create_sky_map(stars, output_file='sky_map.html'):
    """Create interactive sky map HTML with spectrum hover popups."""

    # Extract data arrays
    ra = [s['ra'] for s in stars]
    dec = [s['dec'] for s in stars]
    feh = [s['feh'] for s in stars]
    snr = [s['snr'] for s in stars]
    teff = [s['teff'] for s in stars]
    logg = [s['logg'] for s in stars]

    # Custom data for hover: [filename, plate, mjd, fiber, teff, logg, feh, snr]
    customdata = [
        [get_spectrum_filename(s), s['plate'], s['mjd'], s['fiberid'],
         s['teff'], s['logg'], s['feh'], s['snr']]
        for s in stars
    ]

    # Create scatter plot
    fig = go.Figure(data=go.Scatter(
        x=ra,
        y=dec,
        mode='markers',
        marker=dict(
            size=10,
            color=feh,
            colorscale='Viridis',
            colorbar=dict(
                title='[Fe/H]',
                titleside='right',
                thickness=15,
                len=0.7
            ),
            line=dict(width=1, color='white')
        ),
        customdata=customdata,
        hovertemplate=(
            '<b>RA:</b> %{x:.3f}<br>'
            '<b>Dec:</b> %{y:.3f}<br>'
            '<b>[Fe/H]:</b> %{customdata[6]:.2f}<extra></extra>'
        )
    ))

    # Update layout
    fig.update_layout(
        title=dict(
            text='SDSS F Stars Sky Map',
            font=dict(size=20),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title='RA (deg)',
            autorange='reversed',  # Astronomical convention: East is left
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(128,128,128,0.3)',
            range=[360, 0]
        ),
        yaxis=dict(
            title='Dec (deg)',
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(128,128,128,0.3)',
            range=[-20, 90]
        ),
        plot_bgcolor='#f8f9fa',
        paper_bgcolor='white',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='monospace'
        ),
        margin=dict(l=60, r=60, t=80, b=60)
    )

    # Generate HTML with custom JavaScript for image popup
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>SDSS F Stars Sky Map</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f0f2f5;
        }}
        #sky-plot {{
            width: 100%;
            height: 700px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        #spectrum-popup {{
            display: none;
            position: fixed;
            z-index: 1000;
            background: white;
            border-radius: 8px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.25);
            padding: 10px;
            max-width: 650px;
            pointer-events: none;
        }}
        #spectrum-popup img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
        }}
        #spectrum-popup .info {{
            font-size: 12px;
            color: #333;
            margin-top: 8px;
            padding-top: 8px;
            border-top: 1px solid #eee;
            font-family: monospace;
        }}
        #spectrum-popup .info span {{
            display: inline-block;
            margin-right: 15px;
        }}
        .header {{
            text-align: center;
            margin-bottom: 15px;
            color: #333;
        }}
        .header h1 {{
            margin: 0 0 5px 0;
            font-size: 24px;
        }}
        .header p {{
            margin: 0;
            color: #666;
            font-size: 14px;
        }}
        .legend {{
            text-align: center;
            margin-top: 10px;
            font-size: 13px;
            color: #666;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>SDSS DR18 F Stars Sky Map</h1>
        <p>Hover over a star to see its spectrum. Points colored by metallicity [Fe/H].</p>
    </div>
    <div id="sky-plot"></div>
    <div id="spectrum-popup">
        <img id="popup-img" src="" alt="Spectrum">
        <div class="info" id="popup-info"></div>
    </div>
    <div class="legend">
        {len(stars)} F-type stars (Teff ~ 6500K) from SDSS DR18
    </div>

    <script>
        // Plotly figure data
        var plotData = {fig.to_json()};

        // Render plot
        Plotly.newPlot('sky-plot', plotData.data, plotData.layout, {{
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        }});

        var plot = document.getElementById('sky-plot');
        var popup = document.getElementById('spectrum-popup');
        var popupImg = document.getElementById('popup-img');
        var popupInfo = document.getElementById('popup-info');

        // Handle hover
        plot.on('plotly_hover', function(data) {{
            var point = data.points[0];
            var cd = point.customdata;
            var filename = cd[0];
            var plate = cd[1];
            var mjd = cd[2];
            var fiber = cd[3];
            var teff = cd[4];
            var logg = cd[5];
            var feh = cd[6];
            var snr = cd[7];

            // Set image source
            popupImg.src = 'spectrum_plots/' + filename + '.png';

            // Set info text
            popupInfo.innerHTML =
                '<span><b>Plate-MJD-Fiber:</b> ' + plate + '-' + mjd + '-' + fiber + '</span>' +
                '<span><b>Teff:</b> ' + teff.toFixed(0) + ' K</span>' +
                '<span><b>log g:</b> ' + logg.toFixed(2) + '</span>' +
                '<span><b>[Fe/H]:</b> ' + feh.toFixed(2) + '</span>' +
                '<span><b>S/N:</b> ' + snr.toFixed(1) + '</span>';

            popup.style.display = 'block';

            // Position popup near cursor
            var evt = data.event;
            positionPopup(evt.clientX, evt.clientY);
        }});

        // Update popup position on mouse move over plot
        plot.addEventListener('mousemove', function(evt) {{
            if (popup.style.display === 'block') {{
                positionPopup(evt.clientX, evt.clientY);
            }}
        }});

        function positionPopup(x, y) {{
            var popupWidth = popup.offsetWidth || 650;
            var popupHeight = popup.offsetHeight || 400;
            var offset = 20;

            // Default: position to the right and below cursor
            var left = x + offset;
            var top = y + offset;

            // Adjust if popup would go off right edge
            if (left + popupWidth > window.innerWidth - 10) {{
                left = x - popupWidth - offset;
            }}

            // Adjust if popup would go off bottom edge
            if (top + popupHeight > window.innerHeight - 10) {{
                top = y - popupHeight - offset;
            }}

            // Ensure not off left or top edge
            left = Math.max(10, left);
            top = Math.max(10, top);

            popup.style.left = left + 'px';
            popup.style.top = top + 'px';
        }}

        // Hide popup when not hovering
        plot.on('plotly_unhover', function() {{
            popup.style.display = 'none';
        }});
    </script>
</body>
</html>
'''

    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"Created {output_file}")
    print(f"Plotted {len(stars)} stars")
    print(f"[Fe/H] range: {min(feh):.2f} to {max(feh):.2f}")
    print(f"RA range: {min(ra):.1f} to {max(ra):.1f} deg")
    print(f"Dec range: {min(dec):.1f} to {max(dec):.1f} deg")


def main():
    stars = load_metadata()
    create_sky_map(stars)


if __name__ == '__main__':
    main()
