#!/usr/bin/env python3
"""
Download SDSS DR18 spectra for F stars with Teff near 6500 K.
"""

import os
import requests
from pathlib import Path

# SDSS SkyServer SQL API endpoint
SKYSERVER_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"

# Query for 100 stars with Teff near 6500 K
# Using sppParams table (SSPP - SEGUE Stellar Parameter Pipeline)
SQL_QUERY = """
SELECT TOP 100
    s.specobjid, s.plate, s.mjd, s.fiberid, s.run2d,
    sp.teffadop as teff, sp.loggadop as logg, sp.fehadop as feh,
    s.ra, s.dec, s.snmedian as snr
FROM SpecObjAll s
JOIN sppParams sp ON s.specobjid = sp.specobjid
WHERE sp.teffadop BETWEEN 6300 AND 6700
    AND s.class = 'STAR'
    AND s.snmedian > 10
    AND s.scienceprimary = 1
ORDER BY ABS(sp.teffadop - 6500)
"""

def query_sdss():
    """Query SDSS for F stars with Teff near 6500 K."""
    print("Querying SDSS DR18 for F stars with Teff ~ 6500 K...")

    params = {
        'cmd': SQL_QUERY,
        'format': 'json'
    }

    response = requests.get(SKYSERVER_URL, params=params, timeout=120)
    response.raise_for_status()

    data = response.json()

    if not data or len(data) == 0:
        print("No results returned from query.")
        return []

    # Handle SDSS JSON format: [{"TableName": "Table1", "Rows": [...]}]
    if isinstance(data, list) and len(data) > 0 and "Rows" in data[0]:
        results = data[0]["Rows"]
    elif isinstance(data[0], dict):
        results = data
    else:
        headers = data[0]
        results = [dict(zip(headers, row)) for row in data[1:]]

    print(f"Found {len(results)} matching stars.")
    return results

def download_spectrum(plate, mjd, fiberid, run2d, output_dir):
    """Download a single SDSS spectrum FITS file."""
    # Construct the SAS URL for the spectrum
    # DR18 spectra location format
    plate_str = str(plate).zfill(4)
    fiber_str = str(fiberid).zfill(4)

    # Try the standard BOSS/eBOSS format first
    filename = f"spec-{plate_str}-{mjd}-{fiber_str}.fits"

    # SDSS SAS URL for DR18 spectra
    url = f"https://data.sdss.org/sas/dr18/spectro/sdss/redux/{run2d}/spectra/lite/{plate_str}/{filename}"

    output_path = output_dir / filename

    if output_path.exists():
        print(f"  Already exists: {filename}")
        return True

    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"  Downloaded: {filename}")
            return True
        else:
            # Try alternate URL patterns
            alt_url = f"https://data.sdss.org/sas/dr18/spectro/sdss/redux/{run2d}/spectra/{plate_str}/{filename}"
            response = requests.get(alt_url, timeout=60)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                print(f"  Downloaded: {filename}")
                return True
            print(f"  Failed (HTTP {response.status_code}): {filename}")
            return False
    except Exception as e:
        print(f"  Error downloading {filename}: {e}")
        return False

def main():
    output_dir = Path("/Users/dfink/ccdemo")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Query for stars
    stars = query_sdss()

    if not stars:
        print("No stars found. Exiting.")
        return

    # Save metadata
    import json
    metadata_file = output_dir / "spectra_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(stars, f, indent=2)
    print(f"\nSaved metadata to {metadata_file}")

    # Download spectra
    print(f"\nDownloading {len(stars)} spectra to {output_dir}...")
    success = 0
    failed = 0

    for i, star in enumerate(stars):
        teff = star.get('teff', 0) or 0
        print(f"[{i+1}/{len(stars)}] plate={star['plate']}, mjd={star['mjd']}, fiber={star['fiberid']}, Teff={teff:.0f} K")

        if download_spectrum(
            star['plate'],
            star['mjd'],
            star['fiberid'],
            star['run2d'],
            output_dir
        ):
            success += 1
        else:
            failed += 1

    print(f"\nComplete! Downloaded {success} spectra, {failed} failed.")

if __name__ == "__main__":
    main()
