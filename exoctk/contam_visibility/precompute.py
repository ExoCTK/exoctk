"""
Module to precompute the contamination for targets, store the data to file, and then retrieve results
with per-chunk zero-value masks.

Author: Joe Filippazzo
Date: 01/13/26

Example:
filename = "soss_contam.h5"
data = load_exoplanet(filename, "WASP-39b")

print(data["target_trace"].shape)     # (256, 2048)
print(data["contamination"].shape)    # (360, 256, 2048)
print(len(data["plane_index"]))       # number of non-zero planes
"""

import h5py
import numpy as np
import requests

from . import field_simulator as fs
from ..utils import get_target_data, get_canonical_name


def fetch_exoplanet_names_eaot():
    url = "https://catalogs.mast.stsci.edu/api/v0.1/eaot/search.json"
    params = {"columns": "[Planet_Name]", "pagesize": 20000}

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()
    names = [row[0] for row in data["data"] if row and row[0] is not None]

    return names


def generate_exoplanet_file(filename, test=True):
    """
    Create HDF5 file with empty datasets for exoplanets.

    - target_trace: (256, 2048)
    - contamination: (variable number of planes, 256, 2048), stored only for non-zero planes
    - plane_index: which planes are stored
    """
    planet_names = fetch_exoplanet_names_eaot()
    count = 0

    with h5py.File(filename, "w") as f:
        for name in (planet_names[:3] if test else planet_names):
            if name in f:
                continue

            grp_name = name.strip().replace("/", "_")
            grp = f.create_group(grp_name)
            grp.attrs["name"] = name

            # Target trace
            grp.create_dataset(
                "target_trace",
                shape=(3, 256, 2048),
                dtype="float32",
                compression="gzip",
                compression_opts=4,
                chunks=(1, 256, 2048),
            )

            # Contamination placeholder (0 planes initially)
            grp.create_dataset(
                "contamination",
                shape=(0, 256, 2048),  # variable length along first axis
                maxshape=(None, 256, 2048),
                dtype="float32",
                compression="gzip",
                compression_opts=4,
                chunks=(1, 256, 2048),
            )

            # Plane index placeholder
            grp.create_dataset(
                "plane_index",
                shape=(0,),
                maxshape=(None,),
                dtype="int16",
            )

            count += 1

    print(f"Finished! Saved structure for {count}/{len(planet_names)} exoplanets to {filename}.")
    return planet_names


def save_exoplanet_data(filename, exoplanet_name, target_trace, contamination):
    """
    Save target trace and contamination (only non-zero planes) to HDF5 file.
    """
    grp_name = exoplanet_name.strip().replace("/", "_")

    with h5py.File(filename, "r+") as f:
        if grp_name not in f:
            raise KeyError(f"Exoplanet '{exoplanet_name}' not found in {filename}")

        grp = f[grp_name]

        # --- Target trace (3,256,2048) ---
        grp["target_trace"][:, :, :] = target_trace

        # --- Sparse contamination ---
        plane_index = np.where(contamination.any(axis=(1, 2)))[0]

        if plane_index.size == 0:
            # All-zero contamination
            grp["contamination"].resize((0, 256, 2048))
            grp["plane_index"].resize((0,))
        else:
            nonzero_planes = contamination[plane_index, :, :]

            # Ensure 3D
            if nonzero_planes.ndim == 2:
                nonzero_planes = nonzero_planes[np.newaxis, :, :]

            grp["contamination"].resize(nonzero_planes.shape)
            grp["contamination"][:, :, :] = nonzero_planes

            grp["plane_index"].resize(plane_index.shape)
            grp["plane_index"][:] = plane_index

        grp.attrs["filled"] = True

    print(f"{exoplanet_name} saved ({len(plane_index)} contamination planes)")


def generate_database(filename, aperture='NIS_SUBSTRIP256', test=True):
    """
    Compute the contamination data and save to the file
    """
    # Make the empty file
    planet_names = generate_exoplanet_file(filename, test=True)

    for targname in (planet_names[:3] if test else planet_names):

        # Canonical name and get coordinates
        targname = get_canonical_name(targname)
        data, _ = get_target_data(targname)
        ra_deg = data.get('RA')
        dec_deg = data.get('DEC')

        # Run contamination tool
        target_traces, contamination, _ = fs.field_simulation(ra_deg, dec_deg, aperture, plot=False)

        # Save data to file with mask and plane index
        save_exoplanet_data(filename, targname, target_traces, contamination)
