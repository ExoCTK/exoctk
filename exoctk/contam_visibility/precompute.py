"""
Module to precompute the contamination for targets, store the data to file, and then retrieve results

Author: Joe Filippazzo
Date: 01/13/26

Example Usage:
filename = "soss_contam.h5"

generate_exoplanet_file(filename)
data = load_exoplanet(filename, "WASP-39b")

print(data["orbital_params"].shape)   # (15,)
print(data["target_trace"].shape)     # (256, 2048)
print(data["contamination"].shape)    # (360, 256, 2048)
"""
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


def generate_exoplanet_file(filename):
    """
    Create an HDF5 file with exoplanet data structures but no data... yet.

    Parameters
    ----------
    filename : str
        Path to output HDF5 file
    """
    # Get names from ExoMAST
    planet_names = fetch_exoplanet_names_eaot()
    count = 0

    # Populate the empty file
    with h5py.File(filename, "w") as f:
        for name in planet_names:

            if name not in f:

                # Make the group and store the name
                grp = f.create_group(name)
                grp.attrs["name"] = name

                # Target trace (256, 2048)
                grp.create_dataset("target_trace", shape=(256, 2048), dtype="float32", compression="gzip",
                                   compression_opts=4)

                # Contamination (360, 256, 2048)
                grp.create_dataset( "contamination", shape=(360, 256, 2048), dtype="float32", compression="gzip",
                                    compression_opts=4, chunks=(1, 256, 2048))
                count += 1

    print(f"Finished! Saved data for {count}/{len(planet_names)} exoplanets to {filename}.")

    return planet_names


def save_exoplanet_data(filename, exoplanet_name, target_trace, contamination):
        """
        Save target trace and contamination data for a single exoplanet.

        Parameters
        ----------
        filename : str
            Path to HDF5 file
        exoplanet_name : str
            Exoplanet name (EAOT name)
        target_trace : ndarray
            2D array of shape (256, 2048)
        contamination : ndarray
            3D array of shape (360, 256, 2048)
        """
        # Open the file and save the data
        with h5py.File(filename, "r+") as f:
            if exoplanet_name not in f:
                raise KeyError(f"Exoplanet '{exoplanet_name}' not found in {filename}")

            grp = f[group_name]

            # Write data in-place (no reallocation)
            grp["target_trace"][:, :] = target_trace
            grp["contamination"][:, :, :] = contamination

            # Optional provenance
            grp.attrs["filled"] = True

        print(f"{exoplanet_name} data saved to {filename}")


def generate_database(filename, aperture='NIS_SUBSTRIP256'):
    """
    Compute the contamination data and save to the file

    Parameters
    ----------
    filename: str
        Path to HDF5 file
    aperture:
        The aperture to use, ['NIS_SUBSTRIP96', 'NIS_SUBSTRIP256', 'NRCA5_GRISM256_F444W', 'NRCA5_GRISM256_F322W2', 'MIRI_SLITLESSPRISM']
    """
    # Make the empty file
    planet_names = generate_exoplanet_file(filename)

    with h5py.File(filename, "w") as f:

        for targname in planet_names:

            # Get the data by name
            targname = get_canonical_name(targname)
            data, _ = get_target_data(targname)

            # Update the coordinates
            ra_deg = data.get('RA')
            dec_deg = data.get('DEC')

            # Run contam tool
            target_trace, contamination, _ = fs.field_simulation(ra_deg, dec_deg, aperture, plot=False)

            # Save the data to file
            save_exoplanet_data(filename, exoplanet_name, target_trace, contamination)
