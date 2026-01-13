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
import h5py
import numpy as np
import requests


def fetch_exoplanet_names_eaot():
    url = "https://catalogs.mast.stsci.edu/api/v0.1/eaot/search.json"
    params = {"columns": "[Planet_Name]", "pagesize": 20000}

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()
    names = [row[0] for row in data["data"] if row and row[0] is not None]

    return names


def generate_exoplanet_file(filename: str):
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


def load_exoplanet(filename: str, name: str) -> dict:
    """
    Load exoplanet data by name.

    Parameters
    ----------
    filename : str
        Path to HDF5 file
    name : str
        Exoplanet name (group key)

    Returns
    -------
    dict
        Dictionary containing exoplanet data
    """
    with h5py.File(filename, "r") as f:
        if name not in f:
            raise KeyError(f"Exoplanet '{name}' not found in {filename}")

        # Get the data
        grp = f[name]
        data = {"name": grp.attrs["name"],
                "target_trace": grp["target_trace"][:],
                "contamination": grp["contamination"][:]}

        # TODO: Add slicing for PA ranges to speed up retrieval?

        return data
