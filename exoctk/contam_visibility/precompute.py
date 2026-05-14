"""
Module to precompute the contamination for targets, store the data to file, and then retrieve results
with per-chunk zero-value masks.

Author: Joe Filippazzo
Date: 01/13/26

Example:
from exoctk.contam_visibility import precompute as p
target_list = p.precomputed_target_list()
p.generate_database(target_list, filename='NIS_SUBSTRIP256_db.h5', aperture='NIS_SUBSTRIP256')

Or for a test database:
p.generate_database(['TRAPPIST-1', 'WASP-18'], filename='NIS_SUBSTRIP256_test_db.h5', aperture='NIS_SUBSTRIP256')
"""

import h5py
import numpy as np
import requests

from . import field_simulator as fs
from ..utils import get_target_data, get_canonical_name
from ..pkgdata import resource_filename


def precomputed_target_list():
    """"
    Read in list of targets to precompute
    """
    # Generate list from file of 2117 targets
    target_list = np.genfromtxt(resource_filename('exoctk', 'data/contam_visibility/contam_precompute_targets.txt'), delimiter=', ', dtype=str)

    # Add a 'b' to each target since Exo.MAST is looking for planets, not stars
    target_list = [f"{t} b" for t in target_list]

    return target_list


def save_exoplanet_data(filename, exoplanet_name, target_trace, contamination, goodPA_list=np.arange(360)):
    """
    Save target trace and contamination (only non-zero planes) to HDF5 file.
    """
    grp_name = exoplanet_name.strip().replace("/", "_")

    with h5py.File(filename, "r+") as f:
        if grp_name not in f:
            raise KeyError(f"Exoplanet '{exoplanet_name}' not found in {filename}")

        grp = f[grp_name]

        # --- Target trace ---
        grp["target_trace"][:, :, :] = target_trace
        grp.attrs["goodPA_list"] = goodPA_list

        # --- Sparse contamination ---
        plane_index = np.where(contamination.any(axis=(1, 2)))[0]

        if plane_index.size == 0:
            # All-zero contamination
            grp["contamination"].resize((0, target_trace.shape[0], target_trace.shape[1]))
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


def generate_database(target_names, filename='NIS_SUBSTRIP256_db.h5', aperture='NIS_SUBSTRIP256'):
    """
    Compute the contamination data and save to the file
    """
    # Get the aperture shape
    if aperture == 'NIS_SUBSTRIP256':
        n_traces, nrows, ncols = 3, 256, 2048
    elif aperture == 'NIS_SUBSTRIP96':
        n_traces, nrows, ncols = 3, 96, 2048
    elif aperture in ['NRCA5_40STRIPE1_DHS_F322W2', 'NRCA5_40STRIPE1_DHS_F444W']:
        n_traces, nrows, ncols = 10, 2128, 3192
    else:
        raise NameError(f"Did not recognize the aperture '{aperture}'")

    # Save canonical name
    lookup = {}

    # Make the empty file
    count = 0
    with h5py.File(filename, "w") as f:
        for targname in target_names:
            try:

                # Canonical name and get coordinates
                name = get_canonical_name(targname)
                data, _ = get_target_data(name)
                ra_deg = data.get('RA')
                dec_deg = data.get('DEC')
                lookup[targname] = {'canonical_name': name, 'ra': ra_deg, 'dec': dec_deg}

                # Make the group in the H5 file
                grp_name = name.strip().replace("/", "_")
                if targname != name:
                    print(f"'{targname}', using '{name}'")
                grp = f.create_group(grp_name)
                grp.attrs["name"] = name
                grp.attrs["ra"] = ra_deg
                grp.attrs["dec"] = dec_deg

                # Target trace
                grp.create_dataset("target_trace", shape=(n_traces, nrows, ncols), dtype="float32", compression="gzip",
                                   compression_opts=4, chunks=(1, nrows, ncols))

                # Contamination placeholder (0 planes initially)
                grp.create_dataset("contamination", shape=(0, nrows, ncols), maxshape=(None, nrows, ncols), dtype="float32",
                                   compression="gzip", compression_opts=4, chunks=(1, nrows, ncols))

                # Plane index placeholder
                grp.create_dataset("plane_index", shape=(0,), maxshape=(None,), dtype="int16")

                count += 1

            except Exception as e:
                print(f"Could not add {name}: \n{e}")

    print(f"Saved structure for {count}/{len(target_names)} exoplanets to {filename}.")

    # Generate the contam figures and save to file
    for targname in target_names:

        # Run contamination tool
        target_traces, contamination, goodPA_list = fs.field_simulation(lookup[targname]['ra'], lookup[targname]['dec'], aperture, plot=False)

        # Save data to file with mask and plane index
        save_exoplanet_data(filename, lookup[targname]['canonical_name'], target_traces, contamination, goodPA_list=goodPA_list)

        print(f"Saved '{targname}' contamination results to {filename}")
