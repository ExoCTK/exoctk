# gaia_cache.py

from pathlib import Path

import h5py
from astropy.table import Table


class GaiaCache:
    """
    Local cache for Gaia query results.

    Results are stored as Astropy Tables in an HDF5 file, keyed by
    target name.

    Parameters
    ----------
    filename : str or Path
        Path to the HDF5 cache file.
    """

    def __init__(self, filename="gaia_cache.h5"):
        self.filename = Path(filename)

    def contains(self, target):
        """Return True if a result for `target` is cached."""
        target = str(target)

        if not self.filename.exists():
            return False

        with h5py.File(self.filename, "r") as f:
            return target in f

    def save(self, target, table, overwrite=True):
        """
        Save an Astropy Table to the cache.

        Parameters
        ----------
        target : str
            Target name used as the cache key.
        table : astropy.table.Table
            Gaia query results.
        overwrite : bool
            Replace an existing cached result.
        """
        target = str(target)

        self.filename.parent.mkdir(parents=True, exist_ok=True)

        # Remove an existing entry if requested
        if overwrite and self.filename.exists():
            with h5py.File(self.filename, "a") as f:
                if target in f:
                    del f[target]

        # Astropy handles the conversion of the Table to HDF5
        table.write(
            self.filename,
            path=target,
            format="hdf5",
            append=self.filename.exists(),
            overwrite=overwrite,
        )

        print(f"Saved Gaia results for '{target}' to cache.")

    def load(self, target):
        """
        Load a cached Gaia result.

        Returns
        -------
        astropy.table.Table

        Raises
        ------
        KeyError
            If the target is not in the cache.
        """
        target = str(target)

        if not self.filename.exists():
            raise KeyError(f"No cached Gaia result for '{target}'")

        if not self.contains(target):
            raise KeyError(f"No cached Gaia result for '{target}'")

        return Table.read(
            self.filename,
            path=target,
            format="hdf5",
        )

    def get(self, target):
        """
        Return a cached result if available; otherwise run the query.

        Parameters
        ----------
        target : str
            Target name used as the cache key.

        query_function : callable
            Function that performs the Gaia query and returns an
            Astropy Table. It should take `target` as its only argument.

        Returns
        -------
        astropy.table.Table
            Cached or newly queried result.
        """
        if self.contains(target):
            print(f"Loading Gaia results for '{target}' from cache.")
            return self.load(target)

        else:
            return None

    def remove(self, target):
        """Remove a cached target."""
        target = str(target)

        if not self.filename.exists():
            return

        with h5py.File(self.filename, "a") as f:
            if target in f:
                del f[target]

    def clear(self):
        """Delete the entire cache file."""
        if self.filename.exists():
            self.filename.unlink()

    def targets(self):
        """Return a list of all targets currently in the cache."""
        if not self.filename.exists():
            return []

        with h5py.File(self.filename, "r") as f:
            return list(f.keys())