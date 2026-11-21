# test_gaia_cache.py

import logging

import numpy as np
import pytest
from astropy.table import Table

from exoctk.contam_visibility.gaia_cache import GaiaCache


@pytest.fixture
def cache(tmp_path):
    """Return a GaiaCache using a temporary HDF5 file."""
    return GaiaCache(tmp_path / "gaia_cache.h5")


@pytest.fixture
def sample_table():
    """Return a representative Gaia-like Astropy Table."""
    return Table(
        {
            "source_id": np.array([123456789, 987654321], dtype=np.int64),
            "ra": np.array([10.123, 20.456]),
            "dec": np.array([-5.678, 30.123]),
            "phot_g_mean_mag": np.array([12.34, 15.67]),
        }
    )


# ----------------------------------------------------------------------
# Initialization
# ----------------------------------------------------------------------


def test_init_default_filename():
    """The default cache filename should be gaia_cache.h5."""
    cache = GaiaCache()

    assert cache.filename.name == "gaia_cache.h5"
    assert cache.path == "gaia_cache.h5"


def test_init_custom_filename(tmp_path):
    """A custom filename should be stored as a Path and as a string path."""
    filename = tmp_path / "custom_cache.h5"
    cache = GaiaCache(filename)

    assert cache.filename == filename
    assert cache.path == str(filename)


# ----------------------------------------------------------------------
# contains()
# ----------------------------------------------------------------------


def test_contains_returns_false_when_cache_does_not_exist(cache):
    """contains() should return False when the cache file doesn't exist."""
    assert not cache.filename.exists()
    assert cache.contains("Trappist-1") is False


def test_contains_returns_false_for_uncached_target(cache, sample_table):
    """contains() should return False for a target that hasn't been saved."""
    cache.save("Trappist-1", sample_table)

    assert cache.contains("Other-Target") is False


def test_contains_returns_true_for_cached_target(cache, sample_table):
    """contains() should return True after a target is saved."""
    cache.save("Trappist-1", sample_table)

    assert cache.contains("Trappist-1") is True


def test_contains_converts_target_to_string(cache, sample_table):
    """contains() should accept targets that can be converted to strings."""
    cache.save(12345, sample_table)

    assert cache.contains(12345) is True
    assert cache.contains("12345") is True


# ----------------------------------------------------------------------
# save()
# ----------------------------------------------------------------------


def test_save_creates_cache_file(cache, sample_table):
    """save() should create the HDF5 cache file."""
    assert not cache.filename.exists()

    cache.save("Trappist-1", sample_table)

    assert cache.filename.exists()


def test_save_creates_parent_directories(tmp_path, sample_table):
    """save() should create missing parent directories."""
    filename = tmp_path / "nested" / "directory" / "gaia_cache.h5"
    cache = GaiaCache(filename)

    cache.save("Trappist-1", sample_table)

    assert filename.exists()
    assert cache.contains("Trappist-1")


def test_save_and_load_preserves_table(cache, sample_table):
    """A saved table should be recovered unchanged by load()."""
    cache.save("Trappist-1", sample_table)

    loaded = cache.load("Trappist-1")

    assert isinstance(loaded, Table)
    assert loaded.colnames == sample_table.colnames
    assert len(loaded) == len(sample_table)

    for column in sample_table.colnames:
        np.testing.assert_array_equal(
            loaded[column],
            sample_table[column],
        )


def test_save_multiple_targets(cache, sample_table):
    """Multiple targets should coexist in the same cache."""
    table2 = Table(
        {
            "source_id": np.array([111111111], dtype=np.int64),
            "ra": np.array([1.0]),
            "dec": np.array([2.0]),
        }
    )

    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", table2)

    assert cache.contains("Trappist-1")
    assert cache.contains("Proxima-Centauri")

    loaded1 = cache.load("Trappist-1")
    loaded2 = cache.load("Proxima-Centauri")

    assert len(loaded1) == 2
    assert len(loaded2) == 1


def test_save_overwrites_existing_target(cache, sample_table):
    """Saving an existing target with overwrite=True should replace it."""
    cache.save("Trappist-1", sample_table)

    replacement = Table(
        {
            "source_id": np.array([42], dtype=np.int64),
            "ra": np.array([99.0]),
            "dec": np.array([88.0]),
        }
    )

    cache.save("Trappist-1", replacement, overwrite=True)

    loaded = cache.load("Trappist-1")

    assert len(loaded) == 1
    np.testing.assert_array_equal(
        loaded["source_id"],
        replacement["source_id"],
    )
    np.testing.assert_array_equal(
        loaded["ra"],
        replacement["ra"],
    )
    np.testing.assert_array_equal(
        loaded["dec"],
        replacement["dec"],
    )


def test_save_overwrite_does_not_remove_other_targets(cache, sample_table):
    """Overwriting one target should leave other cached targets intact."""
    table2 = Table(
        {
            "source_id": np.array([42], dtype=np.int64),
            "ra": np.array([99.0]),
            "dec": np.array([88.0]),
        }
    )

    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", table2)

    replacement = Table(
        {
            "source_id": np.array([123], dtype=np.int64),
            "ra": np.array([50.0]),
            "dec": np.array([60.0]),
        }
    )

    cache.save("Trappist-1", replacement, overwrite=True)

    assert cache.contains("Trappist-1")
    assert cache.contains("Proxima-Centauri")

    loaded_other = cache.load("Proxima-Centauri")
    np.testing.assert_array_equal(
        loaded_other["source_id"],
        table2["source_id"],
    )


# ----------------------------------------------------------------------
# load()
# ----------------------------------------------------------------------


def test_load_missing_cache_raises_key_error(cache):
    """load() should raise KeyError if the cache file doesn't exist."""
    with pytest.raises(KeyError, match="No cached Gaia result for 'Trappist-1'"):
        cache.load("Trappist-1")


def test_load_missing_target_raises_key_error(cache, sample_table):
    """load() should raise KeyError for a target that isn't cached."""
    cache.save("Trappist-1", sample_table)

    with pytest.raises(KeyError, match="No cached Gaia result for 'Other-Target'"):
        cache.load("Other-Target")


def test_load_converts_target_to_string(cache, sample_table):
    """load() should accept targets that can be converted to strings."""
    cache.save("12345", sample_table)

    loaded = cache.load(12345)

    assert len(loaded) == len(sample_table)


# ----------------------------------------------------------------------
# get()
# ----------------------------------------------------------------------


def test_get_returns_none_for_uncached_target(cache):
    """get() should return None when the target isn't cached."""
    assert cache.get("Trappist-1") is None


def test_get_returns_cached_table(cache, sample_table):
    """get() should return the cached table when available."""
    cache.save("Trappist-1", sample_table)

    result = cache.get("Trappist-1")

    assert isinstance(result, Table)
    assert len(result) == len(sample_table)

    for column in sample_table.colnames:
        np.testing.assert_array_equal(
            result[column],
            sample_table[column],
        )


def test_get_logs_cache_hit(cache, sample_table, caplog):
    """get() should log when loading a cached target."""
    cache.save("Trappist-1", sample_table)

    with caplog.at_level(logging.INFO):
        result = cache.get("Trappist-1")

    assert result is not None
    assert "Loading Gaia results for 'Trappist-1' from cache" in caplog.text


def test_get_does_not_log_cache_hit_for_missing_target(cache, caplog):
    """get() should not emit the cache-hit message for a missing target."""
    with caplog.at_level(logging.INFO):
        result = cache.get("Trappist-1")

    assert result is None
    assert "Loading Gaia results" not in caplog.text


# ----------------------------------------------------------------------
# remove()
# ----------------------------------------------------------------------


def test_remove_deletes_target(cache, sample_table):
    """remove() should delete the requested target."""
    cache.save("Trappist-1", sample_table)

    assert cache.contains("Trappist-1")

    cache.remove("Trappist-1")

    assert not cache.contains("Trappist-1")


def test_remove_does_not_delete_other_targets(cache, sample_table):
    """remove() should leave other cached targets untouched."""
    table2 = Table(
        {
            "source_id": np.array([42], dtype=np.int64),
            "ra": np.array([1.0]),
            "dec": np.array([2.0]),
        }
    )

    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", table2)

    cache.remove("Trappist-1")

    assert not cache.contains("Trappist-1")
    assert cache.contains("Proxima-Centauri")


def test_remove_missing_target_does_nothing(cache, sample_table):
    """remove() should not fail when the target isn't cached."""
    cache.save("Trappist-1", sample_table)

    cache.remove("Other-Target")

    assert cache.contains("Trappist-1")


def test_remove_when_cache_does_not_exist(cache):
    """remove() should not fail when the cache file doesn't exist."""
    cache.remove("Trappist-1")

    assert not cache.filename.exists()


# ----------------------------------------------------------------------
# clear()
# ----------------------------------------------------------------------


def test_clear_deletes_cache_file(cache, sample_table):
    """clear() should delete the entire cache file."""
    cache.save("Trappist-1", sample_table)

    assert cache.filename.exists()

    cache.clear()

    assert not cache.filename.exists()


def test_clear_when_cache_does_not_exist(cache):
    """clear() should not fail when there is no cache file."""
    cache.clear()

    assert not cache.filename.exists()


def test_clear_removes_all_targets(cache, sample_table):
    """clear() should remove all cached targets."""
    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", sample_table)

    assert len(cache.targets) == 2

    cache.clear()

    assert cache.targets == []


# ----------------------------------------------------------------------
# targets property
# ----------------------------------------------------------------------


def test_targets_empty_when_cache_does_not_exist(cache):
    """targets should be an empty list for a nonexistent cache."""
    assert cache.targets == []


def test_targets_returns_cached_targets(cache, sample_table):
    """targets should contain all cached target names."""
    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", sample_table)
    cache.save("Barnard-Star", sample_table)

    assert set(cache.targets) == {
        "Trappist-1",
        "Proxima-Centauri",
        "Barnard-Star",
    }


def test_targets_updates_after_remove(cache, sample_table):
    """targets should reflect removed entries."""
    cache.save("Trappist-1", sample_table)
    cache.save("Proxima-Centauri", sample_table)

    cache.remove("Trappist-1")

    assert cache.targets == ["Proxima-Centauri"]


# ----------------------------------------------------------------------
# Astropy Table edge cases
# ----------------------------------------------------------------------


def test_save_empty_table(cache):
    """An empty Astropy Table should be cacheable."""
    table = Table(
        {
            "source_id": np.array([], dtype=np.int64),
            "ra": np.array([], dtype=float),
            "dec": np.array([], dtype=float),
        }
    )

    cache.save("Empty-Target", table)

    loaded = cache.load("Empty-Target")

    assert isinstance(loaded, Table)
    assert len(loaded) == 0
    assert loaded.colnames == table.colnames


def test_save_table_with_metadata(cache):
    """Table metadata should survive the cache round trip."""
    table = Table(
        {
            "source_id": [123],
            "ra": [10.0],
            "dec": [20.0],
        }
    )
    table.meta["target"] = "Trappist-1"
    table.meta["query_radius"] = 5.0

    cache.save("Trappist-1", table)

    loaded = cache.load("Trappist-1")

    assert loaded.meta["target"] == "Trappist-1"
    assert loaded.meta["query_radius"] == 5.0
