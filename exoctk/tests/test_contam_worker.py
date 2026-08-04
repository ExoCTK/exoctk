"""Tests for contamination-worker artifact transport."""

import pickle

import numpy as np
import pytest

from exoctk.contam_visibility import field_simulator

app_exoctk = pytest.importorskip(
    "exoctk.exoctk_app.app_exoctk",
    reason="Worker transport tests require the optional web stack")


def _run_task(monkeypatch, params, calculation_result):
    """Run the bound Celery task eagerly without contacting its backend."""

    monkeypatch.setattr(
        app_exoctk.fs, "field_simulation",
        lambda **kwargs: calculation_result)
    monkeypatch.setattr(
        app_exoctk.run_contam_visibility_task,
        "update_state", lambda **kwargs: None)
    task = app_exoctk.run_contam_visibility_task
    task.push_request(id="contam-test")
    try:
        return task.run(params)
    finally:
        task.pop_request()


def test_dhs_worker_persists_only_compact_artifacts(monkeypatch, tmp_path):
    """DHS arrays stay out of Redis and unused target images stay off disk."""

    monkeypatch.setenv("SHARED_DATA_DIR", str(tmp_path))
    targets = [np.ones((3, 4)) for _ in range(10)]
    contamination = field_simulator.DHSContaminationResult(
        order_fractions=tuple(np.zeros((360, 4)) for _ in range(10)),
        position_angles=np.arange(360))
    position_angles = np.array([0, 1, 33])

    returned = _run_task(
        monkeypatch,
        {"aperture": "NRCA5_41STRIPE1_DHS_F322W2"},
        (targets, contamination, position_angles))

    assert returned == "contam-test"
    assert not (tmp_path / "contam-test_targframe.pickle").exists()
    assert not (tmp_path / "contam-test_starcube.pickle").exists()
    artifact = tmp_path / "contam-test_contamination.pickle"
    with artifact.open("rb") as f:
        restored = pickle.load(f)
    assert isinstance(restored, field_simulator.DHSContaminationResult)
    for before, after in zip(contamination, restored):
        np.testing.assert_array_equal(after, before)
    with (tmp_path / "contam-test_results.pickle").open("rb") as f:
        np.testing.assert_array_equal(pickle.load(f), position_angles)
    assert not list(tmp_path.glob(".contam-*"))


def test_non_dhs_worker_preserves_legacy_artifacts(monkeypatch, tmp_path):
    """Established modes retain target-frame and dense-cube artifacts."""

    monkeypatch.setenv("SHARED_DATA_DIR", str(tmp_path))
    targets = [np.ones((3, 4))]
    starcube = np.zeros((360, 3, 4))
    position_angles = np.array([10, 11])

    returned = _run_task(
        monkeypatch,
        {"aperture": "NIS_SUBSTRIP256"},
        (targets, starcube, position_angles))

    assert returned == "contam-test"
    assert (tmp_path / "contam-test_targframe.pickle").exists()
    assert (tmp_path / "contam-test_starcube.pickle").exists()
    assert not (tmp_path / "contam-test_contamination.pickle").exists()
    assert (tmp_path / "contam-test_results.pickle").exists()
    assert not list(tmp_path.glob(".contam-*"))


def test_atomic_pickle_removes_temporary_file_on_failure(
        monkeypatch, tmp_path):
    """A serialization error cannot leave a partial artifact behind."""

    def fail_dump(*args, **kwargs):
        raise RuntimeError("synthetic serialization failure")

    monkeypatch.setattr(app_exoctk.pickle, "dump", fail_dump)
    output = tmp_path / "result.pickle"
    with pytest.raises(RuntimeError, match="synthetic serialization failure"):
        app_exoctk._atomic_pickle_dump(object(), str(output))

    assert not output.exists()
    assert not list(tmp_path.glob(".contam-*"))
