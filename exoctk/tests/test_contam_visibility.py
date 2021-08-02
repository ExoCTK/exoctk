import sys

import pytest

from exoctk.contam_visibility import visibilityPA, field_simulator


def test_field_simulator():
    """Test to see that field_simulator.py works for all modes"""
    for mode, aper in field_simulator.APERTURES.items():
        sim, plt = field_simulator.field_simulation("04 25 29.0162", "-30 36 01.603", mode, nPA=5)


@pytest.mark.skipif(sys.version_info > (3, 9), reason='jwst_gtvt does not currently support python>=3.9.')
def test_using_gtvt():
    """Test to see that gtvt works for all instruments"""
    for instrument in ['NIRISS', 'NIRCam', 'NIRSpec', 'MIRI']:

        # this ra/dec has bad PAs
        ra = "-66"
        dec = "44"
        paMin, paMax, gd, fig, table, grouped_badPAs = visibilityPA.using_gtvt(ra, dec, instrument, targetName="Target", output="bokeh")
        assert grouped_badPAs is not None

        # this ra/dec has 100% coverage (no bad PAs)
        ra = '88'
        dec = '-64'
        output = visibilityPA.using_gtvt(ra, dec, instrument, targetName="Target", output="bokeh")

        assert output is not None
        assert len(output) == 6

        paMin, paMax, gd, fig, table, grouped_badPAs = output

        assert len(grouped_badPAs) == 0
