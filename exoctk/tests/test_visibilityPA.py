import sys

import pytest

from exoctk.contam_visibility import visibilityPA


@pytest.mark.skipif(sys.version_info > (3, 9), reason='jwst_gtvt does not currently support python>=3.9.')
def test_using_gtvt():
    instrument = 'NIRISS'

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
