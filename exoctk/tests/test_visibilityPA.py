from exoctk.contam_visibility import visibilityPA

def test_using_gtvt():
    ra = '88.5897756488149'
    dec = '-64.1297956126618'
    instrument = 'NIRISS'
    output = visibilityPA.using_gtvt(ra, dec, instrument, targetName="Target", output="bokeh")

    assert output is not None
