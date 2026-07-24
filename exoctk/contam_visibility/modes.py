"""Observing modes exposed by the contamination and visibility web form."""


CONTAM_VISIBILITY_MODES = (
    ('NIS_SUBSTRIP256', 'NIRISS - SOSS - SUBSTRIP256'),
    ('NIS_SUBSTRIP96', 'NIRISS - SOSS - SUBSTRIP96'),
    ('NRCA5_41STRIPE1_DHS_F322W2', 'NIRCam - DHS - F322W2'),
    ('NRCA5_41STRIPE1_DHS_F444W', 'NIRCam - DHS - F444W'),
    ('NRCA5_GRISM256_F322W2',
     'NIRCam - Grism Time Series - F322W2 (Visibility Only)'),
    ('NRCA5_GRISM256_F444W',
     'NIRCam - Grism Time Series - F444W (Visibility Only)'),
    ('MIRIM_SLITLESSPRISM_IP', 'MIRI - LRS - SLITLESSPRISM_IP'),
    ('NIRSpec', 'NIRSpec (Visibility Only)'),
)
