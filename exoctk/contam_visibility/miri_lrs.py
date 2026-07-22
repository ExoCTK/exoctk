"""MIRI LRS slitless contamination helpers.

Normal contamination calculations use precomputed Pandeia products.  Pandeia
is intentionally imported only by :mod:`make_miri_lrs_traces`, never here.
Arrays in this module follow Pandeia's displayed detector order:
``(row, x_sci)``.  For MIRI/P750L, increasing array row is opposite to
increasing ``Y_SCI``.
"""

from dataclasses import dataclass
from functools import lru_cache
import glob
import json
import os

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import pysiaf


APERTURE = "MIRIM_SLITLESSPRISM_IP"
SHAPE = (384, 68)
DISPERSION_AXIS = 0
CROSS_DISPERSION_AXIS = 1
SOURCE_SEARCH_RADIUS_ARCSEC = 60.0
TRACE_GLOB = "MIRIM_SLITLESSPRISM_IP_*.fits"


class PositionAngleResults(list):
    """Calculated position angles with explicit visibility bookkeeping.

    Parameters
    ----------
    successful : iterable of int, optional
        Position angles for which contamination was calculated successfully.
    inaccessible : iterable of int, optional
        Position angles that are not observable in the requested epoch.

    Attributes
    ----------
    inaccessible : list of int
        The inaccessible position angles.  The list contents inherited by this
        class are the successfully calculated position angles.
    """

    def __init__(self, successful=(), inaccessible=()):
        """Initialize calculated and inaccessible position angles.

        Parameters
        ----------
        successful : iterable of int, optional
            Successfully calculated position angles stored as list contents.
        inaccessible : iterable of int, optional
            Unobservable position angles stored in :attr:`inaccessible`.
        """
        super().__init__(successful)
        self.inaccessible = list(inaccessible)


@dataclass(frozen=True)
class MiriTraceAsset:
    """One target-centered trace template and its calibration products.

    Attributes
    ----------
    trace : numpy.ndarray
        Source-rate image on the native MIRI/IP detector subarray.
    wavelength : numpy.ndarray
        Wavelength associated with each dispersion-axis detector row.
    extraction_mask : numpy.ndarray
        Pixel weights defining the target extraction aperture.
    valid_wavelength : numpy.ndarray
        Boolean mask identifying rows with calibrated target wavelengths.
    reference_pixel : tuple of int
        Zero-based ``(row, column)`` reference pixel stored in the asset.
    metadata : dict
        Provenance and detector-registration metadata.
    red_extension, blue_extension : numpy.ndarray or None
        Off-detector source-rate tails used when a displaced source spectrum
        moves onto the detector.
    red_extension_y0, blue_extension_y0 : int
        Starting detector-array row of each extension before source shifting.
    """

    trace: np.ndarray
    wavelength: np.ndarray
    extraction_mask: np.ndarray
    valid_wavelength: np.ndarray
    reference_pixel: tuple
    metadata: dict
    red_extension: np.ndarray | None = None
    red_extension_y0: int = SHAPE[0]
    blue_extension: np.ndarray | None = None
    blue_extension_y0: int = 0


def trace_directory(exoctk_data=None):
    """Return the external-data directory containing MIRI/IP templates.

    Parameters
    ----------
    exoctk_data : str or path-like, optional
        ExoCTK data root.  When omitted, the ``EXOCTK_DATA`` environment
        variable is used.

    Returns
    -------
    str
        Path to the MIRI/IP trace directory.

    Raises
    ------
    EnvironmentError
        If neither ``exoctk_data`` nor ``EXOCTK_DATA`` supplies a data root.
    """

    root = exoctk_data or os.environ.get("EXOCTK_DATA")
    if not root:
        raise EnvironmentError("EXOCTK_DATA is required to load MIRI/IP traces")
    return os.path.join(root, "exoctk_contam", "traces", APERTURE)


@lru_cache(maxsize=128)
def load_trace(teff, exoctk_data=None):
    """Load the nearest-temperature MIRI/IP trace and validate its metadata.

    Parameters
    ----------
    teff : float
        Requested source effective temperature in kelvin.
    exoctk_data : str or path-like, optional
        ExoCTK data root.  When omitted, ``EXOCTK_DATA`` is used.

    Returns
    -------
    MiriTraceAsset
        Cached trace asset whose tabulated temperature is closest to ``teff``.

    Raises
    ------
    FileNotFoundError
        If no MIRI/IP trace assets are installed.
    ValueError
        If the selected asset has inconsistent shapes, aperture metadata, or
        off-detector extensions.
    """

    files = glob.glob(os.path.join(trace_directory(exoctk_data), TRACE_GLOB))
    if not files:
        raise FileNotFoundError(
            f"No {APERTURE} trace assets found; run make_miri_lrs_traces.py"
        )
    path = min(files, key=lambda item: abs(
        int(os.path.splitext(os.path.basename(item))[0].rsplit("_", 1)[1])
        - teff))
    with fits.open(path) as hdul:
        trace = np.asarray(hdul["TRACE"].data, dtype=float)
        wavelength = np.asarray(hdul["WAVELENGTH"].data, dtype=float).ravel()
        mask = np.asarray(hdul["EXTRACTION_MASK"].data, dtype=float)
        valid = np.asarray(hdul["VALID_WAVELENGTH"].data, dtype=bool).ravel()
        red_extension = (np.asarray(hdul["RED_EXTENSION"].data, dtype=float)
                         if "RED_EXTENSION" in hdul else None)
        blue_extension = (np.asarray(hdul["BLUE_EXTENSION"].data, dtype=float)
                          if "BLUE_EXTENSION" in hdul else None)
        header = hdul[0].header
        metadata = json.loads(header.get("METAJSON", "{}"))
        reference = (int(header["REFYPIX"]), int(header["REFXPIX"]))

    if trace.shape != SHAPE or mask.shape != SHAPE:
        raise ValueError(f"Invalid {APERTURE} asset shape in {path}: {trace.shape}")
    if wavelength.shape != (SHAPE[DISPERSION_AXIS],) or valid.shape != wavelength.shape:
        raise ValueError(f"Invalid wavelength calibration shape in {path}")
    if header.get("APERTURE") != APERTURE:
        raise ValueError(f"Trace asset aperture mismatch in {path}")

    red_y0 = int(metadata.get("red_extension_y0", SHAPE[0]))
    blue_y0 = int(metadata.get("blue_extension_y0", 0))
    if (red_extension is None or red_extension.ndim != 2 or
            red_extension.shape[1] != SHAPE[1] or red_y0 != SHAPE[0]):
        raise ValueError(f"Invalid red-wavelength extension in {path}")
    if (blue_extension is None or blue_extension.ndim != 2 or
            blue_extension.shape[1] != SHAPE[1] or
            blue_y0 + blue_extension.shape[0] != 0):
        raise ValueError(f"Invalid blue-wavelength extension in {path}")

    return MiriTraceAsset(
        trace=trace, wavelength=wavelength, extraction_mask=mask,
        valid_wavelength=valid, reference_pixel=reference, metadata=metadata,
        red_extension=red_extension, red_extension_y0=red_y0,
        blue_extension=blue_extension, blue_extension_y0=blue_y0)


def load_reference_trace(exoctk_data=None):
    """Load the canonical asset used for shared LRS calibration products.

    MIRI/LRS assets have temperature-dependent source traces but share their
    wavelength solution, extraction mask, and generation provenance. Callers
    that need those shared products should use this helper instead of choosing
    an apparently scientific source temperature themselves.

    Parameters
    ----------
    exoctk_data : str or path-like, optional
        ExoCTK data root.  When omitted, ``EXOCTK_DATA`` is used.

    Returns
    -------
    MiriTraceAsset
        The designated calibration-reference trace asset.
    """

    # These products are duplicated in every temperature-dependent trace
    # asset; the selected temperature has no scientific role here.
    reference_temperature = 4000
    if exoctk_data is None:
        return load_trace(reference_temperature)
    return load_trace(reference_temperature, exoctk_data)


def shift_trace(trace, y_shift, x_shift, shape=SHAPE, source_y0=0,
                source_x0=0):
    """Translate a trace by integer detector-array offsets with clipping.

    Parameters
    ----------
    trace : numpy.ndarray
        Two-dimensional source image to translate.
    y_shift, x_shift : int
        Row and column offsets relative to the target-centered image.
    shape : tuple of int, optional
        ``(rows, columns)`` shape of the returned detector image.
    source_y0, source_x0 : int, optional
        Native row and column origin of ``trace`` relative to the detector.

    Returns
    -------
    numpy.ndarray
        Shifted image clipped to ``shape``.
    """

    source = np.asarray(trace)
    output = np.zeros(shape, dtype=source.dtype)
    y_shift = int(y_shift) + int(source_y0)
    x_shift = int(x_shift) + int(source_x0)
    src_y0, src_x0 = max(0, -y_shift), max(0, -x_shift)
    dst_y0, dst_x0 = max(0, y_shift), max(0, x_shift)
    height = min(source.shape[0] - src_y0, shape[0] - dst_y0)
    width = min(source.shape[1] - src_x0, shape[1] - dst_x0)
    if height > 0 and width > 0:
        output[dst_y0:dst_y0 + height, dst_x0:dst_x0 + width] = source[
            src_y0:src_y0 + height, src_x0:src_x0 + width
        ]
    return output


def render_trace(asset, y_shift, x_shift, shape=SHAPE):
    """Translate a trace asset and its off-detector source-rate tails.

    Parameters
    ----------
    asset : MiriTraceAsset
        Trace asset to render.
    y_shift, x_shift : int
        Row and column offsets from the target reference position.
    shape : tuple of int, optional
        ``(rows, columns)`` shape of the returned detector image.

    Returns
    -------
    numpy.ndarray
        Shifted source-rate image including any detector-crossing extensions.
    """

    rendered = shift_trace(asset.trace, y_shift, x_shift, shape=shape)
    if asset.blue_extension is not None:
        rendered += shift_trace(
            asset.blue_extension, y_shift, x_shift, shape=shape,
            source_y0=asset.blue_extension_y0)
    if asset.red_extension is not None:
        rendered += shift_trace(
            asset.red_extension, y_shift, x_shift, shape=shape,
            source_y0=asset.red_extension_y0)
    return rendered


def source_offsets(stars, aperture, attitude):
    """Populate source coordinates and return rigid array offsets from target.

    ``pysiaf`` SCI coordinates are one-indexed pixel-center coordinates.  Their
    differences provide the source displacement in the aperture.  Pandeia's
    displayed P750L detector product has increasing columns along ``+X_SCI``
    but increasing rows along ``-Y_SCI``, so only the Y difference changes
    sign when converted to a NumPy array shift.

    Parameters
    ----------
    stars : astropy.table.Table
        Sources with ``ra`` and ``dec`` columns and writable detector,
        telescope, and science-coordinate columns.  Row zero is the target.
    aperture : pysiaf.aperture.JwstAperture
        MIRI/IP aperture used for coordinate transformations.
    attitude : numpy.ndarray
        SIAF attitude matrix for the requested V3 position angle.

    Returns
    -------
    list of tuple of int
        ``(row_shift, column_shift)`` for each source, beginning with ``(0, 0)``
        for the target.

    Notes
    -----
    Coordinate columns in ``stars`` are populated in place.
    """

    target_xdet, target_ydet = aperture.reference_point("det")
    target_xtel, target_ytel = aperture.det_to_tel(target_xdet, target_ydet)
    target_xsci, target_ysci = aperture.det_to_sci(target_xdet, target_ydet)
    stars["xdet"][0], stars["ydet"][0] = target_xdet, target_ydet
    stars["xtel"][0], stars["ytel"][0] = target_xtel, target_ytel
    stars["xsci"][0], stars["ysci"][0] = target_xsci, target_ysci

    if len(stars) == 1:
        return [(0, 0)]

    # All pysiaf transforms used here accept arrays.  Keeping the operation
    # vectorized avoids thousands of Python/Quantity calls for crowded fields.
    v2, v3 = pysiaf.utils.rotations.sky_to_tel(
        attitude, np.asarray(stars["ra"][1:], dtype=float),
        np.asarray(stars["dec"][1:], dtype=float))
    stars["xtel"][1:] = v2.to_value(u.arcsec)
    stars["ytel"][1:] = v3.to_value(u.arcsec)
    stars["xdet"][1:], stars["ydet"][1:] = aperture.tel_to_det(
        stars["xtel"][1:], stars["ytel"][1:])
    stars["xsci"][1:], stars["ysci"][1:] = aperture.det_to_sci(
        stars["xdet"][1:], stars["ydet"][1:])
    y_offsets = -np.rint(stars["ysci"][1:] - target_ysci).astype(int)
    x_offsets = np.rint(stars["xsci"][1:] - target_xsci).astype(int)
    return [(0, 0), *zip(y_offsets.tolist(), x_offsets.tolist())]


def nearby_source_indices(stars, radius_arcsec):
    """Select the target and catalog sources within a sky radius.

    Parameters
    ----------
    stars : astropy.table.Table
        Source table whose first row is the target and which contains ``ra``
        and ``dec`` columns in degrees.
    radius_arcsec : float
        Maximum angular separation from the target in arcseconds.

    Returns
    -------
    numpy.ndarray
        Integer row indices, always including target row zero when present.
    """

    if len(stars) == 0:
        return np.array([], dtype=int)
    coordinates = SkyCoord(
        np.asarray(stars["ra"], dtype=float) * u.deg,
        np.asarray(stars["dec"], dtype=float) * u.deg)
    separation = coordinates[0].separation(coordinates).to_value(u.arcsec)
    indices = np.flatnonzero(np.isfinite(separation) &
                             (separation <= radius_arcsec))
    if not len(indices) or indices[0] != 0:
        indices = np.insert(indices, 0, 0)
    return indices


def _finite_positive_scalar(value):
    """Return a finite positive float, or ``None`` for unusable values.

    Parameters
    ----------
    value : object
        Candidate source flux scale, including masked table values.

    Returns
    -------
    float or None
        Positive finite value, or ``None`` when the input is unusable.
    """

    if np.ma.is_masked(value):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def calc_v3pa(v3pa, stars, aperture,
              source_search_radius_arcsec=SOURCE_SEARCH_RADIUS_ARCSEC):
    """Render one MIRI/LRS position angle by rigid source translation.

    Parameters
    ----------
    v3pa : float
        Telescope V3 position angle in degrees.
    stars : astropy.table.Table
        Target and neighboring catalog sources. Row zero must be the target.
    aperture : pysiaf.aperture.JwstAperture
        MIRI slitless-prism aperture used for coordinate transformations.
    source_search_radius_arcsec : float
        Maximum separation of candidate contaminants from the target, in
        arcseconds.

    Returns
    -------
    dict
        Target and contaminant detector images, extracted contamination,
        wavelength calibration, extraction mask, and contributing-source
        metadata for the requested position angle.

    Raises
    ------
    ValueError
        If the target does not have a finite, positive flux scale.
    """

    source_indices = nearby_source_indices(
        stars, source_search_radius_arcsec)
    candidate_stars = stars[source_indices].copy()
    xdet, ydet = aperture.reference_point('det')
    xtel, ytel = aperture.det_to_tel(xdet, ydet)
    attitude = pysiaf.utils.rotations.attitude_matrix(
        xtel, ytel, candidate_stars['ra'][0], candidate_stars['dec'][0],
        v3pa % 360)
    offsets = source_offsets(candidate_stars, aperture, attitude)

    target_fluxscale = _finite_positive_scalar(
        candidate_stars['fluxscale'][0])
    if target_fluxscale is None:
        raise ValueError('The target does not have a valid flux scale.')
    target_asset = load_trace(candidate_stars['Teff'][0])
    target = target_asset.trace * target_fluxscale
    contaminants = np.zeros(SHAPE, dtype=float)
    included = []
    intersecting = []
    for index, star, (y_shift, x_shift) in zip(
            source_indices[1:], candidate_stars[1:], offsets[1:]):
        if star['type'] != 'STAR':
            continue
        fluxscale = _finite_positive_scalar(star['fluxscale'])
        if fluxscale is None:
            continue
        # Padded rows extend beyond both detector edges, so only the unpadded
        # cross-dispersion bound is safe to reject before rendering.
        if abs(x_shift) >= SHAPE[1]:
            continue
        asset = load_trace(star['Teff'])
        shifted = render_trace(asset, y_shift, x_shift)
        if not np.any(shifted):
            continue
        contaminants += shifted * fluxscale
        included.append(int(index))
        overlap = shifted * target_asset.extraction_mask
        if np.any(np.isfinite(overlap) & (overlap != 0)):
            columns = star.colnames
            name = next(
                (str(star[column]) for column in ('name', 'source_id')
                 if column in columns),
                f'Source {index}')
            intersecting.append({
                'index': int(index),
                'name': name,
                'ra': float(star['ra']),
                'dec': float(star['dec']),
                'y_shift': int(y_shift),
                'x_shift': int(x_shift),
                'fluxscale': fluxscale,
                'teff': float(star['Teff']),
            })

    fraction = contamination_fraction(
        target, contaminants, target_asset.extraction_mask)
    return {
        'pa': v3pa,
        'target': target,
        'target_traces': [target],
        'contaminants': contaminants,
        'contamination': fraction,
        'wavelength': target_asset.wavelength,
        'valid_wavelength': target_asset.valid_wavelength,
        'extraction_mask': target_asset.extraction_mask,
        'included_sources': included,
        'intersecting_sources': intersecting,
    }


def contamination_fraction(target, contaminants, mask,
                           collapse_axis=CROSS_DISPERSION_AXIS):
    """Calculate extracted contamination for one MIRI trace.

    Parameters
    ----------
    target : numpy.ndarray
        Target source-rate image.
    contaminants : numpy.ndarray
        Summed contaminating-source rate image with the same shape as
        ``target``.
    mask : numpy.ndarray
        Extraction weights with the same shape as ``target``.
    collapse_axis : int, optional
        Detector axis to average over after applying ``mask``.

    Returns
    -------
    numpy.ndarray
        Mean pixel contamination fraction along the retained spectral axis.
        Channels containing no valid extraction pixels are NaN.
    """

    total = target + contaminants
    fraction = np.divide(
        contaminants,
        total,
        out=np.full_like(total, np.nan, dtype=float),
        where=(total != 0) & np.isfinite(total),
    )
    # Pixels outside the target extraction must not participate in the mean.
    # Multiplying them by zero is insufficient: an unrelated trace can make
    # their fractions finite, causing those zeroes to dilute the in-aperture
    # statistic merely because another source exists elsewhere on the array.
    masked_fraction = np.where(mask > 0, fraction * mask, np.nan)
    valid = np.isfinite(masked_fraction)
    numerator = np.nansum(masked_fraction, axis=collapse_axis)
    denominator = np.sum(valid, axis=collapse_axis)
    return np.divide(
        numerator, denominator,
        out=np.full(np.shape(numerator), np.nan, dtype=float),
        where=denominator > 0)
