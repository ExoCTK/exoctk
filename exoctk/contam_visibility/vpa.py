def calc_v3pa(ra, dec, aperture, pa, plot=True):
    # Converting from decimal degrees
    targetcrd = crd.SkyCoord(ra=ra, dec=dec, unit='deg').to_string('hmsdms')

    # Querying for neighbors with 2MASS IRSA's fp_psc (point-source catalog)
    info = Irsa.query_region(targetcrd, catalog='fp_psc', spatial='Cone', radius=2.5 * u.arcmin)
    # Use Gaia instead?

    # Coordinates of all stars in FOV, including target
    allRA = info['ra'].data.data
    allDEC = info['dec'].data.data

    # Initiating a dictionary to hold all relevant star information
    stars = {}
    stars['RA'], stars['DEC'] = allRA, allDEC

    # Finding the target index which will be used to select the target in our position dictionaries
    sindRA = (ra - stars['RA']) * np.cos(dec)
    cosdRA = dec - stars['DEC']
    stars['distance'] = np.sqrt(sindRA ** 2 + cosdRA ** 2)
    targetIndex = np.argmin(stars['distance'])

    # 2MASS - Teff relations
    jhMod = np.array(
        [0.545, 0.561, 0.565, 0.583, 0.596, 0.611, 0.629, 0.642, 0.66, 0.679, 0.696, 0.71, 0.717, 0.715, 0.706, 0.688,
         0.663, 0.631, 0.601, 0.568, 0.537, 0.51, 0.482, 0.457, 0.433, 0.411, 0.39, 0.37, 0.314, 0.279])
    hkMod = np.array(
        [0.313, 0.299, 0.284, 0.268, 0.257, 0.247, 0.24, 0.236, 0.229, 0.217, 0.203, 0.188, 0.173, 0.159, 0.148, 0.138,
         0.13, 0.123, 0.116, 0.112, 0.107, 0.102, 0.098, 0.094, 0.09, 0.086, 0.083, 0.079, 0.07, 0.067])
    teffMod = np.array(
        [2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500,
         4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5800, 6000])

    # JHK bands of all stars in FOV, including target
    Jmag = info['j_m'].data.data
    Hmag = info['h_m'].data.data
    Kmag = info['k_m'].data.data

    # J-H band, H-K band. This will be used to derive the stellar Temps later
    J_Hobs = Jmag - Hmag
    H_Kobs = Hmag - Kmag

    # Number of stars
    nStars = stars['RA'].size

    # Find/assign Teff of each star
    starsT = np.empty(nStars)
    for j in range(nStars):
        color_separation = (J_Hobs[j] - jhMod) ** 2 + (H_Kobs[j] - hkMod) ** 2
        min_separation_ind = np.argmin(color_separation)
        starsT[j] = teffMod[min_separation_ind]

    # Record keeping
    stars['Temp'] = starsT
    inst, full_aper, xleft, xright, ybot, ytop = APERTURES[aperture]

    # Instantiate SIAF object
    siaf = pysiaf.Siaf(inst)

    aper = siaf.apertures[aperture]
    full = siaf.apertures[full_aper]

    # DET_targ -> TEL_targ -> get attitude matrix for target -> TEL_neighbor -> DET_neighbor -> SCI_neighbor
    xSweet, ySweet = aper.reference_point('det')
    v2targ, v3targ = aper.det_to_tel(xSweet, ySweet)
    attitude = pysiaf.utils.rotations.attitude_matrix(v2targ, v3targ, ra, dec, pa)

    xdet, ydet = [], []
    xsci, ysci = [], []
    for starRA, starDEC in zip(stars['RA'], stars['DEC']):
        # Get the TEL coordinates of each star using the attitude matrix of the target
        V2, V3 = pysiaf.utils.rotations.sky_to_tel(attitude, starRA, starDEC)
        # Convert to arcsec and turn to a float
        V2, V3 = V2.to(u.arcsec).value, V3.to(u.arcsec).value

        XDET, YDET = aper.tel_to_det(V2, V3)
        XSCI, YSCI = aper.det_to_sci(XDET, YDET)

        xdet.append(XDET)
        ydet.append(YDET)
        xsci.append(XSCI)
        ysci.append(YSCI)

    stars['xdet'], stars['ydet'] = np.array(xdet), np.array(ydet)
    stars['xsci'], stars['ysci'] = np.array(xsci), np.array(ysci)

    # Full Frame dimensions
    rows, cols = full.corners('det')
    minrow, maxrow = rows.min(), rows.max()
    mincol, maxcol = cols.min(), cols.max()

    # Just stars in FOV (Should always have at least 1, the target)
    inFOV = []
    for star in range(nStars):

        x, y = stars['xdet'][star], stars['ydet'][star]
        #         print(x, y, mincol, maxcol, minrow, maxrow)
        if (mincol < x) & (x < maxcol) & (minrow < y) & (y < maxrow):
            inFOV.append(star)

    # Only keep stars in FOV
    FOVstars = {k: v[inFOV] for k, v in stars.items()}
    #     FOVstars = stars

    # Figure
    if plot:
        cds = ColumnDataSource(data=FOVstars)
        mapper = linear_cmap(field_name='Temp', palette=Spectral6, low=min(starsT), high=max(starsT))
        fig = figure(title='Generated FOV from 2MASS IRSA fp_psc - RA: {}, DEC: {}'.format(ra, dec), match_aspect=True)
        fullstr = str(full.AperName.replace('_', ' '))
        fig.title = 'The FOV in SCIENCE coordinates at APA {} for {}'.format(pa, aperture)
        #     print('aper corners', aper.corners('sci'))
        #     print('full corners', full.corners('sci'))
        fig.patch(full.corners('sci')[0], full.corners('sci')[1], color="black", alpha=0.1)
        fig.patch(aper.corners('sci')[0], aper.corners('sci')[1], line_color="blue", fill_color='blue', fill_alpha=0.1)

        # Fine-tune trace dims
    trace = get_trace(aperture)
    #     print(trace.shape, xleft, xright, ybot, ytop)
    trace = trace[xleft:-xright, ybot:-ytop]

    # Get subarray dims
    subX, subY = aper.XSciSize, aper.YSciSize

    # Make frame
    frame = np.zeros((subY, subX))

    # Plotting the trace footprints
    for n, (x, y) in enumerate(zip(FOVstars['xsci'], FOVstars['ysci'])):

        if 'NIS' in aperture:
            ptrace = trace.T
            height, width = ptrace.shape
            x0 = x - width + 68
            y0 = y - height

        elif 'F322W2' in aperture:
            ptrace = trace
            height, width = ptrace.shape
            x0 = x - width + 467  # 2048 - 1581
            y0 = y - height / 2

        elif 'F356W' in aperture:
            ptrace = trace
            height, width = ptrace.shape
            x0 = x - width + 467  # 2048 - 1581
            y0 = y - height / 2

        elif 'F277W' in aperture:
            ptrace = trace
            height, width = ptrace.shape
            x0 = x - width - 600
            y0 = y - height / 2

        elif 'F444W' in aperture:
            ptrace = trace.T[:, ::-1]
            height, width = ptrace.shape
            x0 = x - width + 1096  # 2048 - 952
            y0 = y - height / 2

        else:
            ptrace = trace
            height, width = ptrace.shape
            x0 = x - width / 2
            y0 = y - height + 113  # 529 - 416

        if plot:
            fig.image(image=[ptrace], x=x0, dw=width, y=y0, dh=height, alpha=0.5)

    # Stars
    if plot:
        fig.star('xsci', 'ysci', color=mapper, source=cds, size=12)
        fig.circle(stars['xsci'][targetIndex], stars['ysci'][targetIndex], size=15, fill_color=None, line_color='red')
        color_bar = ColorBar(color_mapper=mapper['transform'], width=8, location=(0, 0), title="Teff")
        fig.add_layout(color_bar, 'right')

        show(fig)

    return frame