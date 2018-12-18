def plot_spectral_overlap_F322W2(contaminationF322W2_order1, contaminationF322W2_order2, badPA = [], col, dra, ddec, Decd0, PA, msize, dmag, targetName):
    """This function will use the contamination output for NIRCam's F322W2 GTS
    data to create a bokeh plot that will be outputted in ExoCTK.stsci.edu's
    Contamination Overlap page (when the user opts to plot the contamination
    plot, in addition to the visibility plot).

    Credits
    -------
    Written by Joseph Filipazzo, 201?
    Edited by Jennifer V. Medina, 2018

    Parameters
    ----------

    Returns
    -------
    A bokeh plot
    """
    TOOLS = 'pan, box_zoom, crosshair, reset, hover, save'

    y = np.array([0., 0.])
    y1 = 0.07
    y2 = 0.12
    y3 = 0.17
    y4 = 0.23
    bad_PA_color = '#dddddd'
    bad_PA_alpha = 0.7

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Initializations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    loc = 'data/contam_visibility/lambda_order1-2.txt'
    lam_file = pkg_resources.resource_filename('exoctk', loc)
    ypix, lamO1, lamO2 = np.loadtxt(lam_file, unpack=True)

    #ny = trace2dO1.shape[0] # not necessary for this simulation
    contam01 = contaminationF322W2_order1
    contam02 = contaminationF322W2_order2
    nPA = contam01.shape[0]
    dPA = 360//nPA
    PA = np.arange(nPA)*dPA

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 1 (plots on the far right)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Contam plot (PA vs. wavel)
    xlim0 = lamO1.min()
    xlim1 = lamO1.max()
    ylim0 = PA.min()-0.5*dPA
    ylim1 = PA.max()+0.5*dPA
    color_mapper = LinearColorMapper(palette=inferno(8)[::-1],
                                     low=-4, high=1)
    color_mapper.low_color = 'white'
    color_mapper.high_color = 'black'
    s2 = figure(tools=TOOLS, width=500, height=500, title=None,
                x_range=Range1d(xlim0, xlim1),
                y_range=Range1d(ylim0, ylim1))
    fig_data = np.log10(np.clip(contamO1.T, 1.e-10, 1.))
    s2.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
             color_mapper=color_mapper)
    s2.xaxis.axis_label = 'Wavelength (um)'
    s2.yaxis.axis_label = 'Position Angle (degrees)'

    # Line plot
    s3 = figure(tools=TOOLS, width=150, height=500,
                x_range=Range1d(0, 100), y_range=s2.y_range, title=None)
    s3.line(100*np.sum(contamO1 >= 0.001, axis=0)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s3.line(100*np.sum(contamO1 >= 0.01, axis=0)/ny, PA-dPA/2,
            line_color='green', legend='> 0.01')
    s3.xaxis.axis_label = '% channels contam.'
    s3.yaxis.major_label_text_font_size = '0pt'

    # Add bad PAs
    for ybad0, ybad1 in badPA:
        s2.patch([xlim0, xlim1, xlim1, xlim0],
                 [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha)
        s3.patch([0, 100, 100, 0], [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha, legend='Bad PA')

    # Line list
    s1 = figure(tools=TOOLS, width=500, plot_height=100,
                x_range=s2.x_range, title=None)

    l = np.array([0.89, 0.99])
    s1.line(l, y+y1, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([1.09, 1.2])
    s1.line(l, y+y1, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([1.1, 1.24])
    s1.line(l, y+y2, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                  text='CH4', render_mode='css', text_font_size='8pt'))

    l = np.array([1.3, 1.51])
    s1.line(l, y+y1, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([1.6, 1.8])
    s1.line(l, y+y2, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                  text='CH4', render_mode='css', text_font_size='8pt'))

    l = np.array([1.75, 2.05])
    s1.line(l, y+y1, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([2.3, lamO1.max()])
    s1.line(l, y+y1, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([2.15, 2.5])
    s1.line(l, y+y2, line_color='black', line_width=1.5)
    s1.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                  text='CH4', render_mode='css', text_font_size='8pt'))

    l = np.array([1.1692, 1.1778])
    s1.line(l[0], [y3, y3+0.02], line_color='black')
    s1.line(l[1], [y3, y3+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.2437, 1.2529])
    s1.line(l[0], [y3, y3+0.02], line_color='black')
    s1.line(l[1], [y3, y3+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.5168])
    s1.line(l[0], [y3, y3+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.1384, 1.1409])
    s1.line(l[0], [y4, y4+0.02], line_color='black')
    s1.line(l[1], [y4, y4+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                  y_units='data', text='Na', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.2682])
    s1.line(l[0], [y4, y4+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                  y_units='data', text='Na', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([2.2063, 2.2090])
    s1.line(l[0], [y4, y4+0.02], line_color='black')
    s1.line(l[1], [y4, y4+0.02], line_color='black')
    s1.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                  y_units='data', text='Na', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([2.2935, 2.3227, 2.3525, 2.3830, 2.4141])
    s1.line(l[0], [y3, y3+0.02], line_color='black')
    s1.line(l[1], [y3, y3+0.02], line_color='black')
    s1.line(l[2], [y3, y3+0.02], line_color='black')
    s1.line(l[3], [y3, y3+0.02], line_color='black')
    s1.line(l[4], [y3, y3+0.02], line_color='black')
    s1.line(l[[0, -1]], y+y3+0.02, line_color='black', line_width=1)
    s1.add_layout(Label(x=l[[0, -1]].mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='CO', render_mode='css',
                  text_font_size='8pt'))

    s1.xaxis.major_label_text_font_size = '0pt'
    s1.yaxis.major_label_text_font_size = '0pt'


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Order 2 (plots on the far left)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Contam plot (PA vs. wavelength)
    xlim0 = lamO2.min()
    xlim1 = lamO2.max()
    ylim0 = PA.min()-0.5*dPA
    ylim1 = PA.max()+0.5*dPA

    s5 = figure(tools=TOOLS, width=250, height=500, title=None,
                x_range=Range1d(xlim0, xlim1), y_range=s2.y_range)
    fig_data = np.log10(np.clip(contamO2.T, 1.e-10, 1.))[:, 300:]
    s5.image([fig_data], x=xlim0, y=ylim0, dw=xlim1-xlim0, dh=ylim1-ylim0,
             color_mapper=color_mapper)
    s5.yaxis.major_label_text_font_size = '0pt'
    s5.xaxis.axis_label = 'Wavelength (um)'

    # Line plot
    s6 = figure(tools=TOOLS, width=150, height=500, y_range=s2.y_range,
                x_range=Range1d(100, 0), title=None)
    s6.line(100*np.sum(contamO2 >= 0.001, axis=0)/ny, PA-dPA/2,
            line_color='blue', legend='> 0.001')
    s6.line(100*np.sum(contamO2 >= 0.01, axis=0)/ny, PA-dPA/2,
            line_color='green', legend='> 0.01')
    s6.xaxis.axis_label = '% channels contam.'
    s6.yaxis.major_label_text_font_size = '0pt'

    # Dummy plots for nice spacing
    s0 = figure(tools=TOOLS, width=150, plot_height=100, title=None)
    s0.outline_line_color = "white"
    s7 = figure(tools=TOOLS, width=150, plot_height=100, title=targetName)
    s7.outline_line_color = "white"

    # Add bad PAs
    for ybad0, ybad1 in badPA:
        s5.patch([xlim0, xlim1, xlim1, xlim0],
                 [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha)
        s6.patch([0, 100, 100, 0], [ybad1, ybad1, ybad0, ybad0],
                 color=bad_PA_color, alpha=bad_PA_alpha, legend='Bad PA')

    # Line list
    s4 = figure(tools=TOOLS, width=250, plot_height=100,
                x_range=s5.x_range, title=None)
    l = np.array([0.89, 0.99])
    s4.line(l, y+y1, line_color='black', line_width=1.5)
    s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([1.09, 1.2])
    s4.line(l, y+y1, line_color='black', line_width=1.5)
    s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([1.1, 1.24])
    s4.line(l, y+y2, line_color='black', line_width=1.5)
    s4.add_layout(Label(x=l.mean(), y=y2, x_units='data', y_units='data',
                  text='CH4', render_mode='css', text_font_size='8pt'))

    l = np.array([1.3, lamO2.max()])
    s4.line(l, y+y1, line_color='black', line_width=1.5)
    s4.add_layout(Label(x=l.mean(), y=y1, x_units='data', y_units='data',
                  text='H2O', render_mode='css', text_font_size='8pt'))

    l = np.array([0.7665, 0.7699])
    s4.line(l[0], [y3, y3+0.02], line_color='black')
    s4.line(l[1], [y3, y3+0.02], line_color='black')
    s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.1692, 1.1778])
    s4.line(l[0], [y3, y3+0.02], line_color='black')
    s4.line(l[1], [y3, y3+0.02], line_color='black')
    s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.2437, 1.2529])
    s4.line(l[0], [y3, y3+0.02], line_color='black')
    s4.line(l[1], [y3, y3+0.02], line_color='black')
    s4.add_layout(Label(x=l.mean(), y=y3+0.02, x_units='data',
                  y_units='data', text='K', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.1384, 1.1409])
    s4.line(l[0], [y4, y4+0.02], line_color='black')
    s4.line(l[1], [y4, y4+0.02], line_color='black')
    s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                  y_units='data', text='Na', render_mode='css',
                  text_font_size='8pt'))

    l = np.array([1.2682])
    s4.line(l[0], [y4, y4+0.02], line_color='black')
    s4.add_layout(Label(x=l.mean(), y=y4+0.02, x_units='data',
                  y_units='data', text='Na', render_mode='css',
                  text_font_size='8pt'))

    s4.xaxis.major_label_text_font_size = '0pt'
    s4.yaxis.major_label_text_font_size = '0pt'

    # put all the plots in a grid layout
    fig = gridplot(children=[[s7, s4, s1, s0], [s6, s5, s2, s3]])

return fig



































    title = "%s, NCam LWAR F322W2, V3PA =%3.0f$^\circ$" % (targetName, PA)

    dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
    dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)

    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_title(title)
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)

    ax.set_xlim([-500, 2300])
    ax.set_ylim([-500, 2300])
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    plt.gca().invert_xaxis()  # make X axis increase to left, along +V2 direction
    ax.set_xlabel('Detector pixel (x)')
    ax.set_ylabel('Detector pixel (y)')
    plt.text(2000, 2100, '$\Delta$K<%4.1f mag' % (dmag_limit), color='black', fontsize=10)

    #draw V2 & V3 axes
    ax.arrow(-250, -250, 0, 1500, head_width=75, head_length=125, fc='k', ec='k')
    plt.text(-300, 500, 'V3', color='black', fontsize=14)
    ax.arrow(-250, -250, 1500, 0, head_width=75, head_length=125, fc='k', ec='k')
    plt.text(500, -400, 'V2', color='black', fontsize=14)


    #draw N & E axes
    dx_N = 200 * math.sin(-PA * D2R)
    dy_N = 200 * math.cos(-PA * D2R)
    dx_E = 150 * math.sin((-PA+90) * D2R)
    dy_E = 150 * math.cos((-PA+90) * D2R)
    ax.arrow(1600, 1600, dx_N, dy_N, head_width=50, head_length=80, fc='r', ec='r')
    ax.arrow(1600, 1600, dx_E, dy_E, head_width=50, head_length=80, fc='r', ec='r')
    plt.text(1600+dx_N, 1600+dy_N, 'N', color='black', fontsize=14)
    plt.text(1600+dx_E, 1600+dy_E, 'E', color='black', fontsize=14)

    ax.plot(xbox,ybox) # detector bounds
    #ax.plot([-142, -132, 2254, 2198], [2822, -174, -180, 2800], color='r') # field area that can be imaged on detector
    ax.plot([xmin1_xmax1y, xmin1_miny, xxmax1_miny, xxmax1_xmax1y], [yxmax1_minx, ymin_minx, ymin_xmax1x, yxmax1_xmax1x], color='r') # field area that can be imaged on detector
    ax.plot(X_F322W2, Y_field, color='g', marker='o')
    ax.plot([xmin1_F322W2m1, xmax_F322W2], [Y_field, Y_field], linewidth=2.0, color='g')
    ax.text(xmax_F322W2-200, Y_field+50, 'F322W2 spectrum', color='g')
    #ax.plot(X_F322W2+dx, Y_field+dy, color='black', linestyle='none', marker='o')
    for i in range(nlimit):
        #ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=msize[i])
        ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=(K0-dmag[i])*3)
        if (col[i] == 1 or col[i] == 3):
            ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='red', linestyle='none', marker='o', markersize=msize[i])
            ax.plot([xmin1_F322W2m1+dx[i], xmax_F322W2+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='r')
        if (col[i] == 2 or col[i] == 3):
            ax.plot(X_F322W2+dx[i], Y_field+dy[i],  color='orange', linestyle='none', marker='o', markersize=msize[i])
            ax.plot([xmin1_F322W2m2+dx[i], xmax_F322W2m2+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='orange')
        col[i] = 0 # reset collision value for next filter
    title = "%s_F322W2_PA%03.0f" % (targetName, PA)
#     plt.ylim(-30,30)
    #fig.savefig(title)

#NIRCam source offsets
contaminationF322W2_order1 = np.zeros((nPA, nWaves))
contaminationF322W2_order2 = np.zeros((nPA, nWaves))

for kPA, PA in tqdm(enumerate(setPA), total=nPA):
    #NIRCam source offsets
    dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
    dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)
    col_F322W2= np.zeros(nlimit)

    # Check for F322W2 spectral collisitons
    for i in range(nlimit):
        # min extension for order 1 that the target will illuminate on the detector (unit: pixels)
        xmin1 = X_F322W2+dx[i] + (xmaxWL_F322W2 - wavelen_undev)/disp # m = 1 spectrum limits
        # max extension for order 1 that the target will illuminate on the detector (unit: pixels)
        xmax1 = xmin1 + dx_F322W2m1
        x     = X_F322W2 + dx[i]
        y     = Y_field + dy[i]
        # min extension for order 2 that the target will illuminate on the detector (unit: pixels)
        xmin2 = X_F322W2+dx[i] + (xmaxWL_F322W2 - wavelen_undev)*2.0 /disp # m = 2 spectrum limits
        # max extension for order 2 that the target will illuminate on the detector (unit: pixels)
        xmax2 = xmin1 + dx_F322W2m2

        # First check if 1st order spectrum companion (Background target) steps on 1st order spectrum of target
        # code loops over all targets and calculates dmag for all of them
        # if dmag is 0, that means i = source target.
        if dmag[i] != 0.0:
            #xbuff
            xmin_minusbuff_lt_xmin1 = xmin1 > (xmin1_F322W2m1-X_buff)
            xmax_plusbuff_gt_xmin1 = xmin1 < (xmax_F322W2+X_buff)
            xmin_minusbuff_lt_xmax1 = xmax1 > (xmin1_F322W2m1-X_buff)
            xmax_plusbuff_gt_xmax1 = xmax1 < (xmax_F322W2+X_buff)
            #ybuff
            yfield_minusbuff_lt_y = y > (Y_field-Y_buff)
            yfield_plusbuff_gt_y = y < (Y_field+Y_buff)

            # within x range that we care about (contamination)
            within_xrange1 = xmin_minusbuff_lt_xmin1 and xmax_plusbuff_gt_xmin1
            within_xrange2 = xmin_minusbuff_lt_xmax1 and xmax_plusbuff_gt_xmax1
            within_xrange = within_xrange1 or within_xrange2

            # within y range that we care about (contamination)
            within_yrange = yfield_minusbuff_lt_y and yfield_plusbuff_gt_y

            if within_xrange and within_yrange:
                # Next check to see that contaminating source is withn POM FOV and not the target itself:
                # if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xmax1_miny and math.fabs(dx[i]) > 3.0 and math.fabs(dy[i]) > 3.0): # don't flag the target itself
                col[i] = 1 # boolean 1 means contam exists
                           # boolean 0 means no contam

                # function call is gaussian1D, and (y_gs+Y_field)  is the array being computed over
                # gaussian1D(parameters)(the_array)
                # y_gs (gs = gaussian smooth) - an array ranging from -2048 to +2048
                # adding YField slides this array to be centered to the Y position of the background target
                bgGaussian = gaussian1D(center=y, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)
                contNow = sum(targetGaussian*bgGaussian)
                # the following for-loop adds the 1d gaussian on every column that has contamination
                for kx in range(int(round(xmin1)), int(round(xmax1+1))):
                    # for kx in range(int(round(xmin1_F322W2m1+dx[i])), int(round(xmax_F322W2))):
                    # for kx in range(int(round(xmin1_F322W2m1)), int(round(xmax_F322W2))):
                    if (kx >= 0) and (kx < xmax_F322W2):
                        contaminationF322W2_order1[kPA,kx] += contNow
                # print("m = 1: dx {0:.1f}, dy {1:.1f}, dK = {2:.1f}".format(dx[i], dy[i], dmag[i]))
            # Now check if 2nd order spectrum of companion steps on 1st order spectrum of target
             #xbuff
            xmin_minusbuff_lt_xmin2 = xmin2 > (xmin1_F322W2m2-X_buff)
            xmax_plusbuff_gt_xmin2 = xmin2 < (xmax_F322W2m2+X_buff)
            xmin_minusbuff_lt_xmax2 = xmax2 > (xmin1_F322W2m2-X_buff)
            xmax_plusbuff_gt_xmax2 = xmax2 < (xmax_F322W2m2+X_buff)
            #ybuff
            yfield_minusbuff_lt_y = y > (Y_field-Y_buff)
            yfield_plusbuff_gt_y = y < (Y_field+Y_buff)

            # within x range that we care about (contamination)
            within_xrange1 = xmin_minusbuff_lt_xmin2 and xmax_plusbuff_gt_xmin2
            within_xrange2 = xmin_minusbuff_lt_xmax2 and xmax_plusbuff_gt_xmax2
            within_xrange = within_xrange1 or within_xrange2

            # within y range that we care about (contamination)
            within_yrange = yfield_minusbuff_lt_y and yfield_plusbuff_gt_y


            # if (((xmin2 > (xmin1_F322W2m1-X_buff) and xmin2 < (xmax_F322W2+X_buff)) or (xmax2 > (xmin1_F322W2m1-X_buff) and (xmax2 < (xmax_F322W2+X_buff)))) and y > (Y_field-Y_buff) and y < (Y_field+Y_buff)):
                # if ((x+X_buff) > xmin1_xmax1y and (x-X_buff) < xmax1_miny and math.fabs(dx[i]) > 3.0 and math.fabs(dy[i]) > 3.0): # don't flag the target itself
            if within_xrange and within_yrange:
                # +2 specifies that there could be contamination from order 2 or 1 and 2
                # 0: no contamination
                # 1: contamination from order 1 only
                # 2: contamination from order 2 only
                # 3: contamination from order 1 and 2
                col[i] = col[i] + 2

                bgGaussian = gaussian1D(center=y, width=nyquistWidth, height=normalHeight, offset=0)(y_gs+Y_field)
                contNow = sum(targetGaussian*bgGaussian)
                for kx in range(int(round(xmin1)), int(round(xmax1+1))):
                    # for kx in range(int(round(xmin1_F322W2m1+dx[i])), int(round(xmax_F322W2))):
                    # for kx in range(int(round(xmin1_F322W2m1)), int(round(xmax_F322W2))):
                    if (kx >= xmin1_F322W2m2) and (kx < 2048):
                        contaminationF322W2_order2[kPA,kx] += contNow


plot_spectral_overlap_F322W2(col, dra, ddec, Decd0, 66, msize, dmag, targetName)

xbox = [-0, -0, 2040, 2040, -0]
ybox = [-0, 2040, 2040, -0, -0]

def plot_spectral_overlap_F444W(col, dra, ddec, Decd0, PA, msize, dmag, targetName):

    dx = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.cos((PA-0) * D2R) - ddec * 3600.0 / 0.065 * math.sin((PA-0) *D2R)
    dy = (dra * 3600.0 / 0.065) * math.cos(Decd0 * D2R) * math.sin((PA-0) * D2R) + ddec * 3600.0 / 0.065 * math.cos((PA-0) *D2R)

    title = "%s, NCam LWAR F444W, V3PA =%3.0f$^\circ$" % (targetName, PA)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_title(title)
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)

    ax.set_xlim([-500, 2300])
    ax.set_ylim([-500, 2300])
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    plt.gca().invert_xaxis()  # make X axis increase to left, along +V2 direction
    ax.set_xlabel('Detector pixel (x)')
    ax.set_ylabel('Detector pixel (y)')
    plt.text(2000, 2100, '$\Delta$K<%4.1f mag' % (dmag_limit), color='black', fontsize=10)

    #draw V2 & V3 axes
    ax.arrow(-250, -250, 0, 1500, head_width=75, head_length=125, fc='k', ec='k')
    plt.text(-300, 500, 'V3', color='black', fontsize=14)
    ax.arrow(-250, -250, 1500, 0, head_width=75, head_length=125, fc='k', ec='k')
    plt.text(500, -400, 'V2', color='black', fontsize=14)


    #draw N & E axes
    dx_N = 200 * math.sin(-PA * D2R)
    dy_N = 200 * math.cos(-PA * D2R)
    dx_E = 150 * math.sin((-PA+90) * D2R)
    dy_E = 150 * math.cos((-PA+90) * D2R)
    ax.arrow(1600, 1600, dx_N, dy_N, head_width=50, head_length=80, fc='r', ec='r')
    ax.arrow(1600, 1600, dx_E, dy_E, head_width=50, head_length=80, fc='r', ec='r')
    plt.text(1600+dx_N, 1600+dy_N, 'N', color='black', fontsize=14)
    plt.text(1600+dx_E, 1600+dy_E, 'E', color='black', fontsize=14)

    ax.plot(xbox,ybox) # detector bounds
    #ax.plot([-142, -132, 2254, 2198], [2822, -174, -180, 2800], color='r') # field area that can be imaged on detector
    ax.plot([xmin1_xmax1y, xmin1_miny, xxmax1_miny, xxmax1_xmax1y], [yxmax1_minx, ymin_minx, ymin_xmax1x, yxmax1_xmax1x], color='r') # field area that can be imaged on detector
    ax.plot(X_F444W, Y_field, color='g', marker='o')
    ax.plot([xmin1_F444W, xxmax1_F444W], [Y_field, Y_field], linewidth=2.0, color='g')
    ax.text(xxmax1_F444W-200, Y_field+50, 'F444W spectrum', color='g')
    for i in range(nlimit):
        ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=msize[i])
        #ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='black', linestyle='none', marker='o', markersize=(K0-dmag[i])*3)
        if (col[i] == 1):
            ax.plot(X_F444W+dx[i], Y_field+dy[i],  color='red', linestyle='none', marker='o', markersize=10)
            ax.plot([xmin1_F444W+dx[i], xxmax1_F444W+dx[i]], [Y_field+dy[i], Y_field+dy[i]], linewidth=2.0, color='r')
        col[i] = 0 # reset collision value for next filter

    title = "%s_F444W_PA%03.0f" % (targetName, PA)
    #fig.savefig(title)
