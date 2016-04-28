"""
This program creates finder charts input require is a catlogue with RA and Dec
in (currently in fits format) in addition you can add
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
import aplpy
import os

# =============================================================================
# Define variables
# =============================================================================
# define location of fits file input (if using main program)
workspace = '/local/home/ncook/Projects/Project8/Stage1'
catpath = workspace + '/Data/Catalogues/NJCM_full_e7_forcm.fits'
psave = workspace + '/Plots/finder_charts'
# -----------------------------------------------------------------------------
# define catalogue input columns (None where possible to skip)
idcol = 'WISE_id'
racol = 'ra'
deccol = 'dec'
pmracol = 'pmra'
pmdecol = 'pmde'
pmunit = u.arcsec * u.yr ** -1  # units of pm in catalogue
pmtime = 1000.0  # scale proper motions up (by this many times)
# -----------------------------------------------------------------------------
# define survey (can find info in SkyView.list_surveys())
# SDSSg, SDSSi, SDSSr, SDSSu, SDSSz, 2MASS-J, 2MASS-H, 2MASS-K, UKIDSS-Y,
# UKIDSS-J, UKIDSS-H, UKIDSS-K, WISE 3.4, WISE 4.6, WISE 12, WISE 22
# this is from red --> blue
survey = ['2MASS-K', '2MASS-H', '2MASS-J']
# survey = ['WISE 12', 'WISE 4.6', 'WISE 3.4']
# survey = ['SDSSg', 'SDSSi', 'SDSSz']
# set coord system "J2000", "B1950", "Galactic","E2000", "ICRS"
coords = 'J2000'
# size of image must have units
size_of_image = (300.0 / 3600) * u.deg
# pixels in x and y (will be square)
res_of_image = 1024
# -----------------------------------------------------------------------------
# scale to add
scale = (60.0 / 3600) * u.deg  # must be an astropy quantity
scalelabel = '1 arcmin'
# -----------------------------------------------------------------------------
# plot rgb and if True set stretch
rgb = True
# for rgb define the stretch in red green and blue colours
stretch = ['log', 'log', 'log']
# if rgb is False set cmap for single colour (i.e. 'gist_heat', 'gist_gray')
cmap = 'gist_gray'
# -----------------------------------------------------------------------------
# if true show instead of saving
pshow = True
# colour of target
ctarget = 'r'
# marker size of target
s_target = 300
# colour of compass/scale
cother = 'y'
# size of figure
psize = (10, 10)
# colour of grid
cgrid = '0.75'
# linestyle of grid
lgrid = '--'


# =============================================================================
# Define functions
# =============================================================================
def make_chart(ra=None, dec=None, p=None, name='', pmra=None, pmde=None,
               savename='./finder', show=False, **kwargs):
    """
    Plots the finder chart
    :param ra: array or None, right ascension (floats)
    :param dec: array or None, declination (floats)
    :param p: array or None, SIMBAD or NED name (for resolver)
    :param name: array, names or identifiers (strings)
    :param pmra: array, proper motions in right ascension (floats)
    :param pmde: array, proper motions in declinations (floats)
    :param savename: string, save location of plot
    :param show: bool, if True shows the graph instead of saving to file
    :param kwargs: see below

    kwargs are as follows:

    - coords
                Choose among common equatorial, galactic and ecliptic
                coordinate systems (``"J2000"``, ``"B1950"``, ``"Galactic"``,
                ``"E2000"``, ``"ICRS"``) or pass a custom string.
    - survey: string, surveys used in SkyView.get_images

              define survey (can find info in SkyView.list_surveys())

              e.g.:

                    SDSSg, SDSSi, SDSSr, SDSSu, SDSSz, 2MASS-J, 2MASS-H,
                    2MASS-K, UKIDSS-Y, UKIDSS-J, UKIDSS-H, UKIDSS-K,
                    WISE 3.4, WISE 4.6, WISE 12, WISE 22

              if colour kwarg is True this is a list of three strings
              and is from red --> blue i.e. ['WISE 12', 'WISE 4.6', 'WISE 3.4']

    - radius: float, half width of image
    - pixels: int, size of the image produced in pixels (for x and y)
    - colour: bool, whether to plot in RGB mode or not
    - cmap
                string, matplotlib colormap indicator
                    (used if not rgb), default='gist_gray'
    - stretches
                tuple, defines the rgb stretches to use on rgb plot
                 [stretch in r, stretch in g, stretch in b]
                 default is ['log', 'log', 'log'] can have:
                'linear', 'log', 'sqrt', 'arcsinh', 'power'

    - size: tuple, (x size, y size) in inches same as:
                    fig = plt.figure()
                    fig.set_size_inches(*size)
    - ct: string, colour of the target marker and pm arrows
    - co: string, colour of the compass and scale markers
    - sc: astropy quantity, size of scale to add
    - scl: string, label for scale to add
    - gridc: string, colour of the grid
    - gridl: string, linestyle of the grid

    - pmt: factor proper motions are scaled by

    :return:
    """
    # load params from kwargs
    imres = kwargs.get('pixels', 1024)
    size = kwargs.get('size', (10, 10))
    ct = kwargs.get('ct', 'r')
    co = kwargs.get('co', 'y')
    sc = kwargs.get('sc', 1 * u.deg)
    scl = kwargs.get('scl', '1 deg')
    markersize = kwargs.get('ms', 300)
    pmt = kwargs.get('pmt', 1.0)
    # -------------------------------------------------------------------------
    print ('\n Getting image...')
    gc = plot_image(ra=ra, dec=dec, pos=p, **kwargs)
    # -------------------------------------------------------------------------
    print ('\n Adding markers...')
    gc.show_markers(ra, dec, edgecolor=ct, s=markersize)
    # -------------------------------------------------------------------------
    print ('\n Adding scale bar and compass...')
    gc = add_scale_bar(gc, sc, scl, c=co)
    gc = add_N_E(gc, imres, c=co)
    # -------------------------------------------------------------------------
    # add proper motions if not None
    if pmra is not None and pmde is not None:
        print ('\n Adding proper motion arrows...')
        # plot pm arrows
        gc.show_arrows(ra, dec, pmra.value, pmde.value, color=ct)
    # -------------------------------------------------------------------------
    # add title (with ID) if not None
    plt.title('Target: {0} \npm scaled by factor x{1}'.format(name, pmt))
    # -------------------------------------------------------------------------
    if show:
        plt.show()
        plt.close()
    else:
        gc._figure.set_size_inches(*size)
        plt.savefig(savename + '.png', bbox_inches='tight')
        plt.close()
    return gc


def list_charts():
    """
    Prints user readable form of the surveys available
    :return:
    """
    print ('{0}\n\tList of Surveys\n{0}\n\n'.format('=' * 50))
    d = SkyView.survey_dict
    N = len(d.keys())
    for k, key in enumerate(d.keys()):
        print ('{0}\n\t{1} of {2}\n\t{3}\n{0}'.format('=' * 50, k + 1, N, key))
        for it, entry in enumerate(d[key]):
            if "None" in entry[:4]:
                continue
            print ('\t{0}.\t{1}'.format(it + 1, entry))
    print ('=' * 50 + '\n\n')


def plot_image(ra=None, dec=None, pos=None, **kwargs):
    """
    Use SkyView to get a cutout image

    :param ra: float or string or None, right ascension (deg or hh:mm:ss)
    :param dec: float or string or None, declination (deg or dd:mm:ss)
    :param pos: string or None,
                Determines the center of the field to be retrieved. Both
                coordinates (also equatorial ones) and object names are
                supported. Object names are converted to coordinates via the
                SIMBAD or NED name resolver. See the reference for more info
                on the supported syntax for coordinates.
    :param kwargs: see below

    kwargs are as follows
    - coords
                Choose among common equatorial, galactic and ecliptic
                coordinate systems ("J2000", "B1950", "Galactic",
                "E2000", "ICRS") or pass a custom string.
    - survey string, surveys used in SkyView.get_images

              define survey (can find info in SkyView.list_surveys())

              e.g.:

                    SDSSg, SDSSi, SDSSr, SDSSu, SDSSz, 2MASS-J, 2MASS-H,
                    2MASS-K, UKIDSS-Y, UKIDSS-J, UKIDSS-H, UKIDSS-K,
                    WISE 3.4, WISE 4.6, WISE 12, WISE 22

              if colour kwarg is True this is a list of three strings
              and is from red --> blue i.e. ['WISE 12', 'WISE 4.6', 'WISE 3.4']

    - radius
                float, half width of image
    - pixels
                int, size of the image produced in pixels (for x and y)
    - colour
                bool, whether to plot rgb or not, default=False
    - cmap
                string, matplotlib colormap indicator
                    (used if not rgb), default='gist_gray'
    - stretches
                tuple, defines the rgb stretches to use on rgb plot
                 [stretch in r, stretch in g, stretch in b]
                 default is ['log', 'log', 'log'] can have:
                'linear', 'log', 'sqrt', 'arcsinh', 'power'
    - gridc
                string, colour of the grid
    - gridl
                string, linestyle of the grid

    :return:
    """
    # load params from kwargs
    coords = kwargs.get('coords', 'J2000')
    s = kwargs.get('survey', '2MASS-J')
    radius = kwargs.get('radius', 1.0 * u.deg)
    pixels = kwargs.get('pixels', 1024)
    colour = kwargs.get('colour', False)
    cmap = kwargs.get('cmap', 'gist_gray')
    stretches = kwargs.get('stretches', ('log', 'log', 'log'))
    gridc = kwargs.get('gridc', '0.75')
    gridl = kwargs.get('gridl', '--')
    # get position as ra dec string
    if ra is not None and dec is not None:
        pos = '{0} {1}'.format(ra, dec)
    elif pos is not None:
        pos = pos
    else:
        raise Exception('Error need RA, Dec or Pos to get image')
    # download images
    giargs = dict(position=pos, survey=s, radius=radius, pixels=pixels,
                  coordinates=coords)
    images = SkyView.get_images(**giargs)
    # plot aplpy FITS figure
    apx = aplpy.FITSFigure(images[0][0])
    # if colour expect 3 surveys and make an rgb plot
    if colour:
        files = []
        # loops round each image
        for i, j in enumerate(images):
            # defines a tmp filename for each iamge
            fname = '{0}.fits'.format(i)
            # writes fits to file
            j[0].writeto(fname, clobber=True)
            # appends filename to file list
            files.append(fname)
        # add rgb image to file list
        files.append('rgb.png')
        # create rgb image
        rgbargs = dict(zip(['stretch_r', 'stretch_g', 'stretch_b'], stretches))
        aplpy.make_rgb_image(files[:3], files[3], **rgbargs)
        # show rgb plot on top
        apx.show_rgb('rgb.png')
        for fn in files:
            os.remove(fn)
    # if not colour just show the colour scale
    else:
        apx.show_colorscale(cmap=cmap)
    # add grid
    apx.add_grid()
    apx.grid.set_color(gridc)
    apx.grid.set_linestyle(gridl)
    # return aplpy figure
    return apx


def add_scale_bar(apx, s, sl, c='y'):
    """
    Adds an aplpy scale bar
    :param apx: aplpy fig
    :param s: astropy.units quantity,  size (with units) of the scale bar
    :param sl: string, label for scale bar
    :param c: string, colour of scalebar
    :return:
    """
    apx.add_scalebar(s)
    apx.scalebar.set_label(sl)
    apx.scalebar.set_color(c)
    return apx


def add_N_E(apx, res, posx=0.8, posy=0.8, lenx=0.1, leny=0.1, c='y'):
    """
    Adds a compass to the plot

    :param apx: aplpy fig
    :param res: resolution of the graph (i.e. how many pixels in x and y)
    :param posx: position of compass relative to x axis (0 --> 1)
    :param posy: position of compass relative to y axis (0 --> 1)
    :param lenx: distances away from arrow in x (in axis units 0 --> 1)
    :param leny: distances away from arrow in y (in axis units 0 --> 1)
    :param c: string, colour of compass (text and arrows)
    :return:
    """
    px, py = apx.pixel2world(posx * res, posy * res)
    px1, py1 = apx.pixel2world((posx - lenx) * res, posy * res)
    px2, py2 = apx.pixel2world(posx * res, (posy + leny) * res)
    xlen, ylen = abs(px2 - px1), abs(py2 - py1)
    apx.show_arrows([px, px], [py, py], [xlen, 0.0], [0.0, ylen], color=c)
    apx.add_label(px, py + ylen * 1.1, 'N', color=c)
    apx.add_label(px + xlen * 1.1, py, 'E', color=c)
    return apx


def get_data(path, idc, rac, decc, pmrac, pmdec,
             pmu=u.arcsec * u.yr ** -1, pmtime=1000.0):
    """
    Get data and extract column information
    :param path: str, path to fits file
    :param idc: string or None, key for name or identifier column
    :param rac: string, key for right ascension column
    :param decc: string, key for Declination column
    :param pmrac: string or None, key for proper motion in right ascension
    :param pmdec: string or None, key for proper motion in declination
    :param pmu: astropy quantity, units of raw proper motion
                default is arcsec yr^-1
    :param pmtime: float, number of years to scale proper motion by
                   (1000 default)
    :return:
    """
    print ('\n Loading catalogue...')
    data = fits.getdata(path, ext=1)
    ra, dec = np.array(data[rac]), np.array(data[decc])
    # get the ids if not None
    if idc is not None:
        name = np.array(data[idc])
    else:
        name = np.zeros(len(ra), dtype=str)
    # get the pm info if not None
    if pmdec is not None and pmrac is not None:
        pmra, pmde = np.array(data[pmrac]) * pmu, np.array(data[pmdec]) * pmu
        # convert to degress per 1000 yrs
        uu = u.deg * u.yr ** -1
        pmra, pmde = pmtime * pmra.to(uu), pmtime * pmde.to(uu)
    else:
        pmra, pmde = [None] * len(ra), [None] * len(ra)
    # return
    return ra, dec, name, pmra, pmde


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # get data
    # odata = get_data(catpath, idcol, racol, deccol, pmracol, pmdecol,
    #                  pmunit, pmtime)
    # ras, decs, ids, pmras, pmdes = odata
    etacarinae = SkyCoord.from_name('Eta Carinae').icrs


    ras = [85.24583333, etacarinae.ra.value]
    decs = [-2.4583333, etacarinae.dec.value]
    ids = ['Horse Head Nebula', 'Eta Carinae']
    pmras = [100*pmunit, 100*pmunit]
    pmdes = [100*pmunit, 100*pmunit]
    # -------------------------------------------------------------------------
    # loop round each target
    for row in range(len(ids)):
        # ---------------------------------------------------------------------
        # print progress
        if idcol is not None:
            tstring = 'Target: {0}'.format(ids[row])
        else:
            tstring = 'Row: {0}'.format(row + 1)
        print ('\n{0}\n{1}\n{0}\n'.format('=' * 50, tstring))
        # ---------------------------------------------------------------------
        # plot finder chart
        plotkwargs = dict(ra=ras[row], dec=decs[row], name=ids[row],
                          pmra=pmras[row], pmde=pmdes[row], savename=psave,
                          show=pshow, survey=survey, radius=size_of_image,
                          colour=rgb, size=psize, pixels=res_of_image,
                          ct=ctarget, co=cother, sc=scale, scl=scalelabel,
                          gridc=cgrid, gridl=lgrid, stretches=stretch,
                          cmap=cmap, coord=coords, ms=s_target, pmt=pmtime)
        g = make_chart(**plotkwargs)

# =============================================================================
# End of code
# =============================================================================
