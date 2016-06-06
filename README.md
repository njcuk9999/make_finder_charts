# make_finder_charts

```python
list_charts()
```
Prints user readable form of the surveys available


```python
make_chart(ra=None, dec=None, p=None, name='', pmra=None, pmde=None, 
           savename='./finder', show=False, **kwargs)
```

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

### Example code

```python
# define location of fits file input (if using main program)
workspace = './'
catpath = workspace + 'catalogue.fits'
psave = workspace + '/Plots/'
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
```
