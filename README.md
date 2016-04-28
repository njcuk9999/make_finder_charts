# make_finder_charts

    list_charts():
    """
    Prints user readable form of the surveys available
    :return:
    """

    make_chart(ra=None, dec=None, p=None, name='', pmra=None, pmde=None,
               savename='./finder', show=False, **kwargs)
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
