import numpy as np
from scipy.stats import gamma


def centers(x):
    """
    Compute the centers of the bins based on the input edges.

    Parameters
    ----------
    x : ndarray
        Array of bin edges.

    Returns
    -------
    ndarray
        Array of bin centers.
    """
    return (x[:-1] + x[1:]) * 0.5


def edges(x):
    """
    Compute bin edges with centers that approximately match the input array.

    Returns bin edges with centers that approximately match x. Approximate
    "inverse" of center(x). Note that it is impossible to ensure that
    centers(edges(x)) == x for all x using the functions defined in
    this module as the returned values of centers(e) are subject to
    constraints that do not necessarily exist for arbitrary center
    points.

    Parameters
    ----------
    x : ndarray
        Array of bin centers.

    Returns
    -------
    ndarray
        Array of bin edges.
    """
    c = centers(x)
    return np.concatenate(([2 * x[0] - c[0]], c, [2 * x[-1] - c[-1]]))


def calc_nbins(x):
    """
    Calculate the number of bins for the input data using the Freedman-Diaconis rule.

    Parameters
    ----------
    x : ndarray
        Input data.

    Returns
    -------
    int
        Number of bins.
    """
    n = (np.max(x) - np.min(x)) / (2 * len(x)**(-1. / 3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return 10 if np.isnan(n) else int(n)


def calc_bins(x):
    """
    Calculate the bin edges for the input data.

    Parameters
    ----------
    x : ndarray
        Input data.

    Returns
    -------
    ndarray
        Array of bin edges.
    """
    nbins = calc_nbins(x)
    return np.linspace(np.min(x), np.max(x) + 2, num=nbins + 1)


def cart_to_sphe(x, y, z):
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    x, y, z : float
        Cartesian coordinates.

    Returns
    -------
    tuple
        Spherical coordinates (r, theta, phi).
    """
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return r, theta, phi


def sphe_to_cart(r, theta, phi):
    """
    Convert spherical coordinates to Cartesian coordinates.

    Parameters
    ----------
    r : float
        Radius.
    theta : float
        Polar angle (in radians).
    phi : float
        Azimuthal angle (in radians).

    Returns
    -------
    tuple
        Cartesian coordinates (x, y, z).
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def sphe_to_kent(theta, phi):
    """
    Convert spherical coordinates to the Kent (1982) coordinate system.

    Shift x,y,z->x2,x3,x1 which is the coordinate system used in Kent 1982.
    Returns np.array((x1,x2,x3), ...), a format readable by the fb8 package.


    Parameters
    ----------
    theta : float
        Polar angle (in radians).
    phi : float
        Azimuthal angle (in radians).

    Returns
    -------
    ndarray
        Array of coordinates in Kent system.
    """
    xyz = sphe_to_cart(1, theta, phi)
    xs = np.asarray((xyz[2], xyz[0], xyz[1])).T
    return xs


def center_angle(theta0, phi0, theta1, phi1):
    """
    Calculate the angle between two directions in spherical coordinates.

    Parameters
    ----------
    theta0 : float
        Polar angle (in radians) of the first direction.
    phi0 : float
        Azimuthal angle (in radians) of the first direction.
    theta1 : float
        Polar angle (in radians) of the second direction.
    phi1 : float
        Azimuthal angle (in radians) of the second direction.

    Returns
    -------
    float
        Angle between the two directions (in radians).
    """
    pos1 = sphe_to_cart(1, theta0, phi0)
    pos2 = sphe_to_cart(1, theta1, phi1)
    dotp = np.sum(np.asarray(pos1).T * np.asarray(pos2).T, axis=-1)
    return np.asfarray(np.arccos(np.clip(dotp, -1, 1)))


def vmf_stats(thetas, phis, p=3):
    """
    Calculate statistics for data drawn from a Von-Mises Fisher distribution.

    Return the average of directional coordinates assuming they are
    drawn from a Von-Mises Fisher distribution (gaussian on a sphere)

    Parameters
    ----------
    thetas : ndarray
        Array of polar angles (in radians).
    phis : ndarray
        Array of azimuthal angles (in radians).
    p : int, optional
        Dimensionality of the dataset, by default 3.

    Returns
    -------
    tuple
        Mean direction (theta, phi), mean resultant length (R), 
        concentration parameter (kappa), and standard deviation (sigma).
    """
    norm, theta, phi = mean_ang(thetas, phis)

    R = norm / thetas.size
    # below are approximations
    kappa = R * (p - R**2) / (1 - R**2)
    d = 1 - np.sum(np.dot(np.array(sphe_to_cart(1, theta, phi)), 
                          np.array(sphe_to_cart(1, thetas, phis)))**2) / thetas.size
    sigma = np.sqrt(d / (thetas.size * R**2))
    return theta, phi, R, kappa, sigma


def mean_ang(thetas, phis):
    """
    Calculate the mean direction for given angles.

    Parameters
    ----------
    thetas : ndarray
        Array of polar angles (in radians).
    phis : ndarray
        Array of azimuthal angles (in radians).

    Returns
    -------
    tuple
        Norm, mean polar angle (in radians), and mean azimuthal angle (in radians).
    """
    xsum = np.sum(np.sin(thetas) * np.cos(phis))
    ysum = np.sum(np.sin(thetas) * np.sin(phis))
    zsum = np.sum(np.cos(thetas))

    norm, mean_th, mean_phi = cart_to_sphe(xsum, ysum, zsum)
    return norm, mean_th, mean_phi


def mode_ang(thetas, phis):
    """
    Calculate the mode direction for given angles.

    Parameters
    ----------
    thetas : ndarray
        Array of polar angles (in radians).
    phis : ndarray
        Array of azimuthal angles (in radians).

    Returns
    -------
    tuple
        Norm, mode polar angle (in radians), and mode azimuthal angle (in radians).
    """
    import healpy as hp
    if len(thetas) == 0:
        return np.nan, np.nan, np.nan
    nside = 128
    pixs = hp.ang2pix(nside, thetas, phis)
    counts = np.bincount(pixs)
    mode_pixs = np.flatnonzero(counts == counts.max())
    while nside > 16 and len(mode_pixs) > 1:
        nside /= 2
        pixs = hp.ang2pix(nside, thetas, phis)
        counts = np.bincount(pixs)
        mode_pixs = np.flatnonzero(counts == counts.max())

    norm, mode_th, mode_phi = mean_ang(*hp.pix2ang(nside, mode_pixs))
    return norm, mode_th, mode_phi


def med_ang_res(theta, phi, thetas, phis):
    """
    Calculate the median angular deviation from a reference direction.

    Parameters
    ----------
    theta : float
        Polar angle (in radians) of the reference direction.
    phi : float
        Azimuthal angle (in radians) of the reference direction.
    thetas : ndarray
        Array of polar angles (in radians) of the data points.
    phis : ndarray
        Array of azimuthal angles (in radians) of the data points.

    Returns
    -------
    float
        Median angular deviation (in radians).
    """
    return np.median(center_angle(theta, phi, thetas, phis))


def poisson_llh(hdata, hexp):
    """
    Calculate the Poisson log-likelihood for given data and expectations.

    Parameters
    ----------
    hdata : ndarray
        Observed data.
    hexp : ndarray
        Expected values.

    Returns
    -------
    tuple
        Negative log-likelihood and count of valid entries.
    """
    pllh = np.ma.masked_where(hexp == 0, gamma.logpdf(hexp, hdata + 1))
    return -np.sum(pllh), pllh.count()


def dchi2(hdata, hexp):
    """
    Calculate the unnormalized chi-squared statistic for given data and expectations.

    Parameters
    ----------
    hdata : ndarray
        Observed data.
    hexp : ndarray
        Expected values.

    Returns
    -------
    tuple
        Chi-squared statistic and number of data points.
    """
    return np.sum((hdata - hexp)**2), len(hdata)


def dima_llh(hdata, hexp, sigma=0.1, ns=1, nd=1):
    """
    Calculate Dima's log-likelihood for given data and expectations.

    Parameters
    ----------
    hdata : ndarray
        Observed data.
    hexp : ndarray
        Expected values.
    sigma : float, optional
        Standard deviation, by default 0.1.
    ns : float, optional
        Scaling factor for expected values, by default 1.
    nd : float, optional
        Scaling factor for observed data, by default 1.

    Returns
    -------
    tuple
        Log-likelihood and count of valid entries.
    """
    s = hexp * ns
    d = hdata * nd
    ls = np.ma.masked_where(s + d == 0, (s + d) / (ns + nd))
    ld = ls
    tol = 1e-8

    if sigma > 0:
        while True:
            f = ns * ls - s - np.log(ld / ls) / sigma**2
            df = 1 + (1 / (nd * ld) + 1 / (ns * ls)) / sigma**2
            dls = -f / (ns * df)
            ls += dls
            ld = (s + d - ns * ls) / nd
            if np.all(abs(dls) < tol * ls):
                break

    dllh = s * np.log(hexp / ls) + d * np.log(hdata / ld) + np.log(ld / ls)**2 / (2 * sigma**2)
    return np.sum(dllh), dllh.count()


def most_likely(arr, weights=None):
    """
    Find the densest region in a 1D array of data.

    Parameters
    ----------
    arr : ndarray
        Input data array.
    weights : ndarray, optional
        Weights for the data, by default None.

    Returns
    -------
    float
        Center of the densest region.
    """
    binning = calc_bins(arr)
    harr = np.histogram(arr, binning, weights=weights)[0]
    return centers(binning)[np.argmax(harr)]


def interval(arr, percentile=68.3):
    """
    Calculate the shortest interval around the mode for a given percentile.

    Parameters
    ----------
    arr : ndarray
        Input data array.
    percentile : float, optional
        Percentile for the interval, by default 68.3.

    Returns
    -------
    tuple
        Lower bound, mode, and upper bound of the interval.
    """
    if len(arr) == 0:
        return np.nan, np.nan, np.nan
    center = most_likely(arr)
    sarr = sorted(arr)
    delta = np.abs(sarr - center)
    curr_low = np.argmin(delta)
    curr_up = curr_low
    npoints = len(sarr)
    while curr_up - curr_low < percentile / 100. * npoints:
        if curr_low == 0:
            curr_up += 1
        elif curr_up == npoints - 1:
            curr_low -= 1
        elif sarr[curr_up] - sarr[curr_low - 1] < sarr[curr_up + 1] - sarr[curr_low]:
            curr_low -= 1
        elif sarr[curr_up] - sarr[curr_low - 1] > sarr[curr_up + 1] - sarr[curr_low]:
            curr_up += 1
        elif (curr_up - curr_low) % 2:
            # They are equal, step half of the time up and down
            curr_low -= 1
        else:
            curr_up += 1

    return sarr[curr_low], center, sarr[curr_up]
