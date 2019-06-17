from cosmosis.datablock import names, option_section
from pandas import read_csv
from numpy import log, genfromtxt
from scipy.optimize import curve_fit
import os


# Define some names for identifying files and column names
ini_file = "spectral_index"
spec_file = "wavenumber_spectrum"
twopf_name = "twopf_observables"
tab_names = ["k_phys", "A_s", "A_t"]


# Define a quadratic function to use with curve_fit
def quadratic(x, a, b, c):
    return a*x**2 + b*x + c


# setup function for CosmoSIS - pull in k_pivot & empty list for n_s storage
def setup(options):
    config = {
        'k_pivot': options.get_double(option_section, 'k_pivot'),
    }
    return config


# execute function for CosmoSIS - runs once for each sample
def execute(block, config):
    k_pivot = config['k_pivot']
    # Load in spec_index table from datablock onto a Pandas df
    file = block[spec_file, "spec_index_table"]

    # Pandas implementation of read-in
    df = read_csv(file, sep="\t", header=None, names=tab_names)
    # NumPy implmentation of read-in
    # df = genfromtxt(file, delimiter="\t", names=tab_names)

    # Make ln versions of k_pivot, k, A_s & A_t
    lnKpiv = log(k_pivot)
    lnK = log(df["k_phys"])
    lnAs = log(df["A_s"])
    lnAt = log(df["A_t"])

    # Use curve_fit to find optimised parameters & covariance matrices for quadratic fit
    s_par, s_cov = curve_fit(quadratic, lnK, lnAs)
    t_par, t_cov = curve_fit(quadratic, lnK, lnAt)

    # Use known form for gradient of quadratic at k_pivot
    deriv_s = 2.0*s_par[0]*lnKpiv + s_par[1]
    deriv_t = 2.0*t_par[0]*lnKpiv + t_par[1]

    # Spectral indices are ns=1+deriv_s, nt=deriv_t at k_pivot
    n_s = 1 + deriv_s
    n_t = deriv_t

    # Return the n_s and n_t values to the datablock
    block[twopf_name, "n_s_scipy"] = n_s
    block[twopf_name, "n_t_scipy"] = n_t

    # Delete the temporary data file
    os.remove(file)

    # We tell CosmoSIS that everything went fine by returning zero
    return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
