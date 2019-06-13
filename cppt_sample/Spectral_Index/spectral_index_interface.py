from cosmosis.datablock import names, option_section
import pandas as pd
from numpy import arange, linspace, log
# from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from builtins import str
import uuid

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
        'n_s_storage': []
    }
    return config


# Define matplotlibs subplot for n_s histogram
fig, ax = plt.subplots()


# execute function for CosmoSIS - runs once for each sample
def execute(block, config):
    k_pivot = config['k_pivot']
    # Load in spec_index table from datablock onto a Pandas df
    file = block[spec_file, "spec_index_table"]
    df = pd.read_csv(file, sep="\t", header=None, names=tab_names)

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

    # Spectrum values at k_pivot
    val_s = quadratic(lnKpiv, s_par[0], s_par[1], s_par[2])
    val_t = quadratic(lnKpiv, t_par[0], t_par[1], t_par[2])

    # # Fit a 2nd order spline to the ln(A_s) vs. ln(k) data
    # spline_func_s = splrep(df["k_phys"], df["A_s"], k=2)
    # # Evaluate function @ k_piv
    # val_s = splev(k_pivot, spline_func_s)
    # # Evaluate derivative dln(P_s)/dln(k) @ k=k_piv
    # deriv_s = splev(k_pivot, spline_func_s, der=1)
    # # n_s = deriv_s + 1
    # n_s = 1 + deriv_s * (k_pivot/val_s)

    # # Fit a 2nd order spline to the ln(A_t) vs. ln(k) data
    # spline_func_t = splrep(df["k_phys"], df["A_t"], k=2)
    # # Evaluate function @ k_piv
    # val_t = splev(k_pivot, spline_func_t)
    # # Evaluate derivative n_t = dln(P_t)/dln(k) @ k=k_piv
    # deriv_t = splev(k_pivot, spline_func_t, der=1)
    # n_t = deriv_t * (k_pivot/val_t)

    # create random filename for plots
    filename = uuid.uuid4().hex
    # Folder name
    folder = "Spec_Plots/"

    # Dense region for plotting splines
    k_vals = log(linspace(df["k_phys"].min(), df["k_phys"].max(), 1000))
    k_tang_vals = log(linspace(1.997e-3, 2.003e-3, 20))

    # Now make a plot of the scalar data and fit
    plt.figure()
    plt.plot(lnK, lnAs, 'x', label="data")
    plt.plot(k_vals, quadratic(k_vals, s_par[0], s_par[1], s_par[2]), alpha=0.5, label="quad fit")
    # plt.plot(k_vals, splev(k_vals, spline_func_s), alpha=0.5, label="spline")
    tangent_s = (val_s + deriv_s * (k_tang_vals - lnKpiv))
    plt.plot(k_tang_vals, tangent_s, '--r', label="n_s tangent: " + str(n_s))
    plt.xlabel("ln(k)")
    plt.ylabel("ln(A_s)")
    # plt.xscale('log')
    # plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='best')
    plt.savefig(folder+"scalar"+filename+".pdf")
    plt.close()

    # Now make a plot of the tensor data and fit
    plt.figure()
    plt.plot(lnK, lnAt, 'x', label="data")
    plt.plot(k_vals, quadratic(k_vals, t_par[0], t_par[1], t_par[2]), alpha=0.5, label="quad fit")
    # plt.plot(k_vals, splev(k_vals, spline_func_t), alpha=0.5, label="spline")
    tangent_t = val_t + deriv_t * (k_tang_vals - lnKpiv)
    plt.plot(k_tang_vals, tangent_t, '--r', label="n_t tangent: " + str(n_t))
    plt.xlabel("ln(k)")
    plt.ylabel("ln(A_t)")
    # plt.xscale('log')
    # plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='best')
    plt.savefig(folder+"tensor_"+filename+".pdf")
    plt.close()

    print n_s
    print n_t

    config['n_s_storage'].append(n_s)

    ax.hist(config['n_s_storage'], bins=arange(0.95, 0.98, 1e-4))
    ax.set_xlim([0.95, 0.98])
    fig.savefig(folder+"n_s_histogram.pdf")

    # Return the n_s and n_t values to the datablock
    block[twopf_name, "n_s_scipy"] = n_s
    block[twopf_name, "n_t_scipy"] = n_t

    # We tell CosmoSIS that everything went fine by returning zero
    return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
