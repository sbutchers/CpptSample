from cosmosis.datablock import names, option_section
import pandas as pd
import numpy as np
from scipy.interpolate import splev, splrep
import matplotlib.pyplot as plt
from builtins import str
import uuid

ini_file = "spectral_index"
spec_file = "wavenumber_spectrum"
twopf_name = "twopf_observables"
tab_names = ["k_phys", "A_s", "A_t"]


def setup(options):
    config = {
        'k_pivot': options.get_double(option_section, 'k_pivot')
    }
    return config


def execute(block, config):
    k_pivot = config['k_pivot']
    # Load in spec_index table from datablock onto a Pandas df
    file = block[spec_file, "spec_index_table"]
    df = pd.read_csv(file, sep="\t", header=None, names=tab_names)
    print df
    # transform k->ln k & A->ln(A)
    # df.transform(np.log)
    # ln_k_piv = np.log(k_pivot)
    # Fit a 2nd order spline to the ln(A_s) vs. ln(k) data
    spline_func_s = splrep(df["k_phys"], df["A_s"], k=2)
    # Evaluate function @ k_piv
    val_s = splev(k_pivot, spline_func_s)
    # Evaluate derivative dln(P_s)/dln(k) @ k=k_piv
    deriv_s = splev(k_pivot, spline_func_s, der=1)
    # n_s = deriv_s + 1
    n_s = deriv_s*(k_pivot/val_s) + 1

    # Fit a 2nd order spline to the ln(A_t) vs. ln(k) data
    spline_func_t = splrep(df["k_phys"], df["A_t"], k=2)
    # Evaluate function @ k_piv
    val_t = splev(k_pivot, spline_func_t)
    # Evaluate derivative n_t = dln(P_t)/dln(k) @ k=k_piv
    n_t = splev(k_pivot, spline_func_t, der=1)*(k_pivot/val_t)

    # create random filename for plots
    filename = uuid.uuid4().hex
    # Folder name
    folder = "Spec_Plots/"

    # Now make a plot of the scalar data and fit
    plt.figure()
    plt.plot(df["k_phys"], df["A_s"], 'o', label="data")
    plt.plot(df["k_phys"], splev(df["k_phys"], spline_func_s), label="spline")
    tangent_s = val_s + deriv_s*df["k_phys"]
    plt.plot(df["k_phys"], tangent_s, label="n_s tangent: " + str(n_s))
    plt.xlabel("ln(k)")
    plt.ylabel("ln(A_s)")
    plt.legend(loc='best')
    plt.savefig(folder+"scalar"+filename+".pdf")
    plt.close()

    # Now make a plot of the tensor data and fit
    plt.figure()
    plt.plot(df["k_phys"], df["A_t"], 'o', label="data")
    plt.plot(df["k_phys"], splev(df["k_phys"], spline_func_t), label="spline")
    tangent_t = val_t + n_t*df["k_phys"]
    plt.plot(df["k_phys"], tangent_t, label="n_t tangent: " + str(n_t))
    plt.xlabel("ln(k)")
    plt.ylabel("ln(A_t)")
    plt.legend(loc='best')
    plt.savefig(folder+"tensor_"+filename+".pdf")
    plt.close()

    print n_s
    print n_t

    # Return the n_s and n_t values to the datablock
    block[twopf_name, "n_s_scipy"] = n_s
    block[twopf_name, "n_t_scipy"] = n_t

    # We tell CosmoSIS that everything went fine by returning zero
    return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
