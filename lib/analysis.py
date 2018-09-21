import matplotlib.pyplot as plt
import lib.constants as c
import numpy as np
from scipy.optimize import curve_fit


def find_Cd_line(fit_df_spec, ref):
    cen = [fit_df_spec.loc[0, 'gauss1_cen'], fit_df_spec.loc[0, 'gauss2_cen'],
           fit_df_spec.loc[0, 'gauss4_cen'], fit_df_spec.loc[0, 'gauss5_cen']]
    sig = [fit_df_spec.loc[0, 'gauss1_sig'], fit_df_spec.loc[0, 'gauss2_sig'],
           fit_df_spec.loc[0, 'gauss4_sig'], fit_df_spec.loc[0, 'gauss5_sig']]

    # fit and plot
    def fitfunction(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(fitfunction, ref, cen, sigma=sig)
    a = popt[0]
    b = popt[1]
    a_err = np.sqrt(np.diag(pcov))[0]
    print('rel. error of the fit parameter a: ' + str(abs(a_err / a)))
    b_err = np.sqrt(np.diag(pcov))[1]
    print('rel. error of the fit parameter b: ' + str(abs(b_err / b)))
    print('rel. error of the measured line position: ' + str(
        fit_df_spec.loc[0, 'gauss3_sig'] / fit_df_spec.loc[0, 'gauss3_cen']))

    plt.errorbar(ref, cen, fmt='.', yerr=sig, markersize=c.PLOT_MARKERSIZE)
    plt.plot(ref, fitfunction(ref, *popt))
    plt.hlines(y=fit_df_spec.loc[0, 'gauss3_cen'], xmin=ref[0], xmax=ref[-1], color='black')

    # find intersection, only consider error of measured line position
    cd_line = (fit_df_spec.loc[0, 'gauss3_cen'] - b) / a
    cd_line_err = fit_df_spec.loc[0, 'gauss3_sig'] / a
    print('Cd-line in nm: {}, error in nm: {}'.format(cd_line, cd_line_err))
