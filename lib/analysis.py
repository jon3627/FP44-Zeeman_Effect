import matplotlib.pyplot as plt
import lib.constants as c
import scipy.constants as const
import numpy as np
from scipy.optimize import curve_fit
from lib.util import linear


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

    plt.errorbar(ref, cen, fmt='.', yerr=sig,label='data points', markersize=c.PLOT_MARKERSIZE)
    plt.plot(ref, fitfunction(ref, *popt), label='linear fit', color='r')
    plt.hlines(y=fit_df_spec.loc[0, 'gauss3_cen'], xmin=ref[0], xmax=ref[-1], color='black', label='measurement of the Cd-line')
    plt.ylabel('position [px]')
    plt.xlabel('wavelength [nm]')
    plt.legend()

    # find intersection, only consider error of measured line position
    cd_line = (fit_df_spec.loc[0, 'gauss3_cen'] - b) / a
    cd_line_err = fit_df_spec.loc[0, 'gauss3_sig'] / a
    print('Cd-line in nm: {}, error in nm: {}'.format(cd_line, cd_line_err))
    return cd_line, cd_line_err


def calculate_Bohr_magneton(fit_df, cd_line, cd_line_err, B_field_fit_data):
    df, filename = fit_df
    current = int(filename[4:-1])
    B_field = linear(B_field_fit_data, current)
    orders = df.index.values[::-1]
    pi_lines = df['gauss2_cen'].values
    pi_lines_err = df['gauss2_sig'].values

    # fit and plot
    def fitfunction(x, a, b, c):
        return a*x**2 + b*x + c

    popt, pcov = curve_fit(fitfunction, pi_lines, orders)

    plt.errorbar(pi_lines, orders, fmt='.', xerr=pi_lines_err, label='position of the $\pi$-line')
    plt.plot(pi_lines, fitfunction(pi_lines, *popt), label='quadratic fit', color='r')
    plt.legend()
    plt.xlabel('position [px]')
    plt.ylabel('order')

    daDa = (fitfunction(df['gauss1_cen'].values,*popt)-fitfunction(df['gauss3_cen'].values,*popt))/2
    daDa_err = np.std(daDa)/np.sqrt(len(daDa))
    daDa = np.mean(daDa)
    dlambda = daDa*(cd_line*1e-9)**2/(2*c.THICKNESS*np.sqrt(c.N_LUMMERGEHRCKE**2-1))
    dlambda_err = dlambda*np.sqrt((daDa_err/daDa)**2+(2*cd_line_err/cd_line)**2)

    dEnergy = const.h*const.c/(cd_line*1e-9)-const.h*const.c/(cd_line*1e-9+dlambda)
    dEnergy_err = np.sqrt(((const.h*const.c/(cd_line**2*1e-9)-const.h*const.c*1e-9/(cd_line*1e-9+dlambda)**2)*cd_line_err)**2+
                          (const.h*const.c/(cd_line*1e-9+dlambda)**2*dlambda_err)**2)
    Bohr_magneton = dEnergy/B_field[0]
    Bohr_magneton_err = Bohr_magneton*np.sqrt((dEnergy_err/dEnergy)**2+(B_field[1]/B_field[0])**2)
    print('Bohr magneton in J/T for I='+str(current)+'A: '+str(Bohr_magneton)+', error in J/T: '+str(Bohr_magneton_err))
    return Bohr_magneton, Bohr_magneton_err, dEnergy, dEnergy_err, B_field[0], B_field[1]


def calculate_Bohr_magneton_meth2(arr, arr_err, B_field, B_field_err):
    # fit and plot
    def fitfunction(x, a):
        return a*x

    popt, pcov = curve_fit(fitfunction, B_field, arr, sigma=arr_err)

    plt.errorbar(B_field, arr, fmt='.', yerr=arr_err, xerr=B_field_err, label='data points')
    plt.plot(B_field, fitfunction(B_field, *popt), label='linear fit, intercept = 0', color='r')
    plt.legend()
    plt.ylabel('$\Delta E$ [J]')
    plt.xlabel('magnetic field [T]')
    return popt[0], np.sqrt(np.diag(pcov))[0]


