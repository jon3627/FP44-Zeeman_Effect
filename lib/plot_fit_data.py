from math import ceil
import matplotlib.pyplot as plt
import numpy as np
import lib.constants as c
import pandas as pd
from scipy import sqrt
from lmfit import Parameters, minimize
from lmfit.models import Model, LinearModel, LorentzianModel, GaussianModel
from scipy.optimize import curve_fit


def plot_dfs(dfs, max_column_number, zeeman=None):
    plot_columns = ceil(len(dfs) / max_column_number)
    fig, axis = plt.subplots(plot_columns, max_column_number,
                             figsize=(15, plot_columns * 5), facecolor='w', edgecolor='k')
    if not max_column_number == 1:
        axis = axis.ravel()
    else:
        axis = np.array([axis])
    for i, df in enumerate(dfs):
        if zeeman is not None:
            df, filename = df
            axis[i].set_title(filename)
        axis[i].plot(df.index, df.iloc[:, 0].values, '.', markersize=c.PLOT_MARKERSIZE)
        if df.get('Best fit') is not None:
            axis[i].plot(df.index, df['Best fit'].values, 'r-', markersize=c.PLOT_MARKERSIZE)
    for i in range(len(dfs), max_column_number * plot_columns):
        axis[i].set_visible(False)
    plt.tight_layout()


def plot_B_field_dfs(df, fit_data=None):
    fig, axis = plt.subplots(1, 1, figsize=(15, 5), facecolor='w', edgecolor='k')
    if df.get('Best fit') is not None:
        print(*fit_data.loc[0].values)
        axis.plot(df.loc[:, 'Current I [A]'], df.loc[:, 'Best fit'], color=c.COLORS['2'],
                  linestyle='-', markersize=c.PLOT_MARKERSIZE,
                  label='linear fit: \n slope: ${} \pm {}$ \n intersect: ${} \pm {}$'
                  .format(*np.round(fit_data.loc[0].values, 5)))
        lower = 0
        middle = 9
        upper = 10
    else:
        lower = 1
        middle = 10
        upper = 12
    axis.errorbar(df.loc[lower:middle, 'Current I [A]'],
                  df.loc[lower:middle, 'magnetic field mean [G]'],
                  xerr=df.loc[lower:middle, 'Current error [A]'],
                  yerr=df.loc[lower:middle, 'error magnetic field mean [G] (gauss)'],
                  fmt='.', color=c.COLORS['1'], label='increasing',
                  elinewidth=c.ERRORBAR_LINEWIDTH,
                  capsize=c.ERRORBAR_CAPSIZE)
    axis.errorbar(df.loc[upper:, 'Current I [A]'],
                  df.loc[upper:, 'magnetic field mean [G]'],
                  xerr=df.loc[upper:, 'Current error [A]'],
                  yerr=df.loc[upper:, 'error magnetic field mean [G] (gauss)'],
                  fmt='.', color=c.COLORS['3'], label='decreasing',
                  elinewidth=c.ERRORBAR_LINEWIDTH,
                  capsize=c.ERRORBAR_CAPSIZE)

    axis.set_xlabel('current I [A]')
    axis.set_ylabel('magnetic field B [G]')
    axis.legend()


def fit_b_field_df(df):
    fit_data = {
        'slope': [],
        'slope_err': [],
        'intersect': [],
        'intersect_err': [],
        'redchi': []
    }

    df_clean = df.drop(df.index[[0, 10, 11]])
    x = df_clean.loc[:, 'Current I [A]'].values
    y = df_clean.loc[:, 'magnetic field mean [G]'].values
    xerr = df_clean.loc[:, 'Current error [A]'].values
    yerr = df_clean.loc[:, 'error magnetic field mean [G] (gauss)'].values

    p = Parameters()
    p.add('slope', value=1.)
    p.add('intersect', value=1.)
    mi = minimize(linear_residual, p, args=(x,), kws={'y': y, 'xerr': xerr, 'yerr': yerr}, method='powell')
    mi = minimize(linear_residual, mi.params, args=(x,), kws={'y': y, 'xerr': xerr, 'yerr': yerr}, method='leastsq')

    df_clean['Best fit'] = linear_residual(mi.params, x)
    fit_data = _save_fit_params(mi, fit_data)
    fit_data = pd.DataFrame(data=fit_data)

    return df_clean, fit_data


def fit_dfs(dfs, gaussiannumber, init_params):
    fit_stat = _initialize_fit_data_df(gaussiannumber)
    dfs_fitted = []

    for i, df in enumerate(dfs):
        model = _model_(gaussiannumber)
        params = _get_init_params(gaussiannumber, init_params, model, i)
        fit = model.fit(df.iloc[:, 0].values, x=df.index.values, params=params, method='least_squares')

        fit_stat = _save_fit_params(fit, fit_stat)
        df = _save_fit_in_df(df=df, fit=fit)
        dfs_fitted.append(df)

    fit_df = pd.DataFrame(data=fit_stat)
    return dfs_fitted, fit_df


def _initialize_fit_data_df(gaussiannumber):
    fit_params = {}
    lin_count = 1
    for i in range(1, gaussiannumber + 1):
        fit_params['gauss{}_amp'.format(i)] = []
        fit_params['gauss{}_cen'.format(i)] = []
        fit_params['gauss{}_sig'.format(i)] = []
    for i in range(1, lin_count + 1):
        fit_params['linear{}_intercept'.format(i)] = []
        fit_params['linear{}_slope'.format(i)] = []
    keys = list(fit_params.keys())
    for key in keys:
        fit_params[key + '_err'] = []
    fit_params['redchi'] = []
    return fit_params


def gaussian(x, amp, cen, sig):
    return amp / (np.sqrt(2 * np.pi) * sig) * np.exp(-((x - cen) / sig) ** 2 / 2)


def _model_(gaussiannumber):
    model = LinearModel(independent_vars=['x'], prefix='linear1_')
    for i in range(1, gaussiannumber + 1):
        model += Model(gaussian, independent_vars=['x'], prefix='gauss{}_'.format(i))
    return model


def _save_fit_in_df(df, fit):
    x_init = df.index.values
    df['Best fit'] = fit.eval(x=x_init)
    return df


def _save_fit_params(fit, fit_data):
    fit_params = fit.params
    fit_data['redchi'].append(fit.redchi)
    for key in fit_data.keys():
        if not key[-4:] == '_err' and not key == 'redchi' and not key == 'file':
            fit_data[key].append(fit_params[key].value)
            fit_data[key + '_err'].append(fit_params[key].stderr)
    return fit_data


def _get_init_params(fct, init_params, model, i):
    params = model.make_params()
    if params.keys() == init_params.keys():
        for param in params.keys():
            params[param].set(init_params[param])
    else:
        raise UserWarning('provided and expected parameters do not match')
    return params


def linear_residual(pars, x, y=None, xerr=None, yerr=None):
    vals = pars.valuesdict()
    slope = vals['slope']
    intersect = vals['intersect']
    model = slope * x + intersect
    if y is None:
        return model
    resids = model - y
    if xerr is None:
        return sqrt(resids ** 2 / yerr ** 2)
    if yerr is None:
        return sqrt(resids ** 2 / xerr ** 2)
    if xerr is None and yerr is None:
        return resids

    weighted = resids ** 2 / (xerr ** 2 + yerr ** 2)
    weighted = np.array(weighted, dtype=np.float)
    weighted = sqrt(weighted)
    return weighted
