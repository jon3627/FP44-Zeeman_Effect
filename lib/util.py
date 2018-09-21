import os
import pandas as pd
import numpy as np
import lib.constants as c


def get_txt_csv_xlsx(path):
    files = os.listdir(path)
    file_names = [file for file in files if file.endswith(('.txt', '.TXT', '.CSV', '.csv', '.xlsx', '.XLSX'))]
    file_paths = [os.path.join(path, file) for file in file_names]
    return file_names, file_paths


def make_dfs(path):
    dfs = []
    if path == c.NeCd_path:
        df = pd.read_table(path)
        df.set_index('X', inplace=True)
        df = df.dropna()
        dfs.append(df)
    if path == c.zeeman_path:
        file_names, file_paths = get_txt_csv_xlsx(path)
        for file_path, file_name in zip(file_paths, file_names):
            df = pd.read_table(file_path)
            df.set_index('X', inplace=True)
            df = df.dropna()
            dfs.append((df, file_name[:7]))
    if path == c.B_field_path:
        dfs = pd.read_excel(path, header=1, usecols='A,B,I:K', index_col=None)
    return dfs


def split_trans_dfs(dfs, params):
    arr_of_arr_of_dfs = []
    for j in range(3, len(dfs)):
        df, filename = dfs[j]
        new_dfs = []
        for i in range(0, len(params[j - 3])):
            start, stop = params[j - 3][i]
            new_df = df.loc[start:stop, :]
            new_dfs.append(new_df)
        arr_of_arr_of_dfs.append((new_dfs, filename))
    return arr_of_arr_of_dfs


def linear(fit_data, x, xerr=None):
    value = fit_data.loc[0, 'slope'] * x + fit_data.loc[0, 'intersect']
    if xerr is None:
        error = np.sqrt((fit_data.loc[0, 'slope_err'] * x) ** 2 + fit_data.loc[0, 'intersect_err'] ** 2)
    else:
        error = np.sqrt((fit_data.loc[0, 'slope'] * x * (fit_data.loc[0, 'slope_err'] / xerr)) ** 2
                        + + fit_data.loc[0, 'intersect_err'] ** 2)
    return (value, error)
