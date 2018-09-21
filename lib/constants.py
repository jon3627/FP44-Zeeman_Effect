import os

data_path = '../data'
spectroscopy_path = os.path.join(data_path, 'Spectroscopy')
zeeman_path = os.path.join(data_path, 'Zeeman_without_filter')
NeCd_path = os.path.join(spectroscopy_path, 'NeonCadmium.txt')
B_field_path = os.path.join(data_path, 'execution.xlsx')

PLOT_MARKERSIZE = 4
ERRORBAR_LINEWIDTH = 1
ERRORBAR_CAPSIZE = 1
COLORS = {'1': 'black', '2': 'orangered', '3': 'cornflowerblue'}
N_LUMMERGEHRCKE = 1.4567
THICKNESS = 4.04*1e-3 #m