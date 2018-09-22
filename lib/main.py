from lib.util import make_dfs, split_trans_dfs, linear
from lib.plot_fit_data import plot_dfs, fit_dfs, plot_B_field_dfs, fit_b_field_df
from lib.analysis import find_Cd_line, calculate_Bohr_magneton, calculate_Bohr_magneton_meth2
import lib.constants as c
import matplotlib.pyplot as plt
import numpy as np


def main():
    # part1
    # examination magnetic field for Zeeman effect

    df_B_field = make_dfs(c.B_field_path)
    # plot_B_field_dfs(df_B_field)

    df_B_field, B_field_fit_data = fit_b_field_df(df_B_field)
    plot_B_field_dfs(df_B_field, B_field_fit_data)
    print(B_field_fit_data)
    plt.show()
    # to get certain B-field-value [T] for given current use util.linear(B_field_fit_data, current_value)
    print(linear(B_field_fit_data, 8))

    # Zeeman effect
    dfs_zeeman = make_dfs(c.zeeman_path)
    plot_dfs(dfs_zeeman, 3, zeeman=True)
    plt.show()

    split_params = [[(60, 145), (180, 265), (300, 400), (450, 540), (600, 690), (765, 860), (940, 1050),
                     (1145, dfs_zeeman[0][0].index[-1])],  # trans08A
                    [(40, 130), (165, 255), (285, 390), (425, 530), (570, 680), (730, 860), (900, 1045),
                     (1100, dfs_zeeman[0][0].index[-1])],  # trans11A
                    [(30, 130), (150, 250), (270, 390), (405, 535), (550, 680), (700, 860), (880, 1060),
                     (1100, dfs_zeeman[0][0].index[-1])]]  # trans13A

    dfs_split = split_trans_dfs(dfs_zeeman, split_params)
    params_trans = [[{'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 90, 'gauss2_cen': 110, 'gauss3_cen': 130,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0}, #1
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 205, 'gauss2_cen': 230, 'gauss3_cen': 245,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0}, #2
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 335, 'gauss2_cen': 360, 'gauss3_cen': 285,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0}, #3
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 470, 'gauss2_cen': 505, 'gauss3_cen': 525,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}, #4
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 625, 'gauss2_cen': 650, 'gauss3_cen': 675,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}, #5
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 785, 'gauss2_cen': 810, 'gauss3_cen': 840,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}, #6
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 965, 'gauss2_cen': 995, 'gauss3_cen': 1030,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}, #7
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 1170, 'gauss2_cen': 1210, 'gauss3_cen': 1250,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}], #8

                    [{'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 65, 'gauss2_cen': 90, 'gauss3_cen': 110,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 1
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 185, 'gauss2_cen': 215, 'gauss3_cen': 235,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 2
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 310, 'gauss2_cen': 345, 'gauss3_cen': 370,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 3
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 450, 'gauss2_cen': 485, 'gauss3_cen': 510,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 4
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 595, 'gauss2_cen': 630, 'gauss3_cen': 670,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 5
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 750, 'gauss2_cen': 790, 'gauss3_cen': 835,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 6
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 930, 'gauss2_cen': 975, 'gauss3_cen': 1025,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 7
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.2e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 1135, 'gauss2_cen': 1185, 'gauss3_cen': 1250,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}],  # 8

                    [{'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 50, 'gauss2_cen': 90, 'gauss3_cen': 110,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 1
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 170, 'gauss2_cen': 210, 'gauss3_cen': 235,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 2
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1e8, 'gauss3_amp': 0.8e8,
                      'gauss1_cen': 300, 'gauss2_cen': 332, 'gauss3_cen': 370,
                      'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6,
                      'linear1_intercept': -1.25e7,
                      'linear1_slope': 0},  # 3
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.4e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 435, 'gauss2_cen': 475, 'gauss3_cen': 510,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 4
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.4e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 585, 'gauss2_cen': 620, 'gauss3_cen': 665,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 5
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.4e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 740, 'gauss2_cen': 785, 'gauss3_cen': 835,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 6
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.4e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 920, 'gauss2_cen': 965, 'gauss3_cen': 1025,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0},  # 7
                     {'gauss1_amp': 0.8e8, 'gauss2_amp': 1.4e8, 'gauss3_amp': 1e8,
                      'gauss1_cen': 1125, 'gauss2_cen': 1185, 'gauss3_cen': 1245,
                      'gauss1_sig': 15, 'gauss2_sig': 20, 'gauss3_sig': 15,
                      'linear1_intercept': -1.35e7,
                      'linear1_slope': 0}]]  # 8
    fit_dfs_arr = []
    for j in range(0, len(dfs_split)):
        new_dfs, fit_df = fit_dfs(dfs_split[j][0], gaussiannumber=3, init_params=params_trans[j])
        plot_dfs(new_dfs, 3, zeeman_split=dfs_split[j][1])
        plt.show()
        fit_dfs_arr.append((fit_df, dfs_split[j][1]))

    # part2
    dfs_spec = make_dfs(c.NeCd_path)
    plot_dfs(dfs_spec, 1)
    plt.show()

    params = [{'gauss1_amp': 2.2e8, 'gauss2_amp': 2.9e8, 'gauss3_amp': 2.35e8, 'gauss4_amp': 1.0e8, 'gauss5_amp': 0.35e8,
              'gauss1_cen': 20, 'gauss2_cen': 180, 'gauss3_cen': 480, 'gauss4_cen': 1050, 'gauss5_cen': 1270,
              'gauss1_sig': 6, 'gauss2_sig': 7, 'gauss3_sig': 6, 'gauss4_sig': 5, 'gauss5_sig': 5,
              'linear1_intercept': -1.67e7,
              'linear1_slope': 0}]

    dfs_spec, fit_df_spec = fit_dfs(dfs_spec, gaussiannumber=5, init_params=params)
    plot_dfs(dfs_spec, 1)
    plt.show()

    neon_ref = np.array([638.29914, 640.2248, 650.65277, 653.28824])
    cd_line, cd_line_err = find_Cd_line(fit_df_spec, neon_ref)
    plt.show()

    Bohr_magneton = []
    Bohr_magneton_err = []
    B_field = np.array([])
    B_field_err = []
    dEnergy = []
    dEnergy_err = []
    for j in range(0, len(fit_dfs_arr)):
        B_m, B_m_err, dE, dE_err, B_f, B_f_err = calculate_Bohr_magneton(fit_dfs_arr[j], cd_line, cd_line_err, B_field_fit_data)
        plt.show()
        Bohr_magneton.append(B_m)
        Bohr_magneton_err.append(B_m_err)
        dEnergy.append(dE)
        dEnergy_err.append(dE_err)
        B_field = np.append(B_field, B_f)
        B_field_err.append(B_f_err)

    Bohr_magneton_err = np.std(Bohr_magneton)/np.sqrt(len(Bohr_magneton)) + np.sqrt(np.dot(Bohr_magneton_err, Bohr_magneton_err))
    Bohr_magneton = np.mean(Bohr_magneton)
    print('Bohr magneton in J/T: '+str(Bohr_magneton)+', error in J/T: '+str(Bohr_magneton_err))
    print('difference to the theoretical value: ' + str(abs(Bohr_magneton-9.27400949*1e-24)/Bohr_magneton_err)+' sigma')

    Bohr_magneton, Bohr_magneton_err = calculate_Bohr_magneton_meth2(dEnergy, dEnergy_err, B_field, B_field_err)
    plt.show()
    print('Bohr magneton in J/T (fit method): '+str(Bohr_magneton)+', error in J/T: '+str(Bohr_magneton_err))
    print('difference to the theoretical value: ' + str(abs(Bohr_magneton-9.27400949*1e-24)/Bohr_magneton_err)+' sigma')


if __name__ == '__main__':
    main()
