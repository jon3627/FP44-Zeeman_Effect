from lib.util import make_dfs, split_trans_dfs
from lib.plot_fit_data import plot_dfs, fit_dfs
from lib.analysis import find_Cd_line
import lib.constants as c
import matplotlib.pyplot as plt
import numpy as np

def main():
    # part2
    dfs_spec = make_dfs(c.NeCd_path)
    plot_dfs(dfs_spec, 1)
    plt.show()

    params = {'gauss1_amp': 2.2e8, 'gauss2_amp': 2.9e8, 'gauss3_amp': 2.35e8, 'gauss4_amp': 1.0e8, 'gauss5_amp': 0.35e8,
            'gauss1_cen': 20, 'gauss2_cen': 180, 'gauss3_cen': 480, 'gauss4_cen': 1050, 'gauss5_cen': 1270,
            'gauss1_sig': 6,'gauss2_sig': 7,'gauss3_sig': 6, 'gauss4_sig': 5, 'gauss5_sig': 5,
            'linear1_intercept': -1.67e7,
            'linear1_slope': 0}

    dfs_spec, fit_df_spec = fit_dfs(dfs_spec,gaussiannumber=5, init_params=params)
    plot_dfs(dfs_spec, 1)
    plt.show()

    neon_ref = np.array([638.29914, 640.2248, 650.65277, 653.28824])
    find_Cd_line(fit_df_spec, neon_ref)
    plt.show()

    #part1
    dfs_zeeman = make_dfs(c.zeeman_path)
    plot_dfs(dfs_zeeman, 3, zeeman=True)
    plt.show()

    split_params = [[(60,145),(180,265),(300,400),(450,540),(600,690),(765,860),(940,1050),
                     (1145, dfs_zeeman[0][0].index[-1])], #trans08A
                    [(40, 130), (165, 255), (285, 390), (425, 530), (570, 680), (730, 860), (900, 1045),
                     (1100, dfs_zeeman[0][0].index[-1])], #trans11A
                    [(30, 130), (150, 250), (270, 390), (405, 535), (550, 680), (700, 860), (880, 1060),
                     (1100, dfs_zeeman[0][0].index[-1])]] #trans13A
    dfs_split = split_trans_dfs(dfs_zeeman, split_params)
    for j in range(0, len(dfs_split)):
        plot_dfs(dfs_split[j][0], 3)
        plt.show()


if __name__ == '__main__':
    main()