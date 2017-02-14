import pandas as pd
import numpy as np

from gp import normalize_data, fit_gaussian_process, calculate_kl_divergence

if __name__ == '__main__':

    # import CSV files
    df_control_time = pd.read_csv('python/temp_control_time.csv', delimiter=" ", header=None).dropna()
    df_case_time = pd.read_csv('python/temp_case_time.csv', delimiter=" ", header=None).dropna()
    df_control_volume = pd.read_csv('python/temp_control_volume.csv', delimiter=" ", header=None).dropna()
    df_case_volume = pd.read_csv('python/temp_case_volume.csv', delimiter=" ", header=None).dropna()

    # Replace any zeros in volume size with something that won't error out
    df_control_volume.replace(0, 0.000001, inplace=True)
    df_case_volume.replace(0, 0.000001, inplace=True)

    control_x, control_y, control_y_norm = normalize_data(df_control_time, df_control_volume)
    case_x, case_y, case_y_norm = normalize_data(df_case_time, df_case_volume)

    control_gp = fit_gaussian_process(control_x, control_y_norm, len(df_control_time.columns))
    case_gp = fit_gaussian_process(case_x, case_y_norm, len(df_case_time.columns))

    kl_divergence = calculate_kl_divergence(control_x, case_x, control_gp, case_gp)

    with open("python/temp_statistics_results.csv", "w") as f:
        print("kl_divergence", file=f)
        print(kl_divergence, file=f)
