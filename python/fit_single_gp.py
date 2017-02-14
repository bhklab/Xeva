#!/usr/bin/env python
import pandas as pd
import numpy as np

from gp import normalize_data, fit_gaussian_process

if __name__ == '__main__':

    # import CSV files
    df_time = pd.read_csv('python/temptime.csv', delimiter=" ", header=None).dropna()
    df_volume = pd.read_csv('python/tempvolume.csv', delimiter=" ", header=None).dropna()

    # Replace any zeros in volume size with something that won't error out
    df_volume.replace(0, 0.000001, inplace=True)
    x, y, y_norm = normalize_data(df_time, df_volume)

    fit_gp = fit_gaussian_process(x, y_norm, len(df_time.columns))

    # Print out results of GP
    increment_by = 0.25
    predict_x = np.arange(0, max(x), increment_by)

    with open("python/tempresults.csv", "w") as f:
        print("x", ",", "prediction", ",", "variance", file=f)
        for i in predict_x:
            prediction_at_i = fit_gp.predict(np.asarray([[i]]))
            print(i, ",", prediction_at_i[0][0][0], ",", prediction_at_i[1][0][0], file=f)
