import pandas as pd
import numpy as np


class Values:
    def __init__(self, R, N_blades, P_rated, V_wind_in, V_wind_out, rho, blade_table, FFA_W3_241, FFA_W3_301,
                 FFA_W3_360, FFA_W3_480, FFA_W3_600):
        self.R = R  # Rotor Radius in m
        self.N_blades = N_blades  # number of blades
        self.P_rated = P_rated  # Rated power in kW
        self.V_wind_in = V_wind_in  # Cut in wind speed in m/s
        self.V_wind_out = V_wind_out  # Cut out wind speed in m/s
        self.rho = rho  # air density in kg/m^3
        self.blade_table = blade_table
        self.FFA_W3_241 = FFA_W3_241
        self.FFA_W3_301 = FFA_W3_301
        self.FFA_W3_360 = FFA_W3_360
        self.FFA_W3_480 = FFA_W3_480
        self.FFA_W3_600 = FFA_W3_600

    def interpolate(self, R_needed):
        distances = self.blade_table['r'] - R_needed
        closest_indexes = distances.abs().argsort()[:2]
        upper_index = closest_indexes[1]
        lower_index = closest_indexes[0]
        upper_row = self.blade_table.iloc[upper_index]
        lower_row = self.blade_table.iloc[lower_index]
        position = (R_needed - upper_row['r']) / (lower_row['r'] - upper_row['r'])
        c = upper_row['c'] + (lower_row['c'] - upper_row['c']) * position
        beta = upper_row['beta'] + (lower_row['beta'] - upper_row['beta']) * position
        tc = upper_row['t/c'] + (lower_row['t/c'] - upper_row['t/c']) * position
        return c, beta, tc


Values1 = Values(89.17, 3, 10000, 4, 25, 1.225, pd.read_excel("Assignment1.xlsx"), np.loadtxt('FFA-W3-241.txt'),
                 np.loadtxt('FFA-W3-301.txt'), np.loadtxt('FFA-W3-360.txt'), np.loadtxt('FFA-W3-480.txt'),
                 np.loadtxt('FFA-W3-600.txt'))
c, beta, tc = Values1.interpolate(22)
pass
