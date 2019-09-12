import pandas as pd
import numpy as np


class Values:
    def __init__(self, R, N_blades, P_rated, V_wind_in, V_wind_out, rho, blade_table, FFA_W3_241, FFA_W3_301,
                 FFA_W3_360, FFA_W3_480, FFA_W3_600, cylinder):
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
        self.cylinder = cylinder

    def interpolate_r(self, R_needed):
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

    def interpolate_tc(self, tc, alpha):
        if tc >= 24.1 and tc < 31.0:
            range = 1
            lower = 24.1
            upper = 31.0
            first_table = self.FFA_W3_241
            second_table = self.FFA_W3_301
        elif tc >= 31.0 and tc < 36.0:
            range = 2
            lower = 31.0
            upper = 36.0
            first_table = self.FFA_W3_301
            second_table = self.FFA_W3_360
        elif tc >= 36.0 and tc < 48.0:
            range = 3
            lower = 36.0
            upper = 48.0
            first_table = self.FFA_W3_360
            second_table = self.FFA_W3_480
        elif tc >= 48.0 and tc < 60.0:
            range = 4
            lower = 48.0
            upper = 60.0
            first_table = self.FFA_W3_480
            second_table = self.FFA_W3_600
        elif tc >= 60.0 and tc <= 100:
            range = 5
            lower = 60.0
            upper = 100.
            first_table = self.FFA_W3_600
            second_table = self.cylinder
        else:
            raise ValueError('Value out of range')
        position_tc = (tc - lower) / (upper - lower)
        airfoil_table = pd.DataFrame(data=first_table, columns=['alpha', 'cl', 'cd', 'cm'])
        airfoil_table['cl'] = np.transpose((first_table[:, 1] + position_tc * (
                second_table[:, 1] - first_table[:, 1])))
        airfoil_table['cd'] = np.transpose((first_table[:, 2] + position_tc * (
                second_table[:, 2] - first_table[:, 2])))
        distances_alpha = airfoil_table['alpha'] - alpha
        closest_indexes = distances_alpha.abs().argsort()[:2]
        upper_index = closest_indexes[1]
        lower_index = closest_indexes[0]
        upper_row = airfoil_table.iloc[upper_index]
        lower_row = airfoil_table.iloc[lower_index]
        position_alpha = (alpha - lower_row['alpha']) / (upper_row['alpha'] - lower_row['alpha'])
        cl = lower_row['cl'] + (lower_row['cl'] - upper_row['cl']) * position_alpha
        cd = lower_row['cd'] + (lower_row['cd'] - upper_row['cd']) * position_alpha
        return cl, cd


Values1 = Values(89.17, 3, 10000, 4, 25, 1.225, pd.read_excel("Assignment1.xlsx"), np.loadtxt('FFA-W3-241.txt'),
                 np.loadtxt('FFA-W3-301.txt'), np.loadtxt('FFA-W3-360.txt'), np.loadtxt('FFA-W3-480.txt'),
                 np.loadtxt('FFA-W3-600.txt'), np.loadtxt('cylinder.txt'))
alpha = 4.1
R = 11
c, beta, tc = Values1.interpolate_r(R)
cl, cd = Values1.interpolate_tc(tc, alpha)
pass
