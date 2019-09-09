import pandas as pd
import numpy as np
class Values:
    def __init__(self):
        R = 89.17 # Rotor Radius in m
        N_blades = 3 # number of blades
        P_rated = 10000 # Rated power in kW
        V_wind_in = 4 # Cut in wind speed in m/s
        V_wind_in = 25  # Cut out wind speed in m/s
        rho = 1.225 # air density in kg/m^3
        blade_table = pd.read_excel("Assignment1.xlsx")
        FFA_W3_241 = np.loadtxt('FFA-W3-241.txt')
        FFA_W3_301 = np.loadtxt('FFA-W3-301.txt')
        FFA_W3_360 = np.loadtxt('FFA-W3-360.txt')
        FFA_W3_480 = np.loadtxt('FFA-W3-480.txt')
        FFA_W3_600 = np.loadtxt('FFA-W3-600.txt')
        pass

Values()


