#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:26:23 2019

@author: nickleenders
"""

import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from Assignment_2_Q2 import Vg,Ig,omegaEl,RPM,PlossPhase,Pmech
from Assignment2_1 import Z1,Z2,Zc,Ztotal,Zcap,Efficiency,Pg,S_poc

Apparent_power_loss = np.sqrt(S_poc.real**2+S_poc.imag**2)/(3*Pmech)*100
V_poc = 33*10**3
S_poc_g = 4.2*10**6
P_poc = S_poc_g*Apparent_power_loss
Igrid = P_poc/V_poc