import numpy as np

def sat_pres_arm(T):
   #August-Roche-Magnus formula; from Wikipedia, C-C relation article; DOI:  10.2172/548871
   #T in K, e_w in hPa
   t = T - 273.15
   return 6.1094 * np.exp(17.625 * t / (t + 243.04)) 

def sat_pres_ice(T):
   #Over ice; (2) from Murphy and Koop QJRMS (2005)
   #T in K, e_w in hPa
   return np.exp(28.9074 - 6143.7/T)/100.

def vmr(rh, T, p):
   #return rh*sat_pres(T) / p
   return rh*sat_pres_ice(T) / p

scon = 1368.22
cosz = 0.9
alb = 0.3
emi = 0.99

cpair = 1003. # cp air units of J / K / kg

mdc = 0.658114 # Ratio of mass of dry air to mass of CO2 used by RRTM (see inatm routine)
mdw = 1.607793 # Ratio of mass of dry air to mass of H2O used by RRTM (see inatm routine)
mdo = 0.603428 # Ratio of mass of dry air to mass of O3 used by RRTM (see inatm routine)

zeroC = 273.15 
