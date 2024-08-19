''' Usage: This code contains parameters used in national-canal-cooling.py and will not produce results if run in isolation '''

# Required libraries
from decimal import Decimal
from statistics import mean

''' Main parameters utilised in canal-cooling.py  ''' 

VEC = 3                 # Vertical extinction coefficient
SHADE = 0.2             # Relative light intensity from the shade
G = 9.807               # Gravity (N/kg)
L = 2454.9              # Latent heat of vapourisation (J/g)
RV = 461.5              # Gas constant for water vapour
LRV = (L*1000)/RV       # Ratio of latent heat to water vapour gas constant
RIW = 1.33              # Refractive index of water
RIA = 1                 # Refractive index of air
OSI = 476               # Overhead sun intensity (W/m^2)
E = 0.95                # Emissivity of water
SB = 5.67*10**-8        # Stefan-Boltzmann Constant
TPW = 273.0             # Triple point of water (K)
SHC = 4.186             # Specific heat capacity of water (kJ/kg K)
SHC_WATER = SHC * 1000  # Specific heat capacity of water (J/kg K)
HC_WATER = 0.58         # Heat conductivity of water (W/m K)
CONDUCTIVITY = 0.58     # Thermal conductivity (kc)
CONVECTION = 200        # Convection coefficent (hc)
BASE_PRESSURE = 0.611   # Baseline Pressure at Triple Point (kPa)
WATER_WEIGHT = 999286   # Specific weight of water (g/m^3)
DEW_POINT_B = 17.67     # Constant (b) to calculate dew point (in evaporation)
DEW_POINT_C = 243.5     # Constant (c) to calculate dew point (in evaporation)

''' Chosen canal parameters '''

FLOAT_INTERVAL = 0.2                     # Depth of each canal layer (m)
INTERVAL = Decimal(str(FLOAT_INTERVAL))  # As above, but stored as decimal

''' Concrete parameters  ''' 

C_DEPTH = 0.05          # Concrete section depth (m)
C_EMISSIVITY = 0.85     # Emissivity
C_ALBEDO = 0.15         # Albedo
C_ABSORBIVITY = 0.85    # Absorbivity
C_CONDUCTIVITY = 1.21   # Thermal conductivity (W/m K) 
C_CAPACITY = 921        # Specific heat capacity (J/kg K) 
C_DENSITY = 2238        # Density (kg/m3)

''' Ground parameters  ''' 

G_CONDUCTIVITY = 1      # Thermal conductivity (W/mK)
G_CAPACITY = 1900       # Specific heat capacity (J/kg K) 
G_DENSITY = 1500        # Density (kg/m3)

''' Air parameters '''

A_CAPACITY = 0.716      # Specific heat capacity (kJ/kg*k)	
A_DENSITY = 1.19        # Density (kg/m^3)
A_HEIGHT = 10           # Height of air temperature model domain (m)

# Water temperatures (°C) at 50 - 100 cm depth, for 9 stations for [01/01/2022 00:00:00]
START_CELCIUS = [9.42, 8.82, 9.88, 9.16, 8.2, 9.45, 8.64, 8.72, 8.69]
START_KELVIN = [t + 273.15 for t in START_CELCIUS]

# Average temperature (k)
MEAN_KELVIN = round(mean(START_KELVIN))