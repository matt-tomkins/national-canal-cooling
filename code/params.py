''' Usage: This code contains parameters used in national-canal-cooling.py and will not produce results if run in isolation '''

# Required libraries
from math import e as EULER
from decimal import Decimal
from statistics import mean
from scipy.constants import physical_constants as pyc

# Set reference ("asphalt", "concrete")
REFERENCE_MATERIAL = "asphalt"

''' Main parameters utilised in canal-cooling.py  ''' 

VEC = 3                                     # Vertical extinction coefficient ------------------------ (described https://www.waterontheweb.org/under/lakeecology/04_light.html)
SHADE = 0.2                                 # Relative light intensity from the shade ---------------- (https://doi.org/10.1002/(SICI)1099-1085(199701)11:1%3C79::AID-HYP404%3E3.0.CO;2-N) 
L = 2454.9                                  # Latent heat of vapourisation (J/g) --------------------- (https://doi.org/10.1002/(SICI)1099-1085(199701)11:1%3C79::AID-HYP404%3E3.0.CO;2-N)
RV = 461.5                                  # Gas constant for water vapour -------------------------- (Lide, David R. 1992. CRC Handbook of Chemistry and Physics, 1992-1993 : A Ready-Reference Book of Chemical and Physical Data. 73rd ed. Boca Raton (Florida): CRC Press.)
LRV = (L*1000)/RV                           # Ratio of latent heat to water vapour gas constant ------ (calculated)
RIW = 1.33                                  # Refractive index of water ------------------------------ (Hecht, Eugene. Optics, Global Edition, Pearson Education, Limited, 2016)
RIA = 1                                     # Refractive index of air -------------------------------- (Hecht, Eugene. Optics, Global Edition, Pearson Education, Limited, 2016)
E = 0.95                                    # Emissivity of water ------------------------------------ (https://www.engineeringtoolbox.com/radiation-heat-emissivity-d_432.html)
SB = pyc['Stefan-Boltzmann constant'][0]    # Stefan-Boltzmann Constant ~5.67*10**-8  ---------------- (https://docs.scipy.org/doc/scipy/reference/constants.html)
TPW = 273.16                                # Triple point of water (K) ------------------------------ (https://doi.org/10.1126/science.191.4233.1261, https://doi.org/10.1256/qj.04.94)
SHC = 4.186                                 # Specific heat capacity of water (kJ/kg K) -------------- (common knowledge e.g., https://www.britannica.com/science/specific-heat)
SHC_WATER = SHC * 1000                      # Specific heat capacity of water (J/kg K) --------------- (as above)
HC_WATER = 0.58                             # Heat conductivity of water (W/m K) --------------------- (https://energyeducation.ca/encyclopedia/Thermal_conductivity)
CONDUCTIVITY = 0.58                         # Thermal conductivity (kc) ------------------------------ (https://energyeducation.ca/encyclopedia/Thermal_conductivity)
CONVECTION = 200                            # Convection coefficent (hc) ----------------------------- (source?)
BASE_PRESSURE = 0.611657                    # Baseline Pressure at Triple Point (kPa) ---------------- (https://doi.org/10.1256/qj.04.94)
DEW_POINT_B = 17.67                         # Constant (b) to calculate dew point (in evaporation) --- (https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2)
DEW_POINT_C = 243.5                         # Constant (c) to calculate dew point (in evaporation) --- (https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2)

''' Chosen canal parameters '''

FLOAT_INTERVAL = 0.2                        # Depth of each canal layer (m)
INTERVAL = Decimal(str(FLOAT_INTERVAL))     # As above, but stored as decimal

#============== Concrete parameters ==============# 
if REFERENCE_MATERIAL == "concrete":
    
    R_DEPTH = 0.05                                          # Concrete section depth (m)
    R_EMISSIVITY = 0.85                                     # Emissivity ------------------------ (https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html)
    R_ALBEDO = 0.30                                         # Albedo ---------------------------- (https://www.mdpi.com/2076-3417/12/4/1955 [0.25 - 0.40]), (https://doi.org/10.3390/su132011448 [dark 0.2 - 0.3], [most 0.4 - 0.6])
    R_ABSORBIVITY = 1 - R_ALBEDO                            # Absorbivity ----------------------- (1 - albedo)
    R_CAPACITY = 880                                        # Specific heat capacity (J/kg K) --- (https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html)
    R_DENSITY = 2400                                        # Density (kg/m3) ------------------- (https://doi.org/10.1016/j.jobe.2018.07.002, a 'typical value' might be 2400 kg/m3, but can be made at many densities)
    R_CONDUCTIVITY = 0.0625*EULER**(0.0015 * R_DENSITY)     # Thermal conductivity (W/m K) ------ (https://www.sciencedirect.com/science/article/pii/S2352710218304650, "There is a significant relationship between the unit weight of concrete and the thermal conductivity value")
                                                            # ----------------------------------- Calculate using (1) k = 0.0625e^0.0015p, where k = conductivity, and p = density (https://doi.org/10.1016/j.jobe.2018.07.002) or
                                                            # ----------------------------------- Calculate using (2) k = 0.0865e^0.00125p (https://www.vandidad-co.com/Uploadfiles/CkEditor/Files/%D9%85%D8%AD%D8%B5%D9%88%D9%84%D8%A7%D8%AA/213r_87.pdf)
                                                            # ----------------------------------- Alternatively, could calculate using a subset of Asadi et al. (2018) e.g., k = 0.0415e^0.0019p

#============== Asphalt parameters ==============#  
elif REFERENCE_MATERIAL == "asphalt":

    R_DEPTH = 0.05                      # Asphalt section depth (m)
    R_EMISSIVITY = 0.85                 # Emissivity ------------------------ (https://doi.org/10.1016/j.conbuildmat.2019.05.091, "In general, all the values were in between 0.70 and 0.95")
                                        # ----------------------------------- (Other sources: https://doi.org/10.1061/(ASCE)0899-1561(2007)19:8(683), https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html [0.93])
    R_ALBEDO = 0.1                      # Albedo ---------------------------- (https://doi.org/10.3390/app12041955 [0.1 - 0.2], https://www.eupave.eu/wp-content/uploads/FACT-SHEET-High-albedo-v26062020.pdf [0.05 - 0.15])
    R_ABSORBIVITY = 1 - R_ALBEDO        # Absorbivity ----------------------- (1 - albedo)
    R_DENSITY = 2100                    # Density (kg/m3) ------------------- (https://doi.org/10.3141/1813-09 [Figure 2])
    R_CAPACITY = 1020                   # Specific heat capacity (J/kg K) --- (https://doi.org/10.1016/j.conbuildmat.2021.124757 [~1020 ± 180 = mean of 45 measurements of different asphalt mixtures and ages])
                                        # ----------------------------------- (https://www.engineeringtoolbox.com/specific-heat-solids-d_154.html [920])
    R_CONDUCTIVITY = 1.65               # Thermal conductivity (W/m K) ------ (https://doi.org/10.1016/j.conbuildmat.2021.124757 [~1.65 ± 0.20 = mean of 45 measurements of different asphalt mixtures and ages])
                                        # ----------------------------------- (https://help.iesve.com/ve2021/table_6_thermal_conductivity__specific_heat_capacity_and_density.htm [0.5])
                                        # ----------------------------------- (https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html [0.75])
                                        # ----------------------------------- (https://doi.org/10.1016/j.rineng.2022.100827 [~1.2])
                                        # ----------------------------------- (https://doi.org/10.2514/6.1983-1598 [highly variable])

''' Ground parameters  ''' 

G_CONDUCTIVITY = 1          # Thermal conductivity (W/mK)
G_CAPACITY = 1900           # Specific heat capacity (J/kg K) 
G_DENSITY = 1500            # Density (kg/m3)

''' Air parameters '''

A_CAPACITY = 0.716          # Specific heat capacity (kJ/kg*k) of 0.716 (Specific heat at constant volume) or 1.005 (Specific heat at constant pressure)
A_DENSITY = 1.19            # Density (kg/m^3)
A_HEIGHT = 10               # Height of air temperature model domain (m)
A_R_SPECIFIC = 287.05       # Specific gas constant for dry air (J/kg·K)
A_CONDUCTIVITY = 0.025681   # Conductivity of air (15°C, 1 bar, Carroll et al., 1968)

# Seconds in 15 minutes
MODEL_INTERVAL = 3600 / 4

# Water temperatures (°C) at 50 - 100 cm depth, for 9 stations for [01/01/2022 00:00:00], sourced from the Canal and River Trust
START_CELCIUS = [9.42, 8.82, 9.88, 9.16, 8.2, 9.45, 8.64, 8.72, 8.69]
START_KELVIN = [t + 273.15 for t in START_CELCIUS]

# Average temperature (K)
MEAN_KELVIN = round(mean(START_KELVIN))

# Fixed cloud height parameter, used in thermal_absorption(), based on Goforth et al. (2002) https://doi.org/10.1117/12.459570 
# Cloud height data are not available so we use a fixed value: K=0.34 (cloud height<2km), K=0.18 (2 km < cloud height < 5 km), or K=0.06 (5km < cloud height)
CLOUD_HEIGHT = 0.06