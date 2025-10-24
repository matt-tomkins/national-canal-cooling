''' Usage: This code contains functions used in national-canal-cooling.py and will not produce results if run in isolation. '''

# Force use of Shapely 2.0
from os import environ
environ['USE_PYGEOS'] = '0'

# Required packages
import warnings
from json import load
from datetime import timedelta
from pandas import to_datetime
from numpy import savetxt, array
from datetime import datetime as dt
from pysolar.solar import get_altitude, radiation
from math import radians, cos, sqrt, sin, exp, log

# Import user-defined parameters
from params import * 

# Prevent geopandas CRS warning for pybdshadow (Context: https://github.com/geopandas/geopandas/issues/2606)
# shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)
warnings.filterwarnings("ignore", message="CRS not set for some of the concatenation inputs")

#------------------------- Reference modelling -------------------------#

def solar_intensity(date, longitude, latitude):
    '''
    > This function returns the intensity of the sun (W/m^2) based on its azimuth angle for a specific location and datetime
    > It also returns the solar altitude angle (°), which is used in solar_absorptivity()
    > It replaces the following functions from the original code:
        - sun_intensity(), sun_altitude_d(), sun_altitude_r(), HRA(), sun_azimuth_r(), sun_azimuth_d()
    > Verified using https://www.sunearthtools.com/dp/tools/pos_sun.php for Manchester and Cape Town
    > get_position() returns:
        - Sun altitude in radians 
        - Sun azimuth in radians, measured relative to south (east = negative, west = positive)
    '''
    # Return the solar altitude using Pysolar
    solar_altitude = get_altitude(latitude, longitude, date)

    # Clear sky irradiance on a surface perpendicular to sun (DNI; W/m^2), based on PySolar documentation (https://pysolar.readthedocs.io/en/latest/)
    irradiance = radiation.get_radiation_direct(date, solar_altitude)

    # Adjust for Global Horizontal Irradiance (GHI)
    irradiance_horizontal = irradiance * cos(radians(90 - solar_altitude))

    # Return outputs
    return irradiance_horizontal, solar_altitude

def model_reference_temperature(air_temperature, relative_humidity, cloud_height, cloud_cover, wind_speed, air_pressure, 
                               reference_temp, interval_datetime, longitude, latitude,
                               shading_proportion, canal_area):
    '''
    > A nested function to simplify the code structure
    > Takes an input reference temperature (K) and models a new reference temperature after 15 minutes using:
            - Air temperature (K)
            - Humidity (%)
            - Cloud height and cover (fractional)
            - Wind speed (m/s)
            - Shading proportion [0-1]
    > The following parameters are used:
        - Reference emissivity (R_EMISSIVITY) = 1 - albedo 
    > Similiar in scope to ConcreteNewTemp() from the original code, but individual functions are outlined more fully here
    '''

    # Model the incoming solar (shortwave) radiation (W, J/s), incorporating cloud cover
    csr = reference_solar_radiation(interval_datetime, longitude, latitude, shading_proportion, canal_area, cloud_cover)
    
    # Modify starting reference temperature (K) based on radiation input
    reference_temp = reference_after_sunlight(csr, reference_temp, canal_area)

    # Calculate reference emissions (W, J/s) based on surface concrete temperature
    emissions = reference_emissions(reference_temp, canal_area)
    
    # Thermal absorption (longwave) from sky (W, J/s)
    thermal_abs = thermal_absorption(air_temperature, relative_humidity, cloud_height, cloud_cover, canal_area) 

    # Reference material absorption (W, J/s)
    reference_absorbed = thermal_abs * R_EMISSIVITY
    
    # Net radiation (W, J/s) = total emissions based on surface temperature minus total thermal absorption from sky
    net_radiation =  emissions - reference_absorbed 

    # Calculate energy transfer to air (W, J/s)
    energy_to_air = reference_convection_air(reference_temp,
                                            air_temperature, 
                                            wind_speed, 
                                            canal_area, 
                                            air_pressure)
    
    # Determine heat transfer to layer above (W, J/s)
    heat_transfer = reference_heat_transfer(reference_temp, net_radiation, energy_to_air, canal_area)

    # Total energy change (W), converting from J/s to 15-minute model intervals
    total_energy_change = reference_energy_change(heat_transfer)

    # Corresponding temperature change (K)
    modelled_temperature_change = reference_temperature_change(total_energy_change, canal_area)

    # Final reference temperature for model run, overwriting earlier list
    reference_temp = [x + y for x, y in zip(reference_temp, modelled_temperature_change)]

    # Created dictionary to return variables used for modelling air temperature change
    # Variables measured in W (J/s), accounting for canal area (m/2), multipled by MODEL_INTERVAL (total energy transferred over 15m)
    reference_energy_values = {'net_radiation' : net_radiation * MODEL_INTERVAL,
                   'energy_to_air' : energy_to_air * MODEL_INTERVAL}

    # Return output
    return reference_temp, reference_energy_values

def cloud_cover_effect(cloud_cover):
    '''
    > Function to calculate a constant that represents the fraction of solar radiation that is blocked or diffused by clouds, following Nevins and Apell (2021)
    > Solar irradiance is minimally impacted up to ≈50% cloud cover but decreases by ≈67% at 100% cloud cover.
    > Approximated using a quadratic fit of y = 1 − 0.00243x − (4.24 × 10−5)x^2
    > y is the fraction of clear sky irradiance expected and x is the percentage of cloud cover in the sky
    '''
    return 1 - (0.00243 * (cloud_cover*100)) - (4.24 * 10**-5) * ((cloud_cover*100)**2)

def reference_solar_radiation(date, longitude, latitude, shading_proportion, canal_area, cloud_cover):
    '''
    > Calculates solar radiation for reference material (W/m)
    > Similar to ConcreteSolarRadiation() from original code, renamed from concrete_solar_radiation()
    > Utilises the following parameters:
        - Reference material absorbivity (R_ABSORBIVITY) = 1 - albedo
        - Relative light intensity from the shade (SHADE) = 0.2
    '''

    # Determine shaded and unshaded area (m^2)
    shaded_area = canal_area * shading_proportion
    unshaded_area = canal_area * (1 - shading_proportion)

    # Calculate solar intensity, incorporating cloud_cover
    intensity = solar_intensity(date, longitude, latitude)[0] * cloud_cover_effect(cloud_cover) 

    # If positive, calculate solar radiation
    if intensity > 0: 
        solar_rad = (unshaded_area + shaded_area * SHADE) * intensity * R_ABSORBIVITY 

    # Set to null, if intensity is negative:
    else:
        solar_rad = 0

    # Return output
    return solar_rad

def reference_after_sunlight(ref_solar_radiation, reference_temp, canal_area):
    '''
    > Calculates new reference temperature based on sunlight (°C)
    > Similiar to ConcreteAfterSunlight() from original code
    > Utilises the following parameters:
        - Reference depth (R_DEPTH) = 0.05 m
        - Reference density (R_DENSITY)
        - Reference specific heat capacity (R_CAPACITY)
    '''  

    # Update the surface reference temperature (index = 0)
    reference_temp[0] = ((ref_solar_radiation * MODEL_INTERVAL) / (R_CAPACITY * (canal_area * R_DEPTH * R_DENSITY))) + reference_temp[0] 

    # Return list
    return reference_temp

def reference_emissions(reference_temp, canal_area):
    '''
    > Calculates reference material thermal emissions (W/m)
    > Similiar to ConcreteEmissions() from original code, renamed from concrete_emissions()
    > Utilises the following parameters:
        - Reference emissivity (R_EMISSIVITY) = 1 - albedo
        - Stefan-Boltzmann Constant (SB) = 5.67*10**(-8)
    ''' 
    # Modification of surface reference temperature (K)     
    return(canal_area * R_EMISSIVITY * SB * reference_temp[0]**4)

def thermal_absorption(temperature, relative_humidity, cloud_height, cloud_cover, canal_area):
    '''
    > Calculates thermal absorption (longwave) from sky
    > Similiar to ThermalAbsorption() from original code
    - Canal width and length are excluded
    > Utilises:
        - Air temperature (K)
        - Relative humidity (%)
        - Cloud cover (fractional)
        - Cloud height (fractional)
    > Source: Goforth et al. (2002) https://doi.org/10.1117/12.459570 
    '''   

	# Return output (W/m)
    return ((1 + cloud_height * cloud_cover**2) * 8.78 * (10**-13) * (temperature**5.852) * (relative_humidity**0.07195)) * canal_area

def reference_convection_air(reference_temp, air_temperature, wind_speed, canal_area, air_pressure):
    '''
    > Calculates reference convection to air (W/m)
    > This approach is suitable for a solid surface (e.g., asphalt or concrete) under forced convection, but an alternative approach is needed for liquid surfaces
    > Similiar to ConcreteEnergytoAir() from original code, renamed from concrete_convection_air()
    > Utilises the following parameters:
        - Specific heat capacity of the air (A_CAPACITY)
        - Specific gas constant for dry air (A_R_SPECIFIC)
        - Conductivity of air (A_CONDUCTIVITY)
    > Some related papers:
        - Lee et al. (2009): https://doi.org/10.1016/j.cemconcomp.2008.09.009
        - Yang et al. (2024): https://doi.org/10.1038/s41598-024-64568-6
        - Nguyen et al. (2024): https://doi.org/10.1016/j.tsep.2024.102510
    '''

    # Dynamic viscocity (Sutherlands Law, ~1.81 * 10^-5 at 20°C)
    # Source: https://www.cfd-online.com/Wiki/Sutherland's_law
    u = (1.716 * 10**(-5)) * (air_temperature / 273.15)**(3/2) * ((273.15 + 111) / (air_temperature + 111))

    # Kinematic viscocity
    # Dynamic viscocity / density of air [Ideal Gas Law]
    # Air pressure converted from hPa to Pa
    v = u / ((air_pressure * 100) / (A_R_SPECIFIC * air_temperature))

    # Prandtl number (~0.71 at 20°C)
    # Specific heat * dynamic viscocity / conductivity
    Pr = ((A_CAPACITY * 1000) * u) / A_CONDUCTIVITY

    # Convective heat transfer coefficient
    # Simplified version of the Nusselt number formula
    # 0.664 = empirical constant derived from solving the laminar boundary layer equations for a flat plate.
    hc = 0.664 * A_CONDUCTIVITY * Pr**(1/3) * wind_speed**(1/2) * v**(-1/2)

    # Heat transfer coefficient * area * dT
    energy_to_air = hc * canal_area * (reference_temp[0] - air_temperature)

    # Return output (W, J/s)
    return energy_to_air

def reference_heat_transfer(reference_temp, net_radiation, energy_to_air, canal_area):  
    '''
    > Calculates heat transfer for reference material
    > Similiar to ConcreteHeatTransfer() from original code, renamed from concrete_heat_transfer()
    > Utilises the following parameters:
        - Reference thermal conductivity (R_CONDUCTIVITY)
        - Ground thermal conductivity (G_CONDUCTIVITY) = 1 W/mK
        - Reference depth (R_DEPTH) = 0.05 m
    '''
    # List for storing output, defined length
    transfer = [None] * len(reference_temp)

    # Use enumerate() to iterate through the reference temperatures
    for index, value in enumerate(reference_temp):
    
        # For the surface layer
        if index == 0:
            transfer[index] = net_radiation + energy_to_air

        # For sections [0.05-0.10] to [0.20-0.25] (m), use reference conductivity 
        elif index < 5: 

            # Equal to CONDUCTIVITY * area * difference in temperature (current section - higher section) / section depth
            transfer[index] = R_CONDUCTIVITY * canal_area * (value - reference_temp[index - 1]) / R_DEPTH

        # For sections [0.25-0.50] (m) and below, use ground conductivity 
        else: 

            # As above, but using G_CONDUCTIVITY 
            transfer[index] = G_CONDUCTIVITY * canal_area * (value - reference_temp[index - 1]) / R_DEPTH

    # Return list
    return transfer

def reference_energy_change(heat_transfer):
    '''
    > Calculates change in reference energy
    > Similiar to ConcreteEnergyChange() from original code, renamed from concrete_energy_change()
    '''

    # List for storing output, defined length
    energy_change = [0] * len(heat_transfer)

    # Use enumerate() to iterate through the reference temperatures
    for index, value in enumerate(heat_transfer):
    
        # While there a deeper sections remaining:
        try:
            
            # (Subsequent heat value - current heat value) * (3600 / 4)
            energy_change[index] = (heat_transfer[index + 1] - value) * MODEL_INTERVAL

        # There are no deeper sections remaining:
        except IndexError:
            break
    
    # Return output
    return energy_change

def reference_temperature_change(total_energy_change, canal_area):
    '''
    > Calculates temperature change of reference material
    > Similiar to ConcreteTempRise() from original code, renamed from concrete_temperature_change()
    > Utilises the following parameters:
        - Reference depth (R_DEPTH) = 0.05 m
        - Reference density (R_DENSITY)
        - Reference specific heat capacity (R_CAPACITY)
        - Ground specific heat capacity (G_CAPACITY) = 1900 J/kg K
        - Ground density (G_DENSITY) = 1500 kg/m3
    '''
    
    # List for storing output, defined length
    temperature_change = [None] * len(total_energy_change)

    # Use enumerate() to iterate through the reference temperatures
    for index, value in enumerate(total_energy_change):

        # For shallow sections, use reference characteristics
        if index < 5: 
            temperature_change[index] = value / (R_CAPACITY * (canal_area * R_DEPTH * R_DENSITY))
        
        # For deeper sections, use ground characteristics
        else:
            temperature_change[index] = value / (G_CAPACITY * (canal_area * R_DEPTH * G_DENSITY))

    # Return output
    return temperature_change

#------------------------- Water modelling -------------------------#

def model_water_temperature(date, longitude, latitude, shading, water_energy, water_temp, split_depths, canal_area, cloud_cover):
    '''
    > Takes an input water temperature (K) and energy (J) and produces a new water temperature after 15 minutes based upon:
            - Air temperature (K)
            - Humidity (%)
            - Cloud height and cover (fractional)
            - Wind speed (m/s)
            - Shading proportion [0-1]
    > Similiar in scope to NewTemp() from the original code, but individual functions are outlined more fully here
    '''

    # Return the total absorbed light, absorptivity and solar intensity, incorporating cloud cover
    total_absorbed_light, absorptivity, intensity = absorbed_light(date, longitude, latitude, cloud_cover)

    # Absorbed / reflected radiation due to sunlight (W, J/s)
    absorbed_radiation, _ = absorbed_reflected_radiation(shading, total_absorbed_light, intensity, absorptivity, canal_area)

    # Energy gain (W)
    total_sunlight_energy = sunlight_energy(absorbed_radiation, split_depths)

    # Net energy change (W)
    net_energy_change = net_sunlight_energy(total_sunlight_energy, water_energy)

    # Change in water temperature (K)
    water_temperature_change = temperature_change(net_energy_change, canal_area)

    # Final water temperature for model run (K), overwriting earlier list
    water_temperature = [x + y for x, y in zip(water_temp, water_temperature_change)] 

    # Return output
    return water_temperature, absorbed_radiation
    
def solar_absorptivity(solar_position):
    '''
    > This function returns the absorptivity of water given a direct sunlight angle, utilising the Fresnel Equations
    > Same functionality as the absorptivity() function from the original code
    > Utilises the following parameters:
        - Refractive index of water (RIW) = 1.33
        - Refractive index of air (RIA) = 1
    '''
    
    # Zenith angle (radians)
    zenith = radians(90 - solar_position)

	# Fresnel Equations, implemented by H. Mcdonald and S. Taylor
    x = sqrt(1 - (((RIA / RIW) * sin(zenith))**2))
    y1 = RIA * cos(zenith) - RIW * x
    y2 = RIA * cos(zenith) + RIW * x
    Rs = (y1 / y2)**2
    z1 = RIA * x - RIW * cos(zenith)
    z2 = RIA * x + RIW * cos(zenith)
    Rp = (z1 / z2)**2
    reflectivity = (Rs + Rp) / 2
    absorptivity = 1 - reflectivity
    
    # Return output
    return absorptivity

def absorbed_light(date, longitude, latitude, cloud_cover):
    '''
    > This function returns the total absorbed light (W/m^2), and the solar absorptivity
    > Used later in absorbed_reflected_radiation()
    > Same functionality as the AbsorbedLight() function from the original code, but utilising outputs from earlier functions
    '''

    # Calculate the overhead sun intensity and the solar altitude (°)
    intensity, solar_altitude = solar_intensity(date, longitude, latitude) 

    # Account for cloud cover
    intensity = intensity * cloud_cover_effect(cloud_cover) 

	# if intensity (W/m^2) is > 0 (day) 
    if intensity > 0:
        
        # Calculate the amount of absorbed light (W/m^2) 
        absorptivity = solar_absorptivity(solar_altitude)
        total_absorbed_light = intensity * absorptivity

    # If the intensity is null or negative (night)        
    else:

        # Sets the total amount of absorbed light to 0, and calculates absorptivity
        total_absorbed_light, absorptivity = 0, solar_absorptivity(solar_altitude)

    # Return output
    return total_absorbed_light, absorptivity, intensity

def absorbed_reflected_radiation(shading_proportion, total_absorbed_light, intensity, absorptivity, canal_area):
    '''
    > This function returns the total absorbed radiation (W/m^2) and total reflected radiation (W/m^2)
    > Same functionality as the absorbed_radiation() and reflected_radiation() functions from the original code...
    > ... but combined for simplicity and  utilising outputs from earlier functions
    > Utilises the following parameter:
        - Relative light intensity from the shade (SHADE) = 0.2
    '''  

    # Determine shaded and unshaded area (m^2)
    shaded_area = canal_area * shading_proportion
    unshaded_area = canal_area * (1 - shading_proportion)

    # Absorbed and reflected radiation calculations
    absorbed_radiation = (unshaded_area + shaded_area * SHADE) * total_absorbed_light 
    reflected_radiation = (unshaded_area + shaded_area * SHADE) * intensity * (1 - absorptivity)
    
    # Return the absorbed radiation and reflected radiation (W)
    return absorbed_radiation, reflected_radiation

def sunlight_energy(absorbed_radiation, split_depths):
    '''
    > This function returns the energy gained by sunlight (J)
    > Similar functionality as the SunlightEnergy() function from the original code
    > Utilises the following parameters:
        - Vertical extinction coefficient (VEC) = 3 (% of surface light absorbed or scattered in a 1 m long vertical column of water)
        - Canal depth intervals (INTERVAL) = 0.2 m (decimal)
    ''' 
    # List for storing output, defined length
    total_sunlight_energy = [None] * len(split_depths)

    # Use enumerate() to iterate through the  water depths
    for index, depth in enumerate(split_depths):

        # While there a deeper sections remaining:
        try:

            # Calculate total energy
            total_sunlight_energy[index] = (absorbed_radiation * exp(-VEC * (depth - INTERVAL)) - absorbed_radiation * exp(-VEC * depth)) * MODEL_INTERVAL
        
        # There are no deeper sections remaining:
        except IndexError:
            break

    # Return output
    return total_sunlight_energy

def net_sunlight_energy(total_sunlight_energy, water_energy):
    '''
    > This function returns the net energy gained by sunlight (W, J/s), incorporating the starting energy amount
    > Similar functionality as the NetSunlightEnergy() function from the original code
    ''' 

    # List for storing output, defined length
    net_energy_change = [None] * len(total_sunlight_energy)

    # Use enumerate() to iterate through the energy values
    for index, value in enumerate(total_sunlight_energy):

        # While there a deeper sections remaining:
        try:

            # Updates list value
            net_energy_change[index] = value - water_energy[index] + water_energy[index + 1]

        # There are no deeper sections remaining:
        except IndexError:

             # Updates list value
            net_energy_change[index] = value - water_energy[index]
  
    # Retunrs the output
    return net_energy_change

def temperature_change(net_energy_change, canal_area):
    '''
    > This function returns the temperature change, based on the net sunlight energy and the canal geometry
    > Similar functionality as the TempRise() function from the original code
    > Utilises the following parameters:
        - Specific heat capacity of water (SHC_WATER) = 4186 J/kg K
        - Canal depth intervals (FLOAT_INTERVAL) = 0.2 m (float)
    ''' 

    # List for storing output, defined length
    total_temperature_change = [None] * len(net_energy_change)

    # Iterate through the net energy values
    for index, energy in enumerate(net_energy_change):
            
            # Calculate temperature change = energy / (mass of water (kg) * specific heat capacity)
            total_temperature_change[index] = energy / (((canal_area * FLOAT_INTERVAL) * 1000) * SHC_WATER)

    # Return dictionary (°C)
    return total_temperature_change    

def model_energy_change(air_temperature, relative_humidity, cloud_height, cloud_cover, wind_speed, air_pressure,
                        water_temperature, canal_area, initial_water_k, absorbed_radiation):
    '''
    > Takes an input water temperature (K) and energy (J) and produces a new water energy (J) after 15 minutes
    > Similiar in scope to HeatTransfer() from the original code, but individual functions are outlined more fully here
    '''
     
    # Model net thermal emissions, based on surface water temperature (W, J/s)
    modelled_thermal_emissions = net_thermal_emissions(air_temperature, relative_humidity, cloud_height, cloud_cover, initial_water_k, canal_area) 
    
    # Models the evaporation rate (m/d) and vapour pressures (hPa) after Webb and Zhang (1998)
    modelled_evaporation_rate, modelled_saturation_vapour_pressure, modelled_vapour_pressure = evaporation_rate(air_temperature, wind_speed, relative_humidity, initial_water_k)

    # Modelled energy flux due to evaporation (J/m^2/day) 
    modelled_evaporation_energy = evaporation_energy(modelled_evaporation_rate, air_temperature, initial_water_k)

    # Modelled evaporation power (W, J/s), accounting for canal area
    modelled_evaporation_power = evaporation_power(modelled_evaporation_energy, canal_area)
        
    # Sensible heat calculations
    modelled_sensible_heat = sensible_heat(air_temperature, air_pressure, initial_water_k, wind_speed, canal_area)

    # Energy transfer for the surface layer to / from the air (J / 15 mins), the sum of the sensible and the latent heat fluxes, and net longwave radiation (emissions - thermal absorption)
    modelled_surface_heat_transfer = surface_heat_transfer(modelled_evaporation_power, modelled_thermal_emissions, modelled_sensible_heat) 

    # Energy for all layers (J / 15 mins)
    modelled_energy = canal_heat_transfer_convection(modelled_surface_heat_transfer, water_temperature, canal_area)
    
    # Created dictionary to return key variables (J / 15 mins)
    energy_values = {'sensible_heat' : modelled_sensible_heat * MODEL_INTERVAL,
                   'net_thermal_emissions' : modelled_thermal_emissions * MODEL_INTERVAL}

    # Returns output
    return modelled_energy, energy_values

def model_air_density(air_temperature, air_pressure):
    '''
    > A function to calculate air density (kg/m³) using air temperature (K) and air pressure (Pa)
    > Parameters:
        - Specific gas constant for dry air (A_R_SPECIFIC) = 287.05 J/kg·K
    '''

    # Modelled air density, converting air pressure from hPa to Pa
    return (air_pressure * 100) / (A_R_SPECIFIC * air_temperature)

def model_water_density(water_temperature):
    '''
    > A function to calculate water density (kg/m³) using water temperature (K)
    > Source: Kell (1975) cited in Jones and Harris (1992): https://pmc.ncbi.nlm.nih.gov/articles/PMC4909168/
    '''

    # Convert water K to °C
    wt_celcius = water_temperature - 273.15

    # Return modelled water density, following Kell (1975)
    return (999.83952 
            + 16.945176 * wt_celcius 
            - 7.9870401 * 10**(-3) * wt_celcius**2 
            - 46.170461 * 10**(-6) * wt_celcius**3 
            + 105.56302 * 10**(-9) * wt_celcius**4 
            - 280.54253 * 10**(-12) * wt_celcius**5) / (1 + 16.897850 * 10**(-3) * wt_celcius)

def water_emissions(water_temperature, canal_area):
    '''
    > Function returns the radiative emissions from water (W/m^2), based on the water temperature (K)
    > Similiar functionality as the WaterEmmissions() function from the original code
    > Uses the following parameters:
        - Stefan-Boltzmann Constant (SB) = 5.67*10**(-8)
        - Emissivity of water (E) = 0.95
    '''
    # Return output (W/m) = εσT4 (https://doi.org/10.1016/B978-0-12-083980-3.50013-9)
    return (canal_area * E * SB * water_temperature**4)

def net_thermal_emissions(temperature, relative_humidity, cloud_height, cloud_cover, water_temperature, canal_area):
    '''
    > Net thermal emissions (emissions - thermal absorption)
    > Identical to NetThermalEmissions() from original code
    '''

    # Return output (W)
    return water_emissions(water_temperature, canal_area) - thermal_absorption(temperature, relative_humidity, cloud_height, cloud_cover, canal_area)

def saturation_vapour_pressure(temperature): 
    '''
    > This function returns saturation vapour pressure (hPa), utilising the temperature in K
    > Similar functionality as the Ew() function from the original code, renamed from surface_evaporation_rate()
    > For modelling the evaporation_rate(), this utilises vapour pressure at the water surface, following Webb and Zhang (1998)
    ''' 

    # Saturation vapour pressure (Pa)
    modelled_svp = exp(77.345 + 0.0057 * (temperature) - 7235 / (temperature)) / pow((temperature), 8.2)

    # Return output, converting Pa to hPa (/100)
    return modelled_svp / 100

def evaporation_rate(air_temperature, wind_speed, relative_humidity, initial_water_k):
    '''
    > This function returns the evaporation rate (m/day), plus the saturation vapour pressure and actual vapour pressure
    > Identical functionality as the Ev() function from the original code, but simplified
    '''

    # Saturation and actual vapour pressure (hPa)
    modelled_vapour_pressure = vapour_pressure(air_temperature, relative_humidity)
    modelled_saturation_vapour_pressure = saturation_vapour_pressure(initial_water_k)

    # Convert wind speed from 10 m to 2 m height, assuming neutral conditions (1/7)
    # Source: Touma (1977) https://doi.org/10.1080/00022470.1977.10470503
    wind_speed_2m = wind_speed * pow(2 / 10, 1 / 7)

    # Evaporation rate (mm/day) after Webb and Zhang (1998), Equation 6 
    # An alternative approach would be to use the Penman Equation, which also incorporates radiative inputs
    modelled_evaporation_rate = 0.165 * (0.8 + (0.864 * wind_speed_2m)) * (modelled_saturation_vapour_pressure - modelled_vapour_pressure)

    # Return:
    # [1] Evaporation rate in m/d (/1000)
    # [2] Saturation vapour pressure at the surface water temperature (hPa) 
    # [3] Vapour pressure at air temperature (hPa)
    return modelled_evaporation_rate / 1000, modelled_saturation_vapour_pressure, modelled_vapour_pressure

def vapour_pressure(air_temperature, relative_humidity):   
    '''
    > This function returns air vapour pressure (hPa, mbar), utilising the following climate parameters
        - Air temperature (K)
        - Relative humidity
    > The function also utilises the following parameters:
        - Constant to calculate dew point (DEW_POINT_B) = 17.67
        - Constant to calculate dew point (DEW_POINT_C) = 243.5
        - Ratio of latent heat to water vapour gas constant (LRV) = (L*1000)/RV  
        - Baseline Pressure at Triple Point (BASE_PRESSURE) = 0.611657 kPa
        - Triple point of water (TPW) = 273.16 K
    > Identical functionality as the Ea() function from the original code, but with clearer signposting of constants
    '''

    # Calculate Magnus-Tetens approximation for dew point temperature (converting between to K and °C)
    y = log(relative_humidity * 0.01) + DEW_POINT_B * (air_temperature - 273.15) / (DEW_POINT_C + (air_temperature - 273.15))
    dpt = DEW_POINT_C * y / (DEW_POINT_B - y) + 273.15

    # Vapour pressure at air temperature (kPa, convert to hPa / mbar)
    vapour_pressure = BASE_PRESSURE * (exp(LRV * ((1.0 / TPW) - (1 / dpt)))) * 10

    # Return output
    return vapour_pressure

def evaporation_energy(modelled_evaporation_rate, air_temperature, water_temperature):
    '''
    > This function returns the amount of heat lost by evaporation or gained through condensation for the water body (J/m^2/day) 
    > Source: Webb and Zhang (1998), Equation 7
    > The function also utilises the following parameters:
        - Latent heat of vapourisation (L) = 2454.9 J/g
    > Identical functionality as the EvaporationEnergy() function from the original code, but simplified
    '''

    # Latent heat of vaporisation, as a function of temperature after Webb and Zhang (1998), Equation 8
    modified_L = L - (2.366 * (air_temperature - 273.15))

    # Calculation of water density (kg/m3)
    water_density = model_water_density(water_temperature) 

	# Return evaporative flux (W m/2), using water density in g m^3
    return modelled_evaporation_rate * modified_L * (water_density * 1000)
    
def evaporation_power(modelled_evaporation_energy, canal_area):
    '''
    > This function returns the power loss to evaporation (W/m)
    > Similiar functionality as the EvaporationPower() function from the original code
    '''

    # Power loss to evaporation for entire area (W, J/s) per second [86,400 seconds per day]
    return modelled_evaporation_energy * canal_area * (1 / 86400)  

def sensible_heat(air_temperature, air_pressure, water_temperature, wind_speed, canal_area):
    '''
    > Calculates sensible heat transfer (W, J/s) using water/air temperatures (K), air pressure (Pa) and wind speed (m/s)
    > Similiar functionality to SensibleHeat() from the original code
    > Approach following Moore and Leach (2021) https://doi.org/10.1029/2020WR028712
    '''   

    # Function to calculate air density using temperature and pressure
    air_density = model_air_density(air_temperature, air_pressure)

    # Empirical constants calculated by Moore and Leach (2021), for a reference temperature of 15 ◦C, and using Penman wind function coefficient values
    # "Varying air temperature from 0 to 30 ◦C is associated with a variation of just over ± 5% in the computed coefficients relative to the reference values at 15 ◦C"
    # "...temperature-dependence can be ignored for practical application."
    a = 1.99 * 10**(-3)
    b = 2.13 * 10**(-3)

    # Density of air (kg/m3), Specific heat (J·kg*k), constants, wind speed (m/s), temperatures (K), area (m2)
    Qh = air_density * (A_CAPACITY * 1000) * (a + b * wind_speed) * (water_temperature - air_temperature) * canal_area

    # Return sensible heat flux (W, J/s)
    return Qh

def surface_heat_transfer(latent_flux, radiative_flux, sensible_flux):     
    '''
    > Calculates heat transfer from the water surface, incorporating:
        - latent heat flux, energy used for evaporation or condensation
        - sensible heat flux, heat transfer of heat through conduction and convection, without any change in state, based on temperature difference
        - radiative heat flux, heat transfer via radiation
    > Similiar functionality to SurfaceHeatTransfer() from the original code
    '''   

    # Return output, converting from /s to /15m
    return (latent_flux + radiative_flux + sensible_flux) * MODEL_INTERVAL 

def canal_heat_transfer(modelled_surface_heat_transfer, water_temperature, canal_area):
    '''
    > Calculates heat transfer between depths (split_points), incorporating conductive-convective effects
    > Similiar to HeatTransfer() from original code
    > Uses the following parameters:
        - Heat conductivity of water (HC_WATER) =  0.58 W/m K
        - Convection coefficent (CONVECTION) = 200 hc
        - SPLIT_DEPTH = the depth (m) of each modelled canal section
    '''   

    # List for storing output, defined length
    modelled_heat_transfer = [None] * len(water_temperature)

    # Iterate through the temperature values
    for index, current_temperature in enumerate(water_temperature):

        # For the surface layer:
        if index == 0:
            
            # Update with the modelled surface heat transfer
            modelled_heat_transfer[index] = modelled_surface_heat_transfer

        # For any deeper layers
        else:

            # Obtain the preceding temperature
            higher_depth = water_temperature[index - 1]
      
            # Calculate conductivity
            conductivity = HC_WATER * canal_area * (current_temperature - higher_depth) / FLOAT_INTERVAL

            # If the temperature of the current depth is greater than the temperature above:
            if current_temperature > higher_depth: 

                # Energy value (/15m)
                modelled_heat_transfer[index] = (conductivity + CONVECTION * canal_area * (current_temperature - higher_depth)) * MODEL_INTERVAL
                
            # If temperature <
            else:
                
                # Energy value (/15m)
                modelled_heat_transfer[index] = conductivity * MODEL_INTERVAL

    # Return energy values
    return modelled_heat_transfer


def canal_heat_transfer_convection(modelled_surface_heat_transfer, water_temperature, canal_area):
    '''
    > Calculates heat transfer between depths (split_points), incorporating conductive-convective effects
    > Similiar to HeatTransfer() from original code
    > Adaped from canal_heat_transfer() which is solely based on the difference in temperature between the layers:
        - When temperature in a lower layer is greater than the temperature in the layer above it, there is an energy flow upwards
        - When temperature in a lower layer is less than the layer above, there is no energy flow
        - However, this does not work at or below 4 degrees, where the process is reversed
    > This new function updates this as follows:
        - The convective effect is dependent on absolute distance from 4 degrees (277.15 K)
        - i.e., 2 degrees and 6 degrees would look the same in terms of convection
        - This is a simplistic assumption (density is not directly proportional) but is an improvement on the previous approach
        - Model performance will still degrade once temperature reaches 0 degrees where latent heat and freezing effects would need to be implemented
    > Uses the following parameters:
        - Heat conductivity of water (HC_WATER) =  0.58 W/m K
        - Convection coefficent (CONVECTION) = 200 hc
        - SPLIT_DEPTH = the depth (m) of each modelled canal section
    '''   

    # List for storing output, defined length
    modelled_heat_transfer = [None] * len(water_temperature)

    # Iterate through the temperature values
    for index, current_temperature in enumerate(water_temperature):

        # For the surface layer:
        if index == 0:
            
            # Update with the modelled surface heat transfer
            modelled_heat_transfer[index] = modelled_surface_heat_transfer

        # For any deeper layers
        else:

            # Obtain the preceding temperature
            higher_depth = water_temperature[index - 1]

            # Absolute difference from 4°C (maximum density)
            upper_temp = abs(higher_depth - 277.15)
            lower_temp = abs(current_temperature - 277.15)

            # Calculate conductivity
            conductivity = HC_WATER * canal_area * (current_temperature - higher_depth) / FLOAT_INTERVAL
            
            # If difference from 4°C in lower layer > difference in the layer above
            if lower_temp > upper_temp:

                # If temperature in lower layer > temperature in the layer above
                if current_temperature > higher_depth:

                    # Energy value (/15m)
                    modelled_heat_transfer[index] = (conductivity + CONVECTION * canal_area * (lower_temp - upper_temp)) * MODEL_INTERVAL

                # If temperature in lower layer <= temperature in the layer above, -CONVECTION
                else:
                    # Energy value (/15m)           
                    modelled_heat_transfer[index] = (conductivity + -CONVECTION * canal_area * (lower_temp - upper_temp)) * MODEL_INTERVAL

            # If difference from 4°C in lower layer <= difference in the layer above
            else:

                # Energy value (/15m)
                modelled_heat_transfer[index] = conductivity * MODEL_INTERVAL

    # Return energy values
    return modelled_heat_transfer


def model_spin_up(water_temp, water_energy, reference_temp, 
                  split_points, canal_area, lat, lon, 
                  spin_path, repetitions, record, canal_id):
    '''
    > Function for model spin up, based on a composite climate record for 2021-12-15
    > Using input temperatures (K) and energy values (J), the model is run until equilibrium is reached
    '''

    # Open the spin up record
    with open(spin_path) as path:
            spin_up_climate = load(path)

    # Init lists for storing outputs
    output_water = []
    output_energy = []
    output_reference = []

    # List of keys
    keys = [int(x) for x in spin_up_climate.keys()]

    # Minimum (start date), as unix and datetime
    start = min(keys)
    start_date = dt.fromtimestamp(start)

    # Iterate n times
    for _ in range(repetitions):

        # Sets start date
        local_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
    
        # Sets end date (one day)
        end_datetime = local_datetime + timedelta(hours = 24)

        # Iterate until we've done 24 h
        while local_datetime < end_datetime:

            # Convert to unix, as string
            start_unix = str(int(local_datetime.timestamp()))

            # Starting surface water temperature (K) for this model step
            initial_water_k = water_temp[0]

            # Try and extract current climate
            try:
                current_climate = spin_up_climate[start_unix]

            # Reached the end of the record
            except KeyError:

                # Skip to next iteration
                local_datetime = end_datetime

            # Valid datetime
            else:

                # Model water temperatures
                water_temp, absorbed_radiation = model_water_temperature(local_datetime, #-------- UTC datetime
                                                        lon, lat, #------------------------------- Geographic coordinates (WGS84)
                                                        0, #-------------------------------------- Shading proportion [0-1] 
                                                        water_energy, #--------------------------- Water energy (J)
                                                        water_temp, #----------------------------- Water temperature (K)
                                                        split_points, #--------------------------- Canal depths (m)
                                                        canal_area, #----------------------------- Canal area (m^2)
                                                        current_climate['cloud_cover'] / 100) #--- Cloud cover (fractional)
                
                # Model water energy
                water_energy, _ = model_energy_change(current_climate['air_temperature'], #------- Air temperature (K)
                                    current_climate['humidity'], #-------------------------------- Relative humidity (%)
                                    CLOUD_HEIGHT, #----------------------------------------------- Cloud height (fractional)
                                    current_climate['cloud_cover'] / 100, #----------------------- Cloud cover (fractional),
                                    current_climate['wind_speed'], #------------------------------ Wind speed (m/s)
                                    current_climate['air_pressure'], #---------------------------- Air pressure (hPa)
                                    water_temp, #------------------------------------------------- Water temperature (K)
                                    canal_area, #------------------------------------------------- Canal area (m^2)
                                    initial_water_k, #-------------------------------------------- Starting surface water temperature (K) for this model step 
                                    absorbed_radiation) #----------------------------------------- Radiation (W)
                

                # Model reference temperature
                reference_temp, _ = model_reference_temperature(current_climate['air_temperature'], #--- Air temperature (K)
                                                current_climate['humidity'], #-------------------------- Relative humidity (%)
                                                CLOUD_HEIGHT, #----------------------------------------- Cloud height (fractional)
                                                current_climate['cloud_cover'] / 100, #----------------- Cloud cover (fractional)
                                                current_climate['wind_speed'], #------------------------ Wind speed (m/s)
                                                current_climate['air_pressure'], #---------------------- Air pressure (hPa)
                                                reference_temp, #--------------------------------------- Reference temperature (K)
                                                local_datetime, #--------------------------------------- UTC datetime
                                                lon, lat, #--------------------------------------------- Geographic coordinates (WGS84)
                                                0, #---------------------------------------------------- Shading proportion [0-1] 
                                                canal_area) #------------------------------------------- Area (m^2)

                # Append surface temperatures and energy to output
                output_water.append(water_temp[0])
                output_energy.append(water_energy[0])
                output_reference.append(reference_temp[0])
                
                # Add 15 minutes to the datetime
                local_datetime += timedelta(hours = 0.25)

    # Write to output file
    if record:

        # Zip lists and export
        output = array(list(zip(output_water, output_energy, output_reference)))
        savetxt(f"../outputs/spin_up_output_{canal_id}.csv", output, delimiter = ",", header="water_k, water_j, reference_k", fmt ='%f')

    # Return modified temperatures (K) and energy (J), and time-series values
    return water_temp, water_energy, reference_temp, output_water, output_energy, output_reference

