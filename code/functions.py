''' Usage: This code contains functions used in national-canal-cooling.py and will not produce results if run in isolation.
    Not all functions are utilised but are included here for posterity e.g.,  
    > return_buildings() superseeded by return_buildings_polygon() in auxiliary.py
    > filter_buildings() superseeded by subset_buildings()
    > return_shading() superseeded by teturn_shading() in auxiliary.py
    > determine_shading() superseeded by determine_shading_proportion()'''

# Force use of Shapely 2.0
from os import environ
environ['USE_PYGEOS'] = '0'

# Required packages
import warnings
import requests
import pybdshadow
from pyproj import Geod
from json import loads, load
from pandas import to_datetime
from datetime import timedelta
from suncalc import get_position
from geoalchemy2 import Geometry
from datetime import datetime as dt
from pysolar.solar import get_altitude
from sqlalchemy import create_engine, text
from numpy import savetxt, array, linspace
from geopandas import GeoDataFrame, GeoSeries
from matplotlib.pyplot import subplots, savefig, title
from shapely.geometry import Point, Polygon, LineString
from math import radians, cos, degrees, sqrt, sin, exp, log

# Import user-defined parameters
import config
from params import * 

#---- Key variables ----#

# PostGIS password
PW = config.info['postgis']

# String for connection
DB_STRING = f"postgresql://postgres:{PW}@localhost:5432/canals"

# API key
MT_API = config.info['api_mt']
JH_API = config.info['api_jh']

# Sets ellipsoid model
g_model = Geod(ellps="WGS84")

# Prevent geopandas CRS warning for pybdshadow (Context: https://github.com/geopandas/geopandas/issues/2606)
# shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)
warnings.filterwarnings("ignore", message="CRS not set for some of the concatenation inputs")

#---- Geospatial Functions (buildings, shading, climate data) ----#

def return_buildings(longitude, latitude, distance, plot):
    '''
    > Returns buildings (stored in postGIS) within a specified {distance} from a location {longitude, latitude}
    '''

    # Create the engine
    db_connection = create_engine(DB_STRING) 

    # Convert to GeoSeries and project
    point_geometry_geographic = GeoSeries(Point(longitude, latitude), crs = 4326)
    point_geometry_projected = point_geometry_geographic.to_crs(27700)

    query = text("select * from buildings where ST_DWithin(buildings.geometry, ST_SetSRID(ST_MakePoint(:e, :n), 27700), :d)") \
    .bindparams(n = point_geometry_projected.iloc[0].y, # Northing
                e = point_geometry_projected.iloc[0].x, # Easting
                d = distance) # Threshold distance

    # Extracts geometries using from_postgis, stores as geodataframe
    buildings_via_python = GeoDataFrame.from_postgis(sql = query, 
    con = db_connection,
    geom_col='geometry', # Column name for the geometry
    index_col='os_topo_toid', # Column name for the unique ID
    coerce_float=True)

    # Remove date-time field
    buildings_via_python = buildings_via_python.drop(columns='bha_processdate')

    # Plot the output
    if plot == True:

        # Set up output image
        fig, my_ax = subplots(1, 1, figsize=(16, 10))
        title(f"Buildings-water point test for {round(latitude, 4), round(longitude, 4)}")

        # Plot the building geometries (spatial query)
        buildings_via_python.plot(
            ax = my_ax,
            color = "#343434",
            edgecolor = '#343434',
            linewidth = 0.5,
            )

        # Add north arrow
        x, y, arrow_length = 0.97, 0.99, 0.1
        my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
            arrowprops=dict(facecolor='black', width=5, headwidth=15),
            ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

        # Save the result
        savefig(f'./images/buildings-at-{round(latitude, 4)}_{round(longitude, 4)}.png', bbox_inches='tight')

    # Return geometries of buildings
    return buildings_via_python

def filter_buildings(date, building_data, longitude, latitude, polygon, tolerance, distance, interval):
    '''
    > A new function to improve efficiency
    > return_buildings() returns all buildings within a specified distance (m) of a location
    > Building geometries are then used to calculate shading using return_shading() but this could be more efficient
    > However, only a small number of those buildings could shade the studied location given the sun's azimuth (direction)
    > This function returns those buildings based on:
        - date, latitude, longitude = for calculating the solar azimuth 
        - tolerance (°) = threshold (±) for determining whether buildings are in the current solar direction
    > get_position()['azimuth'] returns:
        - solar azimuth in radians, measured relative to south (east = negative, west = positive)
        - This is converted to the more conventional direction
        - These results are consistent with https://www.sunearthtools.com/dp/tools/pos_sun.php
    > Rather that checking the azmiuth for each building centroid or coordinate sequence...
    > We simply generate the desired area, and then perform an intersects
    '''

    # Returns the solar azimuth in radians
    solar_azimuth = get_position(date, longitude, latitude)['azimuth']

    # Convert azmiuth from relative to south to relative to north (degrees)
    if solar_azimuth < 0: 
        solar_azimuth = degrees(radians(180) - abs(solar_azimuth))
    else:
        solar_azimuth = degrees(radians(180) + abs(solar_azimuth))

    # Line-based analysis
    if polygon == False:

        # Coordiantes of shape, based on azimuth ± tolerance, and specified distance and interval
        forward_long, forward_lat = g_model.fwd([longitude]*2,
                                        [latitude]*2, 
                                        az = [solar_azimuth]*2,
                                        dist = [distance, 0], 
                                        radians = False)[0:2]
        
        forward_linestring = GeoDataFrame(index=[0], crs = 'epsg:4326', geometry = [LineString(zip(forward_long, forward_lat))])\
        .to_crs(building_data.crs)
        
        # Search for buildings using intersects
        filtered_buildings = building_data[building_data.intersects(forward_linestring.geometry.iloc[0])]

        # Return filtered buildings (potential shading)
        return filtered_buildings, forward_linestring
    
    # Polygon-based analysis
    else: 

        # Generate sequence of azimuths and convert to list
        azimuths = linspace(start = solar_azimuth - tolerance, 
                            stop = solar_azimuth + tolerance,
                            num = interval, 
                            endpoint = True).tolist()

        # Sequence of distances
        distances = [distance] * interval

        # Insert (in place) zero at the start of each list, representing the origin location
        azimuths.insert(0,0)
        distances.insert(0,0)

        # Coordiantes of shape, based on azimuth ± tolerance, and specified distance and interval
        forward_long, forward_lat = g_model.fwd([longitude] * (interval + 1),
                                        [latitude] * (interval + 1), 
                                        az = azimuths,
                                        dist = distances, 
                                        radians = False)[0:2]
        

        # Creates polygon from coordinates, converts to GDF, reprojects to BNG
        forward_polygon = GeoDataFrame(index=[0], crs = 'epsg:4326', geometry = [Polygon(zip(forward_long, forward_lat))])\
            .to_crs(building_data.crs)

        # Search for buildings using intersects
        filtered_buildings = building_data[building_data.intersects(forward_polygon.geometry.iloc[0])]
    
    # Return filtered buildings (potential shading)
    return filtered_buildings, forward_polygon

def return_shading(date, filtered_buildings, extracted_buildings, longitude, latitude, forward_shape, plot):
    '''
    > This function returns the shading geometry from buildings, for a specified datetime (UTC)
    > Utilises 'pybdshadow': https://pybdshadow.readthedocs.io/en/latest/
    > Results verified against: https://shadowmap.org/
    > Building geometries (OS) obtained via return_buildings()
    '''

    # Return the solar altitude, converted to degrees
    solar_position = degrees(get_position(date, longitude, latitude)['altitude'])

    # If the sun is above the horizon
    if solar_position < 0:

            return None, solar_position
    
    # The sun is above the horizon, but there are no buildings in the solar direction
    if filtered_buildings.empty:

            return None, solar_position
       
    # The sun is above the horizon AND there are buildings that could shade the location
    else: 
        
        # Keep unique ID, building height, and geometry
        building_data = filtered_buildings[['relhmax', 'geometry']]

        # Rename column to work with pybdshadow
        building_data = building_data.copy().rename(columns={'relhmax': 'height'})

        # Project to geographic coordinates
        building_data_projected = building_data.to_crs('epsg:4326')

        # This part is slow and results in the following warning:

            # > UserWarning: CRS not set for some of the concatenation inputs. \
            # > Setting output's CRS as WGS 84 (the single non-null crs provided)."

        # Preprocessing (remove empty polygons, multipolygons > polygons, generate building_id)
        #building_data_projected = pybdshadow.bd_preprocess(building_data_projected)

        # Explode multi-part geometries into multiple single geometries.
        building_data_projected = building_data_projected.explode(index_parts = False)

        # Add 'building_id' column
        building_data_projected['building_id'] = range(1, len(building_data_projected) + 1)

        # Calculate building shadows
        shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)

        # Assigns a CRS to the output (4326)
        shadows.crs = building_data_projected.crs

        # Reproject to British National Grid
        shadows = shadows.to_crs(building_data.crs)

        # Plotting
        if plot == True:

            # Creates Point from coordinates, converts to GDF, reprojects to BNG
            studied_location = GeoDataFrame(index=[0], crs = 'epsg:4326', geometry = [Point(longitude, latitude)])\
            .to_crs(building_data.crs)

            # Set up output image
            fig, my_ax = subplots(1, 1, figsize=(16, 10))
            title(f"Buildings-shadow test at (UTC) {date}")

            # Plot the building geometries (spatial query)
            shadows.plot(
                ax = my_ax,
                color = "#959595",
                edgecolor = None,
                linewidth = 0.5,
                )
            
            # Plot the extracted building geometries
            extracted_buildings.plot(
                ax = my_ax,
                color = "#FFFFFF",
                edgecolor = "#1189FF",
                linewidth = 0.5,
                )

            # Plot the filtered building geometries
            building_data.plot(
                ax = my_ax,
                color = "#1189FF",
                edgecolor = None,
                linewidth = 0.5,
                )
            
            # Plot the solar azimuth line
            forward_shape.plot(
                ax = my_ax,
                color = "#EA2F2F",
                edgecolor = None,
                linewidth = 0.5,
                )

            # Plot the studied location
            studied_location.plot(
                ax = my_ax,
                color = "#EA2F2F",
                edgecolor = None,
                linewidth = 0.5,
                )

            # Add north arrow
            x, y, arrow_length = 0.97, 0.99, 0.1
            my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

            # Save the result
            savefig(f'./images/shadows-at-{date.date().strftime("%m_%d")}-{date.time().strftime("%H_%M_%S")}.png', bbox_inches='tight')
    
        # Return output
        return shadows, solar_position

def determine_shading(shaded_area, longitude, latitude):
    '''
    > This function determines if the studied location is shaded (True) or not (False) for a specified datetime
    '''

    # If there is a shaded area (store as gdf)
    if isinstance(shaded_area, GeoDataFrame):

        # Creates Point from coordinates, converts to GDF, reprojects to BNG
        studied_location = GeoDataFrame(index=[0], crs = 'epsg:4326', geometry = [Point(longitude, latitude)])\
            .to_crs(shaded_area.crs)

        # Check for intersect between shaded area and current location
        if any(shaded_area.intersects(studied_location.geometry.iloc[0])):
        
            return True # There is a shaded area and it intersects with the studied location
        
        else: 

            return False # There is a shaded area but it does not intersect with the studied location

    else:

        return False # There is no shaded area (sun altitude below horizon)
    
# Obtain climate data from OpenWeatherMap
def get_climate_data(start_datetime, longitude, latitude, API):
    '''
    > This function returns climate data from OpenWeatherMap for a specified datetime (UTC) and location
    > To minimise the number of API calls, the function first checks a postGIS database
    > If the data is absent, the API is used and the dataset is then uploaded to postGIS
    > Function not used, as weather data were bulk downloaded from OpenWeather
    '''

    # First, check for climate data in the postGIS table, establish the connection
    db_connection = create_engine(DB_STRING)

    # Format query, using latitude, longitude and datetime as inputs
    query = text("select * from climate_data where latitude = :a and longitude = :b and datetime = :d") \
    .bindparams(a = latitude, b = longitude, d = start_datetime)

    # Extracts climate data, stores as geodataframe
    weather_gdf = GeoDataFrame.from_postgis(sql = query, 
    con = db_connection,
    geom_col='geometry', # Column name for the geometry
    coerce_float=True)

    # If there is no data, pull the data from the OpenWeatherAPI and add to postGIS
    if weather_gdf.empty: 

        # Sets end datatime (+1 hr)
        end_datetime = start_datetime + timedelta(hours=1)
    
        # Converts dates to Unix timestamp
        start_time = start_datetime.timestamp()
        end_time= end_datetime.timestamp()

        # API call to Open Weather (formatting for readability)
        response = requests.get(f"http://history.openweathermap.org/data/2.5/history/city?&" \
                                f"lat={latitude}&" \
                                f"lon={longitude}&" \
                                f"type=hour&" \
                                f"start={start_time}&" \
                                f"end={end_time}&" \
                                f"appid={API}")
        
        # Loads as json
        weather = loads(response.text)

        # Creating a geodataframe from input latitude and longitude, including those as df columns
        weather_gdf = GeoDataFrame({'latitude': latitude, 'longitude' : longitude, 
            'datetime' : start_datetime,
            'air_temperature' : weather['list'][0]['main']['temp'], 
            'humidity' : weather['list'][0]['main']['humidity'], 
            'cloud_cover' : weather['list'][0]['clouds']['all'], 
            'wind_speed' : weather['list'][0]['wind']['speed'],
            'air_pressure' : weather['list'][0]['main']['pressure'],
            'geometry': [Point(longitude, latitude)]}, crs="epsg:4326").to_crs(27700)
        
        # Add location to table
        weather_gdf.to_postgis('climate_data', db_connection, if_exists='append', dtype={'geom': Geometry('[Point]', srid=27700)})
        
        # Return statement
        return(weather_gdf)

    # The data was previously stored on postGIS
    else:

        # Return gdf
        return weather_gdf

#---- Concrete modelling ----#
def solar_intensity(date, longitude, latitude):
    '''
    > This function returns the intensity of the sun based on its azimuth angle (W/m^2) for a specific location and datetime
    > It also returns the solar altitude angle (°), which is used in solar_absorptivity()
    > It replaces the following functions from the original code:
        - sun_intensity(), sun_altitude_d(), sun_altitude_r(), HRA(), sun_azimuth_r(), sun_azimuth_d()
    > Verified using https://www.sunearthtools.com/dp/tools/pos_sun.php for Manchester and Cape Town
    > get_position() returns:
        - Sun altitude in radians 
        - Sun azimuth in radians, measured relative to south (east = negative, west = positive)
    > Utilises the following parameter:
        - Overhead sun intensity (OSI) = 476 W/m^2
    '''
    # Return the solar altitude using Pysolar (2023-06-27)
    solar_altitude = get_altitude(latitude, longitude, date)

    # Intensity of the sun based on angle (W/m^2) 
    intensity = OSI * cos(radians(90 - solar_altitude)) # 90 - solar altitude = solar zenith

    # Return outputs
    return intensity, solar_altitude

def model_concrete_temperature(air_temperature, relative_humidity, cloud_height, cloud_cover, wind_speed,
                               concrete_temp, interval_datetime, longitude, latitude,
                               shading_proportion, canal_area):
    '''
    > A nested function to simplify the code structure
    > Takes an input concrete temperature (K) and models a new concrete temperature after 15 minutes 
            - Air temperature (K)
            - Humidity (%)
            - Cloud height and cover (fractional)
            - Wind speed (m/s)
            - Shading proportion [0-1]
    > The following parameters are used:
        - Concrete absorbivity (C_ABSORBIVITY) = 0.85 
    > Nested functions are described in full below
    > Similiar in scope to ConcreteNewTemp() from the original code, but individual functions are outlined more fully here
    > All sub-functions checked on 2023-06-28. Results are identical to the spreadsheet approach (when the same inputs are utilised)
    '''

    # Model the incoming solar radiation (W/m)
    csr = concrete_solar_radiation(interval_datetime, longitude, latitude, shading_proportion, canal_area)
    
    # Modify starting concrete temperature (K) based on radiation input
    concrete_temp = concrete_after_sunlight(csr, concrete_temp, canal_area)

    # Calculate concrete emissions (W/m)
    emissions = concrete_emissions(concrete_temp, canal_area)
    
    # Thermal absorption (W/m2)
    thermal_abs = thermal_absorption(air_temperature, relative_humidity, cloud_height, cloud_cover, canal_area) 

    # Concrete absorption (W/m2)
    concrete_absorbed = thermal_abs * C_ABSORBIVITY
    
    # Determine net radiation
    net_radiation =  emissions - concrete_absorbed

    # Calculate energy transfer to air
    energy_to_air = concrete_convection_air(concrete_temp,
                                            air_temperature, 
                                            wind_speed)

    # Determine heat Transfer to layer above (W/m) 
    heat_transfer = concrete_heat_transfer(concrete_temp, net_radiation, energy_to_air, canal_area)

    # Total energy change (J/m)
    total_energy_change = concrete_energy_change(heat_transfer)

    # Corresponding temperature change (K)
    modelled_temperature_change = concrete_temperature_change(total_energy_change, canal_area)

    # Final concrete temperature for model run, overwriting earlier list
    concrete_temp = [x + y for x, y in zip(concrete_temp, modelled_temperature_change)]

    # Created dictionary to return other important values, used for modelling air temperature change
    concrete_energy_values = {'net_radiation' : net_radiation,
                   'energy_to_air' : energy_to_air, 
                   'solar_radiation' : csr}

    # Return output
    return concrete_temp, concrete_energy_values

def concrete_solar_radiation(date, longitude, latitude, shading_proportion, canal_area):
    '''
    > Calculates solar radiation for concrete (W/m)
    > Similar to ConcreteSolarRadiation() from original code
    - Canal width and length are excluded
    > Utilises the following parameters:
        - Concrete absorbivity (C_ABSORBIVITY) = 0.85 
        - Relative light intensity from the shade (SHADE) = 0.2
    '''

    # Determine shaded and unshaded area (m^2)
    shaded_area = canal_area * shading_proportion
    unshaded_area = canal_area * (1 - shading_proportion)

    # Calculate solar intensity
    intensity = solar_intensity(date, longitude, latitude)[0]

    # If positive
    if intensity > 0: 

        # Calculate solar radiation (AN4)
        solar_rad = (unshaded_area + shaded_area * SHADE) * intensity * C_ABSORBIVITY

    # Set to null, if intensity is negative:
    else:
        solar_rad = 0

    # Return output
    return solar_rad

def concrete_after_sunlight(c_solar_radiation, concrete_temp, canal_area):
    '''
    > Calculates new concrete temperature based on sunlight (°C)
    > Similiar to ConcreteAfterSunlight() from original code
        - However, produces values as /h, rather than /15m (* 3600 / 4)
        - Canal width and length are also excluded
    > Takes in the initial concrete temperatures (K)
    > Utilises the following parameters:
        - Concrete depth (C_DEPTH) = 0.05 m
        - Concrete density (C_DENSITY) = 2238 kg/m^3
        - Concrete specific heat capacity (C_CAPACITY) = 921 J/kg K
    '''  

    # Update the surface concrete temperature (index = 0)
    concrete_temp[0] = ((c_solar_radiation * (3600/4)) / (C_CAPACITY * (canal_area * C_DEPTH * C_DENSITY))) + concrete_temp[0] # Surface temperature   

    # Return list
    return concrete_temp

def concrete_emissions(concrete_temp, canal_area):
    '''
    > Calculates concrete thermal radiation (W/m)
    > Similiar to ConcreteEmissions() from original code
    - Canal width and length are included (area)
    > Utilises the following parameters:
        - Concrete emissivity (C_EMISSIVITY) = 0.85
        - Stefan-Boltzmann Constant (SB) = 5.67*10**(-8)
    ''' 
    # Modification of surface concrete temperature (K)     
    return(canal_area * C_EMISSIVITY * SB * concrete_temp[0]**4)

def thermal_absorption(temperature, relative_humidity, cloud_height, cloud_cover, canal_area):
    '''
    > Calculates thermal absorption from sky
    > Similiar to ThermalAbsorption() from original code
    - Canal width and length are excluded
    > Utilises:
        - Air temperature (K)
        - Relative humidity (%)
        - Cloud cover (fractional)
        - Cloud height (fractional)
    '''   

	# Return output (W/m^2) 
    return ((1 + cloud_height * cloud_cover**2) * 8.78 * (10**-13) * (temperature**5.852) * (relative_humidity**0.07195)) * canal_area

def concrete_convection_air(concrete_temp, air_temperature, wind_speed):
    '''
    > Calculates convection to air (W/m)
    > Similiar to ConcreteEnergytoAir() from original code
    - Canal width and length are excluded
    > Wind speed (m/s)

    '''
    # Maths!
    h = 0.664 * 0.025 * (0.713**0.3) * ((15 * (10**(-6)))**(-0.5)) * (wind_speed**0.5)

    # More maths (using surface and concrete temperature)
    energy_to_air = h * (concrete_temp[0] - air_temperature)

    # If the energy is negative, set to 0
    if energy_to_air < 0:
        energy_to_air = 0

    # Return output
    return energy_to_air

def concrete_heat_transfer(concrete_temp, net_radiation, energy_to_air, canal_area):  
    '''
    > Calculates heat transfer for concrete
    > Similiar to ConcreteHeatTransfer() from original code
        - Canal width and length are excluded
    > Utilises the following parameters:
        - Concrete thermal conductivity (C_CONDUCTIVITY) = 1.21 W/mK
        - Ground thermal conductivity (G_CONDUCTIVITY) = 1 W/mK
        - Concrete depth (C_DEPTH) = 0.05 m
    '''
    # List for storing output, defined length
    transfer = [None] * len(concrete_temp)

    # Use enumerate() to iterate through the concrete temperatures
    for index, value in enumerate(concrete_temp):
    
        # For the surface concrete
        if index == 0:

            transfer[index] = net_radiation + energy_to_air

        # For sections [0.05-0.10] to [0.20-0.25] (m), use concrete conductivity 
        elif index < 5: 

            # Equal to C_CONDUCTIVITY * area * difference in temperature (current section - higher section) / section depth
            transfer[index] = C_CONDUCTIVITY * canal_area * (value - concrete_temp[index - 1]) / C_DEPTH

        # For sections [0.25-0.30] (m) and below, use ground conductivity 
        else: 

            # As above, but using G_CONDUCTIVITY 
            transfer[index] = G_CONDUCTIVITY * canal_area * (value - concrete_temp[index - 1]) / C_DEPTH

    # Return list
    return transfer

def concrete_energy_change(heat_transfer):
    '''
    > Calculates change in concrete energy
    > Similiar to ConcreteEnergyChange() from original code
        - However, produces values as /h, rather than /15m (* 3600 / 4)
    > Utilises the following parameters:
        - Concrete depth (C_DEPTH) = 0.05 m
        - Concrete density (C_DENSITY) = 2238 kg/m^3
        - Concrete specific heat capacity (C_CAPACITY) = 921 J/kg K
        - Ground specific heat capacity (G_CAPACITY) = 1900 J/kg K
        - Ground density (G_DENSITY) = 1500 kg/m3
    '''

    # List for storing output, defined length
    energy_change = [0] * len(heat_transfer)

    # Use enumerate() to iterate through the concrete temperatures
    for index, value in enumerate(heat_transfer):
    
        # While there a deeper sections remaining:
        try:
            
            # (Subsequent heat value - current heat value) * 3600
            energy_change[index] = (heat_transfer[index + 1] - value) * (3600/4)

        # There are no deeper sections remaining:
        except IndexError:

            break
    
    # return output
    return energy_change

def concrete_temperature_change(total_energy_change, canal_area):
    '''
    > Calculates temperature change of concrete
    > Similiar to ConcreteTempRise() from original code
        - Canal width and length are now included
    > Utilises the following parameters:
        - Concrete depth (C_DEPTH) = 0.05 m
        - Concrete density (C_DENSITY) = 2238 kg/m^3
        - Concrete specific heat capacity (C_CAPACITY) = 921 J/kg K
        - Ground specific heat capacity (G_CAPACITY) = 1900 J/kg K
        - Ground density (G_DENSITY) = 1500 kg/m3
    '''
    
    # List for storing output, defined length
    temperature_change = [None] * len(total_energy_change)

    # Use enumerate() to iterate through the concrete temperatures
    for index, value in enumerate(total_energy_change):

        # For shallow sections, use concrete characteristics
        if index < 5: 

            temperature_change[index] = value / (C_CAPACITY * (canal_area * C_DEPTH * C_DENSITY))
        
        # For deeper sections, use ground characteristics
        else:

            temperature_change[index] = value / (G_CAPACITY * (canal_area * C_DEPTH * G_DENSITY))

    # Return output
    return temperature_change


#---- Water modelling ----#

def model_water_temperature(date, longitude, latitude, shading, water_energy, water_temp, split_depths, canal_area):
    '''
    > A nested function to simplify the code structure
    > Takes an input water temperature (K) and energy (J) and produces a new water temperature after 15 minutes based upon...
            - Air temperature (K)
            - Humidity (%)
            - Cloud height and cover (fractional)
            - Wind speed (m/s)
            - Shading proportion [0-1]
    > Nested functions are described in full below
    > Similiar in scope to NewTemp() from the original code, but individual functions are outlined more fully here
    > All sub-functions checked on 2023-06-27. Results are identical to the spreadsheet approach (when the same inputs are utilised)
    '''

    # Return the total absorbed light, absorptivity and solar intensity
    total_absorbed_light, absorptivity, intensity = absorbed_light(date, longitude, latitude)

    # Absorbed/reflected radiation due to sunlight 
    absorbed_radiation, reflected_radiation = absorbed_reflected_radiation(shading, total_absorbed_light, intensity, absorptivity, canal_area)

    # Energy gain
    total_sunlight_energy = sunlight_energy(absorbed_radiation, split_depths)

    # Net energy change
    net_energy_change = net_sunlight_energy(total_sunlight_energy, water_energy)

    # Change in water temperature
    water_temperature_change = temperature_change(net_energy_change, canal_area)

    # Final water temperature for model run, overwriting earlier list
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
    > Tested against the Excel spreadsheet: returns an identical value
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

def absorbed_light(date, longitude, latitude):
    '''
    > This function returns the total absorbed light (W/m^2), and the solar absorptivity
    > Used later in absorbed_reflected_radiation()
    > Same functionality as the AbsorbedLight() function from the original code, but utilising outputs from earlier functions
    '''

    # Calculate the overhead sun intensity and the solar altitude (°)
    intensity, solar_altitude = solar_intensity(date, longitude, latitude)

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
    
    # Return the absorbed radiation and reflected radiation
    return absorbed_radiation, reflected_radiation


def sunlight_energy(absorbed_radiation, split_depths):
    '''
    > This function returns the energy gained by sunlight (J)
    > Similar functionality as the SunlightEnergy() function from the original code
    > However, converts to J/h, rather than J/15m (* 3600 / 4)
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

            # Calculates energy: converting from W/m^2 to J/s 
            total_sunlight_energy[index] = (absorbed_radiation * exp(-VEC * (depth - INTERVAL)) - absorbed_radiation * exp(-VEC * depth)) * (3600/4)
        
        # There are no deeper sections remaining:
        except IndexError:
            break

    # Return output
    return total_sunlight_energy

    
def net_sunlight_energy(total_sunlight_energy, water_energy):
    '''
    > This function returns the net energy gained by sunlight (J), incorporating the starting energy amount
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
                        water_temperature, canal_area, initial_water_k):
    '''
    > A nested function to simplify the code structure
    > Takes an input water temperature (K) and energy (J) and produces a new water energy (J) after 15 minutes
    > Nested function are described in full below
    > Similiar in scope to HeatTransfer() from the original code, but individual functions are outlined more fully here
    > All sub-functions checked on 2023-06-27 (some modifications required). Results now identical to the spreadsheet approach (when the same inputs are utilised)
    > Additional fluxes that could be added:
        - Energy loss for bed heat exchange (Qbhf)
        - Energy change, friction at the bed and banks (Qf)
        - Energy change, groundwater discharge, longitudinal advective heat flux (Qa)
    '''
     
    # Model net thermal emissions, based on surface water temperature
    modelled_thermal_emissions = net_thermal_emissions(air_temperature, relative_humidity, cloud_height, cloud_cover, initial_water_k, canal_area) 
    
    # Evaporation rate
    modelled_evaporation_rate, modelled_surface_evaporation_rate, modelled_vapour_pressure = evaporation_rate(air_temperature, wind_speed, relative_humidity, initial_water_k)

    # Modelled evaporation energy
    modelled_evaporation_energy = evaporation_energy(modelled_evaporation_rate)

    # Modelled evaporation power
    modelled_evaporation_power = evaporation_power(modelled_evaporation_energy, canal_area)

    # Cooling, assuming 50% split of evaporation cooling effect between air and water
    modelled_evaporation_cooling = modelled_evaporation_power / 2

    # Sensible heat calculations
    modelled_sensible_heat = sensible_heat(air_temperature, air_pressure, initial_water_k, 
                  modelled_surface_evaporation_rate, modelled_vapour_pressure, modelled_evaporation_power)
    
    # Energy for the surface layer (J/h)
    modelled_surface_heat_transfer = surface_heat_transfer(modelled_evaporation_cooling, modelled_thermal_emissions, modelled_sensible_heat)

    # Energy for all layers (J/h)
    modelled_energy = canal_heat_transfer_convection(modelled_surface_heat_transfer, water_temperature, canal_area)

    # Created dictionary to return other important values, used for modelling air temperature change
    energy_values = {'sensible_heat' : modelled_sensible_heat,
                   'evaporation_cooling' : modelled_evaporation_cooling,
                   'net_thermal_emissions' : modelled_thermal_emissions}

    # Returns output
    return modelled_energy, energy_values

def water_emissions(water_temperature, canal_area):
    '''
    > Function returns the radiative emissions from water (W/m^2), based on the water temperature (K)
    > Similiar functionality as the WaterEmmissions() function from the original code:
        - Rather than normalising to the width of the canal, this operates on a point by point basis
    > Uses the following parameters:
        - Stefan-Boltzmann Constant (SB) = 5.67*10**(-8)
        - Emissivity of water (E) = 0.95
    > Similiar functionality as the WaterEmmissions() function from the original code
    '''
    # Return output (W/m^2) = W = ε σT4 (https://doi.org/10.1016/B978-0-12-083980-3.50013-9)
    return (E * SB * (water_temperature**4)) * canal_area


def net_thermal_emissions(temperature, relative_humidity, cloud_height, cloud_cover, water_temperature, canal_area):
    '''
    > Net thermal emissions (emissions - thermal absorption)
    > Identical to NetThermalEmissions() from original code
    '''
    # Return output
    return water_emissions(water_temperature, canal_area) - thermal_absorption(temperature, relative_humidity, cloud_height, cloud_cover, canal_area)

def surface_evaporation_rate(water_temperature): 
    '''
    > This function returns the evaporation rate (mm/day), utilising the surface water temperature in K
    > Similar functionality as the Ew() function from the original code
    ''' 
    # Saturated vapour pressure at surface water temperature, converting to mbar
    evaporation_rate = exp(77.345 + 0.0057 * (water_temperature) - 7235 / (water_temperature)) / pow((water_temperature), 8.2) / 100   

    # Return output
    return evaporation_rate

def evaporation_rate(air_temperature, wind_speed, relative_humidity, initial_water_k):
    '''
    > This function returns the evaporation rate (m/day), plus the surface rate and vapour pressure
    > Identical functionality as the Ev() function from the original code, but simplified
    '''

    # Evaporation rate (mm/day, converted to m/day)
    vap_pressure = vapour_pressure(air_temperature, relative_humidity)
    surf_evaporation_rate = surface_evaporation_rate(initial_water_k)

    # Return three values 
    return 0.165 * (0.8 + 0.864 * wind_speed) * (surf_evaporation_rate - vap_pressure) / 1000, surf_evaporation_rate, vap_pressure

def vapour_pressure(air_temperature, relative_humidity):   
    '''
    > This function returns air vapour pressure (mbar), utilising the following climate parameters
        - Air temperature (K)
        - Relative humidity
    > The function also utilises the following parameters:
        - Constant to calculate dew point (DEW_POINT_B) = 17.67
        - Constant to calculate dew point (DEW_POINT_C) = 243.5
        - Ratio of latent heat to water vapour gas constant (LRV) = (L*1000)/RV  
        - Baseline Pressure at Triple Point (BASE_PRESSURE) = 0.611 kPa
        - Triple point of water (TPW) = 273.0 K
    > Identical functionality as the Ea() function from the original code, but with clearer signposting of constants
    '''

    # Calculate dew point temperature (converting between to K and °C)
    y = log(relative_humidity * 0.01) + DEW_POINT_B * (air_temperature - 273.15) / (DEW_POINT_C + (air_temperature - 273.15))
    dpt = DEW_POINT_C * y / (DEW_POINT_B - y) + 273.15

    # Vapour pressure at air temperature (K, convert to mbar)
    vapour_pressure = BASE_PRESSURE * (exp(LRV * ((1.0 / TPW) - (1 / dpt)))) * 10

    # Return output
    return vapour_pressure

def evaporation_energy(modelled_evaporation_rate):
    '''
    > This function returns the energy loss to evaporation (J/m^2/day)
    > The function also utilises the following parameters:
        - Latent heat of vapourisation (L) = 2454.9 J/g
        - Specific weight of water (WATER_WEIGHT) = 999286 g/m^3
    > Identical functionality as the EvaporationEnergy() function from the original code, but simplified
    '''

	# Return energy 
    return modelled_evaporation_rate * L * WATER_WEIGHT
    
def evaporation_power(modelled_evaporation_energy, canal_area):
    '''
    > This function returns the power loss to evaporation (W/m)
    > Similiar functionality as the EvaporationPower() function from the original code
    '''

    # Power loss to evaporation per metre area (W/m)
    return modelled_evaporation_energy * canal_area * (1 / 86400)  

def sensible_heat(air_temperature, pressure, water_temperature, 
                  modelled_surface_evaporation_rate, modelled_vapour_pressure, modelled_evaporation_power):
    '''
    > Calculates sensible heat transfer using pressure, and water/air temperatures (K)
    > Similiar functionality to SensibleHeat() from the original code
    '''   
    # Bowen Ratio
    bowen = (0.61 * pressure * (water_temperature - air_temperature) / (modelled_surface_evaporation_rate - modelled_vapour_pressure)) / 1000

    # Return sensible heat (W/m)
    return modelled_evaporation_power * bowen

def surface_heat_transfer(modelled_evaporation_cooling, modelled_thermal_emissions, modelled_sensible_heat):     
    '''
    > Calculates heat transfer between water and air, incorporating evaporation, thermal emissions and sensible heat
    > Similiar functionality to SurfaceHeatTransfer() from the original code
    '''   
    # Return output, converting from /s to /15m
    return (modelled_evaporation_cooling + modelled_thermal_emissions + modelled_sensible_heat) * 3600/4

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

            # Obtaining the preceding temperature
            higher_depth = water_temperature[index - 1]
      
            # Calculate conductivity
            conductivity = HC_WATER * canal_area * (current_temperature - higher_depth) / FLOAT_INTERVAL

            # If the temperature of the current depth is greater than the temperature above:
            if current_temperature > higher_depth: 

                # Energy value (/15m)
                modelled_heat_transfer[index] = (conductivity + CONVECTION * canal_area * (current_temperature - higher_depth)) * 3600/4
                
            # If temperature <
            else:
                
                # Energy value (/15m)
                modelled_heat_transfer[index] = conductivity * 3600/4

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
        - The convective effect dependent on absolute distance from 4 degrees
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

            # Obtaining the preceding temperature
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
                    modelled_heat_transfer[index] = (conductivity + CONVECTION * canal_area * (lower_temp - upper_temp)) * 3600 / 4

                # If temperature in lower layer <= temperature in the layer above, -CONVECTION
                else:
                    # Energy value (/15m)           
                    modelled_heat_transfer[index] = (conductivity + -CONVECTION * canal_area * (lower_temp - upper_temp)) * 3600 / 4

            # If difference from 4°C in lower layer <= difference in the layer above
            else:

                # Energy value (/15m)
                modelled_heat_transfer[index] = conductivity * 3600 / 4

    # Return energy values
    return modelled_heat_transfer

#---- Air temperture modelling ----#

def model_urban_cooling(modelled_sensible_heat, modelled_evaporation_cooling, net_thermal_emissions,
                        absorbed_sunlight, 
                        concrete_net_radiation, concrete_energy_to_air, concrete_solar_radiation,
                        buffered_canal_area):
    '''
    > A nested function to simplify the code structure
    > Models air temperature change in response to changing canal temperatures i.e. urban cooling
    > Nested functions are described in full below
    > Similiar in scope to NetTempChange() from the original code, but rather than calculating many of the variables again...
    > ... these are utilised from earlier functions
    '''



    # Energy per kelvin change (kJ/K) = Thermal capacity * mass of air
    energy_per_kelvin = A_CAPACITY * (A_DENSITY * (buffered_canal_area * A_HEIGHT))

    ''' > Calculations for water '''

    # Energy gain of air from the canal (W/m)
    water_energy_gain = -((-modelled_evaporation_cooling) - modelled_sensible_heat - net_thermal_emissions)

    # Water net energy gain (energy gain - enery removed as absorped sunlight) (W/m)
    water_net_energy_gain = water_energy_gain - absorbed_sunlight

    # Temperature change of air in surrounding area (K) based on water temperature
    air_temperature_change_water = water_net_energy_gain / energy_per_kelvin

    ''' > Calculations for concrete ''' 

    # Energy gain of air from the concrete (W/m)
    concrete_energy_gain = concrete_net_radiation + concrete_energy_to_air
    
    # Concrete net energy gain (energy gain - enery removed as absorped sunlight) (W/m)
    concrete_net_energy_gain = concrete_energy_gain - concrete_solar_radiation

    # Temperature change of air in surrounding area (K) based on water temperature
    air_temperature_change_concrete = concrete_net_energy_gain / energy_per_kelvin

    ''' > Water / concrete differences ''' 
    
    # Energy and temperature difference
    difference_energy = water_net_energy_gain - concrete_net_energy_gain # W/m
    difference_temperature = air_temperature_change_water - air_temperature_change_concrete # K

    # Energy delivered to air over 15 minutes (kJ/(m*h)) 
    energy_delivered = difference_energy * (3600/4) / (buffered_canal_area * A_HEIGHT)

    # Corresponding temperature change
    temperature_change = energy_delivered / energy_per_kelvin

    # Created dictionary to store key output values
    output_values = {'air_water_k' : air_temperature_change_water,
                   'air_concrete_k' : air_temperature_change_concrete,
                   'temperature_difference' : difference_temperature,
                   'energy_difference' : difference_energy,
                   'modelled_temperature_change' : temperature_change,
                   'evaporation_cooling' : modelled_evaporation_cooling,
                   'sensible_heat' : modelled_sensible_heat,
                   'thermal_emissions' : net_thermal_emissions,
                   'water_net_energy_gain' : water_net_energy_gain, 
                   'concrete_net_energy_gain' : concrete_net_energy_gain}

    # Air temperature suppression
    return output_values


def model_spin_up(water_temp, water_energy, concrete_temp, 
                  split_points, canal_area, lat, lon, 
                  spin_path, repetitions, record, canal_id):
    '''
    > Function for model spin up 
    > Using input temperatures (K) and energy values (J), the model is run until equilibrium is reached
    '''

    # Open the spin up record
    with open(spin_path) as path:
            spin_up_climate = load(path)

    # Init lists for storing outputs
    output_water = []
    output_energy = []
    output_concrete = []

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

            # Model water temperatures
            water_temp, _ = model_water_temperature(local_datetime, lon, lat, 
                                                        0, # Presence/absence of shading
                                                        water_energy, # Energy (J)
                                                        water_temp, # Temperature (K)
                                                        split_points, # Canal depths (m)
                                                        canal_area) # Canal area (m^2)
            
            # Try and extract current climate
            try:
                current_climate = spin_up_climate[start_unix]

            # Reached the end of the record
            except KeyError:

                # Skip to next iteration
                local_datetime = end_datetime

            # Valid datetime
            else:
                
                # Model water energy
                water_energy, _ = model_energy_change(current_climate['air_temperature'], # Air temperature (K)
                                    current_climate['humidity'], # Relative humidity (%)
                                    0.06, # Cloud height (fractional)
                                    current_climate['cloud_cover'] / 100, # Cloud cover (fractional),  
                                    current_climate['wind_speed'], # Wind speed (m/s)
                                    current_climate['air_pressure'], # Air pressure
                                    water_temp, # Water temperature (K)
                                    canal_area, # Canal area (m^2)
                                    initial_water_k) # Starting surface water temperature (K) for this model step 
                

                # Model concrete temperature
                concrete_temp, _ = model_concrete_temperature(current_climate['air_temperature'], # Air temperature (K)
                                                current_climate['humidity'],  # Relative humidity (%)
                                                0.06, # Cloud height (fractional)
                                                current_climate['cloud_cover'] / 100, # Cloud cover (fractional)
                                                current_climate['wind_speed'], # Wind speed (m/s)
                                                concrete_temp, # Concrete temperature (K)
                                                local_datetime, # Current datetime (15m interval)
                                                lon, lat, # Location
                                                0, # [0-1]
                                                canal_area) # Area (m^2)

                # Append surface temperatures and energy to output
                output_water.append(water_temp[0])
                output_energy.append(water_energy[0])
                output_concrete.append(concrete_temp[0])
                
                # Add 15 minutes to the datetime
                local_datetime += timedelta(hours = 0.25)

    # Write to output file
    if record:

        # Zip lists and export
        output = array(list(zip(output_water, output_energy, output_concrete)))
        savetxt(f"../outputs/spin_up_output_{canal_id}.csv", output, delimiter = ",", header="water_k, water_j, concrete_k", fmt ='%f')

    # Return modified temperatures (K) and energy (J), and time-series values
    return water_temp, water_energy, concrete_temp, output_water, output_energy, output_concrete


