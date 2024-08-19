'''
Algorithm: Estimating the cooling potential of urban canals 
Github: https://github.com/matt-tomkins/national-canal-cooling
Conference (GISRUK2024): https://zenodo.org/doi/10.5281/zenodo.10927598 
Paper: **Link goes here**

Code contributors: 
- Dr. Matt Tomkins: https://github.com/matt-tomkins
- Dr. Jonny Huck: https://github.com/jonnyhuck

Acknowledgements: 
- Algorithm based upon earlier work by: 
    - Dr. Harry Mcdonald, supervised by Dr. Jo Tippett
    - Sophie Taylor, supervised by Dr. Jonny Huck
- Project funded by the Canal and River Trust

Key packages: 
>---- pybdshadow ----> https://pypi.org/project/pybdshadow/
>---- Pysolar -------> https://pysolar.readthedocs.io/en/latest/
>---- geoalchemy ----> https://geoalchemy-2.readthedocs.io/en/latest/

Datasets:
>---- Building geometries ------> OS MasterMap (https://beta.ordnancesurvey.co.uk/products/os-mastermap-topography-layer)
>---- Canal characteristics ----> Canal and River Trust (https://canalrivertrust.org.uk/)
>---- Climate data -------------> OpenWeather (https://openweathermap.org/)
'''

# Force use of Shapely 2.0
from os import environ
environ['USE_PYGEOS'] = '0'

# Required packages
from json import dump
from os import cpu_count
from numpy import linspace
from time import perf_counter
from pandas import to_datetime
from pyproj import Transformer
from datetime import timedelta
from pyogrio import read_dataframe

# Import user-defined parameters and functions 
from params import * 
from functions import * 
from auxiliary import *

# Prevent datetime warning: https://github.com/pingswept/pysolar/issues/139
warnings.filterwarnings('ignore', message="I don't know about leap seconds after 2021")

#---- Main code ----#
def parallel_model(input_file_path, start_date, days, n_cores, shading):
    ''' > Function to to run the model using multi-processing '''

    # Read the dataframe directly using pyogrio
    canals_gdf = read_dataframe(input_file_path)

    # Split into n clusters, each will be assigned a core
    canals_gdf = array_split(canals_gdf, n_cores)

    # Initialises the Pool for multiprocessing (n cores)
    p = Pool(processes = n_cores)

    # Define arguments, some fixed, some variable
    # Canal data, start date, days of analysis, shading boolean
    args = [(i, start_date, days, shading) for i in canals_gdf]

    # Generate results
    p.starmap(main, args)

def main(canals_gdf, start_date, duration_days, model_shading):

    # Sets the duration of analysis (in hours)
    duration_hours = duration_days * 24

    # Include or exclude shading effects, updates path to file
    if model_shading:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Try and load list of "highly" shaded locations (> threshold)
    try:
        with open(f"../outputs/high-shading-features.json") as shading_subset:
            high_shading = load(shading_subset)

            # Extract feature IDs from dict keys
            high_shading_features = high_shading.keys()

    # If the file is not found, print warning
    except FileNotFoundError:
        print("WARNING. No record of high-shading features. Additional metrics will not be saved. Check file path.")

        # Init empty dict
        high_shading_features = {}

    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # Skip canals with an area < 500 m^2 or > 200_000 m^2
        if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
            continue

        # If the file already exists, skip
        if isfile(f"../results/{folder_path}-shading/model-output-{feature_tuple.code_id}.json"):
            print(f"../results/{folder_path}-shading/model-output-{feature_tuple.code_id}.json exists")
            continue

        # Create dictionary for storing results
        else:
            output = {}

        # Include shading effect
        if model_shading:

            # Try and load shading record for current location
            try:
                with open(f"../../national-canal-cooling-data/shading/shading-{feature_tuple.code_id}.json") as shading_path:
                    shading_record = load(shading_path)

            # If the file is not found
            except FileNotFoundError:

                # Skip to next feature (the model cannot be run without shading data if shading == True)
                print(f"Shading file MISSING for feature {feature_tuple.code_id}")
                continue

        # Try
        try:

            # Load climate record for current location
            with open(f"../../national-canal-cooling-data/interpolated-climate/data-{feature_tuple.code_id}.json") as climate_path:
                climate_record = load(climate_path)

        # If the file is not found
        except FileNotFoundError:

            # Skip to next feature (the model cannot be run without climate data)
            continue

        # Localise niave datetime to aware datetime (UK) and then convert to UTC
        local_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
        
        # Datetime checking
        # print(f"The input datetime is {start_date}. In UTC, this is {local_datetime}")

        # Sets end date
        end_datetime = local_datetime + timedelta(hours = duration_hours)
        
        # Transformer for projected to geographic coordinates (BNG | ING > WGS84)
        transformer = Transformer.from_crs(canals_gdf.crs, 4326)

        # Feature centroid and area (shapely)
        latitude, longitude = transformer.transform(feature_tuple.geometry.centroid.x, feature_tuple.geometry.centroid.y)
        feature_area = feature_tuple.geometry.area

        # Calculate radius and draw circle
        radius = sqrt(feature_area / pi)
        circle = Point(0, 0).buffer(radius, resolution = 10_000)

        # For modelling air temperature change, we can either  use:
        buffered_area = circle.buffer(50).area # Use a circle with the same area as the canal, and buffer by 50 m (minimising shape effects)
        # buffered_area = feature_tuple.geometry.buffer(50).area # Buffer the geometry by 50 m (but this is influenced by feature complexity)

        # Canal depth (m)
        canal_depth = feature_tuple.draught

        # Round the canal draught 'up' to the nearest interval (FLOAT_INTERVAL, 20 cm) e.g. canal depth of 1.7 m becomes 1.8 m
        adjusted_depth = round_up(canal_depth, FLOAT_INTERVAL)

        # Determine the number of layers
        layer_count = int(adjusted_depth / FLOAT_INTERVAL)
        
        # Determine split points (end depths) based on the chosen interval
        depths = linspace(start = 0 + FLOAT_INTERVAL, stop = adjusted_depth, num = layer_count)

        # Convert to Decimal format
        split_points = [Decimal(x).quantize(Decimal("1.0")) for x in depths]

        # Init lists of starting values for the water (K = 282, J = 0)
        water_temp = [MEAN_KELVIN] * len(split_points)
        water_energy = [0] * len(split_points)

        # Init list of starting value for concrete (K = 282), comprising 10 depth layers [0 - .5m]
        concrete_temp = [MEAN_KELVIN] * 10

        # Model spin up (30 iterations), based on a composite climate record for 2021-12-15
        water_temp, water_energy, concrete_temp, _, _, _, = model_spin_up(water_temp, water_energy, concrete_temp, split_points, feature_area,
                                                                latitude, longitude, 
                                                                "../outputs/spin_up_climate.json", 30, False, feature_tuple.code_id)
        
        # Run the full model: iterate until we reach the chosen duration
        while local_datetime < end_datetime:

            # Convert to unix, as string
            unix_time = str(int(local_datetime.timestamp()))

            # If this datetime has already been modelled
            if unix_time in output:

                # Finish the iteration (+ 15 minutes)
                local_datetime += timedelta(hours = 0.25)

            # New datetime encountered
            else:
                
                # Extract current climate
                try:
                    climate = climate_record[unix_time]

                # The datetime does not exist as a key
                except KeyError:
                    
                    # Reached the end datetime
                    break

                # Include the effects of building shading
                if model_shading:
                    
                    # Extract current shading record
                    shading = shading_record[unix_time]

                    # Extracts modelled shading value [0 - 1]
                    shading_proportion = shading['shading']

                # Exclude the effects of building shading
                else:
                    
                    # Set building shading to 0 (no shading)
                    shading_proportion = 0

                # Starting surface water temperature (K) for this model step
                initial_water_k = water_temp[0]

                # Model water temperatures
                water_temp, absorbed_radiation = model_water_temperature(local_datetime, longitude, latitude, 
                                                    shading_proportion, # [0-1]
                                                    water_energy, # Energy (J)
                                                    water_temp, # Temperature (K)
                                                    split_points,
                                                    feature_area) # Area (m^2)
                
                # Model energy change
                water_energy, energy_values = model_energy_change(climate['air_temperature'], # Air temperature (K)
                                    climate['humidity'], # Relative humidity (%)
                                    0.06, # Cloud height (fractional)
                                    climate['cloud_cover'] / 100, # Cloud cover (fractional),  
                                    climate['wind_speed'], # Wind speed (m/s)
                                    climate['air_pressure'], # Air pressure
                                    water_temp, # Modified water temperature (K), after model_water_temperature()
                                    feature_area, # Canal area (m^2)
                                    initial_water_k) # Starting surface water temperature (K) for this model step
                
                # Model concrete temperature
                concrete_temp, concrete_energy_values = model_concrete_temperature(climate['air_temperature'], # Air temperature (K)
                                                climate['humidity'],  # Relative humidity (%)
                                                0.06, # Cloud height (fractional)
                                                climate['cloud_cover'] / 100, # Cloud cover (fractional)
                                                climate['wind_speed'], # Wind speed (m/s)
                                                concrete_temp, # Concrete temperature (K)
                                                local_datetime, # Current datetime (15m interval)
                                                longitude, latitude, # Location
                                                shading_proportion, # [0-1]
                                                feature_area) # Area (m^2)
                
                
                # Model air temperature change
                modelled_air_temp = model_urban_cooling(energy_values['sensible_heat'], 
                                                        energy_values['evaporation_cooling'],
                                                        energy_values['net_thermal_emissions'],
                                                        absorbed_radiation, 
                                                        concrete_energy_values['net_radiation'],
                                                        concrete_energy_values['energy_to_air'],
                                                        concrete_energy_values['solar_radiation'],
                                                        buffered_area) # Canal area, buffered by 50 m
                
                # Output of model run
                # print(f"At {local_datetime}, the new water temperature is {water_temp[0]} K and the concrete temperature is {concrete_temp[0]} K")

                # If CRT water temperature data exists
                if feature_tuple.code_id in ['lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11']:

                    # Add values to output dictionary (surface wT, surface cT, ΔJ, ΔT), plus values at depth (0 - 100 cm), all rounded to 3 dp
                    output.update({unix_time : {'surface_concrete_k' : round(concrete_temp[0], 3), # Concrete temperature at surface
                                            'surface_water_k' : round(water_temp[0], 3), # Temperature at surface (0 - 20 cm)
                                            'depth_water_20' : round(water_temp[1], 3), # Temperature at 20 - 40 cm
                                            'depth_water_40' : round(water_temp[2], 3), # Temperature at 40 - 60 cm
                                            'depth_water_60' : round(water_temp[3], 3), # Temperature at 60 - 80 cm
                                            'depth_water_80' : round(water_temp[4], 3), # Temperature at 80 - 100 cm
                                            'energy_difference' : round(modelled_air_temp['energy_difference'], 3), # Energy difference (J)
                                            'temperature_difference' : round(modelled_air_temp['temperature_difference'], 3)}}) # 'Air' temperature difference (K)
                
                # If feature ID is classed as "highly" shaded (> threshold), extract a wider range of values for testing
                elif feature_tuple.code_id in high_shading_features:

                    # Add values to output dictionary (surface wT, surface cT, ΔJ, ΔT), all rounded to 3 dp
                    output.update({unix_time : {'surface_concrete_k' : round(concrete_temp[0], 3), # Concrete temperature at surface
                                            'surface_water_k' : round(water_temp[0], 3), # Water temperature at surface (0 - 20 cm)
                                            'energy_difference' : round(modelled_air_temp['energy_difference'], 3), # Energy difference (J)
                                            'temperature_difference' : round(modelled_air_temp['temperature_difference'], 3), # 'Air' temperature difference, water-concrete combined (K) 
                                            'water ΔT' : round(modelled_air_temp['air_water_k'], 3), # 'Air' temperature difference from water alone (K) 
                                            'concrete ΔT' : round(modelled_air_temp['air_concrete_k'], 3), # 'Air' temperature difference from concrete alone (K) 
                                            'evaporation_cooling' : modelled_air_temp['evaporation_cooling'], # Flux of energy associated with evaporation
                                            'sensible_heat' : modelled_air_temp['sensible_heat'], # Flux of energy assocaited with sensible heat
                                            'thermal_emissions' : modelled_air_temp['thermal_emissions'], # Flux of energy associated with thermal emissions
                                            'absorbed_radiation' : absorbed_radiation, # Absorbed radiation (W/m^2)
                                            'concrete_energy_to_air' : concrete_energy_values['energy_to_air'],
                                            'concrete_net_radiation' : concrete_energy_values['net_radiation'],
                                            'concrete_solar_radiation' : concrete_energy_values['solar_radiation'],
                                            'water_net_energy_gain' : modelled_air_temp['water_net_energy_gain'], # Total energy exchange from water
                                            'concrete_net_energy_gain' : modelled_air_temp['concrete_net_energy_gain']}}) # Total energy exchange from concrete
                
                # No measured data at this location
                else:


                    # Add values to output dictionary (surface wT, surface cT, ΔJ, ΔT), all rounded to 3 dp
                    output.update({unix_time : {'surface_concrete_k' : round(concrete_temp[0], 3), # Concrete temperature at surface
                                            'surface_water_k' : round(water_temp[0], 3), # Temperature at surface (0 - 20 cm)
                                            'energy_difference' : round(modelled_air_temp['energy_difference'], 3),  # Energy difference (J)
                                            'temperature_difference' : round(modelled_air_temp['temperature_difference'], 3)}}) # 'Air' temperature difference (K)

                # Finish the iteration (+ 15 minutes)
                local_datetime += timedelta(hours = 0.25)
        
        # Save the output to the relevant directory (including- | excluding-shading)
        dump(output, open(f"../results/{folder_path}-shading/model-output-{feature_tuple.code_id}.json", 'w'))

        # Summary statement
        print(f"Completed {feature_tuple.code_id}")

def test(start_date):
    '''
    > Function to compare the model performance with the spreadsheet approach 
    > Comparison is for the first datetime only "2018-08-26 00:00:00" (0.00 tab)
    > Minor differences between the approaches include:
        - Canal depth interval:
            - Code = 0.2 m
            - Spreadsheet = 0.2111r 
        - Location:
            - Code: lat-lon for London (approx.)
            - Spreadsheet = unknown
        - Solar altitude:
            - Code: calculated using Pysolar: https://pysolar.readthedocs.io/en/latest/
            - Spreadsheet: manually
    > Taking these differences into account, the approaches produce identical results
    '''

    # Expected results
    expected_water = [291.445367259280000, 291.739384566412000, 291.874909050089000, 291.740625828866000, 291.590350360768000, 
                        291.520901023302000, 291.489490812955000, 291.474684741435000, 291.467550701115000]
    expected_w_energy = [1345366.4666458, 638436.329968, 294281.1592, -3951.191906, -4421.752829, -2043.499237, 
                        -924.2239483, -435.658524, -209.9142552]
    expected_concrete = [289.5454631, 291.4269587, 292.5376914, 293.0711875, 293.1908854, 
                         293.0334311, 292.7391121, 292.4809406, 292.3390251, 292.3207526]
    
    # Expected air temperature values 
    # Note: modelled_temperature_change = * 1.5 as the Excel approach uses 10-min intervals, rather than 15-m
    expected_air_values = {'air_water_k': 1.7544383240299, 
                           'air_concrete_k': 0.9803777463983330, 
                           'temperature_difference': 0.774060577631567, 
                           'energy_difference': 659.530574565201, 
                           'modelled_temperature_change': 0.464436346578941 * 1.5}  

    # Localise niave datetime to aware datetime (UK) and then convert to UTC
    local_datetime = to_datetime(start_date)\
        .tz_localize('Europe/London')\
        .tz_convert('UTC')
    
    # Initial water temperatures (K) and energy values
    water = [18.33284839, 18.62236339, 18.75216265, 18.59059478, 18.44012332, 18.37079476, 18.33944445, 18.32466334, 18.31753074]
    energy = [1022817.5524852200000, 628659.9028955620000, 281849.2836874990000, -4754.0240826463200, -4427.5197205518100, -2039.9455136730200, -922.4613512158030, -434.9240887584510, -209.8717799875320]
    water = [x + 273.15 for x in water]

    # Initial concrete temperatures (K) and energy values
    concrete_k = [289.7778903933270, 291.5786986757990, 292.6614947829630, 293.1584705321070, 293.2424354091520, 
                  293.0488807801040, 292.7330163436110, 292.4654103379710, 292.3207525655000, 292.3207525655000]

    # Location (approx.)
    lat, lon = 51.51376017851662, -0.09850663133059452

    # Set the canal depth (m)
    canal_depth = 1.8

    # Round the depth 'up' to the nearest 20 cm e.g. canal depth of 1.7 m becomes 1.8 m
    adjusted_depth = round_up(canal_depth, FLOAT_INTERVAL)

    # Number of layers
    layer_count = int(adjusted_depth / FLOAT_INTERVAL)
    
    # Determine split points (start depths) based on the calculated depth
    depths = linspace(start = 0 + FLOAT_INTERVAL, stop = adjusted_depth, num = layer_count)

    # Convert to Decimal format
    split_points = [Decimal(x).quantize(Decimal("1.0")) for x in depths]

    # Climate inputs
    air_temperature = 290.15 # (K)
    relative_humidity = 82 # (%)
    cloud_height = 0.23 # (Fractional)
    cloud_cover = 0.3 # (Fractional)
    wind_speed = 3.0 # (m/s)
    air_pressure = 1015

    # Store starting water temp
    starting_temp = water[0]
    canal_area = 11.9 # (m^2)

    # Model water temperature 
    water, absorbed_radiation = model_water_temperature(local_datetime, lon, lat, 
                                            1, # Presence/absence of shading
                                            energy, # Energy (J)
                                            water, # Temperature (K)
                                            split_points,
                                            canal_area) # Area (m^2)
    
    # Model energy change
    energy, energy_values = model_energy_change(air_temperature, relative_humidity, cloud_height, cloud_cover, wind_speed, air_pressure,
                    water, canal_area, starting_temp)
    
    # Model concrete temperature
    concrete_k, concrete_energy_values = model_concrete_temperature(air_temperature, relative_humidity, cloud_height, cloud_cover, wind_speed,
                                                                    concrete_k, local_datetime, lon, lat,
                                                                    1, canal_area)
    
    # Model air temperature change
    modelled_air_temp = model_urban_cooling(energy_values['sensible_heat'], 
                                            energy_values['evaporation_cooling'], 
                                            energy_values['net_thermal_emissions'],
                                            absorbed_radiation, 
                                            concrete_energy_values['net_radiation'], 
                                            concrete_energy_values['energy_to_air'], 
                                            concrete_energy_values['solar_radiation'],
                                            100, # Area of buffer
                                            canal_area) # Area of canal
    
    # Differences between the modelled value (m) and the expected value (e), at 6 dp
    diff_water = [f"{m - e:.6f}" for m,e in zip(water, expected_water)]
    diff_w_energy = [f"{m - e:.6f}" for m,e in zip(energy, expected_w_energy)]
    diff_concrete = [f"{m - e:.6f}" for m,e in zip(concrete_k, expected_concrete)]
    diff_air = [{k : {'expected' : expected_air_values[k], 
                      'modelled' : v}} for k, v in modelled_air_temp.items()]
      
    # Print outputs
    print(f"Difference in water K: {diff_water}")
    print(f"Difference in water J: {diff_w_energy}")
    print(f"Difference in concrete K: {diff_concrete}")
    print()
    print("Air temperature outputs:")
    for item in diff_air:
        print(item)
    
# The script is meant to be run
if __name__ == '__main__':
    
    # Initialise timer
    start = perf_counter()
    
    # Files for Britain and Ireland
    file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", 
                  "../data/urban-canals-ireland-dissolved-id-tm65.shp"]
    
    # Iterate
    for f in file_names:

        # Run the parallel model (file path, start date, duration, n pools, shading [True | False])
        parallel_model(f, "2022-01-01 00:00:00", 365, 5, True)
    
    
    '''
    # > Code to run model for a single feature or non-parallel
    
    # Read the dataframe directly using pyogrio
    canals_gdf = read_dataframe("../data/canal-geometries-filtered-modified-dissolved-id.shp")
    
    # Select features of interest
    #options = ['lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11'] # Measurement locations
    canals_gdf = canals_gdf.loc[canals_gdf['code_id'].isin(options)]

    # Run main function, start date, duration, shading (True | False)
    main(canals_gdf, "2022-01-01 00:00:00", 365, True)
    
    '''

    # Run the test function, utilising inputs from the Excel spreadsheet
    # test("2018-08-26 00:00:00")

    # Completion
    print(f"Code completion in {perf_counter() - start:.2f} seconds")


