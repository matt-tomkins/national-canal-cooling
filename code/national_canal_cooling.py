'''
Algorithm: Estimating the effects of canals on urban air temperatures
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

# Error with pyroj.Transformer
# Source: https://gis.stackexchange.com/questions/373550/first-call-to-transform-fails-with-inf-all-subsequent-calls-are-ok-what-cou
environ['PROJ_NETWORK'] = 'OFF'

# Required packages
from json import dump
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

# Additional PySolar warnings when using radiation.get_radiation_direct in solar_intensity()
warnings.filterwarnings('ignore', message="invalid value encountered in scalar multiply")
warnings.filterwarnings('ignore', message="overflow encountered in exp")
warnings.filterwarnings('ignore', message="overflow encountered in scalar multiply")

# Function to run the model using multi-processing
def parallel_model(input_file_path, start_date, days, n_cores, shading):

    # Read the dataframe directly using pyogrio
    canals_gdf = read_dataframe(input_file_path)

    ''' Uncomment to enable quicker testing of validation canals '''
    #options = ['nabc_2'] # 'lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11'] 
    #canals_gdf = canals_gdf.loc[canals_gdf['code_id'].isin(options)]

    # Split into n clusters
    canals_gdf = array_split(canals_gdf, n_cores)

    # Initialises the Pool for multiprocessing (n cores)
    p = Pool(processes = n_cores)

    # Define fixed and variable arguments
    # Canal data, start date, days of analysis, shading boolean
    args = [(i, start_date, days, shading) for i in canals_gdf]

    # Generate results
    p.starmap(main, args)

# Main function
def main(canals_gdf, start_date, duration_days, model_shading):

    # Sets the duration of analysis (in hours)
    duration_hours = duration_days * 24

    # Include or exclude shading effects, updates path to file
    if model_shading:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # Skip canals with an area < 500 m^2 or > 200_000 m^2
        if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
            continue

        # If the file already exists, skip
        if isfile(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{feature_tuple.code_id}.json") or isfile(f"../test/results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{feature_tuple.code_id}.json"):
            print(f"model-output-{feature_tuple.code_id}.json exists")
            continue

        # Init dict for storing results
        else:
            output = {}

        # Include shading effect
        if model_shading:

            # Try and load shading record for current location
            try:
                with open(f"../../national-canal-cooling-data/shading/shading-{feature_tuple.code_id}.json") as shading_path:
                    shading_record = load(shading_path)

            # If the file is not found, search for test dataset
            except FileNotFoundError:
                try:
                    with open(f"../test/shading/shading-{feature_tuple.code_id}.json") as shading_path:
                        shading_record = load(shading_path)

                # Skip to next feature (the model cannot be run without shading data if shading == True)
                except FileNotFoundError:
                    print(f"Shading file MISSING for feature {feature_tuple.code_id}")
                    continue

        # Try and load climate record for current location
        try:
            with open(f"../../national-canal-cooling-data/interpolated-climate/data-{feature_tuple.code_id}.json") as climate_path:
                climate_record = load(climate_path)

        # If the file is not found, search for test dataset
        except FileNotFoundError:
            try:
                with open(f"../test/climate/climate-{feature_tuple.code_id}.json") as climate_path:
                        climate_record = load(climate_path)

            #  Skip to next feature (the model cannot be run without climate data)
            except FileNotFoundError:
                continue

        # Localise niave datetime to aware datetime (UK) and then convert to UTC
        local_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')

        # Sets end date
        end_datetime = local_datetime + timedelta(hours = duration_hours)
        
        # Transformer for projected to geographic coordinates (British National Grid | Irish National Grid > WGS84)
        transformer = Transformer.from_crs(canals_gdf.crs, 4326)

        # Feature centroid and area (shapely)
        latitude, longitude = transformer.transform(feature_tuple.geometry.centroid.x, feature_tuple.geometry.centroid.y)
        feature_area = feature_tuple.geometry.area

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

        # Init lists of starting values for the water (K = ~282, J = 0)
        water_temp = [MEAN_KELVIN] * len(split_points)
        water_energy = [0] * len(split_points)

        # Init list of starting value for reference material (K = ~282), comprising 10 depth layers [0 - .5m]
        reference_temp = [MEAN_KELVIN] * 10

        # Model spin up, based on a composite climate record for 2021-12-15
        water_temp, water_energy, reference_temp, _, _, _, = model_spin_up(
                                                                water_temp, #-------------------------- Water temperature (K)
                                                                water_energy, #------------------------ Water energy (J)
                                                                reference_temp, #---------------------- Reference temperature (K)
                                                                split_points, #------------------------ List of split points (m depth)
                                                                feature_area, #------------------------ Modelled area (m^2)
                                                                latitude, longitude, #----------------- Geographic coordinates
                                                                "../outputs/spin_up_climate.json", #--- Path to outfile file
                                                                30, #---------------------------------- Number of iterations
                                                                False, #------------------------------- Save output (Boolean)
                                                                feature_tuple.code_id) #--------------- Feature ID for output file
        
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

                # Return modelled water temperatures (K) and absorbed radiation (W)
                water_temp, absorbed_radiation = model_water_temperature(
                                                    local_datetime, #----------------- UTC datetime
                                                    longitude, latitude, #------------ Geographic coordinates (WGS84)
                                                    shading_proportion, #------------- Shading proportion [0-1]
                                                    water_energy, #------------------- List storing energy (J) for each layer
                                                    water_temp, #--------------------- List storing temperature (K) for each layer
                                                    split_points, #------------------- List of split points (m depth)
                                                    feature_area, #------------------- Modelled area (m^2)
                                                    climate['cloud_cover'] / 100) #--- Cloud cover (fractional)
                
                # Return modelled water energy (J) and related fluxes (W)
                water_energy, energy_values = model_energy_change(
                                                climate['air_temperature'], #--------- Air temperature (K)
                                                climate['humidity'], #---------------- Relative humidity (%)
                                                CLOUD_HEIGHT, #----------------------- Cloud height (fractional, constant)
                                                climate['cloud_cover'] / 100, #------- Cloud cover (fractional) 
                                                climate['wind_speed'], #-------------- Wind speed (m/s)
                                                climate['air_pressure'], #------------ Air pressure (hPa)
                                                water_temp, #------------------------- Modified water temperature (K), after model_water_temperature()
                                                feature_area, #----------------------- Modelled area (m^2)
                                                initial_water_k, #-------------------- Starting surface water temperature (K) for this model step
                                                absorbed_radiation) #----------------- Absorbed radiation (W), after model_water_temperature()
                
                # Model reference material temperature
                reference_temp, reference_energy_values = model_reference_temperature(
                                                    climate['air_temperature'], #----- Air temperature (K)
                                                    climate['humidity'],  #----------- Relative humidity (%)
                                                    CLOUD_HEIGHT, #------------------- Cloud height (fractional, constant)
                                                    climate['cloud_cover'] / 100, #--- Cloud cover (fractional)
                                                    climate['wind_speed'], #---------- Wind speed (m/s)
                                                    climate['air_pressure'], #-------- Air pressure (hPa)
                                                    reference_temp, #----------------- Reference temperature (K)
                                                    local_datetime, #----------------- UTC datetime 
                                                    longitude, latitude, #------------ Geographic coordinates (WGS84)
                                                    shading_proportion, #------------- Shading proportion [0-1]
                                                    feature_area) #------------------- Modelled area (m^2)

                # If modelling for a feature where there is measured data (Source: Canal and River Trust):
                if feature_tuple.code_id in ['lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11']:

                    # Add temperature (K) and energy values (W) to output dictionary, rounded to 3 dp
                    output.update({unix_time : {'surface_reference_k' : round(reference_temp[0], 3), #---------------------------- Reference temperature at surface
                                            'surface_water_k' : round(water_temp[0], 3), #---------------------------------------- Water temperature at surface (0 - 20 cm)
                                            'depth_water_20' : round(water_temp[1], 3), #----------------------------------------- Water temperature at 20 - 40 cm
                                            'depth_water_40' : round(water_temp[2], 3), #----------------------------------------- Water temperature at 40 - 60 cm
                                            'depth_water_60' : round(water_temp[3], 3), #----------------------------------------- Water temperature at 60 - 80 cm
                                            'depth_water_80' : round(water_temp[4], 3), #----------------------------------------- Water temperature at 80 - 100 cm
                                            'water_sensible' : round(energy_values['sensible_heat'], 3), #------------------------ Sensible flux associated with water
                                            'reference_sensible' : round(reference_energy_values['energy_to_air'] , 3)}}) #------- Sensible flux associated with reference material 
                
                # No measured data at this location
                else:

                    # Add temperature (K) and energy values (W) to output dictionary, all rounded to 3 dp
                    output.update({unix_time : {'surface_reference_k' : round(reference_temp[0], 3), #---------------------------- Reference temperature at surface
                                            'surface_water_k' : round(water_temp[0], 3), #---------------------------------------- Water temperature at surface (0 - 20 cm)
                                            'water_sensible' : round(energy_values['sensible_heat'], 3), #------------------------ Sensible flux associated with water
                                            'reference_sensible' : round(reference_energy_values['energy_to_air'] , 3)}}) #------- Sensible flux associated with reference material 

                # Finish the iteration (+ 15 minutes)
                local_datetime += timedelta(hours = 0.25)
        
        # For the test dataset
        if "test_" in feature_tuple.code_id:
            dump(output, open(f"../test/results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{feature_tuple.code_id}.json", 'w'))

        # Save the output to the relevant directory (including- | excluding-shading)
        else:
            dump(output, open(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{feature_tuple.code_id}.json", 'w'))

        # Summary statement
        print(f"Completed {feature_tuple.code_id}")

# The script is meant to be run
if __name__ == '__main__':

    # Select reference material
    user_input = input(f"The reference material is {REFERENCE_MATERIAL.upper()}. Do you want to continue? (yes/no): ")

    # Continue or exit
    if user_input.lower() in ["yes", "y"]:
        print("Continuing...")
    else:   
        exit("Reference material can be updated on Line 10 of params.py. Exiting...")
    
    # Initialise timer
    start = perf_counter()
    
    # Files for Britain and Ireland, and the test dataset
    file_names = ["../data/canal-geometries-filtered-modified-dissolved-id-width.shp", 
                  "../data/urban-canals-ireland-dissolved-id-tm65-width.shp",
                  "../test/test-canals.shp"]

    # Iterate through Britain and Ireland
    for f in file_names:

        # Run the parallel model (file path, start date, duration, n pools, shading [True | False])
        parallel_model(f, "2022-01-01 00:00:00", 365, 5, True)

    # Completion
    print(f"Code completion in {perf_counter() - start:.2f} seconds")

