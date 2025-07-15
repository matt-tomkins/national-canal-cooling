''' Usage: This code is not called by national-canal-cooling.py but does contain important data preparation functions. These are stored in auxiliary.py

# A few relevant sources
>---- gpd memory issues {1} ----> https://gis.stackexchange.com/questions/389027/reading-gdb-file-with-geopandas-throws-memoryerror
>---- gpd memory issues {2} ----> https://stackoverflow.com/questions/70258442/reading-large-gdb-file-with-geopandas
>---- Row loading using gpd ----> https://gis.stackexchange.com/questions/32762/accessing-feature-classes-in-file-geodatabases-using-python-and-gdal
>---- OGR ----------------------> https://gis.stackexchange.com/questions/205861/get-row-counts-of-all-tables-in-file-geodatabase-ideally-from-metadata
>---- Divisors -----------------> https://discuss.codecademy.com/t/how-can-we-get-all-the-divisors-of-a-number/376635
>---- Arange -------------------> https://stackoverflow.com/questions/68721736/divide-number-into-n-equal-parts
>---- Loading to postGIS -------> https://gis.stackexchange.com/questions/239198/adding-geopandas-dataframe-to-postgis-table
>---- postGIS intersects -------> https://geekswithlatitude.readme.io/docs/st_intersects
>---- multiprocessing SQL ------> https://docs.sqlalchemy.org/en/20/core/pooling.html#pooling-multiprocessing
'''

# Force use of Shapely 2.0
from os import environ
environ['USE_PYGEOS'] = '0'

# Required packages
import warnings
from time import perf_counter
from pyogrio import read_dataframe

# Import user-defined functions 
import config 
from auxiliary import * 
from params import *
from functions import *

# Prevent geopandas CRS warning for pybdshadow (Context: https://github.com/geopandas/geopandas/issues/2606)
# shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)
warnings.filterwarnings("ignore", message="CRS not set for some of the concatenation inputs")
warnings.filterwarnings('ignore', message="I don't know about leap seconds after 2021")

''' Load / filter spatial data and add to postGIS '''
LOAD_BUILDING = False #-----------> Load OS buildings
LOAD_LAND_COVER = False #---------> Load Land Cover Map 2021 (LCM2021)
LOAD_WATER = False #--------------> Filter and load OS MasterMap waterways

''' Process / classify canal geometries '''
CANAL_POINTS = False #------------> Generate points along the canal geometries (line)
CLASSIFY_URBAN = False #----------> Return the urban % for each point(LCM2021)
WEIGHT_URBAN = False #------------> Calculate moving average based on neighbouring points 
FILTER_URBAN = False #------------> Return points which exceed the chosen urban proportion
LINE_CONVERSION = False #---------> Convert canal points back to a line, based on a threshold distance

''' Return final canal geometries (polygon) for modelling '''
LOAD_CRT_CANALS = False #---------> Load CRT canal points 
RETURN_CANALS = False #-----------> Return canal polygons from postGIS
CANAL_HEIGHT = False #-------------> Determine the elevation of each canal (to incorporate into building shading calc)
GENERATE_ID = False #-------------> Generate a unique ID for each canal feature

''' Climate extraction and preprocessing '''
PRECOMPUTE_CLIMATE = False #------> Return climate data via API, add to JSON
PROCESS_CLIMATE = False #---------> Using the 40+ year climate record, extract the relevant values and store in JSON
FILTER_CLIMATE = False #----------> Filter the 40+ year record to 2022 only
CLIMATE_WEIGHTING = False #-------> Determine the climate location weighting for each canal feature
EXTRACT_CLIMATE = False #---------> Return unique climate data for each canal feature (based on weighting)
INTERPOLATE_CLIMATE = False #-----> Interpolates climate between each recorded hour (15-m intervals)

''' Building shading '''
PRECOMPUTE_SHADING = False #-------> Calculate shading for each feature, load to JSON
SHADING_COMPARISON = False #-------> Summarise shading for each feature, and return most shaded

''' Canal depth calculation '''
CANAL_DEPTH = False #-------------> Return canal temperature and energy values based on a given depth and interval
SPIN_UP = False #-----------------> Return average climate for 2022-01-01 across all locations - used for model spin up

''' Perform analysis for Irish datasets '''
IRELAND = False

''' File/directory paths '''

# Path to OS building geodatabases
building_path = r'../../canal-cooling/national/data/man-bha-gdb-2021-10/**/*.gdb'

# Path to OS MasterMap geodatabases (containing 'inland water')
water_path = [r"../../canal-cooling/national/data/man-topo-mm-1/**/*.gdb",
              r"../../canal-cooling/national/data/man-topo-mm-2/**/*.gdb",
              r"../../canal-cooling/national/data/man-topo-mm-3/**/*.gdb"]

# Path to LCM2021 shapefile (already filtered to urban areas)
land_cover_path = '.../canal-cooling/national/data/land-cover/lcm-2021-urban-classes.shp'

# Path to canals of interest (.shp)
canal_line_path = "../data/modelled-canals-2023-05-03.shp"
canal_points_path = "../data/modelled-canals-urban-weighted-50-info.shp"
canal_line_output = "../data/modelled-canals-urban-weighted-50-line.shp"

# Database records (JSON)
building_record = "../outputs/building_record.json"
water_record = "../outputs/water_record.json"
land_cover_record = "../outputs/land_cover_record.json"
canal_id_record = "../outputs/canal-id-values.json"
proximity_record = "../outputs/proximity_record.json"
height_record = "../outputs/height_record.json"
spin_climate = "../outputs/spin_up_climate.json"

# PostGIS settings
PW = config.info['postgis']
DB_STRING = f"postgresql://postgres:{PW}@localhost:5432/canals" # Connection string

# API key
MT_API = config.info['api_mt']
JH_API = config.info['api_jh']

def main():
    
    ''' These functions load / filter spatial data and add to postGIS '''

    # Load OS building data (.gdb) to postGIS, filtered based on canal network
    if LOAD_BUILDING:
        load_buildings(building_path, building_record, canal_points_path, calculated_db_size = 93_540, 
                       chunk_size = 10_000, buffer_distance = 1000)

    # Load LCM2021 urban land cover data (.shp) to postGIS ('chunk size' rows at a time)
    if LOAD_LAND_COVER:
        load_land_cover(land_cover_path, land_cover_record, chunk_size = 10_000, epsg = 27700)

    # Filter OS MasterMap (.gdb) to canal geometries and load to postGIS ('chunk size' rows at a time)
    if LOAD_WATER:
        filter_waterways(water_path, water_record, canal_line_output, chunk_size = 10_000, 
                         calculated_db_size = [7_133_916, 6_362_455, 12_771_572], # Max db size in each directory
                         buffer_distance = 1000, directory = 2) 

    ''' These functions process the canal geometries to determine which are to be modelled '''

    # Generate points along canal (line) geometries, every n metres ('interval')
    if CANAL_POINTS:
        generate_canal_points(canal_line_path, "../data/modelled-canals-points.shp", interval = 25, epsg = 27700)

    # Return the urban area within a specified distance from a point, based on LCM 2021
    if CLASSIFY_URBAN:
        urban_proportion("../data/modelled-canals-points.shp", 
                         "../data/modelled-canals-urban.shp", threshold_distance = 500, espg = 27700)

    # Return a weighted urban percentage, based on the points within a threshold distance
    if WEIGHT_URBAN:
        urban_weighting("../data/modelled-canals-urban.shp", 
                        "../data/modelled-canals-urban-weighted.shp", threshold_distance = 500,
                        epsg = 27700)
    
    # Returns points which exceed or equal the specified urban proportion (%)
    if FILTER_URBAN:
        urban_filter("../data/modelled-canals-urban-weighted.shp", threshold_proportion = 50, epsg = 27700)
        
    # Converts canal points back to line
    if LINE_CONVERSION:
        canal_line_conversion(canal_points_path, canal_line_output, threshold_distance = [10, 35])

    ''' These functions extract the canal geometries (polygon) from postGIS '''

    # Load CRT canals file (.shp, Line) to postGIS
    if LOAD_CRT_CANALS:
        load_crt_canals(canal_line_output) # Output of canal_line_conversion() + manual editing

    # Return canal geometries stored on postGIS
    if RETURN_CANALS:
        return_waterways("../data/canal-geometries-filtered.shp")

    # Return the elevation of each canal
    if CANAL_HEIGHT:
        calculate_canal_height("../data/canal-geometries-filtered-modified-dissolved-id.shp", height_record)

    # Calculate shaded area for each feature (buildings + bank), load to postGIS
    if GENERATE_ID:
        generate_unique_id("../data/canal-geometries-filtered-modified-dissolved.shp", canal_id_record)

    ''' This function calculates the shaded area for each canal feature for the desired duration '''

    # Calculate shaded area for each feature (building shadows), save to json
    if PRECOMPUTE_SHADING:
        
        # set start time
        start_time = perf_counter() 

        # Run for Great Britain and Ireland
        for file in ["../data/canal-geometries-filtered-modified-dissolved-id.shp", "../data/urban-canals-ireland-dissolved-id-tm65.shp"]:

            # Sets input file, start date, # cores, duration of analysis and distance (d) to obtain buildings
            shading_parallel(file, # file path
                            "2022-01-01 00:00:00", # start date
                            n_cores = 8, # cores
                            days = 365, # duration of analysis
                            d = 500) # distance threshold for buildings (m)
        
        # Example of non-parallel approach
        #options = ['gc1_22']
        #canals_gdf = canals_gdf.loc[canals_gdf['code_id'].isin(options)]
        #precompute_shading(canals_gdf, "2022-01-01 00:00:00", days = 365, dist = 500)
        
        # Completion time
        print(f"Completed in: {perf_counter() - start_time} seconds")

    if SHADING_COMPARISON:

        # Summarise shading for each feature, and return a list where annual shading >20%
        extract_highly_shaded(0.2)

    ''' These function process climate data from OpenWeatherMap '''

    # Process 40-year climate data from Open Weather Map
    if PROCESS_CLIMATE:

        # Input and output file path
        preprocess_climate_data(r'../../national-canal-cooling-data/raw-climate/*.json', "../../national-canal-cooling-data/processed-climate/")

    # Subset climate to a specific time period (in this case, 2022)
    if FILTER_CLIMATE:

        # Input and output file path, start and end datetime 
        filter_climate_data(r'../../national-canal-cooling-data/processed-climate/*.json', "../../national-canal-cooling-data/filtered-climate/",
                            "2022-01-01 00:00:00", "2023-01-01 00:00:00")
    
    # Relative weighting of climate location data, based on distance
    if CLIMATE_WEIGHTING:

        # Inputs: file path for climate locations and canal, output path, distance (m), write to postGIS (True/False)
        calc_climate_weighting("../data/climate-locations-kmean-100.shp", "../data/canal-geometries-filtered-modified-dissolved-id.shp", 
                               proximity_record, dist = 15_000, write = False)
    
    # Return climate data for each feature based on climate location weighting
    if EXTRACT_CLIMATE:

        extract_climate_data("../data/canal-geometries-filtered-modified-dissolved-id.shp", # Canal geometries
                             proximity_record, # Climate weighting for each geometry
                             "../../national-canal-cooling-data/filtered-climate/", # Climate data for each location
                             "../../national-canal-cooling-data/feature-climate/") # Output directory


    # Return climate data for each feature based on climate location weighting
    if INTERPOLATE_CLIMATE:

        # Analysis hard-coded for 15 minute intervals
        interpolate_climate_data(r'../../national-canal-cooling-data/feature-climate/*.json', # Input files (1-hour resolution)
                                 "../../national-canal-cooling-data/interpolated-climate/") # Output files (15-minute resolution)

    # Return climate data for each location and each datetime, load to JSON, starting "2022-06-10 00:00:00"
    # Progress: as of __16/06/2023__, we have loaded __365__ days of data
    if PRECOMPUTE_CLIMATE:
        '''
        > Approach to obtain climate data via the OpenWeather API
        > No longer required, as funding was obtained for bulk-download of 40-year data
        '''

        # Days loaded so far, total possible per day (40 days, 100k API calls), number of days per iteration
        completed_days = 325
        total = 40
        interval = 20
        
        # List of durations
        durations = range(completed_days + interval, completed_days + interval + total, interval)

        # Iterate through the durations
        for d in durations:

            # For the first 20 days of analysis
            if d <= completed_days + (total / 2):

                # Obtain the climate data (parallel processing) using this input file, 5 cores, d days of analysis and the first API key
                climate_parallel("../data/climate-locations-kmean-100.shp", "2022-06-10 00:00:00", 5, d, MT_API)
            
            # For the final 20 days
            else:
                
                # Obtain the climate data (parallel processing) using this input file, 5 cores, d days of analysis and the second API key
                climate_parallel("../data/climate-locations-kmean-100.shp", "2022-06-10 00:00:00", 5, d, JH_API)
         
    ''' These functions are incomplete '''
    # Process climate data
    if SPIN_UP:

        # Generate a composite climate record for the specified date
        spin_up_climate("2021-12-15 00:00:00", "2021-12-16 00:00:00", "../outputs/spin_up_climate.json")

    ''' Functions for processing Ireland data'''
    if IRELAND:

        # EPSG code for TM65 / Irish Grid
        tm_65 = 29902
        
        # Load land cover data to postGIS
        for f in ["../../national-canal-cooling-data/ireland/lcm-2021-NI-urban-tm65.shp", "../../national-canal-cooling-data/ireland/lcm-2018-ireland-urban-tm65.shp"]:
            load_land_cover(f, land_cover_record, chunk_size = 10_000, epsg = tm_65)

        # Generate canal points 
        generate_canal_points("../../national-canal-cooling-data/ireland/canals-ireland-lines-tm65.shp", 
                              "../../national-canal-cooling-data/ireland/canals-ireland-points-tm65.shp", 
                              interval = 25, epsg = tm_65)

        # Calculate urban proportion
        urban_proportion("../../national-canal-cooling-data/ireland/canals-ireland-points-tm65.shp", 
                         "../../national-canal-cooling-data/ireland/canals-ireland-urban-tm65.shp", 
                         threshold_distance = 500, epsg = tm_65)

        # Weighted average 
        urban_weighting("../../national-canal-cooling-data/ireland/canals-ireland-urban-tm65.shp", 
                        "../../national-canal-cooling-data/ireland/canals-ireland-urban-weighted-tm65.shp", 
                        threshold_distance = 500, epsg = tm_65)

        # Returns points which exceed or equal the specified urban proportion (%)
        urban_filter("../../national-canal-cooling-data/ireland/canals-ireland-urban-weighted-tm65.shp", 
                     threshold_proportion = 50, epsg = tm_65)

        # Irish IDs 
        generate_unique_id_ireland("../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-tm65.shp", 
                                    canal_id_record, 
                                    "../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp")

        # Add canals to postGIS
        db_connection = create_engine(DB_STRING)
        canals = read_dataframe("../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp")
        canals.to_postgis('canals_ireland', db_connection, if_exists='append', dtype={'geom': Geometry('[Polygon]', srid=tm_65)})
        
        # Add building heights to OSM polygons based on Keany et al. (2022)
        add_building_heights("../../national-canal-cooling-data/ireland/buildings-ireland-tm65.shp",
                             "../../national-canal-cooling-data/ireland/building_heights_keany_2022_tm65.tif",
                         "../../national-canal-cooling-data/ireland/buildings-heights-ireland-tm65.shp")
        

        # Add missing building heights to OSM polygons
        missing_building_heights("../../national-canal-cooling-data/ireland/buildings-heights-ireland-tm65.shp",
                         "../../national-canal-cooling-data/ireland/buildings-heights-missing-ireland-tm65.shp")

        
        # Add buildings to postGIS, formatting to match Ordnance Survey 
        load_buildings_ireland("../../national-canal-cooling-data/ireland/buildings-heights-missing-ireland-tm65.shp", chunk_size = 10_000)
        
        # Test to see if return_buildings() still works 
        canals = read_dataframe("../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp")
        options = ["gc1_35"]
        canals = canals.loc[canals['code_id'].isin(options)]
        for feature_tuple in canals.itertuples():
            return_buildings_polygon(feature_tuple, 500, True, canals.crs)
        
        # Subdivide climate data into clusters 
        clusters = ['Cluster-100', 'Cluster-101', 'Cluster-102', 'Cluster-103', 'Cluster-104', 'Cluster-105', 'Cluster-106', 'Cluster-107'] 
        subdivide_climate_data(f"C:/Users/44797/Downloads/c392ebda1af4b2d0c3de8c695c31da52.json", clusters)
        
        # Preprocessing, input and output file path 
        preprocess_climate_data(r'../../national-canal-cooling-data/raw-climate/*.json', "../../national-canal-cooling-data/processed-climate/")
        
        # Filter to 2022, Input and output file path, start and end datetime 
        filter_climate_data(r'../../national-canal-cooling-data/processed-climate/*.json', "../../national-canal-cooling-data/filtered-climate/",
                            "2022-01-01 00:00:00", "2023-01-01 00:00:00")
        

        # Inputs: file path for climate locations and canal, output path, distance (m), write to postGIS (True/False)
        calc_climate_weighting("../../national-canal-cooling-data/ireland/climate-locations-ireland-kmean-8.shp", 
                               "../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp", 
                               "../outputs/proximity_record_ireland.json", dist = 15_000, write = True)
        
        
        # Generate a unique climate record for each canal feature
        extract_climate_data("../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp", # Canal geometries
                                    "../outputs/proximity_record_ireland.json", # Climate weighting for each geometry
                                    "../../national-canal-cooling-data/filtered-climate/", # Climate data for each location
                                    "../../national-canal-cooling-data/feature-climate/") # Output directory
        

        # Analysis hard-coded for 15 minute intervals
        interpolate_climate_data(r'../../national-canal-cooling-data/feature-climate/*.json', # Input files (1-hour resolution)
                                 "../../national-canal-cooling-data/interpolated-climate/") # Output files (15-minute resolution)
        
        # Summary statement
        print("All functions completed on 15-07- and 16-07-2023")
  
# The script is meant to be run
if __name__ == '__main__':

    # Initialise timer and run script
    start_time = perf_counter() 
    main()

    # Code completion
    print(f"Code completion in {perf_counter()-start_time:.2f} seconds!")