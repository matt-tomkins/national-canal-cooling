''' Usage: This code is not called by national-canal-cooling.py, but is used by auxiliary-code.py for data preparation '''

# Required packages
import warnings
import requests
from pyproj import Geod
from glob import glob as g 
from psycopg2 import errors
from rtree.index import Index
from fiona import open as fopen 
from geoalchemy2 import Geometry
from multiprocessing import Pool
from os.path import split, isfile
from pyogrio import read_dataframe
from rasterstats import zonal_stats
from collections import defaultdict
from time import perf_counter, sleep
from sqlalchemy.pool import NullPool
from math import tan, ceil, pi, isnan
from sqlalchemy import create_engine, text
from json import dump, load, loads, decoder
from datetime import datetime, timedelta, date
from math import degrees, radians, hypot, floor
from skimage.draw import line, circle_perimeter
from pysolar.solar import get_altitude, get_azimuth
from rasterio.features import rasterize, Affine, shapes
from rasterio import open as rast_open, uint16, float32
from shapely.geometry import Point, LineString, Polygon
from matplotlib.pyplot import subplots, savefig, title, close
from pandas import to_datetime, to_numeric, DataFrame, merge, concat
from geopandas import GeoDataFrame, GeoSeries, clip, sjoin, read_file as gpd_read
from numpy import zeros, nonzero, argwhere, column_stack, where, sin, cos, intersect1d, arange, array_split, unique, quantile, median, mean

# Import user-defined functions 
import config
from shading import *

# Prevent geopandas CRS warning for pybdshadow (Context: https://github.com/geopandas/geopandas/issues/2606)
# shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)
warnings.filterwarnings("ignore", message="CRS not set for some of the concatenation inputs")

# PostGIS password
PW = config.info['postgis']

# String for connection
DB_STRING = f"postgresql://postgres:{PW}@localhost:5432/canals" # Connection string on local database
# DB_STRING = f"postgresql://matt:{PW}@localhost:5432/canals" # Connection string on server

# API key
MT_API = config.info['api_mt']
JH_API = config.info['api_jh']

# Sets ellipsoid model
g_model = Geod(ellps="WGS84")

#============== Functions ==============#

def partition_v2(database_size, interval_size):
    '''
    > # Function to split a database into chunks based on an interval size
    '''
    # Number of intervals
    n = int(database_size / interval_size)

    # Return list of start and end points using list comprehension
    return [(0 + i * interval_size, interval_size + i * interval_size - 1) for i in range(n)]
    

def round_up(input, value):
    ''' Function to round up inputs to the nearest 'value' (e.g. 100th, 1000th)'''

    return int(ceil(input / value)) * value

def round_down(input, value):
    ''' Function to round down inputs to the nearest 'value' (e.g. 100th, 1000th)'''

    return int(floor(input / value)) * value

def load_buildings(building_path, building_record, canal_file_path, calculated_db_size, chunk_size, buffer_distance):
    '''
    > A function to load buildings from the OS MasterMap (stored as gdb) to a postGIS database
    > 9973 geodatabases, ~23.5 min code run time
    > 3_154_113 polygons stored in 'buildings' postGIS table
    '''

    # Try
    try:
        
        # Load record of geodatabases that have been accessed previously
        with open(building_record) as gd_path:
            database_record = load(gd_path)

    # If the file is not found
    except FileNotFoundError:

        # Create a new dictionary
        database_record = {}

    # File paths for all gdb files within man-bha-gdb-2021-10, 
    building_files = g(building_path, recursive = True)

    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Variable for storing the maximum number of rows
    maximum_database_size = 0

    # Read using gpd and pyogrio, and dissolve points to a single feature
    canal_points = gpd_read(canal_file_path, engine = "pyogrio").dissolve()

    # Buffer the canal geometry (points) by a set distance (1000 m)
    buffered_canal_points = GeoDataFrame(geometry = canal_points.buffer(buffer_distance).geometry, crs = 27700)

    # First iteration of the function
    if calculated_db_size is None: 

        # Iterate through the files
        for b in building_files:

            # Using fiona
            with fopen(b) as features:

                # Return the number of features
                database_size = len(features)

                # Stores if maximum
                if database_size > maximum_database_size:
                    maximum_database_size = database_size

    # If the function has been used before       
    else:

        # The maximum database size has been calculated
        maximum_database_size = calculated_db_size
    
    # Round up database size to nearest chunk size (10,000)
    maximum_database_size = round_up(maximum_database_size, chunk_size)

    # Create list of start-end points
    db_chunks = partition_v2(maximum_database_size, chunk_size)
    
    # Iterate through the geodatabases
    for b in building_files:

        # Get the file name from the full path
        file_name = split(b)[1]

        # If we've accessed this gdb before
        if file_name in database_record.keys():

            # Continue to next file 
            # print(f"Database {file_name} was loaded on {database_record[file_name]}")
            continue
        
        # Otherwise, open the file in chunks
        else: 

            # Iterate through the chunks
            for chunk in db_chunks:
                
                # Use pyogrio
                try:
            
                    # Read using gpd and pyogrio
                    building_gdf = gpd_read(b, 
                                            rows = slice(chunk[0], chunk[1]), # Slice database based on chunk start-end
                                            engine = "pyogrio") # For efficiency
                
                # ValueError: 'skip_features' must be between 0 and dataset length
                except ValueError:

                    # Break this loop (skip to next file)
                    break

                # If the gdf is empty
                if building_gdf.empty:
                    
                    # Break this loop (skip to next file)
                    break

                # If there are buildings, first clip the canal geometry by the total bounds of the building data
                clipped_canal_buffer = clip(buffered_canal_points, mask = building_gdf.total_bounds)

                # Then clip the building geometry by the canal buffer
                filtered_building_gdf = clip(building_gdf, clipped_canal_buffer)

                # If the gdf is empty
                if filtered_building_gdf.empty:

                    # Skip to next chunk
                    continue
                
                # If there is data
                else:
                    
                    # Try to add to postGIS
                    try:

                        # Add to postGIS
                        filtered_building_gdf.to_postgis('buildings', db_connection, if_exists='append', 
                                                        dtype={'geom': Geometry('[Polygon]', srid=27700)})
                    
                    # If multipolygons have been produced:
                    except errors.InvalidParameterValue: 

                        # Multi-part geometries into multiple single geometries
                        exploded_building_gdf = filtered_building_gdf.explode(index_parts = False)

                        # Then, add to postGIS
                        exploded_building_gdf.to_postgis('buildings', db_connection, if_exists='append', 
                                dtype={'geom': Geometry('[Polygon]', srid=27700)})
    
            # If successful, store the file name and date (as string)
            database_record.update({file_name : datetime.today().strftime('%d-%m-%Y')}) 
        
    # When complete, save the buildings record file (json)
    dump(database_record, open(building_record, 'w')) 
    
def load_crt_canals(canal_file_path):
    '''
    > A function to load CRT canal geometries (LineStrings) to postGIS
    '''

    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Load the file using geopandas
    canal_line_gdf = gpd_read(canal_file_path, engine = "pyogrio")
                                
    try: 

        # Add to postGIS database (crt_canals)
        canal_line_gdf.to_postgis('canal_lines', 
                                    db_connection, if_exists='fail', # Raise error if table already exists
                                    dtype={'geom': Geometry('[Line]', srid=27700)})
        
    except ValueError:

        print("Error: The table already exists")

def load_land_cover(file_path, land_cover_record_path, chunk_size, epsg):
    '''
    > A function to load the 2021 Land Cover Map 2021 (LCM2021) to a postGIS database
    > Layer has been filtered to focus on urban areas only (Classes 20 and 21)
    '''

    # Try
    try:
        
        # Load record of geodatabases that have been accessed previously
        with open(land_cover_record_path) as gd_path:
            database_record = load(gd_path)

    # If the file is not found
    except FileNotFoundError:

        # Create a new dictionary
        database_record = {}

    # If we've accessed this gdb before
    if file_path in database_record.keys():

        # Curtail function
        return print(f"Database {file_path} was loaded on {database_record[file_path]}")

    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Using fiona
    with fopen(file_path) as features:

        # Return the number of features
        database_size = len(features)

    # Round up database size to nearest chunk size (10,000)
    maximum_database_size = round_up(database_size, chunk_size)

    # Create list of start-end points
    db_chunks = partition_v2(maximum_database_size, chunk_size)

    # Iterate through the chunks
    for chunk in db_chunks:
    
        # Read using gpd and pyogrio
        land_cover_gdf = gpd_read(file_path, 
                                    rows = slice(chunk[0], chunk[1]), # Slice database based on chunk start-end
                                    engine = "pyogrio") # Use vectorised Pyogrio, rather than Fiona
        
        # If the gdf is empty
        if land_cover_gdf.empty:
            
            # Finish
            break
        
        # For Britain
        if epsg == 27700:

            # Add to postGIS
            land_cover_gdf.to_postgis('land_cover', db_connection, if_exists='append', 
                                            dtype={'geom': Geometry('[Polygon]', srid=epsg)})
        
        # For Ireland
        elif epsg == 29902:

            # Add to postGIS
            land_cover_gdf.to_postgis('land_cover_ireland', db_connection, if_exists='append', 
                                            dtype={'geom': Geometry('[Polygon]', srid=epsg)})
    

    # If successful, store the file name and date (as string)
    database_record.update({file_path : datetime.today().strftime('%d-%m-%Y')}) 

    # When complete, save the land cover record file (json)
    dump(database_record, open(f"../outputs/land_cover_record.json", 'w')) 

def filter_waterways(water_path, water_record, canal_file_path, chunk_size, calculated_db_size, buffer_distance, directory):
    '''
    > A function to filter through the OS MasterMap (stored as gdb) and return inland waterways
    > These are then loaded to postGIS
    > To avoid duplication, database names are stored in a dictionary (JSON) - this is then checked for each file
    > Rather than loading in defined chunks (based on divisible numbers), we sequence into chunks of chunk_size (10,000)
    > As of 16-05-2023, this has been simpified to add everything classed as "inland water", irrespective of canal proximity
    '''

    # List of 'potentially relevant' water classifications, stored in "descriptiveterm" attribute
    # This is then filtered based on the CRT polyline of canals, so we can include more marginal categories for safety
    '''
    categories = ["Aqueduct, Canal", "Canal", "Canal Feeder", "Canal, Aqueduct", "Canal, Reeds",
                  "Conduit", "Ford", "Lock", "Lock, Static Water", "Reservoir", "Reservoir, Reeds", "Static Water",
                  "Watercourse", "Waterfall", "Reeds, Static Water", "Canal, Conduit"]

    # Fields to retain = ['descriptiveterm', 'descriptivegroup', 'geometry', 'fid', 'versiondate']
    fields = ['primary_key', 'featurecode', 'version', 'theme', 'calculatedareavalue', 
              'changedate', 'reasonforchange', 'make', 'physicallevel', 'physicalpresence', 'poly_broken', 
              'session_start_date', 'user_start_date', 'session_end_date', 'user_end_date']

    # Relevant classification, stored in "descriptive group" attribute
    categories = ["Inland Water"]
    '''
    
    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Try
    try:
        
        # Load record of geodatabases that have been accessed previously
        with open(water_record) as gd_path:
            database_record = load(gd_path)

    # If the file is not found
    except FileNotFoundError:

        # Create a new dictionary
        database_record = {}
    
    # File paths for all gdb files within each directory [0, 1, 2]
    water_files = g(water_path[directory], recursive = True)

    # Read the canal points using gpd and pyogrio, and dissolve points to a single feature
    # canal_line = gpd_read(canal_file_path, engine = "pyogrio").dissolve()

    # Buffer the canal geometry (points) by a set distance (1000 m)
    # buffered_canal_line = GeoDataFrame(geometry = canal_line.buffer(buffer_distance).geometry, crs = 27700)
    
    # If this is the first time running the algorithm
    if calculated_db_size is None:

        # Variable for storing the maximum number of rows
        maximum_database_size = 0

        # First iterate through the files
        for w in water_files:

            # Using fiona
            with fopen(w) as features:

                # Return the number of features
                database_size = len(features)
            
                # Stores if maximum
                if database_size > maximum_database_size:
                    maximum_database_size = database_size

    else:

        # Use the previously calculated maximum
        maximum_database_size = calculated_db_size[directory]

    # Round up database size to nearest chunk size (10,000)
    maximum_database_size = round_up(maximum_database_size, chunk_size)

    # Create list of start-end points
    db_chunks = partition_v2(maximum_database_size, chunk_size)

    # Iterate through the files
    for w in water_files:
        
        # Start time
        start = perf_counter() 

        # Get the file name from the full path
        file_name = split(w)[1]
        
        # If we've accessed this gdb before
        if file_name in database_record.keys():

            # Skip to the next
            print(f"Database {file_name} was loaded on {database_record[file_name]}")
            continue
        
        # Otherwise, this is the first time loading this database, iterate through the chunks
        for chunk in db_chunks:

            # Try and use pyogrio
            try:

                # Read directly using pyogrio
                water_gdf = read_dataframe(w, 
                                           layer = 'topographicarea', # Return the polygons
                                           columns = ['descriptiveterm', 'descriptivegroup', 'geometry', 'fid', 'versiondate'],
                                           skip_features = chunk[0],
                                           max_features = 10_000)
            
                
            # ValueError: 'skip_features' must be between 0 and dataset length
            except ValueError:

                # Break this loop (skip to next file)
                break
            
            # If the gdf is empty (probably defunct)
            if water_gdf.empty:
                
                # Skip to the next database
                break
            
            # If the gdf contains data, filter to "inland water" only
            filtered_water_gdf = water_gdf.loc[water_gdf['descriptivegroup'] == 'Inland Water']

            # If there is no inland water within this chunk
            if filtered_water_gdf.empty:
            
                # Skip to the next chunk
                continue

            # If there are waterways, first clip the canal geometry by the total bounds of the inland water data
            # clipped_canal_buffer = clip(buffered_canal_line, mask = filtered_water_gdf.total_bounds)

            # Then clip the water geometry by the canal buffer
            # clipped_water_gdf = clip(filtered_water_gdf, clipped_canal_buffer)

            # Sort by date, remove duplicates based on feature id, keep most recent update
            # clipped_water_gdf = clipped_water_gdf.sort_values('versiondate').drop_duplicates(subset=['fid'], keep='last')
            
            # Add to postGIS
            try: 
                filtered_water_gdf.to_postgis('water', db_connection, if_exists='append', 
                                            dtype={'geom': Geometry('[Polygon]', srid=27700)})
            
            # If multipolygons have been produced:
            except errors.InvalidParameterValue: 

                # Multi-part geometries into multiple single geometries
                exploded_water_gdf = filtered_water_gdf.explode(index_parts = False)

                # Then, add to postGIS
                exploded_water_gdf.to_postgis('water', db_connection, if_exists='append', 
                        dtype={'geom': Geometry('[Polygon]', srid=27700)})
                
            
        # If successful, store the file name and date (as string)
        database_record.update({file_name : str(datetime.now())}) 
        
        # End time and duration
        end = perf_counter() 
        print(f"It took {end - start:.2f} seconds to filter and load database {w}")
        
    # When complete, save the geodatabase record file (json)
    dump(database_record, open(water_record, 'w'))  

def generate_canal_points(input_file_path, output_file_path, interval, epsg):
    '''
    > A function to generate points at a set interval (m) along a line
    '''
    
    # Load the canal geometries using gpd and pyogrio      
    watercourse_line = gpd_read(input_file_path, 
                                     engine = "pyogrio")

    # Dissolve into a single feature
    merged_line = watercourse_line.dissolve()

    # Convert the length of the dissolved line to an integer, and generate corresponding intervals using numpy
    length = int(merged_line.iloc[0].geometry.length)
    distances = arange(0, length, interval)

    # Interpolate points along line using list comprehension
    points = [merged_line.interpolate(distance) for distance in distances]

    # Convert to GeoSeries using list comprehension
    canal_points = GeoSeries([item[0] for item in points], crs=epsg)

    # Output to file
    canal_points.to_file(output_file_path)

def urban_proportion(canal_points_path, canal_points_urban_path, threshold_distance, epsg):
    '''
    > A function to iterate through canal points -
    > Return urban areas (LCM 2021) within a specified threshold distance
    > Calculate the area as a proportion of the threshold area

    > https://postgis.net/docs/ST_GeomFromText.html
    > https://postgis.net/docs/ST_MakePoint.html
    '''

    # Area of a circle based on radius (A=πr2), 785398.1633974483 m2
    potential_area = pi * threshold_distance ** 2

    # Create the engine
    db_connection = create_engine(DB_STRING)

    # Load the canal points using gpd and pyogrio      
    canal_points = gpd_read(canal_points_path, 
                                     engine = "pyogrio")
    
    # Store the values in an output dictionary
    output = {}

    # Iterate efficiently using itertuples
    for point_tuple in canal_points.itertuples():
        
        # Convert to GeoSeries and project
        point_geometry_projected = GeoSeries(point_tuple.geometry, crs = epsg)

        # For Britain
        if epsg == 27700:

            # Return the intersection of land cover geometries and a Point (buffered by threshold distance)
            query = text("select * from land_cover where ST_Intersects(land_cover.geometry,ST_Buffer(ST_SetSRID(ST_MakePoint(:e, :n), :crs), :c))") \
            .bindparams(n = point_geometry_projected.iloc[0].y, # Northing
                        e = point_geometry_projected.iloc[0].x, # Easting
                        d = threshold_distance, # Distance (m)
                        c = epsg) # 27700
        
        # For Ireland
        elif epsg == 29902:
            
            # Return the intersection of land cover geometries and a Point (buffered by threshold distance)
            query = text("select * from land_cover_ireland where ST_Intersects(land_cover_ireland.geometry,ST_Buffer(ST_SetSRID(ST_MakePoint(:e, :n), :c), :d))") \
            .bindparams(n = point_geometry_projected.iloc[0].y, # Northing
                        e = point_geometry_projected.iloc[0].x, # Easting
                        d = threshold_distance, # Distance (m)
                        c = epsg) # 29902

        
        # Extracts geometries using from_postgis, stores as geodataframe
        land_cover_via_python = GeoDataFrame.from_postgis(sql = query, 
        con = db_connection,
        geom_col='geometry', # Column name for the geometry
        index_col='id', # Column name for the unique ID
        coerce_float=True)

        # If there are no urban areas <= threshold distance (intersects)
        if land_cover_via_python.empty:

            # Update the dictionary
            output.update({point_tuple.FID : {'FID' : point_tuple.FID,
                                              'urban_percentage' : 0,
                                              'urban_area' : 0}}) 

        # If there are 'potential' urban areas:
        else:

            # Draw calculation area (circle)
            circle = GeoSeries(Point(point_geometry_projected.iloc[0].x, point_geometry_projected.iloc[0].y).buffer(threshold_distance), crs = epsg)

            # Clip geometries
            land_cover_clipped = land_cover_via_python.clip(circle)

            # If there are no urban areas <= threshold distance (clipped)
            if land_cover_clipped.empty:

                # Update the dictionary
                output.update({point_tuple.FID : {'FID' : point_tuple.FID,
                                                  'urban_percentage' : 0,
                                              'urban_area' : 0}}) 
                
            # There are urban areas that fall within the clipped area
            else:
            
                # Iterate through geometries, store area using list comprehension and sum()
                total_area = sum([cover.geometry.area for cover in land_cover_clipped.itertuples()])
  
                # Percentage urban land cover
                urban_percent = (total_area / potential_area) * 100

                # Update the dictionary, storing the urban area (m^2) and percentage. 
                output.update({point_tuple.FID : {'FID' : point_tuple.FID,
                                                  'urban_percentage' : urban_percent,
                                              'urban_area' : total_area}}) 
    
    # Converts output to dataframe
    df = DataFrame.from_dict(output, 'index', columns=['FID', 'urban_percentage', 'urban_area'])

    # Joins to geodataframe
    canal_points = merge(canal_points, df, on=['FID'])

    # Output to file
    canal_points.to_file(canal_points_urban_path)


def urban_weighting(canal_points_input_path, canal_points_urban_output_path, threshold_distance, epsg):
    '''
    > A function to iterate through the canal points and their urban area %
    > Extract points within a specified threshold distance
    > Calculate the average of that
    > AIM: to avoid short stretches of non-urban canals within urban areas (and vice versa)
    '''

    # Load the canal points using gpd and pyogrio      
    canal_points = gpd_read(canal_points_input_path, 
                                     engine = "pyogrio")
    
    # initialise an rtree Index object
    idx = Index()

    # Store the values in an output dictionary
    output = {}

    # Iterate through the canals gdf and load into the spatial index
    for row in canal_points.itertuples():

        # Add the ID and geometry to the spatial index
        idx.insert(row.FID, (row.geometry.x, row.geometry.y)) 

    # Iterate through the canals gdf when complete
    for row in canal_points.itertuples():

        # Draw calculation area (circle)
        circle = GeoSeries(Point(row.geometry.x, row.geometry.y).buffer(threshold_distance), crs = epsg)

        # Return the points within the threshold distance (bounding geometry)
        potential_points = canal_points.iloc[list(idx.intersection((circle.geometry.bounds.iloc[0])))]

        # Return the points within the threshold distance (actual geometry)
        actual_points = potential_points.loc[potential_points.within(circle.geometry.iloc[0])]

        # Update the dictionary, storing the average urban area
        output.update({row.FID : {'FID' : row.FID,
                                        'avg_ub_per' : actual_points["urban_perc"].mean()}}) 
        
    # Converts output to dataframe
    df = DataFrame.from_dict(output, 'index', columns=['FID', 'avg_ub_per'])

    # Joins to geodataframe
    canal_points = merge(canal_points, df, on=['FID'])

    # Output to file
    canal_points.to_file(canal_points_urban_output_path)

def urban_filter(input_path, threshold_proportion, epsg):
    '''
    > A function to return points which exceed or equal the specified urban proportion (%)
    '''

    # Load the canal points using gpd and pyogrio      
    canal_points = gpd_read(input_path, 
                                     engine = "pyogrio")
    
    # Select rows based on condition
    filtered_canal_points = canal_points[canal_points['avg_ub_per'] >= threshold_proportion]

    # Britain
    if epsg == 27700:

        # Output to file
        filtered_canal_points.to_file(f"../data/modelled-canals-urban-weighted-{threshold_proportion}.shp")

    # Ireland
    elif epsg == 29902:

        # Output to file
        filtered_canal_points.to_file(f"../../national-canal-cooling-data/ireland/canals-ireland-urban-weighted-{threshold_proportion}-tm65.shp")


def canal_line_conversion(canal_points_path, canal_line_output, threshold_distance):
    '''
    > A function to convert canal points to lines
    > This is more complex than a simple Point to LineString conversion
    > We iterate through the point and return the n nearest neighbours
    > If the distance is within the correct threshold (> 50 m and < 100 m), we draw a Line between them
    '''

    # Load the canal points using gpd and pyogrio      
    canal_points = gpd_read(canal_points_path, 
                                     engine = "pyogrio")
    
    # initialise an rtree Index object
    idx = Index()

    # Iterate through the canals gdf and load into the spatial index
    for row in canal_points.itertuples():

        # Add the ID and geometry to the spatial index
        idx.insert(row.FID, (row.geometry.x, row.geometry.y)) 

    # Create empty GeoDataframe
    output_gdf = GeoDataFrame(columns=['FID', 'urban_perc', 'urban_area', 'avg_ub_per',
                                       'Waterway', 'layer', 'geometry'], crs = 27700)
    
    # Create an empty set for storing visited points
    point_list = set()

    # Iterate through the canals gdf
    for row in canal_points.itertuples():
        
        # Return n nearest points using index
        nearest_points = list(idx.nearest(row.geometry.bounds, 4))

        # Remove self intersections
        nearest_points.remove(row.FID)

        # Add the row ID to the point list set()
        point_list.add(row.FID)

        # Iterate through the point indexes
        for point in nearest_points:
            
            # To avoid repeat calculations
            if point in point_list:
            
                # Skip to next
                continue

            # Extract the relevant points
            extracted_point = canal_points.loc[canal_points['FID'] == point]

            # Calculates the distance between two points using Shapely
            distance = row.geometry.distance(extracted_point.geometry)

            # If the distance is between the threshold_distances
            if distance.iloc[0] < threshold_distance[1] and distance.iloc[0] > threshold_distance[0]:

                # Create a LineString based on the two points
                joined_line = LineString([row.geometry, extracted_point.geometry.iloc[0]])

                # Stores attributes
                point_attributes = DataFrame({
                        "FID": [row.FID],
                        "urban_perc": [row.urban_perc],
                        "urban_area": [row.urban_area],
                        "avg_ub_per": [row.avg_ub_per],
                        "Waterway": [row.Waterway],
                        "layer": [row.layer],
                        "draught": [row.draught]})
            
                # Converts to gdf
                row_to_append = GeoDataFrame(point_attributes, geometry = [joined_line], crs = 27700)
                
                # Appends to output table
                output_gdf = concat([output_gdf, row_to_append])

    # When all complete, save to shapefile
    output_gdf.to_file(canal_line_output)        
            
def return_waterways(output_file_name):
    '''
    > A function to return waterways stored in a postGIS table 
    > For effiency, we filter waterways that intersect with the CRT canal network (line)
    > The output is saved as a shapefile
    '''
    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Format query
    query = text("select water.*, canal_lines.layer, canal_lines.draught from water, canal_lines where ST_Intersects(water.geometry,canal_lines.geometry);") # Actual geometry intersects

    # Extracts canal geometries, stores as geodataframe
    canals_gdf = GeoDataFrame.from_postgis(sql = query, 
        con = db_connection,
        geom_col='geometry', # Column name for the geometry
        coerce_float=True)
    
    # Sort by date, remove duplicates based on feature id, keep most recent update
    filtered_canals_gdf = canals_gdf.sort_values('versiondate').drop_duplicates(subset=['fid'], keep='last')

    # Simplify column names (< 10 characters), avoids SetWithCopyWarning
    filtered_canals_gdf = filtered_canals_gdf.rename(columns={'descriptiveterm' : 'term', 
                                        'descriptivegroup' : 'group',
                                        'versiondate' : 'datetime'})
    
    # Convert DateTime to formatted string (vectorised)
    filtered_canals_gdf['datetime'] = filtered_canals_gdf['datetime'].dt.strftime('%Y-%m-%d')
    
    # Save to file
    filtered_canals_gdf.to_file(output_file_name, engine = "pyogrio")


# Return mean value, ignoring 0 
def null_mean(x):
    return x[x!=0].mean()

def add_building_heights(input_polygons, input_raster, output_polygons):
    '''
    > Function to add missing building height values to OSM building polygons...
    > Mean values (zonal statistics) from Keany et al. (2022): https://doi.org/10.5194/asr-19-13-2022
    > An inelegant (slow) solution, but not worth optimising
    '''
    
    # Read directly using pyogrio
    buildings_gdf = read_dataframe(input_polygons)

    # Open using rasterio
    with rast_open(input_raster) as src:

        # Extract affine matrix and band 1
        affine = src.transform
        array = src.read(1)

        # Zonal statistics
        df_zonal_stats = DataFrame(zonal_stats(buildings_gdf, array, affine=affine, stats = ['mean'], 
                                               add_stats={'null_mean':null_mean},
                                               all_touched=True))
        
    # Adding statistics back to original GeoDataFrame
    buildings_heights = concat([buildings_gdf, df_zonal_stats], axis=1) 

    # Write to file
    buildings_heights.to_file(output_polygons, engine = "pyogrio")


def missing_building_heights(input_polygons, output_polygons):
    '''
    > add_building_heights() does not produce heights for all buildings (~10% missing). 
    > This function uses the height of the nearest building
    > All values are also converted to the nearest int
    '''

    # Read directly using pyogrio
    buildings_gdf = read_dataframe(input_polygons)

    # Add an ID column
    buildings_gdf["id"] = buildings_gdf.index + 1

    # Remove unncessary columns, and convert str to float, coercing errors to NaN
    buildings_gdf = buildings_gdf[['id','null_mean', 'geometry']]
    buildings_gdf['null_mean'] = to_numeric(buildings_gdf['null_mean'], errors = 'coerce')

    # Init rtree Index
    idx = Index()

    # Iterate through the buildings gdf and load to spatial index
    for id, building in buildings_gdf.iterrows():
        idx.insert(id, building.geometry.bounds)
    
    # Iterate efficiently using itertuples
    for index, feature_tuple in enumerate(buildings_gdf.itertuples()):

        # If the height is missing
        if isnan(feature_tuple.null_mean): 

            # Init bool and 1 nearest neighbours
            valid_height = True
            neighbours = 1

            # Iterate until we have a height value
            while valid_height:

                # Use the spatial index to get the index of the closest n buildings
                nearest_building_index = list(idx.nearest(feature_tuple.geometry.bounds, neighbours))

                # Iterate through list
                for i in nearest_building_index:

                    # Extract corresponding building and height (m)
                    nearest_building = buildings_gdf.iloc[i]
                    height = nearest_building.null_mean

                    # If the nth neighbour is NaN
                    if isnan(height):
                        
                        # Check next index
                        continue

                    else: 

                        # Add height to gdf and update bool
                        buildings_gdf.loc[index, ['null_mean']] = height
                        valid_height = False

                        print(f"Added a building height of {height}")
                        break
            
                # If all the nth neighbours are 0
                if valid_height:      
                        
                    # Look for n+1 neighbour
                    neighbours += 1

        # If there is a height value, check next building
        else:
            continue
            
    # When complete, round column to nearest integer  
    buildings_gdf["mean_int"] = buildings_gdf['null_mean'].round()      
            
    # Overwrite input file
    buildings_gdf.to_file(output_polygons, engine = "pyogrio") 

def generate_unique_id_ireland(final_canal_polygons, canal_id_path, output_file_path):
    '''
    > A function to generate a unique ID for each canal feature in Ireland, using the name of the canal
    > Utilises the Great Britain dict to start
    '''

    # Read directly using pyogrio
    canals_gdf = read_dataframe(final_canal_polygons)

    # Try
    try:
        
        # Load JSON of canal id values
        with open(canal_id_path) as gd_path:
            canal_id_values = load(gd_path)

    # If the file is not found
    except FileNotFoundError:
        print("File not found")
        exit(1)

    else: 

        # Returns the unique canal names
        canal_names = canals_gdf['layer'].unique()

         # Iterates
        for name in canal_names:

            # Generates canal id string
            letters = "".join([s[0] for s in name.split("-") ])

            # If the id string is already in use
            while letters in canal_id_values.values():

                # If the id string contains numbers
                if any(char.isdigit() for char in letters):

                    # Update the final number (+1)
                    letters = letters[:-1] + str(int(letters[-1]) + 1)

                # If this is the first match
                else:

                    # Add one to the end of the string
                    letters = letters + str(1)
                
            # Add key (canal name) and values (canal id string) to dictionary
            canal_id_values[name] = letters
                
        # When complete, save the geodatabase record file (json)
        dump(canal_id_values, open("../outputs/canal-id-values.json", 'w'))  

    # Map the dataframe canal layers to the dict, and then calculative the cumulative count
    canals_gdf['code_id'] = canals_gdf['layer'].map(canal_id_values)
    canals_gdf['cumulative_count'] = canals_gdf.groupby(['code_id']).cumcount() + 1

    # Combine the string and counnt, and remove the latter
    canals_gdf['code_id']  = canals_gdf['code_id'] + "_" + canals_gdf['cumulative_count'].astype(str)
    canals_gdf.drop(columns = 'cumulative_count', inplace = True)

    # Output to file
    canals_gdf.to_file(output_file_path)


def generate_unique_id(final_canal_polygons, canal_id_path):
    '''
    > A function to generate a unique ID for each canal feature, using the name of the canal
    > This is necessary because OS IDs are no longer valid
    '''

    # Read directly using pyogrio
    canals_gdf = read_dataframe(final_canal_polygons)

    # Try
    try:
        
        # Load JSON of canal id values
        with open(canal_id_path) as gd_path:
            canal_id_values = load(gd_path)

    # If the file is not found
    except FileNotFoundError:

        # Returns the unique canal names
        canal_names = canals_gdf['layer'].unique()

        # Initialises empty dictionary
        canal_id_values = {}

        # Iterates
        for name in canal_names:

            # Generates canal id string
            letters = "".join([s[0] for s in name.split("-") ])

            # If the id string is already in use
            while letters in canal_id_values.values():

                # If the id string contains numbers
                if any(char.isdigit() for char in letters):

                    # Update the final number (+1)
                    letters = letters[:-1] + str(int(letters[-1]) + 1)

                # If this is the first match
                else:

                    # Add one to the end of the string
                    letters = letters + str(1)
                
            # Add key (canal name) and values (canal id string) to dictionary
            canal_id_values[name] = letters
                
        # When complete, save the geodatabase record file (json)
        dump(canal_id_values, open("../outputs/canal-id-values.json", 'w'))  
         
    # Map the dataframe canal layers to the dict, and then calculative the cumulative count
    canals_gdf['code_id'] = canals_gdf['layer'].map(canal_id_values)
    canals_gdf['cumulative_count'] = canals_gdf.groupby(['code_id']).cumcount() + 1

    # Combine the string and counnt, and remove the latter
    canals_gdf['code_id']  = canals_gdf['code_id'] + "_" + canals_gdf['cumulative_count'].astype(str)
    canals_gdf.drop(columns = 'cumulative_count', inplace = True)

    # Output to file
    canals_gdf.to_file("../data/canal-geometries-filtered-modified-dissolved-id.shp")

def climate_parallel(input_file_path, start_date, n_cores, days, key):
    '''
    > Function to pre_compute_climate() using multi-processing
    '''

    # Read the dataframe directly using pyogrio
    climate_gdf = read_dataframe(input_file_path)

    # To exclude the problematic ID. NOTE: this is no longer necessary, as #92 has been replaced (09-06-23)
    # options = [92]
    # climate_gdf = climate_gdf.loc[~climate_gdf['CLUSTER_ID'].isin(options)]

    # Split into n clusters, each will be assigned a CPU
    climate_gdf = array_split(climate_gdf, n_cores)

    # Initialises the Pool for multiprocessing (n cores)
    p = Pool(processes = n_cores)

    # Define arguments, some fixed, some variable
    # Climate locations, start date, API key, days of analysis
    args = [(i, start_date, key, days) for i in climate_gdf]

    # Generate results
    p.starmap(precompute_climate_parallel, args)

def precompute_climate_parallel(climate_locations, start_date, key, days):
    '''
    > A function to iterate through climate locations and for each datetime, obtain climate data from OpenWeatherMap
    > Stored in JSON files for efficiency O(1): https://stackoverflow.com/questions/513882/python-list-vs-dict-for-look-up-table
    > Some unexplained errors with API: occasionally returns no data - these are caught with try-excepts and interpolated
    > Non-parallel version took ~16 mins to return climate data for 24 h for all location (n=100) 
    > Parallel version (5 cores) takes ~3-4 mins to acheive the same
    > When complete, dump everything to postGIS
    '''

    # Hours of analysis
    duration = days * 24

    # Iterate efficiently using itertuples
    for feature_tuple in climate_locations.itertuples():

        # Initialise dictionary to store the missing/interpolated datetimes
        missing_data = {}
        interpolated_data = {}

        # Try
        try:
            
            # Load climate record for current location
            with open(f"../data/climate-data/climate-{feature_tuple.CLUSTER_ID}.json") as json_path:
                climate_record = load(json_path)

        # If the file is not found
        except FileNotFoundError:

            # Create a new dictionary
            climate_record = {}

        # Localise niave datetime to aware datetime (UK) and then convert to UTC
        local_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
        
        # Sets end datetime (+ duration)
        end_datetime = local_datetime + timedelta(hours = duration)

        # Convert to series, project to 4326 and return lon-lat (perhaps could be simplified)
        geographic_geometry = GeoSeries(feature_tuple.geometry, crs = 27700).to_crs(4326)
        longitude, latitude = geographic_geometry.iloc[0].x, geographic_geometry.iloc[0].y
        
        # Continue until we have extracted the required duration
        while local_datetime < end_datetime:

            # Convert current time to UTC
            unix_time = str(local_datetime.timestamp())
    
            # If that datetime (UTC) exists in the dictionary
            if unix_time in climate_record:

                # If that data has been obtained directly from OpenWeatherMap ('measured')
                if ('method', 'measured') in climate_record[unix_time].items():
                    
                    # Add one hour to the datetime
                    local_datetime += timedelta(hours=1)

                # If that data has been estimated based on earlier/later times ('interpolated'), try to obtain the data again
                elif ('method', 'interpolated') in climate_record[unix_time].items():
                    
                    # Sets interval end datatime (current +1 hr)
                    interval_datetime = local_datetime + timedelta(hours=1)

                    # Converts dates to Unix timestamp
                    start_time = local_datetime.timestamp()
                    end_time = interval_datetime.timestamp()

                    # API call to Open Weather (formatting for readability)
                    response = requests.get(f"http://history.openweathermap.org/data/2.5/history/city?&" \
                                            f"lat={latitude}&" \
                                            f"lon={longitude}&" \
                                            f"type=hour&" \
                                            f"start={start_time}&" \
                                            f"end={end_time}&" \
                                            f"appid={key}")

                    # Try to load
                    try:

                        # Loads as json
                        weather = loads(response.text)

                    # Very occasionally there are errors loading from the API
                    except decoder.JSONDecodeError:

                        # Error statement
                        # print(f"API loading error at {local_datetime} for ID {feature_tuple.CLUSTER_ID}")

                        # Add to interpolated dictionary
                        interpolated_data[unix_time] = climate_record[unix_time]

                        # Add one hour to the datetime
                        local_datetime += timedelta(hours=1)

                    # The API returned a valid JSON format
                    else:
                    
                        # Try
                        try:

                            # Add to dictionary
                            climate_record[unix_time] = {'air_temperature' : weather['list'][0]['main']['temp'], 
                                                                'humidity' : weather['list'][0]['main']['humidity'], 
                                                                'cloud_cover' : weather['list'][0]['clouds']['all'], 
                                                                'wind_speed' : weather['list'][0]['wind']['speed'],
                                                                'air_pressure' : weather['list'][0]['main']['pressure'],
                                                                'method' : 'measured',
                                                                'future_hours' : 0}

                        # No data is returned (but the JSON is valid)   
                        except KeyError:

                            # Add to interpolated dictionary
                            interpolated_data[unix_time] = climate_record[unix_time]
                                                            
                            # Add one hour to the datetime
                            local_datetime += timedelta(hours=1)

                        # Data was returned and successfully loaded to the dict
                        else:

                            # Add one hour to the datetime
                            local_datetime += timedelta(hours=1)
                    
                # If that datetime has been checked previously and there are no data 
                elif ('method', 'no-data') in climate_record[unix_time].items():

                    # Add this information to the missing data dict{}
                    missing_data[unix_time] = {'method' : 'no-data'} 

                    # Add one hour to the datetime
                    local_datetime += timedelta(hours=1)

            # This is a new datetime
            else:

                # Sets interval end datatime (current +1 hr)
                interval_datetime = local_datetime + timedelta(hours=1)

                # Converts dates to Unix timestamp
                start_time = local_datetime.timestamp()
                end_time = interval_datetime.timestamp()
        
                # API call to Open Weather (formatting for readability)
                response = requests.get(f"http://history.openweathermap.org/data/2.5/history/city?&" \
                                        f"lat={latitude}&" \
                                        f"lon={longitude}&" \
                                        f"type=hour&" \
                                        f"start={start_time}&" \
                                        f"end={end_time}&" \
                                        f"appid={key}")
                # Try to load
                try:
                
                    # Loads as json
                    weather = loads(response.text)

                # Very occasionally there are errors loading from the API (invalid JSON)
                except decoder.JSONDecodeError:

                    # Error statement
                    # print(f"API loading error at {local_datetime} for ID {feature_tuple.CLUSTER_ID}")

                    # Add this information to the missing data dict{}
                    missing_data[unix_time] = {'method' : 'no-data'}
                    climate_record[unix_time] = {'method' : 'no-data',
                                                 'future_hours' : 0}  

                    # Add one hour to the datetime
                    local_datetime += timedelta(hours=1)

                # The API loaded a valid JSON
                else:

                    # Try
                    try:

                        # Add to dictionary
                        climate_record[unix_time] = {'air_temperature' : weather['list'][0]['main']['temp'], 
                                                            'humidity' : weather['list'][0]['main']['humidity'], 
                                                            'cloud_cover' : weather['list'][0]['clouds']['all'], 
                                                            'wind_speed' : weather['list'][0]['wind']['speed'],
                                                            'air_pressure' : weather['list'][0]['main']['pressure'],
                                                            'method' : 'measured',
                                                            'future_hours' : 0}

                    # No data is returned    
                    except KeyError:

                        # Add to dictionary 'no-data'
                        missing_data[unix_time] = {'method' : 'no-data'} 
                        climate_record[unix_time] = {'method' : 'no-data',
                                                    'future_hours' : 0} 

                        # Add one hour to the output (skip for now)
                        local_datetime += timedelta(hours=1)

                    # The data was correctly added to the dictionary
                    else:

                        # Add one hour to the output
                        local_datetime += timedelta(hours=1)

        ''' When we have iterated through once, collecting the relevant data, we then iterate through the missings and interpolate '''
        
        # Iterate through the missing data
        for missing_datetime in missing_data.keys():

            # Return unix time of current time - 1h (3600s)
            pre_unix_time = str(float(missing_datetime) - 3600)

            # Try
            try:

                # Obtain prior weather, using a copy to avoid modifying the main dict
                prior_weather = climate_record[pre_unix_time].copy()

            # If the first record does not exist
            except KeyError:
                
                # Skip to the next datetime
                continue
            
            # There is a prior value
            else:

                # Remove 'method' key:value
                del prior_weather['method']

                # Try to delete the future key
                try:
                    del prior_weather ['future_hours']

                # If there's an key error, the previous cell is missing data
                except KeyError:
                    
                    # Skip to the next datetime
                    continue
                
                # The future key was successfully deleted
                else:

                    # Initialises objects
                    data_availability = False
                    forward_hours = 1

                    # Continue, until we obtain data
                    while data_availability is False: 

                        # Return unix time of current time + forward hours
                        post_unix_time = str(float(missing_datetime) + (forward_hours * 3600))

                        # Try
                        try:

                            # If there are no measured data at this time step
                            if ('method', 'measured') not in climate_record[post_unix_time].items():

                                # Look forward another hour
                                forward_hours += 1

                            # There are data at this time step
                            else:

                                # Update object, to curtail while statement
                                data_availability = True

                                # Obtain the future weather climate record
                                future_weather = climate_record[post_unix_time].copy()

                                # Remove 'method' key:value
                                del future_weather['method']
                                del future_weather['future_hours']

                                # Determine weights for past/future climate (based on proximity)
                                # The weighting of prior weather increases as the h interval increase (1/2, 2/3, 3/4, 4/5)
                                weights = [forward_hours / (forward_hours + 1), 1 / (forward_hours + 1)]

                                # Convert each dict to dataframe
                                prior_df, future_df  = DataFrame([prior_weather]), DataFrame([future_weather])

                                # Weighted mean of past and future climate (prior_weather, future_weather)
                                weighted_weather = (prior_df * weights[0] + future_df * weights[1])
                            
                                # Round and convert columns to integer
                                integer_cols = ['humidity', 'cloud_cover', 'air_pressure']
                                weighted_weather[integer_cols] = weighted_weather[integer_cols].round(0).astype('int')

                                # Round wind speed and temp to two decimal places
                                decimal_cols = ['air_temperature', 'wind_speed']
                                weighted_weather[decimal_cols] = weighted_weather[decimal_cols].round(2)

                                # Add weighted values to dictionary, plus method/duration for missing datetime
                                climate_record[missing_datetime] = {'air_temperature' : weighted_weather.iloc[0]['air_temperature'], 
                                                            'humidity' :  weighted_weather.iloc[0]['humidity'], 
                                                            'cloud_cover' :  weighted_weather.iloc[0]['cloud_cover'], 
                                                            'wind_speed' :  weighted_weather.iloc[0]['wind_speed'],
                                                            'air_pressure' :  weighted_weather.iloc[0]['air_pressure'],
                                                            'method' : 'interpolated', 
                                                            'future_hours' : forward_hours}
                                
                        # We have reached the end of the table (key does not exist)
                        except KeyError:

                            # Move to the next feature
                            break
        
        ''' Although this seems inefficient (looping through again), it actually ensures the interpolation takes place more often'''

        # Iterate through the interpolated data
        for interpolated_datetime in interpolated_data.keys():

            # Initialise objects
            data_availability = False
            forward_hours = 1

            # If possible, we can improve our interpolation by using data closer in time
            while data_availability is False: 

                # Return unix time of current time + forward hours
                post_unix_time = str(float(interpolated_datetime) + (forward_hours * 3600))

                # Try
                try:

                    # If there are no measured data at this time step
                    if ('method', 'measured') not in climate_record[post_unix_time].items():

                        # Look forward another hour
                        forward_hours += 1

                    # There are data at this time step
                    else:

                        # Update object, to curtail while statement
                        data_availability = True

                        # If we have found reached the same future datetime
                        if forward_hours >= climate_record[interpolated_datetime]['future_hours']:
                            
                            # Skip to the next feature
                            continue

                        # We have found 'actual' data, closer to the timestep
                        else:

                            # Obtain the future weather climate record
                            future_weather = climate_record[post_unix_time].copy()

                            # Remove 'method' key:value
                            del future_weather['method']
                            del future_weather['future_hours']

                            # Return unix time of current time - 1h (3600s)
                            pre_unix_time = str(float(interpolated_datetime) - 3600)

                            # Obtain prior weather, using a copy to avoid modifying the main dict
                            prior_weather = climate_record[pre_unix_time].copy()

                            # Remove 'method' key:value
                            del prior_weather['method']
                            del prior_weather['future_hours']

                            # Determine weights for past/future climate (based on proximity)
                            # The weighting of prior weather increases as the h interval increase (1/2, 2/3, 3/4, 4/5)
                            weights = [forward_hours / (forward_hours + 1), 1 / (forward_hours + 1)]

                            # Convert each dict to dataframe
                            prior_df, future_df  = DataFrame([prior_weather]), DataFrame([future_weather])

                            # Weighted mean of past and future climate (prior_weather, future_weather)
                            weighted_weather = (prior_df * weights[0] + future_df * weights[1])
                        
                            # Round and convert columns to integer
                            integer_cols = ['humidity', 'cloud_cover', 'air_pressure']
                            weighted_weather[integer_cols] = weighted_weather[integer_cols].round(0).astype('int')

                            # Round wind speed and temp to two decimal places
                            decimal_cols = ['air_temperature', 'wind_speed']
                            weighted_weather[decimal_cols] = weighted_weather[decimal_cols].round(2)

                            # Add weighted values to dictionary, plus method/duration for missing datetime
                            climate_record[interpolated_datetime] = {'air_temperature' : weighted_weather.iloc[0]['air_temperature'], 
                                                        'humidity' :  weighted_weather.iloc[0]['humidity'], 
                                                        'cloud_cover' :  weighted_weather.iloc[0]['cloud_cover'], 
                                                        'wind_speed' :  weighted_weather.iloc[0]['wind_speed'],
                                                        'air_pressure' :  weighted_weather.iloc[0]['air_pressure'],
                                                        'method' : 'interpolated', 
                                                        'future_hours' : forward_hours}
                            

                # KeyError (reached the end of the dictionary)
                except KeyError:

                    continue

        
        # Proportion of measured values
        measured_count = sum(value['method'] == 'measured' for key, value in climate_record.items())
        total_count = len(climate_record)          
        measured_proportion = measured_count / total_count  * 100

        # Summary statement
        print(f"For ID {feature_tuple.CLUSTER_ID} and for {total_count/24} days of data, {round(measured_proportion, 1)}% of values are measured ({int(measured_count)}/{int(total_count)})")

        # When complete, save the climate data file (json)
        dump(climate_record, open(f"../data/climate-data/climate-{feature_tuple.CLUSTER_ID}.json", 'w'))

def preprocess_climate_data(input_path, output_path):
    '''
    > A function to process the 40-year climate data from OpenWeatherMap (downloaded in bulk, stored in JSON)
    > Extract relevant metrics (temp, humidity...) and store in a new JSON
    '''
    
    # File paths for all climate JSONs
    climate_files = g(input_path, recursive = True)

    # Iterate through
    for file in climate_files:

        # Create empty dictionary
        climate_dict = {}

        # Open the file
        with open(file) as json_path:
            climate_record = load(json_path)

        # Extract the ID
        cluster_id = climate_record[0]['city_name'][8:]
    
        # If the output file already exists
        if isfile(f"{output_path}climate-{cluster_id}.json"):
            print(f"{output_path}climate-{cluster_id}.json exists")
            continue
        
        else:
        
            # Iterate through each datetime
            for item in climate_record:
                
                # Add to output dictionary
                climate_dict[item['dt']] = {'air_temperature' : item['main']['temp'], 
                                            'humidity' : item['main']['humidity'], 
                                            'cloud_cover' : item['clouds']['all'], 
                                            'wind_speed' : item['wind']['speed'],
                                            'air_pressure' : item['main']['pressure']}

            # When complete, save the climate data file (json)
            dump(climate_dict, open(f"{output_path}climate-{cluster_id}.json", 'w'))


def filter_climate_data(input_path, output_path, start_date, end_date):
    '''
    > A function to process the 40-year climate data from OpenWeatherMap (downloaded in bulk, stored in JSON)
    > Extract relevant metrics (temp, humidity...) and store in a new JSON
    '''
    
    # Localise niave datetime to aware datetime (UK) and then convert to UTC
    start_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
        
    # Sets end datetime 
    end_datetime = to_datetime(end_date)\
        .tz_localize('Europe/London')\
        .tz_convert('UTC')

    # Convert timed to UTC, 
    start_unix, end_unix = start_datetime.timestamp(), end_datetime.timestamp()

    # File paths for all climate JSONs
    climate_files = g(input_path, recursive = True)

    # Iterate through
    for file in climate_files:

        # Get the file name from the full path
        file_name = split(file)[1]

        # If the file already exists
        if isfile(f"{output_path}{file_name}"):
            continue
        
        # Else
        else:

            # Open the file
            with open(file) as json_path:
                climate_record = load(json_path)

                # Filter to defined period only (2022)
                filtered_record = {key:val for key, val in climate_record.items() if float(key) >= start_unix and float(key) < end_unix}

                # When complete, save the climate data file (json), using the same name as input
                dump(filtered_record, open(f"{output_path}{file_name}", 'w'))
            

def calc_climate_weighting(input_climate_path, input_canal_path, proximity_path, dist, write):

    # Read the dataframe directly using pyogrio
    climate_gdf = read_dataframe(input_climate_path)

    # Create the engine
    db_connection = create_engine(DB_STRING)
    
    # Add climate locations to postGIS
    if write:

        # For Britain
        if climate_gdf.crs == 27700:

            # Add to postgis 
            climate_gdf.to_postgis('climate_locations', db_connection, if_exists='fail', dtype={'geom': Geometry('[Point]', srid=27700)})

        # For Ireland
        elif climate_gdf.crs == 29902:

            # Add to postgis 
            climate_gdf.to_postgis('climate_locations_ireland', db_connection, if_exists='fail', dtype={'geom': Geometry('[Point]', srid=29902)})

    # Open the canals data using pyogrio
    canals_gdf = read_dataframe(input_canal_path)

    # Create a dictionary for storing proximity
    proximity_dict = {}

    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # For Britain
        if climate_gdf.crs == 27700:

            # Format postGIS query - extract climate locations within a set distance
            query = text("select climate_locations.* from climate_locations, canal_feature where ST_DWithin(climate_locations.geometry, canal_feature.geometry, :d) and canal_feature.code_id = :f")  \
            .bindparams(d = dist, # Threshold distance
                        f = feature_tuple.code_id) # Feature ID

        # For Ireland
        elif climate_gdf.crs == 29902:  

            # Format postGIS query - extract climate locations within a set distance
            query = text("select climate_locations_ireland.* from climate_locations_ireland, canals_ireland where ST_DWithin(climate_locations_ireland.geometry, canals_ireland.geometry, :d) and canals_ireland.code_id = :f")  \
            .bindparams(d = dist, # Threshold distance
                        f = feature_tuple.code_id) # Feature ID

        # Extracts geometries using from_postgis, stores as geodataframe
        climate_locations = GeoDataFrame.from_postgis(sql = query, 
        con = db_connection,
        geom_col='geometry', # Column name for the geometry
        coerce_float=True)

        # Init lists for storing values
        id_vals = []
        distances = []

        # Iterate through the IDs and calculate distances
        for location in climate_locations.itertuples():

            # Append values to list, distance using Shapely
            id_vals.append(location.CLUSTER_ID)
            distances.append(feature_tuple.geometry.distance(location.geometry))

        try:
            # Determine weighting
            # > https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Inverse-Distance-Weighting/index.html
            inverse = [d**-1 for d in distances]

        # If the distance is 0 (i.e. the climate location is within the feature geometry)
        except ZeroDivisionError:

            # Only use that climate location
            proximity_dict[feature_tuple.code_id] = {'id' : [id_vals[0]],
                                                 'd' : [distances[0]], 
                                                 'w' : [1.0]}

        # If the distances are valid (>0)
        else:

            # > https://stackoverflow.com/questions/22819502/how-to-normalize-data-to-1-such-that-each-value-receives-a-weight-proportional-t
            # Weight, round to three decimal places
            weighting = [round(i * (1 / sum(inverse)), 3) for i in inverse]

            # When complete, add to dictionary
            # Key is the canal feature code, values are ids, distances and weights
            proximity_dict[feature_tuple.code_id] = {'id' : id_vals,
                                                    'd' : distances, 
                                                    'w' : weighting}


    # When complete, save the climate data file (json), using the same name as input
    dump(proximity_dict, open(proximity_path, 'w'))

def extract_climate_data(input_canal_path, proximity_path, climate_path, output_path):
    '''
    > Function to return climate data for each canal geometry, based on weighting of climate location distance
    > Completed 19/06 - took 8414 seconds
    > An inelegent and inefficient (but functional) approach
    '''
    
    # Open the canals data using pyogrio
    canals_gdf = read_dataframe(input_canal_path)

    # Open the proximity record
    with open(proximity_path) as json_path:
                proximity_record = load(json_path)


    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # If the output file already exists
        if isfile(f"{output_path}data-{feature_tuple.code_id}.json"):
            print(f"{output_path}data-{feature_tuple.code_id}.json exists")
            continue

        else:

            # Return the weighting values and ids
            proximity = proximity_record[feature_tuple.code_id]

            # Init dict for storing output
            output_dict = {}

            # If the climate weighting is based on one location only
            if len(proximity['id']) == 1:

                # Open the corresponding climate record
                with open(f"{climate_path}climate-{proximity['id'][0]}.json") as c_path:
                    climate_data = load(c_path)

                # Output to the new directory
                dump(climate_data, open(f"{output_path}data-{feature_tuple.code_id}.json", 'w'))             
            
            # If spatial weighting is required (n > 1)
            else:

                # Iterate through the nested id values and weightings
                for id, wgt in zip(proximity['id'], proximity['w']):

                    # Open the corresponding climate record
                    with open(f"{climate_path}climate-{id}.json") as c_path:
                        climate_data = load(c_path)
        
                    # Iterate through each item
                    for key, val in climate_data.items():    
                        
                        # Update the nested dictionary values
                        climate_data[key] = {'air_temperature' : val['air_temperature'] * wgt, 
                                                    'humidity' : val['humidity'] * wgt, 
                                                    'cloud_cover' : val['cloud_cover'] * wgt, 
                                                    'wind_speed' : val['wind_speed'] * wgt,
                                                    'air_pressure' : val['air_pressure'] * wgt}

                    # If the dictionary contains values
                    if output_dict:

                        # Nested for loop through both dictionaries
                        for key1, val1 in output_dict.items():
                            for key2, val2 in climate_data.items():

                                # If the keys match
                                if key1 == key2:

                                    # Update the output dictionary
                                    output_dict[key1] = {'air_temperature' : round(val1['air_temperature'] + val2['air_temperature'], 2), # 2 dp
                                                    'humidity' : round(val1['humidity'] + val2['humidity'], 0), # Integer
                                                    'cloud_cover' : round(val1['cloud_cover'] + val2['cloud_cover'], 0), # Integer
                                                    'wind_speed' : round(val1['wind_speed'] + val2['wind_speed'], 2), # 2 dp
                                                    'air_pressure' : round(val1['air_pressure'] + val2['air_pressure'], 0)} # Integer
                                    
                                    # Break the inner loop, skip to next outer key
                                    break
                                
                    # The dictionary is empty
                    else:

                        # Fill in the weighted climate data
                        output_dict = climate_data.copy()

                # When complete, save the climate data file (json), using the same name as input
                dump(output_dict, open(f"{output_path}data-{feature_tuple.code_id}.json", 'w'))       

def interpolate_climate_data(input_path, output_path):
    '''
    Function to interpolate 15-m intervals for hourly climate data
    '''

    # File paths for all climate JSONs
    climate_files = g(input_path, recursive = True)

    # Iterate
    for file in climate_files:
        
        # Get the file name from the full path
        file_name = split(file)[1]

        # If the output file already exists
        if isfile(f"{output_path}{file_name}"):

            continue
        
        # Otherwise, interpolate
        else:

            # Init empty dictionary for storing interpolated values
            interpolated_dict = {}

            # Open the corresponding climate record
            with open(file) as c_path:
                climate_data = load(c_path)

            # Iterate
            for key, value in climate_data.items():
                
                # Try
                try:
                
                    # Return the values for the following hour (3600s per hour)
                    post_hour = climate_data[str(int(key) + 3600)]

                # If we've reached the end of the dictionary
                except KeyError:
                    
                    # Exit the current loop
                    break

                # Difference between the values
                difference = {'air_temperature' : post_hour['air_temperature'] - value['air_temperature'], 
                            'humidity' : post_hour['humidity'] - value['humidity'], 
                            'cloud_cover' : post_hour['cloud_cover'] - value['cloud_cover'], 
                            'wind_speed' : post_hour['wind_speed'] - value['wind_speed'],
                            'air_pressure' : post_hour['air_pressure'] - value['air_pressure']}
            
                # Iterate through times (s) and multipliers
                for time, multiplier in zip([900, 1800, 2700], [0.25, 0.50, 0.75]):
                    
                    # Add values to interpolated dictionary, all values rounded to two decimal places
                    interpolated_dict[str(int(key) + time)] = {'air_temperature' : round(value['air_temperature'] + (difference['air_temperature'] * multiplier), 2),
                                                                'humidity' : round(value['humidity'] + (difference['humidity'] * multiplier), 2), 
                                                                'cloud_cover' : round(value['cloud_cover'] + (difference['cloud_cover'] * multiplier), 2),
                                                                'wind_speed' : round(value['wind_speed'] + (difference['wind_speed'] * multiplier), 2), 
                                                                'air_pressure' : round(value['air_pressure'] + (difference['air_pressure'] * multiplier), 2)}

            # Merge the two dictionaries (hourly data | interval data)
            output_dict = climate_data | interpolated_dict

            # When complete, save the climate data file (json) in the new directory, but with the same name
            dump(output_dict, open(f"{output_path}{file_name}", 'w'))        

def spin_up_climate(start_date, end_date, output_path):
    
    '''
    > First filters climate to the desired interval
    > Then interpolates (1h > 15m)
    > Next iterates through climate location (n=100) and returns average values for each climate parameter
    > This is then used as model spin up (to generate starting canal temperatures and energy values)
    '''

    # File paths for all climate JSONs
    climate_files = g(r'../../national-canal-cooling-data/spin-up-climate-interpolated/*.json', recursive = True)

    # Localise niave datetime to aware datetime (UK) and then convert to UTC
    start_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
        
    # Sets end datetime 
    end_datetime = to_datetime(end_date)\
        .tz_localize('Europe/London')\
        .tz_convert('UTC')

    # Convert timed to UTC, 
    start_unix, end_unix = start_datetime.timestamp(), end_datetime.timestamp()

    # Output dict
    output = {}

    # Iterate, beginning at 1
    for count, file in enumerate(climate_files, 1):
        
        # Open the corresponding climate record
        with open(file) as c_path:
            climate_data = load(c_path)

        # Filter to defined period only (1st Jan. 2022)
        filtered_record = {key:val for key, val in climate_data.items() if float(key) >= start_unix and float(key) < end_unix}

        # Init empty lists for storing intermediate and output values
        temps, humid, clouds, wind, pressure = [], [], [], [], []
        output_temp, output_humid, output_clouds, output_wind, output_pressure = [], [], [], [], []

        # Iterate through dict
        for k,v in filtered_record.items():
            # Add values to empty lists
            temps.append(v['air_temperature'])
            humid.append(v['humidity'])
            clouds.append(v['cloud_cover'])
            wind.append(v['wind_speed'])
            pressure.append(v['air_pressure'])

        # If the output list contains data    
        if output_temp:

            # Output is weighted average of inputs
            output_temp = (output_temp * ((1 / count) * (count - 1))) + (temps * (1 / count))

        # The lists are empty
        else: 
            
            # Assign the first values to the output
            output_temp, output_humid, output_clouds, output_wind, output_pressure = temps, humid, clouds, wind, pressure

        print(f"Completed {count} weightings")

    # Now iterate through the keys
    for index, k in enumerate(filtered_record.keys()):

        # Assign to output
        output[k] = {'air_temperature' : output_temp[index], 
                            'humidity' : output_humid[index], 
                            'cloud_cover' : output_clouds[index], 
                            'wind_speed' : output_wind[index],
                            'air_pressure' : output_pressure[index]}
        
    # Init dict for storing 15-m values
    interpolated_dict = {}
        
    # Iterate
    for key, value in output.items():
        
        # Try
        try:
        
            # Return the values for the following hour (3600s per hour)
            post_hour = output[str(int(key) + 3600)]

        # If we've reached the end of the dictionary
        except KeyError:
            
            # Exit the current loop
            break

        # Difference between the values
        difference = {'air_temperature' : post_hour['air_temperature'] - value['air_temperature'], 
                    'humidity' : post_hour['humidity'] - value['humidity'], 
                    'cloud_cover' : post_hour['cloud_cover'] - value['cloud_cover'], 
                    'wind_speed' : post_hour['wind_speed'] - value['wind_speed'],
                    'air_pressure' : post_hour['air_pressure'] - value['air_pressure']}
    
        # Iterate through times (s) and multipliers
        for time, multiplier in zip([900, 1800, 2700], [0.25, 0.50, 0.75]):
            
            # Add values to interpolated dictionary, all values rounded to two decimal places
            interpolated_dict[str(int(key) + time)] = {'air_temperature' : round(value['air_temperature'] + (difference['air_temperature'] * multiplier), 2),
                                                        'humidity' : round(value['humidity'] + (difference['humidity'] * multiplier), 2), 
                                                        'cloud_cover' : round(value['cloud_cover'] + (difference['cloud_cover'] * multiplier), 2),
                                                        'wind_speed' : round(value['wind_speed'] + (difference['wind_speed'] * multiplier), 2), 
                                                        'air_pressure' : round(value['air_pressure'] + (difference['air_pressure'] * multiplier), 2)}

    # Merge the two dictionaries (hourly data | interval data)
    output_dict = output | interpolated_dict
            
    # When complete, save the climate data file (json) in the new directory, but with the same name
    dump(output_dict, open(f"{output_path}", 'w')) 

def check_matches(input_file_path, start_date, days, precision):
    '''
    > Function to check how often solar azimuth-altitude pairs are repeated across a year
    > At two decimal places, this never happens
    > At one decimal places, this happens irregularly
    > At zero decimal places, this happens much more frequently
    > Using 1 d.p. is probably a reasonable approximation 
    '''


    # Try
    try:
        
        # Load shading record for current location
        with open(f"../data/shading-data/shading-satc_1-annual-record.json") as json_path:
            match_record = load(json_path)

    # If the file is not found
    except FileNotFoundError:

        # Read the dataframe directly using pyogrio
        canals_gdf = read_dataframe(input_file_path)

        # For model testing, we'll work on the first one alone
        canals_gdf = canals_gdf.head(1)

        # Create a dictionary to store the results
        match_record = {}

        # Hours of analysis
        duration = days * 24

        # Iterate efficiently using itertuples
        for feature_tuple in canals_gdf.itertuples():

            # Convert to series, project to 4326 and return lon-lat (perhaps could be simplified)
            geographic_geometry = GeoSeries(feature_tuple.geometry.centroid, crs = 27700).to_crs(4326)
            longitude, latitude = geographic_geometry.iloc[0].x, geographic_geometry.iloc[0].y

            # Localise niave datetime to aware datetime (UK) and then convert to UTC
            local_datetime = to_datetime(start_date)\
                .tz_localize('Europe/London')\
                .tz_convert('UTC')
            
            # Sets end datetime (+ duration)
            end_datetime = local_datetime + timedelta(hours = duration)
            
            # Continue until we have modelled the required duration
            while local_datetime < end_datetime:

                # Convert current time to UTC
                unix_time = str(local_datetime.timestamp())

                # Calculate the solar altitude (in degrees), rounded to specified number of decimal places
                solar_altitude = round(get_altitude(latitude, longitude, local_datetime), precision)

                # Calculate the solar azimuth (in degrees), as above
                solar_azimuth = round(get_azimuth(latitude, longitude, local_datetime), precision)

                # Add key (datetime) and previously calculatde shading proportion to dictionary
                match_record[unix_time] = {'altitude' : solar_altitude,
                                            'azimuth' : solar_azimuth}
                
                # Add 15 minutes to the datetime 
                local_datetime += timedelta(hours=0.25)

        # When complete
        dump(match_record, open(f"../data/shading-data/shading-{feature_tuple.code_id}-annual-record.json", 'w'))  

    # Iterate through the dictionary
    else:
                    
        for key, value in match_record.items():
                    
            #print(value)

            # Use dictionary comprehension to check if we have modelled this altitude-azimuth pair befor
            matched_dict = {k:v for (k,v) in match_record.items() if v['altitude'] == value['altitude'] and v['azimuth'] == value['azimuth']}

            print(f"There are {len(matched_dict)} matches for {key}")

def return_buildings_polygon(feature, distance, plot, epsg):
    '''
    > Returns buildings (stored in postGIS) within a specified {distance} from a canal feature (based on ID)
    '''

    # Create the engine, disable pooling
    # Source: https://docs.sqlalchemy.org/en/20/core/pooling.html#pooling-multiprocessing
    db_connection = create_engine(DB_STRING, poolclass=NullPool)

    # If the input data is for Ireland
    if epsg == 29902:
        query = text("select buildings_ireland.* from buildings_ireland, canals_ireland where ST_DWithin(buildings_ireland.geometry, canals_ireland.geometry, :d) and canals_ireland.code_id = :f")  \
        .bindparams(d = distance, # Threshold distance
                    f = feature.code_id) # Feature ID  

    # For Great Britain
    elif epsg == 27700:
        query = text("select buildings.* from buildings, canal_feature where ST_DWithin(buildings.geometry, canal_feature.geometry, :d) and canal_feature.code_id = :f")  \
        .bindparams(d = distance, # Threshold distance
                    f = feature.code_id) # Feature ID  

    # Extracts geometries using from_postgis, stores as geodataframe
    buildings_via_python = GeoDataFrame.from_postgis(sql = query, 
    con = db_connection,
    geom_col='geometry', # Column name for the geometry
    index_col='os_topo_toid', # Column name for the unique ID
    coerce_float=True)
    
    # Remove date-time field
    buildings_via_python = buildings_via_python.drop(columns='bha_processdate')

    # Plot the output
    if plot:

        # If the image already exists
        if isfile(f'../images/buildings/buildings-around-{feature.code_id}.png'):
            pass

        # If not, plot
        else:

            # Store tuple geometry as series
            canal_object = GeoSeries(Polygon(feature.geometry), crs = buildings_via_python.crs)

            # Set up output image
            fig, my_ax = subplots(1, 1, figsize=(16, 10))
            title(f"OS buildings {feature.code_id} (n = {len(buildings_via_python)}) within {distance} m")

            # Plot the canal geometries
            canal_object.plot(
                ax = my_ax,
                color = "#80B8FF",
                edgecolor = '#3A90FE',
                linewidth = 0.5,
                )

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
            savefig(f'../images/buildings/buildings-around-{feature.code_id}.png', bbox_inches='tight')

            # Close the figure
            close()

    # Return geometries of buildings
    return buildings_via_python

def load_buildings_ireland(building_path, chunk_size):
    '''
    > A function to load buildings for Ireland and Northern Ireland (stored as gdb) to a new postGIS database
    > One shapefile, ~300_000 buildings (~10% of GB database)
    '''

    # Establish postGIS database connection
    db_connection = create_engine(DB_STRING)

    # Using fiona
    with fopen(building_path) as features:

        # Return the number of features
        database_size = len(features)

    # Round up file size to nearest chunk size (10,000)
    maximum_database_size = round_up(database_size, chunk_size)

    # Create list of start-end points
    db_chunks = partition_v2(maximum_database_size, chunk_size)

    # Iterate through the chunks
    for chunk in db_chunks:
        
        # Read using pyogrio
        try:
            building_gdf = gpd_read(building_path, 
                                    rows = slice(chunk[0], chunk[1]), # Slice database based on chunk start-end
                                    engine = "pyogrio") # For efficiency
        
        # ValueError: 'skip_features' must be between 0 and dataset length
        except ValueError:

            # Break this loop (skip to next file)
            break

        # If the gdf is empty
        if building_gdf.empty:

            # Break this loop (skip to next file)
            break

        # If buildings are present
        else: 
            
            # Rename columns
            building_gdf.rename(columns={'mean_int': 'relhmax', 'id': 'os_topo_toid'}, inplace=True)

            # Add date column
            building_gdf['bha_processdate'] = date.today()

            # Try to add to postGIS
            try:

                # Add to postGIS
                building_gdf.to_postgis('buildings_ireland', db_connection, if_exists='append', 
                                                dtype={'geom': Geometry('[Polygon]', srid=29902)})
            
            # If multipolygons have been produced:
            except errors.InvalidParameterValue: 

                # Multi-part geometries into multiple single geometries
                exploded_building_gdf = building_gdf.explode(index_parts = False)

                # Then, add to postGIS
                exploded_building_gdf.to_postgis('buildings_ireland', db_connection, if_exists='append', 
                        dtype={'geom': Geometry('[Polygon]', srid=29902)})
    

    print("Buildings added to postGIS")

def subdivide_climate_data(file_path, ids):
    '''
    > Function to subdivide climate data into the cluster locations
    > Required because I accidently selected the bulk download option (one file), rather than a separate file per location
    '''
    
    # Load climate-data
    with open(file_path) as climate_path:
        climate = load(climate_path)

    # Iterate through the cluster IDs
    for id in ids:

        # Init output list
        output = []

        # If the file already exists
        if isfile(f"../../national-canal-cooling-data/raw-climate/{id}.json"):
            print(f"../../national-canal-cooling-data/raw-climate/{id}.json exists")
            
            # Skip to next ID
            continue
        
        # If not, iterate
        else:
            for item in climate:
                
                # If this matches the desired cluster
                if item['city_name'] == id:

                    # Update the list 
                    output.append(item)

            # When we've tested them all, write to file
            dump(output, open(f"../../national-canal-cooling-data/raw-climate/{id}.json", 'w'))

            # Progress statement
            print(f"Completed ../../national-canal-cooling-data/raw-climate/{id}.json")

def split_dataframe(df, chunk_size = 10000): 
    '''
    > Function to split dataframes into chunks of n size
    > Similiar in functionality to array_split, which is currently returning a warning ('GeoDataFrame.swapaxes' is deprecated and will be removed in a future version. Please use 'GeoDataFrame.transpose' instead.)
    > Source: https://stackoverflow.com/questions/17315737/split-a-large-pandas-dataframe
    '''
    # Init empty list
    chunks = list()

    # Determine the number of chunks
    num_chunks = len(df) // chunk_size + 1

    # Iterate, append and return
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

def calculate_canal_height(input_canal_path, height_dict):

    '''
    > Function to return the elevation of each canal 
    > For simplicity, this is via an API call to Open Topo Data (https://www.opentopodata.org/api/)
    > Data source: EU-DEM (25m)
    > Note: not used in the analysis
    '''

    # Open the canals data using pyogrio
    canals_gdf = read_dataframe(input_canal_path)

    # Init count to store # API calls
    count = 0
    
    # Try
    try:
        
        # Load canal height record for current location
        with open(height_dict) as json_path:
            canal_height_record = load(json_path)

    # If the file is not found
    except FileNotFoundError:

        # Create a new dictionary
        canal_height_record = {}


    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # If that canal ID exists in the dictionary
        if feature_tuple.code_id in canal_height_record:

            # Skip to next feature
            continue
        
        # New feature
        else:

            # Maximum 1,000 API calls per day
            if count < 1000:

                # Return centroid location of canal geometry, convert to lat-long
                centroid = GeoSeries(feature_tuple.geometry.centroid, crs = 27700).to_crs(4326)
                
                # API call to Open Topo Data (https://www.opentopodata.org/api/)
                url = "https://api.opentopodata.org/v1/eudem25m"
                data = {"locations": f"{centroid[0].y}, {centroid[0].x}"}
                response = requests.post(url, json=data)

                # Loads as json
                elevation = loads(response.text)

                # Print API output
                print(elevation)

                # Add to dictionary
                canal_height_record[feature_tuple.code_id] = elevation['results'][0]['elevation']

                # Count the API call
                count += 1

                # Add 1 seconds time delay
                sleep(1)

            # If we've reached the maximum number of API calls
            else:

                # Curtail measurement
                break

    # Dump output to JSON
    dump(canal_height_record, open(height_dict, 'w'))  

    # If all the feature IDs are in the dictionary
    if all(k in canal_height_record for k in (canals_gdf['code_id'])):

        # Map the dict to the df
        canals_gdf['ground'] = canals_gdf['code_id'].map(canal_height_record)

        # Save to shapefile
        canals_gdf.to_file("../data/canal-geometries-filtered-modified-dissolved-id-height.shp")

    else:

        print("Not all features have been tested")
        
def shading_parallel(input_file_path, start_date, n_cores, days, d):
    '''
    > Function to precompute_shading() using multi-processing
    '''

    # Read the dataframe directly using pyogrio
    canals_gdf = read_dataframe(input_file_path)

    # Subset to n options
    # options = ['ac1_1', 'ac1_2', 'ac1_3', 'ac1_4', 'ac1_5']
    # canals_gdf = canals_gdf.loc[canals_gdf['code_id'].isin(options)]

    # Split into n clusters, each will be assigned a CPU, using numpy
    #canals_gdf = array_split(canals_gdf, n_cores)

    # Split into n clusters, each will be assigned a CPU, using function
    canals_gdf = split_dataframe(canals_gdf, int(len(canals_gdf) / n_cores))

    # Initialises the Pool for multiprocessing (n cores)
    p = Pool(processes = n_cores)

    # Define arguments, some fixed, some variable
    # Canal data, start date, days of analysis
    args = [(i, start_date, days, d) for i in canals_gdf]

    # Generate results
    p.starmap(precompute_shading, args)

def normalise_azimuth(az):
    '''
    > get_position()['azimuth'] returns:
    - solar azimuth in radians, measured relative to south (east = negative, west = positive)
    - This is converted to the more conventional direction
    - These results are consistent with https://www.sunearthtools.com/dp/tools/pos_sun.php
    '''

    # Convert azmiuth from relative to south to relative to north (degrees) and return
    if az < 0: 
        az = degrees(radians(180) - abs(az))
    else:
        az = degrees(radians(180) + abs(az))
    return az

def filter_buildings_polygon(building_data, canal_data, canal_bounds, date, lon, lat, distance, approach):

    '''
    > A new function to improve efficiency
    > return_buildings() returns all buildings within a specified distance (m) of a location
    > Building geometries are then used to calculate shading using return_shading() but this could be more efficient
    > However, only a small number of those buildings could shade the studied location given the sun's azimuth (direction)
    > This function returns those buildings based on:
        - date, latitude, longitude = for calculating the solar azimuth 
        - tolerance (°) = threshold (±) for determining whether buildings are in the current solar direction
    > Rather that checking the azmiuth for each building centroid or coordinate sequence...
    > We simply generate the desired area, and then perform an intersects
    > Previously used suncalc, shifted to Pysolar e.g.,
        # Returns the solar azimuth in radians (no rounding)
        precise_solar_azimuth = get_position(date, longitude, latitude)['azimuth']
        # Convert azmiuth from relative to south to relative to north (degrees)
        precise_solar_azimuth = normalise_azimuth(precise_solar_azimuth)
    ''' 
    
    # Use the actual geometry
    if approach == "actual":

        # Using PySolar
        precise_solar_azimuth = get_azimuth(lat, lon, date)

        # Return the polygon exterior coordinates
        lon_coords, lat_coords = canal_data.geometry.iloc[0].exterior.coords.xy

        # Generates sequences of distances and azimuths, length corresponding to the number of coordinates
        distances = [distance] * len(lat_coords)
        azimuths = [precise_solar_azimuth] * len(lat_coords)

        # Translate the bounds, by a set distance and direction (solar azimuth) using pyproj (.fwd)
        translated_coordinates = g_model.fwd(lon_coords, lat_coords, az = azimuths, dist = distances, radians = False)[0:2]

        # Init lists
        interior_coords = []
        exterior_coords = []

        # End lon, end lat, start lon, start lat [0,1,2,3]
        for i in zip(translated_coordinates[0], translated_coordinates[1], lon_coords, lat_coords):
            
            # Draw a line
            modelled_line = LineString([Point(i[0], i[1]), Point(i[2], i[3])])

            # If the line crosses the polygon
            if modelled_line.crosses(canal_data.geometry.iloc[0]):
                
                continue
            
            # The line is valid
            else:

                # Append to list
                interior_coords.append(Point(i[2], i[3]))
                exterior_coords.append(Point(i[0], i[1]))

        
        # When complete, reverse one of the lists and combine
        exterior_coords.reverse()
        coords = interior_coords + exterior_coords

        # Convert to Polygon
        p = GeoSeries(Polygon(coords), crs = 4326).to_crs(building_data.crs)

        # Filter buildings: intersect of building data and polygon, creating a copy
        filtered_buildings = building_data[building_data.intersects(p.geometry.iloc[0])].copy()

        # Return 
        return filtered_buildings
    
    # Simplify, using a convex hull
    elif approach == "simplify":

        # Using PySolar
        precise_solar_azimuth = get_azimuth(lat, lon, date) # Correct

        # Return the polygon exterior coordinates
        lon_coords, lat_coords = canal_bounds.iloc[0].exterior.coords.xy

        # Generates sequences of distances and azimuths, length corresponding to the number of coordinates
        distances = [distance] * len(lat_coords)
        azimuths = [precise_solar_azimuth] * len(lat_coords)

        # Translate the bounds, by a set distance and direction (solar azimuth) using pyproj (.fwd)
        translated_coordinates = g_model.fwd(lon_coords, lat_coords, az = azimuths, dist = distances, radians = False)[0:2]

        # Convert to Polygon object
        p = Polygon(list(zip(translated_coordinates[0], translated_coordinates[1])))

        # Combine original and translated bounds into a single GeoSeries, then perform a unary union
        combined_canal_bounds = GeoSeries([p, canal_bounds.geometry.iloc[0]], crs = 4326).to_crs(27700)
        combined_canal_bounds = combined_canal_bounds.unary_union

        # Calculate convex hull 
        bounds_convex_hull = combined_canal_bounds.convex_hull

        # Filter buildings: intersect of building data and convex hull, creating a copy
        filtered_buildings = building_data[building_data.intersects(bounds_convex_hull)].copy()

        # Return 
        return filtered_buildings

    
def line_of_sight(r0, c0, r1, c1, resolution, elevation_data, output, cs):
    """
    * Runs a single ray-trace from one point to another point, returning a list of visible cells
    """

    # Initialise variablesp
    cur_dydx = 0 		  	# Current dydx
    max_dydx = 0 	  		# Maximum dydx
    distance_travelled = 0  # Euclidean distance

    # Pixel coordinates along the line (excluding the first)
    pixels = column_stack(line(r0, c0, r1, c1))[1:]

    # Convert to tuple, excluding the initial canal cell
    modelled_tuple = tuple(map(tuple, pixels))
                           
    # Check whether any set element is present in modelled line using any() 
    if any(i in cs for i in modelled_tuple):
                
        # Return unmodified output
        return output

    # Iterate
    for r, c in pixels:

        # Distance travelled
        distance_travelled = hypot(c0 - c, r0 - r)

        try:

            # Calculates current dy / dx
            cur_dydx = (elevation_data[(r, c)]) / (distance_travelled * resolution)

        except IndexError:
            
            # Location is out-of-bounds, so set to 0
            cur_dydx = 0

        # Update max dy / dx if current max, and set output to 1
        if (cur_dydx > max_dydx):
            max_dydx = cur_dydx
            output[(r, c)] = 1

    # Return updated output surface
    return output

def viewshed(r0, c0, radius_px, resolution, building_raster, cs):
	"""
	* Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
	"""

	# Create output array at the same dimensions as the building raster
	output = zeros(building_raster.shape)

	# Iterate through the circle perimeter, defined by the radius in pxiels
	for r, c in column_stack(circle_perimeter(r0, c0, radius_px)):

		# Calculate line of sight for each 
		output = line_of_sight(r0, c0, r, c, resolution, building_raster, output, cs)

	# return the resulting viewshed
	return output

def subset_buildings(building_data, canal_data, write, cell_size, radius, epsg):
    '''
    # Filtering buildings based on a viewshed algorithm, adapted from: https://github.com/jonnyhuck/green-visibility-index
    '''

    # Convert tuple to GeoSeries
    projected_building = GeoSeries(building_data.geometry, crs = epsg)
    projected_canal = GeoSeries(canal_data.geometry, crs = epsg)

    # Convert canal polygon to line 
    projected_building = projected_building.exterior
    projected_canal = projected_canal.exterior

    # Reset index values as the OSGB values are very long! 
    building_data.reset_index(drop = False, inplace = True) 

    # List of new index values
    index_values = list(building_data.index.values)

    # Return total bounds of multipolygon
    bounds = projected_building.total_bounds

    # Numpy array dimensions based on x- and y-bounds, cell size, as integer
    x_size = int(abs(bounds[0] - bounds[2]) / cell_size)
    y_size = int(abs(bounds[1] - bounds[3]) / cell_size)

    # Create an empty raster, based on the bounds of the building data
    output = zeros(shape = (y_size, x_size))

    # Create an Affine transform object, using the polygon bounds (top left) and cell size
    transform_object = Affine(cell_size, 0, bounds[0], 0, -cell_size, bounds[3])
    
    # create tuples of geometry, value pairs, where value is the attribute value to burn
    building_heights = ((geom,value) for geom, value in zip(projected_building.geometry, building_data.relhmax))
    building_ids = ((geom,value) for geom, value in zip(projected_building.geometry, index_values))

    # Rasterise buildings, with height as the burn-in value
    building_raster = rasterize(building_heights,
                                out_shape = output.shape,
                                fill = 0,
                                out = None,
                                transform = transform_object,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)
    
    # Rasterise buildings, with id as the burn-in value
    building_id_raster = rasterize(building_ids,
                                out_shape = output.shape,
                                fill = 0,
                                out = None,
                                transform = transform_object,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)
    
    # Rasterise canal, with 1-0 as the burn-in value
    canal_raster = rasterize(projected_canal,
                                out_shape = output.shape,
                                fill = 0,
                                out = None,
                                transform = transform_object,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)
    
    # Write to raster
    if write:
    
        with rast_open(f"./buildings-raster.tif", "w",
            driver = "GTiff", crs = projected_building.crs,
            transform = transform_object, 
            dtype = float32,
            count = 1,
            width = x_size,
            height = y_size) as dst:
        
                dst.write(building_raster, indexes = 1)
        
        with rast_open(f"./canal-raster.tif", "w",
            driver = "GTiff", crs = projected_building.crs,
            transform = transform_object, 
            dtype = uint16,
            count = 1,
            width = x_size,
            height = y_size) as dst:
        
                dst.write(canal_raster, indexes = 1)

        with rast_open(f"./buildings-id-raster.tif", "w",
            driver = "GTiff", crs = projected_building.crs,
            transform = transform_object, 
            dtype = uint16,
            count = 1,
            width = x_size,
            height = y_size) as dst:
        
                dst.write(building_id_raster, indexes = 1)
    
    
    # Index of canal locations 
    canal_locations = argwhere(canal_raster)

    # Add canal locations to a canal set (for testing line-of-sight)
    canal_set = set(tuple(map(tuple, canal_locations)))
    
    # Convert the radius (m) to pixels
    radius_px = int(radius / cell_size)

    # Create final output array at the same dimensions as the building raster
    merged_viewshed = zeros(building_raster.shape)

    # Iterate through the canal locations
    for _ , c_coords in enumerate(canal_locations):

        # Return the viewshed for this coordinate
        output = viewshed(c_coords[0], c_coords[1], radius_px, cell_size, building_raster, canal_set)

        # Add output to final visibility raster
        merged_viewshed += output

        # Progress statement
        # print(f"Viewshed completed for {c_coords}, viewshed #{i}")

    # Write to raster
    if write:

        with rast_open(f"./buildings-visibility.tif", "w",
        driver = "GTiff", crs = projected_building.crs,
        transform = transform_object, 
        dtype = uint16,
        count = 1,
        width = x_size,
        height = y_size) as dst:

            dst.write(merged_viewshed, indexes = 1)


    # Normalises (0 = not visible, 1 = visible) and set data type
    merged_viewshed[merged_viewshed > 0] = 1
    merged_viewshed = merged_viewshed.astype(uint16)

    # Polygonise, using rasterio.features.shapes, if the cell is visible (1)
    geoms = list(({'properties': {'value': v}, 'geometry': s} for s, v in shapes(merged_viewshed, transform = transform_object) if v == 1))

    # Convert to GeoDataFrame, setting the crs
    visibility_gdf = GeoDataFrame.from_features(geoms, crs = projected_building.crs)

    # Dissolve for simplicity
    dissolved_gdf = visibility_gdf.dissolve()    

    # Spatial join
    visible_buildings = sjoin(building_data, dissolved_gdf, how='inner', predicate='intersects')
    visible_buildings.drop(columns = ['value', 'index_right'], inplace = True)

    # Return subset buildings, and the bounds
    return visible_buildings, bounds

def return_shading(date, filtered_buildings, all_buildings, canal_feature, plot, end_plot):
    '''
    > This function returns the shading geometry from buildings, for a specified datetime (UTC)
    > Utilises 'pybdshadow': https://pybdshadow.readthedocs.io/en/latest/
    > Results verified against: https://shadowmap.org/
    > Building geometries (OS) obtained via return_buildings()
    '''

    # Count of buildings
    n_buildings = len(filtered_buildings)

    # Keep unique ID, building height, and geometry
    building_data = filtered_buildings[['relhmax', 'geometry']]

    # Rename column to work with pybdshadow
    building_data = building_data.copy().rename(columns={'relhmax': 'height'})

    # Project to geographic coordinates
    building_data_projected = building_data.to_crs(4326)

    # This part is slow and results in the following warning:

        # > UserWarning: CRS not set for some of the concatenation inputs. \
        # > Setting output's CRS as WGS 84 (the single non-null crs provided)."

    # Preprocessing (remove empty polygons, multipolygons > polygons, generate building_id)
    #building_data_projected = pybdshadow.bd_preprocess(building_data_projected)

    # Explode multi-part geometries into multiple single geometries.
    building_data_projected = building_data_projected.explode(index_parts = False)

    # Add 'building_id' column
    building_data_projected['building_id'] = range(1, len(building_data_projected) + 1)

    # Calculate building shadows, directly using pybdshadow
    # shadows = pybdshadow.bdshadow_sunlight(building_data_projected, date)

    # Adapted to use PySolar
    shadows = bdshadow_sunlight(building_data_projected, date)

    # Assigns a CRS to the output (4326)
    shadows.crs = building_data_projected.crs

    # Reproject to British National Grid
    shadows = shadows.to_crs(building_data.crs)

    # If plot = True
    if plot:
                
        # plot for the first day only
        if date < end_plot:             
                    
            # If the image already exists
            if isfile(f'../images/shadows/shadows-for-{canal_feature.code_id}-at-{date.date().strftime("%m_%d")}-{date.time().strftime("%H_%M_%S")}.png'):
                pass

            # Otherwise, do the plotting
            else:

                # Store tuple geometry as series
                canal_object = GeoSeries(Polygon(canal_feature.geometry), crs = 27700)

                # Set up output image
                fig, my_ax = subplots(1, 1, figsize=(16, 10))
                title(f"Buildings-shadow test at (UTC) {date}. Buildings = {n_buildings}")

                # Set axis limits
                my_ax.set_xlim([all_buildings.total_bounds[0],  all_buildings.total_bounds[2]])
                my_ax.set_ylim([all_buildings.total_bounds[1],  all_buildings.total_bounds[3]])

                # Plot the shadow geometries
                shadows.plot(
                    ax = my_ax,
                    color = "#959595",
                    edgecolor = None,
                    linewidth = 0.5,
                    )
                
                # Plot the extracted building geometries
                all_buildings.plot(
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
                
                # Plot the canal geometries
                canal_object.plot(
                    ax = my_ax,
                    color = "#80B8FF",
                    edgecolor = '#3A90FE',
                    linewidth = 0.5,
                    )

                # Add north arrow
                x, y, arrow_length = 0.97, 0.99, 0.1
                my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                    arrowprops=dict(facecolor='black', width=5, headwidth=15),
                    ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)

                # Save the result
                savefig(f'../images/shadows/shadows-for-{canal_feature.code_id}-at-{date.date().strftime("%m_%d")}-{date.time().strftime("%H_%M_%S")}.png', bbox_inches='tight')

                # Close the figure
                close()

    # Return output
    return shadows

def determine_shading_proportion(shaded_feature, canal_feature, canal_area):
    '''
    > This function determines the proportion of shaded area for a feature
    '''

    try:

        # Clip the shaded area by the canal feature 
        clipped_shading = clip(shaded_feature, canal_feature)

    # Bare except
    except:

        # Make geometries valid (self-intersections)
        shaded_feature.geometry = shaded_feature.geometry.buffer(0)

        # Clip the shaded area by the canal feature 
        clipped_shading = clip(shaded_feature, canal_feature)

    # There is a shaded area but it does not intersect with the feature 
    if clipped_shading.empty:

        return 0 

    # There is a shaded area and it does intersect with the feature
    else:

        # Unary union 
        combined_shading = clipped_shading.unary_union

        # Calculate the area of the feature (BNG)
        shaded_area = combined_shading.area

        # Return this area, as a proportion of the total canal feature area (2 d.p)
        return round(shaded_area / canal_area, 2)


def return_distance(building_data, canal_data):
    '''
    > Function to return the distance from each building to the canal
    > Given shadow length = height / tan(solar altitude), we can quickly filter through the buildings
    > i.e. if the shadow length < the building distance, we don't need to calculate the actual shadow
    '''

    # For each building in the input, calculate the distance to the canal
    distances = [x.geometry.distance(canal_data.geometry) for x in building_data.itertuples()]

    # Store as a new attribute and return
    building_data['dist'] = distances
    return building_data


def filter_buildings_raster(building_id_raster, building_data, canal_raster, canal_rows, canal_cols, transformer, date, lon, lat, distance):
    '''
    > Function to filter buildings in the solar direction using numpy arrays
    '''

    # Create an empty raster, based on the canals raster
    solar_output = zeros(shape = canal_raster.shape)

    # Determine solar azimuth using PySolar (radians)
    azimuth_radians = radians(get_azimuth(lat, lon, date))

    # Sines and -cosines (correct y-offset)
    sines, coses = sin(azimuth_radians), -cos(azimuth_radians)

    # Distance in x- and y-dimensions (int), incorporating raster cell size     
    dx = int(sines * (distance / transformer[0]))
    dy = int(coses * (distance / transformer[0]))

    # Translated raster coordinates
    x_n = [dx + c for c in canal_cols]
    y_n = [dy + r for r in canal_rows]

    # Iterate through coordinates
    for start_row, start_col, end_row, end_col in zip(canal_rows, canal_cols, y_n, x_n):

        # Generate line pixel coordinates
        rr, cc = line(start_row, start_col, end_row, end_col)

        # Indexes of non-negative indices that are common to both row-col
        rr_index = nonzero(rr >= 0)
        cc_index = nonzero(cc >= 0)
        common_vals = intersect1d(rr_index, cc_index)

        # Only positive indices
        rr, cc = rr[common_vals], cc[common_vals]

        # Init object to control iteration 
        validity = False

        # Continue until the line is valid
        while validity is False:

            try:

                # Append to output
                solar_output[rr,cc] += 1

                # All coords in image space
                validity = True

            # Coords outside image space
            except IndexError:

                # Remove the last items (row/col) from the end of the lists
                rr, cc = rr[:-1], cc[:-1]

    # Normalises (0 = not visible, 1 = visible) and set data type
    solar_output[solar_output > 0] = 1
    solar_output = solar_output.astype(uint16)
    
    '''
    # Write to file for testing
    with rast_open(f"./buildings-visibility.tif", "w",
        driver = "GTiff", crs = building_data.crs,
        transform = transformer, 
        dtype = uint16,
        count = 1,
        width = solar_output.shape[1],
        height = solar_output.shape[0]) as dst:

            dst.write(solar_output, indexes = 1)
    '''

    # Coordinates of cells in the direction of the sun, and the corresponding building ids
    solar_rows, solar_cols = where(solar_output)
    solar_ids = building_id_raster[solar_rows, solar_cols]

    # Unique building IDs, that are not equal to 0
    unique_values = unique(solar_ids, axis=0)
    actual_values = unique_values[unique_values != 0]

    # Return the relevant rows in the gdf
    output_buildings = building_data[building_data.index.isin(actual_values)].copy()
    return output_buildings

def create_rasters(building_data, canal_data, bounds, cell_size, epsg):
    '''
    > Function to create a higher resolution canal and building id rasters based on the specified bounds and cell size
    > This is utilised for efficient building filtering based on the solar azimuth (filter_buildings_raster())
    '''

    # Convert tuple to GeoSeries
    projected_building = GeoSeries(building_data.geometry, crs = epsg)
    projected_canal = GeoSeries(canal_data.geometry, crs = epsg)

    # Convert polygon to line 
    projected_building = projected_building.exterior
    projected_canal = projected_canal.exterior

    # Reset index values
    building_data.reset_index(drop = False, inplace = True) 

    # Generate list of new index values
    index_values = list(building_data.index.values)

    # Numpy array dimensions based on x- and y-bounds of building data and cell size (integer)
    x_size = int(abs(bounds[0] - bounds[2]) / cell_size)
    y_size = int(abs(bounds[1] - bounds[3]) / cell_size)

    # Create an empty raster, based on the bounds of the building data
    output = zeros(shape = (y_size, x_size))

    # Create an Affine transform object, using the polygon bounds (top left) and cell size
    transform_object = Affine(cell_size, 0, bounds[0], 0, -cell_size, bounds[3])

    # create tuples of geometry, value pairs, where value is the attribute value to burn
    building_ids = ((geom,value) for geom, value in zip(projected_building.geometry, index_values))

    # Rasterise canal, with 1-0 as the burn-in value
    canal_raster = rasterize(projected_canal,
                                out_shape = output.shape,
                                fill = 0,
                                out = None,
                                transform = transform_object,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)
    
    # Rasterise building raster, with id as the burn-in value
    building_id_raster = rasterize(building_ids,
                                out_shape = output.shape,
                                fill = 0,
                                out = None,
                                transform = transform_object,
                                all_touched = False,
                                default_value = 1,
                                dtype = None)

    # Index of canal locations 
    canal_rows, canal_cols = where(canal_raster)

    # Return the rasters, indexes, Affine object, and building data (with reset ID)
    return canal_raster, building_id_raster, canal_rows, canal_cols, transform_object, building_data


def precompute_shading(canals_gdf, start_date, days, dist):
    '''
    > A function to iterate through canal features and return shading proportion
    > For each feature, iterate through datetimes (15m), return and filter buildings and calculate shading
    > Final output: % shaded area for each feature for each datetime
    > Store in JSON
    > Results key:
        > [sbh] = sun below the horizon (no shading)
        > [nb] = no buildings in the solar direction (no shading)
        > [sl] = there are buildings in the solar direction, but shadow length < distance to the canal
        > [m] = modelled shading (using pybdshadow, range = 0 - 1)
        > [pm] = shading proportion based on prior model (matching solar azimuth-altitude, at 1 d.p.)
    '''
    
    # Hours of analysis
    duration = days * 24

    # Add to postgis 
    # canals_gdf.to_postgis('canal_feature', db_connection, if_exists='overwrite', dtype={'geom': Geometry('[Polygon]', srid=27700)})

    # For model testing, we'll work on the first n alone
    #canals_gdf = canals_gdf.head(1)

    # Iterate efficiently using itertuples
    for feature_tuple in canals_gdf.itertuples():

        # Skip canals with an area < 500 m^2 or > 200_000 m^2
        if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
            continue
        
        # Try
        try:
            
            # Load shading record for current location
            with open(f"../../national-canal-cooling-data/shading/shading-{feature_tuple.code_id}.json") as json_path:
                shading_record = load(json_path)

        # If the file is not found
        except FileNotFoundError:

            # Create a new dictionary
            shading_record = {}
        
        # Localise niave datetime to aware datetime (UK) and then convert to UTC
        local_datetime = to_datetime(start_date)\
            .tz_localize('Europe/London')\
            .tz_convert('UTC')
        
        # Sets end datetime (+ duration)
        end_datetime = local_datetime + timedelta(hours = duration)

        # End plot date (+24h)
        end_plot_date = local_datetime + timedelta(hours = 24)
        
        # Convert current time to UTC
        start_unix, end_unix = str(int(local_datetime.timestamp())), str(int(end_datetime.timestamp()))
        
        # If the start and end are already in the dictionary
        if all(k in shading_record for k in (start_unix, end_unix)):

            # Skip to next canal feature
            print(f"All datetimes present for {feature_tuple.code_id}")
            continue

        # Return buildings within a specified distance (500 m) of feature polygon (using id) - this takes 146 seconds for all 2738 features, with 5 cores
        buildings = return_buildings_polygon(feature_tuple, dist, plot = False, epsg = canals_gdf.crs)

        # If there are no buildings within the specified distance
        if buildings.empty:

            # Error statement
            print(f"No buildings within {dist} m for feature {feature_tuple.code_id}")

        # There are buildings:
        else:
   
            # Subset buildings using a viewshed algorithm (line-of-sight), return canal raster
            buildings_los, building_bounds = subset_buildings(buildings, feature_tuple, False, cell_size = 5, radius = 500, epsg = canals_gdf.crs)

            # Add minimum distance column for each building, and the azimuth to each
            buildings_los = return_distance(buildings_los, feature_tuple)

            # Create 'high-resolution' rasters for subsequent filtering
            canal_raster, building_id_raster, canal_rows, canal_cols, transformer, buildings_los = create_rasters(buildings_los, feature_tuple, building_bounds, cell_size = 2, epsg = canals_gdf.crs)

        # Returns shapely area (27700 | 29902)
        feature_area = feature_tuple.geometry.area

        # Convert to series, project to 4326 and return lon-lat (perhaps could be simplified)
        geographic_geometry = GeoSeries(feature_tuple.geometry.centroid, crs = canals_gdf.crs).to_crs(4326)
        longitude, latitude = geographic_geometry.iloc[0].x, geographic_geometry.iloc[0].y

        # Continue until we have modelled the required duration
        while local_datetime <= end_datetime:

            # Convert current time to UTC
            unix_time = str(int(local_datetime.timestamp()))
    
            # If that datetime (UTC) exists in the dictionary
            if unix_time in shading_record:

                # Add 15 minutes to the datetime
                local_datetime += timedelta(hours=0.25)

            # This is a new datetime
            else:

                # If there are no buildings
                if buildings.empty:
                
                    # Add key (datetime) and value (no shading proportion) to dictionary
                    shading_record[unix_time] = {'shading' : 0,
                                                 'altitude' : round(solar_altitude, 1),
                                                 'azimuth' : solar_azimuth,
                                                 'result' : 'nb'}
                             
                    # Add 15 minutes to the datetime 
                    local_datetime += timedelta(hours=0.25)

                # If there are buildings, perform subsequent calculations
                else: 
                
                    # Calculate the solar altitude (in degrees) using PySolar
                    solar_altitude = get_altitude(latitude, longitude, local_datetime)

                    # Calculate the solar azimuth (in degrees) and round to 1 decimal places
                    solar_azimuth = round(get_azimuth(latitude, longitude, local_datetime), 1)

                    # If the sun is below the horizon
                    if solar_altitude <= 0:

                        # Add key (datetime) and value (no shading proportion) to dictionary
                        shading_record[unix_time] = {'shading' : 0,
                                                    'altitude' : round(solar_altitude, 1),
                                                    'azimuth' : solar_azimuth,
                                                    'result' : 'sbh'}
                                
                        # Add 15 minutes to the datetime 
                        local_datetime += timedelta(hours=0.25)

                    # If the sun is above the horizon (and shading is possible)
                    else:
                        
                        # Round the solar altitude to one decimal place
                        solar_altitude = round(solar_altitude, 1)
                        
                        # Use dictionary comprehension to check if we have modelled this altitude-azimuth pair befor
                        matched_dict = {k:v for (k,v) in shading_record.items() if v['altitude'] == solar_altitude and v['azimuth'] == solar_azimuth}
                        
                        # If the dictionary is empty (i.e. this is a new altitude-azimuth pair)
                        if not matched_dict:
                            '''
                            # Filter buildings using a polygon based on the solar azimuth (approach = "simplify" | "actual")
                            filtered_buildings = filter_buildings_polygon(buildings_los, canal_geographic, minimum_bounds, local_datetime, 
                                                                        longitude, latitude, 500, "simplify")
                            '''
                            
                            # Filter buildings based on the solar azimuth     
                            filtered_buildings = filter_buildings_raster(building_id_raster, buildings_los, canal_raster, canal_rows, canal_cols, transformer, 
                                                            local_datetime, longitude, latitude, 500)

                            # If there are no buildings in the solar direction
                            if filtered_buildings.empty:

                                # Add key (datetime) and value (no shading proportion) to dictionary
                                shading_record[unix_time] = {'shading' : 0,
                                                        'altitude' : solar_altitude,
                                                        'azimuth' : solar_azimuth,
                                                        'result' : 'nb'}

                                # Add 15 minutes to the datetime 
                                local_datetime += timedelta(hours=0.25)

                            # There are buildings which could shade the canal
                            else:
                                
                                ''' Using PySolar (rather than suncalc), steps broken up for readability '''

                                # Calculate the solar altitude in degrees, projecting to geographic
                                filtered_buildings['degs'] = [get_altitude(building.geometry.centroid.y, building.geometry.centroid.x, local_datetime) 
                                                            for building in filtered_buildings.to_crs(4326).itertuples()]
                    
                                # Calculate the tan() of the solar altitude for each building (radians)
                                filtered_buildings['solar_tan'] =  [tan(radians(x)) for x in filtered_buildings['degs']]

                                # Calculate shadow length, results now consistent with SunCalc (https://www.suncalc.org/)
                                filtered_buildings['shadow'] = filtered_buildings['relhmax'] / filtered_buildings['solar_tan']

                                # Return the difference between the shadow length and the building distance
                                filtered_buildings['diff'] = filtered_buildings['dist'] - filtered_buildings['shadow']

                                # Filter to buildings where the shadow length > distance to the canal
                                filtered_buildings_shade = filtered_buildings[filtered_buildings['diff'] < 0]

                                # If there are none
                                if filtered_buildings_shade.empty:

                                    # Add key (datetime) and value (no shading proportion) to dictionary
                                    shading_record[unix_time] = {'shading' : 0,
                                                            'altitude' : solar_altitude,
                                                            'azimuth' : solar_azimuth,
                                                            'result' : 'sl'}

                                    # Add 15 minutes to the datetime 
                                    local_datetime += timedelta(hours=0.25)

                                # Perform shading calculations (slow)
                                else:
                                
                                    # Return shading geometry based on the datetime and filtered building geometries (plot = True | False)
                                    shaded_geometry = return_shading(local_datetime, filtered_buildings_shade, buildings_los, feature_tuple, False, end_plot_date)

                                    # Returns shaded proportion (as a % of feature area)
                                    shading_proportion = determine_shading_proportion(shaded_geometry, feature_tuple.geometry, feature_area)

                                    # Add key (datetime) and value (shading proportion) to dictionary
                                    shading_record[unix_time] = {'shading' : shading_proportion,
                                                        'altitude' : solar_altitude,
                                                        'azimuth' : solar_azimuth,
                                                        'result' : 'm'}

                                    # Add 15 minutes to the datetime 
                                    local_datetime += timedelta(hours=0.25)

                        # We have modelled this altitude-azimuth pair before
                        else:
                            
                            # Add key (datetime) and previously calculated shading proportion to dictionary
                            shading_record[unix_time] = {'shading' : list(matched_dict.values())[0]['shading'],
                                                    'altitude' : solar_altitude,
                                                    'azimuth' : solar_azimuth,
                                                    'result' : 'pm'}

                            # Add 15 minutes to the datetime 
                            local_datetime += timedelta(hours=0.25)
   
        # When this feature has been modelled for the desired duration, save dict to JSON
        dump(shading_record, open(f"../../national-canal-cooling-data/shading/shading-{feature_tuple.code_id}.json", 'w'))  

        # Output statement
        print(f"Completed {feature_tuple.code_id}")


def extract_highly_shaded(threshold):
    '''
    > Function to return a list of canal features which exceed a defined shading threshold
    > e.g., annual median daytime shading > threshold (0.2 = 20%)
    > Note to self: the median shading values differ from the outputs of shading_effects_regression() - the latter is focused on May-September only
    '''

    # Check for output
    if isfile("../outputs/feature-shading-summary.json"):
        print("../outputs/feature-shading-summary.json EXISTS. Loading...")

        # Load shading summary file
        with open(f"../outputs/feature-shading-summary.json") as json_path:
            output = load(json_path)

    # If missing:
    else:

        # File paths for Britain and Ireland
        file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", "../data/urban-canals-ireland-dissolved-id-tm65.shp"]

        # Init output dict
        output = {}

        # Iterate and read each file
        for file in file_names:
            canals_gdf = read_dataframe(file)

            # Iterate efficiently using itertuples
            for feature_tuple in canals_gdf.itertuples():

                # Skip canals with an area < 500 m^2 or > 200_000 m^2
                if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
                    continue
                
                # Load shading record for current location
                try:
                    with open(f"../../national-canal-cooling-data/shading/shading-{feature_tuple.code_id}.json") as shading_path:
                        shading_record = load(shading_path)

                # If a file is missing 
                except FileNotFoundError:

                    # Skip to next canal feature
                    print(f"Missing shading files for feature {feature_tuple.code_id}")
                    continue

                # If there is an error with the file
                except decoder.JSONDecodeError:
                    
                    # Skip to next canal feature
                    print(f"Error with feature {feature_tuple.code_id}")
                    continue
                
                # Init dict for storing daytime values
                peak_shading = {}

                # Iterate through the shading dict
                for key, value in shading_record.items():

                    # Convert unix to datetime
                    dat = datetime.fromtimestamp(int(key))

                    # Between 10:00 and 15:45 (5h 45 m)
                    if dat.hour >= 10 and dat.hour < 16:

                        # Add to daytime dictionary
                        peak_shading[key] = value

                # Init lists to store week, fortnight and monthly data
                day_week_list = defaultdict(list)
                day_fortnight_list = defaultdict(list)
                day_month_list = defaultdict(list)

                # Iterate through the daytime values
                for key, value in peak_shading.items():

                    # Convert unix to datetime, extract month integer
                    month = datetime.fromtimestamp(int(key)).month

                    # Extract the day
                    day = datetime.fromtimestamp(int(key)).timetuple().tm_yday

                    # Week of the year
                    week = ceil(day/7)

                    # If the week is even
                    if week % 2 == 0:
                        fortnight = week / 2
                    # The week is odd
                    else:
                        fortnight = (week+1) / 2

                    # Append values to corresponding list with key
                    day_month_list[month].append(value['shading'])

                    # Ignore final week and fortnight 
                    if week < 53:
                        day_week_list[week].append(value['shading'])
                    if fortnight < 27:
                        day_fortnight_list[fortnight].append(value['shading'])

                # Extract all shading values as list
                day_shading = [x['shading'] for x in peak_shading.values()]
                
                # Add summary values to output
                output[feature_tuple.code_id] = {# Annual daytime shading values
                                                'shading values' : {'median' : median(day_shading), 
                                                                'q1' : quantile(day_shading, 0.25),
                                                                'q3' : quantile(day_shading, 0.75),
                                                                'mean' : mean(day_shading), 
                                                                'count' : sum(i > 0.8 for i in day_shading)}, # Count of intervals where shading > 80%
                                                # Week values (daytime)
                                                'weekly values' : {'medians' : [median(day_week_list[x]) for x in day_week_list]}, 
                                                # Fortnight values (daytime)  
                                                'fortnightly values' : {'medians' : [median(day_fortnight_list[x]) for x in day_fortnight_list]},
                                                # Monthly values, daytime and nighttime (median, Q1:Q3)             
                                                'monthly values' : {'medians' : [median(day_month_list[x]) for x in day_month_list], 
                                                                'q1' : [quantile(day_month_list[x], 0.25) for x in day_month_list],
                                                                'q3' : [quantile(day_month_list[x], 0.75) for x in day_month_list]}}

                                                
                # Progress statement
                print(f"Completed feature {feature_tuple.code_id}")  

        # When all features have been summarised, save to outputs
        dump(output, open(f"../outputs/feature-shading-summary.json", 'w'))

    # Filter to highly shaded, and simplify dict
    highly_shaded = {k:v['shading values']['median'] for (k,v) in output.items() if v['shading values']['median'] > threshold}

    # When all features have been summarised, save to outputs
    dump(highly_shaded, open(f"../outputs/high-shading-features.json", 'w'))
