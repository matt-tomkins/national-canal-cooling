''' Functions to analyse and present the model results, produced by national-canal-cooling.py, which are stored in /results 
    > Summary files stored in /outputs
    > Figures stored in /images '''

# Force use of Shapely 2.0
from os import environ, remove
environ['USE_PYGEOS'] = '0'

# Imports
import imageio
from scipy import stats
from copy import deepcopy
from glob import glob as g
from pandas import read_csv
from pywaffle import Waffle
from os.path import basename
from datetime import datetime
from statistics import median
from json import load, decoder
from pyproj import Transformer
import matplotlib.dates as mdates
from geopandas import GeoSeries
from pyogrio import read_dataframe
from collections import defaultdict
import matplotlib.patches as mpatches
from timeit import default_timer as timer
from dateutil.relativedelta import relativedelta
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import Rectangle, subplots, savefig, rcParams, rc, figure, setp
from numpy import arange, quantile, linspace, loadtxt, array, sum, abs, ediff1d, argsort, cumsum, searchsorted

# Import user defined functions
from params import REFERENCE_MATERIAL
from auxiliary import *
from functions import *

# Set font family
rc('font', family='serif')

# Parameters for extracting peak day (10:00 - 15:45) and nighttime periods (22:00 - 03:45)
# This approach provides a consistent basis for comparison across the year (5h 45 m), rather than utilising variable sunrise / sunset times.
DAY_START, DAY_END  = 10, 16
NIGHT_START, NIGHT_END = 22, 4

# 2022 heatwave days ('actual')
HEATWAVE_ACTUAL =  [['166', '167', '168'], #--------------------------------- 15–17 June
                    ['198', '199', '200'], #--------------------------------- 17–19 July      
                    ['221', '222', '223', '224', '225', '226', '227']] #----- 9–15 August
    
# 2022 heatwave days ('extended' ± 2 days)
HEATWAVE_EXTENDED = [['164', '165', '166', '167', '168', '169', '170'], 
                     ['196', '197', '198', '199', '200', '201', '202'], 
                     ['219', '220' '221', '222', '223', '224', '225', '226', '227', '228', '229']] 


# File names for canal shapefiles in Britain and Ireland
FILE_NAMES = ["../data/canal-geometries-filtered-modified-dissolved-id-width.shp", 
             "../data/urban-canals-ireland-dissolved-id-tm65-width.shp"]

#============== Interpretation of sensitivity testing ==============#  

def main():

    # Initialise timer
    start = timer()

    #==================== Analysis ====================#
    
    # [1] Extract summary values for each region
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    #-----> 'region_save' (Boolean): Save (True) intermediate outputs for each region (energy flux, temperature estimate for each datetime)
    #-----> 'buffer_save' (Boolean): Save (True) geometry of buffer areas for each region
    #-----> 'period' (String): For flagging heatwaves, utilise heatwave days only ('actual'), include a longer period ('extended', +/- 2 days), or don't perform analysis ('none', or any other string)
    # model_urban_cooling(True, True, False, "actual")

    # [2] Print summary statistics to console and save to file
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    print_summary_values(True)

    # [3] Analyse temperature effects during heatwaves (15–17 June, 17–19 July and 9–15 August)
    #-----> Requires output of model_urban_cooling(), with 'region_save' = "True" and 'period' = "actual" | "extended"
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    # heatwave_analysis_region(True)

    # [*] Model validation, statistical similarity    
    # Method = ["Spearman", "Pearson", "Kendall"]
    # Not called directly, but utilised in draw_figure_3()
    # time_series_correlation(True, 'lalc_73', depth=40, method="Spearman") # 2024-03-08

    #==================== Figures ====================#

    # Figure 1: Waffle plot and map of urban canal network
    # draw_figure_1()
    
    # Figure 2: Time-series of modelled wT and cT, plus map of feature
    #-----> 'canal_id' (String): Feature name to plot e.g., 'rc2_23'
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    #-----> 'month' (Numeric): Month to highlight [1-12]
    # draw_figure_2("rc2_23", True, 7) 

    # Figure 3: Plot measured vs. modelled water temperatures
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    #-----> 'canal_list' (List): List of canal IDs
    #-----> 'depth' (Numeric): Water depth (cm)
    # draw_figure_3(True, ['lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11'], depth=40) 

    # Figure 4: Fluctuations in median air temperature change at desired interval, plus histogram
    #-----> 'start_date' (String): Start date for plotting ("2022-01-01 10:00:00")
    #-----> 'interval' (String): Interval for aggregration ("week", "fortnight", "month")
    #-----> 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    # draw_figure_4("2022-01-01 10:00:00", "month", True)

    # Create gif showing shadows cast by buildings
    #-----> 'feature_name' (String): Feature ID to model for
    #-----> 'start_date' (String): Start date for plotting ("2022-06-01 10:00:00")
    #-----> 'days' (Numeric): Duration of analysis (days)
    #-----> 'delete' (Boolean): Delete (True) or keep (False) image files after completion
    # shading_gif('ld_16', "2022-06-01 10:00:00", 1, True)

    #==================== Supplementary ====================#

    # [*] Model spin-up, selecting the feature to plot
    #-----> 'canal_id' (String): Feature ID to plot
    # draw_supplementary_figure_A("rc2_23")

    # [*] Plot to showcase the difference in water and reference temperatures when shading is included
    #-----> 'canal_id' (String): Feature ID to plot e.g., 'rc_9' (high shading example) or 'lalc_170' (low shading example)
    # draw_supplementary_figure_B('lalc_170')

    # Completion
    print(f"Code completion in {timer() - start} seconds")

#==================== Analysis ====================#


def model_urban_cooling(shading_boolean, region_save, buffer_save, period):
    '''
    > Originally, air temperature effects were modelled on a feature-wise basis (n = 2,356)
    > However, this approach is sensitive to canal geometry and in particular, the frequency of bridges and locks (etc) which split the canal (arbitrarily) into independent features
    > To account for this, we aggregate the canal features by city / region
    > Modelled temperatures are derived from the net energy flux, which is the difference between the sensible flux associated with water-to-air and the equivalent flux for reference-to-air
    > Parameters:
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
        - 'region_save' (Boolean): Save (True) intermediate outputs for each region (energy flux, temperature estimate for each datetime)
        - 'buffer_save' (Boolean): Save (True) geometry of buffer areas for each region
        - 'period' (String): For flagging heatwaves, utilise heatwave days only ('actual'), include a longer period ('extended', +/- 2 days), or don't perform analysis ('none')
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # If the summary file already exists, exit the function
    if isfile(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json"):
        print(f"Output file: ../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json EXISTS. Exiting...")
        return
    
    # Heatwave days only
    if period == 'actual':
        heatwave_periods =  HEATWAVE_ACTUAL
    
    # Include an extended period (± 2 days)
    elif period == 'extended':
        heatwave_periods =  HEATWAVE_EXTENDED

    # Create a deep copy for nightttimes
    nightwave_periods = deepcopy(heatwave_periods)

    # Init list to store end days
    heatwave_ends = []

    # Iterate through heatwaves
    for sl in nightwave_periods:

        # Store maximum
        heatwave_ends.extend(str([int(max(sl))+1]))

        # Extend by +1 day for nighttime (including the next night)
        sl.extend([int(max(sl))+1])

    # Init output dictionary
    output = {}

    # Read the dataframes (Britain, Ireland) directly using pyogrio
    for f in FILE_NAMES:
        canals_gdf = read_dataframe(f)

        # Region list
        regions = canals_gdf['region'].unique()

        # Iteratively subset to region
        for r in regions:
            subset_canals = canals_gdf.loc[canals_gdf['region'].isin([r])]

            # Init dict for storing values by region, and lists for storing feature-buffer geometries
            region_dict = {}
            region_geoms = []
            region_buffers = []
        
            # Iterate using itertuples
            for feature_tuple in subset_canals.itertuples():

                # Try to load results for this feature
                try:
                    # Load results excluding | including shading
                    with open(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{feature_tuple.code_id}.json") as canal_data_path:
                        canal_data = load(canal_data_path)

                # If the file is missing (e.g, < 500 m^2 | > 200_000 m^2)
                except FileNotFoundError:

                    # Skip to next canal feature
                    # print(f"No file for {feature_tuple.code_id}.json")
                    continue

                # If there is an error with the file
                except decoder.JSONDecodeError:
                    
                    # Skip to next canal feature
                    # print(f"Error with {feature_tuple.code_id}.json")
                    continue

                # Store the geometry
                region_geoms.append(feature_tuple.geometry)

                # Buffer by empirically derived cooling distance (1.5126 * feature_tuple.width + 6.7227)
                # Source: Hathaway and Sharples (2012), Murakawa et al. (1991)
                region_buffers.append(feature_tuple.geometry.buffer(1.5126 * feature_tuple.width + 6.7227))
                
                # If the dict contains data
                if bool(region_dict):
                    
                    # Iterate over the keys, and add the net sensible flux (water-referenec) from new dict
                    for key in region_dict:  
                        region_dict[key]['net_flux'] += canal_data[key]['water_sensible'] - canal_data[key]['reference_sensible']

                # If this is the first feature for this region
                else:
        
                    # Save a simplified dictionary {datetime : {net_flux : value}}
                    region_dict = {outer_key: ({'net_flux': inner_dict['water_sensible'] - inner_dict['reference_sensible']})
                                for outer_key, inner_dict in canal_data.items()
                                }
                             
            # Once the region has been completed, and if the dict contains data      
            if bool(region_dict):

                # Create a geodataframe for this region (no attributes)
                region_gdf = gpd.GeoDataFrame(geometry = region_geoms, crs = subset_canals.crs)

                # Create gdf of buffers and dissolve
                buffer_gdf = gpd.GeoDataFrame(geometry = region_buffers, crs = subset_canals.crs).dissolve()

                # Calculate area
                buffered_area = buffer_gdf.area.iloc[0]

                # Save buffers for region
                if buffer_save:
                    buffer_gdf.to_file(f"../outputs/{REFERENCE_MATERIAL}/region-values/buffered_area_{r}.shp")

                # Heat capacity of the air volume (in J/K), i.e., energy required to change the temperature by 1 K (thermal capacity * mass of air)
                # Utises specific heat capacity (kJ/kg*k) converted to J/kg*k, and fixed air density
                energy_per_kelvin = (A_CAPACITY * 1000) * (A_DENSITY * (buffered_area * A_HEIGHT))

                # Iterate and divide energy flux by energy per kelvin
                for key, value in region_dict.items():
                    region_dict[key]['net_temp'] = value['net_flux'] / energy_per_kelvin

                # When all features for this region have been combined, add heatwave flag
                if region_save and period in ("actual", "extended"):
                    
                    # Iterate through keys
                    for key, value in region_dict.items():

                        # Convert unix to datetime, extract day of year
                        local_datetime = datetime.fromtimestamp(int(key))
                        current_day = str(local_datetime.timetuple().tm_yday)

                        # If a heatwave day and during the daytime, classified as 10:00-15:45 (5h 45 m)
                        if any(current_day in d for d in heatwave_periods) and local_datetime.hour >= DAY_START and local_datetime.hour < DAY_END:

                            # Add flag to output
                            region_dict[key]['heatwave'] = f"{local_datetime.month} day"

                        # Filter to heatwave days + 1 (including the night at the end) and nighttime, classified as 22:00-03:45 (5h 45 m)
                        elif any(current_day in d for d in nightwave_periods) and (local_datetime.hour >= NIGHT_START or local_datetime.hour < NIGHT_END):

                            # Remove the extra night (i.e., after 22:00 the following day)
                            if any(current_day in d for d in heatwave_ends) and local_datetime.hour >= NIGHT_START: 
                                continue
                                
                            # Add flag to output
                            else:
                                region_dict[key]['heatwave'] = f"{local_datetime.month} night"

                    # When complete, save to outputs
                    dump(region_dict, open(f"../outputs/{REFERENCE_MATERIAL}/region-values/region-values-{r}-{folder_path}-shading.json", 'w'))
                
                # Save to outputs, without heatwave flag
                elif region_save: 
                    dump(region_dict, open(f"../outputs/{REFERENCE_MATERIAL}/region-values/region-values-{r}-{folder_path}-shading.json", 'w'))

                # Init dicts for storing values
                peak_day = {}
                peak_night = {}

                # Iterate through the data dict
                for key, value in region_dict.items():

                    # Convert unix to datetime
                    dat = datetime.fromtimestamp(int(key))
                
                    # Between 10:00 and 15:45 (5h 45 m)
                    if dat.hour >= DAY_START and dat.hour < DAY_END:

                        # Add to daytime dictionary
                        peak_day[key] = value

                    # Between 22:00 and 03:45 (5h 45m)
                    elif dat.hour >= NIGHT_START or dat.hour < NIGHT_END:

                        # Add to nightime dictionary
                        peak_night[key] = value
                        
                # Init lists to store week, fortnight and monthly data
                day_week_list = defaultdict(list)
                day_fortnight_list = defaultdict(list)
                day_month_list = defaultdict(list)
                
                # As above, but for night values
                night_week_list = defaultdict(list)
                night_fortnight_list = defaultdict(list)
                night_month_list = defaultdict(list)

                # Iterate through the daytime values
                for key, value in peak_day.items():

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
                    day_month_list[month].append(value['net_temp'])

                    # Ignore final week and fortnight 
                    if week < 53:
                        day_week_list[week].append(value['net_temp'])
                    if fortnight < 27:
                        day_fortnight_list[fortnight].append(value['net_temp'])

                # Iterate through the nighttime values
                for key, value in peak_night.items():

                    # Convert unix to datetime, extract month integer
                    month = datetime.fromtimestamp(int(key)).month

                    # Append values to corresponding key
                    night_month_list[month].append(value['net_temp'])

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

                    # Ignore final week and fortnight 
                    if week < 53:
                        night_week_list[week].append(value['net_temp'])
                    if fortnight < 27:
                        night_fortnight_list[fortnight].append(value['net_temp'])

                # Extract air temperatures
                day_temps = [x['net_temp'] for x in peak_day.values()]
                night_temps = [x['net_temp'] for x in peak_night.values()]

                # Add summary values to output
                output[r] = {# Geometry attributes
                            'geometry' : {'area' :  region_gdf.dissolve().area.iloc[0], 
                                        'buffer area' : buffered_area},
                            # Median (Q1:Q3) daytime values
                            'day values' : {'median' : median(day_temps), 
                                            'q1' : quantile(day_temps, 0.25),
                                            'q3' : quantile(day_temps, 0.75)},
                            # Median (Q1:Q3) nighttime values
                            'night values' : {'median' : median(night_temps), 
                                            'q1' : quantile(night_temps, 0.25),
                                            'q3' : quantile(night_temps, 0.75)},
                            # Week values (day and night)
                            'weekly values' : {'day medians' : [median(day_week_list[x]) for x in day_week_list],
                                            'night medians' : [median(night_week_list[x]) for x in night_week_list]}, 
                            # Fortnight values (day and night)  
                            'fortnightly values' : {'day medians' : [median(day_fortnight_list[x]) for x in day_fortnight_list],
                                                    'night medians' : [median(night_fortnight_list[x]) for x in night_fortnight_list]},
                            # Monthly values, daytime and nighttime (median, Q1:Q3)             
                            'monthly values' : {'day medians' : [median(day_month_list[x]) for x in day_month_list], 
                                            'day q1' : [quantile(day_month_list[x], 0.25) for x in day_month_list],
                                            'day q3' : [quantile(day_month_list[x], 0.75) for x in day_month_list],
                                            'night medians' : [median(night_month_list[x]) for x in night_month_list],
                                            'night q1' : [quantile(night_month_list[x], 0.25) for x in night_month_list],
                                            'night q3' : [quantile(night_month_list[x], 0.75) for x in night_month_list]}}


                # Progress statement
                print(f"Completed {r}")

    # Store output
    dump(output, open(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json", 'w'))


def print_summary_values(shading_boolean):
    '''
    > Utilises summary data for each region produced via model_urban_cooling(), either including or excluding shading effects
    > Returns key summary values e.g., annual daytime and nighttime medians, monthly summaries
    > Results printed to console and saved to file "../outputs/summary-values-{folder_path}-shading.json"
    > Parameters:
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # If the summary file already exists, exit the function
    if isfile(f"../outputs/{REFERENCE_MATERIAL}/summary-values-{folder_path}-shading.json"):
        print(f"Output file: '../outputs/{REFERENCE_MATERIAL}/summary-values-{folder_path}-shading.json' EXISTS. Loading...")
       
        # Load file
        with open(f"../outputs/{REFERENCE_MATERIAL}/summary-values-{folder_path}-shading.json") as file_path:
            output = load(file_path)

            # List of partial keys
            key_list = ["Annual", "January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

            # Title
            print(f"Results {folder_path} shading effects ==================")

            # Iterate, extract day-night data and print
            for k in key_list:
                day = output[f"{k} day"]
                night = output[f"{k} night"]
                print(f"For {k}, median day ΔPT = {day['median']:.2f} (Q1 = {day['q1']:.2f}, Q3 = {day['q3']:.2f}), night ΔPT = {night['median']:.2f} (Q1 = {night['q1']:.2f}, Q3 = {night['q3']:.2f})")

        # Exit function
        return    
    
    # Try and load summary data, including | excluding shading
    try:
        with open(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json") as summary_path:
            summary = load(summary_path)

    # If file does not exist, exit function
    except FileNotFoundError:
        print(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json does not exist")
        return 
    
    # Initialise output dictionary
    output_dict = {}

    # Extract and print median daytime values
    day_medians = array([x['day values']['median'] for x in summary.values()])
    areas = array([x['geometry']['area'] for x in summary.values()])

    # Calculate day median, p25, p75, weighted based on canal area
    wt_day_quantiles = weighted_median(day_medians, areas, [0.25, 0.5, 0.75])

    print(f"Day statistics (annual) {folder_path} shading -----------------")
    print(f"Median day ΔPT of {wt_day_quantiles[1]:.2f}. Q1 = {wt_day_quantiles[0]:.2f}, Q3 = {wt_day_quantiles[2]:.2f}")
    print()

    # Save to output
    output_dict.update({'Annual day' : {'median' : wt_day_quantiles[1],
                                    'q1' : wt_day_quantiles[0],
                                    'q3' : wt_day_quantiles[2]}})
    
    # Extract and print median nightime values
    night_medians = array([x['night values']['median'] for x in summary.values()])

    # Calculate night median, p25, p75, weighted based on canal area
    wt_night_quantiles = weighted_median(night_medians, areas, [0.25, 0.5, 0.75])

    print(f"Night statistics (annual) {folder_path} shading -----------------")
    print(f"Median night ΔPT of {wt_night_quantiles[1]:.2f}. Q1 = {wt_night_quantiles[0]:.2f}, Q3 = {wt_night_quantiles[2]:.2f}")
    print()

    # Save to output
    output_dict.update({'Annual night' : {'median' : wt_night_quantiles[1],
                                    'q1' : wt_night_quantiles[0],
                                    'q3' : wt_night_quantiles[2]}})



    # Extract median values for each month (daytime + nightime)
    mon_day_medians = array([x['monthly values']['day medians'] for x in summary.values()])
    mon_night_medians = array([x['monthly values']['night medians'] for x in summary.values()])

    # Calculate monthly day medians, p25, p75, weighted based on canal area
    month_wt_day_quantiles = [weighted_median(mon_day_medians[:,x], areas, [0.25, 0.5, 0.75]) for x in range(12)]

    # As above, for night 
    month_wt_night_quantiles = [weighted_median(mon_night_medians[:,x], areas, [0.25, 0.5, 0.75]) for x in range(12)]

    # List of months
    months = ["January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

    print(f"Statistics (monthly) {folder_path} shading -----------------")
    # Iterate and print summary statitics
    for index, m in enumerate(months):
        print(f"For {m}, the daytime ΔPT is {month_wt_day_quantiles[index][1]:.2f} ({month_wt_day_quantiles[index][0]:.2f} - {month_wt_day_quantiles[index][2]:.2f}), with nighttime ΔPT of {month_wt_night_quantiles[index][1]:.2f} ({month_wt_night_quantiles[index][0]:.2f} - {month_wt_night_quantiles[index][2]:.2f})")

        # Save to output (daytime)
        output_dict.update({f"{m} day" : {'median' : month_wt_day_quantiles[index][1], 
                                    'q1' : month_wt_day_quantiles[index][0],
                                    'q3' : month_wt_day_quantiles[index][2]}})
        
        # Save to output (nighttime)
        output_dict.update({f"{m} night" : {'median' : month_wt_night_quantiles[index][1], 
                                    'q1' : month_wt_night_quantiles[index][0],
                                    'q3' : month_wt_night_quantiles[index][2]}})
                                 
    # If the output dictionary contains data
    if bool(output_dict):

        # When all features have been summarised, save to outputs
        dump(output_dict, open(f"../outputs/{REFERENCE_MATERIAL}/summary-values-{folder_path}-shading.json", 'w'))

    # No output
    else:
        print(f"The output dictionary is EMPTY. '../outputs/{REFERENCE_MATERIAL}/summary-values-{folder_path}-shading.json' has not been created")
        return
    
def weighted_median(values, weights, quantiles=0.5):
    '''
    Function to calculate the weighted median and other quantiles
    Source: https://stackoverflow.com/questions/20601872/numpy-or-scipy-to-calculate-weighted-median
    '''
    i = argsort(values)
    c = cumsum(weights[i])
    return values[i[searchsorted(c, array(quantiles) * c[-1])]]

def heatwave_analysis_region(shading_boolean):
    '''
    > Function to summarise the temperature effects during heatwaves in June, July and August
    > Requires output of model_urban_cooling(), with 'region_save' = "True" and 'period' = "actual" | "extended"
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # If the heatwave file exists
    if isfile(f"../outputs/{REFERENCE_MATERIAL}/heatwave-values-{folder_path}-shading.json"):
        
        # Inform user and load
        print(f"../outputs/{REFERENCE_MATERIAL}/heatwave-values-{folder_path}-shading.json EXISTS. Loading...")
        with open(f"../outputs/{REFERENCE_MATERIAL}/heatwave-values-{folder_path}-shading.json") as file_path:
            output_dict = load(file_path)
    
    # Summarise by region and by heatwave
    else:

        # Open summary file, containing areas
        with open(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json") as region_path:
            region_dict = load(region_path)

        # File paths, matching the folder path, with wild card for the region name
        files = g(f"../outputs/{REFERENCE_MATERIAL}/region-values/region-values-*-{folder_path}-shading.json", recursive = True)

        # List of months
        months = ["January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

        # Init dict to store flags and datetimes
        flag_dict = {}
        output_dict = {}

        # Iterate and isolate region name
        for f in files:
            region_name = basename(f).removeprefix(f"region-values-").removesuffix(f"-{folder_path}-shading.json")

            # Load file
            with open(f) as file_path:
                data = load(file_path)

                # If the flag dict is empty 
                if not flag_dict:

                    # Extract unique heatwave flags
                    unique_flags = {v["heatwave"] for v in data.values() if "heatwave" in v}

                    # Iterate and store corresponding keys (datetimes) in dict for efficiency
                    for flag in unique_flags:
                        matching_keys = {k for k, v in data.items() if v.get("heatwave") == flag}
                        flag_dict[flag] = {"datetimes" : matching_keys}

                # Iterate through heatwave flags
                for key, values in flag_dict.items():
                    
                    # Subset to relevant keys
                    subset = {k: data[k] for k in values['datetimes'] if k in data}

                    # Extract values
                    net_temp = [x['net_temp'] for x in subset.values()]

                    # Formatted key
                    key_name = f"{months[int(key.split()[0])-1]} {key.split()[1]}"

                    # Add quantiles to output dict
                    try:
                        output_dict[region_name][key_name] = {'median' : median(net_temp),
                                                            'q1' : quantile(net_temp, 0.25),
                                                            'q3' : quantile(net_temp, 0.75)}
                    except KeyError:
                        output_dict[region_name] = {key_name : {'median' : median(net_temp),
                                                                'q1' : quantile(net_temp, 0.25),
                                                                'q3' : quantile(net_temp, 0.75)}}
            
            # Add area to output dict
            output_dict[region_name]['geometry'] = {'area' : region_dict[region_name]['geometry']['area']}

        # When complete, save to output dir
        dump(output_dict, open(f"../outputs/{REFERENCE_MATERIAL}/heatwave-values-{folder_path}-shading.json", 'w'))      

    # Extract areas for weighting
    areas = array([x['geometry']['area'] for x in output_dict.values()])  

    # Area-weighted calculation, for month-period combinations
    for month in ["June", "July", "August"]:
        for period in ['day', 'night']:

            # Extract median daytime values
            day_medians = array([x[f"{month} {period}"]['median'] for x in output_dict.values()])
   
            # Calculate median, p25, p75, weighted based on canal area
            wt_day_quantiles = weighted_median(day_medians, areas, [0.25, 0.5, 0.75])

            # Print summary
            print(f"For {month} ({period}), median ΔPT of {wt_day_quantiles[1]:.2f}. Q1 = {wt_day_quantiles[0]:.2f}, Q3 = {wt_day_quantiles[2]:.2f}")

def time_series_correlation(shading_boolean, canal_id, depth, method):
    '''
    > Function for measuring similarity of time-series (modelled water temperature vs. measured water temperature) using measurement data from the CRT (https://canalrivertrust.org.uk/)
    > Parameters:
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
        - 'canal_id' (String): ID for relevant canal feature
        - 'depth' (Numeric): Water depth (cm) which refers to the upper value:
            - 20 = 20 - 40 cm
            - 40 = 40 - 60 cm
            - 60 = 60 - 80 cm
            - 80 = 80 - 100 cm
        - 'method' (String), choice of statistical method:
            - Pearson
            - Spearman (preferred)
            - Kendall 
    > Sources:
        - https://stats.stackexchange.com/questions/133155/how-to-use-pearson-correlation-correctly-with-time-series
        - https://otexts.com/fpp2/stationarity.html
    > Of the monitoring locations, the following are suitable for comparison:
        - Bridge Pagefield Pipe Crossing (lalc_73) --------- 10 minute frequency
        - Bridge 15 Chain Lane Bridge (nabc_2) ------------- 10 minute frequency
        - New Road Bridge 148 (suc_41) --------------------- 1 hour frequency
        - Bridge 32. Priorswood (batc_14) ------------------ 15 minute frequency (from 03-01-2022)
        - Bridge 4 Wakefield Road Bridge (hc_53) ----------- 1 hour frequency
        - Bridge 73 Anchor Bridge (cc3_11) ----------------- 10 minute frequency
    > The following locations are less suitable for comparison:
        - Anderton Waste Weir Flow (None) ------------------ not directly modelled
        - Bridge 46 Benthouse Bridge (None) ---------------- not directly modelled
        - Hawkesbury Lock (None) --------------------------- not directly modelled
        - Bridge C Lightbody Street (lalc_4) --------------- missing measurement data
        - Bridge 40 Taverners Bridge (cc3_53) -------------- missing measurement data
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Dictionary of CRT water temperature monitoring locations and corresponding canal ids
    location_dict = {'tamc_19' : 'Anderton Waste Weir Flow Water Temperature C (All Data) (Value)', 
                     'lalc_73' : 'Bridge Pagefield Pipe Crossing Mean Temperature C (All Data) (Value)',
                     'nabc_2' : 'Bridge 15 Chain Lane Bridge Mean Temperature C (All Data) (Value)',
                     'suc_41' : 'New Road Bridge 148 Mean Temperature C (All Data) (Value)',
                     'rc_103' : 'Bridge 46 Benthouse Bridge Mean Temperature C (All Data) (Value)',
                     'cc3_60' : 'Hawkesbury Lock Mean Temperature C (All Data) (Value)',
                     'batc_14' : 'Bridge 32. Priorswood Mean Temperature C (All Data) (Value)', 
                     'hc_53' : 'Bridge 4 Wakefield Road Bridge Mean Temperature C (All Data) (Value)',
                     'lalc_4' : 'Bridge C Lightbody Street Mean Temperature C (All Data) (Value)',
                     'cc3_11' : 'Bridge 73 Anchor Bridge Mean Temperature C (All Data) (Value)',
                     'cc3_53' : 'Bridge 40 Taverners Bridge Mean Temperature C (All Data) (Value)'}
    
    # Open csv
    measured_values = read_csv("../data/model_validation.csv")

    # Return value, corresponding to selected id
    location = location_dict[canal_id]
    
    # Extract measured temperatures and datetimes
    measured_water = measured_values[['dt', location]]
    measured_water = measured_water.dropna()

    # Convert from string to datetime
    dt = [datetime.strptime(x, '%d/%m/%Y %H:%M') for x in measured_water['dt']]

    # Convert to unix 
    dt_unix = [str(int(x.timestamp())) for x in dt]

    # Measured dictionary
    measured_dict = dict(zip(dt_unix, measured_water[location]))

    # Load model output record, including | excluding shading
    try:
        with open(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{canal_id}.json") as results_path:
            canal_output = load(results_path)

    # No file present
    except FileNotFoundError:
        print("No modelled data to compare")
        exit()

    # Init dicts to store output values
    measured_dict_filter = {}
    modelled_dict_filter = {}

    # Iterate through the measured values dict
    for key, value in measured_dict.items():

        # If this unix time exists in the modelled dictionary 
        if key in canal_output:
            
            # Add to output
            measured_dict_filter[key] = value
            modelled_dict_filter[key] = canal_output[key]

    # Extract values as numpy array
    measured_wt = array([x for x in measured_dict_filter.values()])

    # Extract modelled values as numpy array, converting to celcius
    try: 
        modelled_wt = array([x[f'depth_water_{depth}'] - 273.15 for x in modelled_dict_filter.values()])

    # EAFP: surface layer
    except KeyError:
        modelled_wt = array([x[f'surface_water_k'] - 273.15 for x in modelled_dict_filter.values()])

    # Return the differences between consecutive array elements
    measured_diff = ediff1d(measured_wt)
    modelled_diff = ediff1d(modelled_wt)

    # Method selection
    if method == "Spearman":
        r, p = stats.spearmanr(measured_diff, modelled_diff)

    elif method == "Pearson":
        r, p = stats.pearsonr(measured_diff, modelled_diff)

    elif method == "Kendall": 
        r, p = stats.kendalltau(measured_diff, modelled_diff)

    else:
        print("Statistical method invalid. Check inputs")

    # Print summary
    print(f"For {canal_id}, {method} correlation = {r:.2f}, p value = {p:.2f}")

    # Return method, test statistic,  p value, array length
    return method, r, p, len(measured_wt)

#==================== Figures ====================#

def draw_figure_1():
    '''
    > Draw Figure 1, which includes a waffle plot and map of the urban canal network
    '''

    # Read in the canal dataframes
    britain_gdf = read_dataframe(FILE_NAMES[0]) 
    ireland_gdf = read_dataframe(FILE_NAMES[1])

    # Read in the border polygon
    coastline = read_dataframe("../../national-canal-cooling-data/plotting/borders.shp")

    # Reproject to BNG
    ireland_gdf = ireland_gdf.to_crs(britain_gdf.crs)
    coastline = coastline.to_crs(britain_gdf.crs)

    # Merge
    merged_df = concat([britain_gdf, ireland_gdf])

    # Filter to canals with areas > 500 m^2 and < 200_000 m^2
    merged_df = merged_df[merged_df['area'].between(500, 200_000)]

    # Produce cluster-statistics (total area)
    stats = merged_df.groupby('region').agg({'area':'sum'})
    
    # Dissolve gdf by region and return centroid
    dissolved_gdf = merged_df.dissolve(by='region').centroid

    # Produce merged df for plotting
    plot_gdf = GeoDataFrame(stats, geometry = dissolved_gdf, crs = britain_gdf.crs)

    # Normalise area [0-1]
    plot_gdf['normalised_area'] = (plot_gdf['area']-plot_gdf['area'].min())/(plot_gdf['area'].max()-plot_gdf['area'].min())

    # Convert to km^2 for readability
    plot_gdf['area_km'] = plot_gdf['area'] * 0.000001

    # Extract the n largest clusters
    cities = plot_gdf.nlargest(6, 'area')

    # Sort by area
    plot_gdf = plot_gdf.sort_values(by=['area'])

    # Update plot params
    rcParams.update({'font.size': 9,
                     "mathtext.fontset" : "dejavuserif",
                     'mathtext.default' : 'regular'})

    # Set up output image
    fig = figure(layout='compressed', figsize=(7, 5))
    gs = fig.add_gridspec(2,1, height_ratios = [0.3, 0.7])
    
    # Add axes
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[0, 0])

    # Colours and labels
    region_colours = ["#754501", "#B76C00", "#FF9700", "#FEAD37", "#FFC26A", "#FFDEAD", "#CBCBCB"]
    region_labels = ['birmingham', 'london', 'liverpool', 'manchester', 'glasgow', 'dublin', 'other']

    # First create a colour column, where all are #CBCBCB
    plot_gdf['color'] = "#CBCBCB"
    plot_gdf['edge_color'] = "#df8344"

    # Iterate and update colours for other regions
    for name, colour in zip(region_labels, region_colours):
        plot_gdf.loc[(plot_gdf.index == name), 'color'] = colour

    # Add region labels
    for name in region_labels:
        plot_gdf.loc[(plot_gdf.index == name), 'edge_color'] = "#242424"

    # List of largest areas, as proportion of total
    urban = cities['area'].tolist()
    urban_prop = [x / sum(plot_gdf['area']) for x in urban]
    urban_prop.append(1 - sum(urban_prop))

    # Waffle chart
    Waffle.make_waffle(
        ax=ax2, 
        rows=5, 
        columns=10, 
        values=urban_prop,
        colors=["#754501", "#B76C00", "#FF9700", "#FEAD37", "#FFC26A", "#FFDEAD", "#CBCBCB"],
        legend={'labels': ['Birmingham', 'London', 'Liverpool', 'Manchester', 'Glasgow', 'Dublin', 'other'],
        'loc' : 'upper center',
        'bbox_to_anchor' : (0.5, 1.28),
        'ncol': 4, 
        'columnspacing' : 0.6,
        'labelspacing' : 0.1,
        'fontsize': 7,
        'framealpha': 0,
        'handlelength' : 0.8}
    )

    # Add y-label
    ax2.text(-0.04, 0,
                s="Area (%)",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=9,
                transform = ax2.transAxes)
    
    # Background coastline (figure ground separation)
    coastline.plot(
        ax = ax1,
        color = None,
        edgecolor = '#b9b9b9',
        linewidth = 1.5,
        ) 

    # Plot the coastline
    coastline.plot(
        ax = ax1,
        color = '#ECECEC',
        edgecolor = '#5E5E5E',
        linewidth = 0.5,
        )   

    # Plot the centroid locations, sized by area
    ax1.scatter(plot_gdf.geometry.x, plot_gdf.geometry.y, 
                     c = "#df8344",
                     s = plot_gdf['area_km'] * 50,
                     edgecolors = plot_gdf['edge_color'],
                     linewidth=0.5)
    
    # Manual legend
    for area in [0.5, 1.5, 2.5]:
           ax1.scatter([], [], c='#df8344', edgecolors = "#df8344", s= area * 50,
                       label= str(area) + ' km$^2$')

    # Add to plot
    ax1.legend(scatterpoints=1, frameon=False,
                  labelspacing=0.75,
                  loc = 'lower left', 
                  fontsize="7")
    
    # Extracts bounds of all canals
    xmin, ymin, xmax, ymax = plot_gdf.geometry.total_bounds

    # Buffer distance (m)
    buffer = 50_000

    # Set axis limits
    ax1.set_xlim(xmin - buffer, xmax + buffer)
    ax1.set_ylim(ymin - buffer, ymax + buffer)

    # Tick params
    ax1.tick_params(axis='x', labelsize=8.5)
    ax1.tick_params(axis='y', labelsize=8.5)
    ax1.tick_params(bottom=False, top=False, left=False, right=False, direction="in")

    # Remove axis labels
    ax1.xaxis.label.set_visible(False)
    ax1.yaxis.label.set_visible(False)

    # Turn off tick labels
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])

    # Add text annotations
    ax1.text(-0.05, 1.05, f"B", transform = ax1.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    ax2.text(-0.05, 1.1, f"A", transform = ax2.transAxes, weight='bold', fontsize = 13, va='top', ha='right')

    # Iterate through the n largest areas
    for geom, label, x_offset, y_offset in zip(
                    cities['geometry'], #----------------------------------------------------------- City region centroids
                    ['Birmingham', 'London', 'Liverpool', 'Manchester', 'Glasgow', 'Dublin'], #----- City labels
                    [0, 0, -117_000, 95_000, 30_000, -30_000], #------------------------------------ Label offset in x-direction (m)
                    [-30_000, -30_000, 14_000, -20_000, -25_000, -25_000]): #----------------------- Label offset in y-direction (y)
            
        # Add labels to map (+ xy offsets)
        ax1.text(geom.x + x_offset, geom.y + y_offset, label, fontsize = 7, va='center', ha='center')

    # Scale bar
    ax1.add_artist(ScaleBar(dx=1, units="m", location="upper right", length_fraction=0.2, font_properties={"size": 8}))

    # North arrow
    x, y, arrow_length = 0.05, 0.99, 0.11
    ax1.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=2, headwidth=8),
	ha='center', va='center', fontsize=12, xycoords=ax1.transAxes)

    # Save to file
    #show()
    savefig(f'../images/figure-1.png', bbox_inches='tight', dpi = 300)

def draw_figure_2(canal_id, shading_boolean, month):
    '''
    > Draw Figure 2, which includes subplots showing:
        - (A) the modelled reference and water temperatures for 2022 for the chosen feature
        - (B) a map showcasing the seleceted canal feature, as well as nearby canals, buildings and buffer areas
        - (C) the modelled reference and water temperatures for a specific month
    > Figure design has been adapted for a canal feature on the Regents Canal, London ('rc2_23') although will work for other inputs
    > Parameters:     
        - 'canal_id' (String): Feature name to plot e.g., 'rc2_23'
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
        - 'month' (Numeric): Month to highlight [1-12]
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Read in the canals dataframe for Britain and extract selected feature
    canals_gdf = read_dataframe(FILE_NAMES[0])
    filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]

    # Code ID is not present
    if len(filtered_canals) == 0:

        # Load dataset for Ireland
        canals_gdf = read_dataframe(FILE_NAMES[1])
        filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]

        # Code ID still not present (LBYL)
        if len(filtered_canals) == 0:

            # Exit the function
            print(f"The chosen canal ID ({canal_id}) is not present, please check input. Exiting function...")
            return

    # Buffer canal network based on empirically derived cooling distance (see model_urban_cooling() function), dissolve
    buffered_canal = GeoDataFrame(geometry=canals_gdf.geometry.buffer(1.5126 * canals_gdf.width + 6.7227)).dissolve()

    # Return buildings within a specified distance (500 m)
    for feature_tuple in filtered_canals.itertuples():
        buildings = return_buildings_polygon(feature_tuple, 500, plot = False, epsg = canals_gdf.crs)

    # Try and load model output record, including | excluding shading
    try:
        with open(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{canal_id}.json") as results_path:
            canal_output = load(results_path)

    # File not found
    except FileNotFoundError:
        print(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{canal_id}.json does NOT exist")

    # Extract lists of values to plot, converting to celcius
    water_temps = [x['surface_water_k'] - 273.15 for x in canal_output.values()]
    reference_temps = [x['surface_reference_k'] - 273.15 for x in canal_output.values()]

    # List of dict keys (unix), convert to datetime format
    dt = list(canal_output.keys())
    dt_datetime = [datetime.utcfromtimestamp(int(x)) for x in dt]

    # Return the indices, corresponding to the chosen month
    idx = [n for n, x in enumerate(dt_datetime) if x.month == month]

    # Filter to indexes
    filter_datetime = [dt_datetime[i] for i in idx]
    filter_water = [water_temps[i] for i in idx]
    filter_reference = [reference_temps[i] for i in idx] 

    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig = figure(layout='compressed', figsize=(7, 5))
    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    # Plot the water and reference temperatures
    ax1.plot(dt_datetime, reference_temps, label = 'Reference', color = "#FFB16F")
    ax1.plot(dt_datetime, water_temps, label ='Water', color = "#707070")
    ax1.set_xlabel("Date")
    ax1.set_ylabel("Temperature (°C)", labelpad=2)
    ax1.xaxis.tick_top()
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in")
    ax1.xaxis.set_label_position('top') 
    ax1.legend(frameon=False)

    # Patch dimensions
    patch_min = min(min(filter_reference), min(filter_water))
    patch_max = max(max(filter_reference), max(filter_water))
    temp_buffer = 2

    # Add rectangle, showing the dimensions of ax3
    ax1.add_patch(Rectangle((filter_datetime[0],patch_min - temp_buffer), filter_datetime[-1] - filter_datetime[0], 
                            (patch_max + temp_buffer) - (patch_min - temp_buffer),
                            fill = False, color=None, # alpha=0.5,
                            edgecolor="#4D4D4D", lw = 1, zorder=2, linestyle = '--'))

    # Plot the buffer area, used for air temperature modelling
    buffered_canal.plot(
        ax = ax2,
        color = "#dde3f1", # "none",
        edgecolor = '#249AC3',
        linewidth = 1,
        linestyle = 'dashed', 
        )

    # Plot the building geometries
    buildings.plot(
        ax = ax2,
        color = "#6D6D6D",
        edgecolor = '#343434',
        linewidth = 0,
        )
            
    # Plot other canal geometries
    canals_gdf.plot(
        ax = ax2,
        color = "#C8F1FF",
        edgecolor = '#249AC3',
        linewidth = 0.5,
        )   
        
    # Plot the selected canal feature
    filtered_canals.plot(
        ax = ax2,
        color = "#C8F1FF",
        edgecolor = '#249AC3',
        linewidth = 1,
        legend = True
        )
    
    # Add custom legends
    canal_polygon = mpatches.Patch(facecolor='#C8F1FF', edgecolor = '#249AC3', linewidth = 1, label='Modelled canal')
    building_polygon = mpatches.Patch(facecolor = "#6D6D6D", edgecolor = '#343434', linewidth = 0, label='Buildings')
    buffer_polygon = mpatches.Patch(facecolor="#dde3f1", edgecolor = '#249AC3', linewidth = 1, linestyle='dashed', label='Buffer area')
    
    # Top legend
    ax2_legend_top = ax2.legend(handles=[canal_polygon, building_polygon], loc = 'upper left',ncol=2)
    ax2_legend_top.get_frame().set_linewidth(0.0)
    ax2.add_artist(ax2_legend_top)

    # Bottom legend
    ax2_legend_bottom = ax2.legend(handles=[buffer_polygon], loc = 'lower right',ncol=1)
    ax2_legend_bottom.get_frame().set_linewidth(0.0)
    ax2.add_artist(ax2_legend_bottom)

    # Scale bar
    ax2.add_artist(ScaleBar(dx=1, units="m", location="lower left", length_fraction=0.2, font_properties={"size": 8},
                            label_loc = "top", box_alpha = 0.9, color = None, frameon = True,
                            scale_loc="top", sep = 2, border_pad = 0.4))
    
    # Extracts bounds of canal feature
    xmin, ymin, xmax, ymax = filtered_canals.geometry.iloc[0].bounds

    # Set desired plot dimension and ratio 
    axis_length = 500
    desired_xy_ratio = 2.5

    # Calculate buffer distances (m)
    x_buffer = ((axis_length - (xmax - xmin)) / 2) * desired_xy_ratio
    y_buffer = (axis_length - (ymax - ymin)) / 2

    # Set axis limits and update plot design
    ax2.set_xlim([filtered_canals.geometry.iloc[0].bounds[0] - x_buffer, filtered_canals.geometry.iloc[0].bounds[2] + x_buffer])
    ax2.set_ylim([filtered_canals.geometry.iloc[0].bounds[1] - y_buffer, filtered_canals.geometry.iloc[0].bounds[3] + y_buffer])
   
    # Remove axis labels
    ax2.xaxis.label.set_visible(False)
    ax2.yaxis.label.set_visible(False)

    # Turn off tick labels
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])

    # Remove ticks
    ax2.tick_params(bottom=False, top=False, left=False, right=False)

    # Add north arrow
    x, y, arrow_length = 0.06, 0.52, 0.21
    ax2.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
	arrowprops=dict(facecolor='black', width=2, headwidth=8),
	ha='center', va='bottom', fontsize=12, xycoords=ax2.transAxes)
        
    # Plot the monthly water and reference temperature record
    ax3.plot(filter_datetime, filter_reference, color = "#FFAD28", label = "Reference (°C)")
    ax3.plot(filter_datetime, filter_water, color = '#707070', label = "Water (°C)")   
    ax3.set_xlabel("Day")
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('-%d'))
    ax3.set_ylabel("Temperature (°C)")
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in")

    # Labelling
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    # Add plot annotations
    ax3.text(0.04, 0.89, f"{months[month-1]} 2022", transform = ax3.transAxes)
    ax1.text(-0.02, 1.08, f"A", transform = ax1.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    ax2.text(-0.05, 1.05, f"B", transform = ax2.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    ax3.text(-0.09, 1.05, f"C", transform = ax3.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    
    # Overview label
    ax1.text(filter_datetime[0]+timedelta(days = -3), patch_max + temp_buffer, f"C", weight='bold', fontsize = 9, va='top', ha='right')

    # Save to file
    #show()
    savefig(f'../images/figure-2-{REFERENCE_MATERIAL}-{canal_id}.png', bbox_inches='tight', dpi = 300)

def draw_figure_3(shading_boolean, canal_list, depth):
    '''
    > Model validation plot, using data from the Canal and River Trust (https://canalrivertrust.org.uk/)
    > Parameters:
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
        - 'canal_list' (List): List of canal IDs.
        - 'depth' (Numeric): Water depth (cm) which refers to the upper value:
            - 20 = 20 - 40 cm
            - 40 = 40 - 60 cm
            - 60 = 60 - 80 cm
            - 80 = 80 - 100 cm
    > Of the monitoring locations, the following are suitable for comparison:
        - Bridge Pagefield Pipe Crossing (lalc_73) --------- 10 minute frequency
        - Bridge 15 Chain Lane Bridge (nabc_2) ------------- 10 minute frequency
        - New Road Bridge 148 (suc_41) --------------------- 1 hour frequency
        - Bridge 32. Priorswood (batc_14) ------------------ 15 minute frequency (from 03-01-2022)
        - Bridge 4 Wakefield Road Bridge (hc_53) ----------- 1 hour frequency
        - Bridge 73 Anchor Bridge (cc3_11) ----------------- 10 minute frequency
    > The following locations are less suitable for comparison:
        - Anderton Waste Weir Flow (None) ------------------ not directly modelled
        - Bridge 46 Benthouse Bridge (None) ---------------- not directly modelled
        - Hawkesbury Lock (None) --------------------------- not directly modelled
        - Bridge C Lightbody Street (lalc_4) --------------- missing measurement data
        - Bridge 40 Taverners Bridge (cc3_53) -------------- missing measurement data
    ''' 

    # Dictionary of CRT wT monitoring locations and corresponding canal ids
    location_dict = {'tamc_19' : 'Anderton Waste Weir Flow Water Temperature C (All Data) (Value)', 
                     'lalc_73' : 'Bridge Pagefield Pipe Crossing Mean Temperature C (All Data) (Value)',
                     'nabc_2' : 'Bridge 15 Chain Lane Bridge Mean Temperature C (All Data) (Value)',
                     'suc_41' : 'New Road Bridge 148 Mean Temperature C (All Data) (Value)',
                     'rc_103' : 'Bridge 46 Benthouse Bridge Mean Temperature C (All Data) (Value)',
                     'cc3_60' : 'Hawkesbury Lock Mean Temperature C (All Data) (Value)',
                     'batc_14' : 'Bridge 32. Priorswood Mean Temperature C (All Data) (Value)', 
                     'hc_53' : 'Bridge 4 Wakefield Road Bridge Mean Temperature C (All Data) (Value)',
                     'lalc_4' : 'Bridge C Lightbody Street Mean Temperature C (All Data) (Value)',
                     'cc3_11' : 'Bridge 73 Anchor Bridge Mean Temperature C (All Data) (Value)',
                     'cc3_53' : 'Bridge 40 Taverners Bridge Mean Temperature C (All Data) (Value)'}
    
    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Open csv
    measured_values = read_csv("../data/model_validation.csv")
    
    # Update plot params
    rcParams.update({'font.size': 9,
                    "mathtext.fontset" : "dejavuserif",
                    'mathtext.default' : 'regular'})
    
    # Set up output image
    fig = figure(layout='compressed', figsize=(8, 6))
    gs = fig.add_gridspec(4,3, height_ratios=[2,1,2,1])

    # List of axis locations
    axes_loc = [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [3, 1], [0, 2], [1, 2], [2, 2], [3, 2]]
    axes_labels = ['A', None,'D', None,'B', None,'E', None,'C', None,'F', None]

    # Iterate through monitoring locations
    for index, canal_id in enumerate(canal_list):

        # Run correlation, using Spearman
        method, r, p, array_length = time_series_correlation(shading_boolean, canal_id, depth, method="Spearman")

        # Return plot locations
        first_plot = axes_loc[2 * index]
        second_plot = axes_loc[2 * index + 1]

        # Add axes
        ax1 = fig.add_subplot(gs[first_plot[0], first_plot[1]])
        ax2 = fig.add_subplot(gs[second_plot[0], second_plot[1]])

        # Return value, corresponding to selected id
        location = location_dict[canal_id]
        
        # Extract measured temperatures and datetimes
        measured_water = measured_values[['dt', location]]
        measured_water = measured_water.dropna()

        # Convert from string to datetime
        dt = [datetime.strptime(x, '%d/%m/%Y %H:%M') for x in measured_water['dt']]

        # Convert to unix 
        dt_unix = [str(int(x.timestamp())) for x in dt]

        # Measured dictionary
        measured_dict = dict(zip(dt_unix, measured_water[location]))

        # Load model output record, including | excluding shading
        try:
            with open(f"../results/{REFERENCE_MATERIAL}/{folder_path}-shading/model-output-{canal_id}.json") as results_path:
                canal_output = load(results_path)

        # Modelled results missing
        except FileNotFoundError:
            print("No modelled data to compare. Exiting function...")
            return

        # Extract lists of values to plot, converting to celcius
        try:
            modelled_water = [x[f'depth_water_{depth}'] - 273.15 for x in canal_output.values()]

        # EAFP: surface layer
        except KeyError:
            modelled_water = [x[f'surface_water_k'] - 273.15 for x in canal_output.values()]

        # List of dict keys (unix), convert to datetime format
        modelled_dt = list(canal_output.keys())
        modelled_datetime = [datetime.utcfromtimestamp(int(x)) for x in modelled_dt]

        # Init dictionary to store residuals
        residual_dict = {}

        # Iterate through the modelled unix times
        for key, value in canal_output.items():

            # If this unix time exists in the measured dictionary
            if key in measured_dict:
                
                # Calculate the residual, converting from K to celcius
                try: 
                    res = (value[f'depth_water_{depth}'] - 273.15) - measured_dict[key]

                # EAFP: surface layer
                except KeyError:
                    res = (value[f'surface_water_k'] - 273.15) - measured_dict[key]

                # Add to residual dictionary
                residual_dict[key] = res

        # Extract datetimes and values
        residual_dt = list(residual_dict.keys())
        residual_datetime = [datetime.utcfromtimestamp(int(x)) for x in residual_dt]
        residual_temp = [x for x in residual_dict.values()]

        # Absolute residuals
        abs_residual = [abs(x) for x in residual_temp]
        print(f"for {canal_id}, the absolute median residual is {median(abs_residual):.2f}, with p25 of {quantile(abs_residual, 0.25):.2f} and p75 of {quantile(abs_residual, 0.75):.2f}")
        print(f"for {canal_id}, the relative median residual is {median(residual_temp):.2f}, with p25 of {quantile(residual_temp, 0.25):.2f} and p75 of {quantile(residual_temp, 0.75):.2f}")

        # Plot the measured-modelled temperatures
        ax1.plot(modelled_datetime, modelled_water, color = '#E28600', label = "Modelled", linewidth = 0.9)
        ax1.plot(dt, measured_water[location], color = '#707070', label = "Measured", linewidth = 0.9)

        # Plot the residuals
        ax2.plot(residual_datetime, residual_temp, color = '#707070', linewidth = 0.9)
        ax2.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)

        # Manually set axis limits
        ax1.set_ylim(-6, 28)
        ax2.set_ylim(-9, 8)

        # Add months only
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('-%m'))
        
        # Change p to string if < 0.01
        if p < 0.01:
            p = "< 0.01"

            # Add correlation statistics
            ax1.text(0.04, 0.05, f"{method} $r_s$ = {r:.2f}\np value = {p}", transform = ax1.transAxes, weight='normal', fontsize = 6.5, va='bottom', ha='left', style='italic')
        
        # Two significant figures
        else:
            ax1.text(0.04, 0.05, f"{method} $r_s$ = {r:.2f}\np value = {p:.2f}", transform = ax1.transAxes, weight='normal', fontsize = 6.5, va='bottom', ha='left', style='italic')

        # Add location and count
        ax1.text(0.95, 0.95, f"Location = {canal_id}", transform = ax1.transAxes, weight='normal', fontsize = 6, va='top', ha='right', 
                 style='italic', color = "#727573", linespacing = 1.5)
        
        # Add absolute residuals
        ax2.text(0.04, 0.07, f"Absolute median residual = {median(abs_residual):.2f} ({quantile(abs_residual, 0.25):.2f}, {quantile(abs_residual, 0.75):.2f})", 
                 transform = ax2.transAxes, weight='normal', fontsize = 6.5, va='bottom', ha='left', style='italic')

        # Modifications based on axis index position
        if index == 0:

            # Axis labelling
            ax1.set_ylabel("Temperature (°C)", labelpad=0)
            ax2.set_ylabel("Residual (°C)", labelpad=0)
            ax2.set_xlabel(None)

            # Tick params
            ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=True)
            ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=True)

        elif index == 1:

            # Axis labelling
            ax1.set_ylabel("Temperature (°C)", labelpad=0)
            ax2.set_ylabel("Residual (°C)", labelpad=0)
            ax2.set_xlabel("Date")

            # Tick params
            ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=True)
            ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=True, labelleft=True)

        elif index in [3,5]:

            # Axis labelling
            ax1.set_ylabel(None)
            ax2.set_ylabel(None)
            ax2.set_xlabel("Date")

            # Tick params
            ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=False)
            ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=True, labelleft=False)

        else: 

            # Axis labelling
            ax1.set_ylabel(None)
            ax2.set_ylabel(None)
            ax2.set_xlabel(None)

            # Tick params       
            ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=False)
            ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False, labelleft=False)

        # Subplot labelling
        ax1.text(0.04, 0.90, f"{axes_labels[2 * index]}", transform = ax1.transAxes, weight='bold', fontsize = 13, va='center', ha='left')

    # Save to file
    # show()
    savefig(f'../images/figure-3-{REFERENCE_MATERIAL}-{depth}-{depth+20}cm.png', bbox_inches='tight', dpi = 300)

def draw_figure_4(start_date, interval, shading_boolean):
    '''
    > Draw Figure 4 which summarises the temporal results for daytime and nighttime, showcasing:
        - (A-C) fluctuations in median air temperature change at weekly, fortnightly, or monthly intervals
        - (B-D) histograms of median air temperature change
    > Parameters:
        - 'start_date' (String): Start date for plotting ("2022-01-01 10:00:00")
        - 'interval' (String): Interval for aggregration ("week", "fortnight", "month")
        - 'shading_boolean' (Boolean): Include (True) or exclude shading effects (False)
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Load summary file, including | excluding shading
    try:
        with open(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json") as path:
                summary_data = load(path)

    # File not found, exit function
    except FileNotFoundError:
        print(f"../outputs/{REFERENCE_MATERIAL}/region-values-{folder_path}-shading.json not found. Exiting function...")
        return

    # Determine averaging interval ('fortnight')
    if interval == "fortnight":
        day_medians = [x['fortnightly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['fortnightly values']['night medians'] for x in summary_data.values()]

    # Extract weekly medians
    elif interval == "week":
        day_medians = [x['weekly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['weekly values']['night medians'] for x in summary_data.values()]

    # Return monthly medians,
    elif interval == "month":
        day_medians = [x['monthly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['monthly values']['night medians'] for x in summary_data.values()]

    # Invalid string
    else:
        print(f"Invalid 'interval' parameter ({interval}). Exiting function...")
        return
    
    # Extract areas for weighting
    areas = array([x['geometry']['area'] for x in summary_data.values()])  

    # Transpose data from by region to by interval
    day_medians = list(map(list, zip(*day_medians)))
    night_medians = list(map(list, zip(*night_medians)))

    # Ignore the final entry
    if interval == "fortnight" or interval == "week":
        day_medians = day_medians[:-1]
        night_medians = night_medians[:-1]

    # Localise niave datetime to aware datetime (UK) and then convert to UTC
    local_datetime = to_datetime(start_date)\
        .tz_localize('Europe/London')\
        .tz_convert('UTC')
    
    # Create list of datetimes for x-axis (month currently unused)
    if interval == "fortnight":
        dt = [local_datetime + timedelta(days = 7 * (2 * n - 1)) for n in range(1, len(day_medians)+1)]
    elif interval == "week":
        dt = [local_datetime + timedelta(days = 3 + (7 * n)) for n in range(0, len(day_medians))]
    elif interval == "month":
        dt = [local_datetime + relativedelta(months=n) for n in range(0, len(day_medians))]

    # Flatten list of aggregated values
    flat_day = [num for sublist in day_medians for num in sublist]
    flat_night = [num for sublist in night_medians for num in sublist]

    # Set global font size
    rcParams.update({'font.size': 9})
    rcParams.update({"mathtext.fontset" : "dejavuserif"})

    # Set up output image
    fig = figure(layout='compressed', figsize=(6, 4))
    gs = fig.add_gridspec(2,3)

    # Histograms
    ax1 = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[1, 2])

    # Boxplots
    ax3 = fig.add_subplot(gs[0, 0:2])
    ax4 = fig.add_subplot(gs[1, 0:2])

    # Iterate through inputs (axes, data, labels, colours)
    for right_ax, left_ax, ax_label, data, flat, label, fill, edge, time_colour in zip(
                                [ax1, ax2], [ax3, ax4], #----------------------------------------------------------------- Matplotlib axes
                                [["B", "A"], ["D", "C"]], #--------------------------------------------------------------- Subplot labels
                                [day_medians, night_medians], #----------------------------------------------------------- Data by aggregation period
                                [flat_day, flat_night], #----------------------------------------------------------------- Flattened data
                                ['Day', 'Night'], #----------------------------------------------------------------------- Labels
                                ["#FFB16F", "#707070"], #------------------------------------------------------------- Fill colour (histogram)
                                ["#707070", "#232323"], #------------------------------------------------------------- Edge colour (histogram)
                                [["#FFD0A6", "#FFA453", "#B05300"],["#B9B9B9", "#656565", "#090909"]]): #---- Colours for shaded quantiles


        # Quantiles (0.25, 0.5, 0.75), no weighting
        #q1 = [quantile(x, 0.25) for x in data]
        #q2 =  [quantile(x, 0.5) for x in data]
        #q3 = [quantile(x, 0.75) for x in data]

        # Quantiles (0.25, 0.5, 0.75), weighted by area
        q1 = [weighted_median(x, areas, 0.25) for x in data]
        q2 = [weighted_median(x, areas, 0.5) for x in data]
        q3 = [weighted_median(x, areas, 0.75) for x in data]
        
        # Tukey whiskers (1.5 * IQR)
        tukey_lower = [a - ((c - a) * 1.5) for a,c in zip(q1, q3)]
        tukey_upper = [c + ((c - a) * 1.5) for a,c in zip(q1, q3)]

        # Init lists
        whisker_lower = []
        whisker_upper = []

        # Iterate, and find the last datum less than tukey upper (and v.v. for tukey lower)
        for val, low_limit, upp_limit in zip(data, tukey_lower, tukey_upper):
            whisker_lower.append(min([i for i in val if i > low_limit]))
            whisker_upper.append(max([i for i in val if i < upp_limit]))
        
        # Set bin min-max, with 0.5°C buffer
        lower_bins = min(whisker_lower) - 0.5
        upper_bins = max(whisker_upper) + 0.5

        # Y axis limits, 1°C buffer
        axis_min, axis_max = lower_bins - 0.5, upper_bins + 0.5

        # Plot histogram
        right_ax.hist(flat, bins=arange(lower_bins, upper_bins + 0.1, 0.1), orientation="horizontal", histtype='stepfilled',
                label = f"{label} " + "$\Delta$$\it{PT}$" , 
                facecolor = fill, edgecolor = edge, linewidth = 0.5, alpha = 1)
        
        # Add common labels 
        right_ax.set_ylabel(None)
        right_ax.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=False, labelbottom=False)
        right_ax.legend(frameon=False)
        right_ax.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)
        right_ax.set_ylim(axis_min, axis_max)

        # Subplot labels
        right_ax.text(-0.05, 1.05, ax_label[0], transform = right_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

        # Boxplot style
        if interval == "month":

            # Bar dimensions, for monthly plot
            h = [up - low for up, low in zip(whisker_upper, whisker_lower)]
            hq =  [up - low for up, low in zip(q3, q1)]

            # X-axis labels
            months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

            # Horizontal line
            left_ax.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)

            # Plot boxplot
            left_ax.bar(x = months, bottom = whisker_lower, width = 1, height = h, facecolor = time_colour[0], edgecolor = None,
                    label = r'$\genfrac{}{}{0}{}{\it{p}75 + 1.5iqr}{\it{p}25 - 1.5iqr}$')
            left_ax.bar(x = months, bottom = q1, width = 1, height = hq, facecolor = time_colour[1], edgecolor = None,
                    label = "$\it{p}$25 – $\it{p}$75")
            left_ax.hlines(y = q2, xmin = arange(-0.5, 11.5, 1), xmax = arange(0.5, 12.5, 1), color = time_colour[2], 
                    label = "$\it{p}$50", linewidth = 1)
            
            # Add legend
            hab, lab = left_ax.get_legend_handles_labels()
            left_ax.legend(handles = [hab[0], hab[2], hab[1]], labels = [lab[0], lab[2], lab[1]], frameon=False, ncol = 3, fontsize="8", loc = 'lower right', labelspacing = 0.2,
                columnspacing=0.5, borderpad=0.1)  

            # Add labels to top subplot
            left_ax.set_xlabel(None)
            left_ax.set_ylabel("$\Delta$$\it{PT}$ (°C)", labelpad=0)
            left_ax.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=False)
            left_ax.set_ylim(axis_min, axis_max)

            # Subplot labels
            left_ax.text(-0.03, 1.05, ax_label[1], transform = left_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

        # Weekly and fortnightly intervals
        else:

            # Horizontal line
            left_ax.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)

            # Time-series style 
            left_ax.fill_between(dt, whisker_lower, whisker_upper, color = time_colour[0], label = r'$\genfrac{}{}{0}{}{\it{p}75 + 1.5iqr}{\it{p}25 - 1.5iqr}$')
            left_ax.fill_between(dt, q1, q3, color = time_colour[1], label = "$\it{p}$25 – $\it{p}$75")
            left_ax.plot(dt, q2, color = time_colour[2], label = "$\it{p}$50")

            # Add legend
            handles, labels = left_ax.get_legend_handles_labels()
            left_ax.legend(handles[::-1], labels[::-1], frameon=False, ncol = 3, fontsize="8", loc = 'upper right', labelspacing = 0.2,
                columnspacing=0.5, borderpad=0.2)  

            # Add labels to top subplot
            left_ax.set_ylabel("$\Delta$$\it{PT}$ (°C)")
            left_ax.set_ylim(axis_min, axis_max)

            # Set date formatter
            left_ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))

            # Axis labels
            if interval == "fortnight":
                left_ax.text(-0.08, 1.05, ax_label[1], transform = left_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
            
            # Axis labels
            elif interval == "week":
                left_ax.text(-0.08, 1.05, ax_label[1], transform = left_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

                # Plot heatwave labels, corresponding to:
                # Week 23: heatwave 15-17 June
                # Week 28: heatwave 17-19 July
                # Week 31: heatwave 9-13 Aug

                # Iterate through registered heatwave weeks
                for week in [23, 28, 31]:
                
                    # Week, Q1-Q2-Q3
                    print(f"For {dt[week]}, weekly median air temperature change ({label}) was {q2[week]:.2f} (Q1: {q1[week]:.2f}, Q3: {q3[week]:.2f})")

                # Heatwave labels for daytime
                if label == "Day":
                    for x, y, label in zip([dt[23], dt[28], dt[31]], [whisker_lower[23], whisker_lower[28], whisker_lower[31]], ['I', 'II', 'III']):
                        left_ax.text(x, y-0.6, label, fontsize = 8, va='bottom', ha='center')

    # Axes specific labelling
    ax1.set_xlabel(None)
    ax2.set_xlabel("Frequency")
    ax4.set_xlabel("Date")
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=False)
    ax4.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=True)
    
    # Remove last xtick label to simplify plot design
    if interval == "week":
        setp(ax4.get_xticklabels()[-1], visible=False)

    #show()
    savefig(f'../images/figure-4-{REFERENCE_MATERIAL}-{interval}_final.png', bbox_inches='tight', dpi = 300)

def shading_gif(feature_name, start_date, days, delete, frames=10, loops=10):
    '''
    > Creates a gif showing shadows cast by buildings for a specific canal feature, datetime and duration
    > Parameters:
        - 'feature_name' (String): Feature ID to model for
        - 'start_date' (String): Start date for plotting ("2022-06-01 10:00:00")
        - 'days' (Numeric): Duration of analysis (days)
        - 'delete' (Boolean): Delete (True) or keep (False) image files (.png) after completion 
        - 'frames' (Numeric): Number of frames per second (fps), defaults to 10
        - 'loops' (Numeric): Nuber of loops before gif end, defaults to 10
    '''

    # Duration of analysis in hours
    duration = days * 24

    # Read the dataframe for Britain
    canals_gdf = read_dataframe(FILE_NAMES[0])

    # Extract chosen canal
    modelled_canal = canals_gdf.loc[canals_gdf['code_id'] == feature_name]

    # If empty
    if modelled_canal.empty:
        
        # Load the dataset for Ireland, and extract canal
        canals_gdf = read_dataframe(FILE_NAMES[1])

        # Extract chosen canal
        modelled_canal = canals_gdf.loc[canals_gdf['code_id'] == feature_name]

        # If feature not present, exit function
        if modelled_canal.empty:
            print(f"Canal feature {feature_name} missing from both datasets. Check feature code. Exiting function...")
            return

    # Convert to series, project to 4326 and return lon-lat
    geographic_geometry = GeoSeries(modelled_canal.geometry.centroid, crs = canals_gdf.crs).to_crs(4326)
    longitude, latitude = geographic_geometry.iloc[0].x, geographic_geometry.iloc[0].y

    # Localise niave datetime to aware datetime (UK) and then convert to UTC
    local_datetime = to_datetime(start_date)\
    .tz_localize('Europe/London')\
    .tz_convert('UTC')

    # Sets end datetime (+ duration)
    end_datetime = local_datetime + timedelta(hours = duration)

    # Create empty lists and copy of datetime
    datetimes = []
    solar_altitudes = []
    local_copy = local_datetime
    filenames = []

    # Iterate, append to list, add 15 minutes
    while local_copy < end_datetime:
        datetimes.append(local_copy)
        solar_altitudes.append(None)
        local_copy += timedelta(hours=0.25)

    # Iterate
    for feature_tuple in modelled_canal.itertuples():

        # Store tuple geometry as series
        canal_object = GeoSeries(Polygon(feature_tuple.geometry), crs = canals_gdf.crs)

        # Return buildings within a specified distance (500 m) of feature polygon (using id)
        buildings = return_buildings_polygon(feature_tuple, 500, plot = False, epsg = canals_gdf.crs)

        # Keep unique ID, building height, and geometry
        building_data = buildings[['relhmax', 'geometry']]

        # Rename column to work with pybdshadow
        building_data = building_data.copy().rename(columns={'relhmax': 'height'})

        # Project to geographic coordinates
        building_data_projected = building_data.to_crs(4326)

        # Explode multi-part geometries into multiple single geometries.
        building_data_projected = building_data_projected.explode(index_parts = False)

        # Add 'building_id' column
        building_data_projected['building_id'] = range(1, len(building_data_projected) + 1)

        # Add custom legends
        canal_polygon = mpatches.Patch(facecolor='#C8F1FF', edgecolor = '#249AC3', linewidth = 1, label='Canals')
        building_polygon = mpatches.Patch(facecolor = "#474747", edgecolor = None, linewidth = 1, label='Buildings')
        shadow_polygon = mpatches.Patch(facecolor = "#B3B3B3", edgecolor = None, linewidth = 0.5, label='Shadows')

        # Continue until we have modelled the required duration
        while local_datetime < end_datetime:
            
            # Filename
            filename = f'../images/shadows/shadows-for-{feature_tuple.code_id}-at-{local_datetime.date().strftime("%m_%d")}-{local_datetime.time().strftime("%H_%M_%S")}.png'

            # If the file exists, append to list and skip
            if isfile(filename):   
                filenames.append(filename)
                local_datetime += timedelta(hours=0.25)

            # Add altitude using dt index, calculate the solar altitude (in degrees) using PySolar
            else:             
                try:
                    solar_altitudes[datetimes.index(local_datetime)] = get_altitude(latitude, longitude, local_datetime)

                # Missing final datetime for list
                except ValueError:
                    
                    # Append missing
                    datetimes.append(local_datetime)
                    solar_altitudes.append(get_altitude(latitude, longitude, local_datetime))

                # Set up output image
                fig, (my_ax, my_ax2) = subplots(2, 1, figsize=(9, 10), height_ratios=[0.8,0.2])

                # Set axis limits
                my_ax.set_xlim([buildings.total_bounds[0],  buildings.total_bounds[2]])
                my_ax.set_ylim([buildings.total_bounds[1],  buildings.total_bounds[3]])

                # Plot the extracted building geometries
                buildings.plot(
                    ax = my_ax,
                    color = "#474747",
                    edgecolor = None,
                    linewidth = 0.5,
                    )

                # Plot the canal geometries
                canals_gdf.plot(
                    ax = my_ax,
                    color = "#C8F1FF",
                    edgecolor = '#249AC3',
                    linewidth = 0.5,
                    )   
                
                # Plot the modelled canal geometry
                canal_object.plot(
                    ax = my_ax,
                    color = "#C8F1FF",
                    edgecolor = '#249AC3',
                    linewidth = 0.5,
                    )


                # Adapted to use PySolar
                try:
                    shadows = bdshadow_sunlight(building_data_projected, local_datetime)

                # Before or after sunrise
                except ValueError:
                    print("Before sunrise or after sunset")
                
                # Valid time
                else:

                    # Assigns a CRS to the output (4326)
                    shadows.crs = building_data_projected.crs

                    # Reproject to British National Grid
                    shadows = shadows.to_crs(building_data.crs)

                    # Plot the shadow geometries
                    shadows.plot(
                        ax = my_ax,
                        color = "#B3B3B3",
                        edgecolor = None,
                        linewidth = 0.5,
                        alpha = 0.5
                        )
                    
  
                # Add a custom legend
                ax_legend_top = my_ax.legend(handles=[canal_polygon, building_polygon, shadow_polygon], loc = 'upper left',ncol=1)        
                my_ax.add_artist(ax_legend_top)    

                # Add north arrow
                x, y, arrow_length = 0.97, 0.99, 0.1
                my_ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                    arrowprops=dict(facecolor='black', width=5, headwidth=15),
                    ha='center', va='center', fontsize=20, xycoords=my_ax.transAxes)
                
                # Add horizontal line
                my_ax2.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)
                
                # Plot the solar altitudes
                my_ax2.plot(datetimes, solar_altitudes, label = 'Solar altitude', color = "#FFB16F")

                # Set the axis limits
                my_ax2.set_xlim([datetimes[0],  datetimes[-1]])
                my_ax2.set_ylim([-30,  70])

                # Update to show times only
                my_ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

                # Axis labels
                my_ax2.set_xlabel("Time")
                my_ax2.set_ylabel("Solar altitude (°)")
                
                # Tight layout
                fig.tight_layout()

                # Save the result
                filenames.append(filename)
                savefig(filename, bbox_inches='tight')

                # Close the figure
                close()

                # Skip to next 15 minutes
                local_datetime += timedelta(hours=0.25)

    # Build gif
    with imageio.get_writer(f'../images/shadows/shadows-for-{feature_tuple.code_id}.gif', mode='I', fps = frames, loop = loops) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

            # Delete file from directory
            if delete:
                remove(filename)

#==================== Supplementary ====================#

def draw_supplementary_figure_A(canal_id):
    '''
    > Function to illustrate the outputs of modelspin up
    > Parameters:
        - 'canal_id' (String): Feature ID to plot. Designed for a feature of the Regents Canal, London (rc2_23) but will work with other inputs (requires model spin-up to be performed)
    '''

    # If the model spin-up file exists
    if isfile(f"../outputs/spin_up_output_{canal_id}.csv"):

        # Load using np
        arr = loadtxt(f"../outputs/spin_up_output_{canal_id}.csv",
                        delimiter=",")
        
        # Extract each column
        water = arr[:, 0]
        energy = arr[:, 1]
        reference = arr[:, 2]

    # Perform model-spin up
    else:

        # Read in the canals dataframe for Britain
        canals_gdf = read_dataframe(FILE_NAMES[0])

        # If the canal exists
        if canal_id in canals_gdf['code_id'].values:

            # Filter to selected
            filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]
        
        # Read in the canals dataframe for Ireland
        else:
            canals_gdf = read_dataframe(FILE_NAMES[1])

            # If the canal exists
            if canal_id in canals_gdf['code_id'].values:

                # Filter to selected
                filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]

            # Exit function
            else:
                print(f"Canal {canal_id} is MISSING from both datasets. Exiting function...")
                return

        # Transformer for projected to geographic coordinates (BNG | ING > WGS84)
        transformer = Transformer.from_crs(canals_gdf.crs, 4326)

        # Feature centroid and area (shapely)
        latitude, longitude = transformer.transform(filtered_canals.geometry.iloc[0].centroid.x, filtered_canals.geometry.iloc[0].centroid.y)

        # Canal area (m^2)
        feature_area = filtered_canals.geometry.iloc[0].area

        # Canal depth (m)
        canal_depth = filtered_canals.iloc[0]['draught']

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

        # Init list of starting value for reference (K = 282), comprising 10 depth layers [0 - .5m]
        reference_temp = [MEAN_KELVIN] * 10

        # Model spin up (30 iterations), based on a composite climate record for 2021-12-15, return full record
        _, _, _, water, energy, reference = model_spin_up(
                                                water_temp, #-------------------------- Water temperature (K)
                                                water_energy, #------------------------ Water energy (J)
                                                reference_temp, #---------------------- Reference temperature (K)
                                                split_points, #------------------------ List of split points (m depth)
                                                feature_area, #------------------------ Modelled area (m^2)
                                                latitude, longitude, #----------------- Geographic coordinates
                                                "../outputs/spin_up_climate.json", #--- Path to outfile file
                                                30, #---------------------------------- Number of iterations
                                                False, #------------------------------- Save output (Boolean)
                                                canal_id) #---------------------------- Feature ID for output file
    
    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig, (ax1, ax2, ax3) = subplots(1, 3, figsize=(8, 4), layout='compressed')

    # Plot the temperature and energy values
    ax1.plot(water, color = '#599FFF')
    ax2.plot(reference, color = '#FFB16F')
    ax3.plot(energy, color = '#707070')

    # Set aspect ratio to 1
    ratio = 1.0
    x_left, x_right = ax1.get_xlim()
    y_low, y_high = ax1.get_ylim()
    ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

    # x-axis labelling
    ax1.set_xlabel("Model time-steps")
    ax2.set_xlabel("Model time-steps")
    ax3.set_xlabel("Model time-steps")

    # y-axis labelling
    ax1.set_ylabel("Temperature (K)")
    ax2.set_ylabel(None)
    ax3.set_ylabel("Energy (J)", labelpad=0)
    ax3.yaxis.offsetText.set_fontsize(7)

    # Tick params
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False)
    ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False)
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=False)

    # Axis labels
    ax1.text(-0.05, 1.05, f"A", transform = ax1.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax2.text(-0.05, 1.05, f"B", transform = ax2.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax3.text(-0.05, 1.05, f"C", transform = ax3.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

    # Save
    #show()
    savefig(f'../images/supp-figure-A.png', bbox_inches='tight', dpi = 300)

def draw_supplementary_figure_B(feature_name):
    '''
    > Figure to compare water and reference temperatures when shading is included | excluded, for a specific feature
    > Parameters:
        - 'canal_id' (String): Feature ID to plot e.g., 'rc_9' (high shading example) or 'lalc_170' (low shading example)
    '''

    # Try and load model output record, including | excluding shading
    try:
        with open(f"../results/{REFERENCE_MATERIAL}/excluding-shading/model-output-{feature_name}.json") as excluding_path:
            excluding_data = load(excluding_path)

        with open(f"../results/{REFERENCE_MATERIAL}/including-shading/model-output-{feature_name}.json") as including_path:
            including_data = load(including_path)

    # File not found
    except FileNotFoundError:
        print(f"One or more results files is MISSING for feature {feature_name}. Exiting function...")
        return

    # Extract water temperature values, excluding | including shading effects
    excluding_water = [x['surface_water_k'] - 273.15 for x in excluding_data.values()]
    including_water = [x['surface_water_k'] - 273.15 for x in including_data.values()]

     # Extract reference temperature values, excluding | including shading effects
    excluding_reference = [x['surface_reference_k'] - 273.15 for x in excluding_data.values()]
    including_reference = [x['surface_reference_k'] - 273.15 for x in including_data.values()]

    # Extract differences
    diff_water = [e - i  for e, i in zip(excluding_water, including_water)]  
    diff_reference = [e - i  for e, i in zip(excluding_reference, including_reference)]  

    # List of dict keys (unix), convert to datetime format
    dt = list(excluding_data.keys())
    dt_datetime = [datetime.utcfromtimestamp(int(x)) for x in dt]

    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig = figure(layout='compressed', figsize=(7, 5))
    gs = fig.add_gridspec(2,2, width_ratios = [0.65, 0.35]) # width_ratio = [0.7, 0.3])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1])
    ax4 = fig.add_subplot(gs[1, 1])

    # Include design common to both subplots
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)
        ax.xaxis.set_label_position('bottom') 
        ax.xaxis.tick_top()

    # Plot the water temperature including | excluding shading
    ax1.plot(dt_datetime, excluding_water, label = 'Excluding shading', color = "#FFB16F", alpha=0.75)
    ax1.plot(dt_datetime, including_water, label = 'Including shading', color = "#707070", alpha=0.75)
    ax1.set_xlabel(None)
    ax1.set_ylabel("Water temperature (°C)", labelpad=-4)
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=False, labeltop=False)
    ax1.legend(frameon=False)

    # Plot the reference temperature including | excluding shading
    ax2.plot(dt_datetime, excluding_reference, label = 'Excluding shading', color = "#FFB16F", alpha=0.75)
    ax2.plot(dt_datetime, including_reference, label = 'Including shading', color = "#707070", alpha=0.75)
    ax2.set_xlabel("Month")
    ax2.set_ylabel("Reference temperature (°C)", labelpad=-4)
    ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=True, labeltop=False)
    ax2.legend(frameon=False)

    # Plot the water temperature difference including | excluding shading
    ax3.plot(dt_datetime, diff_water, label = None, color = "#707070", alpha=0.75)
    ax3.set_xlabel(None)
    ax3.set_ylabel("Temperature difference (°C)")
    ax3.yaxis.set_label_position("right")
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=False, labelbottom=False, labeltop=False, labelright=True)

    # Plot the reference temperature difference including | excluding shading
    ax4.plot(dt_datetime, diff_reference, label = None, color = "#707070", alpha=0.75)
    ax4.set_xlabel("Month")
    ax4.set_ylabel("Temperature difference (°C)")
    ax4.yaxis.set_label_position("right")
    ax4.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=False, labelbottom=True, labeltop=False, labelright=True)

    # Axis labels          
    ax1.text(-0.05, 1.05, f"A", transform = ax1.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax2.text(-0.05, 1.05, f"B", transform = ax2.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax3.text(-0.05, 1.05, f"C", transform = ax3.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax4.text(-0.05, 1.05, f"D", transform = ax4.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

    # Update axis date formatting
    for ax in [ax2, ax4]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m'))

    # show()
    savefig(f'../images/supp-figure-B-{REFERENCE_MATERIAL}-{feature_name}.png', bbox_inches='tight', dpi = 300)


def shading_comparison():
    '''
    > Function to analyse the effects of shading on water and reference materials
    > Rather than comparing the summary for the full models (i.e., "including" vs. "excluding"), we compare the materials independently:
        - Unshaded water vs. shaded water
        - Unshaded reference vs. unshaded reference
    > This is a more valid comparison than comparing the difference between shaded water-shaded reference vs. unshaded water-unshaded reference

    '''
    # File paths for Britain and Ireland
    file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", "../data/urban-canals-ireland-dissolved-id-tm65.shp"]

    # If the one or both summary files already exist, exit the function
    if isfile(f"../outputs/{REFERENCE_MATERIAL}/feature-values-shading-water-comparison.json") or isfile(f"../outputs/{REFERENCE_MATERIAL}/feature-values-shading-water-comparison.json"):
        print(f"Output file: '../outputs/{REFERENCE_MATERIAL}/feature-values-shading-[water/reference]-comparison.json' EXISTS. Exiting function.")
        return    

    # Init output dicts
    output_water, output_reference = {}, {}

    # Iterate and read each file
    for file in file_names:
        canals_gdf = read_dataframe(file)

        # Iterate efficiently using itertuples
        for feature_tuple in canals_gdf.itertuples():
            
            # Try to load  results for this feature
            try:
                # Load results including shading
                with open(f"../results/{REFERENCE_MATERIAL}/including-shading/model-output-{feature_tuple.code_id}.json") as including_path:
                    including_data = load(including_path)
                    
                # Load results excluding shading
                with open(f"../results/{REFERENCE_MATERIAL}/excluding-shading/model-output-{feature_tuple.code_id}.json") as excluding_path:
                    excluding_data = load(excluding_path)

            # If the file is missing (e.g, < 500 m^2 | > 200_000 m^2)
            except FileNotFoundError:

                # Skip to next canal feature
                print(f"One or both files missing for {feature_tuple.code_id}")
                continue

            # If there is an error with the file
            except decoder.JSONDecodeError:
                
                # Skip to next canal feature
                print(f"Error with one or both files for {feature_tuple.code_id}")
                continue

            # Keys for water ΔT, reference ΔT (\u0394 is unicode delta)
            water_key, reference_key = "water \u0394T", "reference \u0394T"

            # Init dicts for daytime | nighttime values for each material
            water_day, water_night, reference_day, reference_night = {}, {}, {}, {}

            # Iterate through one of the data dicts
            for key, value in including_data.items():

                # Convert unix to datetime
                dat = datetime.fromtimestamp(int(key))

                # Between 10:00 and 15:45 (5h 45 m)
                if dat.hour >= DAY_START and dat.hour < DAY_END:

                    # Add to daytime dictionaries, including - excluding
                    water_day[key] = value[water_key] - excluding_data[key][water_key]
                    reference_day[key] = value[reference_key] - excluding_data[key][reference_key]

                # Between 22:00 and 03:45 (5h 45m)
                elif dat.hour >= NIGHT_START or dat.hour < NIGHT_END:

                    # Add to nighttime dictionaries, including - excluding
                    water_night[key] = value[water_key] - excluding_data[key][water_key]
                    reference_night[key] = value[reference_key] - excluding_data[key][reference_key]


            # Init lists to store week, fortnight and monthly data
            water_week_list = defaultdict(list)
            water_fortnight_list = defaultdict(list)
            water_month_list = defaultdict(list)

            # Init lists to store week, fortnight and monthly data
            reference_week_list = defaultdict(list)
            reference_fortnight_list = defaultdict(list)
            reference_month_list = defaultdict(list)

            '''
            # Init lists to store week, fortnight and monthly data
            day_week_list = defaultdict(list)
            day_fortnight_list = defaultdict(list)
            day_month_list = defaultdict(list)
            
            # As above, but for night values
            night_week_list = defaultdict(list)
            night_fortnight_list = defaultdict(list)
            night_month_list = defaultdict(list)
            '''

            '''
            # Iterate through the daytime values
            for key, value in water_day.items():

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
                water_month_list[month].append(value)
                reference_month_list[month].append(reference_day[key])

                # Ignore final week and fortnight 
                if week < 53:
                    water_week_list[week].append(value)
                    reference_week_list[week].append(reference_day[key])

                if fortnight < 27:
                    water_fortnight_list[fortnight].append(value)
                    reference_fortnight_list[fortnight].append(reference_day[key])

            '''

            '''
            # Iterate through the nighttime values
            for key, value in peak_night.items():

                # Convert unix to datetime, extract month integer
                month = datetime.fromtimestamp(int(key)).month

                # Append values to corresponding key
                night_month_list[month].append(value['temperature_difference'])

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

                # Ignore final week and fortnight 
                if week < 53:
                    night_week_list[week].append(value['temperature_difference'])
                if fortnight < 27:
                    night_fortnight_list[fortnight].append(value['temperature_difference'])
            '''

            # Extract air temperatures
            water_dt = [x for x in water_day.values()]
            reference_dt = [x for x in reference_day.values()]


            # Add to output
            output_water[feature_tuple.code_id] = {'day values' : {'median' : median(water_dt), 
                                                            'q1' : quantile(water_dt, 0.25),
                                                            'q3' : quantile(water_dt, 0.75)}}
            
            
            output_reference[feature_tuple.code_id] = {'day values' : {'median' : median(reference_dt), 
                                                            'q1' : quantile(reference_dt, 0.25),
                                                            'q3' : quantile(reference_dt, 0.75)}}

            '''
            # Add summary values to output
            output[feature_tuple.code_id] = {# Geometry attributes
                                            'geometry' : {'area' : feature_tuple.geometry.area, 
                                                        'buffer area' : buffered_area, 
                                                        'depth' : feature_tuple.draught, 
                                                        'easting' : feature_tuple.geometry.centroid.x,
                                                        'northing' : feature_tuple.geometry.centroid.y},
                                            # Median (Q1:Q3) daytime values
                                            'day values' : {'median' : median(day_temps), 
                                                            'q1' : quantile(day_temps, 0.25),
                                                            'q3' : quantile(day_temps, 0.75)},
                                            # Median (Q1:Q3) daytime energy values (J)                
                                            #'energy values' : {'median' : median(day_energy), 
                                            #                'q1' : quantile(day_energy, 0.25),
                                            #                'q3' : quantile(day_energy, 0.75)},
                                            # Median (Q1:Q3) nighttime values
                                            'night values' : {'median' : median(night_temps), 
                                                            'q1' : quantile(night_temps, 0.25),
                                                            'q3' : quantile(night_temps, 0.75)},
                                            # Week values (day and night)
                                            'weekly values' : {'day medians' : [median(day_week_list[x]) for x in day_week_list],
                                                               'night medians' : [median(night_week_list[x]) for x in night_week_list]}, 
                                            # Fortnight values (day and night)  
                                            'fortnightly values' : {'day medians' : [median(day_fortnight_list[x]) for x in day_fortnight_list],
                                                                    'night medians' : [median(night_fortnight_list[x]) for x in night_fortnight_list]},
                                            # Monthly values, daytime and nighttime (median, Q1:Q3)             
                                            'monthly values' : {'day medians' : [median(day_month_list[x]) for x in day_month_list], 
                                                            'day q1' : [quantile(day_month_list[x], 0.25) for x in day_month_list],
                                                            'day q3' : [quantile(day_month_list[x], 0.75) for x in day_month_list],
                                                            'night medians' : [median(night_month_list[x]) for x in night_month_list],
                                                            'night q1' : [quantile(night_month_list[x], 0.25) for x in night_month_list],
                                                            'night q3' : [quantile(night_month_list[x], 0.75) for x in night_month_list]}}

            '''              
            # Progress statement
            # print(f"Completed ../results/excluding-shading/model-output-{feature_tuple.code_id}.json")
   
    # If the output dictionary contains data
    if bool(output_water):

        # When all features have been summarised, save to outputs
        dump(output_water, open(f"../outputs/{REFERENCE_MATERIAL}/feature-values-shading-water-comparison.json", 'w'))
        dump(output_reference, open(f"../outputs/{REFERENCE_MATERIAL}/feature-values-shading-reference-comparison.json", 'w'))

    # No output
    else:
        print(f"The output dictionary is EMPTY. '../outputs/{REFERENCE_MATERIAL}/feature-values-shading-[water / reference]-comparison' has not been created")
        return



# The script is meant to be run
if __name__ == '__main__':

    # Select reference material
    user_input = input(f"The reference material is {REFERENCE_MATERIAL.upper()}. Do you want to continue? (yes/no): ")

    # Continue or exit
    if user_input.lower() in ["yes", "y"]:
        print("Continuing...")
    else:   
        exit("To switch the reference material, update params.py. Exiting...")
    
    # Run main function
    main()

