''' Functions to analyse and present the model results, produced by national-canal-cooling.py, which are stored in /results 
    > Summary files stored in /outputs
    > Figures stored in /images '''

# Force use of Shapely 2.0
from os import environ, remove
environ['USE_PYGEOS'] = '0'

import imageio
from scipy import stats
from copy import deepcopy
from pandas import read_csv
from pywaffle import Waffle
from datetime import datetime
from statistics import median
from json import load, decoder
from pyproj import Transformer
import matplotlib.dates as mdates
from pyogrio import read_dataframe
from collections import defaultdict
import matplotlib.patches as mpatches
from timeit import default_timer as timer
from dateutil.relativedelta import relativedelta
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.pyplot import Rectangle, subplots, savefig, rcParams, rc, figure, setp
from numpy import arange, quantile, sqrt, median as np_median, linspace, loadtxt, array, sum, sqrt, abs, ediff1d

# Import used defined functions
from auxiliary import *
from functions import *

# Set font family
rc('font', family='serif')

# Parameters for extracting peak day (10:00 - 15:45) and nighttime periods (22:00 - 03:45)
# This approach provides a consistent basis for comparison across the year (5h 45 m), rather than utilising variable sunrise / sunset times.
DAY_START, DAY_END  = 10, 16
NIGHT_START, NIGHT_END = 22, 4

def main():

    # Initialise timer
    start = timer()

    #==================== Analysis ====================#

    # [1] Extract summary values for each canal feature, including (True) or excluding (False) shading 
    # extract_summary_data(True) # 2024-03-08

    # [2] Print summary statistics to console and save to file, including (True) or excluding (False) shading 
    # print_summary_values(True) # 2024-03-08

    # [3] Analyse cooling effects during heatwaves (15–17 June, 17–19 July and 9–15 August)
    # Perform analysis for heatwave days only ('actual'), or include a longer period ('extended', +/- 2 days)
    # heatwave_analysis('actual') # 2024-03-08

    # [4] Extract summary values for each canal feature, but using unshaded concrete as the baseline for comparison
    # Analysis on a subset of highly shaded features, so requires output of extract_highly_shaded() from auxiliary.py
    # shading_analysis()

    # [5] Prints summary values produced by shading_analysis()
    # print_shading_summary()

    # [6] Assess impact of shading during heatwaves (15–17 June, 17–19 July and 9–15 August)
    # heatwave_analysis_shading("actual")

    # [*] Model validation, statistical similarity    
    # Method = ["Spearman", "Pearson", "Kendall"]
    # Not called directly, but utilised in draw_figure_3()
    #time_series_correlation(True, 'lalc_73', depth=40, method="Spearman") # 2024-03-08

    #==================== Figures ====================#

    # Figure 1: Waffle plot and map of urban canal network
    #draw_figure_1()
    
    # Figure 2: Time-series of modelled wT and cT, plus map of feature
    # Inputs: Feature name to plot, including (True) or excluding (False) shading, month to highlight [1-12]
    # draw_figure_2("rc2_23", True, 7) # rc2_23

    # Figure 3: Plot measured vs. modelled wT
    # Inputs: Oncluding (True) or excluding (False) shading, list of canal locations, water depth (cm)
    # draw_figure_3(True, ['lalc_73', 'nabc_2', 'suc_41', 'batc_14', 'hc_53', 'cc3_11'], depth=40) 

    # *NEW* Figure 4: Fluctuations in median air temperature change at desired interval, plus histogram
    # Inputs: Start date, interval ("week", "fortnight", "month"), including (True) or excluding (False) shading
    # Updated to a 2x2 plot to illustrate daytime and nighttime change
    # draw_figure_4("2022-01-01 10:00:00", "week", True)

    ''' > Figure 5 was created manually, using the outputs of extract_summary_data(), print_summary_values(), shading_analysis() and print_shading_summary() '''
    
    # Create gif showing shadows cast by buildings
    # Inputs: feature_name, start date, duration of analysis (days), delete (True) or keep (False) image files after completion
    # shading_gif('rc_9', "2022-06-01 01:00:00", 0.5, True)

    #==================== Supplementary ====================#

    # [*] Model spin-up, selecting the feature to plot
    # draw_supplementary_figure_A("rc2_23")

    # [*] Area buffer ratio
    # Inputs: include (True) or exclude shading (False)
    #draw_supplementary_figure_B(True)

    # [*] Canal area vs. geometric buffer and area-normalised buffer
    # draw_supplementary_figure_C()

    # [*] Plot to showcase the difference in water and concrete temperatures when shading is included
    # High shading example ('rc_9'), low shading example ('lalc_170')
    #draw_supplementary_figure_D('rc_9')

    # Completion
    print(f"Code completion in {timer() - start} seconds")

#==================== Analysis ====================#

def extract_summary_data(shading_boolean):
    '''
    > Extract summary data for each canal feature, including and excluding shading effects
    > Outputs stored at f"../outputs/feature-values-{including | excluding}-shading.json"
    '''
    
    # File paths for Britain and Ireland
    file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", "../data/urban-canals-ireland-dissolved-id-tm65.shp"]

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # If the summary file already exists, exit the function
    if isfile(f"../outputs/feature-values-{folder_path}-shading.json"):
        print(f"Output file: '../outputs/feature-values-{folder_path}-shading.json' EXISTS. Exiting function.")
        return    

    # Init output dict
    output = {}

    # Iterate and read each file
    for file in file_names:
        canals_gdf = read_dataframe(file)

        # Iterate efficiently using itertuples
        for feature_tuple in canals_gdf.itertuples():
            
            # Try to load  results for this feature
            try:
                # Load results excluding | including shading
                with open(f"../results/{folder_path}-shading/model-output-{feature_tuple.code_id}.json") as canal_data_path:
                    canal_data = load(canal_data_path)

            # If the file is missing (e.g, < 500 m^2 | > 200_000 m^2)
            except FileNotFoundError:

                # Skip to next canal feature
                print(f"No file for {feature_tuple.code_id}.json")
                continue

            # If there is an error with the file
            except decoder.JSONDecodeError:
                
                # Skip to next canal feature
                print(f"Error with {feature_tuple.code_id}.json")
                continue

            # Init dicts for daytime | nighttime values
            peak_day, peak_night = {}, {}

            # Iterate through the data dict
            for key, value in canal_data.items():

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
                day_month_list[month].append(value['temperature_difference'])

                # Ignore final week and fortnight 
                if week < 53:
                    day_week_list[week].append(value['temperature_difference'])
                if fortnight < 27:
                    day_fortnight_list[fortnight].append(value['temperature_difference'])

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

            # Extract air temperatures
            day_temps = [x['temperature_difference'] for x in peak_day.values()]
            night_temps = [x['temperature_difference'] for x in peak_night.values()]

            # Extract energy values for daytime
            day_energy = [x['energy_difference'] for x in peak_day.values()]
            
            # Calculate radius and draw circle
            radius = sqrt(feature_tuple.geometry.area / pi)
            circle = Point(0, 0).buffer(radius, resolution = 10_000)

            # For modelling air temperature change, we use
            buffered_area = circle.buffer(50).area 

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
                                            'energy values' : {'median' : median(day_energy), 
                                                            'q1' : quantile(day_energy, 0.25),
                                                            'q3' : quantile(day_energy, 0.75)},
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
            # print(f"Completed ../results/excluding-shading/model-output-{feature_tuple.code_id}.json")
   
    # If the output dictionary contains data
    if bool(output):

        # When all features have been summarised, save to outputs
        dump(output, open(f"../outputs/feature-values-{folder_path}-shading.json", 'w'))

    # No output
    else:
        print(f"The output dictionary is EMPTY. '../outputs/feature-values-{folder_path}-shading.json' has not been created")
        return

def print_summary_values(shading_boolean):
    '''
    > Utilises summary data for each canal feature produced via extract_summary_data(), either including or excluding shading effects
    > Returns key summary values e.g., annual daytime and nighttime medians, monthly summaries
    > Results printed to console and saved to file "../outputs/summary-values-{folder_path}-shading.json"
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # If the summary file already exists, exit the function
    if isfile(f"../outputs/summary-values-{folder_path}-shading.json"):
        print(f"Output file: '../outputs/summary-values-{folder_path}-shading.json' EXISTS. Loading...")
       
        # Load file
        with open(f"../outputs/summary-values-{folder_path}-shading.json") as file_path:
            output = load(file_path)

            # Title
            print(f"Results {folder_path} shading effects ==================")

            # Iterate and print output
            for key, value in output.items():
                print(f"For {key}, median ΔPT = {value['median']:.2f} (Q1 = {value['q1']:.2f}, Q3 = {value['q3']:.2f})")

        # Exit function
        return    

    # Try and load summary data, including | excluding shading
    try:
        with open(f"../outputs/feature-values-{folder_path}-shading.json") as summary_path:
            summary = load(summary_path)

    # If file does not exist, exit function
    except FileNotFoundError:
        print(f"../outputs/feature-values-{folder_path}-shading.json does not exist")
        return 
    
    # Initialise output dictionary
    output_dict = {}

    # Extract and print median daytime values
    day_medians = array([x['day values']['median'] for x in summary.values()])

    print(f"Day statistics (annual) {folder_path} shading -----------------")
    print(f"Median day ΔPT of {median(day_medians):.2f}. Q1 = {quantile(day_medians, 0.25)}, Q3 = {quantile(day_medians, 0.75):.2f}")
    print()

    # Save to output
    output_dict.update({'Annual day' : {'median' : median(day_medians),
                                    'q1' : quantile(day_medians, 0.25),
                                    'q3' : quantile(day_medians, 0.75)}})
    
    # Extract and print median nightime values
    night_medians = [x['night values']['median'] for x in summary.values()]

    print(f"Night statistics (annual) {folder_path} shading -----------------")
    print(f"Median night ΔPT of {median(night_medians):.2f}. Q1 = {quantile(night_medians, 0.25):.2f}, Q3 = {quantile(night_medians, 0.75):.2f}")
    print()

    # Save to output
    output_dict.update({'Annual night' : {'median' : median(night_medians),
                                'q1' : quantile(night_medians, 0.25),
                                'q3' : quantile(night_medians, 0.75)}})

    # Extract median values for each month (daytime + nightime)
    mon_day_medians = [x['monthly values']['day medians'] for x in summary.values()]
    mon_night_medians = [x['monthly values']['night medians'] for x in summary.values()]

    # Daytime Median, Q1, Q3
    day_medians = np_median(mon_day_medians,axis=0)
    day_q3 = quantile(mon_day_medians, 0.75, axis=0)
    day_q1 = quantile(mon_day_medians, 0.25, axis=0) 

    # Nighttime Median, Q1, Q3
    night_medians = np_median(mon_night_medians,axis=0)
    night_q3 = quantile(mon_night_medians, 0.75, axis=0)
    night_q1 = quantile(mon_night_medians, 0.25, axis=0) 

    # List of months
    months = ["January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

    print(f"Statistics (monthly) {folder_path} shading -----------------")
    # Iterate and print summary statitics
    for index, m in enumerate(months):
        print(f"For {m}, the daytime ΔPT is {day_medians[index]:.2f} ({day_q1[index]:.2f} - {day_q3[index]:.2f}), with nighttime ΔPT of {night_medians[index]:.2f} ({night_q1[index]:.2f} - {night_q3[index]:.2f})")

        # Save to output (daytime)
        output_dict.update({f"{m} day" : {'median' : day_medians[index], 
                                    'q1' : day_q1[index],
                                    'q3' : day_q3 [index]}})
        
        # Save to output (nighttime)
        output_dict.update({f"{m} night" : {'median' : night_medians[index], 
                                    'q1' : night_q1[index],
                                    'q3' : night_q3[index]}})
                                 
    # If the output dictionary contains data
    if bool(output_dict):

        # When all features have been summarised, save to outputs
        dump(output_dict, open(f"../outputs/summary-values-{folder_path}-shading.json", 'w'))

    # No output
    else:
        print(f"The output dictionary is EMPTY. '../outputs/summary-values-{folder_path}-shading.json' has not been created")
        return
    
def heatwave_analysis(period):
    '''
    > Extracting daytime air temperature change for the three classified heatwaves (15–17 June, 17–19 July and 9–15 August)
    > Printed outputs are based on the full feature set (n = 2,356) and a variable reference material (concrete shading included / excluded) 
    > As a result, the difference due to shading is negligible  
    '''

    # Heatwave days only
    if period == 'actual':
        heatwave_periods =  [['166', '167', '168'], # 15–17 June
                             ['198', '199', '200'], # 17–19 July
                             ['221', '222', '223', '224', '225', '226', '227']] # 9–15 August
    
    # Include an extended period 
    elif period == 'extended':
        heatwave_periods =  [['164', '165', '166', '167', '168', '169', '170'], 
                             ['196', '197', '198', '199', '200', '201', '202'], 
                             ['219', '220' '221', '222', '223', '224', '225', '226', '227', '228', '229']]
        
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

    # For iterating, and printing results
    heatwave_months = ['June', 'July', 'August']

    # If the summary files already exist, load the files
    if isfile(f"../outputs/heatwave-record-day-{period}.json") and isfile(f"../outputs/heatwave-record-night-{period}.json"):
        print(f"Output files: '../outputs/heatwave-record-[day | night]-{period}.json' EXIST. Loading...")

        # Load day file
        with open(f"../outputs/heatwave-record-day-{period}.json") as day_path:
            output = load(day_path)
            
         # Load nighttime file
        with open(f"../outputs/heatwave-record-night-{period}.json") as night_path:
            night_output = load(night_path) 

    # Do the calculations
    else:

        # File paths for Britain and Ireland
        file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", "../data/urban-canals-ireland-dissolved-id-tm65.shp"]

        # Init output dicts
        output = {}
        night_output = {}

        # Iterate and read each file
        for file in file_names:
            canals_gdf = read_dataframe(file)

            # Iterate efficiently using itertuples
            for feature_tuple in canals_gdf.itertuples():

                # Skip canals with an area < 500 m^2 or > 200_000 m^2
                if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
                    continue
            
                # Load results including shading
                try:
                    with open(f"../results/including-shading/model-output-{feature_tuple.code_id}.json") as including_path:
                        including_data = load(including_path)

                    # Excluding shading 
                    with open(f"../results/excluding-shading/model-output-{feature_tuple.code_id}.json") as excluding_path:
                        excluding_data = load(excluding_path)

                # If a file is missing 
                except FileNotFoundError:

                    # Exit function
                    print(f"Missing summary files for feature {feature_tuple.code_id}")
                    return

                # If there is an error with the file
                except decoder.JSONDecodeError:

                    # Exit function
                    print(f"Error with feature {feature_tuple.code_id}")
                    return 
                
                # Add key to output, empty dict
                output[feature_tuple.code_id] = {}
                night_output[feature_tuple.code_id] = {}

                # Init output dicts
                including_filter = {}
                excluding_filter = {}
                including_night_filter = {}
                excluding_night_filter = {}

                # Init dict to store heatwave days and datetimes
                datetimes = {} # dict.fromkeys(heatwave_days, [])#[]) 
                night_datetimes = {}

                # Iterate through dictionary and subset to months of interest
                for key, value in including_data.items():

                    # Convert unix to datetime, extract day of year
                    local_datetime = datetime.fromtimestamp(int(key))
                    current_day = str(local_datetime.timetuple().tm_yday)

                    # Filter to heatwaves only
                    if any(current_day in d for d in heatwave_periods):

                        # Filter to daytime, classified as 10:00-15:45 (5h 45 m)
                        if local_datetime.hour >= DAY_START and local_datetime.hour < DAY_END:

                            # Add to output
                            including_filter[key] = value
                            excluding_filter[key] = excluding_data[key]

                            # Update datetime dict
                            datetimes[key] = current_day
                    
                    # Filter to heatwave days + 1 (including the night at the end)
                    if any(current_day in d for d in nightwave_periods):

                        # Filter to nighttime, classified as 22:00-03:45 (5h 45 m)
                        if local_datetime.hour >= NIGHT_START or local_datetime.hour < NIGHT_END:

                            # Remove the extra night (i.e., after 22:00 the following day)
                            if any(current_day in d for d in heatwave_ends) and local_datetime.hour >= NIGHT_START: 
                                continue
                            
                            # Add to output
                            else:
                                including_night_filter[key] = value
                                excluding_night_filter[key] = excluding_data[key]

                                # Update datetime dict
                                night_datetimes[key] = current_day
                    
                # Iterate through the heatwaves
                for heatwave, nightwave, month in zip(heatwave_periods, nightwave_periods, heatwave_months):

                    # Extract Unix times (keys) that correspond to each
                    keys = [k for k, v in datetimes.items() if v in heatwave]
                    night_keys = [k for k, v in night_datetimes.items() if v in nightwave]

                    # Calculate median air temperature change for daytime
                    output[feature_tuple.code_id].update({f'{month} median day, including shading' : median([including_filter[x]['temperature_difference'] for x in keys]),
                                                            f'{month} median day, excluding shading' : median([excluding_filter[x]['temperature_difference'] for x in keys])})
                    
                    # Calculate median air temperature change for nighttime
                    night_output[feature_tuple.code_id].update({
                                                            f'{month} median night, including shading' : median([including_night_filter[x]['temperature_difference'] for x in night_keys]),
                                                            f'{month} median night, excluding shading' : median([excluding_night_filter[x]['temperature_difference'] for x in night_keys])})
                         
        # When all features have been summarised, save to outputs
        dump(output, open(f"../outputs/heatwave-record-day-{period}.json", 'w'))
        dump(night_output, open(f"../outputs/heatwave-record-night-{period}.json", 'w'))

    # Iterate for printing
    for shading_boolean in ['excluding', 'including']:
        for month in heatwave_months:

            # using item() to extract key value pair as whole
            res = [val[f'{month} median day, {shading_boolean} shading'] for key, val in output.items() if f'{month} median day, {shading_boolean} shading' in val]

            # Summary statement
            print(f"For the {month} heatwave, median daytime ΔT ({shading_boolean} shading) is {median(res):.2f}. Q1 = {quantile(res, 0.25)}, Q3 = {quantile(res, 0.75):.2f}")

            # using item() to extract key value pair as whole
            res = [val[f'{month} median night, {shading_boolean} shading'] for key, val in night_output.items() if f'{month} median night, {shading_boolean} shading' in val]

            # Summary statement
            print(f"For the {month} heatwave, median nighttime ΔT ({shading_boolean} shading) is {median(res):.2f}. Q1 = {quantile(res, 0.25)}, Q3 = {quantile(res, 0.75):.2f}")
    
def shading_analysis():
    '''
    > Simplified shading effects analysis, using unshaded concrete as the baseline for comparison
    > Returns annual median values
    > Returns monthly values
    '''

    # Check for output file
    if isfile(f"../outputs/feature-values-shading-effects.json"):
        print("Output file '../outputs/feature-values-shading-effects.json' EXISTS. Exiting function...")

    # Perform analysis
    else:

        # Try and load list of "highly" shaded locations (> threshold)
        try:
            with open(f"../outputs/high-shading-features.json") as shading_subset:
                high_shading = load(shading_subset)

                # Extract feature IDs from dict keys
                high_shading_features = high_shading.keys()

        # If the file is not found, print warning and exit function
        except FileNotFoundError:
            print("WARNING. No record of high-shading features. Check file path or run extract_highly_shaded()")
            return 0

        # Init output dictionary
        output = {}

        # Iterate through high shading features
        for feature in high_shading_features:
                
            # Try to load results for this feature
            try:
                # Excluding and including shading
                with open(f"../results/excluding-shading/model-output-{feature}.json") as excluding_path:
                    excluding_data = load(excluding_path)
                with open(f"../results/including-shading/model-output-{feature}.json") as including_path:
                    including_data = load(including_path)

            # If the file is missing 
            except FileNotFoundError:

                # Skip to next canal feature
                print(f"No file for {feature}.json")
                continue

            # If there is an error with the file
            except decoder.JSONDecodeError:
                
                # Skip to next canal feature
                print(f"Error with {feature}.json")
                continue

            # Init dicts for daytime | nighttime values
            peak_day, peak_night = {}, {}

            # Iterate through the excluding data dict
            for key, value in excluding_data.items():

                # Convert unix to datetime
                dat = datetime.fromtimestamp(int(key))

                # Between 10:00 and 15:45 (5h 45 m)
                if dat.hour >= DAY_START and dat.hour < DAY_END:

                    # Add to daytime dictionary, storing the water and concrete temperature from BOTH dicts
                    # Also calculating the temperature difference values:
                    # > "excluding-vs-baseline" = difference between water (excluding) and concrete (excluding)
                    # > "including-vs-baseline" = difference between water (including) and concrete (excluding)
                    peak_day[key] = {"water ΔT excluding" : value['water ΔT'],
                                    "water ΔT including" : including_data[key]['water ΔT'],
                                    "concerete ΔT excluding" : value['concrete ΔT'],
                                    "concerete ΔT including" : including_data[key]['concrete ΔT'],
                                    "excluding-vs-baseline" : value['water ΔT'] - value['concrete ΔT'],
                                    "including-vs-baseline" : including_data[key]['water ΔT'] - value['concrete ΔT']}

                # Between 22:00 and 03:45 (5h 45m)
                elif dat.hour >= NIGHT_START or dat.hour < NIGHT_END:

                    try:

                        # Add to nightime dictionary, as above
                        peak_night[key] = {"water ΔT excluding" : value['water ΔT'],
                                        "water ΔT including" : including_data[key]['water ΔT'],
                                        "concerete ΔT excluding" : value['concrete ΔT'],
                                        "concerete ΔT including" : including_data[key]['concrete ΔT'],
                                        "excluding-vs-baseline" : value['water ΔT'] - value['concrete ΔT'],
                                        "including-vs-baseline" : including_data[key]['water ΔT'] - value['concrete ΔT']}
                    
                    except KeyError:
                        print("KeyError: Value missing?")

            # Init lists to store monthly data
            day_month_excluding = defaultdict(list)
            day_month_including = defaultdict(list)
            
            # As above, but for night values
            night_month_excluding = defaultdict(list)
            night_month_including = defaultdict(list)

            # Iterate through the daytime values
            for key, value in peak_day.items():

                # Convert unix to datetime, extract month integer
                month = datetime.fromtimestamp(int(key)).month

                # Append values to corresponding list with key
                day_month_excluding[month].append(value['excluding-vs-baseline'])
                day_month_including[month].append(value['including-vs-baseline'])

            # Iterate through the nighttime values
            for key, value in peak_night.items():

                # Convert unix to datetime, extract month integer
                month = datetime.fromtimestamp(int(key)).month

                # Append values to corresponding key
                night_month_excluding[month].append(value['excluding-vs-baseline'])
                night_month_including[month].append(value['including-vs-baseline'])         
            
            # Extract air temperature difference (ΔT) values
            day_excluding_baseline = [x['excluding-vs-baseline'] for x in peak_day.values()]
            night_excluding_baseline = [x['excluding-vs-baseline'] for x in peak_night.values()]

            day_including_baseline = [x['including-vs-baseline'] for x in peak_day.values()]
            night_including_baseline = [x['including-vs-baseline'] for x in peak_night.values()]

        
            # Add summary values to output
            output[feature] = {
                                # Median (Q1:Q3) daytime values
                                'day excluding values' : {'median' : median(day_excluding_baseline), 
                                                'q1' : quantile(day_excluding_baseline, 0.25),
                                                'q3' : quantile(day_excluding_baseline, 0.75)},
                                # Median (Q1:Q3) daytime values
                                'day including values' : {'median' : median(day_including_baseline), 
                                                'q1' : quantile(day_including_baseline, 0.25),
                                                'q3' : quantile(day_including_baseline, 0.75)},
                                # Median (Q1:Q3) nighttime values
                                'night excluding values' : {'median' : median(night_excluding_baseline), 
                                                'q1' : quantile(night_excluding_baseline, 0.25),
                                                'q3' : quantile(night_excluding_baseline, 0.75)},
                                # Median (Q1:Q3) nighttime values
                                'night including values' : {'median' : median(night_including_baseline), 
                                                'q1' : quantile(night_including_baseline, 0.25),
                                                'q3' : quantile(night_including_baseline, 0.75)},
                                # Monthly values, daytime and nighttime (median, Q1:Q3), excluding shading (in both concrete and water)            
                                'monthly excluding values' : {'day medians' : [median(day_month_excluding[x]) for x in day_month_excluding], 
                                                    'day q1' : [quantile(day_month_excluding[x], 0.25) for x in day_month_excluding],
                                                    'day q3' : [quantile(day_month_excluding[x], 0.75) for x in day_month_excluding],
                                                    'night medians' : [median(night_month_excluding[x]) for x in night_month_excluding],
                                                    'night q1' : [quantile(night_month_excluding[x], 0.25) for x in night_month_excluding],
                                                    'night q3' : [quantile(night_month_excluding[x], 0.75) for x in night_month_excluding]},
                                # Monthly values, daytime and nighttime (median, Q1:Q3), including shading (in water, excluding shading for concrete)             
                                'monthly including values' : {'day medians' : [median(day_month_including[x]) for x in day_month_including], 
                                                    'day q1' : [quantile(day_month_including[x], 0.25) for x in day_month_including],
                                                    'day q3' : [quantile(day_month_including[x], 0.75) for x in day_month_including],
                                                    'night medians' : [median(night_month_including[x]) for x in night_month_including],
                                                    'night q1' : [quantile(night_month_including[x], 0.25) for x in night_month_including],
                                                    'night q3' : [quantile(night_month_including[x], 0.75) for x in night_month_including]}}


        # If the output dictionary contains data
        if bool(output):

            # When all features have been summarised, save to outputs
            dump(output, open(f"../outputs/feature-values-shading-effects.json", 'w'))

        # No output
        else:
            print(f"The output dictionary is EMPTY. '../outputs/feature-values-shading-effects.json' has not been created")
            return

def print_shading_summary():
    '''
    > Utilises summary data for each canal feature produced via shading_analysis(),
    > Returns key summary values e.g., annual daytime and nighttime medians, monthly summaries
    > Results printed to console and saved to file "../outputs/summary-values-{folder_path}-shading.json"
    '''

    # Try and load summary data, including | excluding shading
    try:
        with open(f"../outputs/feature-values-shading-effects.json") as summary_path:
            summary = load(summary_path)

    # If file does not exist
    except FileNotFoundError:
        print(f"File not found: '../outputs/feature-values-shading-effects.json'. Try running shading_analysis(). Exiting function...")
        return

    # Extract and print median daytime values
    day_exc_medians = array([x['day excluding values']['median'] for x in summary.values()])
    day_inc_medians = array([x['day including values']['median'] for x in summary.values()])

    print(f"Day statistics (annual) -----------------")
    print(f"Median day cooling (excluding shading of water) of {median(day_exc_medians):.2f}. Q1 = {quantile(day_exc_medians, 0.25)}, Q3 = {quantile(day_exc_medians, 0.75):.2f}")
    print(f"Median day cooling (including shading of water) of {median(day_inc_medians):.2f}. Q1 = {quantile(day_inc_medians, 0.25)}, Q3 = {quantile(day_inc_medians, 0.75):.2f}")
    print()


    # Extract and print median nightime values
    night_exc_medians = array([x['night excluding values']['median'] for x in summary.values()])
    night_inc_medians = array([x['night including values']['median'] for x in summary.values()])

    print(f"Night statistics (annual)  -----------------")
    print(f"Median night warming (excluding shading of water) of {median(night_exc_medians):.2f}. Q1 = {quantile(night_exc_medians, 0.25):.2f}, Q3 = {quantile(night_exc_medians, 0.75):.2f}")
    print(f"Median night warming (including shading of water) of {median(night_inc_medians):.2f}. Q1 = {quantile(night_inc_medians, 0.25):.2f}, Q3 = {quantile(night_inc_medians, 0.75):.2f}")
    print()


    # Extract median values for each month (daytime)
    mon_day_exc_medians = [x['monthly excluding values']['day medians'] for x in summary.values()]
    mon_day_inc_medians = [x['monthly including values']['day medians'] for x in summary.values()]

    # Extract median values for each month (daytime)
    mon_night_exc_medians = [x['monthly excluding values']['night medians'] for x in summary.values()]
    mon_night_inc_medians = [x['monthly including values']['night medians'] for x in summary.values()]

    # Daytime Median, Q1, Q3, excluding shading effects
    day_exc_medians = np_median(mon_day_exc_medians,axis=0)
    day_exc_q3 = quantile(mon_day_exc_medians, 0.75, axis=0)
    day_exc_q1 = quantile(mon_day_exc_medians, 0.25, axis=0) 

    # Nighttime Median, Q1, Q3, excluding shading effects
    night_exc_medians = np_median(mon_night_exc_medians,axis=0)
    night_exc_q3 = quantile(mon_night_exc_medians, 0.75, axis=0)
    night_exc_q1 = quantile(mon_night_exc_medians, 0.25, axis=0) 

    # Daytime Median, Q1, Q3, including shading effects
    day_inc_medians = np_median(mon_day_inc_medians,axis=0)
    day_inc_q3 = quantile(mon_day_inc_medians, 0.75, axis=0)
    day_inc_q1 = quantile(mon_day_inc_medians, 0.25, axis=0) 

    # Nighttime Median, Q1, Q3,  including shading effects
    night_inc_medians = np_median(mon_night_inc_medians,axis=0)
    night_inc_q3 = quantile(mon_night_inc_medians, 0.75, axis=0)
    night_inc_q1 = quantile(mon_night_inc_medians, 0.25, axis=0) 

    # List of months
    months = ["January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

    print(f"Statistics (monthly) shading -----------------")
    # Iterate and print summary statitics
    for index, m in enumerate(months):
        print(f"For {m}:")
        print(f"- the daytime cooling (excluding shading) is {day_exc_medians[index]:.2f} ({day_exc_q1[index]:.2f} - {day_exc_q3[index]:.2f}), with nighttime warming of {night_exc_medians[index]:.2f} ({night_exc_q1[index]:.2f} - {night_exc_q3[index]:.2f})")
        print(f"- the daytime cooling (including shading) is {day_inc_medians[index]:.2f} ({day_inc_q1[index]:.2f} - {day_inc_q3[index]:.2f}), with nighttime warming of {night_inc_medians[index]:.2f} ({night_inc_q1[index]:.2f} - {night_inc_q3[index]:.2f})")
        print()

def heatwave_analysis_shading(period):
    '''
    > Extracting daytime air temperature change for the three classified heatwaves (15–17 June, 17–19 July and 9–15 August)
    > Unlike heatwave_analysis(), printed outputs are based on the refined set of highly shaded features (n = 305) and a fixed reference material (concrete shading always excluded)   
    > Requires output of extract_highly_shaded() in auxiliary.py
    '''

    # Try and load list of "highly" shaded locations (> threshold)
    try:
        with open(f"../outputs/high-shading-features.json") as shading_subset:
            high_shading = load(shading_subset)

            # Extract feature IDs from dict keys
            high_shading_features = high_shading.keys()

    # If the file is not found, print warning and exit function
    except FileNotFoundError:
        print("WARNING. No record of high-shading features. Check file path or run extract_highly_shaded()")
        return

    # If the summary file already exists, load the file
    if isfile(f"../outputs/heatwave-record-shading-day-{period}.json") and isfile(f"../outputs/heatwave-record-shading-night-{period}.json"):
        print(f"Output files: '../outputs/heatwave-record-shading-[day | night]-{period}.json' EXIST. Loading...")

        # Load day file
        with open(f"../outputs/heatwave-record-shading-day-{period}.json") as file_path:
            output = load(file_path) 
        # Load night file
        with open(f"../outputs/heatwave-record-shading-night-{period}.json") as night_path:
            night_output = load(night_path)
    
    # Do the calculations
    else:

        # Heatwave days only
        if period == 'actual':
            heatwave_periods =  [['166', '167', '168'], # 15–17 June
                                ['198', '199', '200'], # 17–19 July
                                ['221', '222', '223', '224', '225', '226', '227']] # 9–15 August
        
        # Include an extended period 
        elif period == 'extended':
            heatwave_periods =  [['164', '165', '166', '167', '168', '169', '170'], 
                                ['196', '197', '198', '199', '200', '201', '202'], 
                                ['219', '220' '221', '222', '223', '224', '225', '226', '227', '228', '229']]
            
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

        # For iterating, and printing results
        heatwave_months = ['June', 'July', 'August']

        # Init output dicts
        output = {}
        night_output = {}

        # Iterate through high shading features
        for feature in high_shading_features:
                
            # Try to load results for this feature
            try:
                # Excluding and including shading
                with open(f"../results/excluding-shading/model-output-{feature}.json") as excluding_path:
                    excluding_data = load(excluding_path)
                with open(f"../results/including-shading/model-output-{feature}.json") as including_path:
                    including_data = load(including_path)

            # If the file is missing 
            except FileNotFoundError:

                # Skip to next canal feature
                print(f"No file for {feature}.json")
                continue

            # If there is an error with the file
            except decoder.JSONDecodeError:
                
                # Skip to next canal feature
                print(f"Error with {feature}.json")
                continue

            # Add key to output, empty dict
            output[feature] = {}
            night_output[feature] = {}

            # Init output dicts
            including_filter = {}
            excluding_filter = {}
            including_night_filter = {}
            excluding_night_filter = {}

            # Init dict to store heatwave days and datetimes
            datetimes = {} # dict.fromkeys(heatwave_days, [])#[]) 
            night_datetimes = {}

            # Iterate through dictionary and subset to months of interest
            for key, value in excluding_data.items():

                # Convert unix to datetime, extract day of year
                local_datetime = datetime.fromtimestamp(int(key))
                current_day = str(local_datetime.timetuple().tm_yday)

                # Filter to heatwaves only
                if any(current_day in d for d in heatwave_periods):

                    # Filter to daytime, classified as 10:00-15:45 (5h 45 m)
                    if local_datetime.hour >= DAY_START and local_datetime.hour < DAY_END:

                        # Add to output
                        excluding_filter[key] = value['water ΔT'] - value['concrete ΔT']
                        including_filter[key] = including_data[key]['water ΔT'] - value['concrete ΔT']

                        # Update datetime dict
                        datetimes[key] = current_day
                
                # Filter to heatwave days + 1 (including the night at the end)
                if any(current_day in d for d in nightwave_periods):

                    # Filter to nighttime, classified as 22:00-03:45 (5h 45 m)
                    if local_datetime.hour >= NIGHT_START or local_datetime.hour < NIGHT_END:

                        # Remove the extra night (i.e., after 22:00 the following day)
                        if any(current_day in d for d in heatwave_ends) and local_datetime.hour >= NIGHT_START: 
                            continue
                        
                        # Add to output
                        else:
                            excluding_night_filter[key] = value['water ΔT'] - value['concrete ΔT']
                            including_night_filter[key] = including_data[key]['water ΔT'] - value['concrete ΔT']

                            # Update datetime dict
                            night_datetimes[key] = current_day
                     
            # Iterate through the heatwaves
            for heatwave, nightwave, month in zip(heatwave_periods, nightwave_periods, heatwave_months):

                # Extract Unix times (keys) that correspond to each
                keys = [k for k, v in datetimes.items() if v in heatwave]
                night_keys = [k for k, v in night_datetimes.items() if v in nightwave]

                # Calculate median air temperature change for daytime
                output[feature].update({f'{month} median day, including shading' : median([including_filter[x] for x in keys]),
                                        f'{month} median day, excluding shading' : median([excluding_filter[x] for x in keys])})
                
                # Calculate median air temperature change for nighttime
                night_output[feature].update({
                                            f'{month} median night, including shading' : median([including_night_filter[x] for x in night_keys]),
                                            f'{month} median night, excluding shading' : median([excluding_night_filter[x] for x in night_keys])})
                         
        # When all features have been summarised, save to outputs
        dump(output, open(f"../outputs/heatwave-record-shading-day-{period}.json", 'w')) 
        dump(night_output, open(f"../outputs/heatwave-record-shading-night-{period}.json", 'w'))

    # Iterate for printing
    for shading_boolean in ['excluding', 'including']:
        for month in heatwave_months:

            # using item() to extract key value pair as whole
            res = [val[f'{month} median day, {shading_boolean} shading'] for key, val in output.items() if f'{month} median day, {shading_boolean} shading' in val]

            # Summary statement
            print(f"For the {month} heatwave, median daytime ΔT ({shading_boolean} shading) is {median(res):.2f}. Q1 = {quantile(res, 0.25)}, Q3 = {quantile(res, 0.75):.2f}")

            # using item() to extract key value pair as whole
            res = [val[f'{month} median night, {shading_boolean} shading'] for key, val in night_output.items() if f'{month} median night, {shading_boolean} shading' in val]

            # Summary statement
            print(f"For the {month} heatwave, median nighttime ΔT ({shading_boolean} shading) is {median(res):.2f}. Q1 = {quantile(res, 0.25)}, Q3 = {quantile(res, 0.75):.2f}")

def time_series_correlation(shading_boolean, canal_id, depth, method):
    '''
    > Function for measuring similarity of time-series (modelled wT vs. measured wT) using a range of methods including:
        - Pearson
        - Spearman (preferred)
        - Kendall 
    > Sources:
        - https://stats.stackexchange.com/questions/133155/how-to-use-pearson-correlation-correctly-with-time-series
        - https://otexts.com/fpp2/stationarity.html
    > Other options include dynamic time warping (DTW)
    > Measured data from the CRT (https://canalrivertrust.org.uk/)
    > Of the monitoring locations, the following are suitable for comparison:
        - Bridge Pagefield Pipe Crossing (lalc_73): 10 m frequency
        - Bridge 15 Chain Lane Bridge (nabc_2): 10 m frequency
        - New Road Bridge 148 (suc_41): 1 h frequency
        - Bridge 32. Priorswood (batc_14): 15 m frequency (from 03-01-2022)
        - Bridge 4 Wakefield Road Bridge (hc_53): 1 h frequency
        - Bridge 73 Anchor Bridge (cc3_11): 10 m frequency
    > The following locations are less suitable for comparison:
        - Anderton Waste Weir Flow (None): not directly modelled
        - Bridge 46 Benthouse Bridge (None): not directly modelled
        - Hawkesbury Lock (None): not directly modelled
        - Bridge C Lightbody Street (lalc_4): missing measurement data
        - Bridge 40 Taverners Bridge (cc3_53): missing measurement data
    > depth (cm) refers to the upper value
        - 20 = 20 - 40 cm
        - 40 = 40 - 60 cm
        - 60 = 60 - 80 cm
        - 80 = 80 - 100 cm
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

    # Try
    try:

        # Include shading
        if shading_boolean:
        
            # Load model output record, including shading
            with open(f"../results/including-shading/model-output-{canal_id}.json") as results_path:
                canal_output = load(results_path)

        # Exclude shading
        else:

            # Load model output record, excluding shading
            with open(f"../results/excluding-shading/model-output-{canal_id}.json") as results_path:
                canal_output = load(results_path)

    except FileNotFoundError:
        print("No modelled data to compare")
        exit()

    # Init dicts to store output values
    measured_dict_filter = {}
    modelled_dict_filter = {}

    # Iterate through the measured values dictionay (deepcopy)
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
    modeleld_diff = ediff1d(modelled_wt)

    # Method selection
    if method == "Spearman":
        r, p = stats.spearmanr(measured_diff, modeleld_diff)

    elif method == "Pearson":
        r, p = stats.pearsonr(measured_diff, modeleld_diff)

    elif method == "Kendall": 
        r, p = stats.kendalltau(measured_diff, modeleld_diff)

    else:
        print("Statistical method invalid. Check inputs")

    # Print summary
    print(f"For {canal_id}, {method} correlation = {r:.2f}, p value = {p:.2f}")

    # Return method, test statistic,  p value, array length
    return method, r, p, len(measured_wt)

#==================== Figures ====================#

def draw_figure_1():

    # Read in the canal dataframes
    britain_gdf = read_dataframe("../data/canal-geometries-filtered-modified-dissolved-id.shp")
    ireland_gdf = read_dataframe("../data/urban-canals-ireland-dissolved-id-tm65.shp")

    # Read in the border polygon
    coastline = read_dataframe("../../national-canal-cooling-data/plotting/borders.shp")

    # Reproject to BNG
    ireland_gdf = ireland_gdf.to_crs(britain_gdf.crs)
    coastline = coastline.to_crs(britain_gdf.crs)

    # Merge
    merged_df = concat([britain_gdf, ireland_gdf])

    # Filter to canals within 500 and 200_000 m^2
    merged_df = merged_df[merged_df['area'].between(500, 200_000)]

    # Produce cluster-statistics (total area)
    stats = merged_df.groupby('region').agg({'area':'sum'})
    
    # Dissolve gdf by region
    dissolved_gdf = merged_df.dissolve(by='region').centroid

    # Produce merged df for plotting
    plot_gdf = GeoDataFrame(stats, geometry = dissolved_gdf, crs = britain_gdf.crs)

    # Normalise [0-1]
    plot_gdf['normalised_area'] = (plot_gdf['area']-plot_gdf['area'].min())/(plot_gdf['area'].max()-plot_gdf['area'].min())

    # km2
    plot_gdf['area_km'] = plot_gdf['area'] * 0.000001

    # Extract the n largest clusters
    cities = plot_gdf.nlargest(6, 'area')

    # Sort 
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

    # Iterate and update colours for other regions
    for name in region_labels:
        plot_gdf.loc[(plot_gdf.index == name), 'edge_color'] = "#242424"

    # List of largest areas, as proportion of total
    urban = cities['area'].tolist()
    urban_prop = [x / sum(plot_gdf['area']) for x in urban]
    urban_prop.append(1 - sum(urban_prop))

    # Waffle chart
    Waffle.make_waffle(
        ax=ax2,  # pass axis to make_waffle
        rows=5, 
        columns=10, 
        values=urban_prop,
        colors=["#754501", "#B76C00", "#FF9700", "#FEAD37", "#FFC26A", "#FFDEAD", "#CBCBCB"],
        legend={'labels': ['Birmingham', 'London', 'Liverpool', 'Manchester', 'Glasgow', 'Dublin', 'other'],
        'loc' : 'upper center',
        'bbox_to_anchor' : (0.5, 1.28),
        'ncol': 4, # len(urban_prop),
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
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True)

    # Axis labels
    ax1.set_xlabel("Easting", fontsize = 9)
    ax1.set_ylabel("Northing", fontsize = 9)

    # Labelling
    ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax1.yaxis.offsetText.set_fontsize(7)
    ax1.xaxis.offsetText.set_fontsize(0)

    # Add text annotations
    ax1.text(-0.05, 1.05, f"B", transform = ax1.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    ax2.text(-0.05, 1.1, f"A", transform = ax2.transAxes, weight='bold', fontsize = 13, va='top', ha='right')

    # Iterate through the n largest areas
    for loc, geom, label, x_offset, y_offset, label_offset in zip(cities['area'], cities['geometry'], 
                                ['Birmingham', 'London', 'Liverpool', 'Manchester', 'Glasgow', 'Dublin'],
                                [0, 0, -117_000, 95_000, 30_000, -30_000],
                                [-30_000, -30_000, 14_000, -20_000, -25_000, -25_000],
                                [4, 4, 4, 4, 5, 10]):
            
            # Add labels to map (xy offsets)
            ax1.text(geom.x + x_offset, geom.y + y_offset, label, fontsize = 7, va='center', ha='center')

    # Scale bar
    ax1.add_artist(ScaleBar(dx=1, units="m", location="upper right", length_fraction=0.2, font_properties={"size": 8}))


    # Save to file
    #show()
    savefig(f'../images/figure-1.png', bbox_inches='tight', dpi = 300)

def draw_figure_2(canal_id, shading_boolean, month):
    '''
    > Draw Figure 2, which includes subplots showing:
        - (A) the modelled concrete and water temperatures for 2022 for the chosen feature
        - (B) a map showcasing the seleceted canal feature, as well as nearby canals and buildings
        - (C) the modelled concrete and water temperatures for a specific month
    > Figure design has been adapted for a canal feature on the Regents Canal, London (rc2_23) although it will work for other inputs
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Read in the canals dataframe and extract selected feature
    canals_gdf = read_dataframe("../data/canal-geometries-filtered-modified-dissolved-id.shp")
    filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]

    # Return buildings within a specified distance
    for feature_tuple in filtered_canals.itertuples():
        buildings = return_buildings_polygon(feature_tuple, 500, plot = False, epsg = canals_gdf.crs)

    # Try and load model output record, including | excluding shading
    try:
        with open(f"../results/{folder_path}-shading/model-output-{canal_id}.json") as results_path:
            canal_output = load(results_path)

    # File not found
    except FileNotFoundError:
        print(f"../results/{folder_path}-shading/model-output-{canal_id}.json does NOT exist")

    # Extract lists of values to plot, converting to celcius
    water_temps = [x['surface_water_k'] - 273.15 for x in canal_output.values()]
    concrete_temps = [x['surface_concrete_k'] - 273.15 for x in canal_output.values()]

    # List of dict keys (unix), convert to datetime format
    dt = list(canal_output.keys())
    dt_datetime = [datetime.utcfromtimestamp(int(x)) for x in dt]

    # Return the indices, corresponding to the chosen month
    idx = [n for n, x in enumerate(dt_datetime) if x.month == month]

    # Filter to indexes
    filter_datetime = [dt_datetime[i] for i in idx]
    filter_water = [water_temps[i] for i in idx]
    filter_concrete = [concrete_temps[i] for i in idx] 

    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig = figure(layout='compressed', figsize=(7, 5))
    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    # Plot the annual water and concrete temperature
    ax1.plot(dt_datetime, concrete_temps, label = 'Concrete', color = "#FFB16F")
    ax1.plot(dt_datetime, water_temps, label ='Water', color = "#707070")
    ax1.set_xlabel("Date")
    ax1.set_ylabel("Temperature (°C)", labelpad=-4)
    ax1.xaxis.tick_top()
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in")
    ax1.xaxis.set_label_position('top') 
    ax1.legend(frameon=False)

    # Patch dimensions
    patch_min = min(min(filter_concrete), min(filter_water))
    patch_max = max(max(filter_concrete), max(filter_water))
    temp_buffer = 2

    # Add rectangle, showing the dimensions of ax3
    ax1.add_patch(Rectangle((filter_datetime[0],patch_min - temp_buffer), filter_datetime[-1] - filter_datetime[0], 
                            (patch_max + temp_buffer) - (patch_min - temp_buffer),
                            fill = False, color=None, # alpha=0.5,
                            edgecolor="#4D4D4D", lw = 1, zorder=2, linestyle = '--'))

    # For ax2, plot the building geometries
    buildings.plot(
        ax = ax2,
        color = "#6D6D6D",
        edgecolor = '#343434',
        linewidth = 0,
        )
            
    # Plot the canal geometries
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
    ax2_legend_top = ax2.legend(handles=[canal_polygon, building_polygon], loc = 'upper left',ncol=2)
    ax2_legend_top.get_frame().set_linewidth(0.0)
    ax2.add_artist(ax2_legend_top)

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
    ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax2.yaxis.offsetText.set_fontsize(7)
    ax2.xaxis.offsetText.set_fontsize(0)
    ax2.set_xlabel("Easting")
    ax2.set_ylabel("Northing")
    ax2.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='y', labelsize=8)

    # Plot the monthly water and concrete temperature record
    ax3.plot(filter_datetime, filter_concrete, color = "#FFAD28", label = "Concrete (°C)")
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
    ax3.text(-0.05, 1.05, f"C", transform = ax3.transAxes, weight='bold', fontsize = 13, va='top', ha='right')
    
    # Overview label
    ax1.text(filter_datetime[0]+timedelta(days = -3), patch_max + temp_buffer, f"C", weight='bold', fontsize = 9, va='top', ha='right')

    # Save to file
    #show()
    savefig(f'../images/figure-2-{canal_id}.png', bbox_inches='tight', dpi = 300)

def draw_figure_3(shading_boolean, canal_list, depth):
    '''
    > Model validation plot, using data from the CRT (https://canalrivertrust.org.uk/)
    > Of the monitoring locations, the following are suitable for comparison:
        - Bridge Pagefield Pipe Crossing (lalc_73): 10 m frequency
        - Bridge 15 Chain Lane Bridge (nabc_2): 10 m frequency
        - New Road Bridge 148 (suc_41): 1 h frequency
        - Bridge 32. Priorswood (batc_14): 15 m frequency (from 03-01-2022)
        - Bridge 4 Wakefield Road Bridge (hc_53): 1 h frequency
        - Bridge 73 Anchor Bridge (cc3_11): 10 m frequency
    > The following locations are less suitable for comparison:
        - Anderton Waste Weir Flow (None): not directly modelled
        - Bridge 46 Benthouse Bridge (None): not directly modelled
        - Hawkesbury Lock (None): not directly modelled
        - Bridge C Lightbody Street (lalc_4): missing measurement data
        - Bridge 40 Taverners Bridge (cc3_53): missing measurement data
    > depth (cm) refers to the upper value
        - 20 = 20 - 40 cm
        - 40 = 40 - 60 cm
        - 60 = 60 - 80 cm
        - 80 = 80 - 100 cm
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
    #axes_labels = ['A','B','C','D','E','F','G','H','I','J','K','L']
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

        # Try
        try:

            # Include shading
            if shading_boolean:
            
                # Load model output record, including shading
                with open(f"../results/including-shading/model-output-{canal_id}.json") as results_path:
                    canal_output = load(results_path)

            # Exclude shading
            else:

                # Load model output record, excluding shading
                with open(f"../results/excluding-shading/model-output-{canal_id}.json") as results_path:
                    canal_output = load(results_path)

        except FileNotFoundError:
            print("No modelled data to compare")
            exit()

        try: 

            # Extract lists of values to plot, converting to celcius
            modelled_water = [x[f'depth_water_{depth}'] - 273.15 for x in canal_output.values()]

        # EAFP: surface layer
        except KeyError:
        
            # Extract lists of values to plot, converting to celcius
            modelled_water = [x[f'surface_water_k'] - 273.15 for x in canal_output.values()]


        # List of dict keys (unix), convert to datetime format
        modelled_dt = list(canal_output.keys())
        modelled_datetime = [datetime.utcfromtimestamp(int(x)) for x in modelled_dt]

        # Init dictionary to store residuals
        residual_dict = {}

        # Iterate through the modelled unix times
        for key, value in canal_output.items():

            # If this unix time exists in the measured dictionaru
            if key in measured_dict:
                
                try: 

                    # Calculate the residual, converting from K to celcius
                    res = (value[f'depth_water_{depth}'] - 273.15) - measured_dict[key]

                # EAFP: surface layer
                except KeyError:

                    # Calculate the residual, converting from K to celcius
                    res = (value[f'surface_water_k'] - 273.15) - measured_dict[key]

                # Add to residual dictionary
                residual_dict[key] = res

        # Extract datetimes and values
        residual_dt = list(residual_dict.keys())
        residual_datetime = [datetime.utcfromtimestamp(int(x)) for x in residual_dt]
        residual_temp = [x for x in residual_dict.values()]

        # Absolute
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
        ax1.text(0.95, 0.95, f"Location = {canal_id}\nLength = {array_length:,}", transform = ax1.transAxes, weight='normal', fontsize = 6, va='top', ha='right', 
                 style='italic', color = "#727573", linespacing = 1.5)
        
        # Add absolute residuals
        ax2.text(0.04, 0.07, f"Absolute median residual = {median(abs_residual):.2f} ({quantile(abs_residual, 0.25):.2f}, {quantile(abs_residual, 0.75):.2f})", 
                 transform = ax2.transAxes, weight='normal', fontsize = 6.5, va='bottom', ha='left', style='italic')

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
        # ax2.text(0.04, 0.80, f"{axes_labels[2 * index + 1]}", transform = ax2.transAxes, weight='bold', fontsize = 13, va='center', ha='left')

    # Save to file
    # show()
    savefig(f'../images/figure-3-{shading_boolean}_{depth}-{depth+20}cm.png', bbox_inches='tight', dpi = 300)

def draw_figure_4(start_date, interval, shading_boolean):
    '''
    > Draw Figure 4 which summarises the results including | excluding shading effects (daytime vs. nighttime), showcasing:
        - (A-C) fluctuations in median air temperature change at weekly, fortnightly, or monthly intervals
        - (B-D) histograms of median air temperature change
    '''

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Try
    try:

        # Load summary file, including | excluding hading
        with open(f"../outputs/feature-values-{folder_path}-shading.json") as path:
                summary_data = load(path)

    # File not found
    except FileNotFoundError:
        exit("File not found")

    # Determine averaging interval
    if interval == "fortnight":
        day_medians = [x['fortnightly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['fortnightly values']['night medians'] for x in summary_data.values()]

    # Extract weekly medians
    elif interval == "week":
        day_medians = [x['weekly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['weekly values']['night medians'] for x in summary_data.values()]

    # Extract monthly medians
    elif interval == "month":
        day_medians = [x['monthly values']['day medians'] for x in summary_data.values()]
        night_medians = [x['monthly values']['night medians'] for x in summary_data.values()]

    # Invalid string
    else:
        print("Invalid time interval")
        exit(1)

    # Transpose data from by feature to by interval
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
    
    # Create list of datetimes for x-axis
    if interval == "fortnight":
        dt = [local_datetime + timedelta(days = 7 * (2 * n - 1)) for n in range(1, len(day_medians)+1)]

    elif interval == "week":
        dt = [local_datetime + timedelta(days = 3 + (7 * n)) for n in range(0, len(day_medians))]

    # Not currently used
    elif interval == "month":
        dt = [local_datetime + relativedelta(months=n) for n in range(0, len(day_medians))]

    # Flatten list of monthly values
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
    for right_ax, left_ax, ax_label, data, flat, label, fill, edge, time_colour in zip([ax1, ax2], [ax3, ax4], [["B", "A"], ["D", "C"]],
                                     [day_medians, night_medians], [flat_day, flat_night], 
                                     ['Day', 'Night'], ["#FFB16F", "#707070"], ["#707070", "#232323"],
                                    [["#FFD0A6", "#FFA453", "#B05300"], ["#B9B9B9", "#656565", "#090909"]]):

        # Extract quantiles
        q1 = [quantile(x, 0.25) for x in data]
        med =  [quantile(x, 0.5) for x in data]
        q3 = [quantile(x, 0.75) for x in data]

        # Tukey whiskers (1.5*iqr)
        tukey_lower = [a - ((c - a) * 1.5) for a,c in zip(q1, q3)]
        tukey_upper = [c + ((c - a) * 1.5) for a,c in zip(q1, q3)]

        # Init lists
        whisker_lower = []
        whisker_upper = []

        # Iterate, and find the last datum less than tukey upper (and v.v. for tukey lower)
        for val, low_limit, upp_limit in zip(data, tukey_lower, tukey_upper):
            whisker_lower.append(min([i for i in val if i > low_limit]))
            whisker_upper.append(max([i for i in val if i < upp_limit]))
        
        # Set bin min-max, with 0.5C buffer
        lower_bins = min(whisker_lower) - 0.5
        upper_bins = max(whisker_upper) + 0.5

        # Y axis limits, 1C buffer
        axis_min, axis_max = lower_bins - 0.5, upper_bins + 0.5

        # Plot histogram
        right_ax.hist(flat, bins=arange(lower_bins, upper_bins + 0.1, 0.1), orientation="horizontal",histtype='stepfilled',
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
            m = [x - 0.1 for x in med]

            # X-axis labels
            months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

            # Horizontal line
            left_ax.axhline(y = 0, color = '#000000', linestyle = 'dashed', linewidth = 1)

            # Plot boxplot
            left_ax.bar(x = months, bottom = whisker_lower, width = 1, height = h, facecolor = time_colour[0], edgecolor = None,
                    label = r'$\genfrac{}{}{0}{}{\it{p}75 + 1.5iqr}{\it{p}25 - 1.5iqr}$')
            left_ax.bar(x = months, bottom = q1, width = 1, height = hq, facecolor = time_colour[1], edgecolor = None,
                    label = "$\it{p}$25 – $\it{p}$75")
            left_ax.hlines(y = m, xmin = arange(-0.5, 11.5, 1), xmax = arange(0.5, 12.5, 1), color = time_colour[2], 
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
            left_ax.plot(dt, med, color = time_colour[2], label = "$\it{p}$50")

            # Add legend
            handles, labels = left_ax.get_legend_handles_labels()
            left_ax.legend(reversed(handles), reversed(labels), frameon=False, ncol = 3, fontsize="8", loc = 'upper right', labelspacing = 0.2,
                        columnspacing=0.5, borderpad=0.2)  

            # Add labels to top subplot
            left_ax.set_ylabel("$\Delta$$\it{PT}$ (°C)")
            left_ax.set_ylim(axis_min, axis_max)

            # Set date formatter
            left_ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))

            # Axis labels
            if interval == "fortnight":
                left_ax.text(-0.03, 1.05, ax_label[1], transform = left_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
            
            # Axis labels
            elif interval == "week":
                left_ax.text(-0.03, 1.05, ax_label[1], transform = left_ax.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

                # Plot heatwave labels, corresponding to:
                # Week 23: heatwave 15-17 June
                # Week 28: heatwave 17-19 July
                # Week 31: heatwave 9-13 Aug

                # Iterate through registered heatwave weeks
                for week in [23, 28, 31]:
                
                    # Week, Q1-Q2-Q3
                    print(f"For {dt[week]}, weekly median air temperature change ({label}) was {med[week]:.2f} (Q1: {q1[week]:.2f}, Q3: {q3[week]:.2f})")

                # Heatwave labels for daytime
                if label == "Day":
                    for x, y, label in zip([dt[23], dt[28], dt[31]], [whisker_lower[23], whisker_lower[28], whisker_lower[31]], ['I', 'II', 'III']):
                        left_ax.text(x, y-0.8, label, fontsize = 8, va='bottom', ha='center')


    # Axes specific labelling
    ax1.set_xlabel(None)
    ax2.set_xlabel("Frequency")
    ax4.set_xlabel("Date")
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=False)
    ax4.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=True)
    
    # Remove last xtick label to simplify plot design
    if interval == "week":
        setp(ax4.get_xticklabels()[-1], visible=False)

    # show()
    savefig(f'../images/figure-4-{interval}.png', bbox_inches='tight', dpi = 300)

def shading_gif(feature_name, start_date, days, delete):
    '''
    > Creates a gif showing shadows cast by buildings for a specific canal feature and duration
    > Image files (png) deleted upon completion
    > 10 fps, 10 loops 
    '''

    # Duration of analysis in hours
    duration = days * 24

    # Read the dataframe directly using pyogrio
    canals_gdf = read_dataframe("../data/canal-geometries-filtered-modified-dissolved-id.shp")

    # Extract chosen canal
    modelled_canal = canals_gdf.loc[canals_gdf['code_id'] == feature_name]

    # If empty
    if modelled_canal.empty:
        
        # Load the Irish dataset, and extract canal
        canals_gdf = read_dataframe("../data/urban-canals-ireland-dissolved-id-tm65.shp")

        # Extract chosen canal
        modelled_canal = canals_gdf.loc[canals_gdf['code_id'] == feature_name]

        # If feature not present, exit function
        if modelled_canal.empty:
            print(f"Canal feature {feature_name} missing from both datasets. Check feature code. Exiting function...")
            return

    # Convert to series, project to 4326 and return lon-lat (perhaps could be simplified)
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

        # Return buildings within a specified distance (500 m) of feature polygon (using id) - this takes 146 seconds for all 2738 features, with 5 cores
        buildings = return_buildings_polygon(feature_tuple, 500, plot = False, epsg = canals_gdf.crs)

        # Count of buildings
        n_buildings = len(buildings)

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

            # If the file exists
            if isfile(filename):   
                # Append to list and skip
                filenames.append(filename)
                local_datetime += timedelta(hours=0.25)

            # Otherwise
            else:             
                
                # Add altitude using dt index, calculate the solar altitude (in degrees) using PySolar
                try:
                    solar_altitudes[datetimes.index(local_datetime)] = get_altitude(latitude, longitude, local_datetime)

                # Missing final datetime for list
                except ValueError:
                    
                    # Append missing
                    datetimes.append(local_datetime)
                    solar_altitudes.append(get_altitude(latitude, longitude, local_datetime))

                # Set up output image
                fig, (my_ax, my_ax2) = subplots(2, 1, figsize=(9, 10), height_ratios=[0.8,0.2])
                #my_ax.title(f"Building shadows at (UTC) {local_datetime}.")

                # Set axis limits
                my_ax.set_xlim([buildings.total_bounds[0],  buildings.total_bounds[2]])
                my_ax.set_ylim([buildings.total_bounds[1],  buildings.total_bounds[3]])

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
                        )
                    
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
    with imageio.get_writer(f'../images/shadows/shadows-for-{feature_tuple.code_id}.gif', mode='I', fps = 10, loop = 10) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

            # Delete file from directory
            if delete:
                remove(filename)

#==================== Supplementary ====================#

def draw_supplementary_figure_A(canal_id):
    '''
    > Function to illustrate the outputs of model-spin up
    > Designed for a feature of the Regents Canal, London (rc2_23) but will work with other inputs although model spin-up will have to be performed
    '''

    # If the model spin-up file already exists
    if isfile(f"../outputs/spin_up_output_{canal_id}.csv"):

        # Load using np
        arr = loadtxt(f"../outputs/spin_up_output_{canal_id}.csv",
                        delimiter=",")
        
        # Extract each column
        water = arr[:, 0]
        energy = arr[:, 1]
        concrete = arr[:, 2]

    # Else, Pperform model spin-up
    else:

        # Read in the GB canals dataframe 
        canals_gdf = read_dataframe("../data/canal-geometries-filtered-modified-dissolved-id.shp")

        # If the canal exists
        if canal_id in canals_gdf['code_id'].values:

            # Filter to selected
            filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]
        
        # Read in the Irish canals dataframe 
        else:
            canals_gdf = read_dataframe( "../data/urban-canals-ireland-dissolved-id-tm65.shp")

            # If the canal exists
            if canal_id in canals_gdf['code_id'].values:

                # Filter to selected
                filtered_canals = canals_gdf.loc[canals_gdf['code_id'].isin([canal_id])]

            # Exit function
            else:
                print(f"Canal {canal_id} is MISSING from both datasets")
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

        # Init list of starting value for concrete (K = 282), comprising 10 depth layers [0 - .5m]
        concrete_temp = [MEAN_KELVIN] * 10

        # Model spin up (30 iterations), based on a composite climate record for 2021-12-15, return full record
        _, _, _, water, energy, concrete = model_spin_up(water_temp, water_energy, concrete_temp, split_points, feature_area,
                                                                latitude, longitude, 
                                                                "../outputs/spin_up_climate.json", 30, False, canal_id)
    
    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig, (ax1, ax2, ax3) = subplots(1, 3, figsize=(8, 4), layout='compressed')

    # Plot the water temperature
    ax1.plot(water, color = '#599FFF')
    ax2.plot(concrete, color = '#FFB16F')
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

def draw_supplementary_figure_B(shading_boolean):

    # Include or exclude shading effects, updates path to file
    if shading_boolean:
        folder_path = "including"
    else:
        folder_path = "excluding"

    # Load summary file, including | excluding hading
    try:
        with open(f"../outputs/feature-values-{folder_path}-shading.json") as summary_path:
                summary = load(summary_path)

    # File not found
    except FileNotFoundError:
        exit("File not found")

    # Circle areas (m^2) for illustration
    area_one, area_two = 750, 5000

    # Extract values
    day_medians = array([x['day values']['median'] for x in summary.values()])
    area = array([x['geometry']['area'] for x in summary.values()])
    buffer_area = array([x['geometry']['buffer area'] for x in summary.values()])

    # Calculate radius and draw circle
    radius_one = sqrt(area_one / pi)
    circle_one = Point(0, 0).buffer(radius_one, resolution = 10000)
    buffer_one = circle_one.buffer(50)
    buffer_one_gdf = GeoSeries(buffer_one)
    circle_one_gdf = GeoSeries(circle_one)

    # Calculate radius and draw circle
    radius_two = sqrt(area_two / pi)
    circle_two = Point(200, 0).buffer(radius_two, resolution = 10000)
    buffer_two = circle_two.buffer(50)
    buffer_two_gdf = GeoSeries(buffer_two)
    circle_two_gdf = GeoSeries(circle_two)

    # Ratio of canal area to buffer area
    ratio = median(buffer_area/area)

    # Ratio of canal area to buffer area
    area_ratio = buffer_area / area

    # Update params
    rcParams.update({'font.size': 9,
                     "mathtext.fontset" : "dejavuserif",
                     'mathtext.default' : 'regular'})

    # Set up output image
    fig = figure(layout='compressed', figsize=(6, 3))
    gs = fig.add_gridspec(1,2, width_ratios = [1,2])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Plot the wT vs. ratio
    ax1.scatter(area_ratio, day_medians, marker = ".", color = '#707070', edgecolors = 'none', s = 7.5, alpha = 0.5)

    # Set aspect ratio to 1
    ratio = 1.0
    x_left, x_right = ax1.get_xlim()
    y_low, y_high = ax1.get_ylim()
    ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    
    # Add shapely circles and buffer areas
    buffer_two_gdf.plot(ax = ax2, facecolor = '#FFFFFF', edgecolor = '#343434', linewidth = 1, linestyle = '--')
    circle_two_gdf.plot(ax = ax2, color = '#C8F1FF', edgecolor = '#249AC3')
    buffer_one_gdf.plot(ax = ax2, facecolor = '#FFFFFF', edgecolor = '#343434', linewidth = 1, linestyle = '--')
    circle_one_gdf.plot(ax = ax2, color = '#C8F1FF', edgecolor = '#249AC3')

    # Add custom legends
    canal_polygon = mpatches.Patch(facecolor='#C8F1FF', edgecolor = '#249AC3', linewidth = 1, label='Normalised canal area')
    modelled_polygon = mpatches.Patch(facecolor = "none", edgecolor = '#343434', linewidth = 1, linestyle = '--', label='Normalised buffer area ($\Delta$$\it{PT)}$')
    ax2_legend = ax2.legend(handles=[canal_polygon, modelled_polygon], loc = 'lower center', ncols = 1,
            bbox_to_anchor = (0.5, -0.3),
            columnspacing = 0.6,
            labelspacing = 0.6,
            framealpha = 0)
    ax2_legend.get_frame().set_linewidth(0.0)
    ax2.add_artist(ax2_legend)

    # Add scale bar
    ax2.add_artist(ScaleBar(dx=0.8, units="m", location="upper left", length_fraction=0.2, 
                            frameon=False, box_alpha=1, font_properties={"size": 8}))

    # x-axis labelling
    ax1.set_xlabel("Area/buffer ratio")

    # y-axis labelling
    ax1.set_ylabel("$\Delta$$\it{PT}$ (°C)")

    # Tick params
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=True)
    ax2.tick_params(bottom=False, top=False, left=False, right=False, direction="in", labelbottom=False, labelleft=False)

    # Axis labels
    ax1.text(-0.05, 1.05, f"A", transform = ax1.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax2.text(-0.02, 1.05, f"B", transform = ax2.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

    # Add text representing ratios
    ax2.text(0, 25, f"R={(buffer_one.area / circle_one.area):.1f}", fontsize = 8, va='center', ha='center')
    ax2.text(200, 50, f"R={(buffer_two.area / circle_two.area):.1f}", fontsize = 8, va='center', ha='center')

    # Save
    # show()
    savefig(f'../images/supp-figure-B.png', bbox_inches='tight', dpi = 300)

def draw_supplementary_figure_C():

    # If the file already exists
    if isfile("../outputs/area-calculations.json"):

        # Load
        with open("../outputs/area-calculations.json") as results_path:
            results_dict = load(results_path)

    # Perform the analysis
    else:

        # Files for Britain and Ireland
        file_names = ["../data/canal-geometries-filtered-modified-dissolved-id.shp", 
                    "../../national-canal-cooling-data/ireland/urban-canals-ireland-dissolved-id-tm65.shp"]
        
        # Init result dict
        results_dict = {}
        
        # Iterate
        for file in file_names:

            # Read the dataframe directly using pyogrio
            canals_gdf = read_dataframe(file)

            # Iterate through canal features
            for feature_tuple in canals_gdf.itertuples():

                # Skip canals with an area < 500 m^2
                if feature_tuple.geometry.area < 500 or feature_tuple.geometry.area > 200_000:
                    continue

                # For canals > 500 m^2
                else:
                    
                    # Returns shapely area (27700 | 29902)
                    feature_area = feature_tuple.geometry.area

                    # Calculate radius and draw circle
                    radius = sqrt(feature_area / pi)
                    circle = Point(0, 0).buffer(radius, resolution = 10000)

                    # Update results dict
                    results_dict.update({feature_tuple.code_id : {'area' : feature_area,
                                                                '50_m buffer' : feature_tuple.geometry.buffer(50).area,
                                                                '50_m circle buffer' : circle.buffer(50).area}})


        # Write results to file
        dump(results_dict, open(f"../outputs/area-calculations.json", 'w'))

    # Extract lists of values to plot, converting to celcius
    area = [x['area'] for x in results_dict.values()]
    buffer_50m = [x['50_m buffer'] for x in results_dict.values()]
    circle_50m = [x['50_m circle buffer'] for x in results_dict.values()]

    # Combine into dataframe
    df = DataFrame(list(zip(area, buffer_50m, circle_50m)),
               columns =['area', '50_m buffer', '50_m circle buffer'])

    # Sort by area (ascending)
    df.sort_values('area')

    # Set global font size
    rcParams.update({'font.size': 9})

    # Set up output image
    fig, (ax1, ax2) = subplots(1, 2, figsize=(6, 4), layout='compressed')

    # Plot the scatter 
    ax1.scatter(df['area'], df['50_m buffer'], marker = ".",
                color = '#599FFF',
                edgecolors = 'none',
                s = 10,
                alpha = 0.75)
    ax2.scatter(df['area'], df['50_m circle buffer'], marker = ".",
                color = '#FFB16F',
                edgecolors = 'none',
                s = 10,
                alpha = 0.75)
    
    '''
    # Lin reg. 
    b1, a1 = polyfit(area, buffer_50m, deg=1)
    b2, a2 = polyfit(area, circle_50m, deg=1)

    # Create sequence of 100 numbers from 0 to 100 
    xseq = linspace(min(area), max(area), num=1000)

    # Plot lin reg.
    ax1.plot(xseq, a1 + b1 * xseq, color="k", lw=0.5, ls = '--')
    ax2.plot(xseq, a2 + b2 * xseq, color="k", lw=0.5, ls = '--')
    '''

    #set aspect ratio to 1
    ratio = 1.0
    x_left, x_right = ax1.get_xlim()
    y_low, y_high = ax1.get_ylim()
    ax1.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

    '''
    # Add 1:1 lines
    ax1.axline((0, 0), (max(area), max(area)), linewidth=0.5, color='k', linestyle='--')
    ax2.axline((0, 0), (max(area), max(area)), linewidth=0.5, color='k', linestyle='--')
    '''
    
    # x-axis labelling
    ax1.set_xlabel(f"Canal area ($m^2$)")
    ax2.set_xlabel(f"Canal area ($m^2$)")

    # y-axis labelling
    ax1.set_ylabel(f"Modelled area ($m^2$)")
    ax2.set_ylabel(f"Modelled area ($m^2$)")

    # Tick params
    ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=True)
    ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelbottom=True)

    # Axis labels
    ax1.text(-0.05, 1.05, f"A", transform = ax1.transAxes, weight='bold', fontsize = 12, va='top', ha='right')
    ax2.text(-0.05, 1.05, f"B", transform = ax2.transAxes, weight='bold', fontsize = 12, va='top', ha='right')

    # Save
    #show()
    savefig(f'../images/supp-figure-C.png', bbox_inches='tight', dpi = 300)   

def draw_supplementary_figure_D(feature_name):
    '''
    > Figure to compare water and concrete temperatures when shading is included | excluded, for a specific feature
    '''
    # Try and load model output record, including | excluding shading
    try:
        with open(f"../results/excluding-shading/model-output-{feature_name}.json") as excluding_path:
            excluding_data = load(excluding_path)

        with open(f"../results/including-shading/model-output-{feature_name}.json") as including_path:
            including_data = load(including_path)

    # File not found
    except FileNotFoundError:
        exit(f"One or more results files is MISSING for feature {feature_name}. Exiting function...")

    # Extract water temperature values, excluding | including shading effects
    excluding_water = [x['surface_water_k'] - 273.15 for x in excluding_data.values()]
    including_water = [x['surface_water_k'] - 273.15 for x in including_data.values()]

     # Extract water temperature values, excluding | including shading effects
    excluding_concrete = [x['surface_concrete_k'] - 273.15 for x in excluding_data.values()]
    including_concrete = [x['surface_concrete_k'] - 273.15 for x in including_data.values()]

    # Extract differences
    diff_water = [e - i  for e, i in zip(excluding_water, including_water)]  
    diff_concrete = [e - i  for e, i in zip(excluding_concrete, including_concrete)]  

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

    # Plot the concrete temperature including | excluding shading
    ax2.plot(dt_datetime, excluding_concrete, label = 'Excluding shading', color = "#FFB16F", alpha=0.75)
    ax2.plot(dt_datetime, including_concrete, label = 'Including shading', color = "#707070", alpha=0.75)
    ax2.set_xlabel("Month")
    ax2.set_ylabel("Concrete temperature (°C)", labelpad=-4)
    ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=True, labelbottom=True, labeltop=False)
    ax2.legend(frameon=False)

    # Plot the concrete temperature including | excluding shading
    ax3.plot(dt_datetime, diff_water, label = None, color = "#707070", alpha=0.75)
    ax3.set_xlabel(None)
    ax3.set_ylabel("Temperature difference (°C)")
    ax3.yaxis.set_label_position("right")
    ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", labelleft=False, labelbottom=False, labeltop=False, labelright=True)

    # Plot the concrete temperature including | excluding shading
    ax4.plot(dt_datetime, diff_concrete, label = None, color = "#707070", alpha=0.75)
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
    savefig(f'../images/supp-figure-D.png', bbox_inches='tight', dpi = 300)

# The script is meant to be run
if __name__ == '__main__':
    
    # Run main function
    main()

