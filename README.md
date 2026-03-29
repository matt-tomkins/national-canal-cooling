# The cooling power of urban waterways
A University of Manchester project for modelling the cooling effect of urban canals across Great Britain and Ireland.
Project funded by the [Canal and River Trust](https://canalrivertrust.org.uk/).

This project expands upon on an earlier local-scale project, which can be found [here](https://github.com/jonnyhuck/canal-cooling).


## Repository structure:

This repository contains the following directories:

- `/code` which contains the key Python files used to run the model, including:
    - `national_canal_cooling` which can be run directly.
    - `functions` which contains all key functions required by `national_canal_cooling`.
    - `params` which contains key parameters, including designation and selection of the reference material, and which is also called by `national_canal_cooling`.
    - `analysis` for generating results and plots. 
    - `shading` which is not called directly, but is used to generate the shading proportion of each canal section, adapted from [pybdshadow](https://pybdshadow.readthedocs.io/en/latest/index.html).
- `/images` which contains the figures produced for the main text and supplementary information.
- `/outputs` which contains key outputs from `analysis.py`.
- `/tables` which contains Table 1 (*Summary ΔPT values*) and Supplementary Table 1 (*Water temperature data from the Canal and River Trust*)
- `/test` which contains test data and outputs.

## Data availability

**Note** that the model is run using weather data from [OpenWeather](https://openweathermap.org/) and canal and building geometries from the [Ordnance Survey MasterMap Topography layer](https://www.ordnancesurvey.co.uk/products/os-mastermap-topography-layer), which it is not possible to provide in full. 

Instead, the **test** directory contains simulated canal and building geometries and weather data[^1], which allow the user to run the model and evaluate its performance. This is the default behaviour in `national_canal_cooling` as the test data (`../test/test-canals.shp`) is included as one of the input files. 

[^1]: Air temperature (K), humidity (%), cloud cover (%), wind speed (m/s) and air pressure (hPa)

## Project history:

- Mathematical model developed by Harry Mcdonald during his MSc at The University of Manchester, supervised by Dr. Joanne Tippett
- Model adapted to a local-scale spatial model by Sophie Taylor during her MSc at The University of Manchester, supervised by Dr. Jonny Huck
- Model refined and upscaled to consider all urban waterways in the UK by Dr. Matt Tomkins and Dr. Jonny Huck

## Current roles (alphabetical order):

- Dr. Jonny Huck [(Professor in Computational Geography, Department of Geography, The University of Manchester)](https://research.manchester.ac.uk/en/persons/jonathan.huck)
- Dr. Harry Mcdonald [(Research Associate, Department of Electrical & Electronic Engineering, The University of Manchester)](https://research.manchester.ac.uk/en/persons/harry.mcdonald)
- Dr. Joanne Tippett [(Lecturer in Spatial Planning, Department of Planning and Environmental Management, The University of Manchester)](https://research.manchester.ac.uk/en/persons/joanne.tippett)
- Dr. Matt Tomkins [(Lecturer in GIS, Department of Geography, The University of Manchester)](https://research.manchester.ac.uk/en/persons/matthew.tomkins)
