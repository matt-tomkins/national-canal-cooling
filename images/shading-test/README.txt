Shading verification (vs. SunCalc)

Feature: 'ac1_6' (Ashton Canal #6)
Datetime: 2022-01-01 09-45-00
Function: return_shading() in auxiliary-code.py

Output files:
- 'shadows-for-ac1_6-at-01_01-09_45_00.png'
- 'shadows-around-ac1_6-at-01_01-09_45_00.shp'
- 'buildings-around-ac1_6.shp'

Key packages: pybdshadow, pysolar

Results summary ---------------------------------------------------------------------
Most westerly highlight building (~385390, ~398090) has a height of 27.5 m ('relhmax')
At the above datetime, this produces shadows of length ~215 m
Calculation assumes:
 - solar azimuth = 7.3°
 - solar angle = 146.2°

Comparison (SunCalc) ------------------------------------------------------
Source: https://www.suncalc.org/#/53.4795,-2.222,17/2022.01.01/09:45/27.5/3

This gives identical results:
 - solar azimuth = 7.27°
 - solar angle = 146.19°
 - shadow length = 215.49 m 

