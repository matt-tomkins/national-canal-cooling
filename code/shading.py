'''
> Obtained pybdshadow.py, utils.py and preprocess.py from GitHub (23-06-2023: https://github.com/ni1o1/pybdshadow)
> Changelog:
    - The original code uses 'suncalc' for calculating solar altitude-azimuth (https://pypi.org/project/suncalc/)
    - However, while the azimuth values are "correct", the altitude values differ somewhat from other sources (e.g. https://www.suncalc.org/, https://keisan.casio.com/exec/system/1224682277)
    [1] It has now been modified to use 'PySolar' (https://pysolar.readthedocs.io/en/latest/)
    - This returns values which are consistent with the above sources, and produces shadow lengths that match expected shadow lengths i.e. height / tan(rad(angle))
    - Only superficial changes have made to the remainder of the code (structure / style)
'''

"""
BSD 3-Clause License

Copyright (c) 2022, Qing Yu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Required libraries
import math
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
from suncalc import get_position
from shapely.geometry import Polygon, LineString, MultiPolygon

# For use of PySolar
from pysolar.solar import get_altitude, get_azimuth
from math import radians, degrees

def modify_azimuth(az):
    '''
    > This function takes in azimuths relative to north (PySolar) and converts to suncalc version.
    > Confusingly, these are measured relative to south (east = negative, west = positive)
    > Azimuth returned as radians
    '''

    # First convert to degrees for simplicity
    az_degrees = degrees(az)

    # For westerly directions
    if az_degrees > 180:

        # Convert to positive, relative to south
        return radians(az_degrees - 180)
    
    # For easterly directions
    else: 
    
        # COnvert to negative, relative to south
        return radians(-abs(180 - az_degrees))
    

def bd_preprocess(buildings, height=''):
    '''
    Preprocess building data, so that we can perform shadow calculation.
    Remove empty polygons and convert multipolygons into polygons.
    Parameters
    --------------
    buildings : GeoDataFrame
        Buildings.
    height : string
        Column name of building height(meter).
    Return
    ----------
    allbds : GeoDataFrame
        Polygon buildings
    '''
    buildings['geometry'] = buildings.buffer(0)
    buildings = buildings[buildings.is_valid].copy()
    if height!='':
        # 建筑高度筛选
        buildings[height] = pd.to_numeric(buildings[height], errors='coerce')
        buildings = buildings[buildings[height]>0].copy()

    polygon_buildings = buildings[buildings['geometry'].apply(
        lambda r:type(r) == shapely.geometry.polygon.Polygon)]
    multipolygon_buildings = buildings[buildings['geometry'].apply(
        lambda r:type(r) == shapely.geometry.multipolygon.MultiPolygon)]
    allbds = []
    for j in range(len(multipolygon_buildings)):
        r = multipolygon_buildings.iloc[j]
        singlebd = gpd.GeoDataFrame()
        singlebd['geometry'] = list(r['geometry'].geoms)
        for i in r.index:
            if i != 'geometry':
                singlebd[i] = r[i]
        allbds.append(singlebd)
    allbds.append(polygon_buildings)
    allbds = pd.concat(allbds)
    if len(allbds) > 0:
        allbds = gpd.GeoDataFrame(allbds)
        allbds['building_id'] = range(len(allbds))
        allbds['geometry'] = allbds.buffer(0)
    else:
        allbds = gpd.GeoDataFrame()
    return allbds

def gdf_difference(gdf_a,gdf_b,col = 'building_id'):
    '''
    difference gdf_b from gdf_a
    '''
    gdfa = gdf_a.copy()
    gdfb = gdf_b.copy()
    gdfb = gdfb[['geometry']]
    #判断重叠
    from shapely.geometry import  MultiPolygon
    gdfa.crs = gdfb.crs
    gdfb = gpd.sjoin(gdfb,gdfa).groupby([col])['geometry'].apply(
            lambda df: MultiPolygon(list(df)).buffer(0)).reset_index()
    #分割有重叠和无重叠的
    gdfb['tmp'] = 1
    gdfa_1 = pd.merge(gdfa,gdfb[[col,'tmp']],how = 'left')
    gdfa = gdfa_1[gdfa_1['tmp'] == 1].drop('tmp',axis = 1)
    gdfa_notintersected = gdfa_1[gdfa_1['tmp'].isnull()].drop('tmp',axis = 1)
    #对有重叠的进行裁剪
    gdfa = gdfa.sort_values(by = col).set_index(col)
    gdfb = gdfb.sort_values(by = col).set_index(col)
    gdfa.crs = gdfb.crs
    gdfa['geometry'] = gdfa.difference(gdfb).buffer(0)
    gdfa = gdfa.reset_index()
    #拼合
    gdfa = pd.concat([gdfa,gdfa_notintersected])
    return gdfa

def gdf_intersect(gdf_a,gdf_b,col = 'building_id'):
    '''
    intersect gdf_b from gdf_a
    '''
    gdfa = gdf_a.copy()
    gdfb = gdf_b.copy()
    gdfb = gdfb[['geometry']]
    #判断重叠
    from shapely.geometry import  MultiPolygon
    gdfa.crs = gdfb.crs
    gdfb = gpd.sjoin(gdfb,gdfa).groupby([col])['geometry'].apply(
            lambda df: MultiPolygon(list(df)).buffer(0)).reset_index()
    #分割有重叠和无重叠的
    gdfb['tmp'] = 1
    gdfa_1 = pd.merge(gdfa,gdfb[[col,'tmp']],how = 'left')
    gdfa = gdfa_1[gdfa_1['tmp'] == 1].drop('tmp',axis = 1)
    #对有重叠的进行裁剪
    gdfa = gdfa.sort_values(by = col).set_index(col)
    gdfb = gdfb.sort_values(by = col).set_index(col)
    gdfa.crs = gdfb.crs
    gdfa['geometry'] = gdfa.intersection(gdfb).buffer(0)
    gdfa = gdfa.reset_index()

    return gdfa

def lonlat2aeqd(lonlat):
    '''
    Convert longitude and latitude to azimuthal equidistant projection coordinates.
    Parameters
    ----------
    lonlat : numpy.ndarray
        Longitude and latitude in degrees. The shape of the array is (n,m,2), where n and m are the number of pixels in the first and second dimension, respectively. The last dimension is for longitude and latitude.
    Returns
    -------
    proj_coords : numpy.ndarray
        Azimuthal equidistant projection coordinates. The shape of the array is (n,m,2), where n and m are the number of pixels in the first and second dimension, respectively. The last dimension is for x and y coordinates.
    example
    -----------------
    >>> import numpy as np
    >>> from pybdshadow import utils
    >>> lonlat = np.array([[[120,30],[121,31]],[[120,30],[121,31]]])
    >>> proj_coords = utils.lonlat2aeqd(lonlat)
    >>> proj_coords
    array([[[-48243.5939812 , -55322.02388971],
            [ 47752.57582735,  55538.86412435]],
           [[-48243.5939812 , -55322.02388971],
            [ 47752.57582735,  55538.86412435]]])
    '''
    meanlon = lonlat[:,:,0].mean()
    meanlat = lonlat[:,:,1].mean()
    from pyproj import CRS
    epsg = CRS.from_proj4("+proj=aeqd +lat_0="+str(meanlat)+" +lon_0="+str(meanlon)+" +datum=WGS84")
    from pyproj import Transformer
    transformer = Transformer.from_crs("EPSG:4326", epsg,always_xy = True)
    proj_coords = transformer.transform(lonlat[:,:,0], lonlat[:,:,1])
    proj_coords = np.array(proj_coords).transpose([1,2,0])
    return proj_coords


def aeqd2lonlat(proj_coords,meanlon,meanlat):
    '''
    Convert azimuthal equidistant projection coordinates to longitude and latitude.
    Parameters
    ----------
    proj_coords : numpy.ndarray
        Azimuthal equidistant projection coordinates. The shape of the array is (n,m,2), where n and m are the number of pixels in the first and second dimension, respectively. The last dimension is for x and y coordinates.
    meanlon : float
        Longitude of the center of the azimuthal equidistant projection in degrees.
    meanlat : float
        Latitude of the center of the azimuthal equidistant projection in degrees.
    Returns
    -------
    lonlat : numpy.ndarray
        Longitude and latitude in degrees. The shape of the array is (n,m,2), where n and m are the number of pixels in the first and second dimension, respectively. The last dimension is for longitude and latitude.
    Example
    -----------------
    >>> import numpy as np
    >>> from pybdshadow import utils
    >>> proj_coords = proj_coords = np.array(
        [[[-48243.5939812 , -55322.02388971],
          [ 47752.57582735,  55538.86412435]],
         [[-48243.5939812 , -55322.02388971],
          [ 47752.57582735,  55538.86412435]]])
    >>> lonlat = utils.aeqd2lonlat(proj_coords,120.5,30.5)
    >>> lonlat
    array([[[120.,  30.],
            [121.,  31.]],
           [[120.,  30.],
            [121.,  31.]]])
    '''
    from pyproj import CRS
    epsg = CRS.from_proj4("+proj=aeqd +lat_0="+str(meanlat)+" +lon_0="+str(meanlon)+" +datum=WGS84")
    from pyproj import Transformer
    transformer = Transformer.from_crs( epsg,"EPSG:4326",always_xy = True)
    lonlat = transformer.transform(proj_coords[:,:,0], proj_coords[:,:,1])
    lonlat = np.array(lonlat).transpose([1,2,0])
    return lonlat


def calSunShadow_vector(shape, shapeHeight, sunPosition):
    '''
    Calculate the shadow of a building on the ground.

    Parameters
    ----------
    shape : numpy.ndarray
        The shape of the building. The shape of the array is (n,2,2), where n the number of walls, 2 is that each wall has two points, and the last dimension is for longitude and latitude.
    shapeHeight : float
        The height of the building.
    sunPosition : dict
        The position of the sun. The keys are 'azimuth' and 'altitude'.

    Returns
    -------
    shadow : numpy.ndarray
        The shadow of the building on the ground. shape = [n,5,2]
    '''
    # transform coordinate system
    meanlon = shape[:,:,0].mean()
    meanlat = shape[:,:,1].mean()
    shape = lonlat2aeqd(shape)

    azimuth = sunPosition['azimuth']
    altitude = sunPosition['altitude']

    n = np.shape(shape)[0]
    distance = shapeHeight/math.tan(altitude)

    # calculate the offset of the projection position
    lonDistance = distance*math.sin(azimuth) 
    lonDistance = lonDistance.reshape((n, 1))
    latDistance = distance*math.cos(azimuth)
    latDistance = latDistance.reshape((n, 1))

    shadowShape = np.zeros((n, 5, 2)) # n buildings, each building has 5 points, each point has 2 dimensions

    shadowShape[:, 0:2, :] += shape  
    shadowShape[:, 2:4, 0] = shape[:, :, 0] + lonDistance
    shadowShape[:, 2:4, 1] = shape[:, :, 1] + latDistance

    shadowShape[:, [2, 3], :] = shadowShape[:, [3, 2], :]
    shadowShape[:, 4, :] = shadowShape[:, 0, :]

    shadowShape = aeqd2lonlat(shadowShape,meanlon,meanlat)
    return shadowShape


def bdshadow_sunlight(buildings, date,  height='height', roof=False,include_building = True,ground=0):
    '''
    Calculate the sunlight shadow of the buildings.

    Parameters
    ----------
    buildings : GeoDataFrame
        Buildings. coordinate system should be WGS84
    date : datetime
        Datetime
    height : string
        Column name of building height(meter).
    roof : bool
        Whether to calculate the roof shadows.
    include_building : bool
        Whether the shadow include building outline.
    ground : number
        Height of the ground(meter).

    Returns
    ----------
    shadows : GeoDataFrame
        Building shadow
    '''

    # Create a copy of the buildings
    building = buildings.copy()

    # Incorporate topographic elevation
    building[height] -= ground
    building = building[building[height] > 0]

    # Return building average coordinates
    lon1, lat1, lon2, lat2 = list(building.bounds.mean())
    lon = (lon1+lon2)/2
    lat = (lat1+lat2)/2

    # Original approach to obtain solar altitude-azimuth, using suncalc
    #sunPosition = get_position(date, lon, lat)

    #print(sunPosition)

    # Init dict
    sunPosition = {}

    # Obtain sun position using PySolar, converting degrees to radians 
    sunPosition['altitude'] = radians(get_altitude(lat, lon, date))
    sunPosition['azimuth'] = radians(get_azimuth(lat, lon, date))

    # Convert to Suncalc direction
    sunPosition['azimuth'] = modify_azimuth(sunPosition['azimuth'])

    # Print, for checking
    # print(date, get_altitude(lat, lon, date), get_azimuth(lat, lon, date))


    # Checks the solar altitude is above the horizon
    if (sunPosition['altitude'] < 0):
        raise ValueError("Given time before sunrise or after sunset") 
    buildingshadow = building.copy()

    a = buildingshadow['geometry'].apply(lambda r: list(r.exterior.coords))
    buildingshadow['wall'] = a
    buildingshadow = buildingshadow.set_index(['building_id'])
    a = buildingshadow.apply(lambda x: pd.Series(x['wall']), axis=1).unstack()
    walls = a[- a.isnull()].reset_index().sort_values(
        by=['building_id', 'level_0'])
    walls = pd.merge(walls, buildingshadow['height'].reset_index())
    walls['x1'] = walls[0].apply(lambda r: r[0])
    walls['y1'] = walls[0].apply(lambda r: r[1])
    walls['x2'] = walls['x1'].shift(-1)
    walls['y2'] = walls['y1'].shift(-1)
    walls = walls[walls['building_id'] == walls['building_id'].shift(-1)]
    walls = walls[['x1', 'y1', 'x2', 'y2', 'building_id', 'height']]
    walls['wall'] = walls.apply(lambda r: [[r['x1'], r['y1']],
                                           [r['x2'], r['y2']]], axis=1)

    ground_shadow = walls.copy()
    walls_shape = np.array(list(ground_shadow['wall']))

    # calculate shadow for walls
    shadowShape = calSunShadow_vector(
        walls_shape, ground_shadow['height'].values, sunPosition)

    ground_shadow['geometry'] = list(shadowShape)
    ground_shadow['geometry'] = ground_shadow['geometry'].apply(
        lambda r: Polygon(r))
    ground_shadow = gpd.GeoDataFrame(ground_shadow)



    ground_shadow = pd.concat([ground_shadow, building])
    ground_shadow = ground_shadow.groupby(['building_id'])['geometry'].apply(
        lambda df: MultiPolygon(list(df)).buffer(0)).reset_index()
    
    ground_shadow['height'] = 0
    ground_shadow['type'] = 'ground'

    if not roof:
        if not include_building:
            #从地面阴影裁剪建筑轮廓
            ground_shadow = gdf_difference(ground_shadow,buildings)
        return ground_shadow
    else:
        def calwall_shadow(walldata, building):
            walls = walldata.copy()
            walls_shape = np.array(list(walls['wall']))
            # calculate shadow for walls
            shadowShape = calSunShadow_vector(
                walls_shape, walls['height'].values, sunPosition)
            walls['geometry'] = list(shadowShape)
            walls['geometry'] = walls['geometry'].apply(lambda r: Polygon(r))
            walls = gpd.GeoDataFrame(walls)
            walls = pd.concat([walls, building])

            walls = walls.groupby(['building_id'])['geometry'].apply(
                lambda df: MultiPolygon(list(df)).buffer(0)).reset_index()
            return walls

        # 计算屋顶阴影
        roof_shadows = []
        for roof_height in walls[height].drop_duplicates():
            # 高于给定高度的墙
            walls_high = walls[walls[height] > roof_height].copy()
            if len(walls_high) == 0:
                continue
            walls_high[height] -= roof_height
            # 高于给定高度的建筑
            building_high = building[building[height] > roof_height].copy()
            if len(building_high) == 0:
                continue
            building_high[height] -= roof_height
            # 所有建筑在此高度的阴影
            building_shadow_height = calwall_shadow(walls_high, building_high)
            # 在此高度的建筑屋顶
            building_roof = building[building[height] == roof_height].copy()
            building_shadow_height.crs = building_roof.crs
            # 取有遮挡的阴影
            building_shadow_height = gpd.sjoin(
                gpd.GeoDataFrame(building_shadow_height), gpd.GeoDataFrame(building_roof))
            if len(building_shadow_height) == 0:
                continue
            # 与屋顶做交集
            building_roof = gdf_intersect(building_roof,building_shadow_height)

            # 再减去这个高度以上的建筑
            building_higher = building[building[height] > roof_height].copy()
            building_roof = gdf_difference(building_roof,building_higher)
            
            #给出高度信息
            building_roof['height'] = roof_height
            building_roof = building_roof[-building_roof['geometry'].is_empty]

            roof_shadows.append(building_roof)
        if len(roof_shadows) == 0:
            roof_shadow = gpd.GeoDataFrame()
        else:
            roof_shadow = pd.concat(roof_shadows)[
                ['height', 'building_id', 'geometry']]
            roof_shadow['type'] = 'roof'

        if not include_building:
            #从地面阴影裁剪建筑轮廓
            ground_shadow = gdf_difference(ground_shadow,buildings)
        
        shadows = pd.concat([roof_shadow, ground_shadow])
        shadows.crs = None
        shadows['geometry'] = shadows.buffer(0.000001).buffer(-0.000001)
        return shadows


def calPointLightShadow_vector(shape, shapeHeight, pointLight):
    '''
    calculate shadow for a point light
    
    Parameters
    ----------
    shape : numpy.array
        The shape of the building. The shape of the array is (n,2,2), where n the number of walls, 2 is that each wall has two points, and the last dimension is for longitude and latitude.
    shapeHeight : numpy.array
        height of building, shape = [n,1], n is the number of buildings
    pointLight : dict
        point light, pointLight = {'position':[lon,lat,height]}
    
    Returns
    -------
    shadowShape : numpy.array
        shape of shadow, shape = [n,5,2]
    '''
    # 多维数据类型：numpy
    # 输入的shape是一个矩阵（n*2*2) n个建筑物面，每个建筑有2个点，每个点有三个维度
    # shapeHeight(n) 每一栋建筑的高度都是一样的
    n = np.shape(shape)[0]
    pointLightPosition = pointLight['position']  # [lon,lat,height]

    # 高度比
    diff = pointLightPosition[2] - shapeHeight
    scale = np.zeros(n)
    scale[diff != 0] = shapeHeight[diff != 0]/(diff[diff != 0])
    scale[scale <= 0] = 10  # n
    scale = scale.reshape((n, 1))

    shadowShape = np.zeros((n, 5, 2))

    shadowShape[:, 0:2, :] += shape  # 前两个点不变
    vertexToLightVector = shape - pointLightPosition[0:2]  # n,2,2

    shadowShape[:, 2, :] = shape[:, 1, :] + \
        vertexToLightVector[:, 1, :]*scale  # [n,2,2] = [n,2,2]+[n,2,2]*n
    shadowShape[:, 3, :] = shape[:, 0, :] + \
        vertexToLightVector[:, 0, :]*scale

    shadowShape[:, 4, :] = shadowShape[:, 0, :] 

    return shadowShape

def bdshadow_pointlight(buildings,
                        pointlon,
                        pointlat,
                        pointheight,
                        merge=True,
                        height='height',
                        ground=0):
    '''
    Calculate the sunlight shadow of the buildings.

    Parameters
    --------------------
    buildings : GeoDataFrame
        Buildings. coordinate system should be WGS84
    pointlon,pointlat,pointheight : float
        Point light coordinates and height(meter).
    date : datetime
        Datetime
    merge : bool
        Whether to merge the wall shadows into the building shadows
    height : string
        Column name of building height(meter).
    ground : number
        Height of the ground
    
    Returns
    ----------
    shadows : GeoDataFrame
        Building shadow
    '''

    building = buildings.copy()

    building[height] -= ground
    building = building[building[height] > 0]

    if len(building) == 0:
        walls = gpd.GeoDataFrame()
        walls['geometry'] = []
        walls['building_id'] = []
        return walls
    # building to walls
    buildingshadow = building.copy()

    a = buildingshadow['geometry'].apply(lambda r: list(r.exterior.coords))  #裸格式的几何
    buildingshadow['wall'] = a  #
    #print(a[0])
    buildingshadow = buildingshadow.set_index(['building_id']) #设置阴影所对应的id
    a = buildingshadow.apply(lambda x: pd.Series(x['wall']), axis=1).unstack()  #压缩为一个数组
    walls = a[- a.isnull()].reset_index().sort_values(
        by=['building_id', 'level_0'])  #重新排序
    walls = pd.merge(walls, buildingshadow['height'].reset_index())#与高度融合
    
    walls['x1'] = walls[0].apply(lambda r: r[0])  #
    walls['y1'] = walls[0].apply(lambda r: r[1])
    walls['x2'] = walls['x1'].shift(-1)  #向量中的序号全部向前提了一个
    walls['y2'] = walls['y1'].shift(-1)
    walls = walls[walls['building_id'] == walls['building_id'].shift(-1)]
    walls = walls[['x1', 'y1', 'x2', 'y2', 'building_id', 'height']]
    walls['wall'] = walls.apply(lambda r: [[r['x1'], r['y1']],
                                           [r['x2'], r['y2']]], axis=1)
    walls_shape = np.array(list(walls['wall']))

    # Create point light
    pointLightPosition = {'position': [pointlon, pointlat, pointheight]}
    # calculate shadow for walls
    shadowShape = calPointLightShadow_vector(
        walls_shape, walls['height'].values, pointLightPosition)

    walls['geometry'] = list(shadowShape)  #阴影存储
    walls['geometry'] = walls['geometry'].apply(lambda r: Polygon(r))  #将numpy转换成polygon形式
    walls = gpd.GeoDataFrame(walls)  #8列 x1 ， y1 ，x2 ， y2 ， building_id ， height ，  wall ，shadow 
    wallsBuilding = pd.concat([walls, building]) #
    #print(wallsBuilding)
    if merge:
        wallsBuilding = wallsBuilding.groupby(['building_id'])['geometry'].apply(
            lambda df: MultiPolygon(list(df)).buffer(0)).reset_index()
        #print(wallsBuilding)
    shadows=wallsBuilding
    return shadows