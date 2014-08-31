#!python3

## Tommy Carstensen, 2014

import numpy as np
import math
import fractions
import os
import collections
import csv
import json
import shapely
from shapely.geometry import Point


def main():

## http://sedac.ciesin.columbia.edu/data/set/gpw-v3-population-count/data-download
## http://en.wikipedia.org/wiki/Esri_grid

##    lcm = nrows2*n/fractions.gcd(nrows2,n)
##    print(lcm,fractions.gcd(nrows2,n))
##    print(90*8/fractions.gcd(90,8),fractions.gcd(90,8))
##    print(8/fractions.gcd(90,8),90/fractions.gcd(90,8))
##    print(nrows2/fractions.gcd(nrows2,n))
##    stop

    ## argparse
    affix = 'afp00ag'
    affix = 'afp00g'
    affix = 'afds00ag'
    affix = 'afds00g'
    n = nrows = ncols = 5*48*1  # script fast if multiple of 2160
    zero = 0  # add zero values to average?

    layers = 27
    norm = 'log10'
    norm = 'unity'  # feature scaling
    norm = 'log2'
    norm = 'log'

    affix1 = 'out_NASA2LEGO/%s_%ix%i_zero%i' %(affix, n, n, zero)

    if not os.path.isfile('%s.npy' %(affix1)):
##        array_gpw = read_gpw('afp00g.asc')
        array_gpw = read_gpw('%s.asc' %(affix))

        array_LEGO = np.zeros((n, n))
        array_LEGO_cnt = np.zeros((n, n))

        nrows2, ncols2 = np.shape(array_gpw)

        assert np.shape(array_gpw)[0] > n

        den = n/fractions.gcd(ncols2, n)
        num = int(ncols2/fractions.gcd(ncols2, n))

        print('converting NASA array to LEGO array')
        for x1 in range(n):
            print('gpw2LEGO', x1, n)
            for y1 in range(n):
    ##            print(x1,y1)
                for x2 in range(num):
                    for y2 in range(num):
    ##                    x3 = int((nrows2*x1+x2)/n)
    ##                    y3 = int((ncols2*y1+y2)/n)
                        x3 = int((x1*num+x2)/den)
                        y3 = int((y1*num+y2)/den)
                        array_LEGO[x1][y1] += array_gpw[x3][y3]
                        if array_gpw[x3][y3] > 0:
                            array_LEGO_cnt[x1][y1] += 1

        if zero == 0:
            for row in range(n):
                for col in range(n):
                    if array_LEGO[row][col] == 0:
                        continue
                    array_LEGO[row][col] /= array_LEGO_cnt[row][col]

        np.save(affix1, array_LEGO)

        del array_gpw
        del array_LEGO_cnt

    else:
        array_LEGO = np.load('%s.npy' %(affix1))

##    x = []
##    for row in range(240):
##        for col in range(240):
##            if array_LEGO[row][col] == 0:
##                continue
####            if array_LEGO[row][col] >= 200:
####                continue
##            x += [math.log(array_LEGO[row][col])]
####            x += [array_LEGO[row][col]]
##    print(len(x))
##    import matplotlib.pyplot as plt
##    plt.xlabel('Population Density (arbitrary unit)')
##    plt.ylabel('Frequency')
##    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
####    hist, bins = np.histogram(array_LEGO, bins=50)
####    width = 0.7 * (bins[1] - bins[0])
####    center = (bins[:-1] + bins[1:]) / 2
#####    plt.bar(center, hist, width=width)
####    plt.bar(center, hist)
##    plt.show()
##    stop

    array_LEGO = normalize(array_LEGO, layers, norm)

##    x = []
##    for row in range(240):
##        for col in range(240):
##            if array_LEGO[row][col] == 0:
##                continue
####            if array_LEGO[row][col] >= 200:
####                continue
####            x += [math.log(array_LEGO[row][col])]
##            x += [array_LEGO[row][col]]
##    print(len(x))
##    import matplotlib.pyplot as plt
##    plt.xlabel('Population Density (arbitrary unit)')
##    plt.ylabel('Frequency')
##    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
####    hist, bins = np.histogram(array_LEGO, bins=50)
####    width = 0.7 * (bins[1] - bins[0])
####    center = (bins[:-1] + bins[1:]) / 2
#####    plt.bar(center, hist, width=width)
####    plt.bar(center, hist)
##    plt.show()
##    stop

    array_buried = find_buried(array_LEGO, layers)

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header('%s.asc' %(affix))

    array_colors = json2array(
        array_LEGO, nrows, ncols, xllcorner, yllcorner, cellsize)

    lxfml = '%s_y%i_%s.lxfml' %(affix1, layers, norm)
    numpy2lxfml(
        array_LEGO, lxfml, array_colors, array_buried,
        nrows, ncols, xllcorner, yllcorner, cellsize)

##    print(pcount_max)
##    print(np.amin(array_gpw))
##    from collections import Counter
##    print(Counter([float(x) for x in np.nditer(array_gpw)]))

    return


def find_buried(array_LEGO, layers):

    array_buried = np.empty_like(array_LEGO, int)
    n = np.shape(array_LEGO)[0]

    for row in range(n):
        for col in range(n):
            ## no buried plates if only 1 plate
            if array_LEGO[row][col] <= 1:
                continue
            ## nothing buried at the edge (no bricks at the edge anyway)
            if row == 0 or row == n-1 or col == 0 or col == n-1:
                continue
            ## minimum neighbouring height
            z = min([array_LEGO[row+x][col+y]
                     for x in range(-1,2) for y in range(-1,2)])
            ## not buried by neighbouring plates
            if z == 0:
                continue
            array_buried[row][col] = z-1

    return array_buried


def read_gpw_header(file_gpw):

    with open(file_gpw) as f:
        ncols = int(f.readline().strip().split()[1])
        nrows = int(f.readline().strip().split()[1])
        assert ncols%2 == 0 and nrows%2 == 0
        assert ncols > nrows
        assert ncols == 2160  # -30+2160*2.5/60 = 60 = 60E
        assert nrows == 1920  # -40+1920*2.5/60 = 40 = 40N
        xllcorner = int(f.readline().strip().split()[1])
        yllcorner = int(f.readline().strip().split()[1])
        assert xllcorner == -30  # -30 = 30E
        assert yllcorner == -40  # -40 = 40S
        cellsize = float(f.readline().strip().split()[1])
        assert round(cellsize,8) == round(2.5/60,8)

    return ncols, nrows, xllcorner, yllcorner, cellsize


def normalize(a, layers, norm):

    n = np.shape(a)[0]

    assert layers % 3 == 0
    ## https://en.wikipedia.org/wiki/Normalization_(statistics)
    if norm == 'log10':
        amax = np.amax(a)
        a = np.log10(np.divide(a, amax/(10**layers)))
    elif norm == 'log':
        amax = np.amax(a)
        amin = np.amin(a[np.nonzero(a)])
##        cnt = 0
##        sumx = 0
##        sumxx = 0
        for row in range(n):
            for col in range(n):
                if a[row][col] == 0:
                    continue
##                a[row][col] = max(1,math.log(
##                    math.exp(layers)*a[row][col]/amax))
                log = math.log(a[row][col])
##                cnt += 1
##                sumx += log
##                sumxx += log*log
                a[row][col] = max(1,layers*log/math.log(amax))
##                a[row][col] = layers*(log-math.log(amin))/(math.log(amax)-math.log(amin))                                
####        mean = sumx/cnt
####        var = sumxx/cnt-mean**2
##        print(np.percentile(a,100))
##        print(np.percentile(a,99.99))
##        print(np.percentile(a,50))
##        print(np.percentile(a,0.1))
##        print(np.percentile(a,0.01))
##        stop
####        a = np.log(np.divide(a, amax/(math.exp(layers))))
    elif norm == 'unity':
        amin = 0
        amax = np.amax(a)
        a = np.divide(a, amax/layers)
    elif norm == 'log2':
        amax = np.amax(a)
        a = np.log2(
            np.divide(a, amax/(2**layers)))
    else:
        stop

    assert math.ceil(np.amax(a)) == layers

    for row in range(n):
        for col in range(n):
            if a[row][col] < 0:
                a[row][col] = 1
            else:
                a[row][col] = math.ceil(a[row][col])

    return a


def json2array(array_LEGO, nrows, ncols, xllcorner, yllcorner, cellsize):

    ## colors not for sale via Pick-a-Brick:
    ## 2 light gray / grey
    ## 27 dark gray / dark grey

    ## only available via service.lego.com/en-gb/replacementparts/#WhatIndividualBrickBuy/3024
    ## 138 dark tan / sand yellow
    ## 222 bright pink / light purple
    ## 119 lime / bright yellowish green
    ## 102 medium blue / medium blue
    ## 141 dark green / earth green
    ## 330 olive green / olive green
    ## 154 dark red / new dark red

    ## currently out of stock
    ## 268 dark purple / medium lilac

    ## possibly discontinued
    ## 151 sand green / sand green
    ## 321 dark azure / dark azur
    ## 25 brown / earth orange
    ## 135 sand blue

    ## available via PAB but not used
    ## 140 earth blue

    ## Official LDD/PAB color IDs
    d_colors = {
        ## available via PAB
        'Bantu':1, ## White (51380)
        'Bantu / Bantu':1,
        'Semitic: Arab, Bedouin':194, ## Light Bluish Gray (47289)
        'Cushitic':192, ## Reddish Brown (39398)
        'Chadic':5, ## Tan (74477) / Brick Yellow
        ## not available via PAB
        'Saharan':2, ## Light Gray (14325) / Grey
        'Saharan / Nilotic':2,
        'Saharan / Cushitic':2,
        'Chadic / Cushitic':2,
        ## available via PAB
        'Nilotic':106, ## Orange (35644)
        'Nilotic / Bantu':106,
        'Nilotic / Bantoid':106,
        'Chari-Nile / Nilotic':106,
        ## not available via PAB
        'Berber':138, ## Dark Tan (32656) / Sand Yellow
        'Northern Mande':222, ## Bright Pink (34034) / Light Purple
        ## available via PAB
        'Voltaic':28, ## Green (39516) / Dark Green
        'Kru':26, ## Black (35741)
        'Adamawa-Ubangian':23, ## Blue (70095) / Bright Blue
        'Bantoid':21, ## Red (67551) / Bright Red
        'Chari-Nile':199, ## Dark Bluish Gray (43621) / Dark Stone Grey (4210719)
        'Adamawa-Ubangian / Chari-Nile':199,
        'West Atlantic':24, ## Yellow (67832)
        ## not available via PAB
        'Fufulde':119, ## Lime (26976) / Bright Yellowish Green
        'Chadic / Fufulde':119,
        'Fufulde / Adamawa-Ubangia':119,
        'San':151, ## Sand Green (6989) / Sand Green
        'Songhai':102, ## Medium Blue (17453) / Medium Blue
        'Gbaya':141, ## Dark Green (10885) / Earth Green
        'Kwa':321, ## Dark Azure (23275) / Dark Azur (Ivory Coast)
        'Maban':268, ## Dark Purple (21162) / Medium Lilac
        'Southern Mande':330, ## Olive Green (12496) / Olive Green (Ivory Coast)
        'Kordofanian':27, ## Dark Gray (7760)
        'Sandawe':154, ## Dark Red (7042) / New Dark Red
        'Khoi: Nama, Bergdama':25, ## Brown (5452) / Earth Orange
        'Fur':135, ## Sand Blue (4126)
        }

    array_colors = np.empty_like(array_LEGO, int)
    n = np.shape(array_colors)[0]

    with open('etnicity_felix.json') as f:
        js = json.load(f)

    for feature in js['features']:
        polygon = shapely.geometry.shape(feature['geometry'])
        family = feature['properties']['FAMILY']
        if family == 'Malagasy':
            continue
        try:
            color = d_colors[family]
        except KeyError:
            color = 324  # medium lavender (e.g. Afrikaans)
            if family != 'Malagasy':
                print(family, feature['properties']['ETHNICITY'])
        if family in ('',):
##            print(feature['properties'], polygon.bounds, feature['id'], feature['properties']['ID'])
##            continue
            color = 322  # medium azure
##            for key in feature.keys():
##                if key == 'geometry': continue
##                print(key, feature[key])
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        row_min = n*(min_lat-yllcorner+cellsize*(ncols-nrows)/2)/(cellsize*max(nrows,ncols))-0.5
        row_max = n*(max_lat-yllcorner+cellsize*(ncols-nrows)/2)/(cellsize*max(nrows,ncols))-0.5
        col_min = n*(min_lon-xllcorner)/(cellsize*max(nrows,ncols))-0.5
        col_max = n*(max_lon-xllcorner)/(cellsize*max(nrows,ncols))-0.5
        within = False
        ## loop from South to North
        for row in range(math.floor(row_min), math.ceil(row_max)):
            latitude = cellsize*(row+0.5)*max(nrows,ncols)/n+yllcorner-cellsize*(ncols-nrows)/2
            ## Bantu languages in Angola incorrectly assigned to the Kru family
            if family == 'Kru' and latitude < 0:
                family = 'Bantu'
                color = d_colors[family]
            ## loop from West to East
            for col in range(math.floor(col_min), math.ceil(col_max)):
                longitude = cellsize*(col+0.5)*max(nrows,ncols)/n+xllcorner
                ## incorrectly assigned to the Maban family
                if family == 'Maban' and longitude > 30:
                    within = True
                    continue
                point = shapely.geometry.Point(longitude, latitude)
                ## solve the point-in-polygen problem
##                if max_lat < -20 and family != 'Miscellaneous / Unclassified':
##                    print(family)
##                    print(polygon.bounds)
##                    print(row_min, row_max, col_min, col_max)
##                    print('latlon', latitude, longitude)
##                    print(feature['properties']['ETHNICITY'])
##                    print("venda longitude 29'E")
##                    stop1
##                if max_lon < -15 and family != 'Miscellaneous / Unclassified':
##                    print(family)
##                    print(polygon.bounds)
##                    print(
##                        math.floor(row_min), math.ceil(row_max),
##                        math.floor(col_min), math.ceil(col_max))
##                    print('latlon', latitude, longitude)
##                    print(feature['properties']['ETHNICITY'])
##                    print("wolof longitude 16'W")
##                    print('rowcol', row, col)
##                    stop2
##                if feature['properties']['ID'] == 1820:
##                    print(family)
##                    print(polygon.bounds)
##                    print(
##                        math.floor(row_min), math.ceil(row_max),
##                        math.floor(col_min), math.ceil(col_max))
##                    print('latlon', latitude, longitude)
##                    print(feature['properties']['ETHNICITY'])
##                    print('rowcol', row, col)
##                    stop4
                if polygon.contains(point):
##                    if feature['properties']['ID'] == 1820:
##                        print('aaa',array_LEGO[240-row][col])
##                        print('bbb',array_LEGO[row][col])
##                        print(polygon)
##                    if max_lat < -20 or max_lon < -15:
##                        print(family)
##                        print(polygon.bounds)
##                        print(row_min, row_max, col_min, col_max)
##                        print(latitude, longitude)
##                        stop3
                    within = True
                    array_colors[n-row][col] = color
##                    if not array_LEGO[240-row][col]:
##                        array_LEGO[240-row][col] = 1
##                    print(array_LEGO[row][col])
##                    print(family)
##                    print(polygon.bounds)
##                    print(
##                        math.floor(row_min), math.ceil(row_max),
##                        math.floor(col_min), math.ceil(col_max))
##                    print('latlon', latitude, longitude)
##                    print(feature['properties']['ETHNICITY'])
##                    print('rowcol', row, col)
##                    stop
                    pass
                ## polygon might not contain point rounded to nearest grid value
                ## hence add manually
                elif row_max-row_min < 2.5 or col_max-col_min < 2.5:
                    within = True
                    ## don't change color if not within polygon
                    if not array_colors[n-row][col] == 0:
                        continue
                    array_colors[n-row][col] = color
                    pass
                ## continue loop over cols
                continue
            ## continue loop over rows
            continue
##        if family == 'Maban' and longitude > 30:
##            print(row, col, latitude, longitude, family, feature['properties']['ETHNICITY'])
        if within == False:
            print(
                feature['properties']['ETHNICITY'], family, polygon.bounds,
                row_max-row_min, col_max-col_min,
                )
        ## continue loop over features
        continue

    return array_colors


def numpy2lxfml(
    array_LEGO, lxfml, array_colors, array_buried,
    nrows, ncols, xllcorner, yllcorner, cellsize):

    n = np.shape(array_LEGO)[0]
    a = 6378137  # semi major axis
    f = 1/298.257223563  # reciprocal flattening
    e_squared=6.69437999014*10**-3  # first eccentricity squared
    h = 0

    ## initiate refID count
    refID = 0
    ## open file
    with open(lxfml, 'w') as f:
        f.write(head('Africa'))
        ## loop from North to South
        for row in range(n):
            latitude = yllcorner-cellsize*(ncols-nrows)/2+cellsize*(row+0.5)*max(nrows,ncols)/n
            if row % 10 == 0:
                print('numpy2lxfml, row', row, 'of', n-1)
            ## loop from West to East
            for col in range(n):
                longitude = xllcorner+cellsize*(col+0.5)*max(nrows,ncols)/n
##                print(row,col,latitude,longitude)
##                N = a/math.sqrt(1-e_squared*math.sin(latitude)**2)
####                x = (N+h)*math.cos(latitude)*math.cos(longitude)
####                y = (N+h)*math.cos(latitude)*math.sin(longitude)
##                z = (N*(1-e_squared)+h)*math.sin(latitude)
##                z /= 1000  # m to km
##                cm_per_degree = 240*0.8/90
##                km_per_degree = 40076/360
##                bricks_per_cm = 1/0.96
##                z /= km_per_degree  # km to degree
##                z *= cm_per_degree  # degree to cm
##                z *= bricks_per_cm  # cm to bricks
##                z = abs(int(z))
##                if col == 0:
##                    print(row,col,latitude,z)

                materialID = array_colors[row][col]

                ## skip parts of array not covered by GeoJSON array
                if materialID == 0: continue  # tmp!!!

                ## make sure bricks are present
                ## for areas covered by language family polygons
                if materialID and array_LEGO[row][col] == 0:
                    array_LEGO[row][col] = 1

##                print(row,col,latitude,longitude)

                for y in range(int(array_LEGO[row][col])):
                    designID = 3024  # 1x1 plate
                    if array_LEGO[row][col] >= 3*(y//3)+3:
                        if y%3 == 0:
                            designID = 3005  # 1x1 brick
                        else:
                            continue
##                    if array_buried[row][col] >= 3*(y//3)+3:
##                        if y%3 == 0:
##                            designID = 3005  # 1x1 brick
##                        else:
##                            continue

##                    if array_buried[row][col] >= y:
##                        materialID = 1
##                    else:
##                        materialID = array_colors[row][col]

##                    materialID = 119
####                    if row < 192/2 and col < 192/2:
####                        materialID = 1 ## white/NW
####                    elif row < 192/2 and col > 192/2:
####                        materialID = 221 ## violet/NE
####                    elif row > 192/2 and col < 192/2:
####                        materialID = 119 ## green/SW
####                    elif row > 192/2 and col > 192/2:
####                        materialID = 26 ## black/SE
##
####                    ## skip if not Middle East (NE)
####                    if row > 192/2 or col < 192/2:
####                        continue
##
####                    if row > 50 and row < 90:
####                        materialID = row-50
##
##
##                    if int(array_LEGO[row][col]) == 1:
##                        materialID = 1
##                    elif int(array_LEGO[row][col]) == 2:
##                        materialID = 226
##                    elif int(array_LEGO[row][col]) == 3:
##                        materialID = 119
##                    elif int(array_LEGO[row][col]) == 4:
##                        materialID = 28
##                    elif int(array_LEGO[row][col]) == 5:
##                        materialID = 212
##                    elif int(array_LEGO[row][col]) == 6:
##                        materialID = 23
##                    elif int(array_LEGO[row][col]) == 7:
##                        materialID = 222
##                    elif int(array_LEGO[row][col]) == 8:
##                        materialID = 221
##                    elif int(array_LEGO[row][col]) == 9:
##                        materialID = 21
##
##                    if int(array_LEGO[row][col]) <= 3:
##                        materialID = 1
##                    elif int(array_LEGO[row][col]) <= 6:
##                        materialID = 226
##                    elif int(array_LEGO[row][col]) <= 9:
##                        materialID = 119
##                    elif int(array_LEGO[row][col]) <= 12:
##                        materialID = 28
##                    elif int(array_LEGO[row][col]) <= 15:
##                        materialID = 212
##                    elif int(array_LEGO[row][col]) <= 18:
##                        materialID = 23
##                    elif int(array_LEGO[row][col]) <= 21:
##                        materialID = 222
##                    elif int(array_LEGO[row][col]) <= 24:
##                        materialID = 221
##                    elif int(array_LEGO[row][col]) <= 27:
##                        materialID = 21

                    line = '        <Part'
                    line += ' refID="{0}"'.format(refID)
                    line += ' designID="{0}"'.format(designID)
                    line += ' materialID="{0}"'.format(materialID)
                    line += ' itemNos="{0}{1}"'.format(designID, materialID)
                    line += ' angle="0" ax="0" ay="1" az="0"'
                    line += ' tx="{0:.1f}"'.format(row*-0.8)
                    line += ' ty="{0:.2f}"'.format(y*0.32)
##                    line += ' ty="{0}"'.format((y-z)*0.32)
                    line += ' tz="{0:.1f}"'.format(col*0.8)
                    line += '/>\n'
                    f.write(line)
                    refID += 1
                    ## continue loop over ty
                    continue
                ## continue loop over tz
                continue
            ## continue loop over tx
            continue
        f.write(tail())
        ## close file
        pass

    return


def head(name):

    s = '''\
<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<LXFML versionMajor="4" versionMinor="0" name="{0}">
  <Meta>
    <Application name="LEGO Digital Designer" versionMajor="4" versionMinor="3"/>
    <Brand name="LDD"/>
    <BrickSet version="1264"/>
  </Meta>
  <Cameras>
   <Camera refID="0" fieldOfView="80" distance="467.239" angle="73.37" ax="-0.286" ay="-0.937" az="-0.2" tx="-400" ty="191" tz="147"/>
  </Cameras>
  <Scene cameraRefID="0">
    <Model>
      <Group refID="0" angle="0" ax="0" ay="1" az="0" tx="0" ty="0" tz="0">
'''.format(name)

#   <Camera refID="0" fieldOfView="80" distance="400" angle="60" ax="-0.6" ay="-0.6" az="-0.25" tx="0.0" ty="0.0" tz="0.0"/>

    return s


def tail():

#    ty += 0.32 ## per plate in height
    s = '''\
     </Group>
    </Model>
  </Scene>
  <BuildingInstructions>
  </BuildingInstructions>
</LXFML>
'''

    return s


def read_gpw(file_gpw):

    with open(file_gpw) as f:
        ncols = int(f.readline().strip().split()[1])
        nrows = int(f.readline().strip().split()[1])
        assert ncols%2 == 0 and nrows%2 == 0
        assert ncols > nrows
        assert ncols == 2160  # -30+2160*2.5/60 = 60 = 60E
        assert nrows == 1920  # -40+1920*2.5/60 = 40 = 40N
        xllcorner = int(f.readline().strip().split()[1])
        yllcorner = int(f.readline().strip().split()[1])
        assert xllcorner == -30  # -30 = 30E
        assert yllcorner == -40  # -40 = 40S
        cellsize = float(f.readline().strip().split()[1])
        assert round(cellsize,8) == round(2.5/60,8)
        for i in range(1):
            f.readline()

        array_gpw = np.zeros((max(nrows, ncols), max(nrows, ncols)))
        max_gpw = 0
        x1 = (max(nrows, ncols)-nrows)/2
        x2 = (max(nrows, ncols)+nrows)/2
        for row in range(max(nrows, ncols)):
##            print(row)
            latitude = yllcorner+(nrows-row+x1-.5)*2.5/60
            if row >= x1 and row < x2:
                cols = f.readline().rstrip().split()
                pass
            else:
                cols = ['0']*ncols
            N,W = '13.146,43.271'.split(',')
            N = float(N)
            W = float(W)
            Ndelta=0.02
            Wdelta=0.02
            Ndelta=0.2
            Wdelta=0.2
            if latitude < N+Ndelta and latitude > N-Ndelta:
                print()
##                print(latitude,longitude,cols)
            for col,s in enumerate(cols):
                longitude = xllcorner+(col+.5)*2.5/60
##                if False: pass
##                ## https://en.wikipedia.org/wiki/Geography_of_Africa#Extreme_points
##                ## North, Iles des Chiens, Tunisia
##                elif latitude > 37+32/60:
##                    s = 0
##                ## South, Cape Agulhas, South Africa
##                elif latitude < -34-51/60-15/3600:
##                    s = 0
##                ## North, Cape Blanc / Ras ben Sakka, Tunisia
##                elif latitude > 37+21/60:
##                    s = 0
##                ## West, Pointe des Almadies, Senegal
##                elif longitude < -17-53/60:
##                    s = 0
##                ## East, Ras Hafun / Cape Guardafui, Somalia
##                elif longitude > 51+27/60+52/3600:
##                    s = 0
##                ## Mediterranean (Pantelleria, Lampedisa, Malta, Crete,
##                ## Gavdos, Akrotiri, Italy, Greece, Turkey, Cyprus)
##                ## https://en.wikipedia.org/wiki/Extreme_points_of_Europe#Extremes_of_the_European_continent.2C_including_islands
##                ## https://en.wikipedia.org/wiki/Extreme_points_of_Greece   
##                elif col >= 998 and latitude > 30+8/60+43/3600:
##                    s = 0
##                ## Mediterranean
##                elif col >= 675 and col <= 985 and row <= 204+(col-675)*(169-204)/(985-675):
##                    s = 0
##                ## Strait of Gibraltar (incl. Punta de Tarifa)
##                elif col >= 584 and col <= 675 and row <= 217+(col-584)*(204-217)/(675-584):
##                    s = 0
##                elif col == 591 and row <= 213: s = x  # Punta de Europa
##                ## Israel/Egypt Border
##                elif row >= 328 and row <= 372 and col > 1541+(row-328)/((372-328)/(1557-1541)):  # do more accurate slope with more accurate coordinates and then convert to grid
##                    s = 0
##                ## Middle East (from South to North)
##                ## Gulf of Aqaba (Tiran Island ignored)
##                elif row >= 372 and row <= 407 and col > 1558+(row-372)/((406-372)/(1548-1558)):
##                    s = 0
##                elif row == 372 and col >= 1559: s = 0  # Eifat/Aqaba/Jordan/Israel
##                elif row == 410 and col >= 1566: s = 0
##                ## Red Sea (Hanish Islands ignored)
##                elif row <= 730 and row >= 411 and col > 1706+(row-689)/((428-689)/(1555-1706)):
##                    s = 0
##                ## Yemen
##                elif row == 731 and col >= 1750: s = 0
##                elif row == 740 and col >= 1754: s = 0
##                elif row == 745 and col >= 1755: s = 0
##                elif row == 748 and col >= 1757: s = 0
##                elif row == 758 and col >= 1758: s = 0
##                elif row == 764 and col >= 1758: s = 0  # 13.146,43.271
##                elif row == 765 and col >= 1759: s = 0  # 13.104,43.312
##                elif row == 766 and col >= 1760: s = 0  # 13.062,43.354
##                elif row == 767 and col >= 1760: s = 0  # 13.021,43.354
##                elif row == 768 and col >= 1761: s = 0  # 12.979,43.396
##                elif row == 769 and col >= 1764: s = 0  # 12.938,43.521
##                elif row == 770 and col >= 1762: s = 0  # 12.896,43.438
##                elif row == 771 and col >= 1762: s = 0  # 12.854,43.438
##                elif row == 772 and col >= 1763: s = 0  # 12.812,43.479
##                elif row == 773 and col >= 1763: s = 0  # 12.771,43.479
##                elif row == 774 and col >= 1763: s = 0  # 12.729,43.479
##                elif row == 775 and col >= 1762: s = 0  # 12.688,43.438
##                elif row == 776 and col >= 1761: s = 0  # 12.646,43.396 / Perim Island

                if s == '0':
                    array_gpw[row][col] = 0
                elif s == '-9999':
                    array_gpw[row][col] = 0
                else:
                    array_gpw[row][col] = float(s)
                if latitude < N+Ndelta and latitude > N-Ndelta and longitude > W-Wdelta and longitude < W+Wdelta:
                    print(s,row,col,'%.3f,%.3f' %(latitude,longitude))

##        ## first loop to parse values
##        array_gpw = np.resize(
##            np.fromstring(
##                f.read().replace('\n',' '), sep= ' '), (nrows, ncols))

    return array_gpw


if __name__ == '__main__':
    main()
