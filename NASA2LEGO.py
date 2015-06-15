#!python3

## Tommy Carstensen, Aug-Sep2014, Apr-May2015
## Acknowledgments:
## Martin Haspelmath
## William Reno
## Friedrich Riha, Manj Sandhu, Joshua Randall, Mutua Matheka, Jane Walsh
## Wikimedia Foundation, Inc., Glottolog, Ethnologue
## Python Software Foundation
## Free Software Foundation
## GNU Project
## NASA Shuttle Radar Topography Mission (SRTM)
## NASA Socioeconomic Data and Applications Center (SEDAC)
## Felix 2001
## The LEGO company.

## todo: add transparent plats for:
## http://en.wikipedia.org/wiki/List_of_lakes_by_area
## 3) Lake Victoria
## 6) Lake Tanganyika
## 9) Lake Malawi
## 22) Turkana
## 28) Albert
## 29) Mweru
## Lake Chad
## Lake Tana / Blue Nile origin

## http://en.wikipedia.org/wiki/List_of_rivers_of_Africa
## http://en.wikipedia.org/wiki/Geography_of_Africa#Rivers
## !1) Nile River
## 2!) Congo River
## 3) Blue Nile
## 3!) Niger
## 3) Uele
## 3!) Zambezi
## 3) Kasai
## 4!) Senegal???
## 4!) Vaal/Orange
## 4) Benue!
## 5!) Volta River and tributaries: black, white, red
## 6) White Volta
## 6) Sassandra
## 6) Sanaga
## 6) Chari
## 6) Atbara
## 6!) Jubba
## 6) Ruvuma
## 6!) Okavango
## 6!) Limpopo
## http://en.wikipedia.org/wiki/Rufiji_River

d_mountains = {}
d_rivers = {
    'Nile': (30, 10, 31, 6),
    'Congo': (-6, -4, 12, 27),
    'Niger': (13, 60*0.86, -3, -60*0.33),
    'Zambezi': (-18, -34, 36, 28),
    'Senegal': (15, 47, -16, -31),
    'Vaal': (-29, -4, 23, 38),
    'Volta': (5, 46, 0, -41),
    'Jubba': (0, 0.2495*60, 42, 60*0.6307),
    'Okavango': (-18, -57, 22, 29),
    'Limpopo': (-25, -10, 33, 35),
    'WhiteNile': (-2, -17, 29, 20),
    'BlueNile': (12, 0, 37, 15),
    'BlueNile': (15, 38, 32, 32),
    'Rufiji': (-8, 0, 39, 20),
    }
d_lakes = {
    'Lake Victoria': (-1, 0, 33, 0),
    'Lake Tanganyika': (-6, -30, 29, 50),
    'Lake Malawi': (-12, -11, 34, 22),
    'Lake Turkana': (3, 35, 36, 7),
    'Lake Albert': (1, 41, 30, 55),
    'Lake Mweru': (-9, 0, 28, 43),
    'Lake Chad': (13, 0, 14, 0),
    'Lake Tana': (12, 0, 37, 15),
    }

## todo: Slope 30 1 x 1 x 2/3 for mountains or anything above an altitude threshold
## try 3000m or 2000m
## http://en.wikipedia.org/wiki/List_of_highest_mountain_peaks_of_Africa
## http://en.wikipedia.org/wiki/Geography_of_Africa#Plateau_region
## http://en.wikipedia.org/wiki/Geography_of_Africa#Mountains
## https://pypi.python.org/pypi/SRTM.py
## https://github.com/tkrajina/srtm.py

## Flags for: Yoruba, Igbo, Fula, Shona, Zulu

import numpy as np
import math
import fractions
import os
import collections
import csv
import json
import shapely
from shapely.geometry import Point
import itertools
import argparse
import matplotlib.pyplot as plt


## 
d_len2dIDs = {
    ## plates
    1: {
        ## width 1
        1: {1: 3024, 2: 3023, 3: 3623, 4: 3710, 6: 3666, 8: 3460},
        2: {1: 3023, 2: 3022, 3: 3021, 4: 3020, 6: 3795, 8: 3034},
        },
    ## bricks
    3: {
        ## width 1
        1: {1: 3005, 2: 3004, 3: 3622, 4: 3010, 6: 3009, 8: 3008},
        ## width 2
        2: {1: 3004, 2: 3003, 3: 3002, 4: 3001, 6: 2456, 8: 3007},
##        2: {1: 3004, 2: 3003, 3: 3002, 4: 3001, 6: 44237, 8: 93888},
        },
    }
d_dIDs2dim = {}
for h in d_len2dIDs:
    for w in sorted(d_len2dIDs[h].keys()):
        for l, ID in d_len2dIDs[h][w].items():
            if ID in d_dIDs2dim.keys():
                continue
            d_dIDs2dim[ID] = {'dim0': w, 'dim1': l}
## square plates
d_dIDs_sq_plate = {
##    16: 91405,  # 16x16 (GBP2.58) only marginally cheaper than 8x8 (GBP0.67)
    8: 41539,
    6: 3958,  # 6x6 (GBP0.43) more expensive than 4x4 (GBP0.18) per area
    4: 3031,
    }
designID_baseplate = 4186
materialID_grey = 194
for k, v in d_dIDs_sq_plate.items():
    d_dIDs2dim[v] = {'dim0': k, 'dim1': k}
d_dIDs2dim[3001] = {'dim0': 2, 'dim1': 4}
d_dIDs2dim[d_len2dIDs[3][2][6]] = {'dim0': 2, 'dim1': 6}
d_dIDs2dim[d_len2dIDs[3][2][8]] = {'dim0': 2, 'dim1': 8}
## Ubiquitous material IDs; i.e. 1x1 and 2x4 bricks available in PAB.
ubiquitous = {
    1,  # White
    21,  # Red
    23,  # Blue
    24,  # Yellow
    26,  # Black
    28,  # Green
    5,  # Tan
    106,  # Orange
    194,  # Stone Grey
    199,  # Stone Grey
    192,  # Reddish Brown
    119,  # Lime
##    222,  # Light Purple / Light Pink
    }

materialID_buried = 999

## Tuples are sorted from common to rare.
d_colorfamily2color = {
    'red': (
        21,  # red
        154,  # new dark red
##        192,  # reddish brown
        ),  # Khoisan
    'black':(26,),  #  IndoEuropean/Afrikaans
    ## Niger-Congo
    'blue': (
        23, 102, 140,
        1,  # White looks good with bright blue colors and is cheaper.
        323,  # Aqua / Unikitty Blue (1x1 plate out of stock, but got from BL)
##        322, # Medium Azure 322 (1x1 brick - 1x1 plate out of stock)
####        321, # # Dark Azure 321 (1x2 plate - no 1x1 size)
        ),
    ##  AfroAsiatic
    'yellow-orange-brown': (
        24,  # yellow
        5,  # tan
        106,  # orange
        192,  # reddish brown
        138,  # dark tan
        191,  # flame orange
        226,  # cool yellow / bright light yellow
        ),
    ## Bantu
    'green': (
        28,  # green
        119,  # lime
        141,  # earth green / dark green
        330,  # olive
        326,  # Spring Yellowish Green / Yellowish Green (only 1x1 plate available)
##        151,  # / Sand Green (got 500 1x1 plates from BL, but otherwise NA)
        ),
    'purple': (
        222,  # Light Purple / Bright Pink
        124,  # Bright Reddish Violet / Magenta
        324,  # Medium Lavender
##        268,  # Medium Lilac
##        221,  # Bright Purple / Dark Pink (discontinued in 2004)
        ## "The medium lilac element ID 4224857 is out of stock
        ## and we have no plans on bringing it back at this time."
        ),  # NiloSaharan
    'grey': (
        194,  # medium stone grey
        199,  # dark stone grey
        ),  # Other...
##    'white': {1},  # buried
    }
d_color2family = {
    color: family for family, colors in d_colorfamily2color.items()
    for color in colors}

d_max_len = {
    ## plates
    1: {
        ## width 1
        1: {
            1: 12,
            26: 10, 5: 10, 192: 10, 194: 10,
            21: 8, 23: 8, 24: 8, 28: 8, 199: 8,
            106: 6, 102: 6, 141: 6, 140: 6, 154: 6,
            119: 4, 324: 4, 222: 4,
            330: 2,
            },
        2: {
            1: 12, 26: 12, 199: 12, 194: 12,
            21: 10, 5: 10,
            23: 8, 24: 8, 28: 8, 192: 8,
            106: 6, 119: 6, 222: 6,
            154: 4, 141: 4,
            330: 1,
            },
        },
    ## bricks
    3: {
        ## width 1
        1: {
            1: 8, 2: 3004, 3: 3622, 4: 3010, 6: 3009,
            191: 4,
            },
        ## width 2
        2: {
            1: 8,
            191: 1, 303: 1, 324: 1,
            },
        },
    }

## Not available in replacement parts.
d_NA = {
    ## plates
    1: {
        ## width 1
        1: {
            8: {191, 221, 138, 330, 330, 138, 221, 191, 191, 326, 124,},
            6: {330, 141, 221, 268, 323, 326, 151,},
            4: {330, 221, 268, 323, 326, 151,},
            3: {138, 330, 324, 138, 124, 221, 268, 323, 106, 102,},
            2: {323, 326, 322,},  # Aqua, Medium Azure
##            1: {321},  # Dark Azure
            1: {
                322,  # medium azure
##                323,  # aqua / unikitty blue (got 500 1x1 plates from BL)
                268,  # medium lilac (dark purple)
                221,  # bright purple (dark pink)
                321,  # dark azure
                151,  # sand green
                },
            },
        ## width 2
        2: {
            8: {221, 268, 330, 326, 28,},
            6: {221, 268, 326, 141, 124,},
            ## 3020
            4: {221, 268, 323, 326, 222, 151, 222, 124,},
            ## 3021
            3: {191, 124, 330, 221, 268, 323, 326, 154, 119, 102, 151, 326,},
            2: {
                221, 141, 330, 324, 221, 268, 323, 330, 326, 222, 138,
                124, 151, 330, 124, 330, 140, 140, 330, 321,
                },
            1: {},
            },
        },
    ## bricks
    3: {
        ## width 1
        1: {
            8: {
                221, 268, 222, 124, 323, 322, 326, 119, 324, 102, 106,
                140, 141,
                28,  # green available as replacement but 0.42GBP VS two 1x4 = 0.12GBP
                },
            6: {
                191, 221, 140, 324, 221, 268, 140, 323, 324, 222, 151, 102,
                323, 141, 326, 321,
                119,  # lime available as replacement but 0.38GBP
                },
            4: {268, 323, 326, 151,},
            3: {
                221, 268, 138, 124, 323, 322, 326, 141, 324, 140, 154, 191,
                321,
                },
            2: {326,},  # Spring Yellowish Green / Yellowish Green (from 2015)
            1: {
                226, 326,
                321,  # dark azure
                },  # Cool Yellow (226), Spring Yellowish Green (326)
            },
        ## width 2
        2: {
            8: {221, 268, 330, 222, 138, 324, 323, 322, 28, 192, 28,},
            6: {221, 330, 221, 268, 326, 106, 119, 124, 324,},
            4: {140, 124,},
            3: {191, 330, 221, 268, 324, 323, 322, 326,},
            2: {
                191, 330, 324, 323, 326, 141, 330, 151, 141, 330, 124, 321,
                124, 330, 124,
                },
            1: {},
            },
        },
    }
d_NA[1][2][1] = d_NA[1][1][2]
d_NA[3][2][1] = d_NA[3][1][2]

## Most common bricks:
## https://www.bricklink.com/catalogStats.asp?statID=C&itemType=P&inItemType=&catID=5

## Most common plates:



## TODO: FIRST LOOP 3 LAYERS (BRICK) AND THEN 1 LAYER (PLATES)
## I.E. REPLACE CONNECTED PLATES IN 3 LAYERS WITH CHEAPER 12x24 brick...
## e.g. function find_connected_buried
## https://www.bricklink.com/catalogItem.asp?P=30072
## https://service.lego.com/en-gb/replacementparts#WhatIndividualBrickBuy/30072
## 30072 out of producing since 2008 according to TLG


def main():

    args = argparser()

##    lcm = nrows2*n/fractions.gcd(nrows2,n)
##    print(lcm,fractions.gcd(nrows2,n))
##    print(90*8/fractions.gcd(90,8),fractions.gcd(90,8))
##    print(8/fractions.gcd(90,8),90/fractions.gcd(90,8))
##    print(nrows2/fractions.gcd(nrows2,n))
##    stop

    affix1 = 'lxfml/{0}_{1:d}x{1:d}_dens{2}'.format(
        args.affix, args.n, args.density)

    ## 0-args.layers
    if not os.path.isfile('{}.npy'.format(affix1)):
        ## slow
        a_2D_density = asc2np(args, affix1)
    else:
        ## fast
        a_2D_density = np.load('{}.npy'.format(affix1))

    ## Normalise population density.
    a_2D_density = normalize(a_2D_density, args)

    ## slow - convert ethnicity polygons from json file to numpy array
    ## i.e. assign a tentative materialID to each 2D point
    a_2D_mIDs = json2array(args, a_2D_density)

    ## Do a manual fix of zero density in the Eastern Desert in Egypt.
    a_2D_density = fix_zero_density_in_desert(
        a_2D_density, a_2D_mIDs, args)

    ## Do this before finding buried plates.
    a_2D_density, a_2D_mIDs = color_as_nearby(args, a_2D_density, a_2D_mIDs)

    ## fast
    ## Identify how many plates are buried at each grid position.
    a_2D_buried = find_buried(args, a_2D_density, args.layers, a_2D_mIDs)

    ## slow
    ## 0-(args.layers-1)
    a_3D_dIDs = find_connected_buried(
        args, a_2D_buried, a_2D_density, a_2D_mIDs)

    ##
    a_3D_dIDs, a_3D_mIDs, a_3D_angle = find_connected_exposed(
        args, a_2D_density, a_2D_buried,
        a_2D_mIDs, a_3D_dIDs)

##    ## Add rivers and lakes.
##    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
##        '{}.asc'.format(args.affix))
##    for key, coord in itertools.chain(d_rivers.items(), d_lakes.items()):
##        _den = cellsize * max(nrows, ncols)
##        lat = coord[0]+coord[1]/60
##        lon = coord[2]+coord[3]/60
##        row = int(round(args.n*(lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5, 0))
##        col = int(round(args.n*(lon-xllcorner)/_den-0.5, 0))
##        print(key, coord, lat, lon, row, col)
##        assert a_3D_dIDs[a_2D_density[row][col]][row][col] == 0
##        a_3D_dIDs[a_2D_density[row][col]][row][col] = 3024
##        a_3D_mIDs[a_2D_density[row][col]][row][col] = 40

    ## slow
    lxfml = '{}_y{:d}_{}.lxfml'.format(affix1, args.layers, args.norm)
    numpy2lxfml(
        args,
        a_2D_density, lxfml, a_2D_mIDs, a_2D_buried,
        a_3D_dIDs, a_3D_mIDs, a_3D_angle)

##    print(pcount_max)
##    print(np.amin(a_gpw))
##    from collections import Counter
##    print(Counter([float(x) for x in np.nditer(a_gpw)]))

    return


def histogram2(a_2D_density, args):

    path = 'density2_{}.png'.format(args.affix)
    if os.path.isfile(path):
        return

    x = []
    for row in range(args.n):
        for col in range(args.n):
            if a_2D_density[row][col] == 0:
                continue
            x += [a_2D_density[row][col]]

    plt.xlabel('Normalised population density (arbitrary unit)')
    plt.ylabel('Frequency')
    n, bins, patches = plt.hist(
        x, args.layers, normed=1, facecolor='g', alpha=0.75)
    plt.savefig(path)
    plt.close()

    return


def histogram1(a_2D_density, args):

    path = 'density1_{}.png'.format(args.affix)
    if os.path.isfile(path):
        return

    x = []
    for row in range(args.n):
        for col in range(args.n):
            if a_2D_density[row][col] == 0:
                continue
            x += [math.log(a_2D_density[row][col])]
    plt.xlabel('Natural logarithm of population density (arbitrary unit)')
    plt.ylabel('Frequency')
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
    plt.savefig(path)
    plt.close()

    return


def fix_zero_density_in_desert(a_2D_density, a_2D_mIDs, args):

    ## Make sure bricks are present
    ## for areas covered by language family polygons.
    ## Creates low density at coast.
    ## So only do this for the Eastern Desert of Egypt.
    ## https://en.wikipedia.org/wiki/List_of_tripoints
    ## https://en.wikipedia.org/wiki/Geography_of_Egypt#Extreme_points
    ## https://en.wikipedia.org/wiki/22nd_parallel_north
    ## https://en.wikipedia.org/wiki/Bir_Tawil
    ## Egypt, Sudan, Libya tripoint
    tripoint_SW_lat = 22
    tripoint_SW_lon = 25
    ## Egypt, Palestine/Gaza, Israel tripoint
    tripoint_NE_lat = 31.216667
    tripoint_NE_lon = 34.266667
    ## https://en.wikipedia.org/wiki/List_of_countries_by_easternmost_point
    easternmost_lon = 35.75

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))

    ## Loop North to South.
    for row in range(a_2D_density.shape[0]):
        dy = cellsize*(row+0.5)*max(nrows, ncols)/args.n
        dy -= cellsize*(ncols-nrows)/2
        latitude = -yllcorner-dy
        ## Loop West to East.
        for col in range(a_2D_density.shape[1]):
            dx = cellsize*(col+0.5)*max(nrows, ncols)/args.n
            longitude = xllcorner+dx
            ## Language polygon, but no population density.
            if all([
                a_2D_mIDs[row][col] > 0,
                a_2D_density[row][col] == 0]):
                ## todo: generate map with and without!!!
                if all([
                    latitude > tripoint_SW_lat,
                    longitude > tripoint_SW_lon,
                    latitude < tripoint_NE_lat,
                    longitude < easternmost_lon,
                    ]):
                    if args.verbose:
                        print('fix desert', latitude, longitude, row, col)
                    a_2D_density[row][col] = 1
                else:
                    print('tmp skip', row, col, latitude, longitude)
                    continue

    return a_2D_density


def color_as_nearby(args, a_2D_density, a_2D_mIDs):

    assert a_2D_density.shape == a_2D_mIDs.shape

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))

    ## http://en.wikipedia.org/wiki/List_of_countries_by_southernmost_point
    ## http://en.wikipedia.org/wiki/List_of_countries_by_westernmost_point
    Yemen_South = 12+7/60  # 12.66
    Yemen_West = 41+42/60  # 43.12  # 43.42
    Yemen_South1 = 12+(33.75-0.01)/60
    Yemen_West1 = 43+(18.75-0.01)/60
    Yemen_South2 = 13+(41.25-0.01)/60
    Yemen_West2 = 42+(33.75-0.01)/60
    Israel_South = 29.5
    Israel_West = 34+17/60
    Gaza_West = 34+15/60
    ## http://en.wikipedia.org/wiki/Extreme_points_of_Spain
    Spain_South = 36+0/60
    Spain_East = -4.75  # -(3+19/60)
    ## http://en.wikipedia.org/wiki/List_of_countries_by_easternmost_point
    SaudiArabia_South = 16+5/60
    SaudiArabia_South = 17+(3.75-0.01)/60
    SaudiArabia_East = 40+(41.25-0.01)/60
##    Sudan_East = 38+35/60
    ## http://en.wikipedia.org/wiki/Sharm_el-Sheikh
    Sharm_elSheikh_lon = 34+19/60  # Suez Canal East boundary
    # http://en.wikipedia.org/wiki/Ras_Muhammad_National_Park
    Sinai_South_lat = 27+43/60  # Gulf of Aqaba South boundary
    Sinai_South_lon = 34+15/60  # Gulf of Aqaba West boundary
    Egypt_South = 22
    ## http://en.wikipedia.org/wiki/Suez
    Suez_lat = 29+58/60  # Suez Canal North Boundary
    Suez_lon = 32+33/60  # Suez Canal West Boundary
    ## http://en.wikipedia.org/wiki/Safaga
    Safaga_lat = 26+44/60  # Suez Canal South boundary

    ## Canary Islands
    ## http://en.wikipedia.org/wiki/Roque_del_Este
    Canary_Islands_lon = -13-20/60

    ## Strait of Sicily
    Libya_N_lat = 33+10/60
    Tunisia_E_lon = 11+35/60

    ## Setting zero density if no color is not a good solution,
    ## because lakes will randomly appear and disappear,
    ## because the polygons are poorly defined.

    dist = 3
    for x in range(dist, a_2D_density.shape[0]-dist):
        row = x
        lat = -cellsize*(row+0.5)*max(nrows, ncols)/args.n-yllcorner
        lat += cellsize*(ncols-nrows)/2
        for z in range(dist, a_2D_density.shape[1]-dist):
            col = z
            lon = cellsize*(col+0.5)*max(nrows, ncols)/args.n+xllcorner
            ## Popuation density is zero.
            if a_2D_density[x][z] == 0:
                continue
            ## Point is in language polygon.
            if a_2D_mIDs[x][z] > 0:
                continue
            ## If not densely populated and not covered by language polygon
            ## then delete. Because probably near coast or lake.
            ## Including this deletion does not cause problems in Egypt.
            ## todo: generate map with and without!!!
            if a_2D_density[x][z] == 1:
##                if args.verbose:
##                    print(
##                        'near zero density (e.g. outside Africa and coast)',
##                        x, z, lat, lon)
                a_2D_density[x][z] = 0
                continue
##            if lat > Yemen_South and lon > Israel_West:
##                dist = 1
##            else:
##                dist = 2
            ## Count nearby materialIDs.
            cnt = collections.Counter(
                a_2D_mIDs[x+dx][z+dz]
                for dx in range(-dist, dist+1) for dz in range(-dist, dist+1)
                )
            ## Assign greater weight to most nearby.
            cnt += collections.Counter(
                a_2D_mIDs[x+dx][z+dz]
                for dx in range(-1, 1+1) for dz in range(-1, 1+1)
                for i in range(dist**2)
                )
            del cnt[0]
            ## No nearby polygons/colors.
            if len(cnt) == 0:
                continue
            ## No bricks in Yemen.
##            if lat > Yemen_South and lon > Yemen_West:
##                if args.verbose:
##                    print('Yemen',x,z,lat,lon)
##                a_2D_mIDs[x][z] = 221
##                continue
            ## Do not connect colors/polygons across the Red Sea.
            if lat >= Yemen_South and lon >= Israel_West and a_2D_density[x+1][z] == 0:
                a_2D_density[x][z] = 0
                continue
            if lat >= Yemen_South1 and lon >= Yemen_West1:
                if args.verbose:
                    print('Yemen1',x,z,lat,lon)
                a_2D_density[x][z] = 0
                continue
            if lat >= Yemen_South2 and lon >= Yemen_West2:
                if args.verbose:
                    print('Yemen2',x,z,lat,lon)
                a_2D_density[x][z] = 0
                continue
            ## No bricks in Israel and Gaza.
            if lat > Israel_South and lon > Gaza_West:
                if args.verbose:
                    print('Israel',x,z,lat,lon)
                a_2D_density[x][z] = 0
                continue
            ## No bricks in Saudi Arabia across the Red Sea.
            if lat > SaudiArabia_South and lon > SaudiArabia_East:
                if args.verbose:
                    print('Saudi',x,z,lat,lon)
##                a_2D_mIDs[x][z] = 221
                a_2D_density[x][z] = 0
                continue
            ## No bricks in Saudi Arabia across the Gulf of Aqaba.
            if all([
                lat > Sinai_South_lat,
                lon > Sinai_South_lon,
                a_2D_density[x][z] < args.layers/3]):
##                a_2D_mIDs[x][z] = 102  # medium blue
                a_2D_density[x][z] = 0
                continue
            ## No bricks in Saudi Arabia across the Gulf of Aqaba.
            if all([
                lat > Egypt_South,
                lon > Gaza_West,
                any([a_2D_density[x][z-1] == 0, a_2D_density[x][z-2] == 0])]):
##                a_2D_mIDs[x][z] = 221
                a_2D_density[x][z] = 0
                continue
            ## No Bricks in Andalusia.
            if lat > Spain_South and lon < Spain_East:
##                a_2D_mIDs[x][z] = 221
                a_2D_density[x][z] = 0
                continue
            ## No Bricks in Andalusia.
            if lat > Spain_South and a_2D_density[x+1][z] == 0:
##                a_2D_mIDs[x][z] = 221
                a_2D_density[x][z] = 0
                continue
            ## Clear bricks in the Suez canal:
            if all([
                lat > Safaga_lat,
                lat < Suez_lat,
                lon > Suez_lon,
##                lon < Sinai_South_lon,
                lon < Sharm_elSheikh_lon,
                a_2D_density[x][z] < args.layers/3]):
##                a_2D_mIDs[x][z] = 221
                a_2D_density[x][z] = 0
                continue
##            ## Clear brick(s) in the Red Sea south of the Sinai peninsula.
##            if latitutde > Safaga_lat and lat
            ## Skip Canary Islands
            if lon < Canary_Islands_lon and a_2D_mIDs[x][z+1] == 0:
                a_2D_density[x][z] = 0
                continue
            ## No bricks in the strait of Sicily.
            if lat > Libya_N_lat and lon > Tunisia_E_lon:
                a_2D_density[x][z] = 0
                continue
            ## Set materialID to most frequently occuring nearby materialID.
            a_2D_mIDs[x][z] = cnt.most_common(1)[0][0]
##            a_2D_mIDs[x][z] = 5
##            print(cnt.most_common(1)[0][0])

    return a_2D_density, a_2D_mIDs


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--affix', default='afds00g',
        choices=['afp00ag', 'afp00g', 'afds00ag', 'afds00g'])

    parser.add_argument(
        '--plate_size', default=48, type=int, choices=[16, 32, 48],
        help='Size of base plates to build on.')

    parser.add_argument(
        '--plate_cnt', default=5, type=int,
        help='Number of base plates along one dimension.')

    parser.add_argument(
        '--zero', help='add zero values to average?', action='store_true')

    parser.add_argument(
        '--verbose', help='Be verbose?', action='store_true')

    parser.add_argument(
        '--density', help='Take max or mean of density grid?', choices=['max', 'mean'], default='max')

    parser.add_argument('--layers', default=27, type=int)

    parser.add_argument('--norm', default='log', choices=[
        'log10', 'log2', 'unity', 'log'])  # unity is feature scaling

    parser.add_argument('--colors', required=True)

    args = parser.parse_args()

    # script fast if multiple of 2160
    args.n = n = nrows = ncols = args.plate_cnt*args.plate_size

    assert args.layers % 3 == 0

    return args


def asc2np(args, affix1):

    a_gpw = read_gpw('{}.asc'.format(args.affix))
    nrows2, ncols2 = np.shape(a_gpw)

    assert np.shape(a_gpw)[0] > args.n

    den = args.n/fractions.gcd(ncols2, args.n)
    num = int(ncols2/fractions.gcd(ncols2, args.n))

    a_2D_density = np.zeros((args.n, args.n))
    a_2D_density_cnt = np.zeros((args.n, args.n))

    print('converting NASA array to LEGO array')
    for x1 in range(args.n):
        print('gpw2LEGO', x1, args.n)
        for y1 in range(args.n):
            l = []
            for x2 in range(num):
                for y2 in range(num):
                    x3 = int((x1*num+x2)/den)
                    y3 = int((y1*num+y2)/den)
##                    a_2D_density[x1][y1] += a_gpw[x3][y3]
##                    a_2D_density[x1][y1] += a_gpw[x3][y3]
                    l.append(a_gpw[x3][y3])
                    if a_gpw[x3][y3] > 0:
                        a_2D_density_cnt[x1][y1] += 1
            if l.count(0) >= 0 and sum(l) > 0:
                if y1 == 65 and x1 > 95:
                    print(x1, y1, l.count(0), 'mean', sum(l)/len(l), 'max', max(l), sum(l), sum(sorted(l)[1:-1]), l)
            a_2D_density[x1][y1] = max(l)

    np.save(affix1, a_2D_density)

    del a_gpw
    del a_2D_density_cnt

    return a_2D_density


def find_connected_exposed(
    args, a_2D_density, a_2D_buried,
    a_2D_mIDs, a_3D_dIDs):

    print('find connected exposed plates and remaining buried plates')

    a_3D_mIDs = np.zeros(np.shape(a_3D_dIDs), int)
    a_3D_angle = np.zeros(np.shape(a_3D_dIDs), int)

    gap = 'x'
    ## Use symbol for buried bricks instead of materialID,
    ## because surface exposed bricks/plates might have same color as
    ## buried bricks/plates.
    buried = 'o'

    for layer in reversed(range(args.layers)):
        ## build plates horizontally and vertically in each layer
        if layer % 2 == 0:
            irow = 0
            jrow = 1
            icol = 1
            jcol = 0
            pass
        else:
            irow = 1
            jrow = 0
            icol = 0
            jcol = 1
            pass
        if layer % 3 == 2:
            h = 3  # bricks
            layer_insert = layer-2
            layer_remove = (layer, layer-1)
            max_len = 8  # 1x8 (3008) available in few colors
        else:
            h = 1  # plates
            layer_insert = layer
            layer_remove = ()
            max_len = 8

        ## loop baseplates from North to South
        for row1 in range(args.n//args.plate_size):
            ## loop baseplates from West to East
            for col1 in range(args.n//args.plate_size):
                for i in range(args.plate_size):
                    seq = []
                    for j in range(args.plate_size):
                        row = row1*args.plate_size + i*irow+j*jrow
                        col = col1*args.plate_size + i*icol+j*jcol
                        ## Ocean or lake (or desert).
                        ## Or position is above structure.
                        if layer >= a_2D_density[row][col]:
                            seq += [gap]
                        ## Already filled.
                        elif a_3D_dIDs[layer][row][col] != 0:
                            seq += [gap]
                        ## Buried.
                        elif layer < a_2D_buried[row][col]:
                            seq += [buried]
                        ## Not inside Felix2001 GeoJSON polygon
                        ## e.g. Spain and Saudi Arabia
                        elif a_2D_mIDs[row][col] == 0:
                            seq += [gap]
                        else:
                            seq += [int(a_2D_mIDs[row][col])]
                        ## Continue loop over j.
                        continue
                    ## No bricks along line.
                    if seq == args.plate_size*[gap]:
                        continue
                    ## Same color for entire sequence.
                    if seq == args.plate_size*[seq[0]]:
                        pass
                    else:
                        seq = find_consecutive(args, seq, h, gap=gap, buried=buried)
                    append_designID_materialID_main(
                        args, seq, layer, i, irow, jrow, icol, jcol,
                        max_len, row1, col1,
                        a_3D_dIDs, a_3D_mIDs, a_3D_angle,
                        h, layer_insert, layer_remove,
                        gap, buried)
                    ## Continue loop over i.
                    continue
                ## Continue col1 loop.
                continue
            ## Continue row1 loop.
            continue
        ## Continue layer loop.
        continue

    return a_3D_dIDs, a_3D_mIDs, a_3D_angle


def append_designID_materialID_main(
    args, seq, layer, i, irow, jrow, icol, jcol, max_len, row1, col1,
    a_3D_dIDs, a_3D_mIDs, a_3D_angle,
    h, layer_insert, layer_remove, gap, buried):

    width = 1

    pos = 0
    for materialID, g in itertools.groupby(seq):
        ## Get length of iterator.
        len_group = len(list(g))
        if materialID == gap:
            pos += len_group
            continue
        ## Get materialID of buried bricks/plates.
        if materialID == buried:
            materialID = materialID_buried

        ## 1x1 brick not available, but 1x1 plate available.
        if h == 3 and materialID in d_NA[h][1][1] and materialID not in d_NA[1][1][1]:
            ## 1x1 plate if 1x1 brick not available
            ## e.g. cool yellow and spring yellowish green / unikitty green
            h_group = 1
            layer_ins = layer
            layer_rem = ()
        ## 1x1 plate not available, but 1x1 brick available.
        elif layer >= 1 and materialID in d_NA[h][1][1] and not materialID in d_NA[3][1][1]:
            ## 1x1 brick if 1x1 plate not available
            ## e.g. medium azure, medium lilac, dark pink, sand green
            h_group = 3
            layer_ins = layer-2
            layer_rem = (layer, layer-1)
        ## 1x1 plate and 1x1 brick available.
        else:
            h_group = h
            layer_ins = layer_insert
            layer_rem = layer_remove

        ## What is the longest length of this color and size?
        try:
            max_len = min(8, d_max_len[h_group][width][materialID])
        except KeyError:
            for max_len in (8, 6, 4, 3, 2, 1):
                if materialID in d_NA[h_group][width][max_len]:
                    continue
                break

        ## How many plates of max length will fit?
        ## 1st call of sub routine.
        for k in range(len_group//max_len):
            for length in (max_len,):
                designID = d_len2dIDs[h_group][width][length]
                pos = append_designID_materialID_sub(
                    layer, designID, materialID,
                    a_3D_dIDs, a_3D_mIDs, a_3D_angle,
                    pos, i, irow, jrow, icol, jcol, length, width,
                    row1, col1, args,
                    layer_insert=layer_ins, layer_remove=layer_rem,
                    h=h_group)
##                if a_3D_dIDs[3][127][113] == 3022 or a_3D_mIDs[3][127][113] == 330:
##                    stop3
        ## How much space left after filling with plates of max length?
        mod = len_group % max_len
        ## No space left.
        if mod == 0:
            continue
        assert max_len <= 8
        if mod == 7:
            if not any([
                materialID in d_NA[h_group][width][4],
                materialID in d_NA[h_group][width][3],
                ]):
                lengths = (4, 3)
            elif not materialID in d_NA[h_group][width][6]:
                lengths = (6, 1)
            elif not materialID in d_NA[h_group][width][4]:
                lengths = (4, 2, 1)
            elif not materialID in d_NA[h_group][width][3]:
                lengths = (3, 2, 2)
            else:
                lengths = (2, 2, 2, 1)
        elif mod == 5:
            if not any([
                materialID in d_NA[h_group][width][3],
                materialID in d_NA[h_group][width][2],
                ]):
                lengths = (3, 2)
            elif not materialID in d_NA[h_group][width][4]:
                lengths = (4, 1)
            elif not materialID in d_NA[h_group][width][2]:
                lengths = (2, 2, 1)
            else:
                lengths = (1, 1, 1, 1, 1)
        else:
            if mod > 1 and materialID in d_NA[h_group][width][mod]:
                ## 2, 4, 6
                if mod not in (1, 3):
                    assert mod % 2 == 0
                    if not materialID in d_NA[h_group][width][mod]:
                        lengths = int(mod / 2) * (2,)
                    else:
                        lengths = int(mod / 1) * (1,)
                ## 1, 3
                else:
                    lengths = mod * (1,)
            else:
                lengths = (mod,)
        ## 2nd call of sub routine.
        for length in lengths:
            designID = d_len2dIDs[h_group][width][length]
            pos = append_designID_materialID_sub(
                layer, designID, materialID,
                a_3D_dIDs, a_3D_mIDs, a_3D_angle,
                pos, i, irow, jrow, icol, jcol, length, width,
                row1, col1, args,
                layer_ins, layer_rem, h_group)
##            if a_3D_dIDs[3][127][113] == 3022 or a_3D_mIDs[3][127][113] == 330:
##                print(max_len, mod, lengths, length)
##                stop2
        ## Continue loop over materialID groups.
        continue

##    if a_3D_dIDs[3][127][113] == 3022 or a_3D_mIDs[3][127][113] == 330:
##        stop1
##    if a_3D_dIDs[6][103][105] == 3007 or a_3D_mIDs[6][103][105] == 28:
##        stop1

    return


def append_designID_materialID_sub(
    layer, designID, materialID,
    a_3D_dIDs, a_3D_mIDs, a_3D_angle,
    pos, i, irow, jrow, icol, jcol, length, width, row1, col1, args,
    layer_insert, layer_remove, h):

    for j in range(length):
        row = row1*args.plate_size + i*irow+(pos+j)*jrow
        col = col1*args.plate_size + i*icol+(pos+j)*jcol
        ## Odd layer.
        if layer % 2 == 1 and j == length - 1:
            a_3D_dIDs[layer_insert][row][col] = designID
            a_3D_mIDs[layer_insert][row][col] = materialID
            ## Can we replace with double width piece?
            if all([
                ## Not the first row.
                row % args.plate_size > 0,
                ## Previous brick/plate identical to current one.
                a_3D_dIDs[layer_insert][row-1][col] == designID,
                ## Previous color identical or buried.
                any([
                    materialID == materialID_buried,
                    a_3D_mIDs[layer_insert][row-1][col] in (
                        materialID, materialID_buried)]),
                ## Same rotation as brick/plate to be connected with.
                a_3D_angle[layer_insert][row][col] == a_3D_angle[layer_insert][row-1][col],
                ## Color available.
                all([
                    materialID not in d_NA[h][2][length],
                    a_3D_mIDs[layer_insert][row-1][col] not in
                    d_NA[h][2][max(length, width)],
                    ]),
                ]):
                a_3D_dIDs[layer_insert][row][col] = d_len2dIDs[h][2][length]
                a_3D_dIDs[layer_insert][row-1][col] = -1
                ## Do not extend buried color to edge.
                if materialID == materialID_buried:
                    a_3D_mIDs[layer_insert][row][col] = a_3D_mIDs[layer_insert][row-1][col]
                a_3D_mIDs[layer_insert][row-1][col] = 0
                if length == 1:
                    a_3D_angle[layer_insert][row][col] = 1
                if any([
                    a_3D_dIDs[3][127][113] == 3022,
                    a_3D_mIDs[3][127][113] == 330,
                    a_3D_dIDs[6][103][105] == 3007,
                    a_3D_mIDs[6][103][105] == 28,
                    ]):
                    print('h', h, 'l', length, 'w', width)
                    print('mID', materialID)
                    print('NA', d_NA[h][2][max(length, width)])
                    print('max(l, w)', max(length, width))
                    print('rowcol', row, col)
                    print('layer', layer, layer_insert)
                    for x in range(27):
                        print(x, a_3D_mIDs[x][127][113])
                    stop4a
##            if a_3D_dIDs[3][127][113] == 3022 or a_3D_mIDs[3][127][113] == 330:
##                print(h, length)
##                stop4b
        ## Even layer.
        elif layer % 2 == 0 and j == length - 1:
            a_3D_dIDs[layer_insert][row][col] = designID
            a_3D_mIDs[layer_insert][row][col] = materialID
            ## Can we replace with double width piece?
            if all([
                ## Not the first col.
                col % args.plate_size > 0,
                ## Previous brick/plate identical to current one.
                a_3D_dIDs[layer_insert][row][col-1] == designID,
                ## Previous color identical or buried.
                any([
                    materialID == materialID_buried,
                    a_3D_mIDs[layer_insert][row][col-1] in (
                        materialID, materialID_buried)]),
                ## Same rotation as brick/plate to be connected with.
                a_3D_angle[layer_insert][row][col] == a_3D_angle[layer_insert][row][col-1],
                ## Color available.
                all([
                    materialID not in d_NA[h][2][length],
                    a_3D_mIDs[layer_insert][row][col-1] not in
                    d_NA[h][2][max(length, width)],
                    ]),
                ]):
                a_3D_dIDs[layer_insert][row][col] = d_len2dIDs[h][2][length]
                a_3D_dIDs[layer_insert][row][col-1] = -1
                ## Do not extend buried color to edge.
                if materialID == materialID_buried:
                    a_3D_mIDs[layer_insert][row][col] = a_3D_mIDs[layer_insert][row][col-1]
                a_3D_mIDs[layer_insert][row][col-1] = 0
                if length == 1:
                    a_3D_angle[layer_insert][row][col] = 1
                if any([
                    a_3D_dIDs[3][127][113] == 3022,
                    a_3D_mIDs[3][127][113] == 330,
                    a_3D_dIDs[6][103][105] == 3007,
                    a_3D_mIDs[6][103][105] == 28,
                    ]):
                    print('h', h, 'l', length, 'w', width)
                    print('mID', materialID)
                    print('NA', d_NA[h][2][max(length, width)])
                    print('max(l, w)', max(length, width))
                    print('rowcol', row, col)
                    print('layer', layer, layer_insert)
                    for x in range(27):
                        print(x, a_3D_mIDs[x][127][113])
                    stop5a
##            if a_3D_dIDs[3][127][113] == 3022 or a_3D_mIDs[3][127][113] == 330:
##                print(h, length)
##                stop5b
        else:
            a_3D_dIDs[layer_insert][row][col] = -1
            a_3D_mIDs[layer_insert][row][col] = 0
        for layer2 in layer_remove:
            a_3D_dIDs[layer2][row][col] = -1
            a_3D_mIDs[layer2][row][col] = 0

    pos += length

    return pos


def find_consecutive(args, seq, h, gap='x', buried='o'):

    n = len(seq)

    groups = [list(g) for k, g in itertools.groupby(seq)]
    seq2 = []
    for i, g in enumerate(groups):
        ## Append gap.
        if g[0] == gap:
            seq2 += g
        ## Append non-buried.
        elif g[0] != buried:
            seq2 += g
        ## Append buried as buried or non-buried.
        ## elif g[0] == buried
        else:
            ## first group
            if i == 0:
                ## Next group is gap.
                if groups[i+1][0] == gap:
                    seq2 += g
                else:
                    ## Color first group same as next group.
                    if groups[i+1][0] in ubiquitous:
                        seq2 += len(g)*[groups[i+1][0]]
                    ## Color first group as buried.
                    else:
                        seq2 += g
            ## last group
            elif i+1 == len(groups):
                ## previous group is gap
                if groups[i-1][0] == gap:
                    seq2 += g
                ## color last group same as previous group
                else:
                    if groups[i-1][0] in ubiquitous:
                        seq2 += len(g)*[groups[i-1][0]]
                    else:
                        seq2 += g
            ## if end of sequence then continue previous color
            elif len(seq2)+len(g) == n:
                seq2 += len(g)*[groups[i-1][0]]
                print(seq, seq2, g)
                stop7
            ## Gap before and after. Append buried.
            elif groups[i-1][0] == gap and groups[i+1][0] == gap:
                seq2 += g
            ## gap before but not after
            elif groups[i-1][0] == gap and groups[i+1][0] != gap:
                if groups[i+1][0] in ubiquitous:
                    seq2 += len(g)*[groups[i+1][0]]
                else:
                    seq2 += g
            ## gap after but not before
            elif groups[i-1][0] != gap and groups[i+1][0] == gap:
                if groups[i-1][0] in ubiquitous:
                    seq2 += len(g)*[groups[i-1][0]]
                else:
                    seq2 += g
            ## Same color before and after.
            elif groups[i-1][0] == groups[i+1][0]:
                if groups[i-1][0] in ubiquitous:
                    seq2 += len(g)*[groups[i-1][0]]
                else:
                    seq2 += g
            ## Different color before and after.
            else:
                len_prev = len(list(itertools.takewhile(
                    lambda x: x == groups[i-1][0], reversed(seq2))))
                ## Length 2 or longer to be appended.
                if len(g) >= 2:
                    if all([
                        groups[i-1][0] in ubiquitous,
                        groups[i+1][0] not in ubiquitous]):
                        seq2 += len(g)*[groups[i-1][0]]
                    elif all([
                        groups[i-1][0] not in ubiquitous,
                        groups[i+1][0] in ubiquitous]):
                        seq2 += len(g)*[groups[i+1][0]]
                    elif len_prev+len(g)//2 in (1, 2, 3, 4, 6, 8):  # PAB lengths
                        ## append with color from both sides
                        seq2 += (len(g)//2)*[groups[i-1][0]]
                        seq2 += (len(g)-len(g)//2)*[groups[i+1][0]]
                        continue
                    else:
                        seq2 += (len(g)//2+1)*[groups[i-1][0]]
                        seq2 += (len(g)-len(g)//2-1)*[groups[i+1][0]]
                        continue
                else:  # elif len(g) == 1:
                    ## append to previous sequence of length 1
                    try:
                        NA = groups[i-1][0] in d_NA[h][1][len_prev+len(g)]
                    except:
                        NA = False
                    if all([
                        groups[i-1][0] in ubiquitous,
                        groups[i+1][0] not in ubiquitous]):
                        seq2 += len(g)*[groups[i-1][0]]
                    elif all([
                        groups[i-1][0] not in ubiquitous,
                        groups[i+1][0] in ubiquitous]):
                        seq2 += len(g)*[groups[i+1][0]]
                    elif len_prev == 1 and not NA:
                        seq2 += 1*[groups[i-1][0]]
                    elif all([
                        ## before has length greater than 1
                        len_prev > 1,
                        ## after has length equal to 1
                        len(list(itertools.takewhile(
                            lambda x: x == groups[i+1][0], seq[len(seq2)+1:]))) == 1,
                        ]):
                        seq2 += 1*[groups[i+1][0]]
                    elif len_prev+len(g) in (2, 3, 4, 6, 8):  # PAB lengths
                        seq2 += len(g)*[groups[i-1][0]]
                    else:
                        seq2 += len(g)*[groups[i+1][0]]
                    pass
                pass
            pass
        continue

    assert len(seq2) == n

    return seq2


def largest_empty_rectangle(
    args, layer, r1, c1,
    a_2D_buried, a_3D_dIDs, a_2D_density, a_2D_mIDs, zero=0):

    ## http://en.wikipedia.org/wiki/Largest_empty_rectangle

    '''Solve the largest empty rectangle problem in O(mn)'''

    d_rectangle = {'area': 0, 'nrows': 0, 'ncols': 0}
    h = np.zeros((args.plate_size, args.plate_size), int)  # heights
    w = np.zeros((args.plate_size, args.plate_size), int)  # widths
    ## North to South.
    for r2 in range(args.plate_size):
        r3 = r1*args.plate_size+r2
        ## Only consider even numbered grid points.
        ## West to East.
        for c2 in range(args.plate_size):
            c3 = c1*args.plate_size+c2
            ## Skip if no plate in layer.
            if layer >= a_2D_buried[r3][c3]:
                continue
            ## Skip if occupied already (either by material or empty gap);
            ## i.e. a layer i in which a brick was positioned
            ## when looping over layer i+2.
            if a_3D_dIDs[layer][r3][c3] != zero:
                continue
            ## Skip if there is room for ubiquitous non-buried/colored brick above.
            if layer % 3 != 2 and 3*(layer//3) + 2 <= a_2D_density[r3][c3]:
                if a_2D_mIDs[r3][c3] in ubiquitous:
                    continue
            ## first row
            if r2 == 0:
                h[r2][c2] = 1
            ## append to previous row
            else:
                h[r2][c2] = h[r2-1][c2]+1
            ## first col
            if c2 == 0:
                w[r2][c2] = 1
            ## append to previous col
            else:
                w[r2][c2] = w[r2][c2-1]+1

##            ## Only consider even numbered grid points.
##            ## To avoid rotating at edges.
##            if r2 % 2 == 0 or c2 % 2 == 0:
##                continue
            ## Only consider every 4th grid point.
            ## To avoid rotating 2 stud bricks at edges.
            if r2 % 4 != 3 or c2 % 4 != 3:
                continue
##            if layer % 2 == 0:
##                if r2 % 4 != 3 or c2 % 2 != 1:
##                    continue
##            else:
##                if r2 % 2 != 1 or c2 % 4 != 3:
##                    continue

            min_w = w[r2][c2]
            for dh in range(h[r2][c2]):
                min_w = min(min_w, w[r2-dh][c2])
                ## Don't append area if 4x4 brick can't fit inside it.
                if min_w < 4:
                    break
                if dh+1 < 4:
                    continue
                ## Only consider increments of 2.
                if (dh+1) % 2 != 0 or min_w % 2 != 0:
                    continue
##                ## Only consider rectangles of relevant size; e.g. ignore 1x.
##                if dh+1 not in d_dIDs_sq_plate.keys():
##                    continue
##                if min_w not in d_dIDs_sq_plate.keys():
##                    continue
                ## Calculate area.
                area = min_w*(dh+1)
                ## Find largest area.
                if area > d_rectangle['area']:
                    d_rectangle = {
                        'area': area, 'r2': r2, 'c2': c2,
                        'nrows': dh+1, 'ncols': min_w}

            min_h = h[r2][c2]
            for dw in range(w[r2][c2]):
                min_h = min(min_h, h[r2][c2-dw])
                ## Don't append area if 4x4 brick can't fit inside it.
                if min_h < 4:
                    break
                if dw+1 < 4:
                    continue
                ## Only consider increments of 2.
                if (dw+1) % 2 != 0 or min_h % 2 != 0:
                    continue
##                ## Only consider rectangles of relevant size; e.g. ignore 1x.
##                if min_h not in d_dIDs_sq_plate.keys():
##                    continue
##                if dw+1 not in d_dIDs_sq_plate.keys():
##                    continue
                ## Calculate area.
                area = (dw+1)*min_h
                ## Find largest area.
                if area > d_rectangle['area']:
                    d_rectangle = {
                        'area': area, 'r2': r2, 'c2': c2,
                        'nrows': min_h, 'ncols': dw+1}

    return d_rectangle


def find_connected_buried(args, a_2D_buried, a_2D_density, a_2D_mIDs):

    print('''Finding connected 1x1 plates
and replacing them with larger plates.''')

    ## todo: get rid of modulo in this function since only looping every third layer now

    ## 3D array with design IDs at each 3D coordinate
    a_3D_dIDs = np.zeros((args.layers, args.n, args.n), int)
    size_min = min(d_dIDs_sq_plate.keys())

    ## Loop baseplates from North to South
    for r1 in range(int(args.n/args.plate_size)):
        ## Loop baseplates from West to East
        for c1 in range(int(args.n/args.plate_size)):
            if args.verbose:
                print('find_connected_buried', r1, c1)
            ## Loop from top to bottom. Substituting with bricks first
            ## instead of plates if possible.
            for layer in range(args.layers-1, -1, -3):
                while True:
                    d_rectangle = largest_empty_rectangle(
                        args, layer, r1, c1,
                        a_2D_buried, a_3D_dIDs, a_2D_density, a_2D_mIDs)
                    ## No more substitution to be made.
                    if any([
                        d_rectangle['area'] == 0,
                        d_rectangle['nrows'] < size_min,
                        d_rectangle['ncols'] < size_min]):
                        break
##                    if args.verbose:
##                        print('rectangle', r1, c1, layer, d_rectangle)
                    ## Calculate size of area to be filled.
                    for size in reversed(sorted(d_dIDs_sq_plate.keys())):
                        if all([
                            size <= d_rectangle['nrows'],
                            size <= d_rectangle['ncols']]):
                            break
                    if size > d_rectangle['nrows'] or size > d_rectangle['ncols']:
                        stopshouldnothappen
##                    print('size', size)
                    ## Replace with brick(s) instead of plate every third layer.
                    ## Hence calculate modulo.
                    mod = layer % 3
                    ## calculate number of square plates fitting inside area
                    plate_rows = d_rectangle['nrows']//size
                    plate_cols = d_rectangle['ncols']//size
                    r3 = r1*args.plate_size+d_rectangle['r2']
                    c3 = c1*args.plate_size+d_rectangle['c2']
                    ## append plate and empty spaces to 3D array
                    for row_major in range(plate_rows):
                        for col_major in range(plate_cols):
                            ## empty spaces
                            for row_minor in range(size):
                                row = r3-row_major*size-row_minor
                                for col_minor in range(size):
                                    col = c3-col_major*size-col_minor
                                    ## Make empty space for larger plate.
                                    a_3D_dIDs[layer][row][col] = -1
                                    ## Make empty space for bricks
                                    ## for every third layer from bottom.
                                    if mod == 2:
                                        for layer_below in range(
                                            layer-1, layer-3, -1):
                                            a_3D_dIDs[
                                                layer_below][row][col] = -1
                                    continue
                                continue
                            ## plate starts SE and extends to NW
                            ## i.e. from high row to low row
                            ## and from low col to high col
                            row = r3-row_major*size
                            col = c3-col_major*size
                            if mod != 2:
                                ## Insert large plate.
                                designID = d_dIDs_sq_plate[size]
                                a_3D_dIDs[layer][row][col] = designID
                            else:
                                ## Insert brick(s).
                                designID = d_len2dIDs[3][2][size]
                                dx = 2*(layer % 2)
                                dz = 2*(1-layer % 2)
                                for i in range(size//2):
                                    a_3D_dIDs[layer-2][row-i*dx][col-i*dz] = designID
##                                for dx, dz, designID in d_brick_replace[size]:
##                                    a_3D_dIDs[
##                                        layer-2][row-dx][col-dz] = designID
                            ## Continue loop over col_major
                            continue
                        ## Continue loop over row_major.
                        continue
                    ## Continue while loop.
                    continue
                ## Continue layer loop.
                continue
            ## Continue c1 loop.
            continue
        ## Continue r1 loop.
        continue

    return a_3D_dIDs


def find_buried(args, a_2D_density, layers, a_2D_mIDs):

    a_2D_buried = np.zeros(np.shape(a_2D_density), int)

    for row in range(args.n):
        ## nothing buried at the edge (no bricks at the edge anyway)
        if row == 0 or row == args.n-1:
            continue
        for col in range(args.n):
            ## nothing buried at the edge (no bricks at the edge anyway)
            if col == 0 or col == args.n-1:
                continue
            ## no buried plates if only 1 plate
            if a_2D_density[row][col] <= 1:
                continue
            ## minimum neighbouring height
            dens_min = z = min([
                a_2D_density[row+x][col+y]
                for x in range(-1, 1+1) for y in range(-1, 1+1)])
            ## 1x1 plate not available
            ## but 1x1 brick available.
            ## 268 medium lilac (dark purple)
            ## 221 bright purple (dark pink)
            ## 322 medium azure
            ## 151 sand green
            if all([
                ## 1x1 plate not available.
                a_2D_mIDs[row][col] in d_NA[1][1][1],
##                ## 1x1 brick available.
##                a_2D_mIDs[row][col] not in d_NA[3][1][1],
                ]):
                n = a_2D_density[row][col]
                m = dens_min
                dens_min = min(
                    dens_min + 1,
                    n - 2,
##                    n - 3 * ( (n - m + 1 + 1) // 3) + 1,
                    n - 3 * ( (n - (m + 1) + 3) // 3) + 1,
                    )
            ## update array
            a_2D_buried[row][col] = max(0, dens_min-1)

    return a_2D_buried


def read_gpw_header(file_gpw):

    with open(file_gpw) as f:
        ncols = int(f.readline().strip().split()[1])
        nrows = int(f.readline().strip().split()[1])
        assert ncols % 2 == 0 and nrows % 2 == 0
        assert ncols > nrows
        assert ncols == 2160  # -30+2160*2.5/60 = 60 = 60E
        assert nrows == 1920  # -40+1920*2.5/60 = 40 = 40N
        xllcorner = int(f.readline().strip().split()[1])
        yllcorner = int(f.readline().strip().split()[1])
        assert xllcorner == -30  # -30 = 30E
        assert yllcorner == -40  # -40 = 40S
        cellsize = float(f.readline().strip().split()[1])
        assert round(cellsize, 8) == round(2.5/60, 8)

    assert round(yllcorner+cellsize*nrows, 0) == 40
    assert round(xllcorner+cellsize*ncols, 0) == 60

    return ncols, nrows, xllcorner, yllcorner, cellsize


def normalize(a, args):

    n = np.shape(a)[0]

    ## Plot population density histogram prior to normalisation.
    histogram1(a, args)

    ## https://en.wikipedia.org/wiki/Normalization_(statistics)
    if args.norm == 'log10':
        amax = np.amax(a)
        a = np.log10(np.divide(a, amax/(10**layers)))
    elif args.norm == 'log':
        amax = np.amax(a)
        for row in range(n):
            for col in range(n):
                if a[row][col] == 0:
                    continue
                log = math.log(a[row][col])
                a[row][col] = max(1, args.layers*log/math.log(amax))
    elif args.norm == 'unity':
        amin = 0
        amax = np.amax(a)
        a = np.divide(a, amax/args.layers)
    elif args.norm == 'log2':
        amax = np.amax(a)
        a = np.log2(
            np.divide(a, amax/(2**args.layers)))
    else:
        stop

    assert math.ceil(np.amax(a)) == args.layers

    for row in range(n):
        for col in range(n):
            if a[row][col] < 0:
                a[row][col] = 1
            else:
                a[row][col] = math.ceil(a[row][col])

    ## Plot population density histogram after normalisation.
    histogram2(a, args)

    return a


def json2array(args, a_2D_density):

    ## 1) Calculate polygon intersection area for all 1x1 grid points.
    ## 2) Calculate largest intersection area (ID) for each 1x1 grid point.
    ## 3) Assign colors to 1x1 grid points of each ID.

##    AREA = 0

    ## http://en.wikipedia.org/wiki/The_Languages_of_Africa

    print('Converting GeoJSON polygons to array of LEGO color points')

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))

    d_family2color = assign_color_families_to_language_families(args)

    _den = cellsize * max(nrows, ncols)
    factor = _den/args.n

    a_2D_mIDs = np.zeros((args.n, args.n), int)
    a_2D_areas = np.empty((args.n, args.n), dict)
    d_eth2rowcol = {}

    with open('etnicity_felix.json') as f:
        js = json.load(f)

    d_tmp = {}
    ## 1) Calculate polygon intersection areas.
    for i, feature in enumerate(js['features']):

        polygon = shapely.geometry.shape(feature['geometry'])
        FAMILY = feature['properties']['FAMILY']
        ETHNICITY = feature['properties']['ETHNICITY']
        ID = feature['properties']['ID']
        if args.verbose:
            print(i, ID, FAMILY, ETHNICITY)
        if ID == 0 and FAMILY == '' and ETHNICITY == '':
            continue
        ## Do not color/build Madagascar.
        if FAMILY == 'Malagasy':
            continue

        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        ## row is miscalculated...
        ## row = args.n-row
        row_min = args.n*(min_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        row_max = args.n*(max_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        col_min = args.n*(min_lon-xllcorner)/_den-0.5
        col_max = args.n*(max_lon-xllcorner)/_den-0.5

        ## First correct FAMILY of features to allow correct grouping.
        FAMILY, ETHNICITY = fix_family_ethnicity(
            FAMILY, ETHNICITY, ID, max_lat, min_lon)
        ## Do not color multiple separate polygons with the same color,
        ## if their ethnicity is not identical.
        if ETHNICITY == '':
            ETHNICITY = ID

        ## Do floor and ceiling to avoid non-overlapping polygons.
        ## loop from South to North
        assert math.floor(row_min) < math.ceil(row_max)
        for row in range(math.floor(row_min-1), math.ceil(row_max+1)):
            latitude = factor*(row+0.5)+yllcorner
            latitude -= cellsize*(ncols-nrows)/2
            ## loop from West to East
            assert math.floor(col_min-1) < math.ceil(col_max)
            for col in range(math.floor(col_min), math.ceil(col_max+1)):
                longitude = factor*(col+0.5)+xllcorner

                area = calculate_intersection_area(
                    a_2D_areas, longitude, latitude, factor, polygon)

                if area == 0:
                    continue

                if not a_2D_areas[row][col]:
                    a_2D_areas[row][col] = {}
                a_2D_areas[row][col][(FAMILY, ETHNICITY)] = area

    ## 2) Calculate largest intersection areas.
    for row in range(args.n):
        for col in range(args.n):
            ## No polygons intersect.
            if not a_2D_areas[row][col]:
                continue
            ## Avoid grey colors (Misc/Other) if possible.
            for materialID_grey in (194, 199):
                if all([
                    len(a_2D_areas[row][col].keys()) > 1,
                    materialID_grey in a_2D_areas[row][col].keys(),
                    ]):
                    del a_2D_areas[row][col][materialID_grey]
            ## Color by polygon intersection with max area.
            FAMILY, ETHNICITY = max(
                a_2D_areas[row][col], key=a_2D_areas[row][col].get)
            try:
                d_eth2rowcol[(FAMILY, ETHNICITY)].append((row, col))
            except KeyError:
                d_eth2rowcol[(FAMILY, ETHNICITY)] = [(row, col)]

    ## 3) Assign colors.
    for FAMILY, ETHNICITY in sorted(
        d_eth2rowcol, key=lambda k: len(d_eth2rowcol[k]), reverse=True):
        ## Checking heights and extreme points to see if special colors,
        ## for which 1x1 plates are not available, can be used.
        ## Heights within polygon.
        heights = set()
        heights_vicinal = set()
        ##
        extremes = [args.n, 0, args.n, 0]
        color_family = d_color2family[d_family2color[FAMILY]]
##        colors = d_colorfamily2color[color_family]
##        c = collections.Counter(colors)
        c = collections.Counter()
        for row, col in d_eth2rowcol[(FAMILY, ETHNICITY)]:

            ## Heights.
            h = a_2D_density[args.n-row][col]
            heights.add(h)

##            if h == 0:
##                print('zero height', FAMILY, ETHNICITY, args.n-row, col)

            ## Extremes.
            if args.n-row < extremes[0]:
                extremes[0] = args.n-row
            if args.n-row > extremes[1]:
                extremes[1] = args.n-row
            if col < extremes[2]:
                extremes[2] = col
            if col > extremes[3]:
                extremes[3] = col

            ## Count nearby colours.
            dist = 1
            for i in range(-dist, dist+1):
                for j in range(-dist, dist):
                    ## Do not consider diagonally vicinal grid points.
                    heights_vicinal.add(a_2D_density[args.n-row+i][col+j])
                    if i != 0 and j != 0:
                        continue
                    mID = a_2D_mIDs[args.n-row+i][col+j]
                    c[mID] += 1

        ## a) 1x1 plate not available, but 1x2 plate and 1x2 brick available.
        if all([
            any([
                color_family == 'blue' and c[321] == 0,  # dark azure
##                color_family == 'purple' and c[268] == 0,  # medium lilac (dark purple)
##                color_family == 'purple' and c[221] == 0,  # bright purple (dark pink)
                ]),
            ## 1x2 plates available.
            ## 2x2 area can be filled in both directions,
            ## but 1x2 area can only be filled in one direction.
            all([
                ## Rectangular area.
                (
                    (extremes[1]-extremes[0]+1)*
                    (extremes[3]-extremes[2]+1)
                    ) == len(d_eth2rowcol[(FAMILY, ETHNICITY)]),
                ## Length 2 in both directions.
                (extremes[1]-extremes[0]+1) % 2 == 0,
                (extremes[3]-extremes[2]+1) % 2 == 0,
                ## All at same height.
                len(heights) == 1,
                ]),
            ]):
            if color_family == 'blue' and c[321] == 0:
                color = 321  # dark azure
        ## b) 1x1 plate not available, but 1x1 brick available.
        elif all([
            any([
                color_family == 'blue' and c[322] == 0,  # medium azure
                color_family == 'blue' and c[323] == 0,  # aqua (unikitty blue)
                color_family == 'purple' and c[268] == 0,  # medium lilac (dark purple)
                color_family == 'purple' and c[221] == 0,  # bright purple (dark pink)
                color_family == 'green' and c[151] == 0,  # sand green
                ]),
            ## At least 1 brick height.
            all([h >= 3 for h in heights]),
            any([
                ## Not at the edge, so pieces below can be plates.
                ## And vicinal must be able to cover all of a brick.
                all([min(heights_vicinal) > (h % 3) for h in heights]),
                ## Whole stack can be replaced with bricks.
                all([h % 3 == 0 for h in heights]),
                ]),
            any([
                ## 1x1 area.
                all([
                    extremes[1]-extremes[0]+1 == 1,
                    extremes[3]-extremes[2]+1 == 1,
                    ]),
                ## Necessary to avoid identical colors next to each other.
                len(set(d_colorfamily2color[color_family])-set(c)) == 0,
                ]),
            ]):
            if color_family == 'blue' and c[322] == 0:  # medium azure
                color = 322
            elif color_family == 'blue' and c[323] == 0:
                color = 323
            elif color_family == 'purple' and c[268] == 0:  # medium lilac (dark purple)
                color = 268
            elif color_family == 'purple' and c[221] == 0:  # bright purple (dark pink)
                color = 221
            elif color_family == 'green' and c[151] == 0:  # sand green
                color = 151
        ## c) 1x1 plate available.
        else:
            ## Color as least common and preferably use cheap colors.
            for color in d_colorfamily2color[color_family]:
                if c[color] == 0:
                    break
            else:
                color = c.most_common()[-1][0]
                print(
                    'copycolor', color, color_family,
                    FAMILY, ETHNICITY, len(d_eth2rowcol[(FAMILY, ETHNICITY)]), c)

        if color in (194, 199):
            print(
                'grey', color, ID, FAMILY, ETHNICITY,
                len(d_eth2rowcol[(FAMILY, ETHNICITY)]),
                min_lat, min_lon, max_lat, max_lon)

        for row, col in d_eth2rowcol[(FAMILY, ETHNICITY)]:
            a_2D_mIDs[args.n-row][col] = color

        ## Continue loop over family, ethnicity tuples.
        continue

    return a_2D_mIDs


def calculate_intersection_area(
    a_2D_areas, longitude, latitude, factor, polygon):

    polygon_1x1 = shapely.geometry.Polygon([
        (longitude-.5*factor, latitude-.5*factor),
        (longitude+.5*factor, latitude-.5*factor),
        (longitude+.5*factor, latitude+.5*factor),
        (longitude-.5*factor, latitude+.5*factor),
        (longitude-.5*factor, latitude-.5*factor),
        ])
    ## Calculate the area
    ## instead of solving the point-in-polygon problem.
    area = polygon.intersection(
        polygon_1x1).area/(factor*factor)

    return area


def assign_color_families_to_language_families(args):

    d_family2color = {}
    with open(args.colors) as f:
        for line in f:
            if line == '\n':
                continue
            if line[0] == '#':
                continue
            l = line.split('\t')
            d_family2color[l[0]] = int(l[1])
    d_family2color['Semitic'] = d_family2color['Semitic: Arab, Bedouin']
    ## Define super classes.
    d_family2color['Niger-Congo'] = d_family2color['West Atlantic']
    d_family2color['Nilo-Saharan'] = d_family2color['Nilotic']
    ## Make Afrikaans black or dark grey.
    assert 26 not in d_family2color.values()
    assert 194 not in d_family2color.values()
    assert 199 not in d_family2color.values()
    d_family2color['Afrikaans'] = 26  # black
    d_family2color['Other'] = 194  # medium stone grey
    d_family2color['Miscellaneous / Unclassified'] = 199  # dark stone grey

    return d_family2color


def fix_family_ethnicity(family, ethnicity, ID, max_lat, min_lon):

    ## Bantu languages in Angola incorrectly assigned to the Kru family.
    if family == 'Kru' and max_lat < 0:
        family = 'Bantu'
    ## Koman languages incorrectly assigned to the Maban family.
    ## Both Nilo-Saharan though.
    elif family == 'Maban' and min_lon > 30:
        family = 'Nilotic'
        ## http://en.wikipedia.org/wiki/Koman_languages
        ## http://en.wikipedia.org/wiki/Berta_languages
        ## http://en.wikipedia.org/wiki/Eastern_Jebel_languages
        assert ethnicity in (
            ## Koman
            'OPUUO',  # Ethiopia
            'GULE',  # Sudan
            'KOMA',  # Sudan
            ## Chari-Nile > Berta
            'BERTA',
            ## Chari-Nile > Eastern Sudanic
            'GAAM',
            'MUN',
            'MASONGO',
            'BAREA',
            ## Chari-Nile
            'KUNAMA',
            )
    ## http://en.wikipedia.org/wiki/Sandawe_language
    ## http://en.wikipedia.org/wiki/Languages_of_Namibia
    ## http://en.wikipedia.org/wiki/San_languages
    ## http://en.wikipedia.org/wiki/Khoekhoe_language#Dialects
    elif family == 'Sandawe' and max_lat < -10:
        assert ID in (148, 149, 186)
        ## http://en.wikipedia.org/wiki/%C7%82Aakhoe_dialect
        if ID == 186:
            family = 'Khoi: Nama, Bergdama'
            ethnicity = 'Haiǁom'
        ## http://en.wikipedia.org/wiki/Damara_people
        elif ID == 148:
            family = 'San'  # color red
            ethnicity = 'Nama-Damara'
        ## http://en.wikipedia.org/wiki/Nama_people
        else:
            assert ID == 149
            family = 'San'  # color red
            ethnicity = 'Nama-Damara'
    ## http://en.wikipedia.org/wiki/Kwadi_language
    elif ethnicity == 'BERGDAMA' and ID == 147:
        assert family == 'Khoi: Nama, Bergdama'
        ethnicity = 'Kwadi'
    ## https://en.wikipedia.org/wiki/Bajuni_dialect
    elif family == 'Other' and ethnicity == 'TIKUU':
        family = 'Bantu'
    ## http://en.wikipedia.org/wiki/Awngi_language
    elif family == 'Fufulde' and ethnicity == 'AWIYA':
        family = 'Cushitic'
    ## http://en.wikipedia.org/wiki/Dongo_language_(Nilo-Saharan)
    elif family == 'Fufulde' and ethnicity == 'DONGO':
        family = 'Nilo-Saharan'
    ## http://en.wikipedia.org/wiki/Ariaal_people
    elif family == 'Nilotic / Bantoid' and ethnicity == 'ARIAAL':
        family = 'Nilotic'
    ## http://en.wikipedia.org/wiki/Chadian_Arabic
    elif family == 'Saharan / Cushitic' and ethnicity == 'SHUWA':
        family = 'Semitic'
    ## https://www.ethnologue.com/language/pnz
    ## http://en.wikipedia.org/wiki/Pana_language
    elif family == 'Fufulde / Adamawa-Ubangia' and ethnicity == 'PANI':
        family = 'Niger-Congo'
    ## http://en.wikipedia.org/wiki/Yedina_language
    elif family == 'Chadic / Cushitic' and ethnicity == 'BUDUMA':
        family = 'Chadic'
    ## http://en.wikipedia.org/wiki/Fongoro_language
    elif family == 'Saharan / Nilotic' and ethnicity == 'FONGORO':
        family = 'Nilo-Saharan'
    ## http://en.wikipedia.org/wiki/Aja_language_(Nilo-Saharan)
    elif family == 'Adamawa-Ubangian / Chari-Nile' and ethnicity == 'AJA':
        family = 'Chari-Nile'
    ## http://en.wikipedia.org/wiki/Shatt_language
    elif family == 'Chari-Nile / Nilotic' and ethnicity == 'SHATT':
        family = 'Nilo-Saharan'
    ## Afrikaans
    elif family == 'Miscellaneous / Unclassified':
        if ethnicity == 'YEKE':
            family = 'Bantu'
        elif ethnicity in (
            ## http://en.wikipedia.org/wiki/Kambari_languages
            'KAMBARI',
            ## http://en.wikipedia.org/wiki/Kudu-Camo_language
            'KUDU',
            ):
            family = 'Niger-Congo'
        elif max_lat < -20:
            family = 'Afrikaans'
            assert ID in (
                153,  # Johannesburg
                138,  # Bloemfontein
                141,  # Durban
                143,  # East London
                146,  # Cape Town
                )
        elif ethnicity == '' and ID in (
            425,  # Liberia coast
            432,  # Sierra Leone coast
            804,  # Abuja, Nigeria
            1261, 1463,  # Nigeria/Cameroon border
            1496, 1502, 1509, 1616, 1736, 1738, 1769,  # Nigeria
            ):
            pass
        else:
            assert ethnicity == 'DARI'  # Cameroon
    elif family == 'Other':
        ## http://en.wikipedia.org/wiki/Daisu_language
        if ethnicity == 'DAISU':
            family = 'Bantu'
    elif '/' in family and ethnicity != '':
        try:
            colors = [d_family2color[s.strip()] for s in family.split('/')]
        except:
            colors = []
        if len(colors) == 1:
            print(family, ethnicity, colors)
            stop
    else:
        pass

    return family, ethnicity


def numpy2lxfml(
    args,
    a_2D_density, lxfml, a_2D_mIDs, a_2D_buried,
    a_3D_dIDs, a_3D_mIDs, a_3D_angle,
    ):

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))

    plate_size = args.plate_size

    ## initiate refID count
    refID = 0
    ## open file
    with open(lxfml, 'w') as f:
        f.write(head('Africa'))
##        f.write('<<BuildingInstructions>\n')
##        f.write('<BuildingInstruction>\n')
        ## loop from North to South
        for row1 in range(int(args.n/plate_size)):
            ## loop from West to East
            for col1 in range(int(args.n/plate_size)):
                if args.verbose:
                    print('numpy2lxfml, row', row1, col1)
                with open(
                    '{}.{}.{}.lxfml'.format(
                        lxfml[:-len('.lxfml')], row1, col1), 'w') as lxfml_plate:
                    lxfml_plate.write(head('Africa.{}.{}'.format(row1, col1)))
##                    lxfml_plate.write('<BuildingInstructions>\n')
##                    lxfml_plate.write('<BuildingInstruction>\n')
                    ## Add grey 48x48 base plate.
                    assert plate_size == 48
                    line = format_line(
                        refID, designID_baseplate, materialID_grey,
                        (row1+1)*plate_size-1, 0, (col1+1)*plate_size-1,
                        a_3D_angle)
                    f.write(line)
                    lxfml_plate.write(line)
                    ## Write first step.
                    steps = '      <Step>\n'
                    steps += '        <PartRef partRefID="{}"/>\n'.format(refID)
                    steps += '      </Step>\n'
                    ## Increment refID.
                    refID += 1
                    ## Loop over y first to be able do BuildingInstruction Steps.
                    ## This is obviously slower than looping over
                    ## for y in range(int(a_2D_density[row][col])):
                    for y in range(args.layers):
                        ## Building step for layer not yet initiated.
                        bool_init_step = False
##                        lxfml_plate.write('<Step>\n')
##                        f.write('<Step>\n')
                        ## Loop over rows and cols of base plate.
                        for row2 in range(plate_size):
                            row = row1*plate_size+row2
                            for col2 in range(plate_size):
                                col = col1*plate_size+col2
                                ## Generate line for coordinate.
                                line = generate_line(
                                    y, row, col, refID,
                                    a_2D_mIDs, a_2D_density,
                                    a_2D_buried,
                                    a_3D_dIDs, a_3D_mIDs, a_3D_angle,
                                    )
                                if line:
                                    f.write(line)
                                    lxfml_plate.write(line)
                                    if not bool_init_step:
                                        ## Initiate building step for layer.
                                        steps += '      <Step>\n'
                                        bool_init_step = True
                                    steps += '        <PartRef partRefID="{}"/>\n'.format(
                                        refID)
                                    ## Increment refID after appending step.
                                    refID += 1
                                ## continue loop over col2
                                continue
                            ## continue loop over row2
                            continue
##                        lxfml_plate.write('</Step>\n')
##                        f.write('</Step>\n')
                        ## Terminate building step if it was initiated.
                        if bool_init_step:
                            steps += '      </Step>\n'
                        ## continue loop over ty
                        continue
                    ## add tail to plate lxfml file
##                    lxfml_plate.write('</BuildingInstruction>')
##                    lxfml_plate.write('</BuildingInstructions>')
                    lxfml_plate.write(tail(steps))
                    ## close plate lxfml file
                    lxfml_plate.close()
                    pass
                ## continue loop over tz
                continue
            ## continue loop over tx
            continue
##        f.write('</BuildingInstruction>')
##        f.write('</BuildingInstructions>')
        ## add tail to lxfml file
        f.write(tail(''))
        ## close file
        f.close()
        pass

    return


def generate_line(
    y, row, col, refID,
    a_2D_mIDs, a_2D_density, a_2D_buried,
    a_3D_dIDs, a_3D_mIDs, a_3D_angle):

    if y >= a_2D_density[row][col]:
        return None

    ## Skip ocean and land outside Africa.
    if a_2D_mIDs[row][col] == 0:
        return None

    ## Replacement brick/plate.
    if a_3D_mIDs[y][row][col] != 0:
        materialID = a_3D_mIDs[y][row][col]
    ## 1x1 plate, white/buried
    elif y < a_2D_buried[row][col]:
        materialID = materialID_buried
    ## 1x1 plate, exposed/colored
    else:
        materialID = a_2D_mIDs[row][col]

    ## 1x1 plate
    if a_3D_dIDs[y][row][col] == 0:
        designID = 3024  # 1x1 plate
    ## larger plate
    else:
        ## empty space for larger plate
        if a_3D_dIDs[y][row][col] == -1:
            return None
        ## larger plate
        else:
            designID = a_3D_dIDs[y][row][col]
            pass
        pass

    line = format_line(
        refID, designID, materialID, row, y, col, a_3D_angle)

    return line


def format_line(refID, designID, materialID, row, y, col, a_3D_angle):

    if materialID == materialID_buried:
        materialID = 1

    ## LEGO dimensions
    ## http://upload.wikimedia.org/wikipedia/commons/1/1a/Lego_dimensions.svg
    P = 0.8  # cm
    h = 0.32  # cm

    line = '        <Part'
    line += ' refID="{0}"'.format(refID)
    line += ' designID="{0}"'.format(designID)
    line += ' materialID="{0}"'.format(materialID)
    line += ' itemNos="{0}{1}"'.format(designID, materialID)

    ## rotate / rotation
    if y % 2 == 1 and designID != designID_baseplate:
        angle = 90
        dtx = P * (d_dIDs2dim[designID]['dim0']-1)
        dtz = 0
    else:
        angle = 0
        dtx = 0
        dtz = 0
    if a_3D_angle[y][row][col] and designID != 4186:
        angle = 90 - angle
    line += ' angle="{}" ax="0" ay="1" az="0"'.format(angle)
    line += ' tx="{0:.1f}"'.format(row*-P + dtx)
    line += ' ty="{0:.2f}"'.format(y*h)
    line += ' tz="{0:.1f}"'.format(col*P + dtz)
    line += '/>\n'

    return line


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
   <Camera refID="0" fieldOfView="80" distance="0" angle="0" ax="0" ay="0" az="0" tx="0" ty="0" tz="0"/>
  </Cameras>
  <Scene cameraRefID="0">
    <Model>
      <Group refID="0" angle="0" ax="0" ay="1" az="0" tx="0" ty="0" tz="0">
'''.format(name)

    return s


def tail(steps):

#    ty += 0.32 ## per plate in height
    s = '''\
      </Group>
    </Model>
  </Scene>
  <BuildingInstructions>
    <BuildingInstruction>
{}
    </BuildingInstruction>
  </BuildingInstructions>
</LXFML>
'''.format(steps.rstrip())

    return s


def read_gpw(file_gpw):

    with open(file_gpw) as f:
        ncols = int(f.readline().strip().split()[1])
        nrows = int(f.readline().strip().split()[1])
        assert ncols % 2 == 0 and nrows % 2 == 0
        assert ncols > nrows
        assert ncols == 2160  # -30+2160*2.5/60 = 60 = 60E
        assert nrows == 1920  # -40+1920*2.5/60 = 40 = 40N
        xllcorner = int(f.readline().strip().split()[1])
        yllcorner = int(f.readline().strip().split()[1])
        assert xllcorner == -30  # -30 = 30E
        assert yllcorner == -40  # -40 = 40S
        cellsize = float(f.readline().strip().split()[1])
        assert round(cellsize, 8) == round(2.5/60, 8)
        for i in range(1):
            f.readline()

        a_gpw = np.zeros((max(nrows, ncols), max(nrows, ncols)))
        max_gpw = 0
        x1 = (max(nrows, ncols)-nrows)/2
        x2 = (max(nrows, ncols)+nrows)/2
        for row in range(max(nrows, ncols)):
##            print(row)
##            latitude = yllcorner+(nrows-row+x1-.5)*2.5/60
            if row >= x1 and row < x2:
                cols = f.readline().rstrip().split()
                pass
            else:
                cols = ['0']*ncols

            for col, s in enumerate(cols):
##                longitude = xllcorner+(col+.5)*2.5/60

                if s == '0':
                    a_gpw[row][col] = 0
                elif s == '-9999':
                    a_gpw[row][col] = 0
                else:
                    a_gpw[row][col] = float(s)

    return a_gpw


if __name__ == '__main__':
    main()
