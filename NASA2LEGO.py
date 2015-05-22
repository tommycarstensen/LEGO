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
import itertools
import argparse


## 1xlen plates
d_len2designIDs = {1: 3024, 2: 3023, 3: 3623, 4: 3710, 6: 3666, 8: 3460}
d_designIDs2dim = {v:{'dim0':1,'dim1':k} for k,v in d_len2designIDs.items()}
## square plates
d_designIDs = {16: 91405, 8: 41539, 6: 3958, 4: 3031}
for k,v in d_designIDs.items():
    d_designIDs2dim[v] = {'dim0':k,'dim1':k}
designID_baseplate = 4186
materialID_grey = 194


## tmp!!! sep13


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

    affix1 = 'out_NASA2LEGO/{0}_{1:d}x{1:d}_dens{2}'.format(
        args.affix, args.n, args.density)

    if not os.path.isfile('{}.npy'.format(affix1)):
        ## slow
        array_2D_density = asc2np(args, affix1)
    else:
        ## fast
        array_2D_density = np.load('{}.npy'.format(affix1))

##    x = []
##    for row in range(240):
##        for col in range(240):
##            if array_2D_density[row][col] == 0:
##                continue
####            if array_2D_density[row][col] >= 200:
####                continue
##            x += [math.log(array_2D_density[row][col])]
####            x += [array_2D_density[row][col]]
##    print(len(x))
##    import matplotlib.pyplot as plt
##    plt.xlabel('Population Density (arbitrary unit)')
##    plt.ylabel('Frequency')
##    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
####    hist, bins = np.histogram(array_2D_density, bins=50)
####    width = 0.7 * (bins[1] - bins[0])
####    center = (bins[:-1] + bins[1:]) / 2
#####    plt.bar(center, hist, width=width)
####    plt.bar(center, hist)
##    plt.show()
##    stop

    array_2D_density = normalize(array_2D_density, args)

    ## fast
    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))
    assert nrows == 1920 and ncols == 2160
    assert yllcorner == -40 and xllcorner == -30
    assert round(yllcorner+cellsize*nrows,0) == 40
    assert round(xllcorner+cellsize*ncols,0) == 60

##    row = 0
##    col = 0
##    n = 240
##    dy = cellsize*(row+0.5)*max(nrows, ncols)/n
##    dy -= cellsize*(ncols-nrows)/2
##    latitude = yllcorner+dy
##    dx = cellsize*(col+0.5)*max(nrows, ncols)/n
##    longitude = xllcorner+dx
##    print(yllcorner,xllcorner,latitude,longitude,dy,dx,row,col)
##    print(ncols, nrows, cellsize, yllcorner, xllcorner, yllcorner+cellsize*nrows, xllcorner+cellsize*ncols)
##    stop
                
    ## slow - convert ethnicity polygons from json file to numpy array
    ## i.e. assign a tentative materialID to each 2D point
    array_2D_materialIDs = json2array(
        args, nrows, ncols, xllcorner, yllcorner, cellsize)

    ## Set zero density if no color. Not a good solution, because lakes will
    ## randomly appear and disappear, because the polygons are poorly defined.
##    array_2D_density = np.where(
##        array_2D_materialIDs!=0, array_2D_density, 0)
    color_as_nearby(args, array_2D_density, array_2D_materialIDs)

##    x = []
##    for row in range(240):
##        for col in range(240):
##            if array_2D_density[row][col] == 0:
##                continue
####            if array_2D_density[row][col] >= 200:
####                continue
####            x += [math.log(array_2D_density[row][col])]
##            x += [array_2D_density[row][col]]
##    print(len(x))
##    import matplotlib.pyplot as plt
##    plt.xlabel('Population Density (arbitrary unit)')
##    plt.ylabel('Frequency')
##    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='g', alpha=0.75)
####    hist, bins = np.histogram(array_2D_density, bins=50)
####    width = 0.7 * (bins[1] - bins[0])
####    center = (bins[:-1] + bins[1:]) / 2
#####    plt.bar(center, hist, width=width)
####    plt.bar(center, hist)
##    plt.show()
##    stop

    ## fast
    array_2D_buried = find_buried(array_2D_density, args.layers)

    ## fast
    array_3D_designIDs = find_connected_buried(args, array_2D_buried)

    ## 
    array_3D_designIDs, array_3D_materialIDs = find_connected_exposed(
        array_2D_density, array_2D_buried, args.layers,
        array_3D_designIDs, array_2D_materialIDs)

    ## slow
    lxfml = '{}_y{:d}_{}.lxfml'.format(affix1, args.layers, args.norm)
    numpy2lxfml(
        args,
        array_2D_density, lxfml, array_2D_materialIDs, array_2D_buried,
        array_3D_designIDs, array_3D_materialIDs,
        nrows, ncols, xllcorner, yllcorner, cellsize)

##    print(pcount_max)
##    print(np.amin(array_gpw))
##    from collections import Counter
##    print(Counter([float(x) for x in np.nditer(array_gpw)]))

    return


def color_as_nearby(args, array_2D_density, array_2D_materialIDs):

    assert array_2D_density.shape == array_2D_materialIDs.shape

    dist = 2
    for x in range(dist,array_2D_density.shape[0]-dist):
        for y in range(dist,array_2D_density.shape[1]-dist):
            if array_2D_density[x][y] == 0:
                continue
            if array_2D_materialIDs[x][y] > 0:
                continue
            ## If not densely populated and not covered by language polygon
            ## then delete. Because probably near coast or lake.
            ## Including this deletion does not cause problems in Egypt.
            ## todo: generate map with and without!!!
            if array_2D_density[x][y] == 1:
                if args.verbose:
                    print('near zero density', x, y)
                array_2D_density[x][y] = 0
                continue
            ## Count nearby materialIDs.
            cnt = collections.Counter(
                array_2D_materialIDs[x+dx][y+dy]
                for dx in range(-dist, dist+1) for dy in range(-dist, dist+1)
                )
            del cnt[0]
            if len(cnt) == 0:
#                array_2D_materialIDs[x][y] = 119 # tmp!!!

                continue
            ## Set materialID to most frequently occuring nearby materialID.
            array_2D_materialIDs[x][y] = cnt.most_common(1)[0][0]

    return


def argparser():
    
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--affix', default='afds00g',
        choices=['afp00ag', 'afp00g', 'afds00ag', 'afds00g'])

    parser.add_argument(
        '--plate_size', default=48, type=int, choices=[16,32,48],
        help='Size of base plates to build on.')

    parser.add_argument(
        '--plate_cnt', default=5,
        help='Number of base plates along one dimension.')

    parser.add_argument(
        '--zero', help = 'add zero values to average?', action='store_true')

    parser.add_argument(
        '--verbose', help = 'Be verbose?', action='store_true')

    parser.add_argument(
        '--density', help = 'Take max or mean of density grid?', choices=['max','mean'], default='max')

    parser.add_argument('--layers', default = 27, type=int)

    parser.add_argument('--norm', default = 'log', choices=[
        'log10', 'log2', 'unity', 'log']) # unity is feature scaling

    parser.add_argument('--colors', required=True)

    args = parser.parse_args()

    # script fast if multiple of 2160
    args.n = n = nrows = ncols = args.plate_cnt*args.plate_size

    assert args.layers % 3 == 0

    return args


def asc2np(args, affix1):
    
    array_gpw = read_gpw('{}.asc'.format(args.affix))
    nrows2, ncols2 = np.shape(array_gpw)

    assert np.shape(array_gpw)[0] > args.n

    den = args.n/fractions.gcd(ncols2, args.n)
    num = int(ncols2/fractions.gcd(ncols2, args.n))

    array_2D_density = np.zeros((args.n, args.n))
    array_2D_density_cnt = np.zeros((args.n, args.n))

    print('converting NASA array to LEGO array')
    for x1 in range(args.n):
        print('gpw2LEGO', x1, args.n)
        for y1 in range(args.n):
            l = []
            for x2 in range(num):
                for y2 in range(num):
                    x3 = int((x1*num+x2)/den)
                    y3 = int((y1*num+y2)/den)
##                    array_2D_density[x1][y1] += array_gpw[x3][y3]
##                    array_2D_density[x1][y1] += array_gpw[x3][y3]
                    l.append(array_gpw[x3][y3])
                    if array_gpw[x3][y3] > 0:
                        array_2D_density_cnt[x1][y1] += 1
            if l.count(0) >= 0 and sum(l) > 0:
                if y1 == 65 and x1 > 95:
                    print(x1, y1, l.count(0), 'mean', sum(l)/len(l), 'max', max(l), sum(l), sum(sorted(l)[1:-1]), l)
            array_2D_density[x1][y1] = max(l)

    np.save(affix1, array_2D_density)

    del array_gpw
    del array_2D_density_cnt

    return array_2D_density


def find_connected_exposed(
    array_2D_density, array_2D_buried, layers,
    array_3D_designIDs, array_2D_materialIDs):

    print('find connected exposed plates and remaining buried plates')

    n1 = n2 = n = np.shape(array_2D_density)[0]

    array_3D_materialIDs = np.zeros(np.shape(array_3D_designIDs), int)

    max_len = max(d_len2designIDs.keys())

    for layer in range(layers):
        ## build plates horizontally and vertically in each layer
        if layer % 2 == 0:
            irow = 1
            jrow = 0
            icol = 0
            jcol = 1
            continue ## tmp!!! tmp1!!!
        else:
            irow = 0
            jrow = 1
            icol = 1
            jcol = 0
            continue ## tmp!!! tmp2!!! sep13
        for i in range(n):
            seq = []
##            d_colors = {}
            for j in range(n):
                row = i*irow+j*jrow
                col = i*icol+j*jcol
                ## ocean or lake (or desert)
                if array_2D_density[row][col] == 0:
                    seq += ['x']
                ## Position is above structure.
                elif layer >= array_2D_density[row][col]:
                    seq += ['x']
                ## skip filled
                elif array_3D_designIDs[layer][row][col] != 0:
                    seq += ['x']
                ## buried
                elif layer < array_2D_buried[row][col]:
                    seq += ['0']
                ## Not inside Felix2001 GeoJSON polygon
                ## e.g. Spain and Saudi Arabia
                elif array_2D_materialIDs[row][col] == 0:  # tmp!!!
                    seq += ['x']
                    continue  # tmp!!!
                else:
                    seq += [int(array_2D_materialIDs[row][col])]
##                    materialID = array_2D_materialIDs[row][col]
##                    try:
##                        color = d_colors[materialID]
##                    except KeyError:
##                        color = len(d_colors.keys())+1
##                        d_colors[materialID] = color
##                    seq += [color]
                ## Continue loop over j.
                continue
            ## no bricks along line
            if seq == n*['x']:
                continue
            seq = find_consecutive(seq, n)
            append_designID_materialID_main(
                seq, layer, i, irow, jrow, icol, jcol, d_len2designIDs, max_len,
                array_3D_designIDs, array_3D_materialIDs)

    return array_3D_designIDs, array_3D_materialIDs


def append_designID_materialID_main(
    seq, layer, i, irow, jrow, icol, jcol, d_len2designIDs, max_len,
    array_3D_designIDs, array_3D_materialIDs):

    pos = 0
    for materialID, g in itertools.groupby(seq):
        ## Get length of iterator.
        len_group = len(list(g))
##        print('pos', pos, 'len_group', len_group)
        if materialID == 'x':
            pos += len_group
##            print('a', pos, materialID)
            continue
        if materialID == 0:
            print(seq)
            print(materialID)
            print(type(materialID))
            stop
        ## Look up designID of given length.
        try:
            length = len_group
            designID = d_len2designIDs[length]
            if layer == 25 and length == 4 and materialID == 194:
                print('poslen1', pos, length)
            pos = append_designID_materialID(
                layer, designID, materialID,
                array_3D_designIDs, array_3D_materialIDs,
                pos, i, irow, jrow, icol, jcol, length)
            if layer == 25 and length == 4 and materialID == 194:
                print('poslen2', pos, length)
            if layer == 25 and length == 4 and materialID == 194:
                print('b', pos, materialID)
##            print('b', pos, materialID)
        except KeyError:
            ## How many plates of max length will fit?
            for k in range(len_group//max_len):
                length = max_len
                designID = d_len2designIDs[length]
                pos = append_designID_materialID(
                    layer, designID, materialID,
                    array_3D_designIDs, array_3D_materialIDs,
                    pos, i, irow, jrow, icol, jcol, length)
            if layer == 25 and length == 4 and materialID == 194:
                print('c', pos, len_group, len_group//max_len, materialID)
##            print('c', pos, len_group, len_group//max_len, materialID)
            ## How much space left after filling with plates of max length?
            mod = len_group % max_len
            ## No space left.
            if mod == 0:
                continue
            ## Plate length exists.
            try:
                length = mod
                designID = d_len2designIDs[length]
                pos = append_designID_materialID(
                    layer, designID, materialID,
                    array_3D_designIDs, array_3D_materialIDs,
                    pos, i, irow, jrow, icol, jcol, length)
                if layer == 25 and length == 4 and materialID == 194:
                    print('d', pos, materialID)
##                print('d', pos, materialID)
            ## Plate length does not exist.
            except KeyError:
                assert max_len == 8
                if mod == 7:
                    len1 = 4
                    len2 = 3
                elif mod == 5:
                    len1 = 3
                    len2 = 2
                else:
                    print('mod', mod)
                    print(k, list(g), len_group)
                    stop
                for length in (len1, len2):
                    designID = d_len2designIDs[length]
                    pos = append_designID_materialID(
                        layer, designID, materialID,
                        array_3D_designIDs, array_3D_materialIDs,
                        pos, i, irow, jrow, icol, jcol, length)
                    if layer == 25 and length == 4 and materialID == 194:
                        print('e', pos, materialID)
##                    print('e', pos, materialID)
        ## Continue loop over materialID groups.
        continue

    return


def append_designID_materialID(
    layer, designID, materialID,
    array_3D_designIDs, array_3D_materialIDs,
    pos, i, irow, jrow, icol, jcol, length,):

    for j in range(length):
        row = i*irow+(pos+j)*jrow
        col = i*icol+(pos+j)*jcol
        ## replace Southern plate
        ## replace Eastern plate
        if layer % 2 == 1 and j == 0:  ## sep13 tmp!!!
            array_3D_designIDs[layer][row][col] = designID
            array_3D_materialIDs[layer][row][col] = materialID
        elif j == length-1:
            array_3D_designIDs[layer][row][col] = designID
            array_3D_materialIDs[layer][row][col] = materialID
        else:
            array_3D_designIDs[layer][row][col] = -1
            array_3D_materialIDs[layer][row][col] = 0

    pos += length

    return pos


def find_consecutive(seq, repeat):

    gap = 'x'
    buried = '0'

    groups = [list(g) for k, g in itertools.groupby(seq)]
    seq2 = []
    for i, g in enumerate(groups):
##        print('*',i,g)
        if g[0] == gap:
            seq2 += g
        elif g[0] != buried:
            seq2 += g
        ## elif g[0] == '0'
        else:
            ## first group
            if i == 0:
                try:
                    ## next group is gap
                    if groups[i+1][0] == gap:
                        seq2 += g
                    ## color first group same as next group
                    else:
                        seq2 += len(g)*groups[i+1][0]
                except:
                    print(groups)
                    stop
            ## last group
            elif i+1 == len(groups):
                ## previous group is gap
                if groups[i-1][0] == 'x':
                    seq2 += g
                ## color last group same as previous group
                else:
                    seq2 += len(g)*groups[i-1][0]
            ## if end of sequence then continue previous color
            elif len(seq2)+len(g) == repeat:
                seq2 += len(g)*groups[i-1][0]
            ## gap before and after
            elif groups[i-1][0] == gap and groups[i+1][0] == gap:
                seq2 += g
            ## gap before but not after
            elif groups[i-1][0] == gap and groups[i+1][0] != gap:
                seq2 += len(g)*[groups[i+1][0]]
            ## gap after but not before
            elif groups[i-1][0] != gap and groups[i+1][0] == gap:
                seq2 += len(g)*[groups[i-1][0]]
            ## same color before and after
            elif groups[i-1][0] == groups[i+1][0]:
                seq2 += len(g)*[groups[i-1][0]]
            else:
                if len(g) >= 2:
                    len_prev = len(list(itertools.takewhile(
                        lambda x: x == groups[i-1][0], reversed(seq2))))
                    if len_prev+len(g)//2 in (1, 2, 3, 4, 6, 8):  # PAB lengths
                        ## append with color from both sides
                        seq2 += (len(g)//2)*[groups[i-1][0]]
                        seq2 += (len(g)-len(g)//2)*[groups[i+1][0]]
                        continue
                    else:
                        seq2 += (len(g)//2+1)*[groups[i-1][0]]
                        seq2 += (len(g)-len(g)//2-1)*[groups[i+1][0]]
                        continue
                elif len(g) == 1:
                    ## append to previous sequence of length 1
                    len_prev = len(list(itertools.takewhile(
                        lambda x: x == groups[i-1][0], reversed(seq2))))
                    if len_prev == 1:
                        seq2 += 1*[groups[i-1][0]]
                        continue
                    elif all([
                        ## before has length greater than 1
                        len_prev > 1,
                        ## after has length equal to 1
                        len(list(itertools.takewhile(
                            lambda x: x == groups[i+1][0], seq[len(seq2)+1:]))) == 1,
                        ]):
                        seq2 += 1*[groups[i+1][0]]
                        continue
                    elif len_prev+len(g) in (2, 3, 4, 6, 8):  # PAB lengths
                        seq2 += len(g)*[groups[i-1][0]]
                        continue
                    else:
                        seq2 += len(g)*[groups[i+1][0]]
                        continue
##    print(tuple(seq2), t)
    assert len(seq2) == repeat

    return seq2


def find_connected_buried(args, array_2D_buried):

    print('finding connected 1x1 plates and replacing them with larger plates')

    n = args.n
    layers = args.layers
    plate_size = args.plate_size

    ## DO A WHILE LOOP UNTIL area_max IS EMPTY!!!

    ## DO NOT CENTER BRICK IN A DIMENSION
    ## IF FACING OTHER PLATE IN THAT DIRECTION - ie IF AT EDGE
    ## 2015may20 - WTF did I mean by that, when I wrote it???

    ## 3D array with design IDs at each 3D coordinate
    array_3D_designIDs = np.zeros((layers, n, n), int)

    ## loop baseplates from North to South
    for row1 in range(int(n/plate_size)):
        ## loop baseplates from West to East
        for col1 in range(int(n/plate_size)):
            if args.verbose:
                print('find_connected_buried', row1, col1)
            area_max = {layer: {'area': 0} for layer in range(layers)}
            i = -1
            ## Loop while there are areas to be filled.
            while area_max:
                i += 1
                h = np.zeros((layers, plate_size, plate_size), int)  # heights
                w = np.zeros((layers, plate_size, plate_size), int)  # widths
                ## North to South
                for row2 in range(plate_size):
                    row3 = row1*plate_size+row2
                    ## West to East
                    for col2 in range(plate_size):
                        col3 = col1*plate_size+col2
                        ## loop over layers after row and col
                        ## to avoid looping over empty rows and cols in a layer
                        for layer in range(
                            array_2D_buried[row3][col3]):
                            ## skip if layer already filled
                            if not layer in area_max.keys():
                                continue
                            ## skip if filled or blank already (not np.zeros)
                            if array_3D_designIDs[layer][row3][col3] != 0:
                                continue
                            ## first row
                            if row2 == 0:
                                h[layer][row2][col2] = 1
                            ## append to previous row
                            else:
                                h[layer][row2][col2] = h[layer][row2-1][col2]+1
                            ## first col
                            if col2 == 0:
                                w[layer][row2][col2] = 1
                            ## append to previous col
                            else:
                                w[layer][row2][col2] = w[layer][row2][col2-1]+1
                            areas = [(0, None, None)]
                            min_w = w[layer][row2][col2]
                            for dh in range(h[layer][row2][col2]):
                                h1 = dh+1
                                w1 = w[layer][row2-dh][col2]
                                if w1 == 0:
                                    stop
                                    break
                                if w1 < min_w:
                                    min_w = w1
                                ## don't append area if 4x4 brick can't fit inside it
                                if h1 < 4 or min_w < 4:
                                    continue
                                areas.append((h1*min_w, h1, min_w))
                            min_h = h[layer][row2][col2]
                            for dw in range(w[layer][row2][col2]):
                                w1 = dw+1
                                h1 = h[layer][row2][col2-dw]
                                if h1 == 0:
                                    stop
                                    break
                                if h1 < min_h:
                                    min_h = h1
                                ## don't append area if 4x4 brick can't fit inside it
                                if min_h < 4 or w1 < 4:
                                    continue
                                areas.append((min_h*w1, min_h, w1))
                            area, row_span, col_span = sorted(areas)[-1]
                            if area > area_max[layer]['area']:
                                area_max[layer] = {
                                    'area': area,
                                    'row_pos': row3, 'col_pos': col3,
                                    'row_span': row_span, 'col_span': col_span}
                            ## layer loop
                            continue
                        ## col2 loop
                        continue
                    ## row2 loop
                    continue

                for layer in range(layers):
                    if not layer in area_max.keys():
                        continue
                    if area_max[layer]['area'] == 0:
                        del area_max[layer]
                        continue
                    ## get row and col span
                    row_span = area_max[layer]['row_span']
                    col_span = area_max[layer]['col_span']
                    ## calculate size of square plates to fill the area
                    for size in reversed(sorted(d_designIDs.keys())):
                        if size <= row_span and size <= col_span:
                            break
                    ## plate does not fit
                    if size > row_span or size > col_span:
                        del area_max[layer]
                        continue
                    ## get row and col pos and center the square pieces
                    row_pos = area_max[layer]['row_pos']
                    col_pos = area_max[layer]['col_pos']
    ##                ## center the square pieces (todo: see if fewer pieces if cornered)
    ##                ## DONT CENTER IF MOBILE 48x48 ELEMENTS!!!
    ##                row_pos += int(area_max[layer]['row_span']%size)
    ##                col_pos += int(area_max[layer]['col_span']%size)
                    ## calculate number of square plates fitting inside area
                    plate_rows = row_span//size
                    plate_cols = col_span//size
                    ## append plate and empty spaces to 3D array
                    for row_major in range(plate_rows):
                        for col_major in range(plate_cols):
                            ## empty spaces
                            for row_minor in range(size):
                                ## CHECK!!!
                                row = row_pos-row_major*size-row_minor
                                for col_minor in range(size):
                                    col = col_pos-col_major*size-col_minor
                                    ## make empty space for larger plate
                                    array_3D_designIDs[layer][row][col] = -1
                                    continue
                                continue
                            ## plate starts SE and extends to NW
                            ## i.e. from high row to low row
                            ## and from low col to high col
                            row = row_pos-row_major*size
                            col = col_pos-col_major*size
                            ## insert large plate
                            array_3D_designIDs[layer][row][col] = d_designIDs[size]
                            ## Continue loop over col_major
                            continue
                        ## Continue loop over row_major.
                        continue
                    area_max[layer]['area'] = 0
                    ## layer loop
                    continue
                ## while loop
                continue
            ## col1 loop
            continue
        ## row1 loop
        continue

    return array_3D_designIDs


def find_buried(array_2D_density, layers):

    n = np.shape(array_2D_density)[0]

    array_2D_buried = np.zeros(np.shape(array_2D_density), int)

    for row in range(n):
        ## nothing buried at the edge (no bricks at the edge anyway)
        if row == 0 or row == n-1:
            continue
        for col in range(n):
            ## nothing buried at the edge (no bricks at the edge anyway)
            if col == 0 or col == n-1:
                continue
            ## no buried plates if only 1 plate
            if array_2D_density[row][col] <= 1:
                continue
            ## minimum neighbouring height
            z = min([array_2D_density[row+x][col+y]
                     for x in range(-1, 2) for y in range(-1, 2)])
            if z >= 1:
##                if row//48 == 2 and col//48 == 3:
##                    if row == 143 and col == 184:
##                        for x in range(-1, 2):
##                            for y in range(-1, 2):
##                                print(x,y,array_2D_density[row+x][col+y])
##                        stop
                ## update array
                array_2D_buried[row][col] = z-1

    return array_2D_buried


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

    return ncols, nrows, xllcorner, yllcorner, cellsize


def normalize(a, args):

    n = np.shape(a)[0]

    ## https://en.wikipedia.org/wiki/Normalization_(statistics)
    if args.norm == 'log10':
        amax = np.amax(a)
        a = np.log10(np.divide(a, amax/(10**layers)))
    elif args.norm == 'log':
        amax = np.amax(a)
##        amin = np.amin(a[np.nonzero(a)])
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
                a[row][col] = max(1, args.layers*log/math.log(amax))
##                _den = (math.log(amax)-math.log(amin))
##                a[row][col] = layers*(log-math.log(amin))/_den
####        mean = sumx/cnt
####        var = sumxx/cnt-mean**2
##        print(np.percentile(a,100))
##        print(np.percentile(a,99.99))
##        print(np.percentile(a,50))
##        print(np.percentile(a,0.1))
##        print(np.percentile(a,0.01))
##        stop
####        a = np.log(np.divide(a, amax/(math.exp(layers))))
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

    return a


def json2array(args, nrows, ncols, xllcorner, yllcorner, cellsize):

    print('Converting GeoJSON polygons to array of LEGO color points')

    d_colors = {}
    with open(args.colors) as f:
        for line in f:
            if line == '\n':
                continue
            if line[0] == '#':
                continue
            l = line.split('\t')
            d_colors[l[0]] = int(l[1])

    print(d_colors)

    array_2D_materialIDs = np.zeros((args.n, args.n), int)

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
            color = 26 # black
            print(family, feature['properties']['ETHNICITY'])
        if family in ('',):
##            print(
##            feature['properties'], polygon.bounds, feature['id'],
##            feature['properties']['ID'])
##            continue
            color = 322  # medium azure
##            for key in feature.keys():
##                if key == 'geometry': continue
##                print(key, feature[key])
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        _den = cellsize * max(nrows, ncols)
        row_min = args.n*(min_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        row_max = args.n*(max_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        col_min = args.n*(min_lon-xllcorner)/_den-0.5
        col_max = args.n*(max_lon-xllcorner)/_den-0.5
        within = False
        ## Afrikaans
        if family == 'Miscellaneous / Unclassified':
            if max_lat < -20:
                family = 'Afrikaans'
                color = 140  # earth blue
                color = 5  # tan
                color = 119  # lime
            else:
                print(min_lon, min_lat, max_lon, max_lat)
        elif not family in d_colors.keys():
            print(family, min_lon, min_lat, max_lon, max_lat)
        ## loop from South to North
        for row in range(math.floor(row_min), math.ceil(row_max)):
            latitude = cellsize*(row+0.5)*max(nrows, ncols)/args.n+yllcorner
            latitude -= cellsize*(ncols-nrows)/2
            ## Bantu languages in Angola incorrectly assigned to the Kru family
            if family == 'Kru' and latitude < 0:
                family = 'Bantu'
                color = d_colors[family]
            ## loop from West to East
            for col in range(math.floor(col_min), math.ceil(col_max)):
                longitude = cellsize*(col+0.5)*max(nrows, ncols)/args.n+xllcorner
                ## incorrectly assigned to the Maban family
                if family == 'Maban' and longitude > 30:  # min_lon > 30
                    within = True
                    continue  # don't continue but re-assign!!! tmp!!!
                ## don't change from non-miscellaneous to miscellanous
                if (
                    family in ('', 'Miscellaneous / Unclassified') and
                    array_2D_materialIDs[args.n-row][col]):
                    within = True
                    continue
##                ## color with color of largest nearby stack
##                ## i.e. sort by height and filter out non-colored with if
##                if family in ('','Miscellaneous / Unclassified'):
####                if family not in d_colors.keys():
##                    for x in range(3):
##                        for y in range(3):
##                            print(x,y,array_2D_materialIDs[n-row+x-1][col+y-1])
##                            print(array_2D_density[n-row+x-1][col+y-1])
##                    print(i, family, feature['properties']['ETHNICITY'])
##                    stop
                point = shapely.geometry.Point(longitude, latitude)
                ## solve the point-in-polygen problem
                if polygon.contains(point):
                    within = True
                    array_2D_materialIDs[args.n-row][col] = color
                    pass
                ## polygon might not contain point rounded to nearest grid value
                ## hence add manually
                elif row_max-row_min < 2.5 or col_max-col_min < 2.5:
                    within = True
                    ## don't change color if not within polygon
                    if not array_2D_materialIDs[args.n-row][col] == 0:
                        continue
                    array_2D_materialIDs[args.n-row][col] = color
                    pass
                ## continue loop over cols
                continue
            ## continue loop over rows
            continue
##        if family == 'Maban' and longitude > 30:
##            print(
##        row, col, latitude, longitude, family,
##        feature['properties']['ETHNICITY'])
        if within is False:
            print(
                feature['properties']['ETHNICITY'], family, polygon.bounds,
                row_max-row_min, col_max-col_min,
                )
        elif family not in d_colors.keys() and family != 'Afrikaans':
            print(
                'color', color, 'ID', feature['id'],
                'family', family, 'ethnicity', feature['properties']['ETHNICITY'],
                'lonlat', round(longitude, 0), round(latitude, 0),
                'rowcol', row, col)
        ## continue loop over features
        continue

    return array_2D_materialIDs


def numpy2lxfml(
    args,
    array_2D_density, lxfml, array_2D_materialIDs, array_2D_buried,
    array_3D_designIDs, array_3D_materialIDs,
    nrows, ncols, xllcorner, yllcorner, cellsize):

    plate_size = args.plate_size

## todo: minimize indentation level and length of this function

    n = np.shape(array_2D_density)[0]
    a = 6378137  # semi major axis
    f = 1/298.257223563  # reciprocal flattening
    e_squared = 6.69437999014*10**-3  # first eccentricity squared
    h = 0

    ## initiate refID count
    refID = 0
    ## open file
    with open(lxfml, 'w') as f:
        f.write(head('Africa'))
##        f.write('<<BuildingInstructions>\n')
##        f.write('<BuildingInstruction>\n')
        ## loop from North to South
        for row1 in range(int(n/plate_size)):
            ## loop from West to East
            for col1 in range(int(n/plate_size)):
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
                        row1*plate_size, 0, col1*plate_size)
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
                    ## for y in range(int(array_2D_density[row][col])):
                    for y in range(args.layers):
                        ## Building step for layer not yet initiated.
                        bool_init_step = False
##                        lxfml_plate.write('<Step>\n')
##                        f.write('<Step>\n')
                        ## Loop over rows and cols of base plate.
                        for row2 in range(plate_size):
                            row = row1*plate_size+row2
                            ## calculate latitude
                            dy = cellsize*(row+0.5)*max(nrows, ncols)/n
                            dy -= cellsize*(ncols-nrows)/2
                            latitude = -yllcorner-dy
                            for col2 in range(plate_size):
                                col = col1*plate_size+col2
                                ## calculatelate longitude
                                dx = cellsize*(col+0.5)*max(nrows, ncols)/n
                                longitude = xllcorner+dx
                                ## Generate line for coordinate.
                                line = generate_line(
                                    y, row, col, refID, latitude, longitude,
                                    array_2D_materialIDs, array_2D_density,
                                    array_2D_buried,
                                    array_3D_designIDs, array_3D_materialIDs,
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
    y, row, col, refID, latitude, longitude,
    array_2D_materialIDs, array_2D_density, array_2D_buried,
    array_3D_designIDs, array_3D_materialIDs):

    if y >= array_2D_density[row][col]:
        return None

    ## skip parts of array not covered by GeoJSON array
    ## DO THIS BEFORE CONNECTING BRICKS!!!
    ## INSTEAD COLOR AS NEIGBOUR!!
    ## OR USE ALL 4 CORNERS INSTEAD OF JUST CENTER!!!
    ## Skip ocean.
    if array_2D_materialIDs[row][col] == 0:
        return None

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
    ## Language polygon, but no population density.
    if all([
        array_2D_materialIDs[row][col] > 0,
        array_2D_density[row][col] == 0,]):
        ## todo: generate map with and without!!!
        if all([
            latitude > tripoint_SW_lat,
            longitude > tripoint_SW_lon,
            latitude < tripoint_NE_lat,
            longitude < easternmost_lon,
            ]):
            print(latitude, longitude, row, col)
            array_2D_density[row][col] = 1
        else:
            print('tmp skip', row, col, latitude, longitude)
            return None

##                print(row,col,latitude,longitude)

    if array_3D_materialIDs[y][row][col] != 0:
        materialID = array_3D_materialIDs[y][row][col]
    ## white/buried
    elif y < array_2D_buried[row][col]:
        materialID = 1
    ## exposed/colored
    else:
        materialID = array_2D_materialIDs[row][col]
    ## 1x1 plate
    if array_3D_designIDs[y][row][col] == 0:
        designID = 3024  # 1x1 plate
    ## larger plate
    else:
        ## empty space for larger plate
        if array_3D_designIDs[y][row][col] == -1:
            return None
        ## larger plate
        else:
            designID = array_3D_designIDs[y][row][col]
            pass
        pass

    line = format_line(
        refID, designID, materialID, row, y, col)

    return line


def format_line(refID, designID, materialID, row, y, col):

    ## LEGO dimensions
    ## http://upload.wikimedia.org/wikipedia/commons/1/1a/Lego_dimensions.svg
    P = 0.8  # cm
    h = 0.32  # cm

    line = '        <Part'
    line += ' refID="{0}"'.format(refID)
    line += ' designID="{0}"'.format(designID)
    line += ' materialID="{0}"'.format(materialID)
    line += ' itemNos="{0}{1}"'.format(designID, materialID)
    designID_baseplate = 4186
    ## rotate / rotation
    if y % 2 == 1 and designID != designID_baseplate:  ## 2014sep13 2015may20
##        line += ' angle="0" ax="0" ay="1" az="0"'
        line += ' angle="90" ax="0" ay="1" az="0"' ## tmp!!! test!!! East to West (from same point as S/N)
        dtx = P * (d_designIDs2dim[designID]['dim1']-1)
        dtz = 0
##        line += ' angle="90" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
##        line += ' angle="270" ax="0" ay="1" az="0"' ## tmp!!! test!!! East to West (from same point as S/N)
##        line += ' angle="270" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
##        line += ' angle="0" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
    else:
        line += ' angle="0" ax="0" ay="1" az="0"'
        dtx = 0
        dtz = 0
##    line += ' tx="{0:.1f}"'.format(row*-P)
    line += ' tx="{0:.1f}"'.format(row*-P + dtx)
    line += ' ty="{0:.2f}"'.format(y*h)
##                    line += ' ty="{0}"'.format((y-z)*h)
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

        array_gpw = np.zeros((max(nrows, ncols), max(nrows, ncols)))
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
                    array_gpw[row][col] = 0
                elif s == '-9999':
                    array_gpw[row][col] = 0
                else:
                    array_gpw[row][col] = float(s)

    return array_gpw


if __name__ == '__main__':
    main()
