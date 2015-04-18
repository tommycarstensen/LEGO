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


def main():

    args = argparser()

##    lcm = nrows2*n/fractions.gcd(nrows2,n)
##    print(lcm,fractions.gcd(nrows2,n))
##    print(90*8/fractions.gcd(90,8),fractions.gcd(90,8))
##    print(8/fractions.gcd(90,8),90/fractions.gcd(90,8))
##    print(nrows2/fractions.gcd(nrows2,n))
##    stop

    affix1 = 'out_NASA2LEGO/{0}_{1:d}x{1:d}_zero{2:d}'.format(
        args.affix, args.n, args.zero)

    if not os.path.isfile('{}.npy'.format(affix1)):
        array_2D_density = asc2np(args, affix1)
    else:
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

    array_2D_buried = find_buried(array_2D_density, args.layers)

    array_3D_designIDs = find_connected_buried(
        args.n, array_2D_buried, args.layers, args.plate_size)

    ncols, nrows, xllcorner, yllcorner, cellsize = read_gpw_header(
        '{}.asc'.format(args.affix))

    array_2D_materialIDs = json2array(
        args.n, nrows, ncols, xllcorner, yllcorner, cellsize)

    array_3D_designIDs, array_3D_materialIDs = find_connected_exposed(
        array_2D_density, array_2D_buried, args.layers,
        array_3D_designIDs, array_2D_materialIDs)

    lxfml = '{}_y{:d}_{}.lxfml'.format(affix1, args.layers, args.norm)
    numpy2lxfml(
        array_2D_density, lxfml, array_2D_materialIDs, array_2D_buried,
        array_3D_designIDs, array_3D_materialIDs,
        nrows, ncols, xllcorner, yllcorner, cellsize, args.plate_size)

##    print(pcount_max)
##    print(np.amin(array_gpw))
##    from collections import Counter
##    print(Counter([float(x) for x in np.nditer(array_gpw)]))

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

    parser.add_argument('--layers', default = 27, type=int)

    parser.add_argument('--norm', default = 'log', choices=[
        'log10', 'log2', 'unity', 'log']) # unity is feature scaling

    args = parser.parse_args()

    # script fast if multiple of 2160
    args.n = n = nrows = ncols = args.plate_cnt*args.plate_size

    assert args.layers % 3 == 0

    return args


def asc2np(args, affix1):
    
##        array_gpw = read_gpw('afp00g.asc')
        array_gpw = read_gpw('{}.asc'.format(affix))

        array_2D_density = np.zeros((n, n))
        array_2D_density_cnt = np.zeros((n, n))

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
                        array_2D_density[x1][y1] += array_gpw[x3][y3]
                        if array_gpw[x3][y3] > 0:
                            array_2D_density_cnt[x1][y1] += 1

        if not zero:
            for row in range(n):
                for col in range(n):
                    if array_2D_density[row][col] == 0:
                        continue
                    array_2D_density[row][col] /= array_2D_density_cnt[row][col]

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

    d_len2designIDs = {1: 3024, 2: 3023, 3: 3623, 4: 3710, 6: 3666, 8: 3460}
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
            d_colors = {}
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
##            if layer >=15 and col == 182:
##                print(layer)
##                print(seq)
##                print(array_2D_materialIDs[98][col])
##                print(array_2D_buried[98][col])
##                print(array_2D_density[98][col])
##                print(array_3D_designIDs[layer][98][col])
##                stop
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
        if layer % 2 == 1 and j == 0:  ## sep13
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


def find_connected_buried(n, array_2D_buried, layers, plate_size):

    print('finding connected 1x1 plates and replacing them with larger plates')

    ## DO A WHILE LOOP UNTIL area_max IS EMPTY!!!
    ## TODO: FIRST LOOP 3 LAYERS (BRICK) AND THEN 1 LAYER (PLATES) - 12x24 brick...
    ## DO NOT CENTER BRICK IN A DIMENSION
    ## IF FACING OTHER PLATE IN THAT DIRECTION - ie IF AT EDGE

    array_3D_designIDs = np.zeros((layers, n, n), int)

    d_designIDs = {16: 91405, 8: 41539, 6: 3958, 4: 3031}

##    array_2D_density = np.array([
##    [0,0,0,0,0,0,1,1,0],[0,0,0,0,0,1,1,1,1],[0,0,1,1,1,1,1,1,1]])
##    layers = 1
##    array_2D_buried = array_2D_density
##    n1 = np.shape(array_2D_density)[0]
##    n2 = np.shape(array_2D_density)[1]

    ## loop from North to South
    for row1 in range(int(n/plate_size)):
        if False: ##sep20
            continue ## tmp!!! sep13
        ## loop from West to East
        for col1 in range(int(n/plate_size)):
            if row1 != 2 or col1 != 3: continue ## tmp!!! sep13
            print(row1, col1)
            area_max = {layer: {'area': 0} for layer in range(layers)}
            i = -1
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
                            if layer % 2 == 1 and 2+2==5:  ## sep13
                                array_3D_designIDs[layer][row-size+1][col] = d_designIDs[size]
                            else:
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


def json2array(n, nrows, ncols, xllcorner, yllcorner, cellsize):

    print('Converting GeoJSON polygons to array of LEGO color points')

    ## colors not for sale via Pick-a-Brick:
    ## 2 light gray / grey
    ## 27 dark gray / dark grey

    ## only available via service.lego.com/en-gb/replacementparts
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
        'Bantu': 1,  # White (51380)
        'Bantu / Bantu': 1,
        'Semitic: Arab, Bedouin': 194,  # Light Bluish Gray (47289)
        'Cushitic': 192,  # Reddish Brown (39398)
        'Chadic': 5,  # Tan (74477) / Brick Yellow
        ## not available via PAB
        'Saharan': 2,  # Light Gray (14325) / Grey
        'Saharan / Nilotic': 2,
        'Saharan / Cushitic': 2,
        'Chadic / Cushitic': 2,
        ## available via PAB
        'Nilotic': 106,  # Orange (35644)
        'Nilotic / Bantu': 106,
        'Nilotic / Bantoid': 106,
        'Chari-Nile / Nilotic': 106,
        ## not available via PAB
        'Berber': 138,  # Dark Tan (32656) / Sand Yellow
        'Northern Mande': 222,  # Bright Pink (34034) / Light Purple
        ## available via PAB
        'Voltaic': 28,  # Green (39516) / Dark Green
        'Kru': 26,  # Black (35741)
        'Adamawa-Ubangian': 23,  # Blue (70095) / Bright Blue
        'Bantoid': 21,  # Red (67551) / Bright Red
        'Chari-Nile': 199,  # Dark Bluish Gray (43621) / Dark Stone Grey (4210719)
        'Adamawa-Ubangian / Chari-Nile': 199,
        'West Atlantic': 24,  # Yellow (67832)
        ## not available via PAB
        'Fufulde': 119,  # Lime (26976) / Bright Yellowish Green
        'Chadic / Fufulde': 119,
        'Fufulde / Adamawa-Ubangia':  119,
        'San': 151,  # Sand Green (6989) / Sand Green
        'Songhai': 102,  # Medium Blue (17453) / Medium Blue
        'Gbaya': 141,  # Dark Green (10885) / Earth Green
        'Kwa': 321,  # Dark Azure (23275) / Dark Azur (Ivory Coast)
        'Maban': 268,  # Dark Purple (21162) / Medium Lilac
        'Southern Mande': 330,  # Olive Green (12496) / Olive Green (Ivory Coast)
        'Kordofanian': 27,  # Dark Gray (7760)
        'Sandawe': 154,  # Dark Red (7042) / New Dark Red
        'Khoi:  Nama, Bergdama': 25,  # Brown (5452) / Earth Orange
        'Fur': 135,  # Sand Blue (4126)
        }

    d_colors = {
        ## Afro-Asiatic
        ## Yellow
        'Semitic: Arab, Bedouin': 24,
        'Cushitic': 24,
        'Chadic': 24,
        'Berber': 24,
        'Chadic / Cushitic': 24,
        'Chadic / Fufulde': 24,
        ## Nilo-Saharan (Berta, *Fur*, Gumuz, Koman, Kuliak, Kunama, *Maban*,
        ## Saharan, Songhay, Central Sudanic, Eastern Sudanic)
        ## Red
        ## Nilotic (E: Turkana, Maasai; S: Kalenjin, Datooga; W: Luo; Burun)
        'Saharan': 21,
        'Saharan / Nilotic': 21,
        'Saharan / Cushitic': 21,
        'Nilotic': 21,
        'Nilotic / Bantu': 21,
        'Nilotic / Bantoid': 21,
        'Chari-Nile / Nilotic': 21,
        'Maban': 21,
        'Fur': 21,
        ## Niger-Congo A (not Bantu)
        ## Blue
        'Northern Mande': 23,
        'Voltaic': 23,
        'Kru': 23,
        'Adamawa-Ubangian': 23,
        'Bantoid': 23,
        'Chari-Nile': 23,
        'Adamawa-Ubangian / Chari-Nile': 23,
        'West Atlantic': 23,
        'Fufulde': 23,
        'Fufulde / Adamawa-Ubangia':  23,
        'Songhai': 23,
        'Kwa': 23,
        'Southern Mande': 23,
        'Kordofanian': 23,
        ## Niger-Congo B (Bantu)
        ## Green
        'Bantu': 28,
        'Bantu / Bantu': 28,
        'Gbaya': 28,
        ## Khoisan languages / Khoe languages
        ## Orange (do Tan instead?! more common!!!)
        'San': 106,
        'Sandawe': 106,
        'Khoi:  Nama, Bergdama': 106,
        }

    array_2D_materialIDs = np.zeros((n, n), int)

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
##            print(family, feature['properties']['ETHNICITY'])
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
        _den = cellsize*max(nrows,ncols)
        row_min = n*(min_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        row_max = n*(max_lat-yllcorner+cellsize*(ncols-nrows)/2)/_den-0.5
        col_min = n*(min_lon-xllcorner)/_den-0.5
        col_max = n*(max_lon-xllcorner)/_den-0.5
        within = False
        ## Afrikaans
        if family == 'Miscellaneous / Unclassified' and max_lat < -20:
            family = 'Afrikaans'
            color = 140  # earth blue
            color = 5  # tan
            color = 119  # lime
        ## loop from South to North
        for row in range(math.floor(row_min), math.ceil(row_max)):
            latitude = cellsize*(row+0.5)*max(nrows, ncols)/n+yllcorner
            latitude -= cellsize*(ncols-nrows)/2
            ## Bantu languages in Angola incorrectly assigned to the Kru family
            if family == 'Kru' and latitude < 0:
                family = 'Bantu'
                color = d_colors[family]
            ## loop from West to East
            for col in range(math.floor(col_min), math.ceil(col_max)):
                longitude = cellsize*(col+0.5)*max(nrows, ncols)/n+xllcorner
                ## incorrectly assigned to the Maban family
                if family == 'Maban' and longitude > 30:  # min_lon > 30
                    within = True
                    continue  # don't continue but re-assign!!! tmp!!!
                ## don't change from non-miscellaneous to miscellanous
                if (
                    family in ('', 'Miscellaneous / Unclassified') and
                    array_2D_materialIDs[n-row][col]):
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
                    array_2D_materialIDs[n-row][col] = color
                    pass
                ## polygon might not contain point rounded to nearest grid value
                ## hence add manually
                elif row_max-row_min < 2.5 or col_max-col_min < 2.5:
                    within = True
                    ## don't change color if not within polygon
                    if not array_2D_materialIDs[n-row][col] == 0:
                        continue
                    array_2D_materialIDs[n-row][col] = color
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
    array_2D_density, lxfml, array_2D_materialIDs, array_2D_buried,
    array_3D_designIDs, array_3D_materialIDs,
    nrows, ncols, xllcorner, yllcorner, cellsize, plate_size):

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
        ## loop from North to South
        for row1 in range(int(n/plate_size)):
            ## loop from West to East
            for col1 in range(int(n/plate_size)):
                with open(
                    '{}.{}.{}.lxfml'.format(
                        lxfml[:-len('.lxfml')], row1, col1), 'w') as lxfml_plate:
                    lxfml_plate.write(head('Africa.{}.{}'.format(row1, col1)))
                    for row2 in range(plate_size):
                        row = row1*plate_size+row2
                        dy = cellsize*(row+0.5)*max(nrows, ncols)/n
                        dy -= -cellsize*(ncols-nrows)/2
                        latitude = yllcorner+dy
                        if row % 10 == 0:
                            print('numpy2lxfml, row', row, 'of', n-1)
                        for col2 in range(plate_size):
                            col = col1*plate_size+col2
                            dx = cellsize*(col+0.5)*max(nrows, ncols)/n
                            longitude = xllcorner+dx
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

                            ## skip parts of array not covered by GeoJSON array
                            ## DO THIS BEFORE CONNECTING BRICKS!!!
                            ## INSTEAD COLOR AS NEIGBOUR!!
                            ## OR USE ALL 4 CORNERS INSTEAD OF JUST CENTER!!!
                            if array_2D_materialIDs[row][col] == 0:
                                continue  # tmp!!!

                            ## make sure bricks are present
                            ## for areas covered by language family polygons
                            if all([
                                array_2D_materialIDs[row][col] > 0,
                                array_2D_density[row][col] == 0]):
                                array_2D_density[row][col] = 1

            ##                print(row,col,latitude,longitude)

                            for y in range(int(array_2D_density[row][col])):
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
                                        continue
                                    ## larger plate
                                    else:
                                        designID = array_3D_designIDs[y][row][col]
                                        pass
                                    pass

            ##                    ## replace plates with bricks
            ##                    if array_2D_density[row][col] >= 3*(y//3)+3:
            ##                        if y%3 == 0:
            ##                            designID = 3005  # 1x1 brick
            ##                        else:
            ##                            continue
            ##                    if array_2D_buried[row][col] >= 3*(y//3)+3:
            ##                        if y%3 == 0:
            ##                            designID = 3005  # 1x1 brick
            ##                        else:
            ##                            continue

            ##                    if array_2D_buried[row][col] >= y:
            ##                        materialID = 1
            ##                    else:
            ##                        materialID = array_2D_materialIDs[row][col]

                                line = format_line(
                                    refID, designID, materialID, row, y, col)
                                f.write(line)
                                lxfml_plate.write(line)
                                refID += 1
                                ## continue loop over ty
                                continue
                            ## continue loop over col2
                            continue
                        ## continue loop over row2
                        continue
                    ## add grey 48x48 base plate
                    assert plate_size == 48
                    line = format_line(refID, 4186, 194, row, 0, col)
                    f.write(line)
                    lxfml_plate.write(line)
                    refID += 1
                    ## add tail to plate lxfml file
                    lxfml_plate.write(tail())
                    ## close plate lxfml file
                    lxfml_plate.close()
                    pass
                ## continue loop over tz
                continue
            ## continue loop over tx
            continue
        ## add tail to lxfml file
        f.write(tail())
        ## close file
        f.close()
        pass

    return


def format_line(refID, designID, materialID, row, y, col):

    line = '        <Part'
    line += ' refID="{0}"'.format(refID)
    line += ' designID="{0}"'.format(designID)
    line += ' materialID="{0}"'.format(materialID)
    line += ' itemNos="{0}{1}"'.format(designID, materialID)
    if y % 2 == 1 and designID != 4186:  ## sep13
##        line += ' angle="0" ax="0" ay="1" az="0"'
        line += ' angle="90" ax="0" ay="1" az="0"' ## tmp!!! test!!! East to West (from same point as S/N)
##        line += ' angle="90" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
##        line += ' angle="270" ax="0" ay="1" az="0"' ## tmp!!! test!!! East to West (from same point as S/N)
##        line += ' angle="270" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
##        line += ' angle="0" ax="0" ay="-1" az="0"' ## tmp!!! test!!! West to East
    else:
        line += ' angle="0" ax="0" ay="1" az="0"'
    line += ' tx="{0:.1f}"'.format(row*-0.8)
    line += ' ty="{0:.2f}"'.format(y*0.32)
##                    line += ' ty="{0}"'.format((y-z)*0.32)
    line += ' tz="{0:.1f}"'.format(col*0.8)
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
            latitude = yllcorner+(nrows-row+x1-.5)*2.5/60
            if row >= x1 and row < x2:
                cols = f.readline().rstrip().split()
                pass
            else:
                cols = ['0']*ncols
            N, W = '13.146,43.271'.split(',')
            N = float(N)
            W = float(W)
            Ndelta = 0.02
            Wdelta = 0.02
            Ndelta = 0.2
            Wdelta = 0.2
            if latitude < N+Ndelta and latitude > N-Ndelta:
                print()
##                print(latitude,longitude,cols)
            for col, s in enumerate(cols):
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
##                ## https://en.wikipedia.org/wiki/Extreme_points_of_Europe
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
##                ## do more accurate slope with more accurate coordinates
##                ## and then convert to grid
##                elif row >= 328 and row <= 372 and col > 1541+(row-328)/((372-328)/(1557-1541)):
##                    s = 0
##                ## Middle East (from South to North)
##                ## Gulf of Aqaba (Tiran Island ignored)
##                elif row >= 372 and row <= 407 and col > 1558+(row-372)/((406-372)/(1548-1558)):
##                    s = 0
##                elif row == 372 and col >= 1559: s = 0  # Eifat/Aqaba/Jordan/Israel
##                elif row == 410 and col >= 1566: s = 0
##                ## Red Sea (Hanish Islands ignored)
##                elif (
##                    row <= 730 and row >= 411 and
##                    col > 1706+(row-689)/((428-689)/(1555-1706)):
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
                if all(
                    [
                        latitude < N+Ndelta, latitude > N-Ndelta,
                        longitude > W-Wdelta, longitude < W+Wdelta]):
                    print(
                        s, row, col, '{.3f},{.3f}'.format(latitude, longitude))

##        ## first loop to parse values
##        array_gpw = np.resize(
##            np.fromstring(
##                f.read().replace('\n',' '), sep= ' '), (nrows, ncols))

    return array_gpw


if __name__ == '__main__':
    main()
