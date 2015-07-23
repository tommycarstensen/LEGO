#!/usr/bin/env python3

# Tommy Carstensen, 2011, 2015jun, 2015jul

# todo: finish replacements of neigbouring pieces

import os
import math
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def main():

    args = argparser()

    # parse PDB coordinates
    d_coords_PDB = parse_pdb(args)

    # convert Angstrom coordinates to brick coordinates
    d_layers, l_minmax_bricks, d_buried = coordAngstrom2coordBrick(
        args, d_coords_PDB)
    del d_coords_PDB  # save some memory

    plot_layers(args, d_layers, d_buried, l_minmax_bricks)

    # convert brick coordinates to bricks
    coordBrick2coordLDD(
        args, d_layers, l_minmax_bricks, d_buried,
        )
    del d_layers, l_minmax_bricks  # save some memory

    return


def plot_layers(args, d_layers, d_buried, l_minmax_bricks):

    ## N, C, O, H, P, prev layer
    l_prev = []
    l_curr = []
    ## plot previous layer with traingle
    ## plot C, N, O, P, H, buried
    ## C = black circle
    ## N = blue square
    ## O = red triangle
    ## H = white diamond
    ## P = orange circle
    ## buried = green triangle2
    ## support? = green circle
##    ax.set_xlim(l_minmax_bricks[0])
##    ax.set_ylim(l_minmax_bricks[1])
    midx = (l_minmax_bricks[0][1]+l_minmax_bricks[0][0])/2
    midy = (l_minmax_bricks[1][1]+l_minmax_bricks[1][0])/2
    dmax = max(
        l_minmax_bricks[0][1]-l_minmax_bricks[0][0],
        l_minmax_bricks[1][1]-l_minmax_bricks[1][0],
        )
##    plt.xlim(midx-dmax/2, midx+dmax/2)
##    plt.ylim(midy-dmax/2, midy+dmax/2)

##    ax.axis('equal')

##    ax.set_xlim(midx-dmax/2, midx+dmax/2)
##    ax.set_ylim(midy-dmax/2, midy+dmax/2)
##    ax.set_aspect(1)
##        max(l_minmax_bricks[0][0], l_minmax_bricks[1][0]),
##        max(l_minmax_bricks[0][1], l_minmax_bricks[1][1]),
##        )
##    plt.ylim(l_minmax_bricks[1])
##    plt.axis(l_minmax_bricks[0]+l_minmax_bricks[1])

##    plt.axis([midx-dmax/2, midx+dmax/2, midy-dmax/2, midy+dmax/2])
##    plt.axis('equal')
##    plt.tick_params(which='major', length=1)
##    ax.axis('scaled')
    
    d_colors = {
        'O': 'r',
        'N': 'b',
        'C': 'k',
        'S': 'y',
        'H': 'w',
        'P': 'orange',
        }
    d_markers = {'O': 'o', 'N': 's', 'C': 'v', 'S': 'o', 'H': 'D', 'P': 'x'}
    for i, layer in enumerate(sorted(d_layers.keys())):
        print('plt', i)
        fig, ax = plt.subplots()
        for element in d_layers[layer].keys():
            for x in d_layers[layer][element]['x'].keys():
                for z in d_layers[layer][element]['x'][x]:
                    size = 10
                    if d_buried[layer][x][z]['excl_diagonals'] == 7:
                        c = 'g'
                        marker = 'x'
                    else:
                        c = d_colors[element]
                        marker = d_markers[element]
                        if i > 0:
                            for element2 in d_layers[layer-1].keys():
                                try:
                                    l = d_layers[layer-1][element2]['x'][x]
                                except KeyError:
                                    continue
                                if z in l:
                                    bool_supported = True
##                                    ax.scatter(x, z, s=32, c='grey', alpha=0.9, marker='^')
                                    break
                            else:
                                bool_supported = False
                                size = 25
##                            if bool_supported:
##                                ax.scatter(x, z, s=32, c='g', alpha=0.9, marker='^')
                    ax.scatter(x, z, s=size, c=c, alpha=0.9, marker=marker)
##        print(i, layer, d_layers[layer])
##        print(d_buried[layer].keys())
##        print(d_layers[layer]['O'])
##        print(d_buried[layer][8])
##        print(l_minmax_bricks)
##        stop
##        ax.plot(x, y)
##        ax.set_xlim(midx-dmax/2, midx+dmax/2)
##        ax.set_ylim(midy-dmax/2, midy+dmax/2)
        ax.set_xlim(l_minmax_bricks[0])
        ax.set_ylim(l_minmax_bricks[1])
        ax.set_aspect(1)
        ax.xaxis.set_minor_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.grid(b=True, color='k', axis='both', which='major', linestyle='--')
        ax.grid(b=True, color='k', axis='both', which='minor', linestyle=':')
        plt.savefig('png/{}_{:d}_{:d}'.format(args.pdb, args.sizeratio, i))
        plt.clf()
        plt.close()

    return


def calculate_support_position(d_layers, l_minmax_bricks,):

    #
    # dictionary independent of elements
    #
    d_3d = {}
    for z_brick in range(int(l_minmax_bricks[2][0]), int(l_minmax_bricks[2][1])+1):
        d_3d[z_brick] = {}
        for element in d_layers[z_brick].keys():
            for x in d_layers[z_brick][element]['x'].keys():
                if x not in d_3d[z_brick].keys():
                    d_3d[z_brick][x] = []
                d_3d[z_brick][x] += d_layers[z_brick][element]['x'][x]

    #
    # supporters
    #
    d_support = {}
    for z_brick_above in range(
        int(l_minmax_bricks[2][0])+1,
        int(l_minmax_bricks[2][1])+1,
        ):
        d_support[z_brick_above] = {}

        # loop over row/column
        for x in d_3d[z_brick_above].keys():
            l_above = d_3d[z_brick_above][x]

            set_above = set(l_above)
            for z_brick_below in range(int(l_minmax_bricks[2][0]), z_brick_above):
                # nothing in position below
                if x not in d_3d[z_brick_below].keys():
                    d_support[z_brick_above][x] = set_above
                    continue
                l_below = d_3d[z_brick_below][x]
                set_above -= set(l_below)
                # everyting above is also below
                if len(set_above) == 0:
                    break

            # something above was not below
            # use points for support from bottom layer
            if len(set_above) > 0:
                d_support[z_brick_above][x] = set_above

    return d_support


def parse_pdb(args):

    '''this function reads the coordinates of a PDB file and returns:
d_coords_PDB[record][chain][res_no][atom_name] = {
'coord':coord, 'element_size':element_size,
'element_color':element_color,
}
'''

    pdb = args.pdb

    print(('reading and parsing {}.pdb'.format(pdb)))
    fd = open('pdb/{}.pdb'.format(pdb), 'r')
    lines = fd.readlines()
    fd.close()

    d_coords_PDB = {'ATOM': {}, 'HETATM': {}}

    l_minmax_PDB = [
        [1000., -1000.], [1000., -1000.], [1000., -1000.],
        ]

    for line in lines:
        record = line[:6].strip()
        if record == 'ENDMDL':
            break
        if record not in ['ATOM', 'HETATM']:
            continue
        res_name = line[17:20]
        chain = line[21]
        if res_name == 'HOH':
            continue
        element_size = element_color = line[76:78].strip()
        if args.exclude_hydrogens is True and element == 'H':
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        if args.rotate_x90 is True:
            z = float(line[38:46])
            y = -float(line[46:54])
        coord = [x, y, z]
        atom_name = line[12:16].strip()
        res_no = int(line[22:26])
        if chain not in d_coords_PDB[record].keys():
            d_coords_PDB[record][chain] = {}
        if res_no not in d_coords_PDB[record][chain].keys():
            d_coords_PDB[record][chain][res_no] = {}

        #
        # color by residue instead of by atom
        #
        if res_name == 'HEM':
            element_color = 'CL'
        if args.bool_color_by_residue is True:
            d_elements = {
                ' DC': 'CL',  # green
                ' DT': '192',  # reddish brown
                ' DA': '119',  # lime
                ' DG': 'S',  # yellow
                }
            element_color = d_elements[res_name]

        d_coords_PDB[record][chain][res_no][atom_name] = {
            'coord': coord, 'element_size': element_size,
            'element_color': element_color,
            }

        for i in range(3):
            if coord[i] < l_minmax_PDB[i][0]:
                l_minmax_PDB[i][0] = coord[i]
            if coord[i] > l_minmax_PDB[i][1]:
                l_minmax_PDB[i][1] = coord[i]

    return d_coords_PDB


def gcd(a, b):

    '''greatest common divisor using Euclid's algorithm'''

    while b > 1e-14:
        a, b = b, a % b

    return a


def lcm(a, b):

    '''least/lowest/smallest common multiple'''

    multiple = a*b // gcd(a, b)

    return multiple


def calc_brick_per_Angstrom(args):

    #
    # find least common multiple of brick height and length
    #
    height = 9.6  # mm
    width = 8.0  # mm
    multiplier = lcm(9.6, 8.0)

    # slightly slower method for finding least common multiple
    multiplier = 1
    height = 9.6  # mm
    width = 8  # mm
    while True:
        if multiplier*height % width == 0:
            break
        multiplier += 1
    multiplier *= height

    # plates
    bpa_x = bricks_per_angstrom_x = 4.  # 8mm width
    bpa_y = bricks_per_angstrom_y = 4.  # 8mm width
    bpa_z = bricks_per_angstrom_z = 10.  # 3.2mm height
    # bricks (9.6mm height) - 48mm per Angstrom = 480, 000, 000:1
    bpa_x = bricks_per_angstrom_x = multiplier/width  # 6, 8mm width
    bpa_y = bricks_per_angstrom_y = multiplier/width  # 6, 8mm width
    if args.brick_or_plate == 'brick':
        bpa_z = bricks_per_angstrom_z = multiplier/height  # 5, 9.6mm height (brick)
    elif args.brick_or_plate == 'plate':
        bpa_z = bricks_per_angstrom_z = 15.  # 3.2mm height (plate)

    #
    # additional scaling
    #
    bpa_x /= args.sizeratio
    bpa_y /= args.sizeratio
    bpa_z /= args.sizeratio

    # list of scaling factors
    l_bpa = [bpa_x, bpa_y, bpa_z]

    return l_bpa


def coordAngstrom2coordBrick(args, d_coords_PDB):

    pdb = args.pdb

    print('convert PDB coordinates to brick coordinates')

    #
    # set i/o file path
    #
    fn = 'txt/d_layers_{}_{}_{}'.format(
        pdb, args.sizeratio, args.brick_or_plate,)
    if args.bool_hollow is True:
        fn += '_hollow'
    if args.bool_color_by_residue is True:
        fn += '_colorbyresidue'
    if args.bool_grained is True:
        fn += '_flying1x1removal'
    fn += '.txt'

    #
    # set i/o file path
    #
    fn_buried = 'txt/d_buried_{}_{}_{}'.format(
        args.pdb, args.sizeratio, args.brick_or_plate,)
    if args.bool_grained is True:
        fn_buried += '_flying1x1removal'
    fn_buried += '.txt'

    if os.path.isfile(fn) and os.path.isfile(fn_buried):

        #
        # read
        #

        print(('reading', fn))
        fd = open(fn, 'r')
        s = fd.read()
        fd.close()
        d_layers = eval(s)
        del s

        fd = open(fn_buried, 'r')
        s = fd.read()
        fd.close()
        d_buried = eval(s)
        del s

    else:

        #
        # loop over PDB coordinates and check distance to surrounding grid points
        #
        d_vicinal = find_vicinal(args, d_coords_PDB,)

        #
        #
        #
        d_layers = append_bricks_by_atom_vicinity(d_vicinal,)

        #
        # remove pieces before determining what is buried and what is not
        #
        d_layers = remove_bricks_outside(args, d_layers)

        #
        # hollow option
        #
        if args.bool_hollow is True:

            d_layers, d_buried = remove_or_replace_bricks_inside(d_layers,)

        else:

            d_buried = find_buried(d_layers,)

        #
        # write output
        #
        print(('writing', fn))
        fd = open(fn, 'w')
        fd.write(str(d_layers))
        fd.close()

        fd = open(fn_buried, 'w')
        fd.write(str(d_buried))
        fd.close()

    #
    #
    #
    l_minmax_bricks = get_dimensions(args, d_layers,)

    return d_layers, l_minmax_bricks, d_buried


def remove_or_replace_bricks_inside(d_layers,):

    d_buried = find_buried(d_layers,)

    replace_atom = replace_atom

    # either replace buried bricks with black
    # or delete buried bricks
    for z_brick in d_buried.keys():

        if replace_atom not in d_layers[z_brick].keys():
            d_layers[z_brick][replace_atom] = {'x': {}, 'y': {}}

        for element in d_layers[z_brick].keys():

            l_xbricks = list(d_layers[z_brick][element]['x'].keys())
            l_xbricks.sort()
            for x_brick in l_xbricks:
                l_ybricks = list(d_layers[z_brick][element]['x'][x_brick])
                l_ybricks.sort()
                for y_brick in l_ybricks:

                    count_buried_excl_diagonals = (
                        d_buried[z_brick][x_brick][y_brick]['excl_diagonals']
                        )
                    count_buried = d_buried[z_brick][x_brick][y_brick]['all']

                    # delete
                    if args.bool_hollow is True and count_buried == 3*3*3:
                        d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                        d_layers[z_brick][element]['y'][y_brick].remove(x_brick)

    return d_layers, d_buried


def get_dimensions(args, d_layers,):

    if args.verbose is True:
        print('determining brick dimensions')
    l_minmax_bricks = [[1000., -1000.], [1000., -1000.], [1000., -1000.]]
    for z_brick in d_layers.keys():
        for element in d_layers[z_brick].keys():
            for x_brick in d_layers[z_brick][element]['x'].keys():
                for y_brick in d_layers[z_brick][element]['x'][x_brick]:
                    coord_brick = [x_brick, y_brick, z_brick]
                    for i in range(3):
                        if coord_brick[i] < l_minmax_bricks[i][0]:
                            l_minmax_bricks[i][0] = coord_brick[i]
                        if coord_brick[i] > l_minmax_bricks[i][1]:
                            l_minmax_bricks[i][1] = coord_brick[i]

    print(('Brick dimensions (bricks)', l_minmax_bricks))
    print(('Brick dimensions (model, cm)',))

    print(('width', 0.8*(l_minmax_bricks[0][1]-l_minmax_bricks[0][0]),))
    print(('width', 0.8*(l_minmax_bricks[1][1]-l_minmax_bricks[1][0]),))
    if args.brick_or_plate == 'plate':
        print(('height', 0.32*(l_minmax_bricks[2][1]-l_minmax_bricks[2][0])))
    elif args.brick_or_plate == 'brick':
        print(('height', 0.96*(l_minmax_bricks[2][1]-l_minmax_bricks[2][0])))

    return l_minmax_bricks


def remove_bricks_outside(args, d_layers,):

    print('removing "outside" bricks before determining which bricks are buried')
    if args.bool_grained is True:

        l_removal = []

        #
        # loop over sorted layers to remove bricks from the bottom and up
        #
        l_zbricks = list(d_layers.keys())
        l_zbricks.sort()
        for z_brick in l_zbricks:
            # Don't remove the bottom brick no matter what.
            if z_brick == min(d_layers.keys()):
                continue
            for element in d_layers[z_brick].keys():
                for x_brick in d_layers[z_brick][element]['x'].keys():
                    for y_brick in d_layers[z_brick][element]['x'][x_brick]:
                        count_buried = 0
                        count_buried_excl_diagonals = 0
                        count_above_below = 0
                        count_below = 0
                        count_plane = 0
                        for z_add in [-1, 0, 1]:
                            z2 = z_brick+z_add
                            if z_brick+z_add not in d_layers.keys():
                                continue
                            # loop over all neigbouring element types
                            for element2 in d_layers[z_brick+z_add].keys():
                                for x_add in [-1, 0, 1]:
                                    x2 = x_brick+x_add
                                    if (
                                        x_brick+x_add
                                        not in
                                        list(d_layers[z_brick+z_add][element2]['x'].keys())
                                        ):
                                        continue
                                    for y_add in [-1, 0, 1]:
                                        y2 = y_brick+y_add
                                        if (
                                            y_brick+y_add
                                            in
                                            d_layers[z_brick+z_add][element2]['x'][x_brick+x_add]
                                            ):
                                            count_buried += 1
                                            if (
                                                (y_add == 0 and z_add == 0) or
                                                (x_add == 0 and z_add == 0) or
                                                (x_add == 0 and y_add == 0)
                                                ):
                                                count_buried_excl_diagonals += 1
                                                if z_add == 0:
                                                    count_plane += 1
                                                elif z_add != 0:
                                                    count_above_below += 1
                                                    if z_add == -1:
                                                        count_below += 1
                                            continue

                        # "get rid of 1x1 hanging pieces"
                        # nothing *above or below* brick
                        # and only 1 brick in the plane around the brick
                        if count_above_below == 0 and count_plane == 1+1:
                            l_removal += [[z_brick, element, x_brick, y_brick]]
                        # "don't add 1x1s to new layers when building from the bottom and up"
                        # nothing in the plane around the brick
                        # and nothing *below* the brick
                        elif count_below == 0 and count_plane == 1:
                            d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                            d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
                            print(('remove2', z_brick, element))
                        # "get rid of pieces sticking out otherwise troubling the use of 2 knob bricks"
                        # nothing above or below brick
                        # (*after* removal of any bricks below;
                        # thus sorted z_brick values necessary
                        # and immediate removal in previous layer happened)
                        elif count_above_below == 0:
                            d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                            d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
                            print(('remove3', z_brick, element))

                # end of loop over elements

        for z_brick, element, x_brick, y_brick in l_removal:
            d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
            d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
            print(('remove1', z_brick, element))

    return d_layers


def find_vicinal(args, d_coords_PDB,):

    print('find grid points vicinal to atom coordinates')

    l_bpa = calc_brick_per_Angstrom(args)

    d_vicinal = {}

    for record in d_coords_PDB.keys():
        for chain in d_coords_PDB[record].keys():
            for res_no in d_coords_PDB[record][chain].keys():
                if args.verbose:
                    print(('chain', chain, 'residue', res_no))
                for atom_name in d_coords_PDB[record][chain][res_no].keys():
                    if record == 'HETATM':
                        print(('exclude', record, chain, res_no))
                    element_size = (
                        d_coords_PDB[record][chain][res_no][atom_name][
                            'element_size'
                            ]
                        )
                    element_color = (
                        d_coords_PDB[record][chain][res_no][atom_name][
                            'element_color'
                            ]
                        )
                    coord_PDB = (
                        d_coords_PDB[record][chain][res_no][atom_name][
                            'coord'
                            ]
                        )
                    # convert Angstrom coord to grid coord
                    # loop over vicinal grid coords (cube)
                    grid_brick = []
                    for i in range(3):
                        grid_brick += [[
                            int(math.floor(
                                (coord_PDB[i]-args.d_radii[element_size])*l_bpa[i]
                                )),
                            int(math.ceil(
                                (coord_PDB[i]+args.d_radii[element_size])*l_bpa[i]
                                )),
                            ]]
                    sq_radius = args.d_radii[element_size]**2
                    position = 0
                    positionok = 0

                    #
                    # loop over grid bricks and check distance from center of atom
                    #
                    for x_brick in range(grid_brick[0][0], grid_brick[0][1]+1,):
                        x_LEGO_Angstrom = x_brick/l_bpa[0]
                        for y_brick in range(grid_brick[1][0], grid_brick[1][1]+1,):
                            y_LEGO_Angstrom = y_brick/l_bpa[1]
                            for z_brick in range(grid_brick[2][0], grid_brick[2][1]+1,):
                                z_LEGO_Angstrom = z_brick/l_bpa[2]

                                coord_LEGO_Angstrom = [
                                    x_LEGO_Angstrom,
                                    y_LEGO_Angstrom,
                                    z_LEGO_Angstrom,
                                    ]

                                position += 1

                                sq_dist = 0
                                for i in range(3):
                                    sq_dist += (coord_LEGO_Angstrom[i]-coord_PDB[i])**2
                                # grid point too distant from center of sphere
                                if sq_dist > sq_radius:
                                    continue

                                if x_brick not in d_vicinal.keys():
                                    d_vicinal[x_brick] = {}
                                if y_brick not in d_vicinal[x_brick].keys():
                                    d_vicinal[x_brick][y_brick] = {}
                                if z_brick not in d_vicinal[x_brick][y_brick].keys():
                                    d_vicinal[x_brick][y_brick][z_brick] = {}
                                d_vicinal[x_brick][y_brick][z_brick][sq_dist] = element_color

    return d_vicinal


def append_bricks_by_atom_vicinity(d_vicinal,):

    print('appending bricks according to atom nearest to grid point')
    d_layers = {}
    l_x_bricks = list(d_vicinal.keys())
    l_x_bricks.sort()
    for x_brick in l_x_bricks:
        l_y_bricks = list(d_vicinal[x_brick].keys())
        l_y_bricks.sort()
        for y_brick in l_y_bricks:
            for z_brick in d_vicinal[x_brick][y_brick].keys():
                dist_min = min(d_vicinal[x_brick][y_brick][z_brick].keys())
                element = d_vicinal[x_brick][y_brick][z_brick][dist_min]

                #
                # append coordinate to dictionary of coordinates
                #
                if z_brick not in d_layers.keys():
                    d_layers[z_brick] = {}
                if element not in d_layers[z_brick].keys():
                    d_layers[z_brick][element] = {'x': {}, 'y': {}}
                if x_brick not in d_layers[z_brick][element]['x'].keys():
                    d_layers[z_brick][element]['x'][x_brick] = []
                if y_brick not in d_layers[z_brick][element]['y'].keys():
                    d_layers[z_brick][element]['y'][y_brick] = []
    #            d_layers[z_brick][element] += [[x_brick, y_brick]]
                # Why is it not redundant to add for both dimensions?
                d_layers[z_brick][element]['x'][x_brick] += [y_brick]
                d_layers[z_brick][element]['y'][y_brick] += [x_brick]

    return d_layers


def find_buried(d_layers,):

    print('preparing to remove "buried" bricks \
after removing "isolated" bricks on the outside')
    d_buried = {}
    for z_brick in d_layers.keys():
        d_buried[z_brick] = {}
        for element in d_layers[z_brick].keys():
            for x_brick in d_layers[z_brick][element]['x'].keys():
                if x_brick not in d_buried[z_brick].keys():
                    d_buried[z_brick][x_brick] = {}
                for y_brick in d_layers[z_brick][element]['x'][x_brick]:
                    count_buried = 0
                    count_buried_excl_diagonals = 0
                    for z_add in [-1, 0, 1]:
                        z2 = z_brick+z_add
                        if z_brick+z_add not in d_layers.keys():
                            continue
                        # loop over all neigbouring element types
                        for element2 in d_layers[z_brick+z_add].keys():
                            for x_add in [-1, 0, 1]:
                                x2 = x_brick+x_add
                                if x_brick+x_add not in d_layers[z_brick+z_add][element2]['x'].keys():
                                    continue
                                for y_add in [-1, 0, 1]:
                                    y2 = y_brick+y_add
                                    if (
                                        y_brick+y_add
                                        in
                                        d_layers[z_brick+z_add][element2]['x'][x_brick+x_add]
                                        ):
                                        count_buried += 1
                                        if (
                                            (y_add == 0 and z_add == 0) or
                                            (x_add == 0 and z_add == 0) or
                                            (x_add == 0 and y_add == 0)
                                            ):
                                            count_buried_excl_diagonals += 1
                                        continue
                    if count_buried == 3*3*3:
                        bool_buried = True
                    elif count_buried_excl_diagonals == 6+1:
                        bool_buried = True
                    else:
                        bool_buried = False

                    d_buried[z_brick][x_brick][y_brick] = {
                        'all': count_buried,
                        'excl_diagonals': count_buried_excl_diagonals,
                        }

    return d_buried


def check_if_buried(d_layers, z_brick, xy, consecutive, k1, k2):

    l_below_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick-1, k1, xy, consecutive,
        )
    l_above_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick+1, k1, xy, consecutive,
        )
    l_behind_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick, k1, xy-1, consecutive,
        )
    l_front_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick, k1, xy+1, consecutive,
        )
    l_left_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick, k2, consecutive[0]-1, consecutive,
        )
    l_right_all_elements = check_bricks_vicinal_all_elements(
        d_layers, z_brick, k2, consecutive[0]+1, consecutive,
        )
    l_yx = list(range(consecutive[0], consecutive[0]+consecutive[1]))

    if (
        # no bricks below
        len(set(l_yx) & set(l_below_all_elements)) == 0 and
        # no bricks above
        len(set(l_yx) & set(l_above_all_elements)) == 0
        ):
        print(consecutive)
        print(l_below_all_elements)
        print(l_above_all_elements)
        print(l_yx)
        stop_nothing_to_hang_on_to

    # above/below
    set_below_and_above = set(l_below_all_elements) & set(l_above_all_elements)
    set_exposed_vertical = set(l_yx)-set_below_and_above

    # behind/front
    set_behind_and_front = set(l_behind_all_elements) & set(l_front_all_elements)
    set_exposed_horizontal1 = set(l_yx)-set_behind_and_front

    # left/right
    None

    set_exposed = set_exposed_horizontal1 | set_exposed_vertical
    set_buried = set(l_yx)-set_exposed

    l_exposed = list(set_exposed)
    l_buried = list(set_buried)

    l_exposed.sort()
    l_buried.sort()

    return l_exposed, l_buried


def coordBrick2coordLDD(
    args, d_layers, l_minmax_bricks, d_buried,
    ):

    '''convert 1x1x1 brick coordinates to bricks'''

    #
    #
    #
    d_bricks_main = coords_to_lxfml(
        args, d_layers, l_minmax_bricks, d_buried,
        )

    d_bricks_main = check_vicinal_small_bricks(args, d_bricks_main)

    # convert dic of bricks to lines of bricks
    # after replacing small bricks with large bricks and after recoloring
    lines_lxfml_body, refID, steps = dic2lines(args, d_bricks_main)

    #
    # Add 48x48 grey base plate.
    #
    x_min = float('inf')
    x_max = float('-inf')
    y_min = float('inf')
    y_max = float('-inf')
    for layer in d_layers.keys():
        for element in d_layers[layer].keys():
            for x in d_layers[layer][element]['x'].keys():
                x_min = min(x, x_min)
                x_max = max(x, x_max)
                for y in d_layers[layer][element]['x'][x]:
                    y_min = min(y, y_min)
                    y_max = max(y, y_max)
    plate_size = 48
    designID_baseplate = 4186
    materialID_grey = 194
    refID += 1
    s = format_line(
        refID, designID_baseplate, materialID_grey, 0, 0,
        .8*(x_max+x_min-plate_size)/2, 0, .8*(y_max+y_min+plate_size)/2)
    lines_lxfml_body.append(s)

    #
    # Add transparent supports.
    #
    if args.pdb == '2dau':
        d_covered = {}
        for x in range(x_min, x_max+1):
            d_covered[x] = {}
            for y in range(y_min, y_max+1):
                d_covered[x][y] = False
        designID_1x2x5 = 2454
        designID_1x1x5 = 2453
        materialID_transparent = 40
        for i, layer in enumerate(sorted(d_layers.keys())):
            for element in d_layers[layer].keys():
                for x in d_layers[layer][element]['x'].keys():
                    for y in d_layers[layer][element]['x'][x]:
                        if all([
                            i % 5 == 0,
##                            i >= 10,
##                            element in ('O', 'P'),
                            x > (x_max+x_min-plate_size)/2,
                            x < (x_max+x_min+plate_size)/2,
                            y > (y_max+y_min-plate_size)/2,
                            y < (y_max+y_min+plate_size)/2,
                            d_covered[x][y] == False,
                            ]):
                            print(i,x,y,element)
                            for j in range(i // 5):
                                refID += 1
                                s = format_line(
                                    refID, designID_1x1x5,
                                    materialID_transparent,
                                    0, 0, .8*x, 5*.96*j, -.8*y,
                                    )
                                lines_lxfml_body.append(s)
                        d_covered[x][y] = True

    write_lxfml(args, lines_lxfml_body, steps)

    return


def coords_to_lxfml(args, d_layers, l_minmax_bricks, d_buried):

    print((int(l_minmax_bricks[2][1])-int(l_minmax_bricks[2][0])+1, 'layers'))

    volume = 0
    d_bricks_main = {}
    lines_lxfml_body = []

    # loop over layers
    for z_brick in range(
        int(l_minmax_bricks[2][0]),
        int(l_minmax_bricks[2][1])+1,
        ):

        layer = z_brick-l_minmax_bricks[2][0]

        # no bricks in layer
        if z_brick not in d_layers.keys():
            continue

        if z_brick not in d_layers.keys():
            print(('no bricks in layer', layer))
            stop
            continue

        if args.verbose is True:
            print(('layer', layer))

        # set keys
        # odd and even layers rotated relative to each other
        if layer % 2 == 0:
            k1 = k = 'y'
            k2 = 'x'
        else:
            k1 = k = 'x'
            k2 = 'y'

        #
        # support positions
        #
        d_support = {}

        # loop over element types
        for element in d_layers[z_brick].keys():
            materialID = args.d_materialIDs[element]

            l_consecutives_skip = []
            l_coords_skip = []
            l_xys = list(d_layers[z_brick][element][k].keys())
            l_xys.sort()
            # loop over sorted rows/columns in layer
            for xy in l_xys:

                l_coords = d_layers[z_brick][element][k][xy]
                l_consecutives = identify_consecutive_stretches(l_coords)
                l_exposed = []

                #
                # q: bricks of same materialID/element present in
                # the next (xy+1) or previous (xy-1) or next next (xy+2) row?
                # could connect to it, in case it's hanging/flying
                #
                if xy+1 in d_layers[z_brick][element][k].keys():
                    l_coords_next_same_element = d_layers[z_brick][element][k][xy+1]
                    l_consecutives_next_same_element = identify_consecutive_stretches(
                        l_coords_next_same_element,
                        )
                else:
                    l_consecutives_next_same_element = []

                if xy+2 in d_layers[z_brick][element][k].keys():
                    l_coords_nextnext_same_element = d_layers[z_brick][element][k][xy+2]
                    l_consecutives_nextnext_same_element = identify_consecutive_stretches(
                        l_coords_nextnext_same_element,
                        )
                else:
                    l_consecutives_nextnext_same_element = []

                #
                #
                # split consecutive stretches into which bricks?
                #
                #
                for consecutive in l_consecutives:

                    #
                    # was the previous brick a double width brick?
                    #
                    if [xy, consecutive] in l_consecutives_skip:
                        print(l_consecutives_skip)
                        stop1
                        continue

                    width = 1
                    # use next for consecutive_next_same_element loop
                    l = l_consecutives_next_same_element

                    l_consecutives_split = [consecutive]
                    bool_next_rotated = False

                    #
                    # accomodate for rotated 2x1 or 3x1 brick in current row,
                    # if next brick is a 1x1 brick not attached above or below
                    # thus necessary and not just to save bricks, which is done in a later step...
                    # thus not necessary if 1x1 brick connected above or below!!!
                    # splits l_consecutives_split into smaller fragments
                    #
                    if consecutive[1] != 1:

                        for consecutive_next_same_element in l:

                            if (
                                # next brick in layer is a 1x1 brick
                                consecutive_next_same_element[1] == 1 and
                                consecutive_next_same_element[0] in range(
                                    l_consecutives_split[-1][0],
                                    l_consecutives_split[-1][0]+l_consecutives_split[-1][1],
                                    )
                                ):

                                # anything above or below next brick of any element?
                                l_below_all_elements = check_bricks_vicinal_all_elements(
                                    d_layers, z_brick-1, k, xy+width, consecutive_next_same_element,
                                    )
                                l_above_all_elements = check_bricks_vicinal_all_elements(
                                    d_layers, z_brick+1, k, xy+width, consecutive_next_same_element,
                                    )
                                # anything above or below next brick of any element?
                                if len(l_below_all_elements) == 0 and len(l_above_all_elements) == 0:
                                    # split if anything above or below next brick of any element
                                    if len(l_consecutives_split) == 2:
                                        stop_tmp
                                    pos1 = l_consecutives_split[-1][0]
                                    pos2 = l_consecutives_split[-1][0]+consecutive[1]-1
                                    pos_break = consecutive_next_same_element[0]
                                    if pos_break == pos1:
                                        brick1 = [
                                            pos1, 1, 'rotated', width,
                                            ]
                                        brick2 = [
                                            pos_break+1, pos2-pos_break, 'not rotated', width,
                                            ]
                                        brick3 = None
                                        for xy_add in range(1, width+1):
                                            l_consecutives_skip += [[xy+xy_add, consecutive]]
                                    elif pos_break == pos2:
                                        brick1 = [
                                            pos1, pos_break-pos1, 'not rotated', width,
                                            ]
                                        brick2 = [
                                            pos_break, 1, 'rotated', width,
                                            ]
                                        brick3 = None
                                        for xy_add in range(1, width+1):
                                            l_consecutives_skip += [[xy+xy_add, consecutive]]
                                    else:
                                        brick1 = [
                                            pos1, pos_break-pos1, 'not rotated', width,
                                            ]
                                        brick2 = [
                                            pos_break, 1, 'rotated', width,
                                            ]
                                        brick3 = [
                                            pos_break+1, pos2-pos_break, 'not rotated', width,
                                            ]
                                        for xy_add in range(1, width+1):
                                            l_consecutives_skip += [[xy+xy_add, consecutive]]
                                    # even newer method
                                    if brick3:
                                        l_consecutives_split = (
                                            l_consecutives_split[:-1] +
                                            [brick1, brick2, brick3]
                                            )
                                    else:
                                        l_consecutives_split = (
                                            l_consecutives_split[:-1] +
                                            [brick1, brick2]
                                            )

                                bool_next_rotated = True

                                continue

                    # dont skip next line of bricks,
                    # because one or more positions of the current consecutive stretch were rotated
                    if bool_next_rotated is True and width == 2:
                        l_consecutives_skip.remove([xy+1, consecutive])
                        width = 1

                    #
                    # accomodate for rotated 2x1 or 3x1 brick in next row,
                    # if curent brick is a 1x1 brick not attached above or below
                    #
                    elif consecutive[1] == 1:
                        # accomodation happens in function split_consecutive_length_into_bricks
                        # so do nothing
                        pass

                    #
                    #
                    # split consecutive *stretches* into *available* bricks
                    #
                    #
                    for consecutive_split in l_consecutives_split:

                        bool_prev_skipped = False
                        # previous brick was skipped to acomodate for brick previous to it
                        if consecutive_split[1] == 1 and [xy-1, consecutive_split[0]] in l_coords_skip:
                            bool_prev_skipped = True

                        d_bricks, d_layers = split_consecutive_length_into_bricks(
                            args,
                            consecutive_split, materialID, width,
                            d_layers, z_brick, element, k, xy,
                            bool_next_rotated, bool_prev_skipped,
                            )

                        materialID = args.d_materialIDs[element]

                        #
                        # append bricks to lxfml
                        #
                        d_bricks_main = append_bricks(
                            args, d_bricks, materialID, width, k,
                            # position
                            layer, xy,
                            d_bricks_main,
                            )

                        # end of loop over split consecutives

                    # end of loop over consecutives

                # end of loop over rows/columns

            # end of loop over elements

        # end of loop over layers

    return d_bricks_main


def check_vicinal_small_bricks(args, d_bricks_main,):

    '''speed up this function'''

    print('small2large')

    if args.bool_replace_small_bricks_with_larger_bricks is True:

        # also replace 2 * 1x2 with 1x4 !!!
        # also replace 1x2 and 2x2 with 2x3
        # also replace 1x2 and 1x1 with 1x3

        for layer in sorted(d_bricks_main.keys()):

            if layer % 10 == 0:
                print('small2large', layer, len(d_bricks_main.keys()))

            # sorted by cheapest replacement (volume to price ratio)...
            for replacement in [
                # round 1
                # don't replace 1x1 and 1x2 with corner in 1st round, as I prefer 1x3
                # don't replace 1x1 and 1x3 with 1x3, as I prefer 1x4
                '1x8', '1x4', '1x2', '1x6',
                # round 2
                # new 2x2 to 2x4
                # new 1x3 to 2x3
                # new 1x2 to 2x2
                # new 1x2 to 1x4
                # new 1x2 and old 1x1 to 1x3
                # new 1x2 and old 1x1 to corner
                '1x8', '1x4', '1x2', '1x6', '1x3',
                # round 3
                # new 2x2 to 2x4
                # new 1x4 to 2x4
                # new 1x3 to 2x3
                # new 1x3 to 1x6
                '1x8', '1x4', '1x2', '1x6', '1x3',
                # round 4
                # new 1x4 to 2x4
                '2x8', '2x4', '2x6',
                '2x8', '2x3', '2x2',
                '2x8', '2x6', '2x4', '2x3', '2x2',
                # i need to fix the corner function...
                'corner',

                ]:

                #
                # loop over x
                #
                l_tx = list(d_bricks_main[layer]['tx'].keys())
                l_tx.sort()
                for tx in l_tx:

                    #
                    # loop over z
                    #
                    l_tz = list(d_bricks_main[layer]['tx'][tx]['tz'].keys())
                    l_tz.sort()
                    for tz in l_tz:

                        # tz was deleted previously in the loop
                        if tz not in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                            continue
                        # check that tx was not deleted in previous loop over tz2
                        if tx not in d_bricks_main[layer]['tx'].keys():
                            print('txtxtx')
                            stop
                            break

                        designID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID']
                        if designID not in [
                            3005,  # 1x1 (to 1x2, 1x2x2)
                            3004,  # 1x2 (to 1x2x2, 2x2)
                            3003,  # 2x2 (to 2x4)
                            3622,  # 1x3 (to 2x3)
                            3010,  # 1x4 (to 2x4)
                            ]:
                            continue

                        materialID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID']

                        angle = d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
                        ay = d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay']

                        #
                        # loop over vicinal x
                        #
                        for tx2 in range(tx-2, tx+2+1,):

                            # check that there is a brick in the position
                            if tx2 not in d_bricks_main[layer]['tx'].keys():
                                continue
                            # check that tz was not deleted in previous loop over tz2
                            if tz not in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                                break

                            #
                            # loop over vicinal z
                            #
                            for tz2 in range(tz-2, tz+2+1,):

                                # check that there is a brick in the position
                                if tz2 not in d_bricks_main[layer]['tx'][tx2]['tz'].keys():
                                    continue
                                # check that tx2 was not deleted in previous loop over tz2
                                if tx2 not in d_bricks_main[layer]['tx'].keys():
                                    print('tx2tx2tx2')
                                    stop3
                                    break
                                # check that tz was not deleted in previous loop over tz2
                                if tz not in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                                    break
                                # check that tx was not deleted in previous loop over tz2
                                if tx not in d_bricks_main[layer]['tx'].keys():
                                    print('xxx')
                                    STOPSTOP
                                    break

                                designID2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['designID']
                                materialID2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['materialID']
                                angle2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['angle']

                                # obviously bricks have to be of the same color to be connected
                                if materialID != materialID2:
                                    continue

                                # skip if brick1 == brick2
                                if tx == tx2 and tz == tz2:
                                    continue

                                if abs(tx-tx2) < 6 and abs(tz-tz2) < 6:

                                    if (
                                        replacement == '1x3' and
                                        (
                                            # 1x1 (3005) and 1x2 (3004)
                                            # not vice versa 3004 and 3005
                                            (designID == 3005 and designID2 == 3004)
                                            )
                                        ):

                                        if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
                                            continue

                                        d_bricks_main = replace_1x2_and_1x1_w_1x3_main(
                                            args, layer, angle, angle2, tx, tx2, tz, tz2,
                                            d_bricks_main, materialID,
                                            )

                                    elif (
                                        replacement == '2x3' and
                                        (designID == 3622 and designID2 == 3622) and
                                        (tx == tx2 or tz == tz2)
                                        ):

                                        if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
                                            continue

                                        if (
                                            (angle == 90 and angle2 == 0) or
                                            (angle == 0 and angle2 == 90)
                                            ):
                                            continue

                                        d_bricks_main = replace_1x3_and_1x3_w_2x3(
                                            layer, angle, angle2, tx, tx2, tz, tz2,
                                            d_bricks_main, materialID,
                                            )

                                    #
                                    # replace 1x1 (3005) and 1x2 (3004) by a corner
                                    #
                                    elif (
                                        replacement == 'corner' and
                                        (designID == 3005 and designID2 == 3004) and
                                        # corner pieces not in blue and orange
                                        materialID not in [23, 106]
                                        ):
                                        d_bricks_main = (
                                            replace_1x2_and_1x1_w_1x2x2_main(
                                                args, d_bricks_main, layer,
                                                angle, angle2, tx, tx2, tz, tz2,
                                                )
                                            )

                                    # replace 1x1 and 1x1 (3005) by a 1x2 (3004)
                                    elif (
                                        replacement == '1x2' and
                                        designID == 3005 and designID2 == 3005 and
                                        # same row or column
                                        (tx == tx2 or tz == tz2)
                                        ):

                                        if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
                                            continue

                                        d_bricks_main = replace_1x1_and_1x1_w_1x2(
                                            d_bricks_main, layer,
                                            angle, angle2, tx, tx2, tz, tz2
                                            )

                                    # replace 1x2 and 1x2 (3004) by a 2x2 (3003)
                                    elif (
                                        replacement == '2x2' and
                                        designID == 3004 and designID2 == 3004 and
                                        (tx == tx2 or tz == tz2) and
                                        angle == angle2
                                        ):

                                        replace_1x2_w_2x2_main(
                                            args, layer, angle, angle2, tx, tx2, tz, tz2,
                                            d_bricks_main, materialID,
                                            )

                                    # replace 2x2 and 2x2 (3003) by a 2x4
                                    elif (
                                        replacement == '2x4' and
                                        designID == 3003 and designID2 == 3003 and
                                        (tx == tx2 or tz == tz2)
                                        ):

                                        d_bricks_main = replace_2x2_w_2x4_main(
                                            layer, angle, angle2, tx, tx2, tz, tz2,
                                            d_bricks_main, materialID,
                                            )

                                    #
                                    # replace 1x4 and 1x4 (3010) by a 2x4 (3001)
                                    #
                                    elif (
                                        replacement == '2x4' and
                                        (designID == 3010 and designID2 == 3010) and
                                        (tx == tx2 or tz == tz2)
                                        ):

                                        d_bricks_main = replace_1x4_and_1x4_w_2x4_main(
                                            layer, angle, angle2, tx, tx2, tz, tz2,
                                            d_bricks_main, materialID,
                                            )

    print('run loop until no replacements')

    return d_bricks_main


def replace_1x2_w_2x2_main(
    args, layer, angle, angle2, tx, tx2, tz, tz2,
    d_bricks_main, materialID,
    ):

    '''one of those ugly functions with lots of if statements... rewrite!'''

    if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
        return d_bricks_main

    if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
        return d_bricks_main

    # bricks rotated relative to each other
    if (
        (angle == 90 and angle2 == 0) or
        (angle == 0 and angle2 == 90)
        ):
        return d_bricks_main

    if layer % 2 == 0 and angle == 0 and tx == tx2 and tz2 == tz+1:
        d_bricks_main = replace_1x2_w_2x2(
            args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=1, tz_add=1,
            designID=3003, materialID=materialID,
            )
    elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz2 == tz-1:
        d_bricks_main = replace_1x2_w_2x2(
            args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0,
            designID=3003, materialID=materialID,
            )
    elif layer % 2 == 0 and angle == 90 and tz == tz2 and tx == tx2+1:
        if args.verbose is True:
            print('skip 2x2 for now')
    elif layer % 2 == 0 and angle == 90 and tz == tz2 and tx == tx2-1:
        if args.verbose is True:
            print('skip 2x2 for now')
    elif layer % 2 == 1 and angle == 0 and tx == tx2 and tz2 == tz-1:
        if args.verbose is True:
            print('skip 2x2 for now')
    elif layer % 2 == 1 and angle == 0 and tx == tx2 and tz2 == tz+1:
        if args.verbose is True:
            print('skip 2x2 for now')
    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
        if args.verbose is True:
            print('skip 2x2 for now')
    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
        if args.verbose is True:
            print('skip 2x2 for now')
    else:
        print('new 2x2 case')
        print(layer, angle, angle2, tx, tx2, tz, tz2)
        stop

    return d_bricks_main


def replace_2x2_w_2x4_main(
    layer, angle, angle2, tx, tx2, tz, tz2, d_bricks_main, materialID,
    ):

    if layer % 2 == 0 and tx2 == tx and tz2 == tz-2:
        if angle == 90:
            stop
        if angle == 0 and angle2 == 90:
            return d_bricks_main
        print('111', angle, angle2, ay)
        d_bricks_main = replace_2x2_w_2x4(
            d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=1, tz_add=2, materialID=221,
            )  # dark green

    elif layer % 2 == 1 and tx2 == tx+2 and tz2 == tz:
        if angle == 0:
            stop
        d_bricks_main = replace_2x2_w_2x4(
            d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=1, materialID=221,
            )  # bright purple

    elif layer % 2 == 1 and tx2 == tx and tz2 == tz+2:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 1 and tx2 == tx and tz2 == tz-2:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 0 and tx2 == tx-2 and tz2 == tz:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 0 and tx2 == tx+2 and tz2 == tz:
        if verbose is True:
            stop
        return d_bricks_main

    elif layer % 2 == 0 and tx2 == tx and tz2 == tz+2:
        if angle == 0 and ay == 1:
            d_bricks_main = replace_2x2_w_2x4(
                d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=1, tz_add=2,
                )
        # 2x2 is newly added
        elif angle == 0 and angle2 == 90:  # this should be a general exclusion
            None
        # 2x2 is newly added
        elif angle == 90 and angle2 == 0:  # this should be a general exclusion
            None
        else:
            print(layer, angle, angle2, ay)
            stop

    elif layer % 2 == 1 and tx2 == tx-2 and tz2 == tz:
        if angle == 90:
            stop
        else:
            stop

    else:
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = 221
        d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['materialID'] = 221
        print('new 2x4 case')
        print(tx, tz, tx2, tz2)
        print(d_bricks_main[layer]['tx'][tx2]['tz'][tz2])
        print('replace 2 * 2x2 with 1 * 2x4')
        print('layer', layer)
        print('angle', angle)
        stop4

    return d_bricks_main


def replace_1x4_and_1x4_w_2x4_main(
    layer, angle, angle2, tx, tx2, tz, tz2, d_bricks_main, materialID,
    ):

    if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
        return d_bricks_main

    if (
        (angle == 90 and angle2 == 0) or
        (angle == 0 and angle2 == 90)
        ):
        return d_bricks_main

    if layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2-1:
        tz_new = max(tz, tz2)
        tx_new = tx
        pass
    elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2+1:
        tz_new = max(tz, tz2)
        tx_new = tx
        pass
    elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2-1:
        print(layer)
        stop1
        if verbose is True:
            print('skip 2x4 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2+1:
        print(layer)
        stop2
        if verbose is True:
            print('skip 2x4 replacement for now')
        return d_bricks_main

    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
        print(layer)
        stop3
        if verbose is True:
            print('skip 2x4 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
        tz_new = tz
        tx_new = tx
        pass
    elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2-1:
        if verbose is True:
            print('skip 2x4 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2+1:
        print(layer)
        stop5
        if verbose is True:
            print('skip 2x4 replacement for now')
        return d_bricks_main
    else:
        print(layer, angle, angle2)
        print(tx, tx2, tz, tz2)
        stop
    d_bricks_main = replace_1xx_and_1xx_w_2xx(
        d_bricks_main, layer,
        angle, angle2, tx, tx2, tz, tz2,
        tx_new, tz_new,
        materialID,
        designID_new=3001,  # 2x4
        )

    return d_bricks_main


def replace_1x3_and_1x3_w_2x3(
    layer, angle, angle2, tx, tx2, tz, tz2, d_bricks_main, materialID,
    ):

    if layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2-1:
        tz_new = max(tz, tz2)
        tx_new = tx
        pass
    elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2+1:
        tz_new = max(tz, tz2)
        tx_new = tx
        pass
    elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2-1:
        print(layer)
        stop1
        if verbose is True:
            print('skip 2x3 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2+1:
        print(layer)
        stop2
        if verbose is True:
            print('skip 2x3 replacement for now')
        return d_bricks_main

    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
        print(layer)
        stop3
        if verbose is True:
            print('skip 2x3 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
        tz_new = tz
        tx_new = tx
        pass
    elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2-1:
        if verbose is True:
            print('skip 2x3 replacement for now')
        return d_bricks_main
    elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2+1:
        print(layer)
#                                                stop5
        if verbose is True:
            print('skip 2x3 replacement for now')
        return d_bricks_main
    else:
        print(layer, angle, angle2)
        print(tx, tx2, tz, tz2)
        stop
    d_bricks_main = replace_1xx_and_1xx_w_2xx(
        d_bricks_main, layer,
        angle, angle2, tx, tx2, tz, tz2,
        tx_new, tz_new,
        materialID,
        designID_new=3002,  # 2x3
        )

    return d_bricks_main


def replace_1x2_and_1x1_w_1x3_main(
    args, layer, angle, angle2, tx, tx2, tz, tz2, d_bricks_main, materialID,
    ):

    '''# not finished!!!'''
    # not finished!!!

    l_returns = [
        [0, 0, -1, 0, 0],
        [0, -1, -1, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 1, 1, 0, 0],
        [1, 1, 1, 90, 0],
        [1, 0, -1, 90, 0],
        ]

    l_skips = [
        [0, 1, -1, 0, 0],
        [0, 1, 0, 0, 90],
        [1, 1, 0, 90, 90],
        [1, -1, 0, 90, 90],
        [1, -1, -1, 0, 90],
        [0, -1, 1, 0, 0],
        [0, 1, -1, 0, 90],
        [1, -1, -1, 90, 90],
        [1, 1, -1, 90, 90],
        [0, -1, 0, 0, 90],
        [1, -1, 1, 90, 90],
        [1, -1, -1, 90, 0],
        [1, 0, 1, 90, 0],
        [0, -1, -1, 0, 90],
        [0, 1, 1, 0, 90],
        [0, -1, 1, 0, 90],
        [1, 0, 1, 90, 90],
        [0, 1, 0, 0, 0],
        [1, -1, 1, 90, 0],
        [1, 1, -1, 90, 0],
        [1, 1, 1, 90, 90],
        [0, 0, -1, 0, 90],
        ]

    # don't replace (corner)
    if [layer % 2, tx2-tx, tz2-tz, angle, angle2] in l_returns:
        return d_bricks_main

    if [layer % 2, tx2-tx, tz2-tz, angle, angle2] in l_skips:
        if args.verbose is True:
            print('skip 1x3 for now')
        return d_bricks_main

    elif (
        layer % 2 == 0 and tx2 == tx and tz2 == tz+1 and
        angle == 0 and angle2 == 90
        ):
        tx_new = tx
        tz_new = tz+2
        pass
    elif (
        layer % 2 == 1 and tx2 == tx+1 and tz2 == tz and angle == 90 and angle2 == 0
        ):
        tx_new = tx
        tz_new = tz
        pass
    else:
        print('new 1x3 case', layer, 'tx', tx, tx2, tz, tz2, angle, angle2)
        print([layer % 2, tx2-tx, tz2-tz, angle, angle2])
        return d_bricks_main
#                                                stop
    d_bricks_main = replace_1x2_and_1x1_w_1x3(
        d_bricks_main, layer,
        angle, angle2, tx, tx2, tz, tz2,
        tx_new, tz_new,
        )

    return d_bricks_main


def dic2lines(args, d_bricks_main):

    lines_lxfml = []
    refID = 0
    steps = ''
    for layer in sorted(d_bricks_main.keys()):
        steps += '      <Step>\n'
        if args.brick_or_plate == 'brick':
            ty = 0.96*layer
        else:
            ty = 0.32*layer
        for tx in d_bricks_main[layer]['tx'].keys():
            for tz in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                refID += 1
                designID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID']
                materialID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID']
                ay = d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay']
                angle = d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
                s = format_line(
                    refID, designID, materialID, angle, ay, .8*tx, ty, .8*tz)
                lines_lxfml += [s]
                steps += '        <PartRef partRefID="{}"/>\n'.format(refID)
        steps += '      </Step>\n'

    return lines_lxfml, refID, steps


def format_line(refID, designID, materialID, angle, ay, tx, ty, tz):

    itemNo = '{:d}{:d}'.format(designID, materialID)

    # initiate
    s = '        <Part'
    # refID
    s += ' refID="{:d}"'.format(refID)
    # dimension
    s += ' designID="{:d}"'.format(designID)
    # color
    s += ' materialID="{:d}"'.format(materialID)
    # dimension & color
    s += ' itemNos="{:s}"'.format(itemNo)
    # rotation
    s += ' angle="{:d}"'.format(angle)
    s += ' ax="0" ay="{:d}" az="0"'.format(ay)

    # position
    s += ' tx="{:f}"'.format(tx)
    s += ' ty="{:f}"'.format(ty)
    s += ' tz="{:f}"'.format(tz)

    # terminate
    s += '/>\n'

    return s


def replace_1xx_and_1xx_w_2xx(
    d_bricks_main, layer,
    angle, angle2, tx, tx2, tz, tz2,
    tx_new, tz_new,
    materialID,
    designID_new,
    ):

    if tz_new not in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
            d_bricks_main[layer]['tx'][tx]['tz'][tz]
            )
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['materialID'] = (
        materialID
        )
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = (
        designID_new
        )  # 3x2
    # delete bricks from old pos
    if tx != tx_new or tz != tz_new:
        del d_bricks_main[layer]['tx'][tx]['tz'][tz]
    elif tx2 != tx_new or tz2 != tz_new:
        del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
    else:
        stop

    return d_bricks_main


def replace_1x2_and_1x1_w_1x3(

    d_bricks_main, layer,
    angle, angle2, tx, tx2, tz, tz2,
    tx_new, tz_new,
    ):

    '''this function is most likely ver wrong and needs to be expanded.
only tested for 2dau_3 layer 30'''

    if tx_new not in d_bricks_main[layer]['tx'].keys():
        d_bricks_main[layer]['tx'][tx_new] = {'tz': {}}
    if tz_new not in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
            d_bricks_main[layer]['tx'][tx]['tz'][tz]
            )
    # rotation of 3004 1x2
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = angle2
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3622
    # delete bricks from old pos
    if tx != tx_new or tz != tz_new:
        del d_bricks_main[layer]['tx'][tx]['tz'][tz]
    if tx2 != tx_new or tz2 != tz_new:
        del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]

    return d_bricks_main


def replace_1x1_and_1x1_w_1x2(
    d_bricks_main, layer,
    angle, angle2, tx, tx2, tz, tz2):

    return d_bricks_main  # tmp!!!

    if layer % 2 == 0 and tx == tx2 and tz == tz2-1:
        tx_new = tx
        tz_new = tz+1
        if tz_new not in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                d_bricks_main[layer]['tx'][tx]['tz'][tz]
                )
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 90
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
        # delete bricks from old pos
        if tx != tx_new or tz != tz_new:
            del d_bricks_main[layer]['tx'][tx]['tz'][tz]
        if tx2 != tx_new or tz2 != tz_new:
            del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
    elif layer % 2 == 1 and tx == tx2 and tz == tz2-1:
        tx_new = tx
        tz_new = tz
        if tz_new not in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                d_bricks_main[layer]['tx'][tx]['tz'][tz]
                )
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 90
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
        # delete bricks from old pos
        if tx != tx_new or tz != tz_new:
            del d_bricks_main[layer]['tx'][tx]['tz'][tz]
        if tx2 != tx_new or tz2 != tz_new:
            del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
    elif tx == tx2-1 and tz == tz2:
        tx_new = tx
        tz_new = tz
        if tz_new not in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                d_bricks_main[layer]['tx'][tx]['tz'][tz]
                )
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 0
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
        # delete bricks from old pos
        del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
    else:
        print(tx, tz, tx2, tz2)
        print(.8*tx, .8*tz, .8*tx2, .8*tz2)
        print(angle, layer)
        print('replace materialID with 106/221 in lxmfl file')
        stop_compare_with_and_without_replacement

    return d_bricks_main


def replace_1x2_and_1x1_w_1x2x2_main(
    args, d_bricks_main, layer,
    angle, angle2, tx, tx2, tz, tz2,
    ):

    # NB! Replacement of 1x2 and 1x1 by 1x2x2 corner bricks
    # will make the model more stable and reduce the number of bricks
    # but also increase the price

    if layer % 2 == 1 and tx2 == tx-1 and tz2 == tz-1:  # when not rotated
        if angle == 0:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=-1, tz_add=-1, angle=0, materialID=192,
                )  # reddish brown
        elif layer % 2 == 1 and angle == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=180, materialID=119,
                )  # bright yellowish green
        else:
            print(layer, angle)
            stop1a
    elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz-1:  # when not rotated
        if angle == 0:
            pass
        else:
            print(layer, angle)
            stopstop

    elif layer % 2 == 1 and tx2 == tx-1 and tz2 == tz:  # when not rot
        if angle == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=-1, tz_add=1, angle=270,
                )
        elif angle == 0:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=-1, tz_add=1, angle=270, materialID=221,
                )  # bright purple
    elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz:  # when rotated!!!
        if angle == 0:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=-1, tz_add=1, angle=90, materialID=221,
                )  # dark green
        else:
            print(layer, angle, angle2)
            stop

    elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz+1:  # when no rot
        replace_1x2_and_1x1_w_1x2x2(
            args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=270, materialID=221,
            )  # dark green
    elif layer % 2 == 1 and tx2 == tx-1 and tz2 == tz+1:  # no rota
        pass

    elif layer % 2 == 0 and tx2 == tx and tz2 == tz-1:  # when no rotation
        if angle == 0:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=90, materialID=192,
                )  # reddish brown
        else:
            print('replace 1x2, 1x1 with corner - check rotation')
            print(layer, angle)
            stopc
    elif layer % 2 == 0 and tx2 == tx and tz2 == tz+1:  # when no roation
        if angle == 0 and angle2 == 0:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=1, tz_add=1, angle=180, materialID=221,
                )  # bright purple
        elif angle == 0 and angle2 == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=20, angle=90, materialID=221,
                )  # bright purple
        else:
            print('replace 1x2, 1x1 with corner - check rotation')
            print(layer, angle, angle2)
            stopd

    elif layer % 2 == 1 and tx2 == tx+1 and tz2 == tz-1:  # kein rotation of 2x1!
        if angle == 90:
            pass
        else:
            print(layer, angle)
            stope
    elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz-1:  # when no rotatin
        pass

    elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz:
        if angle == 0 and angle2 == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=0, materialID=221,
                )  # bright purple
        elif angle == 0 and angle2 == 0:
            pass  # 1x3
        else:
            print(layer, angle, angle2)
            stop

    elif layer % 2 == 1 and tx2 == tx and tz2 == tz+1:
        if angle == 90:
            pass  # 1x3
        else:
            print(layer, angle)
            stop

    elif layer % 2 == 1 and tx2 == tx and tz2 == tz-1:
        if angle == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=270, materialID=221,
                )  # dark green
        else:
            print(layer, angle)
            stop

    elif layer % 2 == 1 and tx2 == tx+1 and tz2 == tz:
        if angle == 90:
            replace_1x2_and_1x1_w_1x2x2(
                args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add=0, tz_add=0, angle=0,
                )
        else:
            print('replace 1x2, 1x1 with corner - check rotation')
            print(layer, angle)
            stopf

    elif tx2 == tx+1 and tz2 == tz+1:  # when no rota
        pass

    else:
        pass

    return d_bricks_main


def replace_1x2_w_2x2(
    args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add, tz_add, designID, materialID=221,
    ):

    # find consecutive occupants
    i_max = int(min(4, max(args.d_designIDs[args.brick_or_plate][2])))
    for i in range(
        1,
        i_max,
        ):

        bool_break = False
        if tx+i*(tx2-tx) in d_bricks_main[layer]['tx'].keys():
            if tz+i*(tz2-tz) in d_bricks_main[layer]['tx'][tx+i*(tx2-tx)]['tz'].keys():
                if d_bricks_main[layer][
                    'tx'
                    ][tx+i*(tx2-tx)][
                        'tz'
                        ][tz+i*(tz2-tz)][
                            'designID'
                            ] == 3004:  # 1x2
                    del d_bricks_main[layer]['tx'][tx+i*(tx2-tx)]['tz'][tz+i*(tz2-tz)]
                    continue
                else:
                    bool_break = True
                    break
            else:
                bool_break = True
                break
        else:
            bool_break = True
            break

    if bool_break is True:
        i -= 1

    tx_new = tx+tx_add
    tz_new = tz+tz_add

    if d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 90:
        tz_new += (i-1)*(tz2-tz)
    elif (
        layer % 2 == 0 and tx2 == tx and tz2 == tz+1 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
        ):
        tz_new += (i-1)*(tz2-tz)
    elif (
        layer % 2 == 0 and tx2 == tx and tz2 == tz-1 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
        ):
        tz_new += (i-1)*(tz2-tz)
    else:
        tx_new += (i-1)*(tx2-tx)

    # change brick
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = (
        args.d_designIDs[args.brick_or_plate][2][i+1]
        )  # 2x?
#    if recolor is True:
#        d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = materialID
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = abs(
        90-d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
        )
    # copy brick to new pos
    if tx_new not in d_bricks_main[layer]['tx'].keys():
        d_bricks_main[layer]['tx'][tx_new] = {'tz': {}}
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
        d_bricks_main[layer]['tx'][tx]['tz'][tz]
        )
    # delete brick from old pos
    del d_bricks_main[layer]['tx'][tx]['tz'][tz]

    return d_bricks_main


def replace_1x2_and_1x1_w_1x2x2(
    args, d_bricks_main, layer, tx, tz, tx2, tz2, tx_add, tz_add, angle, materialID=221,
    ):

    '''replace 1x2 and 1x1 by a "corner" brick'''

    # delete occupant/neighbour brick
    del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]

    # change brick
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = 2357
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = angle

    # copy brick to new pos
    if tx+tx_add not in d_bricks_main[layer]['tx'].keys():
        d_bricks_main[layer]['tx'][tx+tx_add] = {'tz': {}}  # tmp!!!
    d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add] = dict(
        d_bricks_main[layer]['tx'][tx]['tz'][tz]
        )

    # delete brick from old pos
    if not (tx_add == 0 and tz_add == 0):
        del d_bricks_main[layer]['tx'][tx]['tz'][tz]

    # verbose when finding incorrectly placed bricks
    d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add][
        'materialID'
        ] = materialID

    if args.verbose is True:
        d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add][
            'materialID'
            ] = materialID

    return d_bricks_main


def replace_2x2_w_2x4(
    d_bricks_main, layer, tx, tz, tx2, tz2, tx_add, tz_add, materialID=1,
    ):

    # find consecutive occupants
    i_max = int(max(d_designIDs[args.brick_or_plate][2])/2)
    for i in range(
        1,
        i_max,
        ):

        bool_break = False
        if tx+i*(tx2-tx) in d_bricks_main[layer]['tx'].keys():
            if tz+i*(tz2-tz) in d_bricks_main[layer]['tx'][tx+i*(tx2-tx)]['tz'].keys():
                if d_bricks_main[layer][
                    'tx'][tx+i*(tx2-tx)][
                        'tz'][tz+i*(tz2-tz)][
                            'designID'] == 3003:  # 2x2
                    del d_bricks_main[layer]['tx'][tx+i*(tx2-tx)]['tz'][tz+i*(tz2-tz)]
                    continue
                else:
                    bool_break = True
                    break
            else:
                bool_break = True
                break
        else:
            bool_break = True
            break

    if bool_break is True:
        i -= 1

    tx_new = tx+tx_add
    tz_new = tz+tz_add

    if d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 90:
        tz_new += (i-1)*(tz2-tz)
    elif (
        layer % 2 == 0 and tx2 == tx and tz2 == tz+2 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
        ):
        tz_new += (i-1)*(tz2-tz)
    elif (
        layer % 2 == 0 and tx2 == tx and tz2 == tz-2 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0 and
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
        ):
        tz_new += (i-1)*(tz2-tz)
    else:
        tx_new += (i-1)*(tx2-tx)

    # change brick
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = (
        d_designIDs[args.brick_or_plate][2][2*(i+1)]
        )  # 2x?
#    if recolor is True:
#        d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = materialID
    d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = abs(
        90-d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
        )
    # copy brick to new pos
    if tx_new not in d_bricks_main[layer]['tx'].keys():
        d_bricks_main[layer]['tx'][tx_new] = {'tz': {}}
    d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
        d_bricks_main[layer]['tx'][tx]['tz'][tz]
        )
    # delete brick from old pos
    del d_bricks_main[layer]['tx'][tx]['tz'][tz]

    return d_bricks_main


def write_lxfml(args, lines_lxfml_body, steps):

    # I need to add cameras at correct positions...
    # get those camera positions from brick positions...

    l_head1 = ['''<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<LXFML versionMajor="4" versionMinor="0" name="{}">
  <Meta>
    <Application name="LEGO Digital Designer" versionMajor="4" versionMinor="3"/>
    <Brand name="LDD"/>
    <BrickSet version="1564.2"/>
  </Meta>
'''.format(args.pdb)]

    l_cameras = ['''  <Cameras>
    <Camera refID="0" fieldOfView="80" distance="200" angle="60" \
ax="-0.6" ay="-0.6" az="-0.25" tx="0.0" ty="0.0" tz="0.0"/>
  </Cameras>
''']

    l_head2 = ['''  <Scene cameraRefID="0">
    <Model>
      <Group refID="0" angle="0" ax="0" ay="1" az="0" tx="0" ty="0" tz="0">
''']

    l_tail1 = ['''      </Group>
    </Model>
  </Scene>
  <BuildingInstructions>
    <BuildingInstruction>
{}
    </BuildingInstruction>
  </BuildingInstructions>
</LXFML>
'''.format(steps)]

    lines = l_head1+l_cameras+l_head2+lines_lxfml_body+l_tail1

    fn = os.path.join('lxfml', '{}_{}_{}'.format(
        args.pdb, args.sizeratio, args.brick_or_plate,))
    if args.bool_hollow is True:
        fn += '_hollow'
    if args.bool_color_by_residue is True:
        fn += '_colorbyresidue'
    fn += '.lxfml'
    fd = open(fn, 'w')
    fd.writelines(lines)
    fd.close()

    return


def split_consecutive_length_into_bricks(
    args,
    consecutive, materialID, width,
    d_layers, z_brick, element, k, xy,
    bool_next_rotated, bool_prev_skipped,
    ):

    d_bricks = {
        'rotated_before': [], 'rotated_after': [],
        'normal': [], 'flying': []}

    l_lengths = args.d_brick_sizes[args.brick_or_plate][materialID][width]
    l_lengths.sort()
    max_length = l_lengths[-1]

    # bricks exist
    if (
        consecutive[1] > 1 and
        consecutive[1] in args.d_brick_sizes[args.brick_or_plate][materialID][width]
        ):
        d_bricks['normal'] += [consecutive]
    # odd brick length higher than 1 for which brick exists,
    # if combined with 3x1 brick
    elif (
        consecutive[1] > 1 and
        consecutive[1]-3 in args.d_brick_sizes[args.brick_or_plate][materialID][width]
        ):
        d_bricks['normal'] += [[consecutive[0], 3], [consecutive[0]+3, consecutive[1]-3]]

    #
    # current brick is a 1x1 brick. anything above or below to attach to?
    #
    elif consecutive[1] in [1]:

        if bool_next_rotated is True and consecutive[2] == 'rotated':
            l_below_all_elements = check_bricks_vicinal_all_elements(
                d_layers, z_brick-1, k, xy+consecutive[3], consecutive,
                )
            l_above_all_elements = check_bricks_vicinal_all_elements(
                d_layers, z_brick+1, k, xy+consecutive[3], consecutive,
                )
        else:
            l_below_all_elements = check_bricks_vicinal_all_elements(
                d_layers, z_brick-1, k, xy, consecutive,
                )
            l_above_all_elements = check_bricks_vicinal_all_elements(
                d_layers, z_brick+1, k, xy, consecutive,
                )
        # more options for attachment if width is 2 (xy+1 could also be xy-1...)
            if width == 2:
                l_below_all_elements += check_bricks_vicinal_all_elements(
                    d_layers, z_brick-1, k, xy+1, consecutive,
                    )
                l_above_all_elements += check_bricks_vicinal_all_elements(
                    d_layers, z_brick+1, k, xy+1, consecutive,
                    )

        # nothing above or below current 1x1 brick
        if len(l_below_all_elements) == 0 and len(l_above_all_elements) == 0:

            if args.verbose is True:
                print('use 2x1 brick reversed', end=' ')
                print(', element =', element, ', k =', k, end=' ')
                print(', xy =', xy, consecutive, ', refID =', refID)

            bool_flying = False

            # already accomodated for rotated 2x1 or 3x1 brick in previous row
            if (
                # same element in previous row
                xy-1 in d_layers[z_brick][element][k].keys() and
                # same element in same position in previous row
                consecutive[0] in d_layers[z_brick][element][k][xy-1] and
                # not accomodating in next row after being accomodated itself?
                bool_next_rotated is False
                ):
                before_or_after = 'before'
                pass
            # accomodate for rotated 2x1 or 3x1 brick in next row,
            # if current brick is a 1x1 brick not attached above or below
            elif (
                # same element in next row
                xy+1 in d_layers[z_brick][element][k].keys() and
                # same element in same position in next row
                consecutive[0] in d_layers[z_brick][element][k][xy+1]
                ):
                # remove same element in same position in next row
                d_layers[z_brick][element][k][xy+1].remove(consecutive[0])
                before_or_after = 'after'
            # no brick of any element above or below
            # and no bricks in the same position before or after!
            # 1x1 brick will hang no matter what!
            else:
                bool_flying = True
                if args.verbose is True:
                    print('cant connect the brick to anything! flying brick!')
                    print(xy, 'not below', list(d_layers[z_brick-1][element][k].keys()))
                    print(xy, 'not above', list(d_layers[z_brick+1][element][k].keys()))
                    print(xy-1, 'not before', list(d_layers[z_brick][element][k].keys()))
                    print(xy+1, 'not after', list(d_layers[z_brick][element][k].keys()))
                    print('no consecutive bricks', end=' ')
                    print(d_layers[z_brick][element][k][xy], '-->', consecutive)

            if bool_flying is False:
                if before_or_after == 'before':
                    # previous brick was skipped
                    # to acomodate for brick previous to it
                    if bool_prev_skipped is True:
                        width += 1
                    d_bricks['rotated_before'] += [[consecutive[0], width+1]]
                elif before_or_after == 'after':
                    d_bricks['rotated_after'] += [[consecutive[0], width+1]]
            else:
                d_bricks['flying'] += [[consecutive[0], 1]]

        # something above or below current xy
        # elif len(l_below_all_elements) != 0 or
        # len(l_above_all_elements) != 0:
        else:

            d_bricks['normal'] += [consecutive]
            print('ggg', consecutive)

    # brick doesn't exist
    elif (
        consecutive[1]
        not in
        args.d_brick_sizes[args.brick_or_plate][materialID][width]
        ):
        # even consecutive length
        if consecutive[1] % 2 == 0:
            length = 0
        # odd consecutive length
        else:
            length = 3
            d_bricks['normal'] += [[consecutive[0], 3]]
        max_length = max(
            args.d_brick_sizes[args.brick_or_plate][materialID][width])
        while length < consecutive[1]:
            # avoid using length 2!!!
            # check that next piece is not a 2!!! remainder length...
            if max_length < consecutive[1]-length:
                brick_length = max_length
            elif (
                consecutive[1]-length
                not in
                args.d_brick_sizes[args.brick_or_plate][materialID][width]
                ):
                brick_length = l_lengths[-2]
            else:
                brick_length = consecutive[1]-length
            d_bricks['normal'] += [[consecutive[0]+length, brick_length]]
            length += brick_length
    else:
        print()
        print(consecutive)
        print(materialID, element)
        stop

    return d_bricks, d_layers


def modify_consecutives_exposed(l_consecutives_exposed,):

    # One stretch of consecutive exposed bricks.
    if len(l_consecutives_exposed) == 1:
        pass

    # Two stretches of consecutive exposed bricks.
    elif len(l_consecutives_exposed) == 2:
        # two 1x1 bricks
        if (
            l_consecutives_exposed[0][1] == 1 and
            l_consecutives_exposed[1][1] == 1):
            stop_expected
        # one 1x1 brick
        elif (
            l_consecutives_exposed[0][1] == 1 and
            l_consecutives_exposed[0][1] > 1):
            print(l_consecutives_exposed)
            print(l_consecutives_exposed[0][0]+l_consecutives_exposed[0][1]-1)
            print(l_consecutives_exposed[1][0]-1)
            stop1
        # one 1x1 brick
        elif (
            l_consecutives_exposed[0][1] > 1 and
            l_consecutives_exposed[1][1] == 1):
            print(l_consecutives_exposed)
            # room for a 2x1 brick?
            if (
                (
                    l_consecutives_exposed[0][0] +
                    l_consecutives_exposed[0][1] - 1
                    ) < l_consecutives_exposed[1][0]-1
                ):
                l_consecutives_exposed[1] = [l_consecutives_exposed[1][0], 2]
        # no 1x1 bricks
        else:
            pass

    # Three stretches of consecutive exposed bricks.
    elif len(l_consecutives_exposed) == 3:
        for i in range(len(l_consecutives_exposed)):
            if l_consecutives_exposed[i][1] == 1:
                print(l_consecutives_exposed)
                stop

    # Three or more stretches of consecutive exposed bricks???
    else:
        print(l_consecutives_exposed)
        stop_expected

    return l_consecutives_exposed


def append_bricks(
    args,
    # input
    d_bricks, materialID, width, k,
    layer, xy,
##    # output (append)
##    lines_lxfml_body,
    d_bricks_main,
    ):

    if layer % 2 == 0:
        k1 = k = 'y'
        k2 = 'x'
    else:
        k1 = k = 'x'
        k2 = 'y'

    for k_rot in ['normal', 'rotated_before', 'rotated_after', 'flying']:

        # change xy
        if k_rot == 'rotated':
            xy -= 1

        for brick in d_bricks[k_rot]:
            length = brick_len = brick[1]  # use to determine designID
            brick_pos = brick[0]  # use to determine brick position
            designID = args.d_designIDs[args.brick_or_plate][width][brick_len]

            #
            # bricks aligned with "one" axis
            #
            if k == 'y':
                tx = brick_pos
                tz = -xy
                ay = 1
                angle = 0
                dangle = 90
                if k_rot == 'rotated_before':
                    tz += 1
                    tx -= 1
                    tx += 1
                    tz -= 1
                elif k_rot == 'rotated_after':
                    tz -= 1
                    tx += 2
                    tx += 20

            #
            # bricks aligned with "another" axis
            #
            elif k == 'x':
                tx = xy
                tz = -brick_pos-(brick_len-1)
                ay = -1
                angle = 90
                dangle = -90
                if k_rot == 'rotated_before':
                    tz += (length-1)
                    tx -= (length-1)
                elif k_rot == 'rotated_after':
                    tz += 1

            if k_rot in ['rotated_before', 'rotated_after']:
                angle += dangle
                ay *= -1

            if brick_len == 1 and width == 2:
                angle += dangle

#            volume += width*brick_len

##            ## Redundant?
##            lines_lxfml_body, refID = append_lxfml_line(
##                args, designID, materialID, refID,
##                angle, ay, tx, layer, tz, lines_lxfml_body)

            if layer not in d_bricks_main.keys():
                d_bricks_main[layer] = {'tx': {}}
            if tx not in d_bricks_main[layer]['tx'].keys():
                d_bricks_main[layer]['tx'][tx] = {'tz': {}}
            d_bricks_main[layer]['tx'][tx]['tz'][tz] = {
                'designID': designID, 'materialID': materialID,
                'angle': angle, 'ay': ay,
                }

    return d_bricks_main


##def append_lxfml_line(
##    args, designID, materialID, refID, angle, ay, tx, layer, tz,
##    lines_lxfml_body,
##    ):
##
##    itemNo = '{:d}{:d}'.format(designID, materialID)
##
##    # initiate
##    s = '        <Part'
##    # refID
##    s += ' refID="{:d}"'.format(refID)
##    # dimension
##    s += ' designID="{}"'.format(designID)
##    # color
##    s += ' materialID="{}"'.format(materialID)
##    # dimension & color
##    s += ' itemNos="{}"'.format(itemNo)
##    # rotation
##    s += ' angle="{:d}"'.format(angle)
##    s += ' ax="0" ay="{:d}" az="0"'.format(ay)
##
##    # position
##    s += ' tx="{}"'.format(0.8 * tx)
##
##    if args.brick_or_plate == 'brick':
##        s += ' ty="{}"'.format(0.96 * layer)
##    elif args.brick_or_plate == 'plate':
##        s += ' ty="{}"'.format(0.32 * layer)
##
##    s += ' tz="{}"'.format(0.8 * tz)
##
##    # terminate
##    s += '/>\n'
##
##    # append line
##    lines_lxfml_body += [s]
##
##    refID += 1
##
##    return lines_lxfml_body, refID


def identify_consecutive_stretches(l_coords):

    l_consecutives = []
    j = 0
    i0 = 0
    for i in range(len(l_coords)):
        brick_length = 0
        if i+1 == len(l_coords):
            l_consecutives += [[l_coords[i0], i-i0+1]]
            break
        elif l_coords[i+1] != l_coords[i]+1:
            l_consecutives += [[l_coords[i0], i-i0+1]]
            i0 = i+1
            continue
        else:
            continue

    return l_consecutives


def check_bricks_vicinal_all_elements(
    d_layers, z_brick, k, xy, consecutive):

    if z_brick not in d_layers.keys():
        l_shared = []
    else:
        l_yx = list(range(consecutive[0], consecutive[0] + consecutive[1]))
        l_shared = []
        for element in d_layers[z_brick].keys():
            if xy not in d_layers[z_brick][element][k].keys():
                continue
            set_shared = set(l_yx) & set(d_layers[z_brick][element][k][xy])
            l_shared += list(set_shared)

    return l_shared


def argparser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--affix', default='afds00g',
        choices=['afp00ag', 'afp00g', 'afds00ag', 'afds00g'])

    parser.add_argument('--pdb')

    parser.add_argument('--sizeratio', type=int)

    parser.add_argument(
        '--brick_or_plate', default='brick', choices=['brick', 'plate'])

    parser.add_argument('--exclude_hydrogens', action='store_true')
    parser.add_argument('--bool_hollow', action='store_true')
    parser.add_argument(
        '--bool_grained', action='store_true',
        help='get rid of flying 1x1 pieces')
    parser.add_argument(
        '--rotate_x90', action='store_true', help='rotate around x-axis')
    parser.add_argument('--bool_color_by_residue', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument(
        '--bool_replace_small_bricks_with_larger_bricks',
        action='store_true')

    args = parser.parse_args()

    # http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements
    # van der Waals radii in 1/100 pm
    args.d_radii = {
        'H': 1.20,
        'C': 1.70,
        'N': 1.55,
        'O': 1.52,
        'P': 1.80,
        'S': 1.80,
        'CA': 2.31,
        'NA': 2.27,
        'CL': 1.75,
        'MG': 1.73,
        'FE': 2.00,
        'F': 1.47,
        'ZN': 1.39,
        'CO': 2.00,
        'SE': 1.90,
    }

    # http://en.wikipedia.org/wiki/CPK_coloring
    args.d_colors = {
        'H': 1,  # white
        'C': 26,  # black
        'N': 23,  # bright blue
        'O': 21,  # bright red
        'F': 28,  # (dark) green
        'CL': 28,  # (dark) green
        'BR': 28,
        'I': 28,
        'P': 106,  # bright orange
        'S': 24,  # bright yellow
        'FE': 106,  # bright orange
        'CO': 221,  # bright purple
        'ZN': 221,  # bright purple
        '119': 119,  # bright yellowish green (lime)
        '192': 192,  # reddish brown
        'SE': 106,
        '208': 208,  # grey
    }

    # http://en.wikipedia.org/wiki/CPK_coloring
    args.d_materialIDs = {
            'H': 1,  # white
            'C': 26,  # black
            'N': 23,  # bright blue
            'O': 21,  # bright red
            'F': 28,  # (dark) green
            'CL': 28,  # (dark) green
            'BR': 28,
            'I': 28,
            'P': 106,  # bright orange
            'S': 24,  # bright yellow
            'FE': 106,  # bright orange
            'CO': 221,  # bright purple
            'ZN': 221,  # bright purple
            '119': 119,  # bright yellowish green (lime)
            '192': 192,  # reddish brown
            'SE': 106,
            '208': 208,  # grey
            }

    # Which brick sizes are available for each color?
    with open('dic/d_brick_sizes_PAB.dic') as f:
        args.d_brick_sizes = eval(f.read())

    # What are the designIDs of the different brick sizes?
    args.d_designIDs = {
        'brick':{
            1:{
                1:3005,
                2:3004,
                3:3622,
                4:3010,
                6:3009,
                8:3008,
                10:6111,
                12:6112,
                16:2465,
                },
            2:{
                1:3004,
                2:3003,
                3:3002,
                4:3001,
                6:44237,
                8:3007,
                10:3006,
                },
        },
        'plate':{
            1:{
                1:3024,
                2:3023,
                3:3623,
                4:3710,
                6:3666,
                8:3460,
                10:4477,
                12:60479,
                },
            2:{
                1:3023,
                2:3022,
                3:3021,
                4:3020,
                6:3795,
                8:3034,
                10:3832,
                12:2445,
                16:4282,
                },
        },
    }

    return args


if __name__ == '__main__':
    main()
