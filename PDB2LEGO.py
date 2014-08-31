#!/usr/bin/env python 

## Tommy Carstensen

## todo: finish replacements of neigbouring pieces

## import built-ins
import os, math

class LEGO():

    def main(self,):

        pdb = self.pdb

        ## parse PDB coordinates
        d_coords_PDB = self.parse_pdb(pdb)

        ## convert Angstrom coordinates to brick coordinates
        d_layers, l_minmax_bricks, d_buried = self.coordAngstrom2coordBrick(
            pdb, d_coords_PDB,
            )
        del d_coords_PDB ## save some memory

        ## convert brick coordinates to bricks
        self.coordBrick2coordLDD(
            pdb, d_layers,l_minmax_bricks, d_buried,
            )
        del d_layers, l_minmax_bricks ## save some memory

        return


    def calculate_support_position(self,d_layers,l_minmax_bricks,):

        ##
        ## dictionary independent of elements
        ##
        d_3d = {}
        for z_brick in range(int(l_minmax_bricks[2][0]),int(l_minmax_bricks[2][1])+1):
            d_3d[z_brick] = {}
            for element in d_layers[z_brick].keys():
                for x in d_layers[z_brick][element]['x'].keys():
                    if not x in d_3d[z_brick].keys():
                        d_3d[z_brick][x] = []
                    d_3d[z_brick][x] += d_layers[z_brick][element]['x'][x]

        ##
        ## supporters
        ##
        d_support = {}
        for z_brick_above in range(
            int(l_minmax_bricks[2][0])+1,
            int(l_minmax_bricks[2][1])+1,
            ):
            d_support[z_brick_above] = {}

            ## loop over row/column
            for x in d_3d[z_brick_above].keys():
                l_above = d_3d[z_brick_above][x]

                set_above = set(l_above)
                for z_brick_below in range(int(l_minmax_bricks[2][0]),z_brick_above):
                    ## nothing in position below
                    if not x in d_3d[z_brick_below].keys():
                        d_support[z_brick_above][x] = set_above
                        continue
                    l_below = d_3d[z_brick_below][x]
                    set_above -= set(l_below)
                    ## everyting above is also below
                    if len(set_above) == 0:
                        break

                ## something above was not below
                ## use points for support from bottom layer
                if len(set_above) > 0:
                    d_support[z_brick_above][x] = set_above

        return d_support


    def parse_pdb(self,pdb):

        '''this function reads the coordinates of a PDB file and returns:
d_coords_PDB[record][chain][res_no][atom_name] = {
    'coord':coord,'element_size':element_size,
    'element_color':element_color,
    }
'''

        print 'reading and parsing %s.pdb' %(pdb)
        fd = open('pdb/%s.pdb' %(pdb),'r')
        lines = fd.readlines()
        fd.close()

        d_coords_PDB = {'ATOM':{},'HETATM':{},}

##        l_minmax_PDB = [
##            [1000.,-1000.,],[1000.,-1000.,],[1000.,-1000.,],
##            ]

        for line in lines:
            record = line[:6].strip()
            if record == 'ENDMDL':
                break
            if record not in ['ATOM','HETATM',]:
                continue
            res_name = line[17:20]
            chain = line[21]
            if res_name == 'HOH':
                continue
            element_size = element_color = line[76:78].strip()
            if self.exclude_hydrogens == True and element == 'H':
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if self.rotate_x90 == True:
                z = float(line[38:46])
                y = -float(line[46:54])
            coord = [x,y,z,]
            atom_name = line[12:16].strip()
            res_no = int(line[22:26])
            if not chain in d_coords_PDB[record].keys():
                d_coords_PDB[record][chain] = {}
            if not res_no in d_coords_PDB[record][chain].keys():
                d_coords_PDB[record][chain][res_no] = {}

            ##
            ## color by residue instead of by atom
            ##
            if res_name == 'HEM':
                element_color = 'CL'
            if self.bool_color_by_residue == True:
                d_elements = {
                    ' DC':'CL', ## green
                    ' DT':'192', ## reddish brown
                    ' DA':'119', ## lime
                    ' DG':'S', ## yellow
                    }
                element_color = d_elements[res_name]

            d_coords_PDB[record][chain][res_no][atom_name] = {
                'coord':coord,'element_size':element_size,
                'element_color':element_color,
                }

##            for i in range(3):
##                if coord[i] < l_minmax_PDB[i][0]:
##                    l_minmax_PDB[i][0] = coord[i]
##                if coord[i] > l_minmax_PDB[i][1]:
##                    l_minmax_PDB[i][1] = coord[i]
            
        return d_coords_PDB


    def gcd(self,a,b,):

        '''greatest common divisor using Euclid's algorithm'''

        while b > 1e-14:
            a, b = b, a % b

        return a


    def lcm(self,a,b,):

        '''least/lowest/smallest common multiple'''

        multiple = a*b // self.gcd(a,b,)

        return multiple


    def calc_brick_per_Angstrom(self,):

        ##
        ## find least common multiple of brick height and length
        ##
        height = 9.6 ## mm
        width = 8.0 ## mm
        multiplier = self.lcm(9.6,8.0)

##        ## slightly slower method for finding least common multiple
##        multiplier = 1
##        height = 9.6 ## mm
##        width = 8 ## mm
##        while True:
##            if multiplier*height%width == 0:
##                break
##            multiplier += 1
##        multiplier *= height

##        ## plates
##        bpa_x = bricks_per_angstrom_x = 4. ## 8mm width
##        bpa_y = bricks_per_angstrom_y = 4. ## 8mm width
##        bpa_z = bricks_per_angstrom_z = 10. ## 3.2mm height
        ## bricks (9.6mm height) - 48mm per Angstrom = 480,000,000:1        
        bpa_x = bricks_per_angstrom_x = multiplier/width ## 6, 8mm width
        bpa_y = bricks_per_angstrom_y = multiplier/width ## 6, 8mm width
        if self.brick_or_plate == 'brick':
            bpa_z = bricks_per_angstrom_z = multiplier/height ## 5, 9.6mm height (brick)
        elif self.brick_or_plate == 'plate':
            bpa_z = bricks_per_angstrom_z = 15. ## 3.2mm height (plate)

        ##
        ## additional scaling
        ##
        bpa_x /= self.size_scale
        bpa_y /= self.size_scale
        bpa_z /= self.size_scale

        ## list of scaling factors
        l_bpa = [bpa_x,bpa_y,bpa_z,]

        return l_bpa



    def coordAngstrom2coordBrick(self, pdb, d_coords_PDB):

        print 'convert PDB coordinates to brick coordinates'

        ##
        ## set i/o file path
        ##
        fn = 'txt/d_layers_%s_%s_%s' %(pdb,self.size_scale,self.brick_or_plate,)
        if self.bool_hollow == True:
            fn += '_hollow'
        if self.bool_color_by_residue == True:
            fn += '_colorbyresidue'
        if self.bool_buried_then_black == True:
            fn += '_buriedblack'
        if self.bool_grained == True:
            fn += '_flying1x1removal'
        fn += '.txt'

        ##
        ## set i/o file path
        ##
        fn_buried = 'txt/d_buried_%s_%s_%s' %(pdb,self.size_scale,self.brick_or_plate,)
        if self.bool_buried_then_black == True:
            fn_buried += '_buriedblack'
        if self.bool_grained == True:
            fn_buried += '_flying1x1removal'
        fn_buried += '.txt'

        if os.path.isfile(fn) and os.path.isfile(fn_buried):

            ##
            ## read
            ##

            print 'reading', fn
            fd = open(fn,'r')
            s = fd.read()
            fd.close()
            d_layers = eval(s)
            del s

            fd = open(fn_buried,'r')
            s = fd.read()
            fd.close()
            d_buried = eval(s)
            del s

        else:

            ##
            ## loop over PDB coordinates and check distance to surrounding grid points
            ##
            d_vicinal = self.find_vicinal(d_coords_PDB,)


            ##
            ##
            ##
            d_layers = self.append_bricks_by_atom_vicinity(d_vicinal,)


            ##
            ## remove pieces before determining what is buried and what is not
            ##
            d_layers = self.remove_bricks_outside(d_layers)

            ##
            ## hollow option
            ##
            if self.bool_hollow == True or self.bool_buried_then_black == True:

                d_layers, d_buried = self.remove_or_replace_bricks_inside(d_layers,)

            else:

                d_buried = self.find_buried(d_layers,)

            ##
            ## write output
            ##
            print 'writing', fn
            fd = open(fn,'w')
            fd.write(str(d_layers))
            fd.close()

            fd = open(fn_buried,'w')
            fd.write(str(d_buried))
            fd.close()

        ##
        ##
        ##
        l_minmax_bricks = self.get_dimensions(d_layers,)

        return d_layers, l_minmax_bricks, d_buried


    def remove_or_replace_bricks_inside(self,d_layers,):

        d_buried = self.find_buried(d_layers,)

        replace_atom = self.replace_atom
        
        ## either replace buried bricks with black
        ## or delete buried bricks
        for z_brick in d_buried.keys():

            if replace_atom not in d_layers[z_brick].keys():
                d_layers[z_brick][replace_atom] = {'x':{},'y':{},}
            
            for element in d_layers[z_brick].keys():

                ## don't replace carbon with carbon
                if self.bool_buried_then_black == True and element == replace_atom:
                    continue

                ## don't put carbon in a layer if it's not there in the first place
                if (
                    self.bool_buried_then_black == True
                    and
                    replace_atom not in d_layers[z_brick].keys()
                    ):
                    continue

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

                        ## either delete
                        if self.bool_hollow == True and count_buried == 27:
                            d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                            d_layers[z_brick][element]['y'][y_brick].remove(x_brick)

                        ## or replace
                        elif (
                            (
                                self.bool_buried_then_black == True
                                and
                                count_buried == 27
                                and self.size_scale <= 2
                                )
                            or
                            ## causes trouble in 2dau layer 69 and many other layers
                            ## works if run with graining!
                            ## option 1) allow diagonals to be empty
                            (
                                self.bool_buried_then_black == True
                                and
                                count_buried_excl_diagonals == 7
                                and
                                self.size_scale >= 3)
                            ## option 2) do not allow diagons to be empty
##                                    (
##                                        self.bool_buried_then_black == True
##                                        and
##                                        count_buried == 27
##                                        and self.size_scale >= 3
##                                        )
                            ):

                            ## remove the old element
                            d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                            d_layers[z_brick][element]['y'][y_brick].remove(x_brick)

                            ## replace with the new element
                            if not x_brick in d_layers[z_brick][replace_atom]['x'].keys():
                                d_layers[z_brick][replace_atom]['x'][x_brick] = []
                            d_layers[z_brick][replace_atom]['x'][x_brick] += [y_brick]

                            if not y_brick in d_layers[z_brick][replace_atom]['y'].keys():
                                d_layers[z_brick][replace_atom]['y'][y_brick] = []
                            d_layers[z_brick][replace_atom]['y'][y_brick] += [x_brick]

        return d_layers, d_buried


    def get_dimensions(self,d_layers,):

        if self.verbose == True:
            print 'determining brick dimensions'
        l_minmax_bricks = [[1000.,-1000.,],[1000.,-1000.,],[1000.,-1000.,],]
        for z_brick in d_layers.keys():
            for element in d_layers[z_brick].keys():
                for x_brick in d_layers[z_brick][element]['x'].keys():
                    for y_brick in d_layers[z_brick][element]['x'][x_brick]:
                        coord_brick = [x_brick,y_brick,z_brick,]
                        for i in range(3):
                            if coord_brick[i] < l_minmax_bricks[i][0]:
                                l_minmax_bricks[i][0] = coord_brick[i]
                            if coord_brick[i] > l_minmax_bricks[i][1]:
                                l_minmax_bricks[i][1] = coord_brick[i]

        print 'Brick dimensions (bricks)', l_minmax_bricks
        print 'Brick dimensions (model,cm)',

        print 'width', 0.8*(l_minmax_bricks[0][1]-l_minmax_bricks[0][0]),
        print 'width', 0.8*(l_minmax_bricks[1][1]-l_minmax_bricks[1][0]),
        if self.brick_or_plate == 'plate':
            print 'height', 0.32*(l_minmax_bricks[2][1]-l_minmax_bricks[2][0])
        elif self.brick_or_plate == 'brick':
            print 'height', 0.96*(l_minmax_bricks[2][1]-l_minmax_bricks[2][0])

        return l_minmax_bricks
    

    def remove_bricks_outside(self,d_layers,):

        print 'removing "outside" bricks before determining which bricks are buried'
        if self.bool_grained == True:

            l_removal = []

            ##
            ## loop over sorted layers to remove bricks from the bottom and up
            ##
            l_zbricks = d_layers.keys()
            l_zbricks.sort()
            for z_brick in l_zbricks:
                ## don't remove the bottom brick no matter what
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
                            for z_add in [-1,0,1,]:
                                z2 = z_brick+z_add
                                if not z_brick+z_add in d_layers.keys():
                                    continue
                                ## loop over all neigbouring element types
                                for element2 in d_layers[z_brick+z_add].keys():
                                    for x_add in [-1,0,1,]:
                                        x2 = x_brick+x_add
                                        if not (
                                            x_brick+x_add
                                            in
                                            d_layers[z_brick+z_add][element2]['x'].keys()
                                            ):
                                            continue
                                        for y_add in [-1,0,1,]:
                                            y2 = y_brick+y_add
                                            if (
                                                y_brick+y_add
                                                in
                                                d_layers[z_brick+z_add][element2]['x'][x_brick+x_add]
                                                ):
                                                count_buried += 1
                                                if (
                                                    (y_add == 0 and z_add == 0)
                                                    or
                                                    (x_add == 0 and z_add == 0)
                                                    or
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

                            ## "get rid of 1x1 hanging pieces"
                            ## nothing *above or below* brick
                            ## and only 1 brick in the plane around the brick
                            if count_above_below == 0 and count_plane == 1+1:
                                l_removal += [[z_brick,element,x_brick,y_brick,]]
                            ## "don't add 1x1s to new layers when building from the bottom and up"
                            ## nothing in the plane around the brick
                            ## and nothing *below* the brick
                            elif count_below == 0 and count_plane == 1:
                                d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                                d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
                                print 'remove2', z_brick, element
                            ## "get rid of pieces sticking out otherwise troubling the use of 2 knob bricks"
                            ## nothing above or below brick
                            ## (*after* removal of any bricks below;
                            ## thus sorted z_brick values necessary
                            ## and immediate removal in previous layer happened)
                            elif count_above_below == 0:
                                d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                                d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
                                print 'remove3', z_brick, element

                    ## end of loop over elements

            for z_brick, element, x_brick, y_brick in l_removal:
                d_layers[z_brick][element]['x'][x_brick].remove(y_brick)
                d_layers[z_brick][element]['y'][y_brick].remove(x_brick)
                print 'remove1', z_brick, element

        return d_layers


    def find_vicinal(self,d_coords_PDB,):

        print 'find grid points vicinal to atom coordinates'

        l_bpa = self.calc_brick_per_Angstrom()

        d_vicinal = {}

        for record in d_coords_PDB.keys():
            for chain in d_coords_PDB[record].keys():
                for res_no in d_coords_PDB[record][chain].keys():
                    print chain, res_no
                    for atom_name in d_coords_PDB[record][chain][res_no].keys():
                        if record == 'HETATM':
                            print 'exclude', record, chain, res_no
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
                        ## convert Angstrom coord to grid coord
                        ## loop over vicinal grid coords (cube)
                        grid_brick = []
                        for i in range(3):
                            grid_brick += [[
                                int(math.floor(
                                    (coord_PDB[i]-self.d_radii[element_size])*l_bpa[i]
                                    )),
                                int(math.ceil(
                                    (coord_PDB[i]+self.d_radii[element_size])*l_bpa[i]
                                    )),
                                ]]
                        sq_radius = self.d_radii[element_size]**2
                        position = 0
                        positionok = 0

                        ##
                        ## loop over grid bricks and check distance from center of atom
                        ##
                        for x_brick in range(grid_brick[0][0],grid_brick[0][1]+1,):
                            x_LEGO_Angstrom = x_brick/l_bpa[0]
                            for y_brick in range(grid_brick[1][0],grid_brick[1][1]+1,):
                                y_LEGO_Angstrom = y_brick/l_bpa[1]
                                for z_brick in range(grid_brick[2][0],grid_brick[2][1]+1,):
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
                                    ## grid point too distant from center of sphere
                                    if sq_dist > sq_radius:
                                        continue

                                    if not x_brick in d_vicinal.keys():
                                        d_vicinal[x_brick] = {}
                                    if not y_brick in d_vicinal[x_brick].keys():
                                        d_vicinal[x_brick][y_brick] = {}
                                    if not z_brick in d_vicinal[x_brick][y_brick].keys():
                                        d_vicinal[x_brick][y_brick][z_brick] = {}
                                    d_vicinal[x_brick][y_brick][z_brick][sq_dist] = element_color


        return d_vicinal


    def append_bricks_by_atom_vicinity(self,d_vicinal,):

        print 'appending bricks according to atom nearest to grid point'
        d_layers = {}
        l_x_bricks = d_vicinal.keys()
        l_x_bricks.sort()
        for x_brick in l_x_bricks:
            l_y_bricks = d_vicinal[x_brick].keys()
            l_y_bricks.sort()
            for y_brick in l_y_bricks:
                for z_brick in d_vicinal[x_brick][y_brick].keys():
                    dist_min = min(d_vicinal[x_brick][y_brick][z_brick].keys())
                    element = d_vicinal[x_brick][y_brick][z_brick][dist_min]

                    ##
                    ## append coordinate to dictionary of coordinates
                    ##
                    if not z_brick in d_layers.keys():
                        d_layers[z_brick] = {}
                    if not element in d_layers[z_brick].keys():
                        d_layers[z_brick][element] = {'x':{},'y':{},}
                    if not x_brick in d_layers[z_brick][element]['x'].keys():
                        d_layers[z_brick][element]['x'][x_brick] = []
                    if not y_brick in d_layers[z_brick][element]['y'].keys():
                        d_layers[z_brick][element]['y'][y_brick] = []
        ##            d_layers[z_brick][element] += [[x_brick,y_brick,]]
                    d_layers[z_brick][element]['x'][x_brick] += [y_brick]
                    d_layers[z_brick][element]['y'][y_brick] += [x_brick]

        return d_layers


    def find_buried(self,d_layers,):

        print 'preparing to remove "buried" bricks \
after removing "isolated" bricks on the outside'
        d_buried = {}
        for z_brick in d_layers.keys():
            d_buried[z_brick] = {}
            for element in d_layers[z_brick].keys():
                for x_brick in d_layers[z_brick][element]['x'].keys():
                    if not x_brick in d_buried[z_brick].keys():
                        d_buried[z_brick][x_brick] = {}
                    for y_brick in d_layers[z_brick][element]['x'][x_brick]:
                        count_buried = 0
                        count_buried_excl_diagonals = 0
                        for z_add in [-1,0,1,]:
                            z2 = z_brick+z_add
                            if not z_brick+z_add in d_layers.keys():
                                continue
                            ## loop over all neigbouring element types
                            for element2 in d_layers[z_brick+z_add].keys():
                                for x_add in [-1,0,1,]:
                                    x2 = x_brick+x_add
                                    if not x_brick+x_add in d_layers[z_brick+z_add][element2]['x'].keys():
                                        continue
                                    for y_add in [-1,0,1,]:
                                        y2 = y_brick+y_add
                                        if (
                                            y_brick+y_add
                                            in
                                            d_layers[z_brick+z_add][element2]['x'][x_brick+x_add]
                                            ):
                                            count_buried += 1
                                            if (
                                                (y_add == 0 and z_add == 0)
                                                or
                                                (x_add == 0 and z_add == 0)
                                                or
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
                            'all':count_buried,'excl_diagonals':count_buried_excl_diagonals,
                            }

        return d_buried


    def check_if_buried(self,d_layers,z_brick,xy,consecutive,k1,k2,):

        l_below_all_elements  = self.check_bricks_vicinal_all_elements(
            d_layers,z_brick-1,k1,xy,consecutive,
            )
        l_above_all_elements  = self.check_bricks_vicinal_all_elements(
            d_layers,z_brick+1,k1,xy,consecutive,
            )
        l_behind_all_elements = self.check_bricks_vicinal_all_elements(
            d_layers,z_brick,k1,xy-1,consecutive,
            )
        l_front_all_elements  = self.check_bricks_vicinal_all_elements(
            d_layers,z_brick,k1,xy+1,consecutive,
            )
##        l_left_all_elements   = self.check_bricks_vicinal_all_elements(
##            d_layers,z_brick,k2,consecutive[0]-1,consecutive,
##            )
##        l_right_all_elements  = self.check_bricks_vicinal_all_elements(
##            d_layers,z_brick,k2,consecutive[0]+1,consecutive,
##            )
        l_yx = range(consecutive[0],consecutive[0]+consecutive[1])

        if (
            ## no bricks below
            len(set(l_yx)&set(l_below_all_elements)) == 0
            and
            ## no bricks above
            len(set(l_yx)&set(l_above_all_elements)) == 0
            ):
            print consecutive
            print l_below_all_elements
            print l_above_all_elements
            print l_yx
            stop_nothing_to_hang_on_to

        ## above/below
        set_below_and_above = set(l_below_all_elements) & set(l_above_all_elements)
        set_exposed_vertical = set(l_yx)-set_below_and_above

        ## behind/front
        set_behind_and_front = set(l_behind_all_elements) & set(l_front_all_elements)
        set_exposed_horizontal1 = set(l_yx)-set_behind_and_front

        ## left/right
        None

        set_exposed = set_exposed_horizontal1 | set_exposed_vertical
        set_buried = set(l_yx)-set_exposed

        l_exposed = list(set_exposed)
        l_buried = list(set_buried)

        l_exposed.sort()
        l_buried.sort()

        return l_exposed, l_buried


    def coordBrick2coordLDD(
        self, pdb, d_layers, l_minmax_bricks, d_buried,
        ):

        '''convert 1x1x1 brick coordinates to bricks'''

        ##
        ##
        ##
        d_bricks_main = self.coords_to_lxfml_and_gnuplot(
            d_layers, l_minmax_bricks, d_buried,
            )

        ## clean up
        for element in self.d_radii.keys()+['buried']:
            if os.path.isfile('LEGO_%s_prev.gnuplotdata' %(element)):
                os.remove('LEGO_%s_prev.gnuplotdata' %(element))
        if os.path.isfile('LEGO.gnuplotsettings'):
            os.remove('LEGO.gnuplotsettings')

        d_bricks_main = self.check_vicinal_small_bricks(d_bricks_main)

        ## write shopping list
        ## after replacing small bricks with large bricks
        ## and before recoloring
        self.write_shopping_list_and_calculate_price(pdb,d_bricks_main,)

        ## convert dic of bricks to lines of bricks
        ## after replacing small bricks with large bricks and after recoloring
        lines_lxfml_body = self.dic2lines(d_bricks_main)

        self.write_lxfml(pdb,lines_lxfml_body,)

        if self.bool_split_into_layers == True:

            self.split_into_layers(lines_lxfml_body)

        return


    def coords_to_lxfml_and_gnuplot(
        self,
        d_layers, l_minmax_bricks, d_buried,
        ):

        print int(l_minmax_bricks[2][1])-int(l_minmax_bricks[2][0])+1, 'layers'

        self.refID = 0
        self.volume = 0
        d_bricks_main = {}
        lines_lxfml_body = []

        ## loop over layers
        for z_brick in range(
            int(l_minmax_bricks[2][0]),
            int(l_minmax_bricks[2][1])+1,
            ):

            layer = z_brick-l_minmax_bricks[2][0]

            ## no bricks in layer
            if not z_brick in d_layers.keys():
                continue

            if not z_brick in d_layers.keys():
                print 'no bricks in layer', layer
                stop
                continue

            if self.layer_range:
                if layer not in self.layer_range: continue

            if self.verbose == True:
                print 'layer', layer

            ## set keys
            ## odd and even layers rotated relative to each other
            if layer % 2 == 0:
                k1 = k = 'y'
                k2 = 'x'
            else:
                k1 = k = 'x'
                k2 = 'y'

            ##
            ## support positions
            ##
            lines_buried_gnuplot = []
            lines_support_gnuplot = []
            if self.bool_gnuplot == True:
                d_support = self.calculate_support_position(d_layers,l_minmax_bricks,)
            else:
                d_support = {}

            ## loop over element types
            for element in d_layers[z_brick].keys():
                materialID = self.d_materialIDs[element]
                lines_element_gnuplot = []

                l_consecutives_skip = []
                l_coords_skip = []
                l_xys = d_layers[z_brick][element][k].keys()
                l_xys.sort()
                ## loop over sorted rows/columns in layer
                for xy in l_xys:

                    l_coords = d_layers[z_brick][element][k][xy]
                    l_consecutives = self.identify_consecutive_stretches(l_coords)
                    l_exposed = []

                    ##
                    ## q: bricks of same materialID/element present in
                    ## the next (xy+1) or previous (xy-1) or next next (xy+2) row?
                    ## could connect to it, in case it's hanging/flying
                    ##
                    if xy+1 in d_layers[z_brick][element][k].keys():
                        l_coords_next_same_element = d_layers[z_brick][element][k][xy+1]
                        l_consecutives_next_same_element = self.identify_consecutive_stretches(
                            l_coords_next_same_element,
                            )
                    else:
                        l_consecutives_next_same_element = []

                    if xy+2 in d_layers[z_brick][element][k].keys():
                        l_coords_nextnext_same_element = d_layers[z_brick][element][k][xy+2]
                        l_consecutives_nextnext_same_element = self.identify_consecutive_stretches(
                            l_coords_nextnext_same_element,
                            )
                    else:
                        l_consecutives_nextnext_same_element = []
                        
                    ##
                    ##
                    ## split consecutive stretches into which bricks?
                    ##
                    ##
                    for consecutive in l_consecutives:

                        ##
                        ## was the previous brick a double width brick?
                        ##
                        if [xy,consecutive,] in l_consecutives_skip:
                            print l_consecutives_skip
                            stop1
                            continue

                        ##
                        ## use double or single width for current *consecutive* brick?
                        ## and decide whether to check if nextnext (double width)
                        ## or next (single width) hanging/flying
                        ##
                        if (
                            self.bool_use_double_width == True
                            and
                            consecutive in l_consecutives_next_same_element
                            ):
                            width = 2
                            l_consecutives_skip += [[xy+1,consecutive,]]
                            for yx in range(consecutive[0],consecutive[0]+consecutive[1],):
                                l_coords_skip += [[xy+1,yx,]]
                            ## use nextnext for consecutive_next_same_element loop
                            l = l_consecutives_nextnext_same_element
                        else:
                            width = 1
                            ## use next for consecutive_next_same_element loop
                            l = l_consecutives_next_same_element


                        l_consecutives_split = [consecutive]
                        bool_next_rotated = False

                        ##
                        ## accomodate for rotated 2x1 or 3x1 brick in current row,
                        ## if next brick is a 1x1 brick not attached above or below
                        ## thus necessary and not just to save bricks, which is done in a later step...
                        ## thus not necessary if 1x1 brick connected above or below!!!
                        ## splits l_consecutives_split into smaller fragments
                        ##
                        if consecutive[1] != 1:
                            
                            for consecutive_next_same_element in l:

                                if (
                                    ## next brick in layer is a 1x1 brick
                                    consecutive_next_same_element[1] == 1
                                    and
                                    consecutive_next_same_element[0] in range(
                                        l_consecutives_split[-1][0],
                                        l_consecutives_split[-1][0]+l_consecutives_split[-1][1],
                                        )
                                    ):

                                    ## anything above or below next brick of any element?
                                    l_below_all_elements = self.check_bricks_vicinal_all_elements(
                                        d_layers,z_brick-1,k,xy+width,consecutive_next_same_element,
                                        )
                                    l_above_all_elements = self.check_bricks_vicinal_all_elements(
                                        d_layers,z_brick+1,k,xy+width,consecutive_next_same_element,
                                        )
                                    ## anything above or below next brick of any element?
                                    if len(l_below_all_elements) == 0 and len(l_above_all_elements) == 0:
                                        ## split if anything above or below next brick of any element
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
                                            for xy_add in range(1,width+1):
                                                l_consecutives_skip += [[xy+xy_add,consecutive,]]
                                        elif pos_break == pos2:
                                            brick1 = [
                                                pos1, pos_break-pos1, 'not rotated', width,
                                                ]
                                            brick2 = [
                                                pos_break, 1, 'rotated', width,
                                                ]
                                            brick3 = None
                                            for xy_add in range(1,width+1):
                                                l_consecutives_skip += [[xy+xy_add,consecutive,]]
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
                                            for xy_add in range(1,width+1):
                                                l_consecutives_skip += [[xy+xy_add,consecutive,]]
                                        ## even newer method
                                        if brick3:
                                            l_consecutives_split = (
                                                l_consecutives_split[:-1]
                                                +
                                                [brick1, brick2, brick3,]
                                                )
                                        else:
                                            l_consecutives_split = (
                                                l_consecutives_split[:-1]
                                                +
                                                [brick1, brick2,]
                                                )

                                    bool_next_rotated = True

                                    continue

                        ## dont skip next line of bricks,
                        ## because one or more positions of the current consecutive stretch were rotated
                        if bool_next_rotated == True and width == 2:
                            l_consecutives_skip.remove([xy+1,consecutive,])
                            width = 1

                        ##
                        ## accomodate for rotated 2x1 or 3x1 brick in next row,
                        ## if curent brick is a 1x1 brick not attached above or below
                        ##
                        elif consecutive[1] == 1:
                            ## accomodation happens in function split_consecutive_length_into_bricks
                            ## so do nothing
                            pass
                            

                        ##
                        ##
                        ## split consecutive *stretches* into *available* bricks
                        ##
                        ##
                        for consecutive_split in l_consecutives_split:

                            bool_prev_skipped = False
                            ## previous brick was skipped to acomodate for brick previous to it
                            if consecutive_split[1] == 1 and [xy-1,consecutive_split[0],] in l_coords_skip:
                                bool_prev_skipped = True

                            d_bricks, d_layers = self.split_consecutive_length_into_bricks(
                                consecutive_split, materialID, width,
                                d_layers, z_brick, element, k, xy,
                                bool_next_rotated, bool_prev_skipped,
                                )

                            materialID = self.d_materialIDs[element]

                            ##
                            ## append bricks to lxfml
                            ##
                            (
                                lines_lxfml_body, d_bricks_main,
                                ) = self.append_bricks(
                                d_bricks, materialID, width, k,
                                ## position
                                layer, xy,
                                ## output (append)
                                lines_lxfml_body,
                                d_bricks_main,
                                )

                            ## end of loop over split consecutives

                        ## end of loop over consecutives

                    ## end of loop over rows/columns

                ## append bricks to gnuplot
                lines_support_gnuplot = self.append_gnuplot(
                    d_layers,z_brick,element,d_buried,d_support,
                    lines_element_gnuplot,
                    lines_buried_gnuplot,
                    lines_support_gnuplot,
                    )

                ## write gnuplotdata for element
                if self.bool_gnuplot == True:
                    fd = open('LEGO_%s.gnuplotdata' %(element),'w')
                    fd.writelines(lines_element_gnuplot)
                    fd.close()

                ## end of loop over elements

            if self.bool_gnuplot == True and self.bool_buried_then_black == True:
                fd = open('LEGO_buried.gnuplotdata','w')
                fd.writelines(lines_buried_gnuplot)
                fd.close()

            if self.bool_gnuplot == True:
                fd = open('LEGO_support.gnuplotdata','w')
                fd.writelines(lines_support_gnuplot)
                fd.close()

            ## 1x1 plot can be safely done before brick replacement etc
            if self.bool_gnuplot == True:
                print 'plotting layer', layer
                self.plot(pdb,layer,l_minmax_bricks,)

            ## end of loop over layers

        if self.verbose == True:
            self.write_lxfml(self.pdb+'preSmallToLargeBricks',lines_lxfml_body,)

        return d_bricks_main


    def append_gnuplot(
        self,d_layers,z_brick,element,d_buried,d_support,
        lines_element_gnuplot,
        lines_buried_gnuplot,
        lines_support_gnuplot,
        ):

        for x in d_layers[z_brick][element]['x'].keys():
            for z in d_layers[z_brick][element]['x'][x]:

                bool_buried = False

                if not x in d_buried[z_brick].keys():
                    pass
                elif not z in d_buried[z_brick][x].keys():
                    pass
                elif d_buried[z_brick][x][z]['excl_diagonals'] == 7:
                    bool_buried = True

                if self.bool_buried_then_black == True and bool_buried == True:
                    lines_buried_gnuplot += ['%s %s\n' %(x,z,)]
                else:
                    lines_element_gnuplot += ['%s %s\n' %(x,z,)]

                if z_brick in d_support.keys():
                    if x in d_support[z_brick].keys():
                        for z in d_support[z_brick][x]:
                            lines_support_gnuplot += ['%s %s\n' %(x,z,)]

        return lines_element_gnuplot, lines_buried_gnuplot, lines_support_gnuplot


    def split_into_layers(self,lines_lxfml_body):

        ## this should happen *after* brick replacement!!!

        d_ty = {}
        line_index1 = 0
        for line in lines_lxfml_body:
            k = 'ty'
            index1 = line.index(k)+len(k)+2
            index2 = index1+line[index1:].index('"')
            ty = round(float(line[index1:index2]),2)
            if not ty in d_ty.keys():
                d_ty[ty] = []
            d_ty[ty] += [line]
        l_ty = d_ty.keys()
        l_ty.sort()
        ## loop layers
        for ty in l_ty[:-1]:
            lines = list(d_ty[ty])
            lines += list(d_ty[round(ty+0.96,2)])
            refID = 1
            ## loop lines in layer
            d_htm = {}
            d = {'refID':None,'ty':None,'designID':None,'materialID':None,}
            for i_line in range(len(lines)):
                line = lines[i_line]
                for k in d.keys():
                    index1 = line.index(k)+len(k)+2
                    index2 = index1+line[index1:].index('"')
                    v = line[index1:index2]
                    d[k] = [index1,index2,v,]
                line = '%s%s%s' %(line[:d['refID'][0]], refID, line[d['refID'][1]:])
                lines[i_line] = line

                ## only display bricks for the top layer
                if i_line > len(d_ty[ty])-1:
                    materialID = int(d['materialID'][2])
                    designID = int(d['designID'][2])
                    if not designID in d_htm.keys():
                        d_htm[designID] = {}
                    if not materialID in d_htm[designID].keys():
                        d_htm[designID][materialID] = 0
                    d_htm[designID][materialID] += 1
                    
                refID += 1

##                if not os.path.isdir('lxfml/layers/'+pdb):
##                    os.mkdir('lxfml/layers/'+pdb)
            self.write_lxfml(
                'layers/'+self.pdb+'/'+self.pdb+'_capa%s' %(int(ty/0.96)+1+1),
                lines,
                )
            self.write_htm(
                'layers/'+self.pdb+'_capa%s' %(int(ty/0.96)+1+1),
                d_htm,
                )

        return


    def write_htm(self,pdb,d_htm,):

        lines_htm = []

        fd = open('templates/header.htm','r')
        lines_htm += fd.readlines()
        fd.close()

        ##
        ## get dimensions of each brick
        ##
        d_dimensions = {}
        for width in self.d_designIDs[self.brick_or_plate].keys():
            for length in self.d_designIDs[self.brick_or_plate][width].keys():
                designID = self.d_designIDs[self.brick_or_plate][width][length]
                d_dimensions[designID] = '%s x %s' %(width,length,)
        d_dimensions[2357] = '1x2x2 corner brick'

        ## htm body
        lines_htm += ['<div id="modelParts">\n']
        lines_htm += ['<ul>\n']
        for designID in d_htm.keys():
            for materialID in d_htm[designID].keys():
                count = d_htm[designID][materialID]
                k = [designID,materialID,]
                lines_htm += [
                    '<li class="stepText"><p>%i x (%s)</p></li>\n' %(
                        count,d_dimensions[k[0]],
                        )
                    ]
                lines_htm += [
                    '<li><img height="64" width="64" src="%s_%s.png" /></li>\n' %(
                        k[0],k[1]
                        )
                    ]

        ##
        ## htm tail
        ##
        lines_htm += ['</ul>\n']            
        lines_htm += ['</div>\n'] ## modelParts
        lines_htm += ['</body>\n']            
        lines_htm += ['</html>\n']

        fd = open('htm/%s.htm' %(pdb),'w')
        fd.writelines(lines_htm)
        fd.close()

        return


    def write_shopping_list_and_calculate_price(self,pdb,d_bricks_main,):

        ##
        ## get contents of LEGO item 6177 bucket
        ##
        fd = open('dic/d_bucket6177.dic','r')
        s = fd.read()
        fd.close()
        d_bucket6177 = eval(s)

        ##
        ## get dimensions of each brick
        ##
        d_dimensions = {}
        for width in self.d_designIDs[self.brick_or_plate].keys():
            for length in self.d_designIDs[self.brick_or_plate][width].keys():
                designID = self.d_designIDs[self.brick_or_plate][width][length]
                d_dimensions[designID] = '%s x %s' %(width,length,)
        d_dimensions[2357] = '1x2x2 corner brick'

        ##
        ## parse bricks for shopping list
        ##
        d_shopping_list = {}
        for layer in d_bricks_main.keys():
            for tx in d_bricks_main[layer]['tx'].keys():
                for tz in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                    designID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID']
                    materialID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID']
                    if not designID in d_shopping_list.keys():
                        d_shopping_list[designID] = {}
                    if not materialID in d_shopping_list[designID].keys():
                        d_shopping_list[designID][materialID] = 0
                    d_shopping_list[designID][materialID] += 1

        ##
        ## calculate price and weight before subtraction of inventory
        ##
        price = 0
        mass = weight = 0
        for designID in d_shopping_list.keys():
            for materialID in d_shopping_list[designID].keys():
                v = d_shopping_list[designID][materialID]
                k = [designID,materialID,]
                price += v*self.d_prices_PAB[designID]
                weight += v*self.d_weights[designID]
        print 'price before subtraction of inventory', self.currency, price

        ## subtract inventory
        if self.bool_subtract_inventory == True:
            for designID in d_shopping_list.keys():
                for materialID in d_shopping_list[designID].keys():
                    if not materialID in self.inventory.keys():
                        continue
                    if not designID in self.inventory[materialID].keys():
                        continue
                    d_shopping_list[designID][materialID] = max(
                        0,
                        ## required-inventory
                        d_shopping_list[designID][materialID]-self.inventory[materialID][designID],
                        )

        ##
        ## calculate number of buckets
        ## to replace individually picked bricks (lowest price)
        ##
        print
        print '%-3s %-7s %-7s %-7s' %('n', 'total', 'PAB', '6177')
        for n_buckets in range(1,201):

            if self.bool_buy_buckets == False:
                n_buckets -= 1
                break

            if n_buckets == 200:
                n_buckets -= 1
                break

            d_shopping_list_next = {}

            ## calculate number of bricks in model not in buckets
            for designID in d_shopping_list.keys():
                d_shopping_list_next[designID] = {}
                for materialID in d_shopping_list[designID].keys():
                    d_shopping_list_next[designID][materialID] = (
                        d_shopping_list[designID][materialID]
                        )

                    if not materialID in d_bucket6177.keys():
                        continue
                    if not designID in d_bucket6177[materialID].keys():
                        continue

                    count_model = d_shopping_list[designID][materialID]
                    count_bucket = d_bucket6177[materialID][designID]
                    d_shopping_list_next[designID][materialID] = max(0,count_model-count_bucket)

            ## calculate price
            price_bricks = 0
            for designID in d_shopping_list.keys():
                for materialID in d_shopping_list[designID].keys():
                    v = d_shopping_list_next[designID][materialID]
                    k = [designID,materialID,]
                    price_bricks += v*self.d_prices_PAB[designID]
            ## add price of buckets (EUR 19.95 per bucket)
            if self.currency == 'EUR':
                price_buckets = 19.95*n_buckets
            elif self.currency == 'USD':
                price_buckets = 29.99*n_buckets
            price_next = price_bricks+price_buckets
            print '%3i %7.2f %7.2f %7.2f' %(
                n_buckets, price_next, price_bricks, price_buckets,
                )

            ## continue buing buckets, while it's cheaper than PAB
            if price_next < price:
                price = price_next
                d_shopping_list = d_shopping_list_next
            ## continue buying buckets and end up with excess bricks
##            elif round(price_next-price,2) != 19.95:
            ## don't break if price jump less than 10USD when buying box
            elif price_next-price < 10.:
                price = price_next
                d_shopping_list = d_shopping_list_next
            else:
                n_buckets -= 1
                break
        print n_buckets, price

        ##
        ## htm head
        ##
        lines_htm = []

        fd = open('templates/header.htm','r')
        lines_htm += fd.readlines()
        fd.close()

        lines_htm += ['<div id="modelParts">\n']
        lines_htm += ['<ul>\n']

        ## box 6177
        lines_htm += [
            '<li class="stepText"><p>%i x (650 pcs, set 6177)</p></li>\n' %(
                n_buckets,
                )
            ]
        lines_htm += [
            '<li><img height="64" width="64" src="bucket6177.png" /></li>\n'
            ]

        ##
        ## count bricks (of each type and total)
        ##
        l_counts = []
        count_total = 0
        for designID in d_shopping_list.keys():
            for materialID in d_shopping_list[designID].keys():
                v = d_shopping_list[designID][materialID]
                k = [designID,materialID,]
                l_counts += [[v,k,]]
                count_total += v

        ##
        ## htm body
        ##
        l_counts.sort()
        l_counts.reverse()

        ## sort by count
        for count,k in l_counts:

            if count == 0:
                continue

            lines_htm += [
                '<li class="stepText"><p>%i x (%s)</p></li>\n' %(
                    count,d_dimensions[k[0]],
                    )
                ]
            lines_htm += [
                '<li><img height="64" width="64" src="%s_%s.png" /></li>\n' %(
                    k[0],k[1]
                    )
                ]

        ##
        ## htm tail
        ##
        lines_htm += ['</ul>\n']            
        lines_htm += ['</div>\n'] ## modelParts
        lines_htm += ['</body>\n']            
        lines_htm += ['</html>\n']            

        ##
        ## write the htm shopping list
        ##

        fn = 'htm/%s_%s_%s_%s' %(
            pdb,self.size_scale,self.brick_or_plate,self.currency,
            )
        if self.bool_hollow == True:
            fn += '_hollow'
        if self.bool_buy_buckets == True:
            fn += '_bucket'
        if self.bool_buried_then_black == True:
            fn += '_buriedblack'
        fn += '.htm'
        fd = open(fn,'w')
        fd.writelines(lines_htm)
        fd.close()

        ##
        ## bricklink wish list xml
        ##
        WantedListID = 153019 ## move to init...
        suffix = '%s_%s_%s_%s_%s' %(pdb,self.size_scale,self.brick_or_plate,self.currency,WantedListID,)
        if self.bool_hollow == True:
            suffix += '_hollow'
        if self.bool_buy_buckets == True:
            suffix += '_bucket'

        ## read prices
        fd = open('dic/d_prices_PAB_%s.dic' %('USD'),'r')
        s = fd.read()
        fd.close()
        d_prices_PAB_USD = eval(s)

        self.write_BL_wishlist(suffix,WantedListID,d_shopping_list, d_prices_PAB_USD,)

        print count_total, 'bricks needed to build the model',
        print '(*after* replacement of small bricks with large bricks',
        print 'and *after* replacement of bricks with buckets)'
        print 'price in EUR if PAB (Pick A Brick) =', price
        print 'Brick volume (1x1x1) =', self.volume
        print 'Weight (kg) =', round(weight/1000.,3)

        return


    def write_BL_wishlist(self,suffix,WantedListID,d_shopping_list,d_prices_PAB_USD = None,):

        ## peeron
        d_colors_pab2bricklink = {
            1:1, ## White / White
            5:2, ## Brick yellow / Tan
            6:38, ## Light green / Light Green
            12:29, ## Light Orange Brown / Earth Orange
            18:28, ## Nougat / Flesh
            21:5, ## Bright red / Red
            23:7, ## Bright Blue / Blue
            24:3, ## Bright yellow / Yellow
            25:8, ## Earth Orange / Brown
            26:11, ## Black / Black
            28:6, ## Dark green / Green
            29:37, ## Medium green / Medium Green
            36:96, ## Light Yellowish Orange / Very Light Orange
            37:36, ## Bright Green / Bright Green
            38:68, ## Dark Orange / Dark Orange
            39:44, ## Light bluish violet / Light Violet
            40:12, ## Transparent / Trans-Clear
            41:17, ## Tr. Red / Trans-Red
            42:15, ## Tr. Lg blue / Trans-Light Blue
            43:14, ## Tr. Blue / Trans-Dark Blue
            44:19, ## Tr. Yellow / Trans-Yellow
            45:62, ## Light blue / Light Blue
            47:18, ## Tr. Flu. Reddish orange / Trans-Neon Orange
            48:20, ## Tr. Green / Trans-Green
            49:16, ## Tr. Flu. Green / Trans-Neon Green
            50:46, ## Phosp. White / Glow In Dark Opaque
            100:26, ## Light red / Light Salmon
            102:42, ## Medium blue / Medium Blue
            104:24, ## Bright violet / Purple
            105:31, ## Br. yellowish orange / Medium Orange
            106:4, ## Bright orange / Orange
            107:39, ## Bright bluish green / Dark Turquoise
            110:43, ## Bright Bluish Violet / Violet
            111:13, ## Tr. Brown / Trans-Black
            112:97, ## Medium Bluish Violet / Blue-Violet
            113:50, ## Tr. Medi. reddish violet / Trans-Dark Pink
            115:76, ## Med. yellowish green / Medium Lime
            116:40, ## Med. bluish green / Light Turquoise
            118:41, ## Light bluish green / Aqua
            119:34, ## Br. yellowish green / Lime
            120:35, ## Lig. yellowish green / Light Lime
            124:71, ## Bright reddish violet / Magenta
            125:32, ## Light Orange / Light Orange
            126:51, ## Tr. Bright bluish violet / Trans-Purple
            127:61, ## Gold / Pearl Light Gold
##            128:, ## Dark Nougat / 
            131:66, ## Silver / Pearl Light Gray
            135:55, ## Sand blue / Sand Blue
            136:54, ## Sand violet / Sand Purple
            138:69, ## Sand yellow / Dark Tan
            140:63, ## Earth blue / Dark Blue
            141:80, ## Earth Green / Dark Green
            143:74, ## Tr. Flu. Blue / Trans-Medium Blue
            145:78, ## Sand blue metallic / Metal Blue
            148:77, ## Mettalic Dark Grey / Pearl Dark Gray
            151:48, ## Sand green / Sand Green
            153:58, ## Sand red / Sand Red
            154:59, ## Dark red / Dark Red
##            158:, ## Tr. Flu. Red / 
##            180:, ## Curry / 
            191:110, ## Flame yellowish orange / Bright Light Orange
            192:88, ## Reddish Brown / Reddish Brown
            194:86, ## Medium stone grey / Light Bluish Gray
            198:93, ## Bright Reddish Lilac / Light Purple
            199:85, ## Dark stone grey / Dark Bluish Gray
            208:99, ## Light stone grey / Very Light Bluish Gray
            212:105, ## Light Royal blue / Bright Light Blue
            216:27, ## Rust / Rust
            217:91, ## Brown / Dark Flesh
            222:104, ## Light Purple / Bright Pink
            226:103, ## Cool Yellow / Bright Light Yellow
            232:87, ## Dove blue / Sky Blue
            268:89, ## Medium Lilac / Dark Purple
            }
        d_parts_pab2bricklink = {
            44237:2456,
            }

        ## initiate xml
        lines_xml = ['<INVENTORY>\n']

        for designID in d_shopping_list.keys():
            for materialID in d_shopping_list[designID].keys():

                count = d_shopping_list[designID][materialID]
                if count == 0:
                    continue

                color = d_colors_pab2bricklink[materialID]
                if designID in d_parts_pab2bricklink.keys():
                    itemID = d_parts_pab2bricklink[designID]
                else:
                    itemID = designID

                ## dont subtract inventory again...
##                if materialID in self.inventory.keys():
##                    if designID in self.inventory[materialID].keys():
##                        count -= self.inventory[materialID][designID]
                if count <= 0:
                    continue

                lines_xml += [' <ITEM>\n']
                lines_xml += ['  <ITEMTYPE>P</ITEMTYPE>\n']
                lines_xml += ['  <ITEMID>%s</ITEMID>\n' %(itemID)]
                lines_xml += ['  <COLOR>%s</COLOR>\n' %(color)]
                lines_xml += ['  <MINQTY>%s</MINQTY>\n' %(count)]
                lines_xml += ['  <CONDITION>N</CONDITION>\n'] ## new or used?
                ## append MAXPRICE if PAB prices were read from a dictionary file
                lines_xml += ['  <MAXPRICE>%s</MAXPRICE>\n' %(d_prices_PAB_USD[designID])]
                lines_xml += ['  <REMARKS>%s</REMARKS>\n' %(suffix,)]
                lines_xml += ['  <NOTIFY>N</NOTIFY>\n']
                lines_xml += ['  <WANTEDLISTID>%s</WANTEDLISTID>\n' %(WantedListID)]
                lines_xml += [' </ITEM>\n']
        lines_xml += ['</INVENTORY>\n']
        fp = 'xml/bricklink_%s' %(suffix)
        fp += '.xml'
        fd = open(fp,'w')
        fd.writelines(lines_xml)
        fd.close()

        return


    def check_vicinal_small_bricks(self,d_bricks_main,):

        '''speed up this function'''

        print 'small2large'

        l_layers = d_bricks_main.keys()
        l_layers.sort()

        if self.bool_replace_small_bricks_with_larger_bricks == True:

            ## also replace 2 * 1x2 with 1x4 !!!
            ## also replace 1x2 and 2x2 with 2x3
            ## also replace 1x2 and 1x1 with 1x3

            for layer in l_layers:

                if layer % 10 == 0:
                    print 'small2large', layer

                ## sorted by cheapest replacement (volume to price ratio)...
                for replacement in [
                    ## round 1
                    ## don't replace 1x1 and 1x2 with corner in 1st round, as I prefer 1x3
                    ## don't replace 1x1 and 1x3 with 1x3, as I prefer 1x4
                    '1x8','1x4','1x2','1x6',
                    ## round 2
                    ## new 2x2 to 2x4
                    ## new 1x3 to 2x3
                    ## new 1x2 to 2x2
                    ## new 1x2 to 1x4
                    ## new 1x2 and old 1x1 to 1x3
                    ## new 1x2 and old 1x1 to corner
                    '1x8','1x4','1x2','1x6','1x3',
                    ## round 3
                    ## new 2x2 to 2x4
                    ## new 1x4 to 2x4
                    ## new 1x3 to 2x3
                    ## new 1x3 to 1x6
                    '1x8','1x4','1x2','1x6','1x3',
                    ## round 4
                    ## new 1x4 to 2x4
                    '2x8','2x4','2x6',
                    '2x8','2x3','2x2',
                    '2x8','2x6','2x4','2x3','2x2',
                    ## i need to fix the corner function...                    
##                    'corner',

                    ## if bucket6177, then...
                    ## 2x8 only if red and yellow
                    ## 2x6 only if blue and yellow
                    ## 1x8 only if red and blue
                    ## 1x6 only if rbyw
                    
                    ]:

                    ##
                    ## loop over x
                    ##
                    l_tx = d_bricks_main[layer]['tx'].keys()
                    l_tx.sort()
                    for tx in l_tx:

                        ##
                        ## loop over z
                        ##
                        l_tz = d_bricks_main[layer]['tx'][tx]['tz'].keys()
                        l_tz.sort()
                        for tz in l_tz:

                            ## tz was deleted previously in the loop
                            if tz not in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                                continue
                            ## check that tx was not deleted in previous loop over tz2
                            if not tx in d_bricks_main[layer]['tx'].keys():
                                print 'txtxtx'
                                stop
                                break

                            designID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID']
                            if designID not in [
                                3005, ## 1x1 (to 1x2, 1x2x2)
                                3004, ## 1x2 (to 1x2x2, 2x2)
                                3003, ## 2x2 (to 2x4)
                                3622, ## 1x3 (to 2x3)
                                3010, ## 1x4 (to 2x4)
                                ]:
                                continue

                            materialID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID']

                            angle = d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
                            ay = d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay']

                            ##
                            ## loop over vicinal x
                            ##
                            for tx2 in range(tx-2,tx+2+1,):

                                ## check that there is a brick in the position
                                if not tx2 in d_bricks_main[layer]['tx'].keys():
                                    continue
                                ## check that tz was not deleted in previous loop over tz2
                                if not tz in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                                    break

                                ##
                                ## loop over vicinal z
                                ##
                                for tz2 in range(tz-2,tz+2+1,):

                                    ## check that there is a brick in the position
                                    if not tz2 in d_bricks_main[layer]['tx'][tx2]['tz'].keys():
                                        continue
                                    ## check that tx2 was not deleted in previous loop over tz2
                                    if not tx2 in d_bricks_main[layer]['tx'].keys():
                                        print 'tx2tx2tx2'
                                        stop3
                                        break
                                    ## check that tz was not deleted in previous loop over tz2
                                    if not tz in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                                        break
                                    ## check that tx was not deleted in previous loop over tz2
                                    if not tx in d_bricks_main[layer]['tx'].keys():
                                        print 'xxx'
                                        STOPSTOP
                                        break

                                    designID2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['designID']
                                    materialID2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['materialID']
                                    angle2 = d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['angle']

                                    ## obviously bricks have to be of the same color to be connected
                                    if materialID != materialID2:
                                        continue

                                    ## skip if brick1 == brick2
                                    if tx == tx2 and tz == tz2:
                                        continue

                                    if abs(tx-tx2) < 6 and abs(tz-tz2) < 6:

                                        if (
                                            replacement == '1x3'
                                            and
                                            (
                                                ## 1x1 (3005) and 1x2 (3004)
                                                ## not vice versa 3004 and 3005
                                                (designID == 3005 and designID2 == 3004)
                                                )
                                            ):

                                            if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
                                                continue

                                            d_bricks_main = self.replace_1x2_and_1x1_w_1x3_main(
                                                layer,angle,angle2,tx,tx2,tz,tz2,
                                                d_bricks_main,materialID,
                                                )

                                        elif (
                                            replacement == '2x3'
                                            and
                                            (designID == 3622 and designID2 == 3622)
                                            and
                                            (tx == tx2 or tz == tz2)
                                            ):

                                            if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
                                                continue

                                            if (
                                                (angle == 90 and angle2 == 0)
                                                or
                                                (angle == 0 and angle2 == 90)
                                                ):
                                                continue

                                            d_bricks_main = self.replace_1x3_and_1x3_w_2x3(
                                                layer,angle,angle2,tx,tx2,tz,tz2,
                                                d_bricks_main,materialID,
                                                )

                                        ##
                                        ## replace 1x1 (3005) and 1x2 (3004) by a corner
                                        ##
                                        elif (
                                            replacement == 'corner'
                                            and
                                            (designID == 3005 and designID2 == 3004)
                                            and
                                            ## corner pieces not in blue and orange
                                            materialID not in [23,106,]
                                            ):
                                            d_bricks_main = (
                                                self.replace_1x2_and_1x1_w_1x2x2_main(
                                                    d_bricks_main,layer,
                                                    angle,angle2,tx,tx2,tz,tz2
                                                    )
                                                )

                                        ## replace 1x1 and 1x1 (3005) by a 1x2 (3004)
                                        elif (
                                            replacement == '1x2'
                                            and
                                            designID == 3005 and designID2 == 3005
                                            and
                                            ## same row or column
                                            (tx == tx2 or tz == tz2)
                                            ):

                                            if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
                                                continue

                                            d_bricks_main = self.replace_1x1_and_1x1_w_1x2(
                                                d_bricks_main,layer,
                                                angle,angle2,tx,tx2,tz,tz2
                                                )

                                        ## replace 1x2 and 1x2 (3004) by a 2x2 (3003)
                                        elif (
                                            replacement == '2x2'
                                            and
                                            designID == 3004 and designID2 == 3004
                                            and
                                            (tx == tx2 or tz == tz2)
##                                            and
##                                            angle == angle2
                                            ):

                                            self.replace_1x2_w_2x2_main(
                                                layer,angle,angle2,tx,tx2,tz,tz2,
                                                d_bricks_main,materialID,
                                                )

                                        ## replace 2x2 and 2x2 (3003) by a 2x4
                                        elif (
                                            replacement == '2x4'
                                            and
                                            designID == 3003 and designID2 == 3003
                                            and
                                            (tx == tx2 or tz == tz2)
                                            ):

                                            d_bricks_main = self.replace_2x2_w_2x4_main(
                                                self,layer,angle,angle2,tx,tx2,tz,tz2,
                                                d_bricks_main,materialID,
                                                )

                                        ##
                                        ## replace 1x4 and 1x4 (3010) by a 2x4 (3001)
                                        ##
                                        elif (
                                            replacement == '2x4'
                                            and
                                            (designID == 3010 and designID2 == 3010) ## 1x4
                                            and
                                            (tx == tx2 or tz == tz2)
                                            ):

                                            d_bricks_main = self.replace_1x4_and_1x4_w_2x4_main(
                                                layer,angle,angle2,tx,tx2,tz,tz2,
                                                d_bricks_main,materialID,
                                                )

        print 'run loop until no replacements'   

        return d_bricks_main


    def replace_1x2_w_2x2_main(
        self,layer,angle,angle2,tx,tx2,tz,tz2,
        d_bricks_main,materialID,
        ):

        '''one of those ugly functions with lots of if statements... rewrite!'''

        if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
            return d_bricks_main

        if abs(tx-tx2) > 1 or abs(tz-tz2) > 1:
            return d_bricks_main

        ## bricks rotated relative to each other
        if (
            (angle == 90 and angle2 == 0)
            or
            (angle == 0 and angle2 == 90)
            ):
            return d_bricks_main

        if   layer % 2 == 0 and angle == 0 and tx == tx2 and tz2 == tz+1:
            d_bricks_main = self.replace_1x2_w_2x2(
                d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=1,tz_add=1,
                designID=3003,materialID=materialID,
                )
        elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz2 == tz-1:
            d_bricks_main = self.replace_1x2_w_2x2(
                d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,
                designID=3003,materialID=materialID,
                )
        elif layer % 2 == 0 and angle == 90 and tz == tz2 and tx == tx2+1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        elif layer % 2 == 0 and angle == 90 and tz == tz2 and tx == tx2-1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        elif layer % 2 == 1 and angle == 0 and tx == tx2 and tz2 == tz-1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        elif layer % 2 == 1 and angle == 0 and tx == tx2 and tz2 == tz+1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
            if self.verbose == True:
                print 'skip 2x2 for now'
        else:
            print 'new 2x2 case'
            print layer, angle, angle2, tx, tx2, tz, tz2
##            stop

        return d_bricks_main


    def replace_2x2_w_2x4_main(
        self,layer,angle,angle2,tx,tx2,tz,tz2,d_bricks_main,materialID,
        ):

        if layer % 2 == 0 and tx2 == tx and tz2 == tz-2:
            if angle == 90:
                stop
            if angle == 0 and angle2 == 90:
                return d_bricks_main
            print '111', angle, angle2, ay
            d_bricks_main = self.replace_2x2_w_2x4(
                d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=1,tz_add=2,materialID=221,
                ) ## dark green

        elif layer % 2 == 1 and tx2 == tx+2 and tz2 == tz:
            if angle == 0:
                stop
            d_bricks_main = self.replace_2x2_w_2x4(
                d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=1,materialID=221,
                ) ## bright purple

        elif layer % 2 == 1 and tx2 == tx and tz2 == tz+2:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 1 and tx2 == tx and tz2 == tz-2:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 0 and tx2 == tx-2 and tz2 == tz:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 0 and tx2 == tx+2 and tz2 == tz:
            if self.verbose == True:
                stop
            return d_bricks_main

        elif layer % 2 == 0 and tx2 == tx and tz2 == tz+2:
            if angle == 0 and ay == 1:
##                                            print '333', angle, angle2, ay
                d_bricks_main = self.replace_2x2_w_2x4(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=1,tz_add=2,
                    )
            ## 2x2 is newly added
            elif angle == 0 and angle2 == 90: ## this should be a general exclusion
                None
            ## 2x2 is newly added
            elif angle == 90 and angle2 == 0: ## this should be a general exclusion
                None
            else:
                print layer, angle, angle2, ay
                stop

        elif layer % 2 == 1 and tx2 == tx-2 and tz2 == tz:
            if angle == 90:
                stop
            else:
                stop

        else:
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = 221
            d_bricks_main[layer]['tx'][tx2]['tz'][tz2]['materialID'] = 221
            print 'new 2x4 case'
            print tx, tz, tx2, tz2
            print d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
            print 'replace 2 * 2x2 with 1 * 2x4'
            print 'layer', layer
            print 'angle', angle
##                                            stop4

        return d_bricks_main


    def replace_1x4_and_1x4_w_2x4_main(
        self,layer,angle,angle2,tx,tx2,tz,tz2,d_bricks_main,materialID,
        ):

        if abs(tx-tx2) != 1 and abs(tz-tz2) != 1:
            return d_bricks_main

        if (
            (angle == 90 and angle2 == 0)
            or
            (angle == 0 and angle2 == 90)
            ):
            return d_bricks_main

        if layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2-1:
            tz_new = max(tz,tz2)
            tx_new = tx
            pass
        elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2+1:
            tz_new = max(tz,tz2)
            tx_new = tx
            pass
        elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2-1:
            print layer
            stop1
            if self.verbose == True:
                print 'skip 2x4 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2+1:
            print layer
            stop2
            if self.verbose == True:
                print 'skip 2x4 replacement for now'
            return d_bricks_main

        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
            print layer
            stop3
            if self.verbose == True:
                print 'skip 2x4 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
            tz_new = tz
            tx_new = tx
            pass
        elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2-1:
            if self.verbose == True:
                print 'skip 2x4 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2+1:
            print layer
##                                                stop5
            if self.verbose == True:
                print 'skip 2x4 replacement for now'
            return d_bricks_main
        else:
            print layer, angle, angle2
            print tx, tx2, tz, tz2
            stop
        d_bricks_main = self.replace_1xx_and_1xx_w_2xx(
            d_bricks_main,layer,
            angle,angle2,tx,tx2,tz,tz2,
            tx_new,tz_new,
            materialID,
            designID_new = 3001, ## 2x4
            )

        return d_bricks_main


    def replace_1x3_and_1x3_w_2x3(
        self,layer,angle,angle2,tx,tx2,tz,tz2,d_bricks_main,materialID,
        ):
        
        if layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2-1:
            tz_new = max(tz,tz2)
            tx_new = tx
            pass
        elif layer % 2 == 0 and angle == 0 and tx == tx2 and tz == tz2+1:
            tz_new = max(tz,tz2)
            tx_new = tx
            pass
        elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2-1:
            print layer
            stop1
            if self.verbose == True:
                print 'skip 2x3 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 0 and tz == tz2 and tx == tx2+1:
            print layer
            stop2
            if self.verbose == True:
                print 'skip 2x3 replacement for now'
            return d_bricks_main

        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2+1:
            print layer
            stop3
            if self.verbose == True:
                print 'skip 2x3 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 90 and tz == tz2 and tx == tx2-1:
            tz_new = tz
            tx_new = tx
            pass
        elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2-1:
            if self.verbose == True:
                print 'skip 2x3 replacement for now'
            return d_bricks_main
        elif layer % 2 == 1 and angle == 90 and tx == tx2 and tz == tz2+1:
            print layer
    ##                                                stop5
            if self.verbose == True:
                print 'skip 2x3 replacement for now'
            return d_bricks_main
        else:
            print layer, angle, angle2
            print tx, tx2, tz, tz2
            stop
        d_bricks_main = self.replace_1xx_and_1xx_w_2xx(
            d_bricks_main,layer,
            angle,angle2,tx,tx2,tz,tz2,
            tx_new,tz_new,
            materialID,
            designID_new = 3002, ## 2x3
            )

        return d_bricks_main


    def replace_1x2_and_1x1_w_1x3_main(
        self,layer,angle,angle2,tx,tx2,tz,tz2,d_bricks_main,materialID,
        ):

        '''## not finished!!!'''
        ## not finished!!!

        l_returns = [
            [0,0,-1,0,0,],
            [0,-1,-1,0,0,],
            [0,0,1,0,0,],
            [0,1,1,0,0,],
            [1,1,1,90,0,],
            [1,0,-1,90,0,],
            ]

        l_skips = [
            [0,1,-1,0,0,],
            [0,1,0,0,90,],
            [1,1,0,90,90,],
            [1,-1,0,90,90,],
            [1,-1,-1,0,90,],
            [0,-1,1,0,0,],
            [0,1,-1,0,90,],
            [1,-1,-1,90,90,],
            [1,1,-1,90,90,],
            [0,-1,0,0,90,],
            [1,-1,1,90,90,],
            [1,-1,-1,90,0,],
            [1,0,1,90,0,],
            [0,-1,-1,0,90,],
            [0,1,1,0,90,],
            [0,-1,1,0,90,],
            [1,0,1,90,90,],
            [0,1,0,0,0,],
            [1,-1,1,90,0,],
            [1,1,-1,90,0,],
            [1,1,1,90,90,],
            [0,0,-1,0,90,],
            ]

        ## don't replace (corner)
        if [layer % 2, tx2-tx, tz2-tz, angle, angle2,] in l_returns:
            return d_bricks_main 

        if [layer % 2, tx2-tx, tz2-tz, angle, angle2,] in l_skips:
            if self.verbose == True:
                print 'skip 1x3 for now'
            return d_bricks_main 

        elif (
            layer % 2 == 0 and tx2 == tx and tz2 == tz+1
            and
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
            print 'new 1x3 case', layer, 'tx', tx,tx2,tz,tz2, angle, angle2
            print [layer % 2, tx2-tx, tz2-tz, angle, angle2,]
            return d_bricks_main
    ##                                                stop
        d_bricks_main = self.replace_1x2_and_1x1_w_1x3(
            d_bricks_main,layer,
            angle,angle2,tx,tx2,tz,tz2,
            tx_new, tz_new,
            )

        return d_bricks_main


    def dic2lines(self,d_bricks_main):

        lines_lxfml = []
        refID = 0

        l_layers = d_bricks_main.keys()
        l_layers.sort()
        for layer in l_layers:
            for tx in d_bricks_main[layer]['tx'].keys():
                for tz in d_bricks_main[layer]['tx'][tx]['tz'].keys():
                    refID += 1
                    designID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID']
                    materialID = d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID']
                    itemNo = '%s%s' %(designID,materialID)
                    ay = d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay']
                    angle = d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']


                    ## initiate
                    s = '        <Part'
                    ## refID
                    s += ' refID="%i"' %(refID)
                    ## dimension
                    s += ' designID="%s"' %(designID)
                    ## color
                    s += ' materialID="%s"' %(materialID)
                    ## dimension & color
                    s += ' itemNos="%s"' %(itemNo)
                    ## rotation
                    s += ' angle="%i"' %(angle)
                    s += ' ax="0" ay="%s" az="0"' %(ay)

                    ## position
                    s += ' tx="%f"' %(0.8*tx)

                    if self.brick_or_plate == 'brick':
                        s += ' ty="%f"' %(0.96*layer)
                    elif self.brick_or_plate == 'plate':
                        s += ' ty="%f"' %(0.32*layer)

                    s += ' tz="%f"' %(0.8*tz)

                    ## terminate
                    s += '/>\n'


                    lines_lxfml += [s]

        return lines_lxfml


    def replace_1xx_and_1xx_w_2xx(
        self,
        d_bricks_main,layer,
        angle,angle2,tx,tx2,tz,tz2,
        tx_new,tz_new,
        materialID,
        designID_new,
        ):

        if not tz_new in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                d_bricks_main[layer]['tx'][tx]['tz'][tz]
                )
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['materialID'] = (
            materialID
            )
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = (
            designID_new
            ) ## 3x2
        ## delete bricks from old pos
        if tx != tx_new or tz != tz_new:
            del d_bricks_main[layer]['tx'][tx]['tz'][tz]
        elif tx2 != tx_new or tz2 != tz_new:
            del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
        else:
            stop

        return d_bricks_main


    def replace_1x2_and_1x1_w_1x3(
        self,
        d_bricks_main,layer,
        angle,angle2,tx,tx2,tz,tz2,
        tx_new, tz_new,
        ):

        '''this function is most likely ver wrong and needs to be expanded.
only tested for 2dau_3 layer 30'''

        if not tx_new in d_bricks_main[layer]['tx'].keys():
            d_bricks_main[layer]['tx'][tx_new] = {'tz':{}}
        if not tz_new in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                d_bricks_main[layer]['tx'][tx]['tz'][tz]
                )
        ## rotation of 3004 1x2
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = angle2
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3622
        ## delete bricks from old pos
        if tx != tx_new or tz != tz_new:
            del d_bricks_main[layer]['tx'][tx]['tz'][tz]
        if tx2 != tx_new or tz2 != tz_new:
            del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]

        return d_bricks_main

    def replace_1x1_and_1x1_w_1x2(
        self,
        d_bricks_main,layer,
        angle,angle2,tx,tx2,tz,tz2
        ):

##        return d_bricks_main
        
        if layer % 2 == 0 and tx == tx2 and tz == tz2-1:
            tx_new = tx
            tz_new = tz+1
            if not tz_new in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
                d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                    d_bricks_main[layer]['tx'][tx]['tz'][tz]
                    )
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 90
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
            ## delete bricks from old pos
            if tx != tx_new or tz != tz_new:
                del d_bricks_main[layer]['tx'][tx]['tz'][tz]
            if tx2 != tx_new or tz2 != tz_new:
                del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
        elif layer % 2 == 1 and tx == tx2 and tz == tz2-1:
            tx_new = tx
            tz_new = tz
            if not tz_new in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
                d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                    d_bricks_main[layer]['tx'][tx]['tz'][tz]
                    )
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 90
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
            ## delete bricks from old pos
            if tx != tx_new or tz != tz_new:
                del d_bricks_main[layer]['tx'][tx]['tz'][tz]
            if tx2 != tx_new or tz2 != tz_new:
                del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
        elif tx == tx2-1 and tz == tz2:
            tx_new = tx
            tz_new = tz
            if not tz_new in d_bricks_main[layer]['tx'][tx_new]['tz'].keys():
                d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
                    d_bricks_main[layer]['tx'][tx]['tz'][tz]
                    )
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['angle'] = 0
            d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new]['designID'] = 3004
            ## delete bricks from old pos
            del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]
        else:
            print tx, tz, tx2, tz2
            print .8*tx, .8*tz, .8*tx2, .8*tz2
            print angle, layer
            print 'replace materialID with 106/221 in lxmfl file'
            stop_compare_with_and_without_replacement

        return d_bricks_main


    def replace_1x2_and_1x1_w_1x2x2_main(
        self,d_bricks_main,layer,
        angle,angle2,tx,tx2,tz,tz2,
        ):

        ## NB! Replacement of 1x2 and 1x1 by 1x2x2 corner bricks
        ## will make the model more stable and reduce the number of bricks
        ## but also increase the price

        if layer % 2 == 1 and tx2 == tx-1 and tz2 == tz-1: ## when not rotated
            if angle == 0:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=-1,tz_add=-1,angle=0,materialID=192,
                    ) ## reddish brown
            elif layer % 2 == 1 and angle == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=180,materialID=119,
                    ) ## bright yellowish green
            else:
                print layer, angle
                stop1a
        elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz-1: ## when not rotated
            if angle == 0:
                pass
            else:
                print layer, angle
                stopstop

        elif layer % 2 == 1 and tx2 == tx-1 and tz2 == tz: ## when not rot
            if angle == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=-1,tz_add=1,angle=270,
                    )
            elif angle == 0:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=-1,tz_add=1,angle=270,materialID=221,
                    ) ## bright purple
        elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz: ## when rotated!!!
            if angle == 0:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=-1,tz_add=1,angle=90,materialID=221,
                    ) ## dark green
            else:
                print layer, angle, angle2
                stop

        elif layer % 2 == 0 and tx2 == tx-1 and tz2 == tz+1: ## when no rot
            self.replace_1x2_and_1x1_w_1x2x2(
                d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=270,materialID=221,
                ) ## dark green
        elif layer % 2 == 1 and tx2 == tx-1 and tz2 == tz+1: ## no rota
            pass

        elif layer % 2 == 0 and tx2 == tx and tz2 == tz-1: ## when no rotation
            if angle == 0:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=90,materialID=192,
                    ) ## reddish brown
            else:
                print 'replace 1x2,1x1 with corner - check rotation'
                print layer, angle
                stopc
        elif layer % 2 == 0 and tx2 == tx and tz2 == tz+1: ## when no roation
            if angle == 0 and angle2 == 0:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=1,tz_add=1,angle=180,materialID=221,
                    ) ## bright purple
            elif angle == 0 and angle2 == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=20,angle=90,materialID=221,
                    ) ## bright purple
            else:
                print 'replace 1x2,1x1 with corner - check rotation'
                print layer, angle, angle2
                stopd

        elif layer % 2 == 1 and tx2 == tx+1 and tz2 == tz-1: ## kein rotation of 2x1!
            if angle == 90:
                pass
            else:
                print layer, angle
                stope
        elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz-1: ## when no rotatin
            pass
##
        elif layer % 2 == 0 and tx2 == tx+1 and tz2 == tz:
            if angle == 0 and angle2 == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=0,materialID=221,
                    ) ## bright purple
            elif angle == 0 and angle2 == 0:
                pass ## 1x3
            else:
                print layer, angle, angle2
                stop

        elif layer % 2 == 1 and tx2 == tx and tz2 == tz+1:
            if angle == 90:
                pass ## 1x3
            else:
                print layer, angle
                stop

        elif layer % 2 == 1 and tx2 == tx and tz2 == tz-1:
            if angle == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=270,materialID=221,
                    ) ## dark green
            else:
                print layer, angle
                stop
        
        elif layer % 2 == 1 and tx2 == tx+1 and tz2 == tz:
            if angle == 90:
                self.replace_1x2_and_1x1_w_1x2x2(
                    d_bricks_main,layer,tx,tz,tx2,tz2,tx_add=0,tz_add=0,angle=0,
                    )
            else:
                print 'replace 1x2,1x1 with corner - check rotation'
                print layer, angle
                stopf

        elif tx2 == tx+1 and tz2 == tz+1: ## when no rota
            pass

        else:
            print layer % 2 == 0, tx2 == tx-1, tz2 == tz
            print tx2-tx, tz2-tz
            print tx, tz, tx2, tz2
            print .8*tx, .8*tz, .8*tx2, .8*tz2
##            print materialID
            print 'layer', layer, angle
            print 'replace materialID with 221/106 in lxmfl file'
            stop1

        return d_bricks_main


    def replace_1x2_w_2x2(
        self,d_bricks_main,layer,tx,tz,tx2,tz2,tx_add,tz_add,designID,materialID=221,
        ):

        ## find consecutive occupants
        i_max = int( min( 4, max(self.d_designIDs[self.brick_or_plate][2]) ) )
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
                                ] == 3004: ## 1x2
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

        if bool_break == True:
            i -= 1

        tx_new = tx+tx_add
        tz_new = tz+tz_add

        if d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 90:
            tz_new += (i-1)*(tz2-tz)
        elif (
            layer % 2 == 0 and tx2 == tx and tz2 == tz+1
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
            ):
            tz_new += (i-1)*(tz2-tz)
        elif (
            layer % 2 == 0 and tx2 == tx and tz2 == tz-1
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
            ):
            tz_new += (i-1)*(tz2-tz)
        else:
            tx_new += (i-1)*(tx2-tx)

        ## change brick
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = (
            self.d_designIDs[self.brick_or_plate][2][i+1]
            ) ## 2x?
        if self.recolor == True:
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = materialID
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = abs(
            90-d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
            )
        ## copy brick to new pos
        if not tx_new in d_bricks_main[layer]['tx'].keys():
            d_bricks_main[layer]['tx'][tx_new] = {'tz':{}}
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
            d_bricks_main[layer]['tx'][tx]['tz'][tz]
            )
        ## delete brick from old pos
        del d_bricks_main[layer]['tx'][tx]['tz'][tz]

        return d_bricks_main


    def replace_1x2_and_1x1_w_1x2x2(
        self,d_bricks_main,layer,tx,tz,tx2,tz2,tx_add,tz_add,angle,materialID=221,
        ):

        '''replace 1x2 and 1x1 by a "corner" brick'''

        ## delete occupant/neighbour brick
        del d_bricks_main[layer]['tx'][tx2]['tz'][tz2]

        ## change brick
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = 2357
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = angle

        ## copy brick to new pos
        if not tx+tx_add in d_bricks_main[layer]['tx'].keys():
            d_bricks_main[layer]['tx'][tx+tx_add] = {'tz':{}} ## tmp!!!
        d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add] = dict(
            d_bricks_main[layer]['tx'][tx]['tz'][tz]
            )

        ## delete brick from old pos
        if not (tx_add == 0 and tz_add == 0):
            del d_bricks_main[layer]['tx'][tx]['tz'][tz]

        ## verbose when finding incorrectly placed bricks
        d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add][
            'materialID'
            ] = materialID

        if self.verbose == True:
            d_bricks_main[layer]['tx'][tx+tx_add]['tz'][tz+tz_add][
                'materialID'
                ] = materialID

        return d_bricks_main


    def replace_2x2_w_2x4(
        self,d_bricks_main,layer,tx,tz,tx2,tz2,tx_add,tz_add,materialID=1,
        ):

        ## find consecutive occupants
        i_max = int( max( self.d_designIDs[self.brick_or_plate][2] )/2. )
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
                                'designID'] == 3003: ## 2x2
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

        if bool_break == True:
            i -= 1

        tx_new = tx+tx_add
        tz_new = tz+tz_add

        if d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 90:
            tz_new += (i-1)*(tz2-tz)
        elif (
            layer % 2 == 0 and tx2 == tx and tz2 == tz+2
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
            ):
            tz_new += (i-1)*(tz2-tz)
        elif (
            layer % 2 == 0 and tx2 == tx and tz2 == tz-2
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] == 0
            and
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['ay'] == 1
            ):
            tz_new += (i-1)*(tz2-tz)
        else:
            tx_new += (i-1)*(tx2-tx)

        ## change brick
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['designID'] = (
            self.d_designIDs[self.brick_or_plate][2][2*(i+1)]
            ) ## 2x?
        if self.recolor == True:
            d_bricks_main[layer]['tx'][tx]['tz'][tz]['materialID'] = materialID
        d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle'] = abs(
            90-d_bricks_main[layer]['tx'][tx]['tz'][tz]['angle']
            )
        ## copy brick to new pos
        if not tx_new in d_bricks_main[layer]['tx'].keys():
            d_bricks_main[layer]['tx'][tx_new] = {'tz':{}}
        d_bricks_main[layer]['tx'][tx_new]['tz'][tz_new] = dict(
            d_bricks_main[layer]['tx'][tx]['tz'][tz]
            )
        ## delete brick from old pos
        del d_bricks_main[layer]['tx'][tx]['tz'][tz]

        return d_bricks_main


    def write_lxfml(self,pdb,lines_lxfml_body,):

        ## I need to add cameras at correct positions...
        ## get those camera positions from brick positions...

        l_head1 = [
'''<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<LXFML versionMajor="4" versionMinor="0" name="DNA_building_bricks">
  <Meta>
    <Application name="LEGO Digital Designer" versionMajor="4" versionMinor="0"/>
    <Brand name="Factory"/>
    <BrickSet version="287"/>
  </Meta>
  <Cameras>
'''
]
        l_cameras = [
'''    <Camera refID="0" fieldOfView="80" distance="200" angle="60" \
ax="-0.6" ay="-0.6" az="-0.25" tx="0.0" ty="0.0" tz="0.0"/>
'''
]

        l_head2 = [
'''  </Cameras>
  <Scene cameraRefID="0">
    <Model>
      <Group refID="0" angle="0" ax="0" ay="1" az="0" tx="0" ty="0" tz="0">
'''
]

        l_tail1 = [
'''      </Group>
    </Model>
  </Scene>
  <BuildingInstructions>
    <BuildingInstruction>
'''
]
        l_tail2 = [
'''
    </BuildingInstruction>
  </BuildingInstructions>
</LXFML>
'''
]

        lines_lxfml_instructions = []
        tx_min = None
        tx_max = None
        tz_min = None
        tz_max = None
        for line in lines_lxfml_body:
            d = {'refID':None,'tx':None,'ty':None,'tz':None,}
            for k in d.keys():
                index1 = line.index(k)+len(k)+2
                index2 = index1+line[index1:].index('"')
                v = line[index1:index2]
                d[k] = v
            refID = int(d['refID'])
            tx = float(d['tx'])
            ty = float(d['ty'])
            tz = float(d['tz'])
            if tx_min == None or tx < tx_min: tx_min = tx
            if tx_max == None or tx < tx_min: tx_max = tx
            if tz_min == None or tz < tx_min: tz_min = tz
            if tz_max == None or tz < tx_min: tz_max = tz
        for line in lines_lxfml_body:
            d = {'refID':None,'tx':None,'ty':None,'tz':None,}
            for k in d.keys():
                index1 = line.index(k)+len(k)+2
                index2 = index1+line[index1:].index('"')
                v = line[index1:index2]
                d[k] = v
            refID = int(d['refID'])
            tx = float(d['tx'])+50
            ty = float(d['ty'])+40
            tz = float(d['tz'])+40

            if refID == 1 or (refID) % 10 == 1:

                cameraRefID = (refID+1) / 10
                l_cameras += [
                    '    <Camera refID="%s" fieldOfView="40" distance="40" angle="56.6" \
ax="-0.59" ay="0.77" az="0.2445" tx="%f" ty="%f" tz="%f"/>\n' %(
                        cameraRefID,tx,ty,tz,
                        )
                    ]

                lines_lxfml_instructions += [
                    '    <Step name="Step%i" cameraRefID="%i">\n' %(refID,cameraRefID),
                    ]

            lines_lxfml_instructions += ['      <PartRef partRefID="%i"/>\n' %(refID)]

            if (refID) % 10 == 0:

                lines_lxfml_instructions += ['    </Step>\n']

            tx_prev = tx
            ty_prev = ty
            tz_prev = tz

        if (refID) % 10 != 0:

            lines_lxfml_instructions += ['    </Step>\n']

        lines = l_head1+l_cameras+l_head2+lines_lxfml_body+l_tail1
##        lines += lines_lxfml_instructions
        lines += l_tail2
##        lines = l_head1+l_cameras+l_head2+lines_lxfml_body+l_tail1+l_tail2

        fn = os.path.join('lxfml','%s_%s_%s' %(
            pdb,self.size_scale,self.brick_or_plate,))
        if self.bool_hollow == True:
            fn += '_hollow'
        if self.bool_color_by_residue == True:
            fn += '_colorbyresidue'
        if self.bool_buried_then_black == True:
            fn += '_buriedblack'
        fn += '.lxfml'
        fd = open(fn,'w')
        fd.writelines(lines)
        fd.close()

        return


    def split_consecutive_length_into_bricks(
        self,
        consecutive, materialID, width,
        d_layers, z_brick, element, k, xy,
        bool_next_rotated, bool_prev_skipped,
        ):

        d_bricks = {
            'rotated_before':[],'rotated_after':[],
            'normal':[],'flying':[],
            }

        l_lengths = self.d_brick_sizes[self.brick_or_plate][materialID][width]
        l_lengths.sort()
        max_length = l_lengths[-1]

        ## bricks exist
        if (
            consecutive[1] > 1
            and
            consecutive[1] in self.d_brick_sizes[self.brick_or_plate][materialID][width]
            ):
            d_bricks['normal'] += [consecutive]
        ## odd brick length higher than 1 for which brick exists,
        ## if combined with 3x1 brick
        elif (
            consecutive[1] > 1
            and
            consecutive[1]-3 in self.d_brick_sizes[self.brick_or_plate][materialID][width]
            ):
            d_bricks['normal'] += [[consecutive[0],3],[consecutive[0]+3,consecutive[1]-3],]

        ##
        ## current brick is a 1x1 brick. anything above or below to attach to?
        ##
        elif consecutive[1] in [1,]:

            if bool_next_rotated == True and consecutive[2] == 'rotated':
                l_below_all_elements = self.check_bricks_vicinal_all_elements(
                    d_layers,z_brick-1,k,xy+consecutive[3],consecutive,
                    )
                l_above_all_elements = self.check_bricks_vicinal_all_elements(
                    d_layers,z_brick+1,k,xy+consecutive[3],consecutive,
                    )
            else:
                l_below_all_elements = self.check_bricks_vicinal_all_elements(
                    d_layers,z_brick-1,k,xy,consecutive,
                    )
                l_above_all_elements = self.check_bricks_vicinal_all_elements(
                    d_layers,z_brick+1,k,xy,consecutive,
                    )
            ## more options for attachment if width is 2 (xy+1 could also be xy-1...)
                if width == 2:
                    l_below_all_elements += self.check_bricks_vicinal_all_elements(
                        d_layers,z_brick-1,k,xy+1,consecutive,
                        )
                    l_above_all_elements += self.check_bricks_vicinal_all_elements(
                        d_layers,z_brick+1,k,xy+1,consecutive,
                        )

            ## nothing above or below current 1x1 brick
            if len(l_below_all_elements) == 0 and len(l_above_all_elements) == 0:

                if self.verbose == True:
                    print 'use 2x1 brick reversed',
                    print ', element =', element, ', k =', k,
                    print ', xy =', xy, consecutive, ', refID =', self.refID

                bool_flying = False

                ## already accomodated for rotated 2x1 or 3x1 brick in previous row
                if (
                    ## same element in previous row
                    xy-1 in d_layers[z_brick][element][k].keys()
                    and
                    ## same element in same position in previous row
                    consecutive[0] in d_layers[z_brick][element][k][xy-1]
                    and
                    ## not accomodating in next row after being accomodated itself?
                    bool_next_rotated == False
                    ):
                    before_or_after = 'before'
                    pass
                ## accomodate for rotated 2x1 or 3x1 brick in next row,
                ## if current brick is a 1x1 brick not attached above or below
                elif (
                    ## same element in next row
                    xy+1 in d_layers[z_brick][element][k].keys()
                    and
                    ## same element in same position in next row
                    consecutive[0] in d_layers[z_brick][element][k][xy+1]
                    ):
                    ## remove same element in same position in next row
                    d_layers[z_brick][element][k][xy+1].remove(consecutive[0])
                    before_or_after = 'after'
                ## no brick of any element above or below
                ## and no bricks in the same position before or after!
                ## 1x1 brick will hang no matter what!
                else:
                    bool_flying = True
                    if self.verbose == True:
                        print 'cant connect the brick to anything! flying brick!'
                        print xy, 'not below', d_layers[z_brick-1][element][k].keys()
                        print xy, 'not above', d_layers[z_brick+1][element][k].keys()
                        print xy-1, 'not before', d_layers[z_brick][element][k].keys()
                        print xy+1, 'not after', d_layers[z_brick][element][k].keys()
                        print 'no consecutive bricks',
                        print d_layers[z_brick][element][k][xy], '-->', consecutive

                if bool_flying == False:
                    if before_or_after == 'before':
                        ## previous brick was skipped to acomodate for brick previous to it
                        if bool_prev_skipped == True:
                            width += 1
                        d_bricks['rotated_before'] += [[consecutive[0],width+1,]]
                    elif before_or_after == 'after':
                        d_bricks['rotated_after'] += [[consecutive[0],width+1,]]
                else:
                    d_bricks['flying'] += [[consecutive[0],1,]]

            ## something above or below current xy
            ## elif len(l_below_all_elements) != 0 or len(l_above_all_elements) != 0:
            else:

                d_bricks['normal'] += [consecutive]
##                print 'ggg', consecutive

        ## brick doesn't exist
        elif (
            consecutive[1]
            not in
            self.d_brick_sizes[self.brick_or_plate][materialID][width]
            ):
            ## even consecutive length
            if consecutive[1] % 2 == 0:
                length = 0
            ## odd consecutive length
            else:
                length = 3
                d_bricks['normal'] += [[consecutive[0],3,]]
            max_length = max(self.d_brick_sizes[self.brick_or_plate][materialID][width])
            while length < consecutive[1]:
                ## avoid using length 2!!!
                ## check that next piece is not a 2!!! remainder length...
                if max_length < consecutive[1]-length:
                    brick_length = max_length
                elif (
                    consecutive[1]-length
                    not in
                    self.d_brick_sizes[self.brick_or_plate][materialID][width]
                    ):
                    brick_length = l_lengths[-2]
                else:
                    brick_length = consecutive[1]-length
                d_bricks['normal'] += [[consecutive[0]+length,brick_length,]]
                length += brick_length
        else:
            print
            print consecutive
            print materialID, element
            stop

        return d_bricks, d_layers


    def modify_consecutives_exposed(self,l_consecutives_exposed,):

        ## one stretch of consecutive exposed bricks
        if len(l_consecutives_exposed) == 1:
            pass

        ## two stretches of consecutive exposed bricks
        elif len(l_consecutives_exposed) == 2:
            ## two 1x1 bricks
            if l_consecutives_exposed[0][1] == 1 and l_consecutives_exposed[1][1] == 1:
                stop_expected
            ## one 1x1 brick
            elif l_consecutives_exposed[0][1] == 1 and l_consecutives_exposed[0][1] > 1:
                print l_consecutives_exposed
                print l_consecutives_exposed[0][0]+l_consecutives_exposed[0][1]-1
                print l_consecutives_exposed[1][0]-1
                stop1
            ## one 1x1 brick
            elif l_consecutives_exposed[0][1] > 1 and l_consecutives_exposed[1][1] == 1:
                print l_consecutives_exposed
                ## room for a 2x1 brick?
                if (
                    l_consecutives_exposed[0][0] + l_consecutives_exposed[0][1] - 1
                    <
                    l_consecutives_exposed[1][0]-1
                    ):
                    l_consecutives_exposed[1] = [l_consecutives_exposed[1][0],2,]
            ## no 1x1 bricks
            else:
                pass

        ## three stretches of consecutive exposed bricks
        elif len(l_consecutives_exposed) == 3:
            for i in range(len(l_consecutives_exposed)):
                if l_consecutives_exposed[i][1] == 1:
                    print l_consecutives_exposed
                    stop

        ## three or more stretches of consecutive exposed bricks???
        else:
            print l_consecutives_exposed
            stop_expected

        return l_consecutives_exposed


    def append_bricks(
        self,
        ## input
        d_bricks, materialID, width, k,
        layer, xy,
        ## output (append)
        lines_lxfml_body,
        d_bricks_main,
        ):

        if layer % 2 == 0:
            k1 = k = 'y'
            k2 = 'x'
        else:
            k1 = k = 'x'
            k2 = 'y'

        for k_rot in ['normal','rotated_before','rotated_after','flying',]:

##            ## change xy
##            if k_rot == 'rotated':
##                xy -= 1

            if self.bool_color_flying_bricks == True and k_rot in [
##                'rotated_before','rotated_after',
                'flying',
                ]:
                materialID = 221 ## bright purple

            for brick in d_bricks[k_rot]:
                length = brick_len = brick[1] ## use to determine designID (length and width)
                brick_pos = brick[0] ## use to determine brick position
                designID = self.d_designIDs[self.brick_or_plate][width][brick_len]

                ##
                ## bricks aligned with "one" axis
                ##
                if k == 'y':
                    tx = brick_pos
                    tz = -xy
                    ay = 1
                    angle = 0
                    dangle = 90
    ##                tz *= -1 ## tmp!!!
                    if k_rot == 'rotated_before':
                        tz += 1 ## new!!!
                        tx -= 1 ## new!!!
                        tx += 1
                        tz -= 1
                    elif k_rot == 'rotated_after':
                        tz -= 1 ## new!!!
##                        tx += 2
##                        tx += 20

                ##
                ## bricks aligned with "another" axis
                ##
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
                        tz += 1 ## new!!!


                if k_rot in ['rotated_before','rotated_after',]:
                    angle += dangle
                    ay *= -1

                if brick_len == 1 and width == 2:
                    angle += dangle

                self.volume += width*brick_len

                lines_lxfml_body = self.append_lxfml_line(
                    designID,materialID,angle,ay,tx,layer,tz,lines_lxfml_body,
                    )

                if not layer in d_bricks_main.keys():
                    d_bricks_main[layer] = {
                        'tx':{},
                        }
                if not tx in d_bricks_main[layer]['tx'].keys():
                    d_bricks_main[layer]['tx'][tx] = {'tz':{},}
                d_bricks_main[layer]['tx'][tx]['tz'][tz] = {
                    'designID':designID,'materialID':materialID,
                    'angle':angle,'ay':ay,
                    }

        return lines_lxfml_body, d_bricks_main


    def append_lxfml_line(
        self,designID,materialID,angle,ay,tx,layer,tz,lines_lxfml_body,
        ):

        itemNo = '%s%s' %(designID,materialID)

        ## initiate
        s = '        <Part'
        ## refID
        s += ' refID="%i"' %(self.refID)
        ## dimension
        s += ' designID="%s"' %(designID)
        ## color
        s += ' materialID="%s"' %(materialID)
        ## dimension & color
        s += ' itemNos="%s"' %(itemNo)
        ## rotation
        s += ' angle="%i"' %(angle)
        s += ' ax="0" ay="%s" az="0"' %(ay)

        ## position
        s += ' tx="%f"' %(0.8*tx)

        if self.brick_or_plate == 'brick':
            s += ' ty="%f"' %(0.96*layer)
        elif self.brick_or_plate == 'plate':
            s += ' ty="%f"' %(0.32*layer)

        s += ' tz="%f"' %(0.8*tz)

        ## terminate
        s += '/>\n'

        ## append line
        lines_lxfml_body += [s]

        self.refID += 1

        return lines_lxfml_body


    def identify_consecutive_stretches(self,l_coords):

        l_consecutives = []
        j = 0
        i0 = 0
        for i in range(len(l_coords)):
            brick_length = 0
            if i+1 == len(l_coords):
                l_consecutives += [[l_coords[i0],i-i0+1,]]
                break
            elif l_coords[i+1] != l_coords[i]+1:
                l_consecutives += [[l_coords[i0],i-i0+1,]]
                i0 = i+1
                continue
            else:
                continue

        return l_consecutives


    def check_bricks_vicinal_all_elements(self,d_layers,z_brick,k,xy,consecutive,):

        if not z_brick in d_layers.keys():
            l_shared = []
        else:
            l_yx = range(consecutive[0],consecutive[0]+consecutive[1])
            l_shared = []
            for element in d_layers[z_brick].keys():
                if not xy in d_layers[z_brick][element][k].keys():
                    continue
                set_shared = set(l_yx) & set(d_layers[z_brick][element][k][xy])
                l_shared += list(set_shared)

        return l_shared


    def plot(self,pdb,layer,l_minmax_bricks,):

        lines = []
    ##    lines += ['set terminal postscript eps enhanced color "Helvetica" 48\n']
        lines += ['set terminal png xffffff x000000\n']
    ##    lines += ['set output "%s.ps"\n' %(layer)]
        lines += ['set output "png/%s_%s_%s.png"\n' %(pdb,self.size_scale,layer+1)]
        lines += ['set size 4,4\n']
        lines += ['plot ']
        lines += ['[%s:%s][%s:%s]' %(
            l_minmax_bricks[0][0], l_minmax_bricks[0][1],
            l_minmax_bricks[1][0], l_minmax_bricks[1][1],
            )]

        if self.bool_gnuplot_previous_layer == True:
            if layer != 0:
                ## triangle
                lines += ['"LEGO_C_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
                lines += [',"LEGO_N_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
                lines += [',"LEGO_O_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
                lines += [',"LEGO_P_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
                lines += [',"LEGO_H_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
                lines += [',"LEGO_buried_prev.gnuplotdata" ps 5 pt 10 lc rgb "black" t ""']
            if layer != 0:
                ## black circle
                lines += [',"LEGO_C.gnuplotdata" ps 4 pt 7 lc rgb "black" t ""']
            else:
                ## black circle
                lines += [ '"LEGO_C.gnuplotdata" ps 4 pt 7 lc rgb "black" t ""']
        else:
            ## black circle
            lines += [ '"LEGO_C.gnuplotdata" ps 4 pt 7 lc rgb "black" t ""']

        ## blue square
        lines += [',"LEGO_N.gnuplotdata" ps 4 pt 5 lc 3 t ""']
        ## red triangle1
        lines += [',"LEGO_O.gnuplotdata" ps 4 pt 9 lc 1 t ""']
        ## white diamond
        lines += [',"LEGO_H.gnuplotdata" ps 4 pt 12 lc rgb "black" t ""']
        ## orange circle
        lines += [',"LEGO_P.gnuplotdata" ps 4 pt 7 lc rgb "orange" t ""']

        ## green triangle2
        lines += [',"LEGO_buried.gnuplotdata" ps 3 pt 11 lc rgb "green" t ""']
        ##  green circle
        lines += [',"LEGO_support.gnuplotdata" ps 2 pt 7 lc rgb "green" t ""']

        ## triangle
##        lines += [',"LEGO_ZN.gnuplotdata" ps 3 pt 9 lc 6 t ""']
        ## green square
        lines += [',"LEGO_CL.gnuplotdata" ps 3 pt 5 lc rgb "green" t ""']
        ## green triangle2
        lines += [',"LEGO_119.gnuplotdata" ps 3 pt 11 lc rgb "green" t ""']
        lines += ['\n']
        fd = open('LEGO.gnuplotsettings','w')
        fd.writelines(lines)
        fd.close()

        os.system('"exe\\wgnuplot.exe" LEGO.gnuplotsettings')
    ##    os.system('convert.exe %s.ps %s.png' %(layer,layer,))

        volume = 0
        for element in self.d_materialIDs.keys()+['buried']:

            if os.path.isfile('LEGO_%s_prev.gnuplotdata' %(element)):
                os.remove('LEGO_%s_prev.gnuplotdata' %(element))
            if not os.path.isfile('LEGO_%s.gnuplotdata' %(element)):
                continue

            ## volume calculation (number of data lines)
            fd = open('LEGO_%s.gnuplotdata' %(element),'r')
            lines = fd.readlines()
            fd.close()
            volume += len(lines)

            os.rename(
                'LEGO_%s.gnuplotdata' %(element),
                'LEGO_%s_prev.gnuplotdata' %(element),
                )

        print 'gnuplot volume (number of lines)', volume

        return


    def __init__(self,):

        if os.path.isdir('dic'):

            self.inventory = {}
            fd = open('dic/d_inventory.dic','r')
            s = fd.read()
            fd.close()
            self.inventory = eval(s)

            fd = open('dic/d_radii.dic','r')
            s = fd.read()
            fd.close()
            self.d_radii = eval(s)

            fd = open('dic/d_materials.dic','r')
            s = fd.read()
            fd.close()
            self.d_materialIDs = eval(s)

            fd = open('dic/d_designIDs.dic','r')
            s = fd.read()
            fd.close()
            self.d_designIDs = eval(s)

            fd = open('dic/d_weights.dic','r')
            s = fd.read()
            fd.close()
            self.d_weights = eval(s)

        ## small molecules
        ## lego is made from styrene,  acrylonitrile, (poly)1,3-butadiene
        pdb = 'met' ## Met102 of T4L l2zm
        pdb = 'TPP' ## Vitamin B1
        pdb = 'SRO' ## Serotonin
        pdb = 'CFF' ## Caffeine
        pdb = 'B12' ## Vitamin B12
        pdb = 'his' ## His15 of HEWL NMR 1e8l
        ## DNA
##        pdb = '1d20' ## first NMR from 1990 (TCTATCACCG decamer, straight)
##        pdb = '1d18' ## first NMR from 1990 (CATGCATG octamer, twisted)
        pdb = '1bna' ## first B-DNA x-ray from 1981 (CGCGAATTCGCG dodecamer, B-DNA)
        pdb = '2dau_strand2'
        pdb = '2dau_strand1'
        pdb = '2dau' ## NMR of Dickerson-Drew Dodecamer (straight, NMR, B-DNA)
        ## proteins
        pdb = '2h35' ## only NMR model of Hemoglobin
        pdb = '1y8b' ## longest chain solved by NMR (model is 20m long!!!)
        pdb = '1erp' ## short alpha helix protein
        pdb = '1ed7' ## short beta sheet protein
        pdb = '4INS' ## replaces first insulin structure (1ins)
        pdb = '2LZM_helix' ## helix residues 60-80
        pdb = '1E8L_helix' ## helix residues 88-101
        pdb = '1EMA' ## GFP
        pdb = '1EMA_backbone' ## GFP
        pdb = '1mbn'

        self.pdb = pdb

        ## dimension and content
        self.size_scale = 5 ## default = 3
        self.exclude_hydrogens = False ## default = False
        self.brick_or_plate = 'brick' ## default = brick
        self.rotate_x90 = False ## default = False
        self.bool_replace_small_bricks_with_larger_bricks = True ## default = True
        self.bool_color_by_residue = False
        if self.size_scale >= 3:
            self.bool_grained = True ## get rid of flying 1x1 pieces
            self.exclude_long_bricks_not_in_bucket6177 = True ## price/stability
            self.bool_buried_then_black = False ## default = True
            self.bool_hollow = False ## default = False (use if extra large models!)
        elif self.size_scale <= 2:
            self.bool_grained = False ## get rid of flying 1x1 pieces
            self.exclude_long_bricks_not_in_bucket6177 = True ## price/stability
            self.bool_buried_then_black = False ## default = False
            if self.size_scale == 1:
                self.bool_hollow = True ## default = False (use if extra large models!)
            elif self.size_scale == 2:
                self.bool_hollow = False ## default = False (use if extra large models!)
        self.bool_use_double_width = False
        self.replace_atom = 'CL' ## if bool_buried_then_black == True

        ## verbose in test phase
        self.bool_color_flying_bricks = False ## default = False
        self.recolor = False ## default = False
        self.verbose = False ## default = False
        ## use gnuplot to identify bricks that are not placed
        self.bool_gnuplot = False ## default = False
        self.bool_gnuplot_previous_layer = True
        self.layer_range = None ## default = None

        ## output
        ## lxfml
        self.bool_split_into_layers = False ## default = False
        ## htm shopping list
        self.bool_buy_buckets = True ## default = True, bucket6177_savings.xlsx
        ## htm shopping list
        self.bool_subtract_inventory = False
        ## htm shopping list
        self.currency = 'USD'
        self.currency = 'EUR'

        ##        
        ## non default settings (for testing and other purposes)
        ##

        self.bool_split_into_layers = True ## default = False
##        self.bool_gnuplot = True
        self.bool_gnuplot_previous_layer = True

        ##
        ## input dependent dictionaries
        ##

        if os.path.isdir('dic'):

            ## brick sizes
            if self.exclude_long_bricks_not_in_bucket6177 == False:
                fn = 'dic/d_brick_sizes_PAB.dic'
            elif self.exclude_long_bricks_not_in_bucket6177 == True:
                fn = 'dic/d_brick_sizes_bucket6177.dic'
            fd = open(fn,'r')
            s = fd.read()
            fd.close()
            self.d_brick_sizes = eval(s)

            ## brick and set prices
            fd = open('dic/d_prices_PAB_%s.dic' %(self.currency),'r')
            s = fd.read()
            fd.close()
            self.d_prices_PAB = eval(s)

        return


if __name__ == '__main__':
    instance_LEGO = LEGO()
    instance_LEGO.main()
