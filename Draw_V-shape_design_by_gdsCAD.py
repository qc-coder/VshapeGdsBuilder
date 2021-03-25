'''

Making a design of new V-shape transmons of 'StriplineGeometry'
and 'CircularGeometry'

Code written for Python 2.7.10.
gdsCAD does not work with Python 3

Script made by Vladimir Milchakov [ vladimir_ph@protonmail.com ]
Jul-Oct 2020

::Wafer5 COMMENTS::
- Dose of pads increased from 9.5 to 9.8 C/m2 (comp to wafer3, wafer4)
- dosetests removed
- wires of 4probe-test-structures layer: 4->5. Than dose decreased from 12->10 C/m2
- change width of Ljj to 200 nm back. (we need more resistance)
- change deducted_x_wire_sift = 0.4->0.43 ~!(both in function gen_junction and in globals()). it corresponds to change bridge width: 200nm ->175nm
- maybe the angle will be changed as well (but its still ok with this .gds)
- changed back gap_X_JJ_to_L = 1.8->1.4  (long time ago were same)
- changed squid_Loop_wid  = 1.94 -> 1.3 (to make possible to make small Ljj)

::Wafer6 COMMENTS
- made better symmetry of junctions by shifting it (but this function is disabled now)
- increased width of horizontal microwire in hearts and in 4probe tests
- SQUID size automatically decrease if length of junction in SQUID chain is too small to avoid cut

::Wafer7 COMMENTS
- add list_of_L_JJ_high_sizes to vary same time Length and High of Ljj (big junctions)
- fix the problem with not full crosssection of small JJ Y by shifting X-wires (bigger bridge) "deducted_x_wire_sift"
- increasing size of Big junctions x1,x2,x3,x4 of biggest of last time

::Wafer8 COMMENTS
- same design for two gemini wafers: "Wafer8-1" & "Wafer8-2"
- make more dense JJ parameters (like wafer7ch16)
- improved naming of folder with .njf files
- (for next time) move chip30 structure (heartless) from this horrible spot, and merge two heartless chips to one
- (for NT) !!! Pads_wires_Circulat|Stripline_HEARTLESS.npf does not exist now. (either prepare it, either set it same as normal ones in .njf file)
'''

from __future__ import division
from gdsCAD import *
import numpy as np
import scipy.optimize
import os

################################################################################
#####______________Globals_________________________________________________#####

##__adress_to_save_.gds-file__##
name   = 'Wafer9'
main_folder = 'D:\\Data\\=Exp_Data=\\DESIGNS\\Wafers\\' +name +'\\' +name+'_Structures\\'

def createfolder(folder):
    try:
        os.stat(folder)
    except:
        os.mkdir(folder)
    return True
createfolder(main_folder)


##################################################################
##_______FABRICATION_MASK_PARAMETERS___________________________###
#######################
class BilayerMask:
    low_resist_thick = 0 ##[um]
    top_resist_thick = 0 ##[um]
    evap_angle       = 0
    me1_thick        = 0 ##[um]
    def __init__(self, low_resist_thick, top_resist_thick, evap_angle, me1_thick):
        self.low_resist_thick    = float(low_resist_thick)
        self.top_resist_thick    = float(top_resist_thick)
        self.evap_angle          = float(evap_angle)
        self.me1_thick           = float(me1_thick)
    def as_list(self):
        list_param = [self.low_resist_thick, self.top_resist_thick, self.me1_thick, self.evap_angle ]
        return list_param
    ######____________Predictions____________######
    def coordinates_deposition(self, L_mask):
        [d1,d2,t1,O] = self.as_list()
        tg0 = np.tan(np.deg2rad(35))
        x1 = (d1+d2)*tg0
        x2 = L_mask +d1*tg0
        z1 = -d1*tg0
        z2 = L_mask-t1 -(d1+d2+t1)*tg0
        return [x1,x2,z1,z2]
    #####____big_Junctuion____#####
    def Area_of_Big_junction_from_sizes_of_mask(self, wid_mask, high_mask):
        '''
        high_mask is about 0.2 - 0.4 um usually
        wid_mask is about  5.0 - 10.0 um
        '''
        [x1,x2,z1,z2] = self.coordinates_deposition(wid_mask)
        Ljj = z2-x1
        S = Ljj*high_mask
        return S
    def Lmask_from_Lproj_Big_junction(self, Lproj):
        [d1,d2,t1,O] = self.as_list()
        tg0 = np.tan(np.deg2rad(35))
        Lmask = Lproj +t1*(1+tg0) +2*(d1+d2)*tg0
        return Lmask
    #####____small_Y-junctuion____#####
    def Expected_width_of_small_Y_junction_from_width_of_mask(self, W_mask):
        '''
        Calcilate widthes of real wires after 1-st and 2-nd depositions
        returns list [W1_jj, W2_jj]
        takes list of mask parammeters = [d1,d2,t1,Theta_deg], and width wire in design
        '''
        SQRT2 = np.sqrt(2)
        L_mask = SQRT2*W_mask
        [x1,x2,z1,z2] = self.coordinates_deposition(L_mask)
        W1_jj = (x2-x1)/SQRT2
        W2_jj = (z2-z1)/SQRT2
        return [W1_jj, W2_jj]
    def Expected_AREA_of_small_Y_junction_from_width_of_mask(self, W_mask):
        [W1_jj, W2_jj] = self.Expected_width_of_small_Y_junction_from_width_of_mask(W_mask)
        area = W1_jj * W2_jj
        return area
    def Wmask_from_Area_JJ_Y_junction_analitical(self, S):
        '''
        Calculate what width should have wires of small Y-jj on mask
         to make junction of given area with given bi-layer mask parameters
        '''
        [d1,d2,t1,O] = self.as_list()
        tg0 = np.tan(np.deg2rad(35))
        D = d2*tg0
        T = t1*(1+tg0)
        radical = np.sqrt( 8*S +T*(T+1) )
        L_mask = 0.5 *( 2*D + radical)
        W_mask = L_mask/np.sqrt(2)
        return W_mask
    def Wmask_from_Area_JJ_Y_junction_iterative(self, WANTED_AREA):
        def minimization_func(x, wanted_real_area=WANTED_AREA):
            expext_Area = self.Expected_AREA_of_small_Y_junction_from_width_of_mask(x)
            return  abs(expext_Area- wanted_real_area)
        func = minimization_func
        init_guess = np.sqrt(WANTED_AREA)
        result = scipy.optimize.minimize(func, init_guess)
        x0 = result.x[0]
        return x0
#######################
##_____________________________________________________________###
'''
SET THE MASK
'''
# Target_Mask = BilayerMask(0.725, 0.250, 35, 0.020)
# Remy_Mask   = BilayerMask(0.745, 0.260, 35, 0.020)
# Wafer1_Mask = BilayerMask(0.704, 0.227, 35, 0.020)
# Wafer2_Mask = BilayerMask(0.693, 0.230, 35, 0.020)
# Wafer3_Mask = BilayerMask(0.692, 0.235, 35, 0.020)
# Wafer4_Mask = BilayerMask(0.685, 0.228, 35, 0.0202)
# Wafer5_Mask = BilayerMask(0.674, 0.226, 35, 0.0202)
# Wafer7_Mask = BilayerMask(0.697, 0.223, 35, 0.0202)
Wafer8_Mask = BilayerMask(0.710, 0.220, 35, 0.0202)
ACTUAL_MASK = Wafer8_Mask


'''
CHOOSE THE PROJECTED AREAS OF JUNCTIONS
SET THE AREAS
'''
## Wafer_5 ##
WITH_DOSETESTS = False
N_SQUIDS        = 10
# IND_J_HIG       = 0.2 ## scince wafer7 we define high of Ljj separately for each
################ _ reset manually same stuff (just avoiding negative width proj _ #####) ##direct length and width


####__from_wafer7__####
# list_of_JJ_seizes       = [ 0.180, 0.180,  0.180,  0.180,
#                             0.240, 0.240,  0.240,  0.240,
#                             0.300, 0.300,  0.300,  0.300 ]
# list_of_L_JJ_seizes     = [ 3.800, 6.300, 6.300, 11.300,   ##took biggest from w6, and x1,x2,x3,x4
#                             3.800, 6.300, 6.300, 11.300,
#                             3.800, 6.300, 6.300, 11.300 ]
# list_of_L_JJ_high_sizes = [ 0.20,  0.20,  0.30,   0.30,
#                             0.20,  0.20,  0.30,   0.30,
#                             0.20,  0.20,  0.30,   0.30  ]

###__wafer7__###
list_of_JJ_seizes       = [ 0.220, 0.220,  0.220,  0.220,
                            0.230, 0.230,  0.230,  0.230,
                            0.245, 0.245,  0.245,  0.245 ]
list_of_L_JJ_seizes     = [ 5.800, 6.000, 6.200, 6.800,
                            5.800, 6.000, 6.200, 6.800,
                            5.800, 6.000, 6.200, 6.800 ]
list_of_L_JJ_high_sizes = [ 0.210, 0.240, 0.280, 0.320,
                            0.210, 0.240, 0.280, 0.320,
                            0.210, 0.240, 0.280, 0.320  ]


#####_________SIZE_OF_small_JJ_Transmon-test_________#####
test_transmon_JJ_size_list = [ 0.220, 0.230, 0.245 ]

################################################################################
################################################################################
#
#
#
#


##_____________________________________________________________###

### real doses ####
dose_dict_structures = {
    1 : 2.2, ##subwires
    2 : 0.6, ##undrct
}

dose_dict_pads = {
    4 : 11,     #wires (arms)
    5 : 10,    #pads for wafer2
    6 : 10      #big wires (in test_structures)
}

##_____test structures id's plan_____###
id_test   =  26
id_trans1 =  1
id_trans2 =  6
id_trans3 = 25
id_NoHeart_Circle = 30
id_NoHeart_Stripl = 29

####______________####

##_______Main_Sizes_______##
##################################################################
NUM_CHIPS_WAF_HORIZ = 6
NUM_CHIPS_WAF_VERTC = 5

CHIP_X_SIZE = 5e3       ##[um]
CHIP_Y_SIZE = 6.8e3     ##[um]

WAFER_DIAMETER = 50.8e3 #um
WAFER_EDGE = 15e3 #um
ORIGIN_LEFT_LOW_MARKER = ( -( 400 +3*CHIP_X_SIZE - 0.5*WAFER_EDGE) +400, 7000) ##coordinates of marker from left low wafer corner


HEART_SIZE =(40,100)

SUB_WIRS_WIDTH = 0.35
SUB_WIRE_HORIZONTAL_TO_JJs = 0.45
WIRS_UNDCT_WIDTH = 0.70
Y_JJ_UNDCT_WIDTH = 0.67
L_JJ_UNDCT_WIDTH = 0.73
SUB_WIRS_PENET = 2 ## make wires penetrate other shapes, to have for sure galvanic contact [um]
WIDWIRE_STPL_BIG = 20
WIDWIRE_STPL_MED = 4
WIDWIRE_CIRC_MED = 8
WIDWIRE_TRANSMON = 8
WIDWIRE_TEST_CHAINS = 8

## Holes
HOLES_RAD = 10
HOLES_DIST = 40
NO_HOLES_BORDER = 15

## Test structures 4-probe ##
TEST_WR2_WID = 15 ## big wires of 4-probe test structures width
# GAPS_FROM_THE_EDGE = (237, -26)
GAPS_FROM_THE_EDGE = (420, -80) ## avoid to overlap markers
SIZE_PADS = (140,120)

## HEARTLESS SAMPLES ##
WIDTH_OF_SHORTCUT_L_WIRE = 1.0 ##[um]

### Technical details for heart ###
gap_X_JJ_to_L           = 1.4  #1.8    ### X-shift both JJ from edge of Inductance #0.4
horiz_wire_shift        = 1.0          ### Y-shift between Inductance and horis wire
gap_JJ_to_horiz_wire    = 8.0          ### Y-shift how far JJ is shifted from the edge

### inductance ###
squid_Loop_hig  = 1.8
squid_Loop_wid  = 2.0 #1.3 ##wafer5 #1.9 ##allprev
shapo_gap_hig   = 0.8
squids_shift_X  = 0.27

### single junctions ###
error_X_move    = 1.2 ##[um]
# deducted_x_wire_sift = 0.4
# deducted_x_wire_sift = 0.42
deducted_x_wire_sift = 0 ## crossection of Yjj - size of bridge
JJ_coupler_heigh = 0.5 ##trapezoid high
cap_wid = 0.05 #shapo undercut

## Test transmon ## (all three will be same)
(Transmon_pad_SX, Transmon_pad_SX) = (500,500)
Transmon_Gap = 200


##################################################################
##################################################################

##_______Layers_______##
##################################################################
LAYER_WAFER     = 102
LAYER_MARKS     = 101
LAYER_SUB_WIRES = 1
LAYER_SUB_UNDCT = 2
LAYER_MED_WIRES = 4
LAYER_PADS      = 5
LAYER_BIG_WIRES = 6
LAYER_SIDE_CIRCULAR_PADS = LAYER_PADS
LAYER_LABELS_TEST = LAYER_PADS
LAYER_HOLES_PADS = LAYER_PADS

### DOSE TESTS ###
if WITH_DOSETESTS:
    DOSES_FOR_DOSE_TEST  = [8.25,   8.50,   8.75,   9.00,   9.25,   9.50,   9.75]
    LAYERS_FOR_DOSE_TEST = [121,     122,    123,    124,   125,     126,    127]
else:
    DOSES_FOR_DOSE_TEST  = []
    LAYERS_FOR_DOSE_TEST = []
##################################################################


##################################################################
##--------------------------------------------------------------##
##################################################################
##_______constants_______##
COS45 = np.cos(np.radians(45))
SQRT2 = np.sqrt(2)

##################################################################
################################################################################

#####______________LOG-file___________________________________####
f = open(main_folder + "log_GDS_build_" +name +".txt", "w")
f.write('Wafer name: '+name+'. V-shape transmons by Vladimir Aug 2020 '+'\n')
# f.write('Ic_JJ_projected = ' +str(Ic_want_JJ__nA) +'[nA]\n\n')
f.close()
def add_line_to_log(line):
    f = open(main_folder + "log_GDS_build_" +name +".txt", "a")
    f.write(line +'\n')
    f.close()
add_line_to_log('chip[#]\tType\tJJ_wid_mask[um]\tL_len_mask[um]\tL_hig[um]\tJJ_project_Area[um2]\tLjj_project_Area[um2]\tproj_JJ_w1[um]\tproj_JJ_w2[um]\tproj_Len_big_JJ[um]')
##################################################################

#####______________NAMES_OF_FILES_DICT________________________####
dict_filenames = {}
dict_types = {}

##################################################################

################################################################################
####____________NECESSARY_SHAPE_FUNCTIONS______#################################
################################################################################

#######________Returns points sequences____#####################################

def draw_sector(N, x0,y0, rad, theta1, theta2, direction=+1):
    points = []
    for i in np.linspace(theta1, theta2, N):
        x = x0+ direction* rad * np.cos( i )
        y = y0+ rad * np.sin( i )
        points.append( (x,y) )
    return points

def draw_sector_ring( N, rad1 = 5, rad2=8, theta1=np.pi/2, theta2=np.pi):
    '''
    Making a sector of ring wirh rounded corners
    '''
    points = []
    # small rad sector
    points += draw_sector(N, 0,0, rad1, theta1, theta2)
    # rounded corner1
    small_rad = (rad2 - rad1)/2.
    mean_rad = np.mean( [rad1,rad2] )
    x0 = mean_rad * np.cos( theta2 )
    y0 = mean_rad * np.sin( theta2 )
    t1 = theta2+np.pi
    t2 = theta2
    points += draw_sector(N, x0,y0, small_rad, t1, t2, direction=+1)
    # big rad sector
    points += draw_sector(N, 0,0, rad2, theta2, theta1)
    # rounded corner2
    x0 = mean_rad * np.cos( theta1 )
    y0 = mean_rad * np.sin( theta1 )
    t1 = -theta1 - np.pi
    t2 = -theta1
    points = points + draw_sector(N, x0,y0, small_rad, t1, t2, direction=-1)
    return points

def draw_sector_wire( N, rad1 = 5, rad2=8, theta1=np.pi/2, theta2=np.pi):
    '''
    Making a sector of ring (simple one)
    '''
    points = []
    # small rad sector
    points += draw_sector(N, 0,0, rad1, theta1, theta2)
    # big rad sector
    points += draw_sector(N, 0,0, rad2, theta2, theta1)
    return points

def draw_trapezoid_coupler_wire(Wid_low, Wid_top, Heigh, posit_low_cen=(0,0)):
    ( x0 , y0 ) = posit_low_cen
    A = ( x0-Wid_low/2,  y0 )
    B = ( x0+Wid_low/2,  y0 )
    C = ( x0+Wid_top/2,  y0+Heigh )
    D = ( x0-Wid_top/2,  y0+Heigh )
    points = [A,B,C,D]
    return points

def draw_trapezoid_inclined_wire(Wid, Len, angle, posit_low_cen=(0,0)):
    sin0 = np.sin( np.radians(angle) )
    cos0 = np.cos( np.radians(angle) )
    base_width = Wid / cos0
    ( x0 , y0 ) = posit_low_cen
    A = ( x0-base_width/2,               y0 )
    B = ( A[0] +base_width, A[1] )
    C = ( B[0] +Len*sin0,   B[1] +Len*cos0 )
    D = ( C[0] -Wid*cos0,   C[1] +Wid*sin0 )
    points = [A,B,C,D]
    return points

def draw_rectangle_round_corners(X_size, Y_size, rad, x0=0, y0=0):
    points = []
    ## set corners
    [ax,ay] = [ x0          , y0          ]
    [bx,by] = [ x0 + X_size , y0          ]
    [cx,cy] = [ x0 + X_size , y0 + Y_size ]
    [dx,dy] = [ x0          , y0 + Y_size ]
    ## to check positions of corners
    if ax > bx:
        [ax,ay, bx,by] = [bx,by, ax,ay]
    if dx > cx:
        [cx,cy, dx,dy] = [dx,dy, cx,cy]
    if ay > dy:
        [ax,ay, dx,dy] = [dx,dy, ax,ay]
    if by > cy:
        [bx,by, cx,cy] = [cx,cy, bx,by]
    ## sector of round corner ( 1/4 pi )
    sec = 0.25
    ###  draw corners one by one ###
    ## A ##
    ca_x = ax + rad
    ca_y = ay + rad
    t_mid = -0.75
    points += draw_sector(50, ca_x,ca_y, rad, (t_mid - sec)*np.pi, (t_mid + sec)*np.pi)
    ## B ##
    cb_x = bx - rad
    cb_y = by + rad
    t_mid = -0.25
    points += draw_sector(50, cb_x,cb_y, rad, (t_mid - sec)*np.pi, (t_mid + sec)*np.pi)
    ## C ##
    cc_x = cx - rad
    cc_y = cy - rad
    t_mid = 0.25
    points += draw_sector(50, cc_x,cc_y, rad, (t_mid - sec)*np.pi, (t_mid + sec)*np.pi)
    ## D ##
    cd_x = dx + rad
    cd_y = dy - rad
    t_mid = 0.75
    points += draw_sector(50, cd_x,cd_y, rad, (t_mid - sec)*np.pi, (t_mid + sec)*np.pi)
    return points

def draw_turning_wire(Wid_wire=20, Len_wire=140, Angle=45, Wid_low=80, Hgt_cplr=60, posit_low_cen=(0,0)):
    # def draw_trapezoid_coupler_wire(Wid_low, Wid_top, Heigh, posit_low_cen=(0,0)):
    #     ( x0 , y0 ) = posit_low_cen
    #     A = ( x0-Wid_low/2,  y0 )
    #     B = ( x0+Wid_low/2,  y0 )
    #     C = ( x0+Wid_top/2,  y0+Heigh )
    #     D = ( x0-Wid_top/2,  y0+Heigh )
    #     points = [A,B,C,D]
    #     return points
    # def draw_trapezoid_inclined_wire(Wid, Len, angle, posit_low_cen=(0,0)):
    #     sin0 = np.sin( np.radians(angle) )
    #     cos0 = np.cos( np.radians(angle) )
    #     base_width = Wid / cos0
    #     ( x0 , y0 ) = posit_low_cen
    #     A = ( x0-base_width/2,               y0 )
    #     B = ( A[0] +base_width, A[1] )
    #     C = ( B[0] +Len*sin0,   B[1] +Len*cos0 )
    #     D = ( C[0] -Wid*cos0,   C[1] +Wid*sin0 )
    #     points = [A,B,C,D]
    #     return points
    # ########
    ( x0 , y0 ) = posit_low_cen
    [F,C,D,E] = draw_trapezoid_inclined_wire(Wid_wire, Len_wire, Angle, posit_low_cen=( x0 , y0+Hgt_cplr ))
    Wid_top = C[0] - F[0]
    [A,B,C,F] = draw_trapezoid_coupler_wire(Wid_low, Wid_top, Hgt_cplr, posit_low_cen=( x0 , y0 ))
    return [A,B,C,D,E,F]

def segment_crossection_point(AB, CD):
    '''
    Find and return a crossection point between segments AB and CD
    each segment is a tuple of two points, each point is a tuple of X,Y coordinates
    '''
    def get_line_equation_from_two_points(A,B):
        (xa, ya) = A
        (xb, yb) = B
        if xa == xb:
            k = np.inf
        else:
            k = (yb-ya) / (xb-xa)
        b = ya - k*xa
        return (k,b)

    (A, B) = AB
    (C, D) = CD
    ## if one of segments is a point - no crossection
    if A==B or C==D:
        return None
    ## if segments own same point - no crossection
    if A==C or B==C or A==D or B==D:
        return None

    ### if both are not vertical
    is_between = lambda x, (a, b): (a <= x <= b) or (a >= x >= b)
    is_vertical = lambda AB: AB[0][0] == AB[1][0]
    if not is_vertical(AB) and not is_vertical(CD):
        (k1, b1) = get_line_equation_from_two_points(A,B)
        (k2, b2) = get_line_equation_from_two_points(C,D)
        ### find a lines crossection point (if exist)
        if k1 != k2:
            zx = (b2-b1) / (k1 - k2)
            zy = k1*zx + b1

            ### if this point belongs to both line segments
            if is_between(zx, (A[0],B[0]))  and  is_between(zx, (C[0],D[0])):
                return (zx,zy)
            else:
                return None
        else:
            # print 'lines are parallel'
            return None

    ## if both are vertical - no crossection
    elif is_vertical(AB) and is_vertical(CD):
        return None

    ## one is vertical - check if it is crossection
    else:
        if is_vertical(AB):
            (vertical, normal) = (AB, CD)
        else:
            (vertical, normal) = (CD, AB)
        (k, b) = get_line_equation_from_two_points(normal[0],normal[1])
        zx = vertical[0][0]
        zy = k*zx + b
        if is_between(zy, (vertical[0][1], vertical[1][1])) and is_between(zx, (normal[0][0],normal[1][0])):
            return (zx,zy)
        else:
            return None


#######________Returns sahpes______________#####################################

gen_rectangle = lambda (x0,y0), (SX,SY), layr: shapes.Rectangle( (x0,y0), (x0+SX,y0+SY), layer=layr)


#######________Returns cells_______________#####################################
def gen_arrow_cell(layer, origin=(0,0), direction='down', size=1.0):
    '''
    Making an arrow of given direction
    ( direction = 'up', 'down', 'left', 'right',  'up-left', 'up-right', 'down-left', 'down-right' )
    size = 200 x 250 um
    returns cell
    '''
    direction = str(direction)
    if direction not in ['up', 'down', 'left', 'right',  'up-left', 'up-right', 'down-left', 'down-right']:
        print('Error of gen_arrow_cell() wrong parameter direction')
        return None

    points = [ (0,0),(-100,100),(-100,150),(-30,100),(-30,250),(30,250),(30,100),(100,150),(100,100) ]

    points_sized = []
    for p in points:
        p_sized = (p[0]*size, p[1]*size)
        points_sized.append(p_sized)

    bdy = core.Boundary(points_sized, layer=layer)
    if direction=='down':
        pass
    elif direction=='up':
        bdy = bdy.rotate(180)
    elif direction=='left':
        bdy = bdy.rotate(-90)
    elif direction=='right':
        bdy = bdy.rotate(90)

    elif direction=='up-left':
        bdy = bdy.rotate(180+45)
    elif direction=='up-right':
        bdy = bdy.rotate(180-45)
    elif direction=='down-left':
        bdy = bdy.rotate(-45)
    elif direction=='down-right':
        bdy = bdy.rotate(45)


    bdy.translate(origin)

    cell = core.Cell('Arrow_'+direction)
    cell.add(bdy)
    return cell

def gen_rect_pad_cell(x0, y0, xsize, ysize, rad, name='Pad', layer=LAYER_PADS):
    '''
    returns cell with rectangle with rounded corners
    '''
    cell = core.Cell('Pad')
    points = draw_rectangle_round_corners(xsize, ysize, rad, x0=x0, y0=y0)
    bdy = core.Boundary(points, layer=layer)
    cell.add(bdy)
    return cell

def gen_circle_pad_cell(x0, y0, rad):
    disk = shapes.Disk((x0,y0), rad, layer=LAYER_PADS)
    cell = core.Cell('Pad')
    cell.add(disk)
    return cell

def gen_arc_pads_cell(R_int, R_out, Sector, PushSect, X_c,Y_c):
    '''
    Sector= 2 == full circle, 1 == half
    PushSect - rotation angle
    '''
    cell = core.Cell('Pad')
    N_sectors = 30
    points = draw_sector_ring(N_sectors, rad1=R_int, rad2=R_out, theta1=np.pi*(1-Sector)/2, theta2=np.pi*(1+Sector)/2)
    Pad1 = core.Boundary(points, layer=LAYER_SIDE_CIRCULAR_PADS)
    Pad3 = Pad1.copy()

    Pad1.rotate( PushSect )
    Pad3.rotate( -PushSect )

    Pad3.rotate(180).translate((X_c,Y_c))
    Pad1.translate( (X_c,Y_c) )
    cell.add(Pad1)
    cell.add(Pad3)
    return cell

def gen_holes_net_in_RoundedCornerRectangle( (X0, Y0), (SIZE_X, SIZE_Y), round_rad, rad=HOLES_RAD, dist=HOLES_DIST, layer=LAYER_HOLES_PADS ):
    '''
    Check if point is in rounded corners rectangle
    (X0, Y0) - coordinates of left low point of rectangle
    round_rad - redius of rounded corner
    '''
    point_is_in_rectagle = lambda (x,y), (x0,y0,x1,y1): x0 <= x <= x1 and y0 <= y <= y1
    point_is_in_circle   = lambda (x,y), (x_c,y_c, circl_rad): (x-x_c)**2 +(y-y_c)**2 <= circl_rad**2

    ## define corner squares and circles
    R = round_rad
    (X1, Y1) = (X0+SIZE_X, Y0+SIZE_Y)
    [sqre1, sqre2, sqre3, sqre4] = [ (X0, Y0, X0+R, Y0+R), (X0, Y1-R, X0+R, Y1), (X1-R, Y1-R, X1, Y1), (X1-R, Y0, X1, Y0+R) ]
    [circ1, circ2, circ3, circ4] = [ (X0+R,Y0+R,R), (X0+R,Y1-R,R), (X1-R,Y1-R,R), (X1-R,Y0+R,R) ]

    def point_is_in_RoundedCornerRectangle( p ):
        '''
        return False IF a) point is outside a big square; b)point in inside corner square AND NOT inside appropriate corner circle
        otherwise returns True
        '''
        if not point_is_in_rectagle( p, (X0,Y0,X1,Y1) ):
            return False

        if point_is_in_rectagle(p,sqre1) and not point_is_in_circle(p,circ1):
            return False
        if point_is_in_rectagle(p,sqre2) and not point_is_in_circle(p,circ2):
            return False
        if point_is_in_rectagle(p,sqre3) and not point_is_in_circle(p,circ3):
            return False
        if point_is_in_rectagle(p,sqre4) and not point_is_in_circle(p,circ4):
            return False
        return True

    ###############################
    dist_y = 2*dist*np.sin(np.pi/3)

    (mean_x, mean_y)    = ( X0 +SIZE_X/2, Y0 +SIZE_Y/2 )
    (delta_x, delta_y)  = ( dist, dist_y )
    (remain_x_space, remain_y_space)  = ( SIZE_X %delta_x, SIZE_Y %delta_y )  ## needs to place holes symmetrically
    (n_columns, n_rows) = ( int(SIZE_X //delta_x ), int(SIZE_Y //delta_y) )

    ( X_shift , Y_shift ) = (  dist*np.cos(np.pi/3) , dist*np.sin(np.pi/3)  )
    cell = core.Cell('Holes')
    for j in range( n_rows+1 ):
        for i in range( n_columns+1 ):
            X = X0 +i*delta_x +remain_x_space/2
            Y = Y0 +j*delta_y +remain_y_space/2
            hole_0 = shapes.Disk( (X,Y), rad, layer=layer )
            if point_is_in_RoundedCornerRectangle( (X,Y) ):
                cell.add(hole_0)
            if point_is_in_RoundedCornerRectangle( (X +X_shift, Y +Y_shift) ):
                hole_1 = hole_0.copy()
                hole_1.translate(( X_shift, Y_shift ))
                cell.add(hole_1)
    return cell

def gen_holes_net_in_circle( (X0,Y0), CIRCL_RAD, rad=HOLES_RAD, dist=HOLES_DIST, layer=LAYER_HOLES_PADS ):
    point_is_in_circle = lambda (x,y), (x_c,y_c, circl_rad): (x-x_c)**2 +(y-y_c)**2 < circl_rad**2

    dist_y = 2*dist*np.sin(np.pi/3)
    NX = int( 2*CIRCL_RAD/dist   ) +1
    NY = int( 2*CIRCL_RAD/dist_y ) +1
    ( X_shift , Y_shift ) = (  dist*np.cos(np.pi/3) , dist_y/2  )
    cell = core.Cell('Holes')
    for y in np.arange(  Y0 -dist_y*NY/2, Y0 +dist_y*NY/2, dist_y):
        for x in np.arange(X0 -dist*NX/2, X0 +dist*NX/2, dist):
            hole_0 = shapes.Disk( (x,y), rad, layer=layer )
            if point_is_in_circle( (x,y), (X0,Y0,CIRCL_RAD) ):
                cell.add(hole_0)
            if point_is_in_circle( (x +X_shift, y +Y_shift), (X0,Y0,CIRCL_RAD) ):
                hole_1 = hole_0.copy()
                hole_1.translate(( X_shift, Y_shift ))
                cell.add(hole_1)
    return cell

def gen_holes_net_in_ringsector( (X0,Y0), sect_rad_int, sect_rad_out, sector_angle, rad=HOLES_RAD, dist=HOLES_DIST, layer=LAYER_HOLES_PADS ):
    cell = core.Cell('Holes')
    ## choose radiuses for the rows of holes (simmetrical):
    mean_rad = np.mean([ sect_rad_out, sect_rad_int ])
    delta_rad = np.sqrt( 0.75*dist**2 )
    interval = sect_rad_out - sect_rad_int
    num_of_rows = int( interval // delta_rad )
    remain_rad  = interval % delta_rad
    row_radius_list = []
    for i in range(num_of_rows+1):
        radius = sect_rad_int +i*delta_rad +remain_rad/2
        row_radius_list.append( radius )
    flag = False
    for rad_row in row_radius_list:
        ## plot a given radius circle of holes
        flag = not flag
        theta_step = 2*np.arcsin( dist/2/rad_row )
        theta1 = np.pi*(1-sector_angle)/2
        theta2 = np.pi*(1+sector_angle)/2
        #######################################
        ## top arc
        for theta in np.arange(theta2 -flag*theta_step/2, theta1 -flag*theta_step/2, -theta_step):
            x = X0  + rad_row*np.cos(theta)
            y = Y0  + rad_row*np.sin(theta)
            hole_0 = shapes.Disk( (x,y), rad, layer=layer )
            cell.add(hole_0)
        ## down arc
        theta1 += np.pi
        theta2 += np.pi
        for theta in np.arange(theta1 +flag*theta_step/2, theta2 +flag*theta_step/2, theta_step):
            x = X0  + rad_row*np.cos(theta)
            y = Y0  + rad_row*np.sin(theta)
            hole_0 = shapes.Disk( (x,y), rad, layer=layer )
            cell.add(hole_0)
    return cell

def gen_dosetest_rectangular(x0, y0, xsize, ysize, round_rad, layer):
    hole_brd = HOLES_RAD+NO_HOLES_BORDER
    cell = core.Cell('DosePad')
    cell.add( gen_rect_pad_cell(x0, y0, xsize, ysize, round_rad, layer=layer))
    cell.add( gen_holes_net_in_RoundedCornerRectangle( (x0+hole_brd, y0+hole_brd), (xsize-2*hole_brd, ysize-2*hole_brd), round_rad-hole_brd , layer=layer ) )
    return cell


#####################
####  Examples      #
#####################
def shape_tests():
    '''
    neeeds only to test the above functions
    '''
    points = draw_sector_ring(20, rad1=6, rad2=8, theta1=(-0.3)*np.pi, theta2=(0.3)*np.pi)
    bdy = core.Boundary(points, layer=101)
    bdy.show()

    points = draw_rectangle_round_corners(8, 6, 0.8, x0=0, y0=0)
    bdy = core.Boundary(points, layer=101)
    bdy.show()

    points = draw_turning_wire(Wid_wire=20, Len_wire=140, Angle=-45, Wid_low=80, Hgt_cplr=60)
    bdy = core.Boundary(points, layer=101)
    bdy.reflect('x')
    bdy.reflect('y')
    bdy.show()

    points = draw_turning_wire(Wid_wire=20, Len_wire=140, Angle=45, Wid_low=80, Hgt_cplr=60)
    bdy = core.Boundary(points, layer=101)
    bdy.show()


### Other necessary functions ###
################################################################################
def my_date_yyyymmdd():
    def my_inttostr_twodigits(num):
        num = int(num)
        string = str(num)
        if num < 10:
            string = '0' + string
        return string
    import datetime
    dt = datetime.datetime.now()
    day   = my_inttostr_twodigits( dt.day )
    month = my_inttostr_twodigits( dt.month )
    year  = my_inttostr_twodigits( dt.year )
    return year+month+day

def copy_and_translate_cell( cell_init, (X_shift, Y_shift) ):
    cell = cell_init.copy()
    for el in cell.elements:
        el.translate((X_shift, Y_shift))
    return cell

def copy_and_rotate_cell( cell_init, angle, center=(0,0)):
    cell = cell_init.copy()
    for el in cell.elements:
        el.rotate(angle, center=center)
    return cell

def copy_and_reflect_cell( cell_init, axis='x' ):
    cell = cell_init.copy()
    bb = cell.bounding_box
    (low, top, right, left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
    X_c = (right +left)/2
    Y_c = (top +low)/2

    for el in cell.elements:
        el.reflect(axis, origin=(X_c,Y_c))
    return cell

def get_chip_position_from_number( num_of_chip ):
    '''
    num_of_chip = { 1, 2, <...>, 30 } (for 5x6 wafer)
    update 2020-08-26 {left-low must be (0,0)}
    '''
    ### checking ###
    num_of_chip = abs(int(num_of_chip))
    if num_of_chip<=0:
        print('\n!____num_of_chip must be > 0')
        return (None, None)
    print('num_of_chip=',num_of_chip)

    num_of_chip -= 1

    ### working ###
        ## chipposition
    X_num = num_of_chip % 6
    Y_num = num_of_chip // 6

        ## coordinates ##
    X =  X_num*CHIP_X_SIZE
    Y = -Y_num*CHIP_Y_SIZE

    Y += (NUM_CHIPS_WAF_VERTC-1)*CHIP_Y_SIZE ## update

    return (X,Y)

################################################################################
################################################################################
####____________GENERATING_CELLS_FUNCTIONS_____#################################
################################################################################

def draw_wafer_2inch():
    Wafer_edge_cell = core.Cell('Wafer_edge')
    Wafer_Radius   = 25.4 *1e3 #[um]
    Wafer_center_x = 3*CHIP_X_SIZE
    Wafer_center_y = 2.5*CHIP_Y_SIZE
    wafer_base_len = 19.1 *1e3 #[um]
    wafer_theta_deg = np.degrees( np.arcsin(0.5*wafer_base_len/Wafer_Radius) )
    Wafer_edge_cell.add(shapes.Circle((Wafer_center_x,Wafer_center_y), Wafer_Radius, 1, layer=0, initial_angle=270-wafer_theta_deg, final_angle=-90+wafer_theta_deg) )
    Wafer_edge_cell = copy_and_translate_cell( Wafer_edge_cell, (-5450, 6536) )
    return Wafer_edge_cell

def gen_wafer_universal( x0=0,y0=0, wafer_diameter=WAFER_DIAMETER, edge_length=WAFER_EDGE, LAYER=LAYER_WAFER):
    '''
    Generating cell with wafer of given diameter and given edge length
    takes:
        coordinates of left low corner,
        diameter,
        legth of edge
        and layout layer
        ( all dimensions are in [um])
    returns:
        cell with wafer layout
    '''
    ##_calc_vars__##
    R_WAFER = wafer_diameter/2
    xc = x0 +0.5*edge_length
    yc = y0 +np.sqrt(R_WAFER**2 -(0.5*edge_length)**2)
    alpha = np.arcsin( (0.5*edge_length)/R_WAFER )
    fraction = 1.0 -( 2*alpha / (2*np.pi) )
    ###########_define_technical_parameters_####
    lw = 0.002*R_WAFER #um
    N_sectors = 100
    [R_int, R_out, Sector] = [R_WAFER-lw, R_WAFER, 2*fraction]
    #############_do_the_things_###
    cell = core.Cell('wafer')
    points = draw_sector_ring(N_sectors, rad1=R_int, rad2=R_out, theta1=np.pi*(1-Sector)/2, theta2=np.pi*(1+Sector)/2)
    sector1 = core.Boundary(points, layer=LAYER)
    cell.add(sector1)
    cell = copy_and_translate_cell(cell, (xc,yc))
    cell.add(gen_rectangle( (x0,y0), (edge_length,lw), LAYER))
    return cell

def gen_wafer_markers(sqr_mk_size=8):
    Wafer_marks = core.Cell('Wafer_marks')
    Wafer_marks.add( gen_rectangle( (0 -400, 0 -200),        (-sqr_mk_size,-sqr_mk_size),  LAYER_MARKS) ) #left low
    Wafer_marks.add( gen_rectangle( (0 -400, 34000 +200),    (-sqr_mk_size, sqr_mk_size),  LAYER_MARKS) ) #left top
    Wafer_marks.add( gen_rectangle( (30000 +400, 0 -200),    ( sqr_mk_size,-sqr_mk_size),  LAYER_MARKS) ) #rigth low
    Wafer_marks.add( gen_rectangle( (30000 +400, 34000 +200),( sqr_mk_size, sqr_mk_size),  LAYER_MARKS) ) #rigth top

    #######__arrows to markers__#######
    gap_btw_ar = 200
    ar_length = 250
    ### vertical traces ###
    for i in range(5):
        # left low marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='up', origin=(-400-sqr_mk_size/2,  -200-sqr_mk_size -(i+1)*gap_btw_ar -i*ar_length )))
        # right low marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='up', origin=(30400+sqr_mk_size/2, -200-sqr_mk_size -(i+1)*gap_btw_ar -i*ar_length )))
    for i in range(6):
        # left up marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', origin=(-400-sqr_mk_size/2,  34000+sqr_mk_size +(i+2)*gap_btw_ar +(i+0)*ar_length )))
        # right up marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', origin=(30400+sqr_mk_size/2, 34000+sqr_mk_size +(i+2)*gap_btw_ar +(i+0)*ar_length )))
    ### horizontal traces ###
    for i in range(6):
        # left low marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='right', origin=(-400-sqr_mk_size/2 -(i+1)*gap_btw_ar -i*ar_length,  -200-sqr_mk_size/2 )))
        # right low marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='left',  origin=(30400+sqr_mk_size/2 +(i+1)*gap_btw_ar +i*ar_length,  -200-sqr_mk_size/2 )))
    for i in range(6):
        # left up marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='right', origin=(-400-sqr_mk_size/2 -(i+1)*gap_btw_ar -i*ar_length,   34200+sqr_mk_size/2 )))
        # right up marker
        Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='left',  origin=(30400+sqr_mk_size/2 +(i+1)*gap_btw_ar +i*ar_length,  34200+sqr_mk_size/2 )))

    #######__arrows_aligment_wafer_big__#######
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='up', size=2.5, origin=(15000, 42500 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='up', size=2.5, origin=(15000, 42500-700 )))

    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800+800 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800+800*2 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800+800*3 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800+800*4 )))
    Wafer_marks.add(gen_arrow_cell(LAYER_MARKS, direction='down', size=2.5, origin=(15000, -6800+800*5 )))

    return Wafer_marks

def gen_mark_cell( chip_size_x, chip_size_y, size=25, position=(0,0), sqr_mk_size=8, sqr_mk_dist=200 ):
    (x0, y0) = position
    ####################################
    ####_CORNER_BIG_MARKERS_############
    ####################################
    cell_bg_mkrs = core.Cell('big markers')
    dicing_marker=[(0,0), (0,3*size), (size,4*size), (size,size), (4*size,size), (3*size,0)]
    dicing_marker_origin = core.Boundary(dicing_marker)
    dicing_marker_origin.layer = LAYER_MARKS

    dicing_marker_LB = dicing_marker_origin.copy()
    dicing_marker_LB.translate((x0,y0))
    cell_bg_mkrs.add(dicing_marker_LB)

    dicing_marker_LT = dicing_marker_origin.copy()
    dicing_marker_LT.rotate(270).translate((0,chip_size_y))
    dicing_marker_LT.translate((x0,y0))
    cell_bg_mkrs.add(dicing_marker_LT)

    dicing_marker_RT = dicing_marker_origin.copy()
    dicing_marker_RT.rotate(180).translate((chip_size_x,chip_size_y))
    dicing_marker_RT.translate((x0,y0))
    cell_bg_mkrs.add(dicing_marker_RT)

    dicing_marker_RB = dicing_marker_origin.copy()
    dicing_marker_RB.rotate(90).translate((chip_size_x,0))
    dicing_marker_RB.translate((x0,y0))
    cell_bg_mkrs.add(dicing_marker_RB)

    ####################################
    ####_SMALL_SQUARE MARKERS_##########
    ####################################
    ####__LEFT_LOW_CORNER___####
    cell_corner_left_low =core.Cell('corner_small_markers')
    orig_marker = gen_rectangle( (0 -sqr_mk_size/2, 0 -sqr_mk_size/2), (sqr_mk_size, sqr_mk_size), LAYER_MARKS)
    ###############
    sm_marker_corner = orig_marker.copy()
    sm_marker_corner.translate(( sqr_mk_dist, sqr_mk_dist ))
    cell_corner_left_low.add(sm_marker_corner)
    ###
    sm_marker_upper = sm_marker_corner.copy()
    sm_marker_upper.translate(( 0, sqr_mk_dist ))
    cell_corner_left_low.add(sm_marker_upper)
    ##
    sm_marker_uppest = sm_marker_upper.copy()
    sm_marker_uppest.translate(( 0, sqr_mk_dist ))
    cell_corner_left_low.add(sm_marker_uppest)
    ###
    sm_marker_righter = sm_marker_corner.copy()
    sm_marker_righter.translate(( sqr_mk_dist, 0 ))
    cell_corner_left_low.add(sm_marker_righter)
    ##
    sm_marker_rightest = sm_marker_righter.copy()
    sm_marker_rightest.translate(( sqr_mk_dist, 0 ))
    cell_corner_left_low.add(sm_marker_rightest)
    ###############
    cell_corner_right_low = copy_and_translate_cell( copy_and_rotate_cell(cell_corner_left_low,  90, center=(sqr_mk_dist, sqr_mk_dist)), (chip_size_x-2*sqr_mk_dist, 0                         ) )
    cell_corner_left_top  = copy_and_translate_cell( copy_and_rotate_cell(cell_corner_left_low, 270, center=(sqr_mk_dist, sqr_mk_dist)), (0                        , chip_size_y-2*sqr_mk_dist ) )
    cell_corner_right_top = copy_and_translate_cell( copy_and_rotate_cell(cell_corner_left_low, 180, center=(sqr_mk_dist, sqr_mk_dist)), (chip_size_x-2*sqr_mk_dist, chip_size_y-2*sqr_mk_dist ) )

    cell_sm_mkrs = core.Cell('small markers')
    cell_sm_mkrs.add(cell_corner_left_low)
    cell_sm_mkrs.add(cell_corner_right_low)
    cell_sm_mkrs.add(cell_corner_left_top)
    cell_sm_mkrs.add(cell_corner_right_top)
    ####################################
    ####################################
    cell_markers=core.Cell('Marker')
    cell_markers.add(cell_bg_mkrs)
    cell_markers.add(cell_sm_mkrs)
    return cell_markers

def gen_strip_pads_cell(param_stripline):
    '''
     __   _______   __
    |  | |       | |  |
    |__| |_______| |__|
    '''
    [Gap, MidLen, SideLen, MidWid, SideWid, Shift, RoundRad, Strc_Gap, X_c, Y_c] = param_stripline

    Pad1 = gen_rect_pad_cell( X_c-SideWid/2,       Y_c+MidLen/2+Gap, SideWid, SideLen, RoundRad, name='PadA' )
    Pad2 = gen_rect_pad_cell( X_c-MidWid/2+Shift,  Y_c-MidLen/2,  MidWid,  MidLen,  RoundRad, name='PadB' )
    Pad3 = gen_rect_pad_cell( X_c-SideWid/2,       Y_c-MidLen/2-SideLen-Gap, SideWid, SideLen, RoundRad, name='PadC' )
    cell = core.Cell('Pads')
    cell.add(Pad1)
    cell.add(Pad2)
    cell.add(Pad3)
    return cell

def gen_round_pads_cell(param_circular):
    '''
    (O)
    '''
    [R3, Gap, RingWidth, Sector, PushSect, X_c, Y_c] = param_circular
    Pad2   = gen_circle_pad_cell(X_c, Y_c, R3-RingWidth-Gap)
    Pad1_3 = gen_arc_pads_cell(R3-RingWidth, R3, Sector, PushSect, X_c,Y_c)
    cell = core.Cell('Pads')
    cell.add(Pad2)
    cell.add(Pad1_3)
    return cell

################################################################################

def gen_wires_for_stripline_cell(param_stripline, w_big=WIDWIRE_STPL_BIG, w_med=WIDWIRE_STPL_MED, l_med=100, heart_size=(40,100)):
    [Gap, MidLen, SideLen, MidWid, SideWid, Shift, RoundRad, Strc_Gap, X_c, Y_c] = param_stripline
    (Heart_SX, Heart_SY) = heart_size

    ### Calc Geometry ###
    wire2_big_len = Strc_Gap -l_med -Heart_SX +Shift
    wire2_big_strartX = X_c -MidWid/2 +Shift
    wire2_big_strartY = Y_c -w_big/2

    wire2_med_strartX = wire2_big_strartX -wire2_big_len
    wire2_med_strartY = wire2_big_strartY +w_big/2 -w_med/2

    wire1_rad_out = Strc_Gap +w_big/2
    wire1_rad_in  = wire1_rad_out -w_big
    wire1_rnd_startX = X_c -SideWid/2
    wire1_rnd_startY = Y_c +MidLen/2 +Gap +RoundRad +w_big -wire1_rad_out

    wire1_strhgt_startX = X_c -MidWid/2 -wire1_rad_out
    wire1_strhgt_startY = Y_c +l_med +Heart_SY/2
    wire1_strght_len = wire1_rnd_startY -Y_c -l_med -Heart_SY/2

    wire1_med_startX = X_c -MidWid/2 -wire1_rad_out +w_big/2 -w_med/2
    wire1_med_startY = Y_c +Heart_SY/2
    #####################

    ### Draw objects ###
    ### Wire2 ###
    wire2_big = shapes.Rectangle((0,0),(-wire2_big_len, w_big), layer=LAYER_BIG_WIRES)
    wire2_big.translate(( wire2_big_strartX, wire2_big_strartY ))
    Wire2_big_cell = core.Cell('Wire2_big')
    Wire2_big_cell.add(wire2_big)

    wire2_med = shapes.Rectangle((0,0), ( -l_med, w_med ), layer=LAYER_MED_WIRES)
    wire2_med.translate(( wire2_med_strartX, wire2_med_strartY ))
    Wire2_med_cell = core.Cell('Wire2_med')
    Wire2_med_cell.add(wire2_med)

    Wire2_cell = core.Cell('Wire2')
    Wire2_cell.add( Wire2_big_cell )
    Wire2_cell.add( Wire2_med_cell )

    ### Wire1 ###
    wire1_rnd = shapes.Disk((wire1_rnd_startX, wire1_rnd_startY), wire1_rad_out, inner_radius=wire1_rad_in, initial_angle=90, final_angle=180, layer=LAYER_BIG_WIRES)
    Wire1_rnd_cell = core.Cell('Wire1_rnd')
    Wire1_rnd_cell.add( wire1_rnd )

    wire1_strght = shapes.Rectangle( ( 0,0 ), ( w_big, wire1_strght_len ), layer=LAYER_BIG_WIRES)
    wire1_strght.translate(( wire1_strhgt_startX, wire1_strhgt_startY ))
    Wire1_strght_cell = core.Cell('Wire1_strght')
    Wire1_strght_cell.add( wire1_strght )

    wire1_med = shapes.Rectangle( (0,0), (w_med,l_med), layer=LAYER_MED_WIRES)
    wire1_med.translate(( wire1_med_startX, wire1_med_startY ))
    Wire1_med = core.Cell('Wire1_med')
    Wire1_med.add( wire1_med )

    Wire1_cell = core.Cell('Wire1')
    Wire1_cell.add( Wire1_rnd_cell )
    Wire1_cell.add( Wire1_strght_cell )
    Wire1_cell.add( Wire1_med )

    ### Wire3 ###
    wire3_rnd = wire1_rnd.copy()
    wire3_rnd = wire3_rnd.reflect('x', origin=(X_c,Y_c))
    Wire3_rnd_cell = core.Cell('Wire3_rnd')
    Wire3_rnd_cell.add( wire3_rnd )

    wire3_strght = wire1_strght.copy()
    wire3_strght = wire3_strght.reflect('x', origin=(X_c,Y_c))
    Wire3_strght_cell = core.Cell('Wire3_strght')
    Wire3_strght_cell.add( wire3_strght )

    wire3_med = wire1_med.copy()
    wire3_med.reflect('x', origin=(X_c,Y_c))
    Wire3_med = core.Cell('Wire3_med')
    Wire3_med.add( wire3_med )

    Wire3_cell = core.Cell('Wire3')
    Wire3_cell.add( Wire3_rnd_cell )
    Wire3_cell.add( Wire3_strght_cell )
    Wire3_cell.add( Wire3_med )

    ### Assemble all ###
    cell = core.Cell('Wires')
    cell.add(Wire1_cell)
    cell.add(Wire2_cell)
    cell.add(Wire3_cell)
    return cell

def gen_wires_for_circular_cell(param_circular, w_big=WIDWIRE_CIRC_MED, heart_size=(40,100), heart_position=(0,0)):
    '''
    Make wiring for circular design
    '''
    ######__parameters___######
    N_sectors = 30
    overlap_big = 2

    [R3, Gap, RingWidth, Sector, PushSect, X_c, Y_c] = param_circular
    R2 = R3 - RingWidth
    R1 = R2 - Gap
    R_wire_mid = np.mean([ R2, R3 ])
    R_int = R_wire_mid - w_big/2
    R_out = R_wire_mid + w_big/2
    [heart_SX, heart_SY] = HEART_SIZE
    (heart_X, heart_Y) = heart_position

    ### Making precice wires overlap with Pad1 & Pad3 for any wire_width ###
    ###____###
    rad_small_round = (R3 - R2)/2
    len_add = np.sqrt( rad_small_round**2 -(w_big/2)**2 )  ### linear approx[um]
    theta_add = np.arcsin( ( len_add - overlap_big ) / ( R_wire_mid+w_big/2 ) )  ### positive in [rad]
    ###____###
    theta1 =np.pi*(1+Sector)/2 + theta_add ##end of Pad
    theta2 = np.arccos( heart_SY/(2*R_wire_mid) ) + np.pi/2 ##end of Heart


    l_middle_wire = abs(heart_X -(X_c -R1) ) -heart_SX +overlap_big

    ######__draw_objects__#######
    points = draw_sector_wire(N_sectors, rad1=R_int, rad2=R_out, theta1=theta1, theta2=theta2)
    wire1 = core.Boundary(points, layer=LAYER_MED_WIRES)
    wire1.translate((X_c, Y_c))
    Wire1_cell = core.Cell('Wire1')
    Wire1_cell.add(wire1)

    wire3 = wire1.copy()
    wire3.reflect('x', origin=(X_c,Y_c))
    Wire3_cell = core.Cell('Wire3')
    Wire3_cell.add(wire3)

    wire2 = shapes.Rectangle( (0,0), (-l_middle_wire, w_big), layer=LAYER_MED_WIRES)
    wire2.translate((X_c -R1 +overlap_big, Y_c -w_big/2))
    Wire2_cell = core.Cell('Wire2')
    Wire2_cell.add(wire2)

    ######__assembling__#######
    cell = core.Cell('Wires')
    cell.add(Wire1_cell)
    cell.add(Wire2_cell)
    cell.add(Wire3_cell)
    return cell

################################################################################
def gen_junction( position_low_wire, JJ_wid, error_X_move=1.2, deducted_x_wire_sift=0.42, JJ_coupler_heigh=0.5, cap_wid=0.05, upsdown=False ):
    '''
    Making Y-junctions
    Takes position of center of the wire which is low and left (connected to L)
    Area of junction (make it square)
    and shift of wire after entering JJ
    '''
    sub_wire_wid = SUB_WIRS_WIDTH
    sub_wire_undct_wid = Y_JJ_UNDCT_WIDTH
    #
    if deducted_x_wire_sift >= sub_wire_undct_wid:
        print('\n!____deducted_x_wire_sift   must be less than   sub_wire_undct_wid. Set to 0')
        deducted_x_wire_sift = 0
    elif deducted_x_wire_sift < 0:
        print('\n!____deducted_x_wire_sift must be positive or zero. Set to 0')
        deducted_x_wire_sift = 0

    ## values calc ##
    limit_X_move = error_X_move / 2
    add_tail = limit_X_move / ( 2*SQRT2 )
    # print '\n\n\n\n\n\n\n\n\nadd_tail=', add_tail

    Y_ovrlp_JJ = SQRT2*JJ_wid + 2*add_tail*COS45
    J_SY = 2*JJ_coupler_heigh + Y_ovrlp_JJ
    (JJ1_X, JJ1_y) = position_low_wire
    JJ_len_w = Y_ovrlp_JJ*COS45 +add_tail
    wire_shift = SQRT2*JJ_wid  +sub_wire_undct_wid  +SQRT2*add_tail
    wire_shift -= deducted_x_wire_sift
    #####################

    ## Low part of low JJ
    points = draw_turning_wire(Wid_wire=JJ_wid, Len_wire=JJ_len_w,
                            Wid_low=sub_wire_wid, Hgt_cplr=JJ_coupler_heigh)
    Low_JJ_bdy = core.Boundary(points, layer=LAYER_SUB_WIRES)
    Low_JJ_bdy.translate(( JJ1_X, JJ1_y ))

    ## Low part of low JJ - undercut
    edge_pnts = points[1:4]
    undct_pnts = edge_pnts[::-1]
    for i in range(len(edge_pnts)):
        undct_pnts.append(( edge_pnts[i][0] + sub_wire_undct_wid,  edge_pnts[i][1] ))
    undct_pnts[-1] = ( undct_pnts[-1][0] - 0.5*sub_wire_undct_wid ,
                        undct_pnts[-1][1] - 0.5*sub_wire_undct_wid )
    Low_JJ_undct_bdy = core.Boundary( undct_pnts, layer=LAYER_SUB_UNDCT )
    Low_JJ_undct_bdy.translate(( JJ1_X, JJ1_y ))
    [A_1, B_1, C_1, D_1, E_1, F_1] = undct_pnts

    ## + (( JJ1_X, JJ1_y ))

    ## Top part of low JJ
    Top_JJ_bdy = Low_JJ_bdy.copy()
    Top_JJ_bdy.reflect('x', origin=(JJ1_X, JJ1_y))
    Top_JJ_bdy.translate(( wire_shift, J_SY ))

    ## Top part of low JJ - undercut
    edge_pnts = Top_JJ_bdy.points[4:7].tolist()
    undct_pnts = edge_pnts[::-1]
    for i in range(len(edge_pnts)):
        undct_pnts.append( (edge_pnts[i][0] - sub_wire_undct_wid,  edge_pnts[i][1]) )
    ## Top undercut cap on low JJ
       ###avoiding the overlap
    [A,B,C,D,E,F] = undct_pnts

    EX_shift = deducted_x_wire_sift
    C   = ( C[0], C[1] )
    D2  = ( E[0]+EX_shift, E[1] ) #done
    D1  = ( F_1[0]+JJ1_X, F_1[1]+JJ1_y ) #same as another undercut
    D0  = ( D1[0]-EX_shift/2, D1[1]-EX_shift/2 )
    D3  = ( E[0]+EX_shift/2, E[1]-EX_shift/2 )

    Z_c = segment_crossection_point( (D1,D0), (D,C) )
    if Z_c is not None:         #### exception for short shifts
        #print '______________CROSSECTION!!!!_______'
        undct_pnts    = [A, B, C, Z_c, D1, D2, D3, E, F] ### and even this way (full points clamp)
        (top_x_left, top_y_left) = Z_c
    else:
        #print '_____NOOOO_CROSSECTION!!!!_______'
        undct_pnts    = [A, B, C, D, D0, D1, D2, D3, E, F] ### and even this way (full points clamp) #as before
        (top_x_left, top_y_left) = D

    top_x_right  = C[0]
    top_y_right  = C[1] -cap_wid
    Top_JJ_undct_bdy = core.Boundary(undct_pnts, layer=LAYER_SUB_UNDCT)
    Top_JJ_undct_cap = shapes.Rectangle( (top_x_left, top_y_left), (top_x_right, top_y_right), layer=LAYER_SUB_UNDCT )


    #####################
    if upsdown:
        Low_JJ_bdy.reflect('x', origin=(JJ1_X, JJ1_y) )
        Low_JJ_undct_bdy.reflect('x', origin=(JJ1_X, JJ1_y) )
        Top_JJ_bdy.reflect('x', origin=(JJ1_X, JJ1_y) )
        Top_JJ_undct_bdy.reflect('x', origin=(JJ1_X, JJ1_y) )
        Top_JJ_undct_cap.reflect('x', origin=(JJ1_X, JJ1_y) )
    ## Assembling cell ##
    cell = core.Cell('JJ')
    cell.add(Low_JJ_bdy)
    cell.add(Low_JJ_undct_bdy)
    cell.add(Top_JJ_undct_cap)
    cell.add(Top_JJ_bdy)
    cell.add(Top_JJ_undct_bdy)
    param = [wire_shift]
    return [param, cell]

def gen_inductance( position, N_squids, J_wid, J_hig, Loop_hig, Loop_wid, gap_hig, squids_shift_X):
    '''
    Takes central postion of L
    '''

    sub_wire_wid = SUB_WIRS_WIDTH
    sub_wire_undct_wid = WIRS_UNDCT_WIDTH
    big_junction_undct = L_JJ_UNDCT_WIDTH
    ### checks ####
    N_squids = abs(int(N_squids))
    if N_squids %2 !=0:
        print '\n!____N_squids must be even number'
        N_squids +=1

    if Loop_wid > (J_wid -2*sub_wire_wid ):
        Loop_wid = J_wid -2*sub_wire_wid

    ### values ####
    (L_x0, L_y0) = position
    L_len = N_squids*(Loop_hig+2*J_hig) + (N_squids-1)*gap_hig
    L_wid = J_wid + 2*sub_wire_undct_wid
    SQ1_x = L_x0 + squids_shift_X
    SQ1_y = L_y0 -L_len/2
    Gap_wr_x = L_x0 -sub_wire_wid/2

    ### Drawing ###
    def gen_squid( position, J_wid, J_hig, Loop_hig, Loop_wid ):
        ###   variables   ###
        Undrct = big_junction_undct
        Undrct_wir = sub_wire_undct_wid
        Wir_wid = sub_wire_wid
        (x0,y0) = position
        DJJ_x = x0 -J_wid/2
        DJJ_y = y0
        UJJ_x = DJJ_x
        UJJ_y = DJJ_y +Loop_hig +J_hig
        LC_x  = x0 -Loop_wid/2 -Wir_wid
        LC_y  = y0 +J_hig
        RC_x  = x0 +Loop_wid/2
        RC_y  = y0 +J_hig

        ####   Shaping   #####
        ### Down JJ (DJJ)###
        DJJ = shapes.Rectangle(        (DJJ_x, DJJ_y),       (DJJ_x+J_wid, DJJ_y+J_hig),       layer=LAYER_SUB_WIRES)
        DJJ_r_udct = shapes.Rectangle( (DJJ_x+J_wid, DJJ_y), (DJJ_x+J_wid+Undrct, DJJ_y+J_hig),layer=LAYER_SUB_UNDCT)
        DJJ_l_udct = shapes.Rectangle( (DJJ_x, DJJ_y),       (DJJ_x-Undrct, DJJ_y+J_hig),      layer=LAYER_SUB_UNDCT )

        ### Up   JJ (UJJ)###
        UJJ = shapes.Rectangle(        (UJJ_x, UJJ_y),       (UJJ_x+J_wid, UJJ_y+J_hig),       layer=LAYER_SUB_WIRES)
        UJJ_r_udct = shapes.Rectangle( (UJJ_x+J_wid, UJJ_y), (UJJ_x+J_wid+Undrct, UJJ_y+J_hig),layer=LAYER_SUB_UNDCT)
        UJJ_l_udct = shapes.Rectangle( (UJJ_x, UJJ_y),       (UJJ_x-Undrct, UJJ_y+J_hig),      layer=LAYER_SUB_UNDCT )

        ### Left Column (LC) ###
        LC = shapes.Rectangle(        ( LC_x, LC_y), ( LC_x+Wir_wid, LC_y+Loop_hig), layer=LAYER_SUB_WIRES )
        LC_l_udct = shapes.Rectangle( ( LC_x, LC_y), ( LC_x-Undrct_wir,  LC_y+Loop_hig), layer=LAYER_SUB_UNDCT )

        ### Right Column (RC)###
        RC = shapes.Rectangle( ( RC_x, RC_y), ( RC_x+Wir_wid, RC_y+Loop_hig), layer=LAYER_SUB_WIRES )
        RC_r_udct = shapes.Rectangle( ( RC_x+Wir_wid, RC_y), ( RC_x+Wir_wid+Undrct_wir,  RC_y+Loop_hig), layer=LAYER_SUB_UNDCT )

        cell = core.Cell('squid')
        cell.add(DJJ)
        cell.add(DJJ_r_udct)
        cell.add(DJJ_l_udct)
        cell.add(UJJ)
        cell.add(UJJ_r_udct)
        cell.add(UJJ_l_udct)
        cell.add(LC)
        cell.add(LC_l_udct)
        cell.add(RC)
        cell.add(RC_r_udct)
        return cell

    def gen_gapwire(position, gap_hig, orient=0):
        (x0,y0) = position
        wr_wid = sub_wire_wid
        un_wid = sub_wire_undct_wid
        wr_x = un_x = x0
        wr_y = un_y = y0
        if   orient==1:
            un_x += wr_wid
        elif orient==0:
            wr_x += un_wid
        wire  = shapes.Rectangle( (wr_x, wr_y), (wr_x+wr_wid, wr_y+gap_hig), layer=LAYER_SUB_WIRES)
        udrct = shapes.Rectangle( (un_x, un_y), (un_x+un_wid, un_y+gap_hig), layer=LAYER_SUB_UNDCT)
        cell = core.Cell('gapwire')
        cell.add(wire)
        cell.add(udrct)
        return cell

    cell = core.Cell('L')
    ##
    for i in range(N_squids):
        (SQx, SQy) = ( SQ1_x, SQ1_y +i*( Loop_hig +2*J_hig +gap_hig ) )
        squid = gen_squid( (SQx, SQy), J_wid, J_hig, Loop_hig, Loop_wid )
        cell.add( squid )
        if i < N_squids-1:
            gap_wr = gen_gapwire( (Gap_wr_x, SQy+Loop_hig+2*J_hig), gap_hig, orient=i%2 )
            cell.add( gap_wr )

    ##############
    bb = cell.bounding_box
    L_wid = bb[1,0] - bb[0,0]
    L_len = bb[1,1] - bb[0,1]
    return cell

def gen_heart_cell(position, JJ_wid_mask=0.2, ind_J_wid_mask=1.0, ind_J_hig_mask=0.2, N_SQUIDS=N_SQUIDS, heart_size=(40,100), HEARTLESS=False, BROKENHEART=False ):
    ### VARIABLES ###
    sub_wire_wid=SUB_WIRS_WIDTH
    sub_wire_undct_wid=WIRS_UNDCT_WIDTH

    ####################################
    ## parameters ##
    (x0, y0) = position ## low left point
    (Heart_SX, Heart_SY) = heart_size

    # JJ_wid_mask     = round( Wmask_from_Area_JJ_Y_junction(mask_params, JJ_area)        , 2)
    # ind_J_wid_mask  = round( Lmask_from_Lproj_Big_junction(mask_params, ind_J_wid_proj) , 2)
    ######################################################################
    ######################################################################
    ###   Making Inductance   ###
    L_position = (x0, y0+Heart_SY/2)
    Inductance_cell = gen_inductance( L_position, N_SQUIDS, ind_J_wid_mask, ind_J_hig_mask, squid_Loop_hig, squid_Loop_wid, shapo_gap_hig, squids_shift_X )
    ## extracting parameters ##
    bb = Inductance_cell.bounding_box
    L_wid = bb[1,0] - bb[0,0]
    L_len = bb[1,1] - bb[0,1]
    L_right_lim = bb[1,0]
    knee_jj_wid = JJ_wid_mask *SQRT2
    JJ_X = L_right_lim +np.max([ knee_jj_wid, sub_wire_wid ])/2 +gap_X_JJ_to_L
    JJ1_pos_wire_c = (JJ_X,  y0 +Heart_SY/2 -L_len/2 -horiz_wire_shift +gap_JJ_to_horiz_wire)
    JJ2_pos_wire_c = (JJ_X,  y0 +Heart_SY/2 +L_len/2 +horiz_wire_shift -gap_JJ_to_horiz_wire)
    horiz_wire_len = JJ_X - x0
    vert_wire_len = (Heart_SY -L_len)/2
    ##            ##

    ######################################################################
    #####__JJ1__####

    [[JJ_wire_shift], JJ_1_cell] = gen_junction( JJ1_pos_wire_c, JJ_wid_mask, error_X_move=error_X_move, deducted_x_wire_sift=deducted_x_wire_sift, JJ_coupler_heigh=JJ_coupler_heigh, cap_wid=cap_wid)
    bb = JJ_1_cell.bounding_box
    JJ1_wid = bb[1,0] - bb[0,0]
    JJ1_len = bb[1,1] - bb[0,1]

    #####__JJ2__####
    [[JJ_wire_shift], JJ_2_cell] = gen_junction( JJ2_pos_wire_c, JJ_wid_mask, error_X_move=error_X_move, deducted_x_wire_sift=deducted_x_wire_sift, JJ_coupler_heigh=JJ_coupler_heigh, cap_wid=cap_wid, upsdown=True)
    bb = JJ_2_cell.bounding_box
    len_btwn_jj_wire = L_len +2*horiz_wire_shift -2*JJ1_len -2*gap_JJ_to_horiz_wire
    len_low_vert_undct = vert_wire_len +SUB_WIRS_PENET -horiz_wire_shift -sub_wire_wid
    len_right_wire = Heart_SX -horiz_wire_len -JJ_wire_shift -sub_wire_wid/2

    ######################################################################
    ###   Making Wiring   ################################################
    ######################################################################
    ######################################################################
    SubWires_cell = core.Cell('SubWires')
    SubUndct_cell = core.Cell('SubWires_Undercut')

    ###################################
    ## Low Part ##
      ##vertical wire from bottom to inductance
    low_vert_wire = shapes.Rectangle( (0,0),(sub_wire_wid, vert_wire_len +SUB_WIRS_PENET), layer=LAYER_SUB_WIRES)
    low_vert_wire.translate(( x0 -sub_wire_wid/2, y0 -SUB_WIRS_PENET ))
    SubWires_cell.add(low_vert_wire)

    low_vert_undct1 = shapes.Rectangle( (0,0), (sub_wire_undct_wid, len_low_vert_undct), layer=LAYER_SUB_UNDCT )
    low_vert_undct1.translate(( x0 +sub_wire_wid/2, y0 -SUB_WIRS_PENET))
    SubUndct_cell.add(low_vert_undct1)

    if not BROKENHEART:
        low_vert_undct2 = shapes.Rectangle( (0,0), (sub_wire_undct_wid, horiz_wire_shift), layer=LAYER_SUB_UNDCT )
        low_vert_undct2.translate(( x0 +sub_wire_wid/2, y0 +len_low_vert_undct -SUB_WIRS_PENET +sub_wire_wid))
    else:
        low_vert_undct2 = shapes.Rectangle( (0,0), (sub_wire_undct_wid, horiz_wire_shift+sub_wire_wid), layer=LAYER_SUB_UNDCT )
        low_vert_undct2.translate(( x0 +sub_wire_wid/2, y0 +len_low_vert_undct -SUB_WIRS_PENET))
    SubUndct_cell.add(low_vert_undct2)

      ##horizontal low wire under inductance
    if not BROKENHEART:
        low_horz_wire = shapes.Rectangle( (0,0),(horiz_wire_len, sub_wire_wid), layer=LAYER_SUB_WIRES)
        low_horz_wire.translate(( x0 +sub_wire_wid/2, y0 +vert_wire_len -sub_wire_wid -horiz_wire_shift ))
    else:
        ## HERE is smth strange
        XXX1    = x0 +Heart_SX
        XXX_0   = x0 +horiz_wire_len -sub_wire_wid/2
        # XXX_len = XXX1 -XXX_0 +SUB_WIRS_PENET
        XXX_len = XXX1 -XXX_0
        low_horz_wire = shapes.Rectangle( (0,0), (XXX_len, sub_wire_wid), layer=LAYER_SUB_WIRES)
        low_horz_wire.translate(( XXX_0,  y0 +vert_wire_len -sub_wire_wid -horiz_wire_shift ))
    SubWires_cell.add(low_horz_wire)

      ##vertical wire from edge to JJ low
    low_vert_wire_to_JJ = shapes.Rectangle( (0,0),(sub_wire_wid, gap_JJ_to_horiz_wire), layer=LAYER_SUB_WIRES)
    low_vert_wire_to_JJ.translate(( x0 +horiz_wire_len -sub_wire_wid/2,  y0 +vert_wire_len -horiz_wire_shift ))
    SubWires_cell.add(low_vert_wire_to_JJ)

      ##vertical undercut from edge to JJ low
    if not BROKENHEART:
        low_vert_undct_to_JJ = shapes.Rectangle( (0,0),(sub_wire_undct_wid, gap_JJ_to_horiz_wire+sub_wire_wid), layer=LAYER_SUB_UNDCT)
        low_vert_undct_to_JJ.translate(( x0 +horiz_wire_len +sub_wire_wid/2,  y0 +vert_wire_len -horiz_wire_shift-sub_wire_wid ))
    else:
        low_vert_undct_to_JJ = shapes.Rectangle( (0,0),(sub_wire_undct_wid, gap_JJ_to_horiz_wire), layer=LAYER_SUB_UNDCT)
        low_vert_undct_to_JJ.translate(( x0 +horiz_wire_len +sub_wire_wid/2,  y0 +vert_wire_len -horiz_wire_shift ))
    SubUndct_cell.add(low_vert_undct_to_JJ)


    ###################################

    ## Top Part (reflected) ##
      ##vertical wire from top to inductance
    top_vert_wire = low_vert_wire.copy()
    top_vert_wire.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubWires_cell.add(top_vert_wire)

      ##undercut for vertical wire from top to inductance
    top_vert_undct1 = low_vert_undct1.copy()
    top_vert_undct1.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubUndct_cell.add(top_vert_undct1)

    top_vert_undct2 = low_vert_undct2.copy()
    top_vert_undct2.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubUndct_cell.add(top_vert_undct2)

      ##horizontal top wire over inductance
    top_horz_wire = low_horz_wire.copy()
    top_horz_wire.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubWires_cell.add(top_horz_wire)

      ##vertical wire from edge to JJ top
    top_vert_wire_to_JJ = low_vert_wire_to_JJ.copy()
    top_vert_wire_to_JJ.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubWires_cell.add(top_vert_wire_to_JJ)
      ##vertical undercut from edge to JJ low
    top_vert_undct_to_JJ = low_vert_undct_to_JJ.copy()
    top_vert_undct_to_JJ.reflect('x', ( x0, y0 +Heart_SY/2 ))
    SubUndct_cell.add(top_vert_undct_to_JJ)
    ###################################
    ## Right Part ##

      ## horizontal wire from right to JJs
    rigth_wire = shapes.Rectangle( (SUB_WIRS_PENET,0),(-len_right_wire, sub_wire_wid), layer=LAYER_SUB_WIRES)
    rigth_wire.translate(( x0 +Heart_SX,  y0 +Heart_SY/2 -sub_wire_wid/2 ))
    SubWires_cell.add(rigth_wire)
      ## vertical wire between JJ
    btwn_JJ_wire = shapes.Rectangle( (0,0),(sub_wire_wid, len_btwn_jj_wire), layer=LAYER_SUB_WIRES)
    btwn_JJ_wire.translate(( x0 +Heart_SX -len_right_wire -sub_wire_wid,  y0 +Heart_SY/2 -len_btwn_jj_wire/2 ))
    SubWires_cell.add(btwn_JJ_wire)
      ## vertical undercut between JJ
    btwn_JJ_undct = shapes.Rectangle( (0,0),(sub_wire_undct_wid, len_btwn_jj_wire), layer=LAYER_SUB_UNDCT)
    btwn_JJ_undct.translate(( x0 +Heart_SX -len_right_wire -sub_wire_wid -sub_wire_undct_wid,  y0 +Heart_SY/2 -len_btwn_jj_wire/2 ))
    SubUndct_cell.add(btwn_JJ_undct)

    ######################################################################

    if BROKENHEART:
        ## low part wire
        vertical_edge_wire = shapes.Rectangle( (0,0), (sub_wire_wid, len_low_vert_undct), layer=LAYER_SUB_WIRES )
        # vertical_edge_wire.translate(( x0 +Heart_SX +SUB_WIRS_PENET -sub_wire_wid, y0 -SUB_WIRS_PENET ))
        vertical_edge_wire.translate(( x0 +Heart_SX -sub_wire_wid, y0 -SUB_WIRS_PENET ))
        SubWires_cell.add(vertical_edge_wire)

        ## low part undrct
        vertical_edge_wire_undct = shapes.Rectangle( (0,0), (sub_wire_undct_wid, len_low_vert_undct+sub_wire_wid), layer=LAYER_SUB_UNDCT )
        vertical_edge_wire_undct.translate(( x0 +Heart_SX, y0 -SUB_WIRS_PENET ))
        SubUndct_cell.add(vertical_edge_wire_undct)

        ## top part
        vertical_edge_wire2 = vertical_edge_wire.copy()
        vertical_edge_wire2.reflect('x', ( x0, y0 +Heart_SY/2 ))
        SubUndct_cell.add(vertical_edge_wire2)

        vertical_edge_wire_undct2 = vertical_edge_wire_undct.copy()
        vertical_edge_wire_undct2.reflect('x', ( x0, y0 +Heart_SY/2 ))
        SubUndct_cell.add(vertical_edge_wire_undct2)

    ######################################################################
    ######################################################################

    ###   Assembling all heart  ###
    cell = core.Cell('Heart')
    if not HEARTLESS:
        cell.add( Inductance_cell )
        cell.add( JJ_1_cell )
        cell.add( JJ_2_cell )
    else:
        ### shortcut of JJ ###
        J1_bb = JJ_1_cell.bounding_box
        J2_bb = JJ_2_cell.bounding_box
        cell.add( shapes.Rectangle( (J1_bb[0,0], J1_bb[0,1]),  (J1_bb[1,0]+.5, J1_bb[1,1]), layer=LAYER_SUB_WIRES) )  ## !V KOSTIL
        cell.add( shapes.Rectangle( (J2_bb[0,0], J2_bb[0,1]),  (J2_bb[1,0]+.5, J2_bb[1,1]), layer=LAYER_SUB_WIRES) )  ## !V KOSTIL

        ### shortcut of L ###
        L_bb  = Inductance_cell.bounding_box
        [x0,y0,x1,y1] = [ L_bb[0,0], L_bb[0,1], L_bb[1,0], L_bb[1,1] ]
        x_center = (x0+x1)/2
        [x0,y0,x1,y1] = [x_center-WIDTH_OF_SHORTCUT_L_WIRE ,y0, x_center+WIDTH_OF_SHORTCUT_L_WIRE, y1]
        cell.add( shapes.Rectangle( (x0        ,         y0),  (x1        , y1        ), layer=LAYER_SUB_WIRES) )
        ###_______________###

    cell.add( SubWires_cell )
    cell.add( SubUndct_cell )
    return cell


################################################################################

def gen_holes_for_stripline_cell( param_stripline ):
    '''
    Making holes in proper place for Stripline design
    '''
    #########################################
    [Gap, MidLen, SideLen, MidWid, SideWid, Shift, RoundRad,
                        Strc_Gap, X_c, Y_c] = param_stripline
    holes_border = HOLES_RAD+NO_HOLES_BORDER
    RoundRad -= holes_border
    cell_0 = core.Cell('Holes_Strip')
    cell_0.add( gen_holes_net_in_RoundedCornerRectangle( (X_c -MidWid/2  +holes_border,   Y_c -MidLen/2               +holes_border), (MidWid -2*holes_border,  MidLen-2*holes_border), RoundRad ) )
    cell_0.add( gen_holes_net_in_RoundedCornerRectangle( (X_c -SideWid/2 +holes_border,   Y_c -MidLen/2 -Gap -SideLen +holes_border), (SideWid-2*holes_border, SideLen-2*holes_border), RoundRad ) )
    cell_0.add( gen_holes_net_in_RoundedCornerRectangle( (X_c -SideWid/2 +holes_border,   Y_c +MidLen/2 +Gap          +holes_border), (SideWid-2*holes_border, SideLen-2*holes_border), RoundRad ) )
    return cell_0

def gen_holes_for_circular_cell( param_circular, rad=HOLES_RAD, dist=HOLES_DIST ):
    '''
    Making holes in proper place for Circular design
    '''
    #################################
    [R3, Gap, RingWidth, Sector, PushSect, X_c, Y_c] = param_circular
    R2 = R3 -RingWidth
    R1 = R2 -Gap
    holes_border = HOLES_RAD +NO_HOLES_BORDER

    cell_0 = core.Cell('Holes_Circ')
    cell_0.add( gen_holes_net_in_circle(     (X_c,Y_c), R1-holes_border   ) )
    cell_0.add( gen_holes_net_in_ringsector( (X_c,Y_c), R2+holes_border, R3-holes_border, Sector  ) )
    return cell_0

################################################################################
def gen_test_structures_4probe(heart_cell, chip_size=(0,0), size_pads=SIZE_PADS, TEST_WR_WID=WIDWIRE_STPL_MED, TEST_WR2_WID=TEST_WR2_WID, GAPS_FROM_THE_EDGE=GAPS_FROM_THE_EDGE ):
    '''
    Making 4-probe test structure as heart with all wiring in corners of sample
    '''
    ####__Parameters__####

    (SX_pad,  SY_pad)  = size_pads
    (chip_SX, chip_SY) = chip_size
    (X_GAP_EDGE, Y_GAP_EDGE) = GAPS_FROM_THE_EDGE

    (X_0, Y_0) = (X_GAP_EDGE, Y_GAP_EDGE)
    (X_hrt_cell, Y_hrt_cell) = (X_0+4*SX_pad, Y_0+1*SY_pad) ##left_low corner of square like size_pads
    (X_hrt,   Y_hrt)   = (X_hrt_cell+SX_pad/2, Y_hrt_cell +SY_pad/2) ##center of inductance
    (SX_hrt,  SY_hrt)  = (HEART_SIZE[0], HEART_SIZE[1])
    (SX_ful_strc, SY_ful_strc) = ( 11*SX_pad, 3*SY_pad  ) ##size of full structure


    #######___DOING_WIRES____#####
    def gen_test_heart_wiring_cell():

        #################################################
        make_reflection = lambda orig_shape : orig_shape.copy().reflect('x', origin=(X_hrt,Y_hrt))
        make_reflection_v_axis = lambda orig_shape : orig_shape.copy().reflect('y', origin=(X_hrt,Y_hrt))
        #################################################

        cell_0 = core.Cell('Test_structure')

        #################################################
        ##____________________WIRES____________________##
        #################################################
        ### Tiny wires ###
        cell_sm_wr = core.Cell('Test_sm_wires')
        X = X_hrt +SX_hrt
        Y = Y_hrt +TEST_WR_WID/4
        SY = SY_hrt/2 +SUB_WIRS_PENET

        ### SMALL WIRES ################
        #### LOW part ###

        ########___to_Inductance__######
        # wr2_low = gen_rectangle( (X_hrt +SUB_WIRS_PENET, Y_hrt +SY_hrt/2),  (-SX_pad/2 -SUB_WIRS_PENET, TEST_WR_WID),  LAYER_MED_WIRES)
        wr2_low = gen_rectangle( (X_hrt +SUB_WIRS_PENET +10, Y_hrt +SY_hrt/2),  (-SX_pad/2 -SUB_WIRS_PENET -10, TEST_WR_WID),  LAYER_MED_WIRES)  ### !V !V KOSTIL KOSTIL KOSTIL
        cell_sm_wr.add(wr2_low)
        ########___to_JJs__######
        # wr3_top = gen_rectangle( (X_hrt +SX_hrt -TEST_WR_WID/2, Y_hrt +SY_hrt/2),  (SX_pad/2 -SX_hrt +TEST_WR_WID/2, TEST_WR_WID),  LAYER_MED_WIRES) ## !V
        wr3_top = gen_rectangle( (X_hrt +SX_hrt -TEST_WR_WID/2 -5, Y_hrt +SY_hrt/2),  (SX_pad/2 -SX_hrt +TEST_WR_WID/2 +5, TEST_WR_WID),  LAYER_MED_WIRES) ## !V ### !V !V KOSTIL KOSTIL KOSTIL
        # wr5_mid = gen_rectangle( (X_hrt +SX_hrt, Y_hrt -TEST_WR_WID/2),  (SX_pad/2 -SX_hrt, TEST_WR_WID),  LAYER_MED_WIRES )
        wr5_mid = gen_rectangle( (X_hrt +SX_hrt -SUB_WIRS_PENET, Y_hrt -TEST_WR_WID/2),  (SX_pad/2 -SX_hrt +SUB_WIRS_PENET, TEST_WR_WID),  LAYER_MED_WIRES )
        cell_sm_wr.add(wr3_top)
        cell_sm_wr.add(wr5_mid)
        #### TOP part ###
        ########___to_Inductance__######
        cell_sm_wr.add( make_reflection( wr2_low ))
        ########___to_JJs__######
        cell_sm_wr.add( make_reflection( wr3_top ))
        cell_0.add(cell_sm_wr)
        #################################################

        ### BIG WIRES ################
        cell_bg_wr = core.Cell('Test_bg_wires')
        ##### Right part #####
        ######################
        ######## Middle right wire ###
        wr_b_mid_right_srtght = gen_rectangle( (X_hrt_cell +SX_pad -SUB_WIRS_PENET,  Y_hrt -TEST_WR2_WID/2),   ( 5*SX_pad +2*SUB_WIRS_PENET, TEST_WR2_WID), LAYER_BIG_WIRES)
        cell_bg_wr.add( wr_b_mid_right_srtght )
        ######## Lower right wire ###
        wr_b_lower_right_h    = gen_rectangle( (X_hrt_cell +SX_pad -SUB_WIRS_PENET,  Y_hrt -SY_pad/4 -TEST_WR2_WID/2),   ( 3.5*SX_pad +SUB_WIRS_PENET +TEST_WR2_WID/2, TEST_WR2_WID), LAYER_BIG_WIRES)
        wr_b_lower_right_v    = gen_rectangle( (X_hrt_cell +4.5*SX_pad -TEST_WR2_WID/2,  Y_hrt -SY_pad/4 -TEST_WR2_WID/2),   ( TEST_WR2_WID, -SY_pad/4), LAYER_BIG_WIRES)
        cell_bg_wr.add( wr_b_lower_right_h )
        cell_bg_wr.add( wr_b_lower_right_v )
        ######## Lowest right wire ###
        wr_b_low_right_h      = gen_rectangle( (X_hrt_cell +SX_pad -SUB_WIRS_PENET, Y_hrt_cell -0.5*SY_pad +TEST_WR2_WID/2),   (SX_pad +2*SUB_WIRS_PENET, TEST_WR2_WID),  LAYER_BIG_WIRES)
        wr_b_low_right_v      = gen_rectangle( (X_hrt_cell +SX_pad -SUB_WIRS_PENET, Y_hrt_cell -0.5*SY_pad +1.5*TEST_WR2_WID),   (TEST_WR2_WID, SY_pad/2), LAYER_BIG_WIRES )
        cell_bg_wr.add( wr_b_low_right_h )
        cell_bg_wr.add( wr_b_low_right_v )
        ######## Upper right wire ###
        wr_b_upper_right_h    = make_reflection(wr_b_lower_right_h)
        wr_b_upper_right_v    = make_reflection(wr_b_lower_right_v)
        cell_bg_wr.add( wr_b_upper_right_h )
        cell_bg_wr.add( wr_b_upper_right_v )
        ######## Up right wire ###
        wr_b_up_right_h       = make_reflection(wr_b_low_right_h)
        wr_b_up_right_v       = make_reflection(wr_b_low_right_v)
        cell_bg_wr.add( wr_b_up_right_h )
        cell_bg_wr.add( wr_b_up_right_v )
        ##### Left part #####
        ######################
        ######## Up left wire ###
        wr_b_up_left_h = make_reflection_v_axis( wr_b_up_right_h )
        wr_b_up_left_v = make_reflection_v_axis( wr_b_up_right_v )
        cell_bg_wr.add( wr_b_up_left_h )
        cell_bg_wr.add( wr_b_up_left_v )
        ######## Upper left wire ###
        wr_b_upper_left_h = make_reflection_v_axis( wr_b_upper_right_h )
        wr_b_upper_left_v = make_reflection_v_axis( wr_b_upper_right_v )
        cell_bg_wr.add( wr_b_upper_left_h )
        cell_bg_wr.add( wr_b_upper_left_v )
        ######## Lowest left wire ###
        wr_b_low_left_h = make_reflection_v_axis( wr_b_low_right_h )
        wr_b_low_left_v = make_reflection_v_axis( wr_b_low_right_v )
        cell_bg_wr.add( wr_b_low_left_h )
        cell_bg_wr.add( wr_b_low_left_v )
        ######## Lower left wire ###
        wr_b_lower_left_h = make_reflection_v_axis( wr_b_lower_right_h )
        wr_b_lower_left_v = make_reflection_v_axis( wr_b_lower_right_v )
        cell_bg_wr.add( wr_b_lower_left_h )
        cell_bg_wr.add( wr_b_lower_left_v )
        ###____###
        cell_0.add(cell_bg_wr)

        #################################################
        ##____________________PADS_____________________##
        #################################################
        ###### left part ##
        cell_pads = core.Cell('Test_pads')
        cell_pads.add(gen_rectangle(  (X_hrt_cell -4*SX_pad, Y_hrt_cell +SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        cell_pads.add(gen_rectangle(  (X_hrt_cell -4*SX_pad, Y_hrt_cell -SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        ##
        cell_pads.add(gen_rectangle(  (X_hrt_cell -2*SX_pad, Y_hrt_cell +SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        cell_pads.add(gen_rectangle(  (X_hrt_cell -2*SX_pad, Y_hrt_cell -SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        ##
        ###### right part ##
        cell_pads.add(gen_rectangle(  (X_hrt_cell +2*SX_pad, Y_hrt_cell +SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        cell_pads.add(gen_rectangle(  (X_hrt_cell +2*SX_pad, Y_hrt_cell -SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        ##
        cell_pads.add(gen_rectangle(  (X_hrt_cell +4*SX_pad, Y_hrt_cell +SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        cell_pads.add(gen_rectangle(  (X_hrt_cell +4*SX_pad, Y_hrt_cell -SY_pad), (SX_pad,SY_pad), LAYER_PADS  ))
        ##
        cell_pads.add(gen_rectangle(  (X_hrt_cell +6*SX_pad, Y_hrt_cell        ), (SX_pad,SY_pad), LAYER_PADS  ))
        ############################
        cell_0.add(cell_pads)


        #################################################
        ##______________GRID_OF_MARKERS________________##
        #################################################
        # cell_grid = core.Cell('Marker_grid')
        # for i in np.arange(-4, 7, 1):
        #     cell_grid.add( gen_mark_cell(SX_pad,SY_pad, size=4, position=(X_hrt_cell +i*SX_pad, Y_hrt_cell -SY_pad) ))
        #     cell_grid.add( gen_mark_cell(SX_pad,SY_pad, size=4, position=(X_hrt_cell +i*SX_pad, Y_hrt_cell        ) ))
        #     cell_grid.add( gen_mark_cell(SX_pad,SY_pad, size=4, position=(X_hrt_cell +i*SX_pad, Y_hrt_cell +SY_pad) ))
        # cell_0.add( cell_grid )

        #################################################
        return cell_0


    #######################__EXECUTION__###################
    def move_heart_to_right_place(heart_cell, X_hrt, Y_hrt):
        magic_X_shift = -5 ### the problem - we don't know the width of Inductance here, but want to have a galvanic connection
        bb = heart_cell.bounding_box
        (low, top, right, left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
        X_shift = -left +X_hrt +magic_X_shift
        Y_shift = -low  +Y_hrt
        return copy_and_translate_cell( heart_cell, (X_shift, Y_shift) )

    cell_ts1 = core.Cell('4probe')
    heart_cell = move_heart_to_right_place( heart_cell, X_hrt, Y_hrt -SY_hrt/2 -TEST_WR_WID/2 )
    cell_ts1.add( heart_cell )

    cell_ts1.add( gen_test_heart_wiring_cell() )
    cell_ts2 = copy_and_translate_cell( cell_ts1, ( chip_SX -SX_ful_strc -2*X_GAP_EDGE, 0                                         ) )
    cell_ts3 = copy_and_translate_cell( cell_ts1, ( 0                                 , chip_SY -SY_ful_strc -2*Y_GAP_EDGE ) )
    cell_ts4 = copy_and_translate_cell( cell_ts1, ( chip_SX -SX_ful_strc -2*X_GAP_EDGE, chip_SY -SY_ful_strc -2*Y_GAP_EDGE ) )
    ####
    cell = core.Cell('4probes')
    cell.add(cell_ts1)  ## low left corner
    cell.add(cell_ts2)  ## low right corner
    # cell.add(cell_ts3)  ## top left corner
    # cell.add(cell_ts4)  ## top right corner
    return cell

def gen_transmon_cell( param_transmon, JJ_wid_mask=0.05, ind_j_hig=0.2, w_big=WIDWIRE_TRANSMON, num_of_jj=1, num_of_L_SQUIDS=0, sub_wire_len=36, PRECISE_WIR_LEN=False ):
    num_of_jj =int(abs(num_of_jj))
    num_of_L_SQUIDS = int(abs(num_of_L_SQUIDS))

    if num_of_jj * num_of_L_SQUIDS  !=0:
        print('\n!____gen_transmon_cell() error: Either num_of_jj either num_of_L must be zero!')
        return core.Cell('None_cell')

    sub_wire_wid = SUB_WIRS_WIDTH
    sub_wire_undct_wid = WIRS_UNDCT_WIDTH
    [Gap, SX_pad, SY_pad, X_c, Y_c] = param_transmon
    position_low_wire = (X_c+sub_wire_wid/2, Y_c)

    sub_wire_len += SUB_WIRS_PENET

    ####################

    ###_____Junction____###
    Len_of_structures = 0
    if num_of_L_SQUIDS ==0: ## normal situation - no inductance instead of Junction
        if num_of_jj == 1: ## normal situation (one JJ), all others used only for test-structures
            cell_junction = core.Cell('Junction')
            [[JJ_wire_shift], JJ_cell] = gen_junction( position_low_wire, JJ_wid_mask )
            cell_junction.add(JJ_cell)
            bb = JJ_cell.bounding_box
            (jj_low, jj_top, jj_right, jj_left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
            JJ1_wid = jj_right - jj_left
            JJ1_len = jj_top   - jj_low
            Len_of_structures += JJ1_len
        elif num_of_jj==0:
            cell_junction = core.Cell('no Junction')
            JJ_wire_shift = 0
            jj_low = position_low_wire[1]
            jj_top = jj_low
            JJ1_len = 0
            Len_of_structures += JJ1_len
        elif num_of_jj>1:
            cell_junction = core.Cell('Junctions_chain')
            for i in range(num_of_jj):
                [[JJ_wire_shift], JJ_cell] = gen_junction( position_low_wire, JJ_wid_mask, JJ_coupler_heigh=1.5 )
                bb = JJ_cell.bounding_box
                (jj_low, jj_top, jj_right, jj_left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
                JJ1_wid = jj_right - jj_left
                JJ1_len = jj_top   - jj_low
                JJ_cell = copy_and_translate_cell( JJ_cell, (0, i*JJ1_len) )
                if (i+1)%2 ==0:
                    JJ_cell = copy_and_reflect_cell(JJ_cell)
                cell_junction.add(JJ_cell)

            ####___###
            JJ_wire_shift = (num_of_jj%2) *JJ_wire_shift
            bb = cell_junction.bounding_box
            (jj_low, jj_top, jj_right, jj_left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
            JJ1_wid = jj_right - jj_left
            JJ1_len = jj_top   - jj_low
            Len_of_structures += JJ1_len
            ####___############################
    else: ## (if we make an Inductance instead of JJ for test structures)
        cell_junction = core.Cell('Junction')
        # !V !V !V Troubles
        JJ_cell = gen_inductance( position_low_wire, num_of_L_SQUIDS, JJ_wid_mask, ind_j_hig, squid_Loop_hig, squid_Loop_wid, shapo_gap_hig, squids_shift_X )
        ### shift ##
        bb = JJ_cell.bounding_box
        (jj_low, jj_top, jj_right, jj_left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
        JJ_mid_X = (jj_right + jj_left) /2
        JJ_cell = copy_and_translate_cell(JJ_cell, (X_c -JJ_mid_X, 0) )

        cell_junction.add(JJ_cell)
        bb = JJ_cell.bounding_box
        (jj_low, jj_top, jj_right, jj_left) = ( bb[0,1], bb[1,1], bb[1,0], bb[0,0])
        JJ1_wid = jj_right - jj_left
        JJ1_len = jj_top   - jj_low
        Len_of_structures += JJ1_len/2
        # print 'Len_of_structures', Len_of_structures
        JJ_wire_shift = 0

    ####################

    ###_____Pads____###
    if PRECISE_WIR_LEN:
        Gap    += 2*Len_of_structures
        SY_pad -= Len_of_structures
    cell_pads = core.Cell('pads_transmon')
    Pad1 = shapes.Rectangle((0,0),(SX_pad,SY_pad), layer=LAYER_PADS)
    Pad1.translate(( X_c, Y_c ))
    Pad1.translate(( -SX_pad/2, Gap/2  ))
    cell_pads.add(Pad1)
    Pad2 = Pad1.copy()
    Pad2.translate(( 0, -Gap -SY_pad ))
    cell_pads.add(Pad2)

    ####################

    ########_____Wires____########

    ######_SUB_wires_###

    ########_LOW_subwire1_##
    wire_low_small = shapes.Rectangle((0,0),(sub_wire_wid, -sub_wire_len), layer=LAYER_SUB_WIRES)
    wire_low_small.translate(( X_c, jj_low ))


    ########_TOP_subwire2_##
    wire_top_small = shapes.Rectangle((0,0),(sub_wire_wid, sub_wire_len), layer=LAYER_SUB_WIRES)
    wire_top_small.translate(( X_c +JJ_wire_shift, jj_top ))


    ########_LOW&TOP_subwire1_undct##
    wire_low_small_undct = shapes.Rectangle((0,0),(sub_wire_undct_wid, -sub_wire_len), layer=LAYER_SUB_UNDCT)
    wire_top_small_undct = shapes.Rectangle((0,0),(sub_wire_undct_wid, sub_wire_len), layer=LAYER_SUB_UNDCT)
    ##
    if   num_of_jj ==1: #normal situation of 1 JJ
        wire_low_small_undct.translate(( X_c +sub_wire_wid, jj_low ))
        wire_top_small_undct.translate(( X_c +JJ_wire_shift -sub_wire_undct_wid, jj_top ))
    elif num_of_jj ==0: ## just wire or squid-chain
        wire_low_small_undct.translate(( X_c +sub_wire_wid, jj_low ))
        wire_top_small_undct.translate(( X_c +sub_wire_wid, jj_top ))
    elif num_of_jj %2==1:
        wire_low_small_undct.translate(( X_c +sub_wire_wid, jj_low ))
        wire_top_small_undct.translate(( X_c +JJ_wire_shift -sub_wire_undct_wid, jj_top ))
    elif num_of_jj %2==0:
        wire_low_small_undct.translate(( X_c +sub_wire_wid, jj_low ))
        wire_top_small_undct.translate(( X_c +JJ_wire_shift +sub_wire_wid, jj_top ))

    cell_sm_wires = core.Cell('Small_Wires')
    cell_sm_wires.add( wire_low_small )
    cell_sm_wires.add( wire_top_small )
    cell_sm_wires.add( wire_low_small_undct )
    cell_sm_wires.add( wire_top_small_undct )

    # if num_of_L_SQUIDS != 0:
    #     pass
    #     cell_sm_wires = copy_and_translate_cell( cell_sm_wires, (0.7,0) )

    ####################

    ######_BIG_wires_###

    cell_bg_wires = core.Cell('Big_Wires')
    ########low_wire####
    low_wire = shapes.Rectangle((0,0),( w_big, (jj_low - (Y_c -Gap/2)) -sub_wire_len +2*SUB_WIRS_PENET ), layer=LAYER_MED_WIRES)
    low_wire.translate(( X_c-w_big/2,  Y_c -Gap/2 -SUB_WIRS_PENET))
    cell_bg_wires.add( low_wire )
    ########top_wire####
    top_wire = shapes.Rectangle((0,0),(    w_big,  -((Y_c +Gap/2) -jj_top -sub_wire_len +2*SUB_WIRS_PENET)    ), layer=LAYER_MED_WIRES)
    top_wire.translate(( X_c-w_big/2,  Y_c +Gap/2 +SUB_WIRS_PENET))
    cell_bg_wires.add( top_wire )


    ###__Assembling__###
    cell_0 = core.Cell('Transmon')
    cell_0.add( cell_pads )
    cell_0.add( cell_junction )
    cell_0.add( cell_sm_wires )
    cell_0.add( cell_bg_wires )
    return cell_0


################################################################################
################################################################################
####_________GENERATING_FULL_SAMPLES_FUNCTIONS_____#############################
################################################################################

def gen_v_shape_stripline_chip_cell(JJ_wid_mask=.4, ind_J_wid_mask=1.0, ind_J_hig_mask=0.2,  Gap=200, MidLen=1800, SideLen=150, MidWid=600, SideWid=600, RoundRad=60, Strc_Gap=250, name='V_shape_Stripline', HEARTLESS=False, with_holes=False):
    '''
    JJ_area and L_JJ_area here are real areas of the structure. To calculate the size of mask we use knowledge about given mask
    '''
    Shift=0 ### add assimetry
    ################################
    param_stripline = [Gap, MidLen, SideLen, MidWid, SideWid, Shift, RoundRad,
                        Strc_Gap, CHIP_X_SIZE/2, CHIP_Y_SIZE/2]
    heart_position = ( CHIP_X_SIZE/2 -MidWid/2 -Strc_Gap,
                                CHIP_Y_SIZE/2 -HEART_SIZE[1]/2 )
    ### Assembling ###
    cell_0 = core.Cell(name)
    cell_0.add( gen_mark_cell( CHIP_X_SIZE, CHIP_Y_SIZE ) )
    cell_0.add( gen_strip_pads_cell( param_stripline ) )
    if with_holes:
        cell_0.add( gen_holes_for_stripline_cell( param_stripline ) )
    cell_0.add( gen_wires_for_stripline_cell( param_stripline, heart_size=HEART_SIZE ))

    heart_cell = gen_heart_cell( heart_position, JJ_wid_mask=JJ_wid_mask, ind_J_wid_mask=ind_J_wid_mask, ind_J_hig_mask=ind_J_hig_mask, heart_size=HEART_SIZE, HEARTLESS=HEARTLESS )
    cell_0.add( heart_cell)

    broken_heart_cell = gen_heart_cell( heart_position, JJ_wid_mask=JJ_wid_mask, ind_J_wid_mask=ind_J_wid_mask, ind_J_hig_mask=ind_J_hig_mask, heart_size=HEART_SIZE, HEARTLESS=HEARTLESS, BROKENHEART=True )
    cell_0.add( gen_test_structures_4probe( broken_heart_cell, chip_size=(CHIP_X_SIZE, CHIP_Y_SIZE) ))


    if WITH_DOSETESTS:
        # ###---------------------------------------------------###
        if HEARTLESS: ## !V Kostil to add dosetesters only for Heartless stripline
            dose_tests_x_left = 3800
            dose_tests_y_low  = 800
            dose_tests_y_gap  = 200
            dose_test_sx = MidLen/2
            dose_test_sy = MidWid
            for i in range(len(LAYERS_FOR_DOSE_TEST)):
                layer = LAYERS_FOR_DOSE_TEST[i]
                x0 = dose_tests_x_left
                y0 = dose_tests_y_low + i*dose_tests_y_gap +i*dose_test_sy
                cell_0.add( gen_dosetest_rectangular(x0, y0, dose_test_sx, dose_test_sy, RoundRad, layer) )
        ###---------------------------------------------------###
    return cell_0

def gen_v_shape_circular_chip_cell( JJ_wid_mask=.4, ind_J_wid_mask=1.0, ind_J_hig_mask=0.2,  R3=710, Gap=140, RingWidth=70, Sector=0.75,                                            name='V_shape_Circle'   , HEARTLESS=False, with_holes=False):
    PushSect  = 0 ###shift in degrees (must be zero if symmetrical) ### add assimetry
    param_circular = [R3, Gap, RingWidth, Sector, PushSect, CHIP_X_SIZE/2, CHIP_Y_SIZE/2]
    R_wire_mid = R3 - RingWidth/2
    heart_position = ( CHIP_X_SIZE/2 -np.sqrt( R_wire_mid**2 -(HEART_SIZE[1]/2)**2 ),  CHIP_Y_SIZE/2 -HEART_SIZE[1]/2 )

    ### Assembling ###
    cell_0 = core.Cell(name)
    cell_0.add( gen_mark_cell( CHIP_X_SIZE, CHIP_Y_SIZE ) )
    cell_0.add( gen_round_pads_cell(param_circular ))
    if with_holes:
        cell_0.add( gen_holes_for_circular_cell(param_circular) )
    cell_0.add( gen_wires_for_circular_cell( param_circular, heart_size=HEART_SIZE, heart_position=heart_position ))

    heart_cell = gen_heart_cell( heart_position, JJ_wid_mask=JJ_wid_mask, ind_J_wid_mask=ind_J_wid_mask, ind_J_hig_mask=ind_J_hig_mask, heart_size=HEART_SIZE, HEARTLESS=HEARTLESS )
    cell_0.add( heart_cell )

    broken_heart_cell = gen_heart_cell( heart_position, JJ_wid_mask=JJ_wid_mask, ind_J_wid_mask=ind_J_wid_mask, ind_J_hig_mask=ind_J_hig_mask, heart_size=HEART_SIZE, HEARTLESS=HEARTLESS, BROKENHEART=True )
    cell_0.add( gen_test_structures_4probe( broken_heart_cell, chip_size=(CHIP_X_SIZE, CHIP_Y_SIZE)))

    return cell_0

def gen_transmon_chip_cell( JJ_wid_mask=0.1, Gap=Transmon_Gap, SX_pad=Transmon_pad_SX, SY_pad=Transmon_pad_SX, name='Transmon'):
    param_transmon = [Gap, SX_pad, SY_pad, CHIP_X_SIZE/2, CHIP_Y_SIZE/2]
    cell_0 = core.Cell(name)
    cell_0.add( gen_mark_cell( CHIP_X_SIZE, CHIP_Y_SIZE ) )
    cell_0.add( gen_transmon_cell(param_transmon, JJ_wid_mask=JJ_wid_mask)  )
    return cell_0

def gen_test_junctions_chip_cell( name='Test_junctions_chains', list_chips=None):
    print '\n\n\n____function enter!!!'
    ###__PARAMETERS___###
    ( X_shift, Y_shift ) = (400,400) ## gap from the edge of chip
    (SX_pad, SY_pad) = (150, 75)
    Distance = 200 ## distance between neighbor structures
    num_of_L_vary = 4   ##X-axis on chip
    num_of_jj_vary = 8  ##X-axis on chip
    num_of_struc = num_of_jj_vary + num_of_L_vary ##X-axis on chip
    ###_values_###
    num_of_tests = len(list_of_JJ_seizes) ## number_of_junctions_to_test ##Y-axis on chip
    Gap = ( CHIP_Y_SIZE -2*Y_shift -2*SY_pad )  / num_of_tests   -2*SY_pad

    #####################
    def gen_basics_for_chain(LABEL=False):
        pads_cell = core.Cell('basic_cell')
        pad1 = gen_rectangle( (0,0), (SX_pad,SY_pad), LAYER_PADS )
        pads_cell.add(pad1)

        basic_cell = core.Cell('basic')
        pads_cell_new = pads_cell.copy()
        for i in range(num_of_struc):
            basic_cell.add( pads_cell_new )
            pads_cell_new = copy_and_translate_cell(pads_cell_new, (Distance +SX_pad, 0) )

        ## Labels ##
        if LABEL:
            list_of_top_labels = ['J0', 'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'S2', 'S4', 'S6', 'S8' ]
            cell_labels = core.Cell('Labels_top')
            text_size = SY_pad
            for i in range(num_of_struc):
                text = list_of_top_labels[i]
                X_text = i *(Distance +SX_pad)
                cell_labels.add( shapes.Label(text, text_size, (X_text, SY_pad+text_size), layer=LAYER_LABELS_TEST) )
            basic_cell.add( cell_labels )
        ##________##

        return basic_cell

    def gen_test_junctions_chain_cell(num_of_chain, JJ_wid_mask=.4, ind_J_wid_mask=1.0, ind_J_hig_mask=0.2, label=''):
        ###__Execution__###
        chain_cell = core.Cell('JJ_wid_mask:'+str(JJ_wid_mask)+'; ind_J_wid_mask:'+str(ind_J_wid_mask))
        ####____________###
        TEST_SRT_SUBWIR_LEN = 110
        DELTA_LENGTH_TEST_LINE = -10

        delta_techn = DELTA_LENGTH_TEST_LINE/2
        subwir_len_techn = TEST_SRT_SUBWIR_LEN/2 ## two wires become one with twice length
        ###### JJ TESTS ######
        for i in range(num_of_jj_vary):
            struc_cell = core.Cell('test_struc'+str(i))
            if i == 0:
                ## if it is simple line - vary length
                sub_wire_len = subwir_len_techn + num_of_chain * delta_techn
                Gap_i    = Gap    + num_of_chain * DELTA_LENGTH_TEST_LINE
                SY_pad_i = SY_pad - (num_of_chain * DELTA_LENGTH_TEST_LINE)/2
                struc_cell.add( gen_transmon_cell([Gap_i, SX_pad, SY_pad_i, SX_pad/2, Gap_i/2 +SY_pad_i], num_of_jj=i,                         JJ_wid_mask=JJ_wid_mask,   ind_j_hig=ind_J_hig_mask,  sub_wire_len= sub_wire_len,    w_big=WIDWIRE_TEST_CHAINS, PRECISE_WIR_LEN=True ))
            else:
                struc_cell.add( gen_transmon_cell([Gap,   SX_pad, SY_pad,   SX_pad/2, Gap/2   +SY_pad  ], num_of_jj = i, JJ_wid_mask=JJ_wid_mask, sub_wire_len=subwir_len_techn, w_big=WIDWIRE_TEST_CHAINS, PRECISE_WIR_LEN=True ))
            struc_cell = copy_and_translate_cell(struc_cell, ( i*(Distance +SX_pad), 0) )
            chain_cell.add( struc_cell )

        ###### SQUIDS TESTS ######
        for i in range(num_of_L_vary):
            struc_cell = core.Cell('test_struc'+str(i))
            struc_cell.add( gen_transmon_cell([Gap, SX_pad, SY_pad, SX_pad/2, Gap/2 +SY_pad],            num_of_jj=0, num_of_L_SQUIDS=2*(i+1), JJ_wid_mask=ind_J_wid_mask, ind_j_hig=ind_J_hig_mask, sub_wire_len= subwir_len_techn, w_big=WIDWIRE_TEST_CHAINS, PRECISE_WIR_LEN=True ))
            struc_cell = copy_and_translate_cell(struc_cell, ( (num_of_jj_vary +i)*(Distance +SX_pad), 0) )
            chain_cell.add( struc_cell )

        ### Right label individual ###
        text = label
        text_size = SY_pad
        X_text = -SX_pad
        Y_text = SY_pad
        cell_label = core.Cell('Label')
        cell_label.add( shapes.Label(text, text_size, (X_text, Y_text), layer=LAYER_LABELS_TEST) )

        chain_cell.add(cell_label)

        return chain_cell

    #####################
    ###__EXECUTION___####
    cell_0 = core.Cell(name)
    cell_0.add( gen_mark_cell( CHIP_X_SIZE, CHIP_Y_SIZE ) )
    ### add boarders
    cell_0.add(copy_and_translate_cell( gen_basics_for_chain(), (X_shift, Y_shift) ))

    ### KOSTIL !V!V!V !V !V KOSTIL
    ## top labels ##
    # cell_0.add(copy_and_translate_cell( gen_basics_for_chain(LABEL=True), (X_shift, CHIP_Y_SIZE -SY_pad -Y_shift ) ))
    cell_0.add(copy_and_translate_cell( gen_basics_for_chain(LABEL=True), (X_shift, 3399 ) ))
    ## bottom labels ##
    cell_0.add(copy_and_translate_cell( gen_basics_for_chain(LABEL=True), (X_shift, 50 ) ))

    ######
    done_already = lambda x, y, done_x_list, done_y_list : (x in done_x_list) and (y in done_y_list)
    good_type = lambda type : type=='Stripline' or type=='Circular'
    ##
    counter = 0
    done_jj_areas  = []
    done_Ljj_areas = []

    # print('\n\n\n\n\n')
    # print 'LOOK HERE'
    # for chip in list_chips:
    #     print chip.num,'\t' ,chip.chiptype,'\t', chip.ind_J_hig,'\t', chip.ind_J_wid_mask,'\t'
    # print('\n\n\n\n\n')

    for chip in list_chips:
        if not good_type(chip.chiptype):
            continue
        if done_already( chip.jj_area, chip.l_jj_area, done_jj_areas, done_Ljj_areas ):
            continue
        cell_chain = gen_test_junctions_chain_cell( counter, JJ_wid_mask=chip.JJ_wid_mask, ind_J_wid_mask=chip.ind_J_wid_mask, ind_J_hig_mask=chip.ind_J_hig, label=str(chip.num) )
        cell_chain = copy_and_translate_cell(cell_chain, (X_shift, counter *(Gap +2*SY_pad) +Y_shift +SY_pad))
        cell_0.add( cell_chain )
        counter +=1
        done_jj_areas.append(chip.jj_area)
        done_Ljj_areas.append(chip.l_jj_area)
    ######
    return cell_0


################################################################################
####_________class CHIP_____#############################
################################################################################
class Chip:
    chiptype        =None
    filename        = ''
    jj_area         =None
    l_jj_area       =None
    mask_params     =None
    JJ_wid_mask     =None
    ind_J_hig       =None
    ind_J_wid_proj  =None
    ind_J_wid_mask  =None

    def gen_cell(self):
        pass
        return None

    def __init__(self, num, JJ_wid_mask, L_jj_sizes, name=None, saveandlog=True):

        self.num = num
        if name is None:
            self.name = str(self.num) +'_' +str(self.chiptype)
        else:
            self.name = name


        ### define jj_mask_sizes from areas ##
        self.mask_params    = ACTUAL_MASK.as_list()
        ## small Y-junction ##
        self.JJ_wid_mask    = JJ_wid_mask
        ## big JJ of SQUID ##

        ########################################## w6->w7 (two parameters of Ljj kostil !V)
        [L_jj_len_mask, L_jj_high_mask]  = L_jj_sizes
        self.ind_J_hig      = L_jj_high_mask
        self.ind_J_wid_mask = L_jj_len_mask

        # if type(L_jj_sizes) == type([1,2]): ## !V KOSTIL FOR TRANSMON (TEST 2PROBE)
        #     [L_jj_len_mask, L_jj_high_mask]  = L_jj_sizes
        # else:
        #     L_jj_len_mask = L_jj_sizes
        #     L_jj_high_mask = 0.2

        # self.ind_J_hig      = IND_J_HIG
        ##########################################

        ### Projected areas:
        ## small jj
        [JJ_wid1, JJ_wid2] = ACTUAL_MASK.Expected_width_of_small_Y_junction_from_width_of_mask(self.JJ_wid_mask)
        self.jj_area    = JJ_wid1*JJ_wid2
        ## big Ljj
        self.l_jj_area  = ACTUAL_MASK.Area_of_Big_junction_from_sizes_of_mask(self.ind_J_wid_mask, self.ind_J_hig)
        self.ind_J_wid_proj = self.l_jj_area /self.ind_J_hig



        self.cell = self.gen_cell()
        if saveandlog:
            self.savegds()
            self.add_to_log()

    def __str__(self):
        string = ''
        string += '\nChip_num:\t' +str(self.num)
        string += '\nType:\t' +self.chiptype
        string += '\nName:\t' + self.name
        string += '\nfilename:\t' +self.filename
        string += '\njj_area:\t'  +str(self.jj_area)
        string += '\nl_jj_area\t' +str(self.l_jj_area)
        return string

    def __repr__(self):
        return self.__str__()

    def savegds(self, folder=main_folder+'chips_gds'):
        createfolder(folder)
        self.filename = self.name
        if self.cell is not None:
            layout = core.Layout('LIBRARY')
            layout.add( self.cell )
            layout.save( folder +'\\' +self.filename +'.gds')

    def add_to_log(self, Ndigits=4):
        logline = ''
        logline += str(self.num)
        logline += '\t'+ self.chiptype

        if self.JJ_wid_mask is not None:
            pass
            logline += '\t'+str(round(self.JJ_wid_mask, 4)) ### width of JJ mask
        else:
            logline += '\t'+'None'

        if self.ind_J_wid_mask is not None:
            pass
            logline += '\t'+str(round(self.ind_J_wid_mask ,4)) ## length of Ljj mask
        else:
            logline += '\t'+'None'

        if self.ind_J_wid_mask is not None:
            pass
            logline += '\t'+str(round( self.ind_J_hig, Ndigits )) ## width of Ljj (mask = real)
        else:
            logline += '\t'+'None'

        if self.jj_area is not None:
            [JJ_w1, JJ_w2] = ACTUAL_MASK.Expected_width_of_small_Y_junction_from_width_of_mask(self.JJ_wid_mask )
            expected_small_Y_jj_area = JJ_w1 * JJ_w2
            logline += '\t'+str(round( expected_small_Y_jj_area, Ndigits ))
        else:
            logline += '\t'+'None'

        if self.l_jj_area is not None:
            expected_bigJJ_area = ACTUAL_MASK.Area_of_Big_junction_from_sizes_of_mask( self.ind_J_wid_mask, self.ind_J_hig)
            logline += '\t'+str(round( expected_bigJJ_area, Ndigits ))
        else:
            logline += '\t'+'None'

        if self.jj_area is not None:
            logline += '\t'+str(round( JJ_w1, Ndigits ))
            logline += '\t'+str(round( JJ_w2, Ndigits ))
        else:
            logline += '\t'+'None'

        if self.l_jj_area is not None:
            logline += '\t'+str(round( expected_bigJJ_area / self.ind_J_hig, Ndigits ))
        else:
            logline += '\t'+'None'



        add_line_to_log(logline)

class StriplineChip(Chip):
    chiptype = 'Stripline'
    def gen_cell(self):
        return gen_v_shape_stripline_chip_cell(JJ_wid_mask=self.JJ_wid_mask, ind_J_wid_mask=self.ind_J_wid_mask, ind_J_hig_mask=self.ind_J_hig, name=str(self.num) +'_' +self.chiptype, with_holes=True)

class CircularChip(Chip):
    chiptype = 'Circular'
    def gen_cell(self):
        return gen_v_shape_circular_chip_cell(JJ_wid_mask=self.JJ_wid_mask, ind_J_wid_mask=self.ind_J_wid_mask, ind_J_hig_mask=self.ind_J_hig, name=str(self.num) +'_' +self.chiptype, with_holes=True)

class HearlessStriplineChip(StriplineChip):
    chiptype = 'Stripline_heartless'
    def gen_cell(self):
        return gen_v_shape_stripline_chip_cell(JJ_wid_mask=self.JJ_wid_mask, ind_J_wid_mask=self.ind_J_wid_mask, ind_J_hig_mask=self.ind_J_hig, name=str(self.num) +'_' +self.chiptype, HEARTLESS=True, with_holes=True)

class HearlessCircularChip(CircularChip):
    chiptype = 'Circular_heartless'
    def gen_cell(self):
        return gen_v_shape_circular_chip_cell(JJ_wid_mask=self.JJ_wid_mask, ind_J_wid_mask=self.ind_J_wid_mask, ind_J_hig_mask=self.ind_J_hig, name=str(self.num) +'_' +self.chiptype, HEARTLESS=True, with_holes=True)

class TransmonChip(Chip):
    chiptype = 'Transmon'
    def gen_cell(self):
        self.l_jj_area = None
        self.ind_J_hig = None
        self.ind_J_wid_proj = None
        self.ind_J_wid_mask = None
        return gen_transmon_chip_cell(JJ_wid_mask=self.JJ_wid_mask, name=str(self.num)+'_' +self.chiptype )

class TestJunctionsChip(Chip):
    chiptype = 'TestJunctions'
    def __init__(self, num, all_chip_list, saveandlog=True):
        self.num = num
        self.name = str(self.num) +'_' +str(self.chiptype)
        self.jj_area = None
        self.l_jj_area = None
        self.cell = self.gen_cell( all_chip_list )
        if saveandlog:
            self.savegds()
            self.add_to_log()

    def gen_cell(self, all_chip_list):
        return gen_test_junctions_chip_cell( list_chips = all_chip_list )


################################################################################
####_________GENERATING_FULL_WaferS_FUNCTIONS_____#############################
def build_wafer():
    '''
    Make main_folder with each .gds for each chip + .gds with all wafer in one for user
    '''
    def chipposition_from_chipnumber(chip_num):
        '''
        returns chip position in the grid by its number
        number:
              |31|32|
        |25|26|27|28|29|30|
     |33|19|20|21|22|23|24|36|
  |39|34|13|14|15|16|17|18|37|40|
     |35|07|08|09|10|11|12|38|
        |01|02|03|04|05|06|

        postion: [X,Y]
              |25|35|
        |04|14|24|34|44|54|
    |-13|03|13|23|33|43|53|63|
|-22|-12|02|12|22|32|42|52|62|72|
    |-11|01|11|21|31|41|51|61|
        |00|10|20|30|40|50|
        '''
         #CHECK
        if chip_num<1 or chip_num>40:
            print 'ERROR OF SAMPLE NUMBER'
            return (None, None)

        chip_num = chip_num -1  #to start chips from "one" not "zero"

        if 0<=chip_num<=29:     ## simple case inside wafer
            X_pos = chip_num % 6
            Y_pos = chip_num // 6
        ## manual cases additional samples
            ## top couple
        elif chip_num==31-1:
                (X_pos, Y_pos) = (2, 5)
        elif chip_num==32-1:
                (X_pos, Y_pos) = (3, 5)
            ## left triple
        elif chip_num==33-1:
                (X_pos, Y_pos) = (-1, 3)
        elif chip_num==34-1:
                (X_pos, Y_pos) = (-1, 2)
        elif chip_num==35-1:
                (X_pos, Y_pos) = (-1, 1)
            ## right triple
        elif chip_num==36-1:
                (X_pos, Y_pos) = (6, 3)
        elif chip_num==37-1:
                (X_pos, Y_pos) = (6, 2)
        elif chip_num==38-1:
                (X_pos, Y_pos) = (6, 1)
            ## edge left
        elif chip_num==39-1:
                (X_pos, Y_pos) = (-2, 2)
            ## edge right
        elif chip_num==40-1:
                (X_pos, Y_pos) = (7, 2)
        else:
            (X_pos, Y_pos) = (None, None)


        return (X_pos, Y_pos)

    def coordinates_from_chipposition(X_pos, Y_pos):
            '''
            takes chipposition in the grid (can be given by chipposition_from_chipnumber() )
            and returns the real coordinates in [nm] for build a wafer and .njf file
            '''
            X =  X_pos*CHIP_X_SIZE
            Y =  Y_pos*CHIP_Y_SIZE
            return (X,Y)

    def coordinates_from_chipnumber(chip_num):
        '''
        uses functions:
        chipposition_from_chipnumber()
        coordinates_from_chipposition()
        to return coordinates of left low corner of the chip by it's number
        '''
        (X_pos, Y_pos) = chipposition_from_chipnumber(chip_num)
        print '(X_pos, Y_pos)=',(X_pos, Y_pos)
        (X,Y) = coordinates_from_chipposition(X_pos, Y_pos)
        return (X,Y)
    ############################################

    def move_to_position(cell_0, chip_num):
        '''
        function for moving chip-cell onto wafer, to put rhem all in one .gds-file
        '''
        ( X_chip, Y_chip ) = coordinates_from_chipnumber(chip_num)
        if X_chip is None or Y_chip is None:
            return core.Cell('None')
        print '_____STRIPLINE IS DONE______\n'
        return copy_and_translate_cell( cell_0, (X_chip, Y_chip) )

    Wafer = core.Cell('Wafer')

    all_chips = []

    #################
    def remove_chip_from_list(chipnum, list_of_chips):
        if chipnum in list_of_chips:
            list_of_chips.remove(chipnum)
            return list_of_chips
        else:
            print '\n\n\n ATTENTION! ERROR OF CHIPNUM FOR REMOVE \n check id_test, id_trans etc'
            return list_of_chips
    # total_number_of_chips = NUM_CHIPS_WAF_HORIZ *NUM_CHIPS_WAF_VERTC
    total_number_of_chips = 30 #!V !V magic_number
    list_all_chip_ids = range(1, total_number_of_chips+1, 1)
    list_all_chip_ids = remove_chip_from_list(id_test           , list_all_chip_ids)
    list_all_chip_ids = remove_chip_from_list(id_trans1         , list_all_chip_ids)
    list_all_chip_ids = remove_chip_from_list(id_trans2         , list_all_chip_ids)
    list_all_chip_ids = remove_chip_from_list(id_trans3         , list_all_chip_ids)
    list_all_chip_ids = remove_chip_from_list(id_NoHeart_Stripl , list_all_chip_ids)
    list_all_chip_ids = remove_chip_from_list(id_NoHeart_Circle , list_all_chip_ids)

    ##______________________________###

    ###__Generate_samples____####
    j = 0
    for i in range(len(list_of_JJ_seizes)):
        ###
        jj_size     = list_of_JJ_seizes[i]
        l_jj_length = list_of_L_JJ_seizes[i]
        l_jj_high   = list_of_L_JJ_high_sizes[i]
        l_jj_size   = [l_jj_length, l_jj_high]
        ######
        chip = StriplineChip( list_all_chip_ids[j]  , jj_size, l_jj_size)
        all_chips.append(chip)
        Wafer.add( move_to_position(chip.cell, chip.num) )
        ###
        chip  = CircularChip( list_all_chip_ids[j+1], jj_size, l_jj_size)
        all_chips.append(chip)
        Wafer.add( move_to_position(chip.cell, chip.num) )
        ######
        j+=2
    ###_______________________####


    ####___Generate_test_structures___####
    chip = TestJunctionsChip(id_test, all_chips)
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    ##__________## Transmons _____________##
    chip = TransmonChip(id_trans1, test_transmon_JJ_size_list[0], [0,0])
    # chip = TransmonChip(id_trans1,0.4,[2.1,0.3])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    chip = TransmonChip(id_trans2, test_transmon_JJ_size_list[1], [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    chip = TransmonChip(id_trans3, test_transmon_JJ_size_list[2], [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    ##__________## Heartless samples______##
    chip = HearlessStriplineChip( id_NoHeart_Stripl, 0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    chip = HearlessCircularChip( id_NoHeart_Circle,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, chip.num) )

    ##____Additional markers, arrows, etc_##
    Wafer.add( gen_wafer_markers() )
    # Wafer.add( gen_wafer_universal(x0=5396, y0=-6004) ) #it was this!
    Wafer.add( gen_wafer_universal( x0=-ORIGIN_LEFT_LOW_MARKER[0], y0=-ORIGIN_LEFT_LOW_MARKER[1] ))

    ### try new slots
    ####___Generate_test_structures___####
    chip = TestJunctionsChip(31, all_chips)
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 31) )

    chip = HearlessCircularChip( 32,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 32) )

    chip = HearlessCircularChip( 33,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 33) )

    chip = HearlessCircularChip( 34,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 34) )

    chip = HearlessCircularChip( 35,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 35) )

    chip = HearlessCircularChip( 36,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 36) )

    chip = HearlessCircularChip( 37,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 37) )

    chip = HearlessCircularChip( 38,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 38) )

    chip = HearlessCircularChip( 39,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 39) )

    chip = HearlessCircularChip( 40,  0, [0,0])
    all_chips.append(chip)
    Wafer.add( move_to_position(chip.cell, 40) )


    #################
    layout = core.Layout('LIBRARY')
    layout.add( Wafer )
    layout.save(main_folder +'Overview_' +name  +'_forUser.gds')
    return all_chips

################################################################################
####_________MAKING .NJF-FILE______________________#############################
def prepare_NJF_file(NJF_FILE_NAME, FOLDER_PATTERNS, PATTERNS, ID_FILES_DICT, DOSE_DICT, BASE_DOSE=1, FOLDER_SAVE_NJF="forLitho\\", SQR_MK_DIST=200, ZERO_POINT=(400,200), COMMENT=name ):

    COMMENT += '  ' +my_date_yyyymmdd()
    ####___TEXT_BLOCKS___#####
    boarder1 = "###------------------------------------------------------------###\n"
    boarder2 = "###____________"
    boarder3 = "\n##-----------------------------------------\n"

    ###__FINCTIONS____####

    um2nm = lambda x: int(1e3*x)

    def myinttostr(number):
        number = round(number)
        number = int(number)
        string = str(number)
        ls = len(string)
        new_string = ''
        for i in np.arange(ls-1,-1,-1):
            new_string = string[i] +new_string
            if (ls-i) %3 == 0:
                new_string = ' ' +new_string
        return new_string

    def origin(chipnum):
        position = get_chip_position_from_number(chipnum)
        X = ZERO_POINT[0] + position[0]
        Y = ZERO_POINT[1] + position[1]
        origin_string = '(' +myinttostr(um2nm(X)) +',\t' +myinttostr(um2nm(Y)) +')'
        return origin_string

    def registration_top():
        lefttop  = myinttostr( um2nm(SQR_MK_DIST))          +',\t' + myinttostr( um2nm(CHIP_Y_SIZE-SQR_MK_DIST))
        righttop = myinttostr( um2nm(CHIP_X_SIZE-SQR_MK_DIST)) +',\t' + myinttostr( um2nm(CHIP_Y_SIZE-SQR_MK_DIST))
        reg_string = '(' +lefttop +')\t(' +righttop +')'
        return reg_string

    def registration_low():
        leftlow  = myinttostr( um2nm(SQR_MK_DIST))             +',\t' + myinttostr( um2nm(SQR_MK_DIST))
        rightlow = myinttostr( um2nm(CHIP_X_SIZE-SQR_MK_DIST)) +',\t' + myinttostr( um2nm(SQR_MK_DIST))
        reg_string = '(' +leftlow +')\t(' +rightlow +')'
        return reg_string

    def stepsize():
        SX_chip = myinttostr( um2nm(CHIP_X_SIZE) )
        SY_chip = myinttostr( um2nm(CHIP_Y_SIZE) )
        chipsize_str = '(' +SX_chip +',\t' +SY_chip +')'
        return chipsize_str

    def write_title(f, title):
        endline = ''
        for i in range( len(boarder1)-len(title)-len(boarder2)-4 ) :
            endline += '_'
        f.write('\n'+boarder1)
        f.write(boarder2+title+endline+'###\n')
        f.write(boarder1+'\n')

    def write_header(f):
        f.write('.global\n')
        f.write('#registration\t(0, 34 400 000)\t(30 800 000, 34 400 000)\n') ##why?
        f.write('registration\t(0,0)\t\t\t(30 800 000, 0)\n')
        f.write('marktype		sqr8n\n')
        f.write('focus			auto\n')
        f.write('.end\n\n')
        return

    def write_blocks(f, basedose, patterns_list):
        for chipnum in np.arange(1,31,1):
            f.write('#Chip:'+str(chipnum)+' ('+patterns_list[chipnum-1]+')\n')
            f.write('.block\n')
            f.write('origin'      +'\t\t\t'+origin(chipnum)+'\n')
            f.write('registration'+'\t'    +registration_top()+'\n')
            f.write('registration'+'\t'    +registration_low()+'\n')
            f.write('marktype'    +'\t\t'     +'sqr8N'+'\n')
            f.write('focus'       +'\t\t\t'   +'map1' +'\n')
            f.write('stepsize'    +'\t\t'     +stepsize() +'\n')
            f.write('grid'        +'\t\t\t'   +'(1, 1)' +'\n')
            f.write('base_dose'   +'\t\t'     +str(basedose) +'\n')
            f.write('pattern'     +'\t\t\t'   +patterns_list[chipnum-1] +'\t\t(0,0)' +'\n')
            f.write('.end\n')
            f.write(boarder3+'\n')

    def write_patterns(f, folder, id_filenames_dict, dose_dict):
        list_keys = list( id_filenames_dict.keys() )
        list_keys.sort()    ### sort id_filenames_dict
        for key in list_keys:
            f.write('.pattern\n')
            f.write('id' +'\t\t' +key +'\n')
            f.write('filename' +'\t\t' +str(folder +id_filenames_dict[key]) +'\n')
            for d_key in dose_dict.keys():
                f.write('dose' +'\t\t' +str(d_key) +'\t\t' +str(dose_dict[d_key]) +'\n')
            f.write('.end\n')
            f.write(boarder3+'\n')
        return

    def write_ender(f):
        f.write('\n')
        f.write('.write\n')
        f.write('current \t\tauto\n')
        f.write('.end\n')

    ########################

    ########################
    #### EXECUTE ####
    directory = main_folder +FOLDER_SAVE_NJF
    createfolder( directory )

    Pads_file = open( directory +'\\' +NJF_FILE_NAME +".njf", "w")
    Pads_file.write('### '+COMMENT+'\n')

    write_title(Pads_file, 'START')
    write_header(Pads_file)
    write_title(Pads_file, 'BLOCKS')
    write_blocks(Pads_file, BASE_DOSE, PATTERNS)
    write_title(Pads_file, 'PATTERNS')
    write_patterns(Pads_file, FOLDER_PATTERNS, ID_FILES_DICT, DOSE_DICT)
    write_title(Pads_file, 'ENDER')
    write_ender(Pads_file)

    Pads_file.close()

    return True

################################################################################
############_________MAIN____________###########################################
all_chips = build_wafer()

local_litho_dir = my_date_yyyymmdd()[2:] #FORMAT YYMMDD
local_litho_dir=local_litho_dir+'_Vladimir_' +name +'_struc'
folder_with_npf = 'nicolasroch/' +local_litho_dir+'/'

'''
# [Dose] = [C/m2]
# real_dose = base_dose * coeff.
# base_dose = 1 (for pads)
# base_dose = 5 (for wires)
# here are coefficients:
'''



#### add doses for dose tests
dose_dict_heartless_dose_test = dose_dict_pads
for i in range(len(LAYERS_FOR_DOSE_TEST)):
    dose_dict_heartless_dose_test[ LAYERS_FOR_DOSE_TEST[i] ] = DOSES_FOR_DOSE_TEST[i]



id_files_dict_pads = {
    'Stripline'     : 'PadsWires_Stripline.npf',
    'Circular'      : 'PadsWires_Circular.npf',
    'Transmon'      : 'PadsWires_Transmon.npf',
    'TestJunctions' : 'PadsWires_TestJunctions.npf',
    'Stripline_heartless'   :   'PadsWires_Stripline_heartless.npf', ### with holes and dose tests
    'Circular_heartless'    :   'PadsWires_Circular_heartless.npf'   ### with holes and dose tests
}

##
patterns_for_PADS = [] #v ## patterns_for_PADS[chip.num] = chip.chiptype
patterns_for_HRTS = []
id_files_dict_structures = {}
for counter in range(len(all_chips) +1):
    for chip in all_chips:
        if chip.num == counter:
            patterns_for_PADS.append(chip.chiptype) #v
            if chip.num < 10:
                pattern = 'Hrts_0'+str(chip.num)
            else:
                pattern = 'Hrts_'+str(chip.num)
            patterns_for_HRTS.append(pattern)
            id_files_dict_structures[pattern] = chip.filename +'.npf'


########################################################

#####__________________________________####
prepare_NJF_file("Hearts_write"   , folder_with_npf, patterns_for_HRTS, id_files_dict_structures, dose_dict_structures, BASE_DOSE=5, COMMENT=name+' Hearts'   , FOLDER_SAVE_NJF = local_litho_dir)
########################################################

#####__________________________________####
prepare_NJF_file("PadsWires_write", folder_with_npf, patterns_for_PADS,    id_files_dict_pads,    dose_dict_pads      , BASE_DOSE=1, COMMENT=name+' PadsWires', FOLDER_SAVE_NJF = local_litho_dir)
########################################################




################################################################################
################################################################################
################################################################################
################################################################################
