import sys
import os

from ROOT import *
import numpy as np
from sys import exit
from numpy.random import uniform, normal
from random import choice
from pprint import pprint 
import copy

# Check out git clone https://github.com/mdj/NexDet.git and point sys.path to that place
sys.path.append("/home/philip/Documents/bachelor/kode/")
from NexDet.KalmanFilters import UncentedKalmanFilter,UKF
from NexDet.PhysicsObjects import TruthParticle, TruthTrack, Detector, Measurement, Event, RecoTrack


np.set_printoptions(edgeitems=3,infstr='Inf',linewidth=200, nanstr='NaN', precision=8, suppress=False, threshold=1000)

## NUmber of particles to simulate
n_parts = 8
doHough = False

## Information on systematics and uncertainties
pixelSize =  0.0002    # (mm)
resolution =  0.00006  #measurement resolution (mm)
tailAmplitude =   0.1  #probability of badly measured hit
tailWidth =    0.00018   #resolution of badly measured hit (mm)

# The multiple scattering times momentum per plane is estimated as follows
# Rl 50mu Si = 5.3e-4, Rl 50cm He Rl 8.8e-5
# multScattAngle=0.0175*sqrt(Rl)*(1+0.125*log(10Rl)/log(10))/sqrt(2)

multScattAngle = 0.0002  # effective theta0*E(GeV) (mult scatt) per plane
thetaxz =        0.0     # incident track angle in the xz plane

# There is an adjustable threshold with which we can get the noise occupancy
# as low is 10^-7, at a cost in the hit efficiency

noiseOccupancy = 0.00001  # noise probability in readout time window
hitEfficiency =  0.97     # probability for a track to make a hit
                          #  assume: noise=0.00001  eff=0.97
                          #                0.0001   eff=0.98
                          #                0.001    eff=0.995
                          #                0.000001 eff=0.93)

DIM_EXPERIMENT_X = 100.0
DIM_EXPERIMENT_Y = 100.0
DIM_EXPERIMENT_Z = 1000.0


DIM_MIMOSA_X = 21.5
DIM_MIMOSA_Y = 13.8
DIM_MIMOSA_Z = 0.245

NPIXELS_MIMOSA_X = 1152
NPIXELS_MIMOSA_Y = 576
NPIXELS_MIMOSA_Z = 1

MIMOSA_PITCH_X = 0.0184
MIMOSA_PITCH_Y = 0.0184
MIMOSA_PITCH_Z = 0.0

STRIP_SIZE_X = 0.0184
STRIP_SIZE_Y = 0.0184
STRIP_SIZE_Z = 0.020

POS_MIMOSA_ONE_X = 0
POS_MIMOSA_ONE_Y = 0
POS_MIMOSA_ONE_Z = 336.200

POS_MIMOSA_TWO_X = 0
POS_MIMOSA_TWO_Y = 0
POS_MIMOSA_TWO_Z = 622.200

POS_MIMOSA_THREE_X = 0
POS_MIMOSA_THREE_Y = 0
POS_MIMOSA_THREE_Z = DIM_MIMOSA_Z ## set to surface to avoid 

POS_MIMOSA_FOUR_X = 0
POS_MIMOSA_FOUR_Y = 0
POS_MIMOSA_FOUR_Z = 258.000

RANGE_MIMOSA_ONE_X = [POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_ONE_X + DIM_MIMOSA_X/2.0]
RANGE_MIMOSA_ONE_Y = [POS_MIMOSA_ONE_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_ONE_Y + DIM_MIMOSA_Y/2.0]
RANGE_MIMOSA_ONE_Z = [POS_MIMOSA_ONE_Z - DIM_MIMOSA_Z/2.0, POS_MIMOSA_ONE_Z + DIM_MIMOSA_Z/2.0]
        
RANGE_MIMOSA_TWO_X = [POS_MIMOSA_TWO_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_TWO_X + DIM_MIMOSA_X/2.0]
RANGE_MIMOSA_TWO_Y = [POS_MIMOSA_TWO_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_TWO_Y + DIM_MIMOSA_Y/2.0]
RANGE_MIMOSA_TWO_Z = [POS_MIMOSA_TWO_Z - DIM_MIMOSA_Z/2.0, POS_MIMOSA_TWO_Z + DIM_MIMOSA_Z/2.0]

RANGE_MIMOSA_THREE_X = [POS_MIMOSA_THREE_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_THREE_X + DIM_MIMOSA_X/2.0]
RANGE_MIMOSA_THREE_Y = [POS_MIMOSA_THREE_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_THREE_Y + DIM_MIMOSA_Y/2.0]
RANGE_MIMOSA_THREE_Z = [POS_MIMOSA_THREE_Z - DIM_MIMOSA_Z/2.0, POS_MIMOSA_THREE_Z + DIM_MIMOSA_Z/2.0]

RANGE_MIMOSA_FOUR_X = [POS_MIMOSA_FOUR_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_FOUR_X + DIM_MIMOSA_X/2.0]
RANGE_MIMOSA_FOUR_Y = [POS_MIMOSA_FOUR_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_FOUR_Y + DIM_MIMOSA_Y/2.0]
RANGE_MIMOSA_FOUR_Z = [POS_MIMOSA_FOUR_Z - DIM_MIMOSA_Z/2.0, POS_MIMOSA_FOUR_Z + DIM_MIMOSA_Z/2.0]

POS_MAGNET_X = 0.0
POS_MAGNET_Y = 0.0
POS_MAGNET_Z = (POS_MIMOSA_FOUR_Z + (POS_MIMOSA_ONE_Z-POS_MIMOSA_FOUR_Z)/2.0) 

DIM_MAGNET_X = 50.0
DIM_MAGNET_Y = 16.0
DIM_MAGNET_Z = 150.0

x_diamond = [0, 0, POS_MIMOSA_THREE_Z-144.0] # scattering point    

RANGE_EXPERIMENT_X = [-(DIM_EXPERIMENT_X/2), DIM_EXPERIMENT_X-x_diamond[0]]
RANGE_EXPERIMENT_Y = [-(DIM_EXPERIMENT_Y/2), DIM_EXPERIMENT_Y-x_diamond[1]]
RANGE_EXPERIMENT_Z = [x_diamond[2]-2.0, DIM_EXPERIMENT_Z-x_diamond[2]]


# c = 299792458.0 # m/s speed of light
FIELD_MAGNET_INTEGRATED = 193.89 # kG*cm (l=150 mm=0.15m)
FIELD_MAGNET_STRENGTH_Z = 1.2926 # kG
kappa = 2.99792458e-4 # GeV/c kG^-1 cm^-1



def hit_to_position(hit, det_id):
    """docstring for hit_to_position"""
    if det_id == 1:
        pos_x = RANGE_MIMOSA_ONE_X[0] + hit[0]*(DIM_MIMOSA_X/NPIXELS_MIMOSA_X)
        pos_y = RANGE_MIMOSA_ONE_Y[0] + hit[1]*(DIM_MIMOSA_Y/NPIXELS_MIMOSA_Y)
        pos_z = RANGE_MIMOSA_ONE_Z[0] + hit[2]*(DIM_MIMOSA_Z/NPIXELS_MIMOSA_Z)
        

    if det_id == 2:
        pos_x = RANGE_MIMOSA_TWO_X[0] + hit[0]*(DIM_MIMOSA_X/NPIXELS_MIMOSA_X)
        pos_y = RANGE_MIMOSA_TWO_Y[0] + hit[1]*(DIM_MIMOSA_Y/NPIXELS_MIMOSA_Y)
        pos_z = RANGE_MIMOSA_TWO_Z[0] + hit[2]*(DIM_MIMOSA_Z/NPIXELS_MIMOSA_Z)


    if det_id == 3:
        pos_x = RANGE_MIMOSA_THREE_X[0] + hit[0]*(DIM_MIMOSA_X/NPIXELS_MIMOSA_X)
        pos_y = RANGE_MIMOSA_THREE_Y[0] + hit[1]*(DIM_MIMOSA_Y/NPIXELS_MIMOSA_Y)
        pos_z = RANGE_MIMOSA_THREE_Z[0] + hit[2]*(DIM_MIMOSA_Z/NPIXELS_MIMOSA_Z)


    if det_id == 4:
        pos_x = RANGE_MIMOSA_FOUR_X[0] + hit[0]*(DIM_MIMOSA_X/NPIXELS_MIMOSA_X)
        pos_y = RANGE_MIMOSA_FOUR_Y[0] + hit[1]*(DIM_MIMOSA_Y/NPIXELS_MIMOSA_Y)
        pos_z = RANGE_MIMOSA_FOUR_Z[0] + hit[2]*(DIM_MIMOSA_Z/NPIXELS_MIMOSA_Z)
    return [pos_x, pos_y, pos_z, det_id]
    


def bounds(position):
    """docstring for bounds"""
    if (position[0] >= RANGE_EXPERIMENT_X[0] and position[0] <= RANGE_EXPERIMENT_X[1]) and (position[1] >= RANGE_EXPERIMENT_Y[0] and position[1] <= RANGE_EXPERIMENT_Y[1]) and position[2] >= RANGE_EXPERIMENT_Z[0] and position[2] <= RANGE_EXPERIMENT_Z[1]/2:
        return True 
    return False
    
def inside_mimosa_xy(position):
    """docstring for inside_mimosa_one"""
    if position:
        if (position[0] >= RANGE_MIMOSA_ONE_X[0] and position[0] <= RANGE_MIMOSA_ONE_X[1]) and (position[1] >= RANGE_MIMOSA_ONE_Y[0] and position[1] <= RANGE_MIMOSA_ONE_Y[1]):# and position[2] >= RANGE_MIMOSA_ONE_Z[0] and position[2] <= RANGE_MIMOSA_ONE_Z[1]:
            return 1
        elif (position[0] >= RANGE_MIMOSA_TWO_X[0] and position[0] <= RANGE_MIMOSA_TWO_X[1]) and(position[1] >= RANGE_MIMOSA_TWO_Y[0] and position[1] <= RANGE_MIMOSA_TWO_Y[1]):# and position[2] >= RANGE_MIMOSA_TWO_Z[0] and position[2] <= RANGE_MIMOSA_TWO_Z[1]:
            return 2
        elif (position[0] >= RANGE_MIMOSA_THREE_X[0] and position[0] <= RANGE_MIMOSA_THREE_X[1]) and (position[1] >= RANGE_MIMOSA_THREE_Y[0] and position[1] <= RANGE_MIMOSA_THREE_Y[1]):# and position[2] >= RANGE_MIMOSA_THREE_Z[0] and position[2] <= RANGE_MIMOSA_THREE_Z[1]:
            return 3
        elif (position[0] >= RANGE_MIMOSA_FOUR_X[0] and position[0] <= RANGE_MIMOSA_FOUR_X[1]) and (position[1] >= RANGE_MIMOSA_FOUR_Y[0] and position[1] <= RANGE_MIMOSA_FOUR_Y[1]):# and position[2] >= RANGE_MIMOSA_FOUR_Z[0] and position[2] <= RANGE_MIMOSA_FOUR_Z[1]:
            return 4




def inside_mimosa(position):
    """docstring for inside_mimosa_one"""
    if (position[0] >= RANGE_MIMOSA_ONE_X[0] and position[0] <= RANGE_MIMOSA_ONE_X[1]) and (position[1] >= RANGE_MIMOSA_ONE_Y[0] and position[1] <= RANGE_MIMOSA_ONE_Y[1]) and position[2] >= RANGE_MIMOSA_ONE_Z[0] and position[2] <= RANGE_MIMOSA_ONE_Z[1]:
        return 1
    elif (position[0] >= RANGE_MIMOSA_TWO_X[0] and position[0] <= RANGE_MIMOSA_TWO_X[1]) and(position[1] >= RANGE_MIMOSA_TWO_Y[0] and position[1] <= RANGE_MIMOSA_TWO_Y[1]) and position[2] >= RANGE_MIMOSA_TWO_Z[0] and position[2] <= RANGE_MIMOSA_TWO_Z[1]:
        return 2
    elif (position[0] >= RANGE_MIMOSA_THREE_X[0] and position[0] <= RANGE_MIMOSA_THREE_X[1]) and (position[1] >= RANGE_MIMOSA_THREE_Y[0] and position[1] <= RANGE_MIMOSA_THREE_Y[1]) and position[2] >= RANGE_MIMOSA_THREE_Z[0] and position[2] <= RANGE_MIMOSA_THREE_Z[1]:
        return 3
    elif (position[0] >= RANGE_MIMOSA_FOUR_X[0] and position[0] <= RANGE_MIMOSA_FOUR_X[1]) and (position[1] >= RANGE_MIMOSA_FOUR_Y[0] and position[1] <= RANGE_MIMOSA_FOUR_Y[1]) and position[2] >= RANGE_MIMOSA_FOUR_Z[0] and position[2] <= RANGE_MIMOSA_FOUR_Z[1]:
        return 4
    else:
        return -1

def pixel_hit(pos, det_id=None):
    """Estimate which pixel where hit in an ideal case"""
    if not det_id:
        det_id = inside_mimosa(pos)
    hit_x = -1 
    hit_y = -1 
    hit_z = -1
    
    
    if det_id == 1:
        hit_x = ((pos[0] - RANGE_MIMOSA_ONE_X[0])/DIM_MIMOSA_X) * NPIXELS_MIMOSA_X
        hit_y = ((pos[1] - RANGE_MIMOSA_ONE_Y[0])/DIM_MIMOSA_Y) * NPIXELS_MIMOSA_Y
        hit_z = ((pos[2] - RANGE_MIMOSA_ONE_Z[0])/DIM_MIMOSA_Z) * NPIXELS_MIMOSA_Z


    if det_id == 2:
        hit_x = ((pos[0] - RANGE_MIMOSA_TWO_X[0])/DIM_MIMOSA_X) * NPIXELS_MIMOSA_X
        hit_y = ((pos[1] - RANGE_MIMOSA_TWO_Y[0])/DIM_MIMOSA_Y) * NPIXELS_MIMOSA_Y
        hit_z = ((pos[2] - RANGE_MIMOSA_TWO_Z[0])/DIM_MIMOSA_Z) * NPIXELS_MIMOSA_Z


    if det_id == 3:
        hit_x = ((pos[0] - RANGE_MIMOSA_THREE_X[0])/DIM_MIMOSA_X) * NPIXELS_MIMOSA_X
        hit_y = ((pos[1] - RANGE_MIMOSA_THREE_Y[0])/DIM_MIMOSA_Y) * NPIXELS_MIMOSA_Y
        hit_z = ((pos[2] - RANGE_MIMOSA_THREE_Z[0])/DIM_MIMOSA_Z) * NPIXELS_MIMOSA_Z


    if det_id == 4:
        hit_x = ((pos[0] - RANGE_MIMOSA_FOUR_X[0])/DIM_MIMOSA_X) * NPIXELS_MIMOSA_X
        hit_y = ((pos[1] - RANGE_MIMOSA_FOUR_Y[0])/DIM_MIMOSA_Y) * NPIXELS_MIMOSA_Y
        hit_z = ((pos[2] - RANGE_MIMOSA_FOUR_Z[0])/DIM_MIMOSA_Z) * NPIXELS_MIMOSA_Z
    return [np.floor(hit_x), np.floor(hit_y), np.floor(hit_z), det_id] # pixelize

    


def point_plane_intersection(line0, line1, plane1, plane2, plane3):
    """From MAthworld
    http://mathworld.wolfram.com/Line-PlaneIntersection.html
    """
    
    x1 = plane1[0]
    y1 = plane1[1]
    z1 = plane1[2]

    x2 = plane2[0]
    y2 = plane2[1]
    z2 = plane2[2]

    x3 = plane3[0]
    y3 = plane3[1]
    z3 = plane3[2]

    x4 = line0[0]
    y4 = line0[1]
    z4 = line0[2]

    x5 = line1[0]
    y5 = line1[1]
    z5 = line1[2]
    
    
    tup = TMatrixD(4,4)
    tup[0][0] = 1
    tup[0][1] = 1
    tup[0][2] = 1
    tup[0][3] = 1
    
    tup[1][0] = x1
    tup[1][1] = x2
    tup[1][2] = x3
    tup[1][3] = x4
    
    tup[2][0] = y1
    tup[2][1] = y2
    tup[2][2] = y3
    tup[2][3] = y4
    
    tup[3][0] = z1
    tup[3][1] = z2
    tup[3][2] = z3
    tup[3][3] = z4


    tlow = TMatrixD(4,4)
    tlow[0][0] = 1
    tlow[0][1] = 1
    tlow[0][2] = 1
    tlow[0][3] = 0
    
    tlow[1][0] = x1
    tlow[1][1] = x2
    tlow[1][2] = x3
    tlow[1][3] = x5 - x4
    
    tlow[2][0] = y1
    tlow[2][1] = y2
    tlow[2][2] = y3
    tlow[2][3] = y5 - y4
    
    tlow[3][0] = z1
    tlow[3][1] = z2
    tlow[3][2] = z3
    tlow[3][3] = z5 - z4

    t = - tup.Determinant() / tlow.Determinant()
    

    x = x4 + (x5 - x4) * t
    y = y4 + (y5 - y4) * t
    z = z4 + (z5 - z4) * t

    if t > 0 and t <= 1:
        return [x,y,z]
    
def detect_entry_point(x0, x1, detector_id, exit_point=False):
    """
    This method calculates a timestep where a particle might 
    have entered the detector and returns the dt, to redo the 
    integration in smaller steps from x0 and how many times
    dt*dt0 should be done to catch up with the previous time
    """

    if detector_id == 1:
        px1 = POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0
        py1 = POS_MIMOSA_ONE_Y - DIM_MIMOSA_Y/2.0 
        pz1 = POS_MIMOSA_ONE_Z - DIM_MIMOSA_Z/2.0 
        plane1 = [px1, py1, pz1]

        px2 = POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py2 = POS_MIMOSA_ONE_Y - DIM_MIMOSA_Y/2.0 
        pz2 = POS_MIMOSA_ONE_Z - DIM_MIMOSA_Z/2.0 
        plane2 = [px2, py2, pz2]

        px3 = POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py3 = POS_MIMOSA_ONE_Y - DIM_MIMOSA_Y/2.0 + DIM_MIMOSA_Y
        pz3 = POS_MIMOSA_ONE_Z - DIM_MIMOSA_Z/2.0 
        plane3 = [px3, py3, pz3]
        
    if detector_id == 2:

        px1 = POS_MIMOSA_TWO_X - DIM_MIMOSA_X/2.0
        py1 = POS_MIMOSA_TWO_Y - DIM_MIMOSA_Y/2.0 
        pz1 = POS_MIMOSA_TWO_Z - DIM_MIMOSA_Z/2.0 
        plane1 = [px1, py1, pz1]

        px2 = POS_MIMOSA_TWO_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py2 = POS_MIMOSA_TWO_Y - DIM_MIMOSA_Y/2.0 
        pz2 = POS_MIMOSA_TWO_Z - DIM_MIMOSA_Z/2.0 
        plane2 = [px2, py2, pz2]

        px3 = POS_MIMOSA_TWO_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py3 = POS_MIMOSA_TWO_Y - DIM_MIMOSA_Y/2.0 + DIM_MIMOSA_Y
        pz3 = POS_MIMOSA_TWO_Z - DIM_MIMOSA_Z/2.0 
        plane3 = [px3, py3, pz3]

    if detector_id == 3:

        px1 = POS_MIMOSA_THREE_X - DIM_MIMOSA_X/2.0
        py1 = POS_MIMOSA_THREE_Y - DIM_MIMOSA_Y/2.0 
        pz1 = POS_MIMOSA_THREE_Z - DIM_MIMOSA_Z/2.0 
        plane1 = [px1, py1, pz1]

        px2 = POS_MIMOSA_THREE_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py2 = POS_MIMOSA_THREE_Y - DIM_MIMOSA_Y/2.0 
        pz2 = POS_MIMOSA_THREE_Z - DIM_MIMOSA_Z/2.0 
        plane2 = [px2, py2, pz2]

        px3 = POS_MIMOSA_THREE_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py3 = POS_MIMOSA_THREE_Y - DIM_MIMOSA_Y/2.0 + DIM_MIMOSA_Y
        pz3 = POS_MIMOSA_THREE_Z - DIM_MIMOSA_Z/2.0 
        plane3 = [px3, py3, pz3]
        
    if detector_id == 4:

        px1 = POS_MIMOSA_FOUR_X - DIM_MIMOSA_X/2.0
        py1 = POS_MIMOSA_FOUR_Y - DIM_MIMOSA_Y/2.0 
        pz1 = POS_MIMOSA_FOUR_Z - DIM_MIMOSA_Z/2.0 
        plane1 = [px1, py1, pz1]

        px2 = POS_MIMOSA_FOUR_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py2 = POS_MIMOSA_FOUR_Y - DIM_MIMOSA_Y/2.0 
        pz2 = POS_MIMOSA_FOUR_Z - DIM_MIMOSA_Z/2.0 
        plane2 = [px2, py2, pz2]

        px3 = POS_MIMOSA_FOUR_X - DIM_MIMOSA_X/2.0 + DIM_MIMOSA_X
        py3 = POS_MIMOSA_FOUR_Y - DIM_MIMOSA_Y/2.0 + DIM_MIMOSA_Y
        pz3 = POS_MIMOSA_FOUR_Z - DIM_MIMOSA_Z/2.0 
        plane3 = [px3, py3, pz3]
    
    if exit_point:
        plane1 = [px1, py1, pz1 + DIM_MIMOSA_Z]
        plane2 = [px2, py2, pz2 + DIM_MIMOSA_Z]
        plane3 = [px3, py3, pz3 + DIM_MIMOSA_Z]

    hit = point_plane_intersection(x0, x1, plane1, plane2, plane3)
    if inside_mimosa_xy(hit):
        return hit


def getBfield(pos):
    """Return B-field in Tesla for a given position"""

    
    if (pos[2] > (POS_MAGNET_Z - DIM_MAGNET_Z/4.0)  and pos[2] < (POS_MAGNET_Z + DIM_MAGNET_Z/4.0 )) and (pos[1] > (POS_MAGNET_Y - DIM_MAGNET_Y/2)  and pos[1] < (POS_MAGNET_Y + DIM_MAGNET_Y/2)) and (pos[0] > (POS_MAGNET_X - DIM_MAGNET_X/2)  and pos[0] < (POS_MAGNET_X + DIM_MAGNET_X/2)):
        a = FIELD_MAGNET_STRENGTH_Z#1.10 #kiloGauss
        # c = 2.245 # sigma kGauss
        # 
        #     z= pos[2] - POS_MAGNET_Z
        #     if z > 6.0 and z < 20.0: fieldz = a*exp(-(z-6.0)**2/(2*c*c))
        #     elif z < -6.0 and z > -20.0: fieldz = a*exp(-(z+6.0)**2/(2*c*c))
        #     elif z > -6.0 and 0 < 6.0: fieldz = a
        #     else: fieldz = 0.0
        return np.array([ 0.0, a, 0.0]).T
    else:
        return np.array([0.0, 0.0, 0.0]).T

    



def main():
    """docstring for main"""

    
    eman = TEveManager.Create()
    frm =  TEveFrameBox();
    frm.SetAABox(RANGE_EXPERIMENT_X[0], RANGE_EXPERIMENT_Y[0], RANGE_EXPERIMENT_Z[0], RANGE_EXPERIMENT_X[1], RANGE_EXPERIMENT_Y[1], RANGE_EXPERIMENT_Z[1])
    frm.SetFrameColor(kCyan);
    frm.SetBackColorRGBA(120,120,120,20);
    frm.SetDrawBack(kTRUE);

    gStyle.SetPalette(1);

    pal = TEveRGBAPalette(0, 130);
    detbox = TEveBoxSet("BoxSet");
    detbox.UseSingleColor();
    detbox.SetMainColor(kCyan-2);
    detbox.SetMainTransparency(70);

    detbox.SetPalette(pal);
    detbox.SetFrame(frm);
    detbox.Reset(TEveBoxSet.kBT_AABox, kFALSE, 64);

    # detbox.AddBox(0, 0, 50, 10, 10, 30); #x, y, z, dx, dy, dz

    # Mimosas
    detbox.AddBox(POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_ONE_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_ONE_Z - DIM_MIMOSA_Z/2.0, DIM_MIMOSA_X, DIM_MIMOSA_Y, DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
    detbox.AddBox(POS_MIMOSA_TWO_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_TWO_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_TWO_Z - DIM_MIMOSA_Z/2.0, DIM_MIMOSA_X, DIM_MIMOSA_Y, DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
    detbox.AddBox(POS_MIMOSA_THREE_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_THREE_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_THREE_Z - DIM_MIMOSA_Z/2.0, DIM_MIMOSA_X, DIM_MIMOSA_Y, DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
    detbox.AddBox(POS_MIMOSA_FOUR_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_FOUR_Y - DIM_MIMOSA_Y/2.0, POS_MIMOSA_FOUR_Z - DIM_MIMOSA_Z/2.0, DIM_MIMOSA_X, DIM_MIMOSA_Y, DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
    # Magnet box

    def detect_entry_points(x0, x1, exit_point=False):
        planes = [1,2,3,4]
        hits = []
        for plane in planes:
            hit = detect_entry_point(x0, x1, plane, exit_point)
            if hit:
                hit.append(plane) # add plane id to list
                hits.append(hit)

        return hits


    magbox = TEveBoxSet("MagSet");
    magbox.UseSingleColor();
    magbox.SetMainColor(kGreen);
    magbox.SetMainTransparency(90);
    magbox.SetPalette(pal);
    magbox.SetFrame(frm);
    magbox.Reset(TEveBoxSet.kBT_AABox, kFALSE, 64);

    magbox.AddBox(POS_MAGNET_X - DIM_MAGNET_X/2.0, POS_MAGNET_Y - DIM_MAGNET_Y/2.0, POS_MAGNET_Z - DIM_MAGNET_Z/4.0, DIM_MAGNET_X, DIM_MAGNET_Y, DIM_MAGNET_Z/2.0); #x, y, z, dx, dy, dz

    detbox.RefitPlex();
    gEve.AddElement(detbox);
    gEve.AddElement(magbox);
    gEve.Redraw3D(kTRUE);


    phist = TH1D("p", "momentum;[GeV];[Entries]", 100, 0, 10)
    pthist = TH1D("pt", "transverse momentum p_{T};[GeV];[Entries]", 100, 0, 10)

    detdist = TH1D("detector_distance", "Lengths traversed in detectors;[m];[Entries]", 100, DIM_MIMOSA_Z*0.9, DIM_MIMOSA_Z*1.1)

    hit_xz_hough_in = TH2D("hit_xz_hough_in","hit_xz_hough_in",2*int(RANGE_EXPERIMENT_Z[1]), x_diamond[2]-10.0, 700.0, NPIXELS_MIMOSA_X, POS_MIMOSA_ONE_X - DIM_MIMOSA_X/2.0, POS_MIMOSA_ONE_X + DIM_MIMOSA_X/2.0)

    hit_planes = [TH2D("hit_plane%d" % ((1+i)), "Hits in Plane %d" % ((1+i)), NPIXELS_MIMOSA_X, 0, NPIXELS_MIMOSA_X, NPIXELS_MIMOSA_Y, 0, NPIXELS_MIMOSA_Y) for i in xrange(4)]
    hit_planes_pos = [TH2D("hit_plane_pos%d" % ((1+i)), "Hits in Plane %d positrons" % ((1+i)), NPIXELS_MIMOSA_X, 0, NPIXELS_MIMOSA_X, NPIXELS_MIMOSA_Y, 0, NPIXELS_MIMOSA_Y) for i in xrange(4)]
    hit_planes_neg = [TH2D("hit_plane_neg%d" % ((1+i)), "Hits in Plane %d electrons" % ((1+i)), NPIXELS_MIMOSA_X, 0, NPIXELS_MIMOSA_X, NPIXELS_MIMOSA_Y, 0, NPIXELS_MIMOSA_Y) for i in xrange(4)]
    
    param_hits = [TH1F("x_truth", "x truth", 400, -0.1, 0.1), TH1F("y_truth", "y truth", 400, -0.1, 0.1), TH1F("tx_truth", "tx truth", 400, -0.1, 0.1), TH1F("ty_truth", "ty truth", 400, -0.1, 0.1), TH1F("qop_truth", "qop truth", 400, -1.0, 1.0)]
    
    
    force_time = TGraph()
    p_time = TGraph()
    dist_time = TGraph()
    bfield_time = TGraph()
    det_time = TGraph()


    
    

    # def multiple_scattering(hit_in, hit_out, p0, p, q):
    #     """docstring for multiple_scattering"""
    #     z = q #abs(sign(self.particle.type["id"]))
    #     X0 = 304.2 #meters, D. Green table 1.2 Radiation length   change!!!!!!! - needs to be atmospheric dependable
    #     x = sqrt((hit_out[0] - hit_in[0])**2 + (hit_out[1] - hit_in[1])**2 + (hit_out[2] - hit_in[2])**2)
    #     
    #     #PDG2006 eq:  27.12
    #     theta0 = 13.6/(beta2(self.particle.p1[0,0],self.particle.type["mass"])*momentum(self.particle.p1[0,0],self.particle.type["mass"]))*z*sqrt(x/X0)*(1+0.038*log(x/X0))
    # 
    #     z1 = random.normal(0,1) #these are independent according to PDG
    #     z2 = random.normal(0,1)
    #     rhoytheta = sqrt(3.0)/2.0 #Correlation coefficient
    #     yplane = z1 * x * theta0 * sqrt(1-rhoytheta**2) / sqrt(3) + z2 * rhoytheta *x*theta0/sqrt(3)
    #     thetaplane= z2 * theta0
    #     print yplane
    # 
    #     multiple = mat(array([[sqrt(2)*yplane],[sqrt(2)*yplane],[0.0]],'d'))
    #     print multiple
    # 
    # 
    #     ray = self.particle.x1[1:]
    #     print ray
    #     rho, theta, phi = toSperical(ray[0,0],ray[1,0],ray[2,0])
    #     plane = (Ry(theta)*Rz(phi))*mat(ray) #rotates 
    #     plane = plane + multiple #adds
    #     self.particle.x1[1:] = (Rz(-phi)*Ry(-theta))*plane #rotates back
    #     print self.particle.x1[1:]

    # def Fmag(x, p, q, m):
    #     """docstring for Bfield"""
    #     return (kappa*q) * np.cross(p/m, getBfield(x))
    # 
    # def Vmag(x,p, q, m):
    #     """docstring for Vmag"""
    #     return p/m

    # def rk4(x, p, dt):
    #     """Runge-Kutta 4'th order integration"""
    #     k1p = dt * Fmag(x,p, q, m)
    #     k1x = dt * Vmag(x,p, q, m)
    # 
    #     k2p = dt * Fmag(x + dt/2.0,p + 0.5*k1p, q, m)
    #     k2x = dt * Vmag(x + 0.5*k1x, p + dt/2.0, q, m)
    # 
    #     k3p = dt * Fmag(x + dt/2.0,p + 0.5*k2p, q, m)
    #     k3x = dt * Vmag(x + 0.5*k2x, p + dt/2.0, q, m)
    # 
    #     k4p = dt * Fmag(x + dt,p + k3p, q, m)
    #     k4x = dt * Vmag(x + k3x, p, q, m)
    # 
    #     p = p + 1.0/6.0 * (k1p + 2*k2p + 2*k3p + k4p) # RK4        
    #     x = x + 1.0/6.0 * (k1x + 2*k2x + 2*k3x + k4x) # RK4
    #     
    #     return x,p
    # 
    # def euler_step(r, dz, z):
    #     """Euler integration but in proper 5-elment format"""
    #     x = r[0]
    #     y = r[1]
    #     tx = r[2]
    #     ty = r[3]
    #     qoverp = r[4]
    # 
    #     B = getBfield([x,y,z])
    #     Bx = B[0]
    #     By = B[1]
    #     Bz = B[2]
    # 
    #     dx = dz * tx
    #     dy = dz * ty
    #     
    #     ds = dz * kappa * qoverp * sqrt(1+tx**2 + ty**2)
    #     dtx = ds * (tx*ty*Bx - (1+tx**2)*By + ty*Bz)
    #     dty = ds * ((1+ty**2)*Bx -tx*ty*By - tx*Bz)
    #     
    #     dqoverp = dz 
    # 
    #     return r + np.array([dx,dy,dtx,dty,dqoverp])


    def rk_45_step(z, r, dz=0.001, direction=1):
        """
        Calculate the Runge-Kutta 4'th degree step of the 5-element formatted version
        """

        def f(z, r):
            """Euler integration but in proper 5-elment format"""

            # print 10*"-"            
            # print "In function term"
            # print z, r
            x = r[0]
            y = r[1]
            tx = r[2]
            ty = r[3]
            qoverp = r[4]

            B = getBfield([x,y,z])
            Bx = B[0]
            By = B[1]
            Bz = B[2]

            dx = tx
            dy = ty
            # print kappa*qoverp, B, x,y,z
            dtx = kappa * qoverp * sqrt(1+tx**2 + ty**2) * (tx*ty*Bx - (1+tx**2)*By + ty*Bz)
            dty = kappa * qoverp * sqrt(1+tx**2 + ty**2) * ((1+ty**2)*Bx -tx*ty*By - tx*Bz)
            dqoverp = 0

            # print np.array([dx,dy,dtx,dty,dqoverp])
            # print 10*"-"
            return np.array([dx,dy,dtx,dty,dqoverp])
        
 
        
        # Runge-Kutta-Fehlberg method
        # Returns t,x, and the single step trunctation error, eps
        # coefficients
        c20=0.25
        c21=0.25
        c30=0.375
        c31=0.09375
        c32=0.28125
        c40=12/13.0
        c41=1932/2197.0
        c42= -7200/2197.0
        c43=7296/2197.0
        c51=439/216.0
        c52= -8.0
        c53=3680/513.0
        c54=-845/4104.0
        c60=0.5
        c61= -8/27.0
        c62=2.0
        c63= -3544/2565.0
        c64=1859/4104.0
        c65= -0.275
        a1=25/216.0
        a2=0.0
        a3=1408/2565.0
        a4=2197/4104.0
        a5= -0.2
        b1=16/135.0
        b2=0.0
        b3=6656/12825.0
        b4=28561/56430.0
        b5= -0.18
        b6=2/55.0
        # K values
        dz = direction * dz
        
        K1 = dz*f(z,r);
        K2 = dz*f(z+c20*dz, r+c21*K1);
        K3 = dz*f(z+c30*dz, r+c31*K1+c32*K2);
        K4 = dz*f(z+c40*dz, r+c41*K1+c42*K2+c43*K3);
        K5 = dz*f(z+dz,r+c51*K1+c52*K2+c53*K3+c54*K4);
        K6 = dz*f(z+c60*dz,r+c61*K1+c62*K2+c63*K3+c64*K4+c65*K5);
        r4 = r + a1*K1 + a3*K3 + a4*K4 + a5*K5; #RK4
        r = r + b1*K1 + b3*K3 + b4*K4 + b5*K5 + b6*K6; #RK5
        z = z+dz;
        eps = abs(r - r4)
        
        return z, r, np.amax(eps)

    def adpative_runge_kutta(z,zf,r):
        """Integrate from z to zf"""
        

        emax = 1.0e-3
        emin = 1.0e-8
        dzmin = 1.0e-3
        dzmax = 1.0
        iflag = 1
        delta = 0.5e-5
        dz = dzmax

        direction = 1
        if zf < z:
            direction = -1
        
        r_sol = []
        z_sol = []
        nitr = 0
        while 1: # Set iteration max here potentially
            nitr += 1
            r0 = r
            z0 = z

            if abs(dz) < dzmin:
                dz = np.sign(dz)*dzmin
                # print "dz up to dzmin"
            if abs(dz) > dzmax:
                dz = np.sign(dz)*dzmax
                # print "dz down to dzmax"

            d = abs(zf - z)
            if d <= abs(dz):
                iflag = 0 # terminate at next round
                if d <= delta*max(abs(zf), abs(z)): break
                # print "breaking up is hard to do.."

            # Call integrator
            z, r, eps = rk_45_step(z, r, dz, direction)

            if eps < emin:
                dz = 2.0*dz
                # print "Increasing dz by 2", dz
            if eps > emax:
                dz = dz/2.0
                z = z0 # Redo the previous step with the finer resolution
                r = r0
                # print "descreasing dz by 0.5 and retrying..", dz
                continue
            
            r_sol.append(r)
            z_sol.append(z)
            if iflag == 0: break
            
        return z_sol, r_sol, nitr

    reco_hits = [[],[],[],[]]
    
    
    evt = Event("Electron on Diamond scattering Simulation")
    
    evt.detector_fz = [3, 4, 1, 2]
    evt.ip0 = x_diamond

    for ipart in xrange(n_parts):
        # x = np.array([uniform(-0.002, 0.002), uniform(-0.002, 0.002), uniform(-DIM_MAGNET_Z/2.0 - 0.05, DIM_MAGNET_Z/2.0 + 0.05)]) # Position 
        # x = np.array([x_diamond[0]*normal(0,0.05),x_diamond[1]*normal(0,0.05),x_diamond[2]*normal(0,0.1)]) # Position 
        x = np.array(x_diamond) # Position 
        
        ## FIXME add normal distribution as initial position to simulate scattering point uncertainty

        p = np.array([uniform(-0.07, 0.07), uniform(-0.01, 0.01), uniform(0.1, 10.0)]) # momentum GeV
        m = 5.109989e-4 # GeV electron mass
        q = choice([-1,1]) # e
        pdgId = q*11
        t = 0.0
        spin=0.5

        color_charge = [kRed, kGreen]
        xes =  TEveLine("q/p = %d/%2.2f" % (int(q),np.linalg.norm(p)))
        xes.SetLineColor(color_charge[q > 0])
        xes.SetMarkerSize(0.2)
        xes.SetMarkerStyle(4);

        mag_points =  TEveLine("mag_points")
        mag_points.SetLineColor(kViolet)
        mag_points.SetMarkerSize(0.4)
        mag_points.SetMarkerStyle(4);


        intersects_entry = []
        for i in xrange(4):
            
            intersects1 =  TEvePointSet("intersects_entry%d" % (i + 1))
            intersects1.SetMarkerColor(kOrange)
            intersects1.SetMarkerSize(0.8)
            intersects1.SetMarkerStyle(4);
            intersects_entry.append(intersects1)

        intersects_exit = []
        for i in xrange(4):

            intersects1 =  TEvePointSet("intersects_exit%d" % (i + 1))
            intersects1.SetMarkerColor(kBlue)
            intersects1.SetMarkerSize(0.8)
            intersects1.SetMarkerStyle(4);
            intersects_exit.append(intersects1)


        detector_hits = TGraph2D()

        print("Particle %d (q/p = %d/%2.2f GeV) running" % (ipart, q, np.linalg.norm(p)))

        # Initial track parameters
        e = p / np.linalg.norm(p)
        dx = e[0]
        dy = e[1]
        dz = e[2]
        tx = dx/dz
        ty = dy/dz
        r = np.array([x[0],x[1],tx,ty,(q/np.linalg.norm(p))])
        z = x[2]
        print("Initial r\n", r)
        print("Initial z", z)
        tru_part = TruthParticle(pdgId, x, p, q, spin )
        tru_track = TruthTrack()
        tru_track.r = r
        tru_track.z = z
        tru_part.tru_track = tru_track
        evt.true_particles.append(tru_part)
        
        
        print(x_diamond)
        hit_xz_hough_in.Fill(x_diamond[2],x_diamond[0])

        
        # configure integrator
        emax = 1.0e-3
        emin = 1.0e-8
        dzmin = 1.0e-3
        dzmax = 2.0
        iflag = 1
        delta = 0.5e-5
        dz = dzmax

        zf = DIM_EXPERIMENT_Z - abs(x_diamond[2])
        
        magnet_integrated = 0.0
        magfield_length = 0.0
        while 1: # Set iteration max here potentially
            # if dz < dzmax:
                # print "dz %12.12f" %dz
        
            # print r, z
            r0 = r
            z0 = z

            for i,para in enumerate(r):
                param_hits[i].Fill(para)
                
            if abs(dz) < dzmin: # Adjust the proposed dz value if it is too small
                dz = np.sign(dz)*dzmin
            if abs(dz) > dzmax:
                dz = np.sign(dz)*dzmax # Adjust the proposed dz value if it is too large
            
            d = abs(zf - z) # stop of we have propagated to the edge
            if d <= abs(dz):
                iflag = 0 # terminate at next round
                if d <= delta*max(abs(zf), abs(z)): # terminate now if we are at the edge
                    break
                
            # Call integrator
            z, r, eps = rk_45_step(z, r, dz)
                                
            if eps < emin: # if dz is too small to be efficicent increase a bit
                dz = 2.0*dz

            if eps > emax: # if dz is too large, decrease the value by half.
                dz = dz/2.0
                z = z0 # Redo the previous step with the finer resolution
                r = r0
                continue # avoid division by zero in the hit finder below
                
        
            # Detector crossing of the detector planes
            ex = [r[0], r[1], z]
            ex0 = [r0[0], r0[1], z0]
            detector_plane_hits = detect_entry_points(ex, ex0)
            detector_plane_hits_exit = detect_entry_points(ex, ex0, True)

                
            for hit in detector_plane_hits:
                detector_hits.SetPoint(detector_hits.GetN(), hit[0], hit[1], hit[2])
                dpoint = pixel_hit(hit[0:3], hit[3])
                
                hit_planes[hit[3]-1].Fill(dpoint[0], dpoint[1])
                
                if q > 0:
                    hit_planes_pos[hit[3]-1].Fill(dpoint[0], dpoint[1])
                else:
                    hit_planes_neg[hit[3]-1].Fill(dpoint[0], dpoint[1])
                    
                intersects_entry[hit[3]-1].SetNextPoint(hit[0], hit[1], hit[2])
                reco_hits[hit[3]-1].append(hit) # Save hits for reconstruction
                
            
                # Save hit information
                # print hit[3]
                if not hit[3] in evt.detectors:
                    evt.detectors[hit[3]] = Detector()
                    evt.detectors[hit[3]].z = hit[2]
                    evt.detectors[hit[3]].name = "Plane %d" % hit[3]

                pixelized_hit = hit_to_position(dpoint, hit[3])
                m = Measurement(evt.detectors[hit[3]])
                m.position = pixelized_hit #hit convert to global position again
                m.local_coordinate = dpoint
                m.true_particle = tru_part
                m.true_track = tru_track
                m.true_parameters_at_surface = copy.deepcopy(r),z
                tru_track.measurements.append(m) # link back to measurement from truth track
                evt.detectors[hit[3]].hits.append(m)
                evt.hits.append(pixelized_hit)
                print(pixelized_hit)
                hit_xz_hough_in.Fill(pixelized_hit[2],pixelized_hit[0])

            for hitf in detector_plane_hits_exit:
                intersects_exit[hitf[3]-1].SetNextPoint(hitf[0], hitf[1], hitf[2])

            
            # Check some stats on the B-field            

            bsum = sum(getBfield([r[0],r[1],z]))
            if abs(bsum) > 0.001:
                # print bsum
                mag_points.SetMarkerSize(bsum)
                mag_points.SetNextPoint(r[0],r[1],z)
                # print "field", bsum, dz
                magnet_integrated += dz * bsum
                magfield_length += dz
                
                
            xes.SetNextPoint(r[0], r[1], z)
            
            if iflag == 0:  # if we are told to break by the above stuff, do so..
                break
        
        print("Integrated field: %f kG*cm ~   (nominal: %f kG*cm) (l=%f cm)" % (magnet_integrated, FIELD_MAGNET_INTEGRATED, magfield_length))
        # Draw Event
        gEve.AddElement(xes);
        for intersecter in intersects_entry:
            gEve.AddElement(intersecter)

        for intersecter in intersects_exit:
            gEve.AddElement(intersecter)

        gEve.AddElement(mag_points)
        gEve.Redraw3D()
    # raw_input("DO")
    ###########################################################################
    if doHough:
        
        xz = np.zeros([hit_xz_hough_in.GetNbinsX(), hit_xz_hough_in.GetNbinsY()])
        print(xz.shape)
        # z_bounds = 

        for x in xrange(1,hit_xz_hough_in.GetNbinsX()):
            for y in xrange(1,hit_xz_hough_in.GetNbinsY()):
                xz[x,y] = hit_xz_hough_in.GetBinContent(x,y)
        rho, theta, H = hough_transform(xz)
        print(rho)
        print(theta)
        print(H)

        hit_xz_hough_out = TH2D("meg","meg", len(theta), theta[0], theta[-1], len(rho), rho[0], rho[-1])
        for x in xrange(1,hit_xz_hough_out.GetNbinsX()):
            for y in xrange(1,hit_xz_hough_out.GetNbinsY()):
                hit_xz_hough_out.SetBinContent(x,y,H[y,x])

        ctho = TCanvas("hough")
        ctho.Divide(2,1)
        ctho.cd(1)
        hit_xz_hough_in.Draw("COLZ")
        ctho.cd(2)
        hit_xz_hough_out.Draw("COLZ")
        raw_input("m")


        # sys.exit(0)

    
    def predictor_func(r0, aux):
        """docstring for predictor"""
        z0 = aux["z"]
        zf = aux["zf"]
        z, r, nitr = adpative_runge_kutta(z0,zf,r0)
        # print "predicted after |z|=%f, " % (z), r
        # print "Number of hits in precition range (%f-%f m): %d" % (z0,zf,len(r)) 
        if "pred" in aux:
            for i in xrange(len(r)):
                # if r[i][4] > 0.0:
                #     aux["pred"].SetMarkerColor(kGreen)
                # else:
                #     aux["pred"].SetMarkerColor(kRed)
                    
                aux["pred"].SetNextPoint(r[i][0], r[i][1], z[i])

        return r[-1]
        
    def measurement_func(r, aux):
        """
        Transform prediction state into measurement (resolution, efficiency etc..)
        """
        return np.array([r[0], r[1]])
        
    


    c = TCanvas("c")
    detdist.Draw()

    cc = TCanvas("hits", "hit planes", 1000, 1000)
    cc.Divide(2,2)
    for i,h in enumerate(hit_planes_neg):
        cc.cd(i+1)
        h.Draw("")
        h.SetMarkerColor(kRed)
        h.SetFillColor(kRed)
        h.SetMarkerSize(1.0)
        h.SetMarkerStyle(21)
        
    for i,h in enumerate(hit_planes_pos):
        cc.cd(i+1)
        h.Draw("same")
        h.SetMarkerColor(kGreen)
        h.SetFillColor(kGreen)
        h.SetMarkerSize(1.0)
        h.SetMarkerStyle(21)


    c2 = TCanvas("para_hits_g")
    c2.Divide(1,5)
    for i,h in enumerate(param_hits):
        c2.cd(1+i)
        h.Draw()
        



    def add_to_param(trk, param, z):
        """docstring for add_to_param"""
        for k,p in enumerate(param):
            trk.param_graph[k].SetPoint(trk.param_graph[k].GetN(), z, p[0])
            

    top_view = TCanvas("topview")
    top_view.Divide(1,2)
    top_view.cd(1)
    h = TH2F("Top_view", "Top view (x-z plane)", 200, RANGE_EXPERIMENT_Z[0]-10, RANGE_EXPERIMENT_Z[1], 300, -DIM_MIMOSA_X/2, DIM_MIMOSA_X/2)
    h.DrawCopy()

    # Loop over the last two detector layers
    det0 = evt.detectors[evt.detector_fz[0]]
    det1 = evt.detectors[evt.detector_fz[1]]
    det2 = evt.detectors[evt.detector_fz[2]]
    det3 = evt.detectors[evt.detector_fz[3]]
    trk_top_lines = []
    f = TF1("pol1", "pol1")

    ip_hit = TGraph()
    ip_hit.SetMarkerSize(1)
    ip_hit.SetMarkerStyle(5)
    ip_hit.SetMarkerColor(kBlack)
    ip_hit.SetPoint(0, evt.ip0[2], evt.ip0[0])
    ip_hit.Draw("SAME P")
    d0_hits = TGraph()
    d0_hits.SetMarkerSize(1)
    d0_hits.SetMarkerStyle(4)
    d0_hits.SetMarkerColor(kRed)
    for hl0 in det0.hits:
        d0_hits.SetPoint(d0_hits.GetN(), hl0.position[2], hl0.position[0])
    d0_hits.Draw("SAME P")

    d1_hits = TGraph()
    d1_hits.SetMarkerSize(1)
    d1_hits.SetMarkerStyle(4)
    d1_hits.SetMarkerColor(kGreen)
    for hl0 in det1.hits:
        d1_hits.SetPoint(d1_hits.GetN(), hl0.position[2], hl0.position[0])
    d1_hits.Draw("SAME P")

    d2_hits = TGraph()
    d2_hits.SetMarkerSize(1)
    d2_hits.SetMarkerStyle(4)
    d2_hits.SetMarkerColor(kBlue)
    for hl0 in det2.hits:
        d2_hits.SetPoint(d2_hits.GetN(), hl0.position[2], hl0.position[0])
    d2_hits.Draw("SAME P")

    d3_hits = TGraph()
    d3_hits.SetMarkerSize(1)
    d3_hits.SetMarkerStyle(4)
    d3_hits.SetMarkerColor(kOrange)
    for hl0 in det3.hits:
        d3_hits.SetPoint(d3_hits.GetN(), hl0.position[2], hl0.position[0])
    d3_hits.Draw("SAME P")


        
    class TrackSeedHit(object):
            """docstring for TrackSeedHit"""
            def __init__(self, hit):
                super(TrackSeedHit, self).__init__()
                self.hit = hit
                self.fit_with_parent_prob = None
                self.fit_with_parent_chi2 = None
                self.fit_with_parent_ndof = None
                self.parent = None
                self.children = []

            def self_parent(self):
                out = [self]
                if self.parent:
                    out += self.parent.self_parent()
                return out
            def get_all_hit_chains_from_self(self):
                '''
                Returns an exploded list of all combinations of hits starting from this node
                '''
                chain = []



            def __repr__(self):
                return "Hit at z=%2.4f" % self.hit.position[2]


    track_seeds_outside_in = []
    track_seeds_inside_out= []
    # Iterative fitter
    for hl3 in det3.hits:
        x3 = hl3.position[0]
        z3 = hl3.position[2]

        trk_seed = TrackSeedHit(hl3)
        track_seeds_outside_in.append(trk_seed)

        for hl2 in det2.hits:
            x2 = hl2.position[0]
            z2 = hl2.position[2]

            tmp = TGraph()
            tmp.SetPoint(tmp.GetN(), z3, x3)
            tmp.SetPoint(tmp.GetN(), z2, x2)
            tmp.SetPoint(tmp.GetN(), evt.ip0[2], evt.ip0[0])

            f = TF1("pol1", "pol1")
            tmp.Fit(f,"Q")
            hl2chi = f.GetChisquare()
            if f.GetProb() > 0.95 and f.GetChisquare() < 0.0001:
                trk_seed_l2 = TrackSeedHit(hl2)
                trk_seed_l2.fit_with_parent_prob = f.GetProb()
                trk_seed_l2.fit_with_parent_chi2 = f.GetChisquare()
                trk_seed_l2.fit_with_parent_ndof = f.GetNDF()

                trk_seed_l2.parent = trk_seed
                trk_seed.children.append(trk_seed_l2)

                for hl1 in det1.hits:
                    x1 = hl1.position[0]
                    z1 = hl1.position[2]
                    tmp = TGraph()
                    tmp.SetPoint(tmp.GetN(), z3, x3)
                    tmp.SetPoint(tmp.GetN(), z2, x2)
                    tmp.SetPoint(tmp.GetN(), z1, x1)
                    tmp.SetPoint(tmp.GetN(), evt.ip0[2], evt.ip0[0])

                    tmp.Fit(f,"Q")
                    hl1chi = f.GetChisquare()
                    
                    # print "prob, ",f.GetProb(), "chi2", f.GetChisquare()
                    if f.GetProb() > 0.99 and f.GetChisquare() < 0.0001:
                        # print "selected"
                        trk_seed_l1 = TrackSeedHit(hl1)
                        trk_seed_l1.parent = trk_seed_l2 # current hit
                        trk_seed_l1.fit_with_parent_prob = f.GetProb()
                        trk_seed_l1.fit_with_parent_chi2 = f.GetChisquare()
                        trk_seed_l1.fit_with_parent_ndof = f.GetNDF()

                        trk_seed_l2.children.append(trk_seed_l1) # previous hit


                        for hl0 in det0.hits:
                            x0 = hl0.position[0]
                            z0 = hl0.position[2]
                            tmp = TGraph()
                            tmp.SetPoint(tmp.GetN(), z3, x3)
                            tmp.SetPoint(tmp.GetN(), z2, x2)
                            tmp.SetPoint(tmp.GetN(), z1, x1)
                            tmp.SetPoint(tmp.GetN(), z0, x0)
                            tmp.SetPoint(tmp.GetN(), evt.ip0[2], evt.ip0[0])
                            tmp.Fit(f,"Q")
                            hl0chi = f.GetChisquare()
                            print("prob, ",f.GetProb(), "chi2", f.GetChisquare())
                            
                            if f.GetProb() > 0.999 and f.GetChisquare() < 0.0004:
                                print("selecteded")
                                trk_seed_l0 = TrackSeedHit(hl0)
                                trk_seed_l0.parent = trk_seed_l1 # current hit
                                trk_seed_l0.fit_with_parent_prob = f.GetProb()
                                trk_seed_l0.fit_with_parent_chi2 = f.GetChisquare()
                                trk_seed_l0.fit_with_parent_ndof = f.GetNDF()

                                trk_seed_l1.children.append(trk_seed_l0) # previous hit
                                track_seeds_inside_out.append(trk_seed_l0)
                                # tmp.Draw("L SAME")
                                # trk["hits"].append(hl0)
                                # raw_input("Perfect track")


    top_view.cd(2)
    hx = TH2F("Side_view", "Side view (y-z plane)", 100, RANGE_EXPERIMENT_Z[0]-10, RANGE_EXPERIMENT_Z[1], 300, -DIM_MIMOSA_Y/2,DIM_MIMOSA_Y/2)
    hx.DrawCopy()
    vip_hit = TGraph()
    vip_hit.SetMarkerSize(1)
    vip_hit.SetMarkerStyle(5)
    vip_hit.SetMarkerColor(kBlack)
    vip_hit.SetPoint(0, evt.ip0[2], evt.ip0[1])
    vip_hit.Draw("SAME P")
    vd0_hits = TGraph()
    vd0_hits.SetMarkerSize(1)
    vd0_hits.SetMarkerStyle(4)
    vd0_hits.SetMarkerColor(kRed)
    for hl0 in det0.hits:
        vd0_hits.SetPoint(vd0_hits.GetN(), hl0.position[2], hl0.position[1])
    vd0_hits.Draw("SAME P")

    vd1_hits = TGraph()
    vd1_hits.SetMarkerSize(1)
    vd1_hits.SetMarkerStyle(4)
    vd1_hits.SetMarkerColor(kGreen)
    for hl0 in det1.hits:
        vd1_hits.SetPoint(vd1_hits.GetN(), hl0.position[2], hl0.position[1])
    vd1_hits.Draw("SAME P")

    vd2_hits = TGraph()
    vd2_hits.SetMarkerSize(1)
    vd2_hits.SetMarkerStyle(4)
    vd2_hits.SetMarkerColor(kBlue)
    for hl0 in det2.hits:
        vd2_hits.SetPoint(vd2_hits.GetN(), hl0.position[2], hl0.position[1])
    vd2_hits.Draw("SAME P")

    vd3_hits = TGraph()
    vd3_hits.SetMarkerSize(1)
    vd3_hits.SetMarkerStyle(4)
    vd3_hits.SetMarkerColor(kOrange)
    for hl0 in det3.hits:
        vd3_hits.SetPoint(vd3_hits.GetN(), hl0.position[2], hl0.position[1])
    vd3_hits.Draw("SAME P")


    track_draw_stack  = []
    trklet_extrapolated_stack = []
    trklet_r0extrapolated_stack =[]
    
    dRdp = TProfile("dRdp", "ratio as a function of momentum", 100, 0, 10)
    pBias = TH1D("pBias", "diff", 200, -10, 10)
    
    
    def propagate_gfx(r, z0, z1):
        """docstring for propagate_gfx"""
        trklet =  TEveLine("Extrapolated r0  trk %d" % iii)
        # trklet.SetDefaultSmooth(True)
        trklet.SetMarkerColor(kOrange)
        trklet.SetLineColor(kBlue)
        # trklet.SetMarkerSize(0.9)
        # trklet.SetMarkerStyle(4);
        predictor_func(r, {"z" : z0, "zf" : z1, "pred" : trklet})
        gEve.AddElement(trklet)
        return trklet
        
    for iii,track in enumerate(track_seeds_inside_out):
        proto_tracks =  track.self_parent()
        print(proto_tracks)
        gtx = TGraph()
        gty = TGraph()

        if len(proto_tracks) == 4:
            print(20*"+")

            # Before magnet
            p0 = np.array(proto_tracks[0].hit.position[0:3])
            p1 = np.array(proto_tracks[1].hit.position[0:3])
            e = ((p1-p0) / np.linalg.norm(p1-p0))
            dx = e[0]
            dy = e[1]
            dz = e[2]
            txi = dx/dz
            tyi = dy/dz

            # After magnet
            p2 = np.array(proto_tracks[2].hit.position[0:3])
            p3 = np.array(proto_tracks[3].hit.position[0:3])
            e = ((p3-p2) / np.linalg.norm(p3-p2))
            dx = e[0]
            dy = e[1]
            dz = e[2]
            txf = dx/dz
            tyf = dy/dz

            int_field_tesla_m = 0.019389 # 1.5 is a fudgefactor, don;t know why... se Bos. E,
            poverq = (int_field_tesla_m*1.5)/(tyf/ np.sqrt(1.0 + txf*txf + tyf*tyf) - tyi/np.sqrt(1.0 + txi*txi + tyi*tyi))
            # source
            # 1.  Bos, E. Reconstruction of charged particles in the LHCb experiment
            # . nikhef.nl at <http://www.nikhef.nl/pub/services/biblio/theses_pdf/thesis_E_Bos.pdf>
            # page 70.
            
            print("pg: %10.10f" % poverq)
            ptrue = np.linalg.norm(proto_tracks[0].hit.true_particle.p)
            qtrue = proto_tracks[0].hit.true_particle.charge
            print("tru p %10.10f, q = %d" % (ptrue, qtrue)) 
            dRdp.Fill(ptrue, abs(poverq)/ptrue)
            pBias.Fill(ptrue - abs(poverq))

            track.match = False
            if (ptrue - abs(poverq)) < 0.1:
                track.match = True
                
            e = ((p0-evt.ip0) / np.linalg.norm(p0-evt.ip0))
            dx = e[0]
            dy = e[1]
            dz = e[2]
            tx = dx/dz
            ty = dy/dz
            track.r = np.array([[evt.ip0[0],evt.ip0[1],tx,ty,1.0/poverq]]).T


            # track.r = np.array([[p0[0],p0[1],txi,tyi,np.sign(pq)/abs(pq)]]).T
            track.r0 = copy.deepcopy(track.r)
            track.z = evt.ip0[2] #inital z point

            track.P = np.diag([1,1,1,1,1])
            track.Q = np.diag([0.1,0.1,0.1,0.1,0.1])
            track.R = np.diag([0.0000001, 0.0000001])# Create a noise matrix for each measurement

                        
            n_det_hits = len(proto_tracks)
            n_param = 5
            MM = np.zeros((n_param, n_det_hits+1))
            PP = np.zeros((n_param, n_param, n_det_hits+1))
            ZZ = np.zeros((1, n_det_hits+1))
            MM[:,0] = track.r.T
            PP[:,:,0] = track.P
            ZZ[:,0] = track.z
            f_param_smo = []
            
            ukf = UKF() # Unscented Kalamn filter instance
            for det_idx in xrange(n_det_hits):
                    ## Kalman stuff below
                    print("Propagating from %f to %f" %(track.z, proto_tracks[det_idx].hit.position[2]))
                    
                    # Prediction
                    f_param = {"z" : track.z, "zf" : proto_tracks[det_idx].hit.position[2]}
                    m,P,D = ukf.predict1(copy.deepcopy(track.r),copy.deepcopy(track.P),predictor_func,track.Q, f_param)   
                    
                    # Convert to measurement frame
                    m,P,K0,MU0,S0,LH0 = ukf.update1(copy.deepcopy(m),copy.deepcopy(P),np.array([[proto_tracks[det_idx].hit.position[0], proto_tracks[det_idx].hit.position[1]]]).T, measurement_func,track.R)
                    print("True r", proto_tracks[det_idx].hit.true_parameters_at_surface)
                    print("Estimated", m.T)
                    # Use this for track rejection at some point
                    print("\nLH0 %f %f\n" % (LH0[0], LH0[1]), "match", track.match)
                    
                    # Save for smoothing
                    MM[:,det_idx+1] = m.T
                    PP[:,:,det_idx+1] = P
                    ZZ[:,det_idx+1] = f_param["zf"]
                    
                    # Update state
                    f_param_smo.append({"zf" : f_param["z"], "z" : f_param["zf"]}) # for smoothing reversed
                    track.z = f_param["zf"]
                    track.r = m
                    track.P = P
            
            # Smooth track
            f_param_smo.reverse()
            # track.Q = np.diag([0,0,0,0,0])
            for par in f_param_smo:
                print(par)
            SM,SP,D = ukf.smooth1(copy.deepcopy(MM),copy.deepcopy(PP),predictor_func,track.Q, f_param_smo,same_p=False)

            # Smoothed vs fitted but unsmoothed
            print("MM\n",MM)
            print("SM\n",SM)
            print("ZZ\n", ZZ)
            track.r = np.array([SM[:,0]]).T
                        
            print("Initial track guess")
            print(track.r0)
            print("Final track fit")
            print(track.r[:,0])
            
            for ii in SP:
                print(ii)
            
            color_charge = [kRed, kGreen]
                        
            # trklet_extrapolated_stack.append(propagate_gfx(SM[:,0], x_diamond[2], RANGE_EXPERIMENT_Z[1]))
            mm = propagate_gfx(MM[:,0], x_diamond[2], RANGE_EXPERIMENT_Z[1])
            mm.SetLineStyle(2)
            mm.SetLineColor(color_charge[MM[4,0] > 0]-2)
            trklet_extrapolated_stack.append(mm)
            

            
            gEve.Redraw3D()
            print("final momentum estimate after smoothing %f GeV" % (abs(1.0/track.r[4,0])))
        for seedhit in proto_tracks:
            gty.SetPoint(gty.GetN(), seedhit.hit.position[2], seedhit.hit.position[1])            
            gtx.SetPoint(gtx.GetN(), seedhit.hit.position[2], seedhit.hit.position[0])            

        top_view.cd(2)

        gty.SetLineColor(kRed)
        gty.Draw("L SAME")
        # gty.Fit("pol3","Q")
        
        top_view.cd(1)
        gtx.Draw("L SAME")
        # gtx.Fit("pol1","Q")
        track_draw_stack.append(gtx)
        track_draw_stack.append(gty)



    cccc = TCanvas("meh")
    dRdp.Draw()
    ccccc =  TCanvas("dfdsf")
    
    pBias.Draw()
    raw_input("Exit [enter]: ")
    
def hough_transform(img_bin, theta_res=1, rho_res=1):
    nR,nC = img_bin.shape
    theta = np.linspace(-90.0, 0.0, np.ceil(90.0/theta_res) + 1.0)
    theta = np.concatenate((theta, -theta[len(theta)-2::-1]))

    D = np.sqrt((nR - 1)**2 + (nC - 1)**2)
    q = np.ceil(D/rho_res)
    nrho = 2*q + 1
    rho = np.linspace(-q*rho_res, q*rho_res, nrho)
    H = np.zeros((len(rho), len(theta)))
    for rowIdx in range(nR):
        for colIdx in range(nC):
            if img_bin[rowIdx, colIdx]:
                for thIdx in range(len(theta)):
                    rhoVal = colIdx*np.cos(theta[thIdx]*np.pi/180.0) + \
                        rowIdx*np.sin(theta[thIdx]*np.pi/180)
                    rhoIdx = np.nonzero(np.abs(rho-rhoVal) == np.min(np.abs(rho-rhoVal)))[0]
                    H[rhoIdx[0], thIdx] += 1
    return rho, theta, H
        
if __name__ == '__main__':
    main()


#################### 
#
# TODO
#
#  - Rescale to MM
#  - Make kalman filter work
#  - port to c++
#  - restructure everything
#  - look at alignment 