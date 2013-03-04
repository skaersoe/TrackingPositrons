from ROOT import *
import numpy as np
import sys, os


# Check out git clone https://github.com/mdj/NexDet.git and point sys.path to that place
sys.path.append("/Users/mdj/Dropbox/Research/NexDet")
from NexDet.KalmanFilters import UncentedKalmanFilter,UKF
from NexDet.PhysicsObjects import TruthParticle, TruthTrack, Detector, Measurement, Event, RecoTrack


class NA63Detector(object):
    """docstring for NA63Detector"""
    def __init__(self):
        super(NA63Detector, self).__init__()
    
        self.pixelSize =  0.0002    # (mm)
        self.resolution =  0.00006  #measurement resolution (mm)
        self.tailAmplitude =   0.1  #probability of badly measured hit
        self.tailWidth =    0.00018   #resolution of badly measured hit (mm)

        # The multiple scattering times momentum per plane is estimated as follows
        # Rl 50mu Si = 5.3e-4, Rl 50cm He Rl 8.8e-5
        # multScattAngle=0.0175*sqrt(Rl)*(1+0.125*log(10Rl)/log(10))/sqrt(2)

        self.multScattAngle = 0.0002  # effective theta0*E(GeV) (mult scatt) per plane
        self.thetaxz =        0.0     # incident track angle in the xz plane

        # There is an adjustable threshold with which we can get the noise occupancy
        # as low is 10^-7, at a cost in the hit efficiency

        self.noiseOccupancy = 0.00001  # noise probability in readout time window
        self.hitEfficiency =  0.97     # probability for a track to make a hit
                                  #  assume: noise=0.00001  eff=0.97
                                  #                0.0001   eff=0.98
                                  #                0.001    eff=0.995
                                  #                0.000001 eff=0.93)



        self.DIM_MIMOSA_X = 21.5
        self.DIM_MIMOSA_Y = 13.8
        self.DIM_MIMOSA_Z = 0.245

        self.NPIXELS_MIMOSA_X = 1152
        self.NPIXELS_MIMOSA_Y = 576
        self.NPIXELS_MIMOSA_Z = 1

        self.MIMOSA_PITCH_X = 0.0184
        self.MIMOSA_PITCH_Y = 0.0184
        self.MIMOSA_PITCH_Z = 0.0

        self.STRIP_SIZE_X = 0.0184
        self.STRIP_SIZE_Y = 0.0184
        self.STRIP_SIZE_Z = 0.020

        self.POS_MIMOSA_ONE_X = 1.516
        self.POS_MIMOSA_ONE_Y = -1.193
        self.POS_MIMOSA_ONE_Z = 336.200
        self.TILT_MIMOSA_ONE_X = 0.181
        self.TILT_MIMOSA_ONE_Y = 0.0
        self.TILT_MIMOSA_ONE_Z = 0.0
        
        
        self.POS_MIMOSA_TWO_X = 2.640
        self.POS_MIMOSA_TWO_Y = -1.805
        self.POS_MIMOSA_TWO_Z = 622.200
        self.TILT_MIMOSA_TWO_X = 0.533
        self.TILT_MIMOSA_TWO_Y = 0.0
        self.TILT_MIMOSA_TWO_Z = 0.0



        self.POS_MIMOSA_THREE_X = 0
        self.POS_MIMOSA_THREE_Y = 0
        self.POS_MIMOSA_THREE_Z = 0.0 ## set to surface to avoid 
        self.TILT_MIMOSA_THREE_X = 1.5
        self.TILT_MIMOSA_THREE_Y = 0.0
        self.TILT_MIMOSA_THREE_Z = 0.0

        self.POS_MIMOSA_FOUR_X = 1.133
        self.POS_MIMOSA_FOUR_Y = -0.757
        self.POS_MIMOSA_FOUR_Z = 258.000
        self.TILT_MIMOSA_FOUR_X = 1.734
        self.TILT_MIMOSA_FOUR_Y = 0.0
        self.TILT_MIMOSA_FOUR_Z = 0.0

        self.RANGE_MIMOSA_ONE_X = [self.POS_MIMOSA_ONE_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_ONE_X + self.DIM_MIMOSA_X/2.0]
        self.RANGE_MIMOSA_ONE_Y = [self.POS_MIMOSA_ONE_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_ONE_Y + self.DIM_MIMOSA_Y/2.0]
        self.RANGE_MIMOSA_ONE_Z = [self.POS_MIMOSA_ONE_Z - self.DIM_MIMOSA_Z/2.0, self.POS_MIMOSA_ONE_Z + self.DIM_MIMOSA_Z/2.0]

        self.RANGE_MIMOSA_TWO_X = [self.POS_MIMOSA_TWO_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_TWO_X + self.DIM_MIMOSA_X/2.0]
        self.RANGE_MIMOSA_TWO_Y = [self.POS_MIMOSA_TWO_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_TWO_Y + self.DIM_MIMOSA_Y/2.0]
        self.RANGE_MIMOSA_TWO_Z = [self.POS_MIMOSA_TWO_Z - self.DIM_MIMOSA_Z/2.0, self.POS_MIMOSA_TWO_Z + self.DIM_MIMOSA_Z/2.0]

        self.RANGE_MIMOSA_THREE_X = [self.POS_MIMOSA_THREE_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_THREE_X + self.DIM_MIMOSA_X/2.0]
        self.RANGE_MIMOSA_THREE_Y = [self.POS_MIMOSA_THREE_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_THREE_Y + self.DIM_MIMOSA_Y/2.0]
        self.RANGE_MIMOSA_THREE_Z = [self.POS_MIMOSA_THREE_Z - self.DIM_MIMOSA_Z/2.0, self.POS_MIMOSA_THREE_Z + self.DIM_MIMOSA_Z/2.0]

        self.RANGE_MIMOSA_FOUR_X = [self.POS_MIMOSA_FOUR_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_FOUR_X + self.DIM_MIMOSA_X/2.0]
        self.RANGE_MIMOSA_FOUR_Y = [self.POS_MIMOSA_FOUR_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_FOUR_Y + self.DIM_MIMOSA_Y/2.0]
        self.RANGE_MIMOSA_FOUR_Z = [self.POS_MIMOSA_FOUR_Z - self.DIM_MIMOSA_Z/2.0, self.POS_MIMOSA_FOUR_Z + self.DIM_MIMOSA_Z/2.0]

        self.POS_MAGNET_X = 0.0
        self.POS_MAGNET_Y = 0.0
        self.POS_MAGNET_Z = (self.POS_MIMOSA_FOUR_Z + (self.POS_MIMOSA_ONE_Z-self.POS_MIMOSA_FOUR_Z)/2.0) 

        
        self.DIM_MAGNET_X = 50.0
        self.DIM_MAGNET_Y = 16.0
        self.DIM_MAGNET_Z = 150.0

        
        self.x_diamond = [0, 0, self.POS_MIMOSA_THREE_Z-144.0] # scattering point    

        self.DIM_EXPERIMENT_X = 100.0
        self.DIM_EXPERIMENT_Y = 100.0
        self.DIM_EXPERIMENT_Z = 1000.0

        self.RANGE_EXPERIMENT_X = [-(self.DIM_EXPERIMENT_X/2), self.DIM_EXPERIMENT_X-self.x_diamond[0]]
        self.RANGE_EXPERIMENT_Y = [-(self.DIM_EXPERIMENT_Y/2), self.DIM_EXPERIMENT_Y-self.x_diamond[1]]
        self.RANGE_EXPERIMENT_Z = [self.x_diamond[2]-2.0, self.DIM_EXPERIMENT_Z-self.x_diamond[2]]


        # c = 299792458.0 # m/s speed of light
        self.FIELD_MAGNET_INTEGRATED = 193.89 # kG*cm (l=150 mm=0.15m)
        self.FIELD_MAGNET_STRENGTH_Z = 1.2926 # kG
        self.kappa = 2.99792458e-4 # GeV/c kG^-1 cm^-1



    def get_center_of_detector(self, det_id):
        """docstring for get_center_of_detector"""
        
        if det_id == 1:
            pos_x = self.POS_MIMOSA_ONE_X 
            pos_y = self.POS_MIMOSA_ONE_Y 
            pos_z = self.POS_MIMOSA_ONE_Z 


        if det_id == 2:
            pos_x = self.POS_MIMOSA_TWO_X 
            pos_y = self.POS_MIMOSA_TWO_Y 
            pos_z = self.POS_MIMOSA_TWO_Z 


        if det_id == 3:
            pos_x = self.POS_MIMOSA_THREE_X 
            pos_y = self.POS_MIMOSA_THREE_Y 
            pos_z = self.POS_MIMOSA_THREE_Z 


        if det_id == 4:
            pos_x = self.POS_MIMOSA_FOUR_X 
            pos_y = self.POS_MIMOSA_FOUR_Y 
            pos_z = self.POS_MIMOSA_FOUR_Z 
            
        return [pos_x, pos_y, pos_z]

    def hit_to_position(self, hit, det_id):
        """docstring for hit_to_position"""
        if det_id == 1:
            pos_x = self.RANGE_MIMOSA_ONE_X[0] + hit[0]*(self.DIM_MIMOSA_X/self.NPIXELS_MIMOSA_X)
            pos_y = self.RANGE_MIMOSA_ONE_Y[0] + hit[1]*(self.DIM_MIMOSA_Y/self.NPIXELS_MIMOSA_Y)
            pos_z = self.RANGE_MIMOSA_ONE_Z[0] + hit[2]*(self.DIM_MIMOSA_Z/self.NPIXELS_MIMOSA_Z)


        if det_id == 2:
            pos_x = self.RANGE_MIMOSA_TWO_X[0] + hit[0]*(self.DIM_MIMOSA_X/self.NPIXELS_MIMOSA_X)
            pos_y = self.RANGE_MIMOSA_TWO_Y[0] + hit[1]*(self.DIM_MIMOSA_Y/self.NPIXELS_MIMOSA_Y)
            pos_z = self.RANGE_MIMOSA_TWO_Z[0] + hit[2]*(self.DIM_MIMOSA_Z/self.NPIXELS_MIMOSA_Z)


        if det_id == 3:
            pos_x = self.RANGE_MIMOSA_THREE_X[0] + hit[0]*(self.DIM_MIMOSA_X/self.NPIXELS_MIMOSA_X)
            pos_y = self.RANGE_MIMOSA_THREE_Y[0] + hit[1]*(self.DIM_MIMOSA_Y/self.NPIXELS_MIMOSA_Y)
            pos_z = self.RANGE_MIMOSA_THREE_Z[0] + hit[2]*(self.DIM_MIMOSA_Z/self.NPIXELS_MIMOSA_Z)


        if det_id == 4:
            pos_x = self.RANGE_MIMOSA_FOUR_X[0] + hit[0]*(self.DIM_MIMOSA_X/self.NPIXELS_MIMOSA_X)
            pos_y = self.RANGE_MIMOSA_FOUR_Y[0] + hit[1]*(self.DIM_MIMOSA_Y/self.NPIXELS_MIMOSA_Y)
            pos_z = self.RANGE_MIMOSA_FOUR_Z[0] + hit[2]*(self.DIM_MIMOSA_Z/self.NPIXELS_MIMOSA_Z)
        return [pos_x, pos_y, pos_z, det_id]



    def bounds(self, position):
        """docstring for bounds"""
        if (position[0] >= self.RANGE_EXPERIMENT_X[0] and position[0] <= self.RANGE_EXPERIMENT_X[1]) and (position[1] >= self.RANGE_EXPERIMENT_Y[0] and position[1] <= self.RANGE_EXPERIMENT_Y[1]) and position[2] >= self.RANGE_EXPERIMENT_Z[0] and position[2] <= self.RANGE_EXPERIMENT_Z[1]/2:
            return True 
        return False

    def inside_mimosa_xy(self, position):
        """docstring for inside_mimosa_one"""
        if position:
            if (position[0] >= self.RANGE_MIMOSA_ONE_X[0] and position[0] <= self.RANGE_MIMOSA_ONE_X[1]) and (position[1] >= self.RANGE_MIMOSA_ONE_Y[0] and position[1] <= self.RANGE_MIMOSA_ONE_Y[1]):# and position[2] >= self.RANGE_MIMOSA_ONE_Z[0] and position[2] <= self.RANGE_MIMOSA_ONE_Z[1]:
                return 1
            elif (position[0] >= self.RANGE_MIMOSA_TWO_X[0] and position[0] <= self.RANGE_MIMOSA_TWO_X[1]) and(position[1] >= self.RANGE_MIMOSA_TWO_Y[0] and position[1] <= self.RANGE_MIMOSA_TWO_Y[1]):# and position[2] >= self.RANGE_MIMOSA_TWO_Z[0] and position[2] <= self.RANGE_MIMOSA_TWO_Z[1]:
                return 2
            elif (position[0] >= self.RANGE_MIMOSA_THREE_X[0] and position[0] <= self.RANGE_MIMOSA_THREE_X[1]) and (position[1] >= self.RANGE_MIMOSA_THREE_Y[0] and position[1] <= self.RANGE_MIMOSA_THREE_Y[1]):# and position[2] >= self.RANGE_MIMOSA_THREE_Z[0] and position[2] <= self.RANGE_MIMOSA_THREE_Z[1]:
                return 3
            elif (position[0] >= self.RANGE_MIMOSA_FOUR_X[0] and position[0] <= self.RANGE_MIMOSA_FOUR_X[1]) and (position[1] >= self.RANGE_MIMOSA_FOUR_Y[0] and position[1] <= self.RANGE_MIMOSA_FOUR_Y[1]):# and position[2] >= self.RANGE_MIMOSA_FOUR_Z[0] and position[2] <= self.RANGE_MIMOSA_FOUR_Z[1]:
                return 4




    def inside_mimosa(self, position):
        """docstring for inside_mimosa_one"""
        if (position[0] >= self.RANGE_MIMOSA_ONE_X[0] and position[0] <= self.RANGE_MIMOSA_ONE_X[1]) and (position[1] >= self.RANGE_MIMOSA_ONE_Y[0] and position[1] <= self.RANGE_MIMOSA_ONE_Y[1]) and position[2] >= self.RANGE_MIMOSA_ONE_Z[0] and position[2] <= self.RANGE_MIMOSA_ONE_Z[1]:
            return 1
        elif (position[0] >= self.RANGE_MIMOSA_TWO_X[0] and position[0] <= self.RANGE_MIMOSA_TWO_X[1]) and(position[1] >= self.RANGE_MIMOSA_TWO_Y[0] and position[1] <= self.RANGE_MIMOSA_TWO_Y[1]) and position[2] >= self.RANGE_MIMOSA_TWO_Z[0] and position[2] <= self.RANGE_MIMOSA_TWO_Z[1]:
            return 2
        elif (position[0] >= self.RANGE_MIMOSA_THREE_X[0] and position[0] <= self.RANGE_MIMOSA_THREE_X[1]) and (position[1] >= self.RANGE_MIMOSA_THREE_Y[0] and position[1] <= self.RANGE_MIMOSA_THREE_Y[1]) and position[2] >= self.RANGE_MIMOSA_THREE_Z[0] and position[2] <= self.RANGE_MIMOSA_THREE_Z[1]:
            return 3
        elif (position[0] >= self.RANGE_MIMOSA_FOUR_X[0] and position[0] <= self.RANGE_MIMOSA_FOUR_X[1]) and (position[1] >= self.RANGE_MIMOSA_FOUR_Y[0] and position[1] <= self.RANGE_MIMOSA_FOUR_Y[1]) and position[2] >= self.RANGE_MIMOSA_FOUR_Z[0] and position[2] <= self.RANGE_MIMOSA_FOUR_Z[1]:
            return 4
        else:
            return -1

    def pixel_hit(self, pos, det_id=None):
        """Estimate which pixel where hit in an ideal case"""
        if not det_id:
            det_id = inside_mimosa(pos)
        hit_x = -1 
        hit_y = -1 
        hit_z = -1


        if det_id == 1:
            hit_x = ((pos[0] - self.RANGE_MIMOSA_ONE_X[0])/self.DIM_MIMOSA_X) * self.NPIXELS_MIMOSA_X
            hit_y = ((pos[1] - self.RANGE_MIMOSA_ONE_Y[0])/self.DIM_MIMOSA_Y) * self.NPIXELS_MIMOSA_Y
            hit_z = ((pos[2] - self.RANGE_MIMOSA_ONE_Z[0])/self.DIM_MIMOSA_Z) * self.NPIXELS_MIMOSA_Z


        if det_id == 2:
            hit_x = ((pos[0] - self.RANGE_MIMOSA_TWO_X[0])/self.DIM_MIMOSA_X) * self.NPIXELS_MIMOSA_X
            hit_y = ((pos[1] - self.RANGE_MIMOSA_TWO_Y[0])/self.DIM_MIMOSA_Y) * self.NPIXELS_MIMOSA_Y
            hit_z = ((pos[2] - self.RANGE_MIMOSA_TWO_Z[0])/self.DIM_MIMOSA_Z) * self.NPIXELS_MIMOSA_Z


        if det_id == 3:
            hit_x = ((pos[0] - self.RANGE_MIMOSA_THREE_X[0])/self.DIM_MIMOSA_X) * self.NPIXELS_MIMOSA_X
            hit_y = ((pos[1] - self.RANGE_MIMOSA_THREE_Y[0])/self.DIM_MIMOSA_Y) * self.NPIXELS_MIMOSA_Y
            hit_z = ((pos[2] - self.RANGE_MIMOSA_THREE_Z[0])/self.DIM_MIMOSA_Z) * self.NPIXELS_MIMOSA_Z


        if det_id == 4:
            hit_x = ((pos[0] - self.RANGE_MIMOSA_FOUR_X[0])/self.DIM_MIMOSA_X) * self.NPIXELS_MIMOSA_X
            hit_y = ((pos[1] - self.RANGE_MIMOSA_FOUR_Y[0])/self.DIM_MIMOSA_Y) * self.NPIXELS_MIMOSA_Y
            hit_z = ((pos[2] - self.RANGE_MIMOSA_FOUR_Z[0])/self.DIM_MIMOSA_Z) * self.NPIXELS_MIMOSA_Z
        return [np.floor(hit_x), np.floor(hit_y), np.floor(hit_z), det_id] # pixelize




    def point_plane_intersection(self, line0, line1, plane1, plane2, plane3):
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

    def detect_entry_point(self, x0, x1, detector_id, exit_point=False):
        """
        This method calculates a timestep where a particle might 
        have entered the detector and returns the dt, to redo the 
        integration in smaller steps from x0 and how many times
        dt*dt0 should be done to catch up with the previous time
        """

        if detector_id == 1:
            px1 = self.POS_MIMOSA_ONE_X - self.DIM_MIMOSA_X/2.0
            py1 = self.POS_MIMOSA_ONE_Y - self.DIM_MIMOSA_Y/2.0 
            pz1 = self.POS_MIMOSA_ONE_Z - self.DIM_MIMOSA_Z/2.0 
            plane1 = [px1, py1, pz1]

            px2 = self.POS_MIMOSA_ONE_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py2 = self.POS_MIMOSA_ONE_Y - self.DIM_MIMOSA_Y/2.0 
            pz2 = self.POS_MIMOSA_ONE_Z - self.DIM_MIMOSA_Z/2.0 
            plane2 = [px2, py2, pz2]

            px3 = self.POS_MIMOSA_ONE_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py3 = self.POS_MIMOSA_ONE_Y - self.DIM_MIMOSA_Y/2.0 + self.DIM_MIMOSA_Y
            pz3 = self.POS_MIMOSA_ONE_Z - self.DIM_MIMOSA_Z/2.0 
            plane3 = [px3, py3, pz3]

        if detector_id == 2:

            px1 = self.POS_MIMOSA_TWO_X - self.DIM_MIMOSA_X/2.0
            py1 = self.POS_MIMOSA_TWO_Y - self.DIM_MIMOSA_Y/2.0 
            pz1 = self.POS_MIMOSA_TWO_Z - self.DIM_MIMOSA_Z/2.0 
            plane1 = [px1, py1, pz1]

            px2 = self.POS_MIMOSA_TWO_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py2 = self.POS_MIMOSA_TWO_Y - self.DIM_MIMOSA_Y/2.0 
            pz2 = self.POS_MIMOSA_TWO_Z - self.DIM_MIMOSA_Z/2.0 
            plane2 = [px2, py2, pz2]

            px3 = self.POS_MIMOSA_TWO_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py3 = self.POS_MIMOSA_TWO_Y - self.DIM_MIMOSA_Y/2.0 + self.DIM_MIMOSA_Y
            pz3 = self.POS_MIMOSA_TWO_Z - self.DIM_MIMOSA_Z/2.0 
            plane3 = [px3, py3, pz3]

        if detector_id == 3:

            px1 = self.POS_MIMOSA_THREE_X - self.DIM_MIMOSA_X/2.0
            py1 = self.POS_MIMOSA_THREE_Y - self.DIM_MIMOSA_Y/2.0 
            pz1 = self.POS_MIMOSA_THREE_Z - self.DIM_MIMOSA_Z/2.0 
            plane1 = [px1, py1, pz1]

            px2 = self.POS_MIMOSA_THREE_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py2 = self.POS_MIMOSA_THREE_Y - self.DIM_MIMOSA_Y/2.0 
            pz2 = self.POS_MIMOSA_THREE_Z - self.DIM_MIMOSA_Z/2.0 
            plane2 = [px2, py2, pz2]

            px3 = self.POS_MIMOSA_THREE_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py3 = self.POS_MIMOSA_THREE_Y - self.DIM_MIMOSA_Y/2.0 + self.DIM_MIMOSA_Y
            pz3 = self.POS_MIMOSA_THREE_Z - self.DIM_MIMOSA_Z/2.0 
            plane3 = [px3, py3, pz3]

        if detector_id == 4:

            px1 = self.POS_MIMOSA_FOUR_X - self.DIM_MIMOSA_X/2.0
            py1 = self.POS_MIMOSA_FOUR_Y - self.DIM_MIMOSA_Y/2.0 
            pz1 = self.POS_MIMOSA_FOUR_Z - self.DIM_MIMOSA_Z/2.0 
            plane1 = [px1, py1, pz1]

            px2 = self.POS_MIMOSA_FOUR_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py2 = self.POS_MIMOSA_FOUR_Y - self.DIM_MIMOSA_Y/2.0 
            pz2 = self.POS_MIMOSA_FOUR_Z - self.DIM_MIMOSA_Z/2.0 
            plane2 = [px2, py2, pz2]

            px3 = self.POS_MIMOSA_FOUR_X - self.DIM_MIMOSA_X/2.0 + self.DIM_MIMOSA_X
            py3 = self.POS_MIMOSA_FOUR_Y - self.DIM_MIMOSA_Y/2.0 + self.DIM_MIMOSA_Y
            pz3 = self.POS_MIMOSA_FOUR_Z - self.DIM_MIMOSA_Z/2.0 
            plane3 = [px3, py3, pz3]

        if exit_point:
            plane1 = [px1, py1, pz1 + self.DIM_MIMOSA_Z]
            plane2 = [px2, py2, pz2 + self.DIM_MIMOSA_Z]
            plane3 = [px3, py3, pz3 + self.DIM_MIMOSA_Z]

        hit = point_plane_intersection(x0, x1, plane1, plane2, plane3)
        if inside_mimosa_xy(hit):
            return hit


    def getBfield(self, pos):
        """Return B-field in Tesla for a given position"""


        if (pos[2] > (self.POS_MAGNET_Z - self.DIM_MAGNET_Z/4.0)  and pos[2] < (self.POS_MAGNET_Z + self.DIM_MAGNET_Z/4.0 )) and (pos[1] > (self.POS_MAGNET_Y - self.DIM_MAGNET_Y/2)  and pos[1] < (self.POS_MAGNET_Y + self.DIM_MAGNET_Y/2)) and (pos[0] > (self.POS_MAGNET_X - self.DIM_MAGNET_X/2)  and pos[0] < (self.POS_MAGNET_X + self.DIM_MAGNET_X/2)):
            a = self.FIELD_MAGNET_STRENGTH_Z#1.10 #kiloGauss
            # c = 2.245 # sigma kGauss
            # 
            #     z= pos[2] - self.POS_MAGNET_Z
            #     if z > 6.0 and z < 20.0: fieldz = a*exp(-(z-6.0)**2/(2*c*c))
            #     elif z < -6.0 and z > -20.0: fieldz = a*exp(-(z+6.0)**2/(2*c*c))
            #     elif z > -6.0 and 0 < 6.0: fieldz = a
            #     else: fieldz = 0.0
            return np.array([0, a, 0.0]).T
        else:
            return np.array([0.0, 0.0, 0.0]).T


    def detect_entry_points(self, x0, x1, exit_point=False):
        planes = [1,2,3,4]
        hits = []
        for plane in planes:
            hit = self.detect_entry_point(x0, x1, plane, exit_point)
            if hit:
                hit.append(plane) # add plane id to list
                hits.append(hit)

        return hits

    def createEve(self):
        """docstring for createEve"""
        self.eman = TEveManager.Create()
        frm =  TEveFrameBox();
        frm.SetAABox(self.RANGE_EXPERIMENT_X[0], self.RANGE_EXPERIMENT_Y[0], self.RANGE_EXPERIMENT_Z[0], self.RANGE_EXPERIMENT_X[1], self.RANGE_EXPERIMENT_Y[1], self.RANGE_EXPERIMENT_Z[1])
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
        detbox.AddBox(self.POS_MIMOSA_ONE_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_ONE_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_ONE_Z - self.DIM_MIMOSA_Z/2.0, self.DIM_MIMOSA_X, self.DIM_MIMOSA_Y, self.DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
        detbox.AddBox(self.POS_MIMOSA_TWO_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_TWO_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_TWO_Z - self.DIM_MIMOSA_Z/2.0, self.DIM_MIMOSA_X, self.DIM_MIMOSA_Y, self.DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
        detbox.AddBox(self.POS_MIMOSA_THREE_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_THREE_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_THREE_Z - self.DIM_MIMOSA_Z/2.0, self.DIM_MIMOSA_X, self.DIM_MIMOSA_Y, self.DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
        detbox.AddBox(self.POS_MIMOSA_FOUR_X - self.DIM_MIMOSA_X/2.0, self.POS_MIMOSA_FOUR_Y - self.DIM_MIMOSA_Y/2.0, self.POS_MIMOSA_FOUR_Z - self.DIM_MIMOSA_Z/2.0, self.DIM_MIMOSA_X, self.DIM_MIMOSA_Y, self.DIM_MIMOSA_Z); #x, y, z, dx, dy, dz
        # Magnet box


        magbox = TEveBoxSet("MagSet");
        magbox.UseSingleColor();
        magbox.SetMainColor(kGreen);
        magbox.SetMainTransparency(90);
        magbox.SetPalette(pal);
        magbox.SetFrame(frm);
        magbox.Reset(TEveBoxSet.kBT_AABox, kFALSE, 64);

        magbox.AddBox(self.POS_MAGNET_X - self.DIM_MAGNET_X/2.0, self.POS_MAGNET_Y - self.DIM_MAGNET_Y/2.0, self.POS_MAGNET_Z - self.DIM_MAGNET_Z/4.0, self.DIM_MAGNET_X, self.DIM_MAGNET_Y, self.DIM_MAGNET_Z/2.0); #x, y, z, dx, dy, dz

        detbox.RefitPlex();
        gEve.AddElement(detbox);
        gEve.AddElement(magbox);
        gEve.Redraw3D(kTRUE);
        
    

# f = TFile("/Users/mdj/Dropbox/Research/AUTestBeamGPU/taf.dev/datDSF/run7_05.root", "read")
f = TFile("/Users/mdj/Dropbox/Research/AUTestBeamGPU/taf.dev/datDSF/run7_07.root", "read")
# f = TFile("/Users/mdj/Dropbox/Research/AUTestBeamGPU/taf.dev/datDSF/run18_02.root", "read")

t = f.Get("T") # TTeee


class DetectorHit(object):
    """docstring for DetectorHit"""
    def __init__(self, det, v, u):
        super(DetectorHit, self).__init__()
        self.det = det
        self.v = v # vertical
        self.u = u # horisontal
                
def getPlanes(evt=1):
    """docstring for getPlanes"""
    sel = "Hu:Hv:Hpk"
    #           Hpk Hu          Hv
    blacklist = [[2, 10239.599, 2235.5998],
                 [3, 9117.2002, -3247.599],
                 [3, 9190.7998, -3128.0]]
                
                
    # print "!(Hpk==2 && abs(Hu-10239.599) < 0.1 && abs(Hv-2235.5998) < 0.1) && !(Hpk==3 && abs(Hu-9117.2002) < 0.1 && abs(Hv-(-3247.599)) < 0.1) && !(Hpk==3 && abs(Hu-9190.7998) < 0.1 && abs(Hv-(-3128)) < 0.1)"

    cuts = []
    for blk in blacklist:
        cuts.append("!(Hpk==%(plane)d && abs(Hu-(%(Hu)10.4f)) < %(eps)1.10f && abs(Hv-(%(Hv)10.4f)) < %(eps)1.10f )" % {"plane" : blk[0], "Hu" : blk[1], "Hv" : blk[2], "eps" : 0.01})
    
    
    acceptance_cut_1 = "(Hpk==1 && (Hu > -9000))"
    acceptance_cut_2 = "(Hpk==2 && (Hu > -8000))"
    acceptance_cut_3 = "(Hpk==3 && (Hu > -12000))"
    acceptance_cut_4 = "(Hpk==4 && (Hu > -7000))"
    
    # cut = "Hevt==%d && (%s) && (%s)" % (evt,  " && ".join(cuts),  " || ".join([acceptance_cut_1, acceptance_cut_2, acceptance_cut_3, acceptance_cut_4]))
    cut = "Hevt==%d" % evt
    # print cut

    t.SetEstimate(t.GetEntries())
    t.Draw(sel,cut,"goff")
    N = t.GetSelectedRows()
    Hu = t.GetV1()
    Hv = t.GetV2()
    plane = t.GetV3()

    # t.Scan(sel,cut)
    
    allhits = []
    perplanehits = [[],[],[],[]]
    for i in xrange(N):
        dhit = DetectorHit(int(plane[i]), Hv[i], Hu[i])
        allhits.append(dhit)
        perplanehits[int(plane[i])-1].append(dhit)
        
    return allhits, perplanehits
    
Nevt = 500
Sevt = 5000

na63 = NA63Detector()
na63.createEve()

intersects_entry = []
for i in xrange(4):
    
    intersects1 =  TEvePointSet("intersects_entry%d" % (i + 1))
    intersects1.SetMarkerColor(kOrange)
    intersects1.SetMarkerSize(0.8)
    intersects1.SetMarkerStyle(4);
    intersects_entry.append(intersects1)


hit_planes = [TH2D("hit_plane%d" % ((1+i)), "Hits in Plane %d" % ((1+i)), na63.NPIXELS_MIMOSA_X, -na63.NPIXELS_MIMOSA_X, na63.NPIXELS_MIMOSA_X, na63.NPIXELS_MIMOSA_Y, -na63.NPIXELS_MIMOSA_Y, na63.NPIXELS_MIMOSA_Y) for i in xrange(4)]

for evt_i in xrange(Sevt,Sevt+Nevt):
    # print na63.POS_MIMOSA_ONE_Z
    hits, perplanehits = getPlanes(evt_i)
    # print hits

    if evt_i % 100 == 0:
        print "%2.0f%% done..." % (float(evt_i-Sevt)/float(Nevt)*100.0)

        
    for hit in hits:
        # print "det", hit.det
        det_pos = na63.get_center_of_detector(hit.det)
        intersects_entry[hit.det-1].SetNextPoint(det_pos[0]+(hit.u/(11520.0)*na63.DIM_MIMOSA_X/2),det_pos[1]+(hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2),det_pos[2])
        
        hit_planes[hit.det-1].Fill(hit.u/10.0, hit.v/10.0)

    # if len(perplanehits[3]) > 0 and len(perplanehits[2]) > 0:
    #     for p4hit in perplanehits[3]:
    #         for p3hit in perplanehits[2]:
    #             xes =  TEveLine("34line")
    #             xes.SetLineColor(kWhite)    
    #             
    #             xes.SetNextPoint(p4hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p4hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p4hit.det)[2])
    #             xes.SetNextPoint(p3hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p3hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p3hit.det)[2])
    #             gEve.AddElement(xes)
    # 
    # 
    # 
    # if len(perplanehits[3]) > 0 and len(perplanehits[0]) > 0:
    #     for p4hit in perplanehits[3]:
    #         for p3hit in perplanehits[0]:
    #             xes =  TEveLine("31line")
    #             xes.SetLineColor(kBlue)    
    # 
    #             xes.SetNextPoint(p4hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p4hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p4hit.det)[2])
    #             xes.SetNextPoint(p3hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p3hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p3hit.det)[2])
    #             gEve.AddElement(xes)
    # 
    # 
    # if len(perplanehits[0]) > 0 and len(perplanehits[1]) > 0:
    #     for p4hit in perplanehits[0]:
    #         for p3hit in perplanehits[1]:
    #             xes =  TEveLine("01line")
    #             xes.SetLineColor(kGreen)    
    # 
    #             xes.SetNextPoint(p4hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p4hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p4hit.det)[2])
    #             xes.SetNextPoint(p3hit.u/(11520.0)*na63.DIM_MIMOSA_X/2, p3hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2, na63.get_center_of_detector(p3hit.det)[2])
    #             gEve.AddElement(xes)
    #             
    if len(perplanehits[0]) > 0 and len(perplanehits[1]) and len(perplanehits[2]) > 0 and len(perplanehits[3]) > 0:
        for p1hit in perplanehits[0]:
            for p2hit in perplanehits[1]:
                for p3hit in perplanehits[2]:
                    for p4hit in perplanehits[3]:
                    
                        xes =  TEveLine("all")
                        xes.SetLineColor(kGreen)    
                        
                        
                        dp3 = na63.get_center_of_detector(p3hit.det)
                        dp4 = na63.get_center_of_detector(p4hit.det)
                        dp1 = na63.get_center_of_detector(p1hit.det)
                        dp2 = na63.get_center_of_detector(p2hit.det)
                        
                        
                        xes.SetNextPoint(dp3[0]+(p3hit.u/(11520.0)*na63.DIM_MIMOSA_X/2), dp3[1] + (p3hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2), dp3[2])
                        xes.SetNextPoint(dp4[0]+(p4hit.u/(11520.0)*na63.DIM_MIMOSA_X/2), dp4[1] + (p4hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2), dp4[2])
                        xes.SetNextPoint(dp1[0]+(p1hit.u/(11520.0)*na63.DIM_MIMOSA_X/2), dp1[1] + (p1hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2), dp1[2])
                        xes.SetNextPoint(dp2[0]+(p2hit.u/(11520.0)*na63.DIM_MIMOSA_X/2), dp2[1] + (p2hit.v/(5760.0)*na63.DIM_MIMOSA_Y/2), dp2[2])
                        gEve.AddElement(xes)
    
    
    
    # raw_input("meh")
    # for inter in intersects_entry: gEve.RemoveElement(inter, 0)
    
    # xes =  TEveLine("q/p = %d/%2.2f" % (int(q),np.linalg.norm(p)))
    # xes.SetLineColor(color_charge[q > 0])
    # xes.SetMarkerSize(0.2)
    # xes.SetMarkerStyle(4);
    # 
    # 
    # intersects_entry[hit[3]-1].SetNextPoint(hit[0], hit[1], hit[2])
    
    
    
    # 
    # ca = TCanvas("meh", "meh", 1000, 1000)
    # ca.Divide(2,2)
    # 
    # if p1.GetEntries() > 0:
    #     ca.cd(1)
    #     p1.Draw("Text")
    # 
    # if p2.GetEntries() > 0:
    #     ca.cd(2)
    #     p2.Draw("Text")
    # 
    # if p3.GetEntries() > 0:
    #     ca.cd(3)
    #     p3.Draw("Text")
    #     
    # if p4.GetEntries() > 0:
    #     ca.cd(4)
    #     p4.Draw("Text")
    #     
for inter in intersects_entry: gEve.AddElement(inter)
gEve.Redraw3D()


cc = TCanvas("hits", "hit planes", 1000, 1000)
cc.Divide(2,2)
for i,h in enumerate(hit_planes):
    cc.cd(i+1)
    h.Draw("TEXT")
    h.SetContour(888)
    
    # h.SetMarkerColor(kRed)
    # h.SetFillColor(kRed)
    # h.SetMarkerSize(1.0)
    # h.SetMarkerStyle(21)


raw_input("done")






















