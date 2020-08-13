#Nicolaas van der Voort
#August 13, 2020
#AG Seidel, HHU Dusseldorf

import matplotlib.pyplot as plt
import numpy as np
import rmsd
import copy
import developmental_functions as df

################################################################################
class fourPointSet:
    """Coarsely align a four-point Origami structure
    See help of Orient function for usecase example
    currently only four-points implemented, but is is written such that it does 
    not take so much time to expand it to the general case.
    The trick is to have a structure symmetry to exploit for coarsely aligning
    the constructs. For the Origamis this is the assumption that all are flipped
    upwards. (no jelly-sandwiches here).
    """
    keys = ['D1', 'D2', 'A1', 'A2']
    
    def rotate_via_numpy(self, xy, radians):
        """Use numpy to build a rotation matrix and take the dot product.
        adopted from: https://gist.github.com/LyleScott/e36e08bfb23b1f87af68c9051f985302"""
        x, y = xy
        c, s = np.cos(radians), np.sin(radians)
        j = np.matrix([[c, s], [-s, c]])
        m = np.dot(j, [x, y])
        return float(m.T[0]), float(m.T[1])
    
    def __init__(self, loc):
        self.points = dict.fromkeys(self.keys, 0)
        self.points['D1'] = np.array([ loc['G'].spotLst[0].posx, loc['G'].spotLst[0].posy ])
        self.points['D2'] = np.array([ loc['G'].spotLst[1].posx, loc['G'].spotLst[1].posy ])
        self.points['A1'] = np.array([ loc['Y'].spotLst[0].posx, loc['Y'].spotLst[0].posy ])
        self.points['A2'] = np.array([ loc['Y'].spotLst[1].posx, loc['Y'].spotLst[1].posy ])
                                     
    def RepositionToPoint(self, pointname):
        """pointname is e.g. 'A1'"""
        dpos = copy.deepcopy(self.points[pointname])
        for point in self.points.values():
            point -= dpos
                                     
    def calcAngle(self, point1 = 'A1', point2 = 'A2'):
        Delta = self.points[point1] - self.points[point2]
        Dx, Dy = Delta
        self.angle = np.arctan(Dy / Dx)
                                     
    def rotate(self, angle = None):
        if angle is None: angle = self.angle
        #using self.points.values creates a bug where the rotation is no 
        #longer executed correctly, perhaps overwriting itself prematurely?
        for name in self.points.keys():
            self.points[name] = \
                 np.array(self.rotate_via_numpy(self.points[name], angle))
    def mirrorInX(self):
        for name in self.points.keys():
                self.points[name] *= [1,-1]

    def relabelLeftToRight(self):
        """relabels the spots such that the spots that have lowest x are '1'"""
        for color in ['D', 'A']:
            #0 position is for x dimension
            if self.points[color + '1'][0] > self.points[color + '2'][0]:                    
                point_holder = self.points[color +'1']
                self.points[color + '1'] = self.points[color + '2']
                self.points[color + '2'] = point_holder
                
    def isYAverageBelowZero(self, color):
        #1 position is for y dimension
        return self.points[color + '1'][1] + self.points[color + '2'][1] < 0
    def calcDisplacement(self, nameset1, nameset2):
        """
        calculates vector displacement between two sets of points
        E.g. for displacement between all acceptors and all Donors:
            nameset1 = ['D1'..'Dn']
            nameset2 = ['A1'..'An']
        """
        Nentries = len(nameset1) * len(nameset2)
        displacements = np.zeros((Nentries, 2))
        i=0
        for name1 in nameset1:
            for name2 in nameset2:
                displacements[i] = self.points[name1] - self.points[name2]
                i+=1
        return displacements
    def applyChannelShift(self, nameset, shift):
        for name in nameset:
            self.points[name] -= shift
    def Orient(self):
        """
        utility function for the most common case of aligning
        First: A1 point is chosen to align to origin
        second: the angle between points A1 and A2 is calculated
        third: the structure is rotated by this angle, the two acceptors now lie
            on the x-axis
        fourth: If the Donors are below the acceptors, the structure is rotated
            180 degrees
        fifth: Assuming that all origamis are labelled from the top: all HF 
            pairs are left and NF are right. The points are re-labelled left-
            to-right.
        sixth: the poinsets are repositioned such that the new A1 is in the 
            origin.
        """
        self.RepositionToPoint('A1')
        self.calcAngle()
        self.rotate()
        if self.isYAverageBelowZero('D'):
            self.rotate(angle = np.pi)
        self.relabelLeftToRight()
        self.RepositionToPoint('A1')
        
    def getPoints(self):
        points_out = []
        for point in self.points.values():
            points_out.append(point)
        return np.array(points_out)
    def setPoints(self, pointsetList):
        for key, pointset in zip(self.keys, pointsetList):
            self.points[key] = pointset
        return
    def plotPoints(self, pxSize, c = ('orange', 'orange', 'r', 'r')):
        plt.scatter(*self.points['D1'] * pxSize, c = c[0], marker = '.')
        plt.scatter(*self.points['D2'] * pxSize, c = c[1], marker = '.')
        plt.scatter(*self.points['A1'] * pxSize, c = c[2], marker = '.')
        plt.scatter(*self.points['A2'] * pxSize, c = c[3], marker = '.')
        return
        
        
################################################################################
class ensemblePointSet:
    """operates on an ensemble of point structures, E.g. an ensemble of origami 
    structures.
    Usage example:
    
    DADApointset = ensemblePointSet(DADAlocLst, pxSize) # initialize
    #calculate dispalcement of Donor Channel wrt Acceptor Channel
    disp = DADApointset.calcAvgDisplacement(['D1', 'D2'], ['A1', 'A2'])
    #apply this shift (other shifts can be taken as well) to all points
    #applyChannelShift is a function on the pointset level, therefore the 
    #function callBatchFun is used
    DADApointset.callBatchFun('applyChannelShift', ('D1', 'D2'), disp)
    #call the Orient functions, which consists of a series of displacements
    #and rotations to coarsely align all structures and make the acceptors lie
    #on the x axis
    DADApointset.callBatchFun('Orient')
    #generate an anchor. Here it is generated from the design model
    DADApointset.genAnchorFromModel()
    #alternative the average structure can be taken
    DADApointset.genAnchorFromMean()
    #run minimize the RMSD to the anchor for all structures
    DADApointset.batchRmsd()
    #for visualization
    _ = plt.hist(DADApointset.scores, bins = 100, range = (0,10))
    #prune (kick-out bad fits
    DADApointset.pruneByScore(2)
    #plot
    DADApointset.plotPointsets()
    """
    def __init__(self, locLst, pxSize):
        self.pointsets = []
        self.locLst = locLst # keeping the raw data for now
        for loc in locLst:
            self.pointsets.append(fourPointSet(loc))
        self.Nentries = len(self.pointsets)
        self.pxSize = pxSize
        return
    def calcAvgDisplacement(self, nameset1, nameset2, maxval = 10):
        #first axis: structs
        #second axis: pair sets
        #third axis: distances
        Npermutations = len(nameset1) * len(nameset2)
        disp = np.zeros((self.Nentries, Npermutations, 2))
        for i, ps in enumerate(self.pointsets):
            disp[i] = ps.calcDisplacement(nameset1, nameset2)
        disp = disp.reshape(self.Nentries * Npermutations, 2) # throw away struct specificity
        disp = df.kickvector(disp, maxval)
        return np.mean(disp, axis = (0))
    def callBatchFun(self, functionName, *args, **kwargs):
        """batch run an operation on the constituent pointset objects"""
        out = []
        for ps in self.pointsets:
            out.append(getattr(ps,functionName)(*args, **kwargs))
        return out #out can be anything
    def batchRmsd(self):
        """uses an external rmsd library to align all structures to an anchor"""
        self.anchor -= rmsd.centroid(self.anchor)
        self.scores = []
        for ps in self.pointsets:
            #convert to numpy array for rmsd lib
            psArr = ps.getPoints()
            #potentially the fun rmsd.kabsch_rmsd does multiple steps in one
            #but I don't understand what it does exatly
            psArr -= rmsd.centroid(psArr)
            U = rmsd.kabsch(psArr, self.anchor)
            psArr = np.dot(psArr, U)
            self.scores.append(rmsd.rmsd(psArr, self.anchor))
            ps.setPoints(psArr)
        return
    def pruneByScore(self, minscore):
        prunedIds = []
        for i in range(len(self.scores)-1, -1, -1):
            if self.scores[i] > minscore:
                self.scores.pop(i)
                self.pointsets.pop(i)
                prunedIds.append(i)
        print('pruned %i elements' %len(prunedIds))    
        return prunedIds
    def genAnchorFromMean(self):
        pointsetArr = np.array(self.callBatchFun('getPoints'))
        self.anchor = np.mean(pointsetArr, axis = 0)
        return
    def genAnchorFromModel(self, H2H = 2.7, pxSize = 10):
        """distances based on calculations made by Anders based on the number 
        of strands and helices separating the dyes"""
        A1 = np.array([0,0])
        A2 = A1 + np.array([74.59, - 5*H2H])
        D1 = A1 + np.array([0, 2 * H2H])
        D2 = A2 + np.array([0, 6 * H2H])
        self.anchor = np.array([A1, A2, D1, D2]) / pxSize
        return
    
    def plotPointsets(self, c = ['orange', 'orange', 'r', 'r'], outfile = None, \
        title = None, xlim = [-60, 60]):
        fontsize = 24
        plt.figure(figsize=(10,10))
        self.callBatchFun('plotPoints', self.pxSize)
        if hasattr(self, 'anchor'):
            anchor = copy.deepcopy(self.anchor) #localbound var
            anchor *= self.pxSize
            plt.scatter(anchor[:,0], anchor[:,1], c='k', marker = '+', s = 100)
        plt.axis('equal')
        plt.grid(True)
        plt.tick_params(labelsize= fontsize)
        plt.xlabel('x (nm)', size = fontsize)
        plt.ylabel('y (nm)', size = fontsize)
        plt.xlim(xlim)
        plt.ylim((-25, 25))
        plt.title(title, size = fontsize)
        if outfile:
            plt.savefig(outfile, bbox_inches = 'tight', dpi = 300)
        plt.show()
        return