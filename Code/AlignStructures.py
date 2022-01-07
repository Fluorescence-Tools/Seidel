#Nicolaas van der Voort
#August 13, 2020
#AG Seidel, HHU Dusseldorf

import matplotlib.pyplot as plt
import numpy as np
import rmsd
import copy
import developmental_functions as df
from matplotlib import patches

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
    

        
    def __init__(self, pointArr, loc = None):
        #save loc information for later filtering
        self.loc = loc
        assert(len(self.keys) == pointArr.shape[0])
        self.points = {}
        for key, point in zip(self.keys, pointArr):
            self.points[key] = point
        
    @classmethod
    def fromLoc(cls, loc):
        D1 = np.array([ loc['G'].spotLst[0].posx, loc['G'].spotLst[0].posy ])
        D2 = np.array([ loc['G'].spotLst[1].posx, loc['G'].spotLst[1].posy ])
        A1 = np.array([ loc['Y'].spotLst[0].posx, loc['Y'].spotLst[0].posy ])
        A2 = np.array([ loc['Y'].spotLst[1].posx, loc['Y'].spotLst[1].posy ])
        pointArr = np.array([D1, D2, A1, A2])
        return cls(pointArr, loc)
                                     
    def RepositionToPoint(self, pointname):
        """pointname is e.g. 'A1'"""
        dpos = copy.deepcopy(self.points[pointname])
        for point in self.points.values():
            point -= dpos
                                     
    def calcAngle(self, point1 = 'A1', point2 = 'A2'):
        Delta = self.points[point2] - self.points[point1]
        Dx, Dy = Delta
        self.angle = np.arctan2(Dy, Dx)
                                     
    def rotate(self, angle = None):
        if angle is None: angle = self.angle
        #using self.points.values creates a bug where the rotation is no 
        #longer executed correctly, perhaps overwriting itself prematurely?
        for name in self.points.keys():
            self.points[name] = \
                 np.array(self.rotate_via_numpy(self.points[name], angle))
    def mirror(self, mirrorAxis = 'horizontal'):
        if mirrorAxis == 'horizontal':
            for name in self.points.keys():
                self.points[name] *= [1, -1]
        elif mirrorAxis == 'vertical':
            for name in self.points.keys():
                self.points[name] *= [-1, 1]
    
    def swapPoints(self, point1, point2):
        point_holder = self.points[point1]
        self.points[point1] = self.points[point2]
        self.points[point2] = point_holder
    def relabelLeftToRight(self):
        """relabels the spots such that the spots that have lowest x are '1'"""
        for color in ['D', 'A']:
            #0 position is for x dimension
            if self.points[color + '1'][0] > self.points[color + '2'][0]:
                self.swapPoints(color + '1', color + '2')
                
    def relabelOnFRETIndicator(self, FRETind):
        """FRETind can be e.g., 'tauG' or 'proxRatio'. 
        The spot pair with the higher value obtains index 2.
        E.g. FRETind = 'proxRatio', then the highFRET pair is labelled 2.
        The whole set of FRETindicators is called 'FRETind', confusing."""
        FRETindicators = [getattr(spot, FRETind) for spot in self.loc['FRETind']]
        maxid = np.argmax(FRETindicators)
        if maxid == 0: #the higher FRET indicator has label '1', swap labels
            self.swapPoints('D1', 'D2')
            self.swapPoints('A1', 'A2')
        elif maxid == 1: #the higher FRET indicator has label '2', do nothing
            pass
        
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
    def Orient(self, isMirrorInY = False):
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
        if isMirrorInY:
            self.mirror(mirrorAxis = 'vertical')
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
    def getDistance(self, name1, name2):
        return np.linalg.norm(self.getDisplacement(name1, name2))
    def getDisplacement(self, name1, name2):
        return self.points[name2] - self.points[name1]
    def getFourPointAngle(self, name1, name2, name3, name4):
        """
        calculates angle between AB and CD vectors where points ABCD are denoted
        by name1, name2, name3, name 4 respectively.
        Rotates the whole structure such that AB has angle 0. Then the angle of 
        CD denotes the angle of AB and CD.
        """
        self.calcAngle(name1, name2)
        self.rotate()
        CD = self.getDisplacement(name3, name4)
        #arctan2 takes first y, then x.
        return np.arctan2(CD[1], CD[0])
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
        #self.locLst = locLst # keeping the raw data for now
        for loc in locLst:
            self.pointsets.append(fourPointSet.fromLoc(loc))
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
        anchor = self.anchor.getPoints() #convert to np arr
        anchor -= rmsd.centroid(anchor)
        self.scores = []
        for ps in self.pointsets:
            #convert to numpy array for rmsd lib
            psArr = ps.getPoints()
            #potentially the fun rmsd.kabsch_rmsd does multiple steps in one
            #but I don't understand what it does exatly
            psArr -= rmsd.centroid(psArr)
            U = rmsd.kabsch(psArr, anchor)
            psArr = np.dot(psArr, U)
            self.scores.append(rmsd.rmsd(psArr, anchor))
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
        pointArr = np.mean(pointsetArr, axis = 0)
        self.anchor = fourPointSet(pointArr)
        return
    def genAnchorFromModel(self, H2H = 2.7, alongHelix = 74.8):
        """distances based on calculations made by Anders based on the number 
        of strands and helices separating the dyes"""
        A1 = np.array([0,0])
        #A2 = A1 + np.array([74.59, - 5*H2H])
        #based on A-A any No. of G spot as calibration A-A: 75.8 corrected for 
        #12nm y displacement: 74.8
        A2 = A1 + np.array([alongHelix, - 5*H2H])
        D1 = A1 + np.array([0, 2 * H2H])
        D2 = A2 + np.array([0, 6 * H2H])
        pointArr = np.array([D1, D2, A1, A2]) / self.pxSize
        self.anchor = fourPointSet(pointArr)
        return
    def genAnchorFromPointset(self, pointsetID):
        self.anchor = self.pointsets[pointsetID]
        return
    def selectOnPointsetIDs(self, pointsetIDs):
        self.pointsets = [self.pointsets[ID] for ID in pointsetIDs]
        return
    def plotPointsets(self, c = ['orange', 'orange', 'r', 'r'], outfile = None, \
        title = None, xlim = [-60, 60], ylim = [-25,25]):
        fontsize = 24
        plt.figure(figsize=(10,10))
        self.callBatchFun('plotPoints', self.pxSize)
        if hasattr(self, 'anchor'):
            anchor = self.anchor.getPoints() #localbound var
            anchor *= self.pxSize
            plt.scatter(anchor[:,0], anchor[:,1], c='k', marker = '+', s = 100)
        # plot center for each point
        points = np.array(self.callBatchFun('getPoints')) * self.pxSize
        meanpoints = np.mean(points, axis = 0)
        plt.scatter(meanpoints[:,0], meanpoints[:,1], marker = 's', c = c, \
        edgecolor = 'k')
        plt.axis('equal')
        plt.grid(True)
        plt.tick_params(labelsize= fontsize)
        plt.xlabel('x (nm)', size = fontsize)
        plt.ylabel('y (nm)', size = fontsize)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.title(title, size = fontsize)
        if outfile:
            plt.savefig(outfile, bbox_inches = 'tight', dpi = 300)
        #plt.show()
        return plt.gcf()
        
    def plotPairDistance(self, pointname1, pointname2, 
        bins = np.arange(0,100,2), alpha = 1):
        dist = np.array(self.callBatchFun('getDistance', pointname1, pointname2))
        dist *= self.pxSize
        plt.hist(dist, bins = bins, alpha = alpha)
        return dist
    def plotAngle(self, pointname1, pointname2, pointname3, pointname4, **axKwargs):
        #get all data
        d, d_mean, d_std, d_meanstd, d_model, polar, polar_mean, polar_std, polar_meanstd, polar_model= \
            self.getPointPairStats(pointname1, pointname2, pointname3, pointname4)
        #plot points
        r, phi = polar.T
        p = plt.polar(phi, r, '.', **axKwargs)
        #plot mean
        r_mean, phi_mean = polar_mean
        r_meanstd, phi_meanstd = polar_meanstd
        darkc = lighten_color(p[0].get_color(), 1.2)
        plt.polar(phi_mean, r_mean, '+', 
            c = darkc, 
            markersize = 10, markeredgewidth = 2)
        
        #plot model expectation value
        r_model, phi_model = polar_model
        plt.polar(phi_model, r_model, '*',
            markersize = 10, markeredgewidth = 1,
            c=darkc )
            
        return
    def getPointPairStats(self, axisPoint1, axisPoint2, pairPoint1, pairPoint2,
        verbose = False, title = 'some distance', addangle = 0):
        #isse:checked again and it is unclear which approach is best. Leave as is for now.
            #issue: it was discovered that it is better to use overall RMSD alignment rather
            #then aligning on A1-A2 (I guess it uses more statistics).
            #So the axispoints can go out alltogether.
        #orient the structure according to two axisPoints on the x axis
        self.callBatchFun('calcAngle', axisPoint1, axisPoint2)
        #I have to mangle my code to make a pretty figure
        #self.callBatchFun('rotate')
        for pointset in self.pointsets:
            pointset.rotate(angle = pointset.angle + addangle)
        #get distances
        d = self.callBatchFun('getDisplacement', pairPoint1, pairPoint2)
        d = np.array(d) * self.pxSize
        d_mean = np.mean(d, axis = 0)
        d_std = np.std(d, axis = 0)
        d_meanstd = d_std / np.sqrt(d.shape[0])
        #unsure what the correct treatment in polar coordinates is. Have asked Oleg.
        polar = np.array([[np.linalg.norm( el ), np.arctan2(el[1], el[0])]
                     for el in d])
        # polar_mean = [np.linalg.norm( d_mean ), 
            # np.rad2deg(np.arctan2(d_mean[1], d_mean[0]))]
        # polar_std = np.std(polar, axis = 0)
        # polar_meanstd = polar_std / np.sqrt(d.shape[0])
        polar_mean = cart2pol(*d_mean)
        #take perpendicular component of std, divide by r_mean.
        perpendicular_std = np.linalg.norm([d_std[0] * np.sin(polar_mean[1]), d_std[1] * np.cos(polar_mean[1])])
        phi_std = np.arctan2(perpendicular_std , polar_mean[0])
        r_std = np.linalg.norm([d_std[0] * np.cos(polar_mean[1]), d_std[1] * np.sin(polar_mean[1])])
        polar_std = [r_std, phi_std]
        polar_meanstd = polar_std / np.sqrt(d.shape[0])
        
        #calc model points
        self.anchor.calcAngle(axisPoint1, axisPoint2)
        self.anchor.rotate()
        d_model = self.anchor.getDisplacement(pairPoint1, pairPoint2) * self.pxSize
        polar_model = cart2pol(*d_model)
        if verbose:
            print(title)
            print('(x,y)')
            print('mean is (%.2f, %.2f)' % (d_mean[0], d_mean[1]))
            print('uncertainty of mean is (%.2f, %.2f)'% (d_meanstd[0], d_meanstd[1]))
            print('model is (%.2f, %.2f)' % (d_model[0], d_model[1]))
            print('spread is (%.2f, %.2f)\n'% (d_std[0], d_std[1]))
            #print('!!! all polar statistics transformed from cartesian space!!!')
            #print('!!!approach unverified!!!')
            #print('(r,phi)')
            #print('mean is (%.2f, %.2f)' % (polar_mean[0], np.rad2deg(polar_mean[1])))
            #print('uncertainty of mean is (%.2f, %.2f)'% (polar_meanstd[0], np.rad2deg(polar_meanstd[1])))
            #print('model is (%.2f, %.2f)' % (polar_model[0], np.rad2deg(polar_model[1])))
            #print('spread is (%.2f, %.2f)\n'% (polar_std[0], np.rad2deg(polar_std[1])))
        #might want to return the variables later
        return d, d_mean, d_std, d_meanstd, d_model, polar, polar_mean, polar_std, polar_meanstd, polar_model
        
    def savePointsetsToCsv(self, ffile):
        """saves all coordinates, distances and model data"""
        #probably want to move definition of keys to ensemble level
        keys = self.pointsets[0].keys
        header = ''
        for key in keys:
            header += key + 'x'
            header += ','
            header += key + 'y'
            header += ','
        pointout = np.array(self.callBatchFun('getPoints'))
        shape = pointout.shape
        pointout = pointout.reshape([shape[0], shape[1]* shape[2]])
        pointout *= self.pxSize
        np.savetxt(ffile, pointout, delimiter = ',', header = header, fmt = '%.2e')

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])