# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 14:29:22 2021

@author: voort
"""

from PyQt5 import QtWidgets, QtCore, uic
from PyQt5.QtCore import Qt, QUrl
import pyqtgraph as pg
import sys
import os
import numpy as np
from random import randint

#further wishes:
#have like a terminal where outputs are given at the bottom, like Thomas has

#to recompile the resources qrc file, go to the folder containing the .qrc file and do
# > pyrcc5 resources.qrc -o resources.py
#see also: https://www.pythonguis.com/tutorials/qresource-system/
import resources

#scatter <- histogram <- Fit <- Support plane
#xy data is transformed into a distance distribution histogram
#the histogram is fitted
#the support plane of the fit may be investigated
#at each step one may try different parameters etc. such that multiple histograms
#may belong to one scatter dataset, but each histogram belongs to a single scatter
class Scatterdata():
    def __init__(self, GUI, fname):
        self.GUI = GUI
        self.scatterfname = fname
        self.data = np.genfromtxt(fname).astype(np.float)
        assert(self.data.shape[1] == 2) #two column data only

        #draw a scatterplot
        self.updateScatter()
        
        
    def makehist(self, *args):
        #here the scatter data is centered and then corrected.
        #also the binning is given
        #Next an histogram object is made by making an instance of the child 
        #class and copying over all attributes. Could this be done in an
        #automated way?
        raise NotImplementedError()
    def drawscatterplot(self):
        #draw a QT scatterplot
        #want to change such that it only updates the data after the first try
        self.plotdrawn = False
        if not self.plotdrawn:
            self.GUI.scatterWidget.setBackground('k')
            #there is a problem in that points keep getting added
            #but the old points are not being removed.
            #need to understand pyqtgraph.scatterWidget API better.
            #stopped here on 17 Sept 2021
            #pen = pg.mkPen(color = (255,0,0))
            self.data_line = self.GUI.scatterWidget.plot(self.shiftdata[:,0], self.shiftdata[:,1], pen = None, symbol = 'o')
        else:
            self.data_line.setData(self.shiftdata[:,0], self.shiftdata[:,1])
    def updateScatter(self):
        #get the shift set by the user
        shift = np.array([float(self.GUI.shiftx.text()), float(self.GUI.shifty.text())])
        self.shiftdata = self.data - shift
        self.drawscatterplot()
        self.GUI.scCmx.setText('in X: %.2f' % np.mean(self.shiftdata[:,0]))
        self.GUI.scCmy.setText('in Y: %.2f' % np.mean(self.shiftdata[:,1]))
        
class Histdata(Scatterdata):
    #this __init__ overloading doesn't work in python, look for a solution that does work
    def __init__(self, fname):
        #if a filename is given, then the data is loaded
        #add to list
        pass
    def __init__(self, numpyarray):
        #if a numpy array is given, build from that
        #add to list
        pass
    def addfit(self, *args):
        #return a Fit child object with the all properties copied over
        pass

class Fit(Histdata):
    def __init__(self):
        #draw PyQtgraph floatable object
        #add name to list
        #set initial data and fit variables
        #call some update function that re-draws the fit with those parameters
        pass
    def buildmodel(self, GUI, dist_type, ndist):
        #makes a new parameter object (from lmfit) with the appropriate
        #number of distribution of the correct type. E.g. 3 gaussians.
        #also update the number of inputs on the gui to match the p
        #parameters
        #need to work out how to interface with the GUI
        pass
    def fit(self):
        #optimizes the fit
        #redraws fit
        pass
    def drawFit(self):
        #draws the fit, potentially only redraws the data in the fit
        pass
    def save(self):
        #save all parameters
        pass
    
class SupportPlane(Fit):
    def __init__(self):
        #copy over stuff from child
        #add name to list
        #set buttons for which parameter to scan
        pass
    def calculateSupport(self):
        #calculates the support plane
        pass
    def abortCalculation(self):
        #aborts the running calculation, because those can take a long time
        pass
    def save():
        #save all parameters, overrides parent child function
        #also saves support plane
        pass
        
def printhello():
    print('hello')
#this class has been adapted from:
#https://learndataanalysis.org/implement-files-and-urls-to-listbox-widget-drag-and-drop-function-pyqt5-tutorial/
class dataListWidget(QtWidgets.QListWidget):
    def __init__(self, parent = None):#not sure what the parent kwarg does
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.addItem('my first list item!')
        #self.itemDoubleClicked.connect(self.currentItem().data.updateScatter())
        print(self.itemDoubleClicked.connect(self.activateCurrentItem))
    def activateCurrentItem(self):
        self.currentItem().data.updateScatter()
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()
    
    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()
        
    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
            #links = []
            for url in event.mimeData().urls():
                # https://doc.qt.io/qt-5/qurl.html
                if url.isLocalFile():
                   fname = str(url.toLocalFile())
                else:
                    fname = str(url.toString())
                #add to list
                item = QtWidgets.QListWidgetItem(fname)
                item.data = Scatterdata(self.window(), fname)
                self.addItem(item)
        else:
            event.ignore()
            
            
class MainWindow(QtWidgets.QMainWindow):
    
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        
        uic.loadUi('mainwindow.ui', self)
        
        # #init the data containing lists
        # self.scatterdatas = []
        # self.histdatas = []
        # self.fits = []
        # self.supportplanes = []
        
        #self.shiftx.textChanged.connect(#find relevant object here.updateScatter)
        
    #     # demo moving plot
    #     self.x = list(range(100))
    #     self.y = [randint(0,100) for _ in range(100)]
    #     self.scatterWidget.setBackground('w')
    #     pen = pg.mkPen(color = (255,0,0))
    #     self.data_line = self.scatterWidget.plot(self.x, self.y, pen = pen)
    #     self.timer = QtCore.QTimer()
    #     self.timer.setInterval(50)
    #     self.timer.timeout.connect(self.update_plot_data)
    #     self.timer.start()
        
    # def update_plot_data(self):
    #     self.x = self.x[1:] # remove first x element
    #     self.x.append(self.x[-1] + 1) #add a value 1 higher than the last
        
    #     self.y = self.y[1:] #remove the first y element
    #     self.y.append(randint(0,100))
        
    #     self.data_line.setData(self.x, self.y)
        
def main():        
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    #app.exec()
    sys.exit(app.exec())
    #app.quit()

if __name__ == '__main__':
    main()