#!/usr/bin/env python2
import sys
import time
import copy
import signal
import json
import numpy as np
import matplotlib.cm as cm
import rat
import argparse
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import ctypes as ct
from PyQt4 import QtGui, QtCore

class MainWindow(QtGui.QWidget):
    '''
    Order of operations:
    1. Create the main window
    2. Create Grid layout
    3. Open the ratroot file to get information and store event info
    4. 
    '''
    App = None
    plot = None
    colorMask = None
    isLog = False

    def __init__(self, fname, app=None):
        self.fname = fname
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        super(MainWindow,self).__init__()
        self.beginRat()
        # Set the layout
        self.buildLayout()
        self.show()

    def buildLayout(self):
        #self.layout = QtGui.QFormLayout()
        formGroupBox = QtGui.QGroupBox()
        formLayout = QtGui.QFormLayout()
        # Main layout
        self.layout = QtGui.QGridLayout()
        # Events
        self.total_events_line = QtGui.QLineEdit("0")
        self.total_events_line.setReadOnly(True)
        self.total_events_line.setText(str(self.entries))
        #self.layout.addRow(QtGui.QLabel("Total Events"), self.total_events_line)
        formLayout.addRow(QtGui.QLabel("Total Events"), self.total_events_line)
        ## Event switcher buttons
        left_button = QtGui.QPushButton("<")
        left_button.clicked.connect(self.leftEvent)
        right_button = QtGui.QPushButton(">")
        right_button.clicked.connect(self.rightEvent)
        self.event = QtGui.QLineEdit("0")
        self.event.textChanged.connect(self.newEvent)
        self.hboxchooser = QtGui.QHBoxLayout()
        self.hboxchooser.addWidget(left_button)
        self.hboxchooser.addWidget(right_button)
        self.hboxchooser.addWidget(self.event)
        formLayout.addRow(self.hboxchooser)
        # Color masks on QHL, QHS, TAC
        QHL_button = QtGui.QPushButton("QHL")
        QHL_button.clicked.connect(self.QHLMask)
        QHS_button = QtGui.QPushButton("QHS")
        QHS_button.clicked.connect(self.QHSMask)
        TAC_button = QtGui.QPushButton("TAC")
        TAC_button.clicked.connect(self.TACMask)
        LOG_button = QtGui.QPushButton("LOG")
        LOG_button.clicked.connect(self.LOGMask)
        MUON_button = QtGui.QPushButton("Muon Fitter")
        MUON_button.clicked.connect(self.findMuon)
        MATH_button = QtGui.QPushButton("Math!")
        MATH_button.clicked.connect(self.timingCorrection)
        hboxcolor = QtGui.QHBoxLayout()
        hboxcolor.addWidget(QHL_button)
        hboxcolor.addWidget(QHS_button)
        hboxcolor.addWidget(TAC_button)
        hboxcolor.addWidget(MUON_button)
        hboxcolor.addWidget(MATH_button)
        #hboxcolor.addWidget(LOG_button)
        formLayout.addRow(hboxcolor)
        # Muon Fitter Button
        ##
        formGroupBox.setLayout(formLayout)
        self.layout.addWidget(formGroupBox,0,0)
        # Main graphics here
        #self.glWin = gl.GLViewWidget()
        #self.newEvent()
        self.plot = gl.GLScatterPlotItem()
        self.newEvent()
        self.glWin = PlotObject(self.plot, self.App)
        self.layout.addWidget(self.glWin,1,0,4,1)
        #self.hboxPlotter = QtGui.QHBoxLayout()
        #self.hboxPlotter.addWidget(self.glWin)
        #self.layout.addRow(self.hboxPlotter)
        #self.layout.addWidget(self.glWin)
        # Manipulation buttons here
        self.setLayout(self.layout)
        # Startup

    def leftEvent(self):
        min_event = 0
        cur_event = int(self.event.text())
        if cur_event > (min_event):
            cur_event -= 1
        self.event.setText(str(cur_event))

    def rightEvent(self):
        max_event = int(self.total_events_line.text())
        cur_event = int(self.event.text())
        if cur_event < (max_event-1):
            cur_event += 1
        self.event.setText(str(cur_event))

    def newEvent(self):
        min_event, max_event = 0, int(self.total_events_line.text())
        try:
            cur_event = int(self.event.text())
        except ValueError:
            self.event.setText("0")
            cur_event = int(self.event.text())
        if (cur_event >= min_event) and (cur_event < max_event):
            self.updateEvent()
        else:
            self.event.setText("0")

    def beginRat(self):
        # Load in the rat root file dsreader
        self.ds = rat.RAT.DU.DSReader(self.fname)
        self.pmt_info = rat.utility().GetPMTInfo()
        self.entries = self.ds.GetEntryCount()

    def readRat(self):
        event = int(self.event.text())
        self.ev = self.ds.GetEntry(event).GetEV(0)
        calpmt = self.ev.GetCalPMTs()
        npmt = calpmt.GetNormalCount()
        allpmtcount = self.pmt_info.GetCount()
        # Loop through the pmts for this event and get their coordinates
        # as well as their qhl, qhs, and timing tac
        # store these in numpy arrays
        self.posArray = np.zeros((allpmtcount, 3))
        self.pmtidArray = np.zeros((allpmtcount))
        self.QHL = np.zeros((allpmtcount)) - 100
        self.QHS = np.zeros((allpmtcount)) - 100
        self.TAC = np.zeros((allpmtcount)) - 100
        for pmt_id in range(allpmtcount):
            pmtx = self.pmt_info.GetPosition(pmt_id).X()
            pmty = self.pmt_info.GetPosition(pmt_id).Y()
            pmtz = self.pmt_info.GetPosition(pmt_id).Z()
            self.pmtidArray[pmt_id] = pmt_id
            self.posArray[pmt_id] = np.array([pmtx, pmty, pmtz])

        self.plWeights = np.zeros((allpmtcount)) + 3 

        for pmt in range(npmt):
            thispmt = calpmt.GetPMT(pmt)
            pmt_id = thispmt.GetID()
            self.QHL[pmt_id] = thispmt.GetQHL()
            self.QHS[pmt_id] = thispmt.GetQHS()
            self.TAC[pmt_id] = thispmt.GetTime()
            self.plWeights[pmt_id] = 10

    def readRatOld(self):
        event = int(self.event.text())
        self.ev = self.ds.GetEntry(event).GetEV(0)
        calpmt = self.ev.GetCalPMTs()
        npmt = calpmt.GetNormalCount()
        # Loop through the pmts for this event and get their coordinates
        # as well as their qhl, qhs, and timing tac
        # store these in numpy arrays
        self.posArray = np.zeros((npmt, 3))
        self.pmtidArray = np.zeros((npmt))
        self.QHL = np.zeros((npmt))
        self.QHS = np.zeros((npmt))
        self.TAC = np.zeros((npmt))
        for pmt in range(npmt):
            thispmt = calpmt.GetPMT(pmt)
            pmt_id = thispmt.GetID()
            pmtx = self.pmt_info.GetPosition(pmt_id).X()
            pmty = self.pmt_info.GetPosition(pmt_id).Y()
            pmtz = self.pmt_info.GetPosition(pmt_id).Z()
            self.pmtidArray[pmt] = pmt_id
            self.posArray[pmt] = np.array([pmtx, pmty, pmtz])
            self.QHL[pmt] = thispmt.GetQHL()
            self.QHS[pmt] = thispmt.GetQHS()
            self.TAC[pmt] = thispmt.GetTime()

    def addPMT(self, pmtpos):
        self.posArray = np.append(self.posArray, np.array([pmtpos]), axis=0)

    def drawEvent(self):
        self.colorize()
        # If this is the first draw, then create the object
        if self.plot is None:
            self.plot = gl.GLScatterPlotItem(pos=self.posArray, color=self.colorArray, pxMode=False)
        else:
            self.plot.pos = self.posArray
            self.plot.color=self.colorArray
            self.plot.size = self.plWeights
            self.plot.update()

    def QHSMask(self):
        self.colorMask = 'QHS'
        self.drawEvent()

    def QHLMask(self):
        self.colorMask = 'QHL'
        self.drawEvent()

    def TACMask(self):
        self.colorMask = 'TAC'
        self.drawEvent()

    def LOGMask(self):
        if self.isLog:
            self.isLog = False
        else:
            self.isLog = True
        self.drawEvent()

    def colorize(self):
        if self.colorMask == 'QHS':
            variable = self.QHS
        elif self.colorMask == 'QHL':
            variable = self.QHL
        elif self.colorMask == 'TAC':
            variable = self.TAC
        elif self.colorMask == 'MUON':
            variable = self.fOfq
        else:
            variable = self.QHL #default
        if self.isLog:
            variable = np.log(variable)
        true_min = 0.01
        min_var = max(true_min, min(variable))
        max_var = max(variable)
        var_range = max_var - min_var
        colorList = []
        for val in variable:
            if val < true_min:
                val = true_min
            try:
                colorList.append(cm.jet(int((val-min_var)*(256/var_range))))
            except ValueError:
                print 'value,min_value,value_range:',val, min_var, var_range
        self.colorArray = np.array(colorList)

    def updateEvent(self):
        # read the event
        self.readRat()
        # plot the event
        self.drawEvent()

    def findMuon(self):
        '''
        Requires an exit point to already be clicked, if there is not one
        simply do nothing. Given an exit, do the following.
        1. Loop over ALL other (hit) pmts
            - Project cone from pmt to exit
            - f(q) = Sum charge of all PMT in cone
        2. Take PMT with highest f(q) as entrance and print its ID and coords
        3. Highlight that PMT and draw a line through the detector
            - Bonus, also draw the cone?
        '''
        useSlow = False
        print 'Begin looking for muons at %s' % time.asctime()
        stime = time.time()
        spotvecpy = np.array([0,0,0])
        if useSlow:
            bestpmt, bestpmtidx = self.badReplaceWithCtypes()
        else:
            ## Using ctypes
            npmt = len(self.posArray)
            libmu = ct.CDLL('libmuviewer.so')
            pmtfinder = libmu.findBestPMT
            pdub = ct.POINTER(ct.c_double)
            cint = ct.c_int
            cdub = ct.c_double
            pmtfinder.argtypes = [pdub, pdub, pdub, pdub, cdub, cdub, cdub, cint]
            #pmtfinder.restype = np.ctypeslib.ndpointer(dtype=ct.c_double, shape=(npmt,))
            pmtfinder.restype = pdub
            cx, cy, cz = self.glWin.Candidates[0]
            #px, py, pz = self.posArray.transpose()
            #print "type of qhs: %s" % type(self.QHS[0])
            #for i in range(4):
            #    print "python: charge %f, px %f, py %f, pz %f" % (self.QHS[i], self.posArray[i][0], self.posArray[i][1], self.posArray[i][2])
            #print self.posArray[0:4]
            #px, py, pz = self.posArray.transpose()
            px = np.array([val[0] for val in self.posArray])
            py = np.array([val[1] for val in self.posArray])
            pz = np.array([val[2] for val in self.posArray])
            chrg = np.array([val for val in self.QHS])
            #for i in range(4):
            #    print "python2: charge %f, px %f, py %f, pz %f" % (chrg[i], px[i], py[i], pz[i])
            #chrg = np.array([val for val in self.QHS])
            #print "len px,py,pz,chrg: %i, %i, %i, %i" % (len(px), len(py), len(pz), len(chrg))
            ## Lets test other charges
            px_p = px.ctypes.data_as(pdub)
            py_p = py.ctypes.data_as(pdub)
            pz_p = pz.ctypes.data_as(pdub)
            chrg_p = chrg.ctypes.data_as(pdub)
            indata = (ct.c_double *npmt)()

            indata = pmtfinder(chrg_p, px_p, py_p, pz_p, cx, cy, cz, npmt)
            array_holder = np.array(np.fromiter(indata, dtype=np.float64, count=npmt))
            #for i in range(10):
            #    print array_holder[i]
            self.fOfq = []
            max_fOfq = np.max(array_holder)
            for f in array_holder:
                self.fOfq.append((f/max_fOfq)**20)
            #for i in range(10):
            #    print self.fOfq[i]
            ## Test the new program, which finds a better entrance
            pointfinder = libmu.findBestSpot
            #pointfinder = libmu.findNormBestSpot
            pointfinder.argtypes = [pdub, pdub, pdub, pdub, cdub, cdub, cdub, cint]
            pointfinder.restype = pdub
            spotvector = (ct.c_double * 3)()
            spotvector = pointfinder(chrg_p, px_p, py_p, pz_p, cx, cy, cz, npmt)
            spotvecpy = np.array(np.fromiter(spotvector, dtype=np.float64, count=3))
            print 'High precision entry:', spotvecpy
        #print self.fOfq
        self.exit_point = self.glWin.Candidates[0]
        self.entry_point = spotvecpy
        self.colorMask = 'MUON'
        self.drawEvent()
        bestpmtidx = np.argmax(self.fOfq)
        #print bestpmt
        etime = time.time()
        print 'Done looking at %s' % time.asctime()
        print 'The search took %i seconds' % (etime-stime)
        #self.plot.color[bestpmtidx] = np.array([1,1,1,1])
        print "%i = %i ?" % (len(self.colorArray), len(self.posArray))
        self.addPMT( spotvecpy )
        #self.plot.color = np.append( self.plot.color, np.array([1,1,1,1]) )
        self.colorArray = np.append( self.colorArray, np.array([np.array([1,1,1,1])]), axis=0)
        self.plWeights = np.append( self.plWeights, np.array([10]), axis=0 )
        print "%i = %i ?" % (len(self.colorArray), len(self.posArray))
        self.plot.color = self.colorArray
        self.plot.size = self.plWeights
        self.plot.pos = self.posArray
        #self.plot.color = np.append( self.colorArray, np.array([np.array([1,1,1,1])]), axis=0)
        #self.plot.pos = np.append( self.posArray, np.array([spotvecpy]))
        self.plot.update()

    def timingCorrection(self):
        muIn = self.entry_point
        muOut = self.exit_point
        # Create a table which stores timing corrections
        # Define some linear algebra
        cross = np.cross
        dot = np.vdot
        mag = np.linalg.norm
        c = 299792458.0 * 1e3 # mm / s
        n = 1.33
        theta = np.arccos(1/n)
        # Fill a database file with: key: gitd, value: {pmt_id: correction factor}
        for idx,pmt in enumerate(self.posArray):
            muVec = (muOut - muIn)
            muHat = muVec / mag(muVec)
            costheta = dot(pmt,muHat) / mag(pmt)
            if costheta > np.cos(theta):
                conVec = cross(pmt,muHat)-cross(muIn,muHat)
                t2 = (n*mag(conVec))/(c*np.sin(theta))
                t1 = dot(pmt,muHat)/c - dot(muIn,muHat)/c - t2*np.cos(theta)/n
                real_time = t2 + t1
                print real_time*1e9, self.TAC[idx], ' --- ', self.TAC[idx]-real_time*1e9
        

    def badReplaceWithCtypes(self):
        QHS = self.QHS
        posArray = self.posArray
        self.fOfq = []
        pout = self.glWin.Candidates[0]
        costheta = 1/1.33
        for pin in posArray:
            charge = 0
            v1 = pout - pin
            if np.dot(v1,v1) > 0:
                for pi,q in zip(posArray, QHS):
                    v2 = pi - pin
                    if np.dot(v2,v2) > 0:
                        comp = np.dot(v1,v2)/(np.dot(v1,v1)**0.5 * np.dot(v2,v2)**0.5)
                        if comp > costheta:
                            charge += q
            self.fOfq.append(charge)
        #self.colorMask = 'MUON'
        #self.drawEvent()
        return posArray[np.argmax(self.fOfq)], np.argmax(self.fOfq)



class PlotObject(gl.GLViewWidget):
    App = None
    can_idx = []
    can_clr = []
    def __init__(self, plot, app=None):
        if self.App is None:
            if app is not None:
                self.App = app
            else:
                self.App = QtGui.QApplication([])
        super(PlotObject,self).__init__()
        self.Poss = []
        self.Plot = plot
        self.addItem(self.Plot)
        self._downpos = []
        # init camera
        self.setCameraPosition(pos=np.array([0, 0, 0]))

    # All keyboard and mouse event handlers
    def keyPressEvent(self, ev):
        super(PlotObject,self).keyPressEvent(ev)
        if ev.key() == ord('O'):
            self.orbit(45, 45)

    def mousePressEvent(self, ev):
        ''' Store the position of the mouse press for later use '''
        super(PlotObject, self).mousePressEvent(ev)
        self._downpos = self.mousePos

    def mouseReleaseEvent(self, ev):
        ''' Allow for single click to move and right click for context menu '''
        super(PlotObject, self).mouseReleaseEvent(ev)
        if self._downpos == ev.pos():
            x = ev.pos().x()
            y = ev.pos().y()
            if ev.button() == 2:
                self.mPosition()
            elif ev.button() == 1:
                x = x - self.width() / 2
                y = y - self.height() / 2
                print self.opts['center']
                print x,y
        self._prev_zoom_pos = None
        self._prev_pan_pos = None

    def mouseMoveEvent(self, ev):
        ''' Allow shift to move and ctrl to pan '''
        shift = ev.modifiers() & QtCore.Qt.ShiftModifier
        ctrl = ev.modifiers() & QtCore.Qt.ControlModifier
        if shift:
            y = ev.pos().y()
            if not hasattr(self, '_prev_zoom_pos') or not self._prev_zoom_pos:
                self._prev_zoom_pos = y
                return
            dy = y - self._prev_zoom_pos
            def delta():
                return -dy * 5
            ev.delta = delta
            self._prev_zoom_pos = y
            self.wheelEvent(ev)
        elif ctrl:
            pos = ev.pos().x(), ev.pos().y()
            if not hasattr(self, '_prev_pan_pos') or not self._prev_pan_pos:
                self._prev_pan_pos = pos
                return
            dx = pos[0] - self._prev_pan_pos[0]
            dy = pos[1] - self._prev_pan_pos[1]
            self.pan(dx, dy, 0, relative=True)
            self._prev_pan_pos = pos
        else:
            super(PlotObject, self).mouseMoveEvent(ev)

    #def plotGLPlot(self, objs):
    #    poss = np.array([0, 0, 0])
    #    self.Poss = []
    #    self.GlobalInds = []
    #    weights = np.array(0, dtype=float)
    #    def pswc (x) : return 10 * x**0.25
    #    for obj in objs:
    #        for i,cogs in enumerate(obj.CoG):
    #            for cog in cogs:
    #                if obj.PieceWeight[i]:
    #                    poss = np.vstack([poss,np.asarray(cog.T)])
    #                    self.Poss.append(np.matrix(cog.T))
    #                    self.GlobalInds.append(obj.Index[i])
    #                    pw = pswc(obj.PieceWeight[i])
    #                    weights = np.append(weights, pw)
    #    maxw = max(weights)
    #    threshold = np.mean(weights)
    #    self.Colors = np.empy([len(weights),4])
    #    for i, pw in enumerate(weights):
    #        if pw <= threshold:
    #            c = pw / maxw
    #            self.Colors[i] = np.array([c,1,0,0.7])
    #        else:
    #            c = 1 - pw/maxw
    #            self.Colors[i] = np.array([1,c,0,0.7])
    #    self.removeItem(self.Plot)
    #    self.Plot = gl.GLScatterPlotItem()
    #    self.Plot.setData(pos=poss, size=weights, color=self.Colors, pxMode=False)
    #    self.Sizes = weights
    #    self.addItem(self.Plot)
    #    self.show()

    def mPosition(self):
        # See: 
        # Step 0: 2d viewport coordinates
        # these are x,y from -1 to 1
        mx = self._downpos.x()
        my = self._downpos.y()
        self.Candidates = []
        view_w = self.width()
        view_h = self.height()
        x = 2.0 * mx / view_w - 1.0
        y = 1.0 - (2.0 * my / view_h)
        PMi = self.projectionMatrix().inverted()[0]
        VMi = self.viewMatrix().inverted()[0]
        ray_clip = QtGui.QVector4D(x, y, -1.0, 1.0)
        ray_eye = PMi * ray_clip
        ray_eye.setW(0)
        # Convert to world coordinates
        ray_world = VMi * ray_eye
        ray_world = QtGui.QVector3D(ray_world.x(), ray_world.y(), ray_world.z())
        ray_world.normalize()
        origin = np.matrix(self.cameraPosition())
        ray_world = np.matrix([ray_world.x(), ray_world.y(), ray_world.z()])
        r = 125 # millimeters
        # Restore old colors before selecting new colors
        for idx,clr in zip(self.can_idx, self.can_clr):
            self.Plot.color[idx] = clr

        self.can_idx, self.can_clr = [], []
        for i, C in enumerate(self.Plot.pos):
            OC = origin - C
            b = np.inner(ray_world, OC)
            b = b.item(0)
            c = np.inner(OC, OC)
            c = c.item(0) - r**2
            bsqr = np.square(b)
            if (bsqr - c) >= 0:
                self.Candidates.append(self.Plot.pos[i])
                self.can_idx.append(i)
        print self.Candidates
        for i in self.can_idx:
            self.can_clr.append(copy.deepcopy(self.Plot.color[i]))
            self.Plot.color[i] = np.array([1, 1, 1, 1])
        self.Plot.update()



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('ratfile', type=str, help=('rat root file'))
    return parser.parse_args()

if __name__ == '__main__':
    signal.signal(signal.SIGINT, lambda *args: QtGui.QApplication.quit())
    args = get_args()
    app = pg.QtGui.QApplication([])
    timer = QtCore.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)
    gui = MainWindow(args.ratfile)
    sys.exit(app.exec_())
