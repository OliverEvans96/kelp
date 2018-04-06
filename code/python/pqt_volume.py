# pqt_volume.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Fri 08 Jul 2016 11:09:33 PM EDT
# Last Edited: Fri 08 Jul 2016 11:45:39 PM EDT

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
from numpy import *

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.opts['distance'] = 200
w.show()
w.setWindowTitle('pyqtgraph example: GLVolumeItem')

g = gl.GLGridItem()
g.scale(10, 10, 1)
w.addItem(g)

##
A = linspace(0,10,101)
L = linspace(0,10,101)
t = linspace(0,10,101)
c = [255,0,0,100]
A,L,t,c = ix_(A,L,t,c)
G = A*exp(L+t)
G = G/G.max()
G1 = G*c

print(G1)
print(G1.shape)

v = gl.GLVolumeItem(G1)
#v.translate(-50,-50,-50)
w.addItem(v)

ax = gl.GLAxisItem()
w.addItem(ax)

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
