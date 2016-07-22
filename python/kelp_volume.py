# kelp_volume.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Mon 11 Jul 2016 07:55:11 PM EDT
# Last Edited: Mon 11 Jul 2016 07:57:29 PM EDT

from numpy import *
from vispy_volume import volume
import pickle
with open('kelp_pickle/xx_3d.pickle','rb') as inFile:
    xx = pickle.load(inFile)

with open('kelp_pickle/yy_3d.pickle','rb') as inFile:
    yy = pickle.load(inFile)

with open('kelp_pickle/zz_3d.pickle','rb') as inFile:
    zz = pickle.load(inFile)

with open('kelp_pickle/PP_3d.pickle','rb') as inFile:
    PP = pickle.load(inFile)

PP1 = flipud(rollaxis(PP,2))
volume(xx,yy,zz,PP1)


