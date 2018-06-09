# coding: utf-8

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
