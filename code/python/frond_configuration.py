# frond_configuration.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Fri 17 Jun 2016 03:50:07 PM CEST
# Last Edited: Fri 17 Jun 2016 03:55:02 PM CEST

# Come up with a simple function which can approximate frond curvature
# as a function of base and tip locations.

from numpy import *
from matplotlib.pyplot import *

## Kelp parameters ##

# Frond length
LL = 10

x = linspace(0,LL,101)
y = cosh(x)

plot(x,y)
ion()
show()

