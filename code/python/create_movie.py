# create_movie.py
# Oliver Evans
# Clarkson University REU 2016
# Created: Thu 28 Jul 2016 12:18:30 AM EDT
# Last Edited: Thu 28 Jul 2016 12:25:26 AM EDT

# Create movie from working_visualization.py of constant rotation about z-axis

from working_visualization import *
from vispy import io
from vispy_volume import Canvas

c = Canvas(xlim,ylim,zlim,PP_3d,clim=None)
c.run()

c._view.camera.set_state({'elevation':20})

for i in range(40):
    c._view.camera.set_state({'azimuth':i*9})
    img=c.render()
    io.write_png('volume_img/img/{:03d}.png'.format(i),img)

