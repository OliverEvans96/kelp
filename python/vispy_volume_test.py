#vispy_volume.py
#Oliver Evans
#4-21-2016

#Plot function of 2 variables in 3D using vispy

import numpy as np
from vispy import app,scene,visuals,gloo
from vispy.color import BaseColormap,colormap

# Ignore IPython ShimWarning
import warnings
from IPython.utils.shimmodule import ShimWarning
warnings.filterwarnings("ignore",category=ShimWarning)

# create colormaps that work well for translucent and additive volume rendering
class TransFire(BaseColormap):
    glsl_map = """
    vec4 translucent_fire(float t) {
        return vec4(pow(t, 0.5), t, t*t, max(0, t*1.05 - 0.05));
    }
    """


class TransGrays(BaseColormap):
    glsl_map = """
    vec4 translucent_grays(float t) {
        return vec4(t, t, t, t*0.05);
    }
    """

# xlim,ylim,zlim are axes limits (2-vector: [xmin,xmax])
# vol_array is a 3d array of values
class Canvas(scene.SceneCanvas):
    def __init__(self,xlim,ylim,zlim,vol_array,color=None,*args,**kwargs):
        scene.SceneCanvas.__init__(self,'3D Volume Plot',keys='interactive',bgcolor=(0.9,0.9,0.9),*args,**kwargs)

        #Allow addition of new attributes
        self.unfreeze()

        #Initialize bounds
        self._xlim=xlim
        self._ylim=ylim
        self._zlim=zlim

        #Initialize axes positions
        self._axpos=np.empty([3,3],int)
        self._last_axpos=np.empty([3,3],int)
        self._axpos_changed=True

        #Initialize plane position
        self._plane_pos=[0,0,0]

        #Keywords to pass to plot
        self._plot_kwargs={}
        if color is not None:
            self._plot_kwargs['color']=color

        #Calculate everything
        self.set_values(vol_array)
        self.calc_bounds()
        self.setup_camera()
        self.draw_axes()
        self.draw_planes()
        self.redraw_planes()

    def set_values(self,vol_array):

        #Calculate center and range
        self.calc_bounds()

        # Flip axes on volume array
        vol_flipped = np.flipud(np.swapaxes(vol_array,0,2))

        # Draw volumetric plot
        self._volume=scene.visuals.Volume(vol_flipped)
        self._nx,self._ny,self._nz = vol_array.shape
        self._volume.transform = scene.STTransform(
                #translate=(-self._nx/2,-self._ny/2,-self._nz/2),
                translate=((1-self._nx)/2,(1-self._nz)/2,0),
                scale=(self._xr/self._nx,self._yr/self._ny,self._zr/self._nz))
        self._volume.cmap = colormap.CubeHelixColormap()
        self._volume.method = 'mip'
        #self._volume.visible = False

    #Calculate bounds info
    def calc_bounds(self):
        #Center
        self._xc=np.mean(self._xlim)
        self._yc=np.mean(self._ylim)
        self._zc=np.mean(self._zlim)

        #Range
        self._xr=self._xlim[1]-self._xlim[0]
        self._yr=self._ylim[1]-self._ylim[0]
        self._zr=self._zlim[1]-self._zlim[0]

    def setup_camera(self):
        #Initial camera distance
        self._cam_dist=2*np.mean([self._xr,self._yr,self._zr])

        # Distance between axes and axes labels
        self._label_dist = self._cam_dist * 0.15

        #Create scene and camera
        self._view = self.central_widget.add_view()
        self._view.camera = scene.TurntableCamera(up='z',
                fov=60,azimuth=45,elevation=30,
                center=(self._xc,self._yc,self._zc),distance=self._cam_dist)

    # Draw axes
    def draw_axes(self):
        #First, calculate positions
        self.position_axes()

        #Draw x axis
        self._xax=scene.Axis(
            pos=[[self._xlim[0],
                self._ylim[self._axpos[0,1]],
                self._zlim[self._axpos[0,2]]],
                [self._xlim[1],
                self._ylim[self._axpos[0,1]],
                self._zlim[self._axpos[0,2]]]],
            domain=self._xlim,
            tick_direction=(0,-1,0),
            major_tick_length=10*self._cam_dist,
            font_size=10*self._cam_dist,
            font_face='FreeSerif',
            axis_color='k',
            tick_color='k',
            text_color='k',
            parent=self._view.scene)

        self._xlabel = scene.Text("x",
            pos=[(self._xlim[0]+self._xlim[1])/2,
                self._ylim[self._axpos[0,1]]-self._label_dist,
                self._zlim[self._axpos[0,2]]],
            font_size=10*self._cam_dist,
            face='FreeSerif',
            color='k',
            parent=self._view.scene)
                
        #Draw y axis
        self._yax=scene.Axis(
            pos=[[self._xlim[self._axpos[1,0]],
                self._ylim[0],
                self._zlim[self._axpos[1,2]]],
                [self._xlim[self._axpos[1,0]],
                self._ylim[1],
                self._zlim[self._axpos[1,2]]]],
            domain=self._ylim,
            tick_direction=(-1,0,0),
            major_tick_length=10*self._cam_dist,
            font_size=10*self._cam_dist,
            font_face='FreeSerif',
            axis_color='k',
            tick_color='k',
            text_color='k',
            parent=self._view.scene)

        self._ylabel = scene.Text("y",
            pos=[self._xlim[self._axpos[1,0]]-self._label_dist,
                (self._ylim[0]+self._ylim[1])/2,
                self._zlim[self._axpos[1,2]]],
            font_size=10*self._cam_dist,
            face='FreeSerif',
            color='k',
            parent=self._view.scene)

        #Draw z axis
        self._zax=scene.Axis(
            pos=[[self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]],
                self._zlim[0]],
                [self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]],
                self._zlim[1]]],
            domain=self._zlim,
            tick_direction=(0,-1,0),
            major_tick_length=10*self._cam_dist,
            font_face='FreeSerif',
            font_size=10*self._cam_dist,
            axis_color='k',
            tick_color='k',
            text_color='k',
            parent=self._view.scene)

        self._zlabel = scene.Text("z",
            pos=[self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]]-self._label_dist,
                (self._zlim[0]+self._zlim[1])/2],
            font_size=10*self._cam_dist,
            face='FreeSerif',
            color='k',
            parent=self._view.scene)

        self._xax._update_subvisuals()
        self._yax._update_subvisuals()
        self._zax._update_subvisuals()

        #Add plot to scene
        self._view.add(self._volume)

    def draw_planes(self):
        color=(245/255,245/255,245/255)

        xticks=np.unique(np.r_[self._xax._ticks.pos[:,0],self._xlim])
        yticks=np.unique(np.r_[self._yax._ticks.pos[:,1],self._ylim])
        zticks=np.unique(np.r_[self._zax._ticks.pos[:,2],self._zlim])

        self._px=scene.Plane(width=self._yr,height=self._zr,
                direction='+x',
                width_segments=yticks,height_segments=-zticks,
                color=None,edge_color='k')
        self._py=scene.Plane(width=self._zr,height=self._xr,direction='+y',
                width_segments=zticks,height_segments=-xticks,
                color=None,edge_color='k')
        self._pz=scene.Plane(width=self._xr,height=self._yr,direction='+z',
                width_segments=xticks,height_segments=-yticks,
                color=None,edge_color='k')

        self._view.add(self._px)
        self._view.add(self._py)
        self._view.add(self._pz)


    # Change position of axes without recreating them
    def redraw_axes(self):
        # Change axes
        self._xax.pos=[[self._xlim[0],
                self._ylim[self._axpos[0,1]],
                self._zlim[self._axpos[0,2]]],
                [self._xlim[1],
                self._ylim[self._axpos[0,1]],
                self._zlim[self._axpos[0,2]]]]

        self._yax.pos=[[self._xlim[self._axpos[1,0]],
                self._ylim[0],
                self._zlim[self._axpos[1,2]]],
                [self._xlim[self._axpos[1,0]],
                self._ylim[1],
                self._zlim[self._axpos[1,2]]]]

        self._zax.pos=[[self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]],
                self._zlim[0]],
                [self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]],
                self._zlim[1]]]

        # Change axes labels
        self._xlabel.pos = [(self._xlim[0]+self._xlim[1])/2,
                self._ylim[self._axpos[0,1]] + self._label_dist * (2*self._axpos[0,1]-1),
                self._zlim[self._axpos[0,2]]]

        self._ylabel.pos = [self._xlim[self._axpos[1,0]] + self._label_dist * (2*self._axpos[1,0]-1),
                (self._ylim[0]+self._ylim[1])/2,
                self._zlim[self._axpos[1,2]]]

        self._zlabel.pos = [self._xlim[self._axpos[2,0]],
                self._ylim[self._axpos[2,1]] + self._label_dist * (2*self._axpos[2,1]-1),
                (self._zlim[0]+self._zlim[1])/2]

        # Change tick directions accordingly
        self._xax.tick_direction = np.array([0, 2 * self._axpos[0,1] - 1, 0],dtype='float64')
        self._yax.tick_direction = np.array([2 * self._axpos[1,0] - 1, 0, 0],dtype='float64')
        self._zax.tick_direction = np.array([0, 2 * self._axpos[0,1] - 1, 0],dtype='float64')

    # Determine position of axes based on camera position to keep them out of the way
    # Axes positions - one row for each axis
    # One column per direction (x,y,z)
    # 0 means lower edge, 1 means upper edge
    # Recalculated every draw based on view position (like matplotlib's Axes 3D)
    # Main diagonal is meaningless
    def position_axes(self):
        azim=self._view.camera.get_state()['azimuth']
        elev=self._view.camera.get_state()['elevation']
        
        if elev >= 0:
            self._axpos[:,2]=0
        else:
            self._axpos[:,2]=1
        
        if 0 < azim:
            self._axpos[:,0]=[0,1,0]
        else:
            self._axpos[:,0]=[0,0,1]

        if -90<=azim<90:
            self._axpos[:,1]=0
        else:
            self._axpos[:,1]=1

        # Reset diagonal to 0
        np.fill_diagonal(self._axpos,0)

        if (self._axpos != self._last_axpos).any():
            self._axpos_changed=True

        self._last_axpos=np.copy(self._axpos)

    # Position of planes
    def position_planes(self):
        # x plane should be opposite y axis
        self._plane_pos[0]=1-self._axpos[1,0]

        # y plane should be opposite x axis
        self._plane_pos[1]=1-self._axpos[0,1]

        # z plane should be the same as the x and y axes
        self._plane_pos[2]=self._axpos[0,2] #(or self._axpos[1,2])

    def redraw_planes(self):
        self._px.transform = visuals.transforms.STTransform(
                translate=(
                    self._xlim[0]+self._plane_pos[0]*self._xr,
                    0,
                    0))

        self._py.transform = visuals.transforms.STTransform(
                translate=(
                    0,
                    self._ylim[0]+self._plane_pos[1]*self._yr,
                    0))

        self._pz.transform = visuals.transforms.STTransform(
                translate=(
                    0,
                    0,
                    self._zlim[0]+self._plane_pos[2]*self._zr))

    def on_draw(self,event):
        #self.context.set_clear_color((0, 0, 0, 1))
        #self.context.set_viewport(0, 0, *self.physical_size)
        #self.context.clear()

        self.position_axes()
        if(self._axpos_changed):
            self.position_planes()
            self.redraw_axes()
            self.redraw_planes()
            self._axpos_changed=False
        self._draw_scene()

    #Run!
    def run(self):
        self.show()
        with warnings.catch_warnings():
            app.run()

def volume(xlim,ylim,zlim,vol_array,**kwargs):
    c=Canvas(xlim,ylim,zlim,vol_array,**kwargs)
    c.run()
