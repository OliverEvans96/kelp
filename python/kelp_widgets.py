import numpy as np
import bqplot as bq
import bqplot.interacts as bqi
import ipywidgets as ipw
import ipyvolume as ipv
import traitlets as tr

class HandDrawFigure(bq.Figure):
    def __init__(self, traitful, trait_name, xdim=None, ylim=None, labels=None, color='blue'):
        """
        traitful: inherits from traitlets.HasTraits
        trait_name: string - name of trait of traitful, 1d numpy array of y values
        xdim: Traitful dimension with traits:
            'vals' as x coordinates for plot
            'minval' as minimum x coordinate
            'maxval' as maximum x coordinate
        ylim: length-2 list or tuple of y bounds: [ymin, ymax]
        labels: dict with keys: title, xlabel, ylabel
        """

        self.traitful = traitful
        self.trait_name = trait_name
        self.xdim = xdim
        self.ylim = ylim
        self.labels = labels
        self.color = color

        self.init_vals()
        self.init_elements()
        self.init_logic()
        self.init_style()

    def init_vals(self):
        self.yvals = getattr(
            self.traitful,
            self.trait_name
        )

        if self.xdim is None:
            self.xvals = None
        else:
            self.xvals = self.xdim.vals

        if self.xvals is None:
            self.xvals = np.arange(
                len(self.yvals),
                dtype=float
            )
        if self.ylim is None:
            self.ylim = (np.min(self.yvals), np.max(self.yvals))
        if self.labels is None:
            self.labels={}
        if 'ylabel' not in self.labels.keys():
            self.labels['ylabel'] = self.trait_name
        if 'xlabel' not in self.labels.keys():
            self.labels['xlabel'] = 'x'
        if 'title' not in self.labels.keys():
            self.labels['title'] = ''

        self.title = self.labels['title']

    def init_elements(self):
        self.xscale = bq.LinearScale(
            min=np.min(self.xvals),
            max=np.max(self.xvals)
        )
        self.yscale = bq.LinearScale(
            min=self.ylim[0],
            max=self.ylim[1]
        )
        self.xax = bq.Axis(
            scale=self.xscale,
            label=self.labels['xlabel'],
            grid_lines='none'
        )
        self.yax = bq.Axis(
            scale=self.yscale,
            label=self.labels['ylabel'],
            orientation='vertical',
            grid_lines='none'
        )
        self.line = bq.Lines(
            x=self.xvals,
            y=self.yvals,
            scales={'x': self.xscale, 'y': self.yscale},
            colors=[self.color],
            interpolation='cardinal'
        )
        self.handdraw = bqi.HandDraw(lines=self.line)

        super().__init__(
            marks=[self.line],
            axes=[self.xax, self.yax],
            interaction=self.handdraw
        )

    def init_logic(self):
        tr.link(
            (self.line, 'y'),
            (self.traitful, self.trait_name)
        )
        if self.xdim is not None:
            tr.link(
                (self.xscale, 'min'),
                (self.xdim, 'minval')
            )
            tr.link(
                (self.xscale, 'max'),
                (self.xdim, 'maxval')
            )
            tr.link(
                (self.line, 'x'),
                (self.xdim, 'vals')
            )

    def init_style(self):
        self.layout = ipw.Layout(height=u'300px')

class GridWidget(ipw.VBox):
    grid = tr.Any()

    def __init__(self, grid):
        super().__init__()
        self.grid = grid
        self.init_elements()
        self.init_style()
        self.init_layout()
        self.init_values()
        self.init_logic()


    def init_elements(self):

        self.title = ipw.HTML('<h3>Space/Angle Grid</h3>')

        self.xwidth_slider = ipw.FloatSlider(
            min=1,
			max=10,
			description='x width',
            continuous_update=False
        )
        self.ywidth_slider = ipw.FloatSlider(
            min=1,
			max=10,
			description='y width',
            continuous_update=False
		)
        self.zdepth_slider = ipw.FloatSlider(
            min=1,
			max=10,
			description='z depth',
            continuous_update=False
		)
        self.nx_spinner = ipw.BoundedIntText(
            min=1,
			max=1000,
			description='nx',
		)
        self.ny_spinner = ipw.BoundedIntText(
            min=1,
			max=1000,
			description='ny',
		)
        self.nz_spinner = ipw.BoundedIntText(
            min=1,
			max=1000,
			description='nz',
		)
        self.ntheta_spinner = ipw.BoundedIntText(
            min=1,
			max=1000,
			description='ntheta',
		)
        self.nphi_spinner = ipw.BoundedIntText(
            min=1,
			max=1000,
			description='nphi',
		)

    def init_layout(self):
        self.children = [
            self.title,
            ipw.HBox([
                ipw.VBox([
                    self.nx_spinner,
                    self.ny_spinner,
                    self.nz_spinner,
                    self.ntheta_spinner,
                    self.nphi_spinner
                ]),
                ipw.VBox([
                    self.xwidth_slider,
                    self.ywidth_slider,
                    self.zdepth_slider,
                ])
            ])
        ]


    def init_values(self):
        self.nx_spinner.value = 10
        self.ny_spinner.value = 10
        self.nz_spinner.value = 10
        self.ntheta_spinner.value = 10
        self.nphi_spinner.value = 10
        self.zdepth_slider.value = 10

    def init_style(self):
        self.nx_spinner.layout.width='150px'
        self.ny_spinner.layout.width='150px'
        self.nz_spinner.layout.width='150px'
        self.ntheta_spinner.layout.width='150px'
        self.nphi_spinner.layout.width='150px'

    def init_logic(self):
        tr.link(
            (self.nx_spinner, 'value'),
            (self.grid.x, 'num')),
        tr.link(
            (self.ny_spinner, 'value'),
            (self.grid.y, 'num')),
        tr.link(
            (self.nz_spinner, 'value'),
            (self.grid.z, 'num')),
        tr.link(
            (self.ntheta_spinner, 'value'),
            (self.grid.theta, 'num')),
        tr.link(
            (self.nphi_spinner, 'value'),
            (self.grid.phi, 'num')),
        tr.link(
            (self.zdepth_slider, 'value'),
            (self.grid.z, 'maxval'))

        self.xwidth_slider.dim = self.grid.x
        self.ywidth_slider.dim = self.grid.y
        self.zdepth_slider.dim = self.grid.z

        self.xwidth_slider.observe(self.set_width, names='value')
        self.ywidth_slider.observe(self.set_width, names='value')

    def set_width(self, change):
        space_dim = change['owner'].dim
        halfwidth = change['new'] / 2
        space_dim.minval = -halfwidth
        space_dim.maxval = halfwidth

# Define variables over depth
class RopeWidget(ipw.VBox):

    rope = tr.Any()

    def __init__(self, rope):

        super().__init__()

        title = ipw.HTML("<h3>Rope Parameters</h3>")

        self.rope = rope
        grid = rope.grid

        z_quants = ['water_speeds', 'water_angles', 'frond_lengths', 'frond_stds']
        z_scale = bq.LinearScale(min=grid.z.minval, max=grid.z.maxval)
        z_ax = bq.Axis(scale=z_scale, label='Depth (z)', grid_lines='none')

        mins = {
            'water_speeds': 0,
            'water_angles': 0,
            'frond_lengths': 0,
            'frond_stds': 0
        }

        maxs = {
            'water_speeds': 10,
            'water_angles': 2*np.pi,
            'frond_lengths': 1,
            'frond_stds': 1
        }

        colors = {
            'water_speeds': 'red',
            'water_angles': 'green',
            'frond_lengths': 'blue',
            'frond_stds': 'yellow'
        }

        ylabels = {
            'water_speeds': 'Water current velocity',
            'water_angles': 'Water current angle',
            'frond_lengths': 'Frond length mean',
            'frond_stds': 'Frond length std. dev.'
        }


        values = {}
        figs = {}

        for quant in z_quants:
            values[quant] = getattr(rope, quant)
            ylim = np.array([mins[quant], maxs[quant]], dtype=float)
            labels = {'xlabel': 'z', 'ylabel': '', 'title': ylabels[quant]}
            figs[quant] = HandDrawFigure(rope, quant, grid.z, ylim, labels, color=colors[quant])

        self.children = [
            title,
            ipw.HBox([
                figs['water_speeds'],
                figs['water_angles'],
            ]),
            ipw.HBox([
                figs['frond_lengths'],
                figs['frond_stds']
            ])
        ]

class BCWidget(ipw.VBox):
    bc = tr.Any()

    def __init__(self, bc):
        super().__init__()
        self.bc = bc
        self.init_elements()
        self.init_layout()
        self.init_logic()

    def init_elements(self):
        self.title = ipw.HTML("<h3>Sunlight Boundary Condition</h3>")
        self.theta_s_slider = ipw.FloatSlider(
            min=0,
            max=2*np.pi,
            step=np.pi/16,
            value=self.bc.theta_s,
            description=r'$\theta_s$'
        )
        self.phi_s_slider = ipw.FloatSlider(
            min=0,
            max=np.pi/2,
            step=np.pi/64,
            value=self.bc.phi_s,
            description=r'$\phi_s$'
        )
        self.max_rad_slider = ipw.FloatSlider(
            min=0,
            max=10,
            value = self.bc.max_rad,
            description = r'$L_{max}$'
        )
        self.decay_slider = ipw.FloatSlider(
            min=0,
            max=10,
            value = self.bc.decay,
            description = 'Decay'
        )

    def init_layout(self):
        self.children = [
            self.title,
            self.theta_s_slider,
            self.phi_s_slider,
            self.max_rad_slider,
            self.decay_slider
        ]

    def init_logic(self):
        tr.link((self.theta_s_slider, 'value'), (self.bc, 'theta_s'))
        tr.link((self.phi_s_slider, 'value'), (self.bc, 'phi_s'))
        tr.link((self.max_rad_slider, 'value'), (self.bc, 'max_rad'))
        tr.link((self.decay_slider, 'value'), (self.bc, 'decay'))

class ParamWidget(ipw.VBox):
    params = tr.Any()

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.init_elements()
        self.init_layout()
        self.init_logic()

    def init_elements(self):
        self.title = ipw.HTML("<h3>Numerical Parameters</h3>")
        self.maxiter_inner_slider = ipw.IntSlider(
            min=0,
            max=100,
            value=self.params.maxiter_inner,
            description=r'maxiter_inner'
        )
        self.maxiter_outer_slider = ipw.IntSlider(
            min=0,
            max=100,
            value=self.params.maxiter_outer,
            description=r'maxiter_outer'
        )
        self.tol_abs_picker = ipw.BoundedFloatText(
            min=0,
            max=10,
            value = self.params.tol_abs,
            description = 'Tol_Abs'
        )
        self.tol_rel_picker = ipw.BoundedFloatText(
            min=0,
            max=10,
            value = self.params.tol_rel,
            description = 'Tol_Rel'
        )

    def init_layout(self):
        self.children = [
            self.title,
            self.maxiter_inner_slider,
            self.maxiter_outer_slider,
            self.tol_abs_picker,
            self.tol_rel_picker
        ]

    def init_logic(self):
        tr.link((self.maxiter_inner_slider, 'value'), (self.params, 'maxiter_inner'))
        tr.link((self.maxiter_outer_slider, 'value'), (self.params, 'maxiter_outer'))
        tr.link((self.tol_abs_picker, 'value'), (self.params, 'tol_abs'))
        tr.link((self.tol_rel_picker, 'value'), (self.params, 'tol_rel'))


class IOPWidget(ipw.VBox):
    iops = tr.Any()
    vsf_vals = tr.Any()

    def __init__(self, iops):
        super().__init__()
        self.iops = iops
        self.init_vals()
        self.init_elements()
        self.init_layout()
        self.init_logic()

    def init_elements(self):
        self.title = ipw.HTML("<h3>Optical Properties</h3>")
        self.aw_slider = ipw.FloatSlider(min=0,max=10, description='$a_w$')
        self.bw_slider = ipw.FloatSlider(min=0,max=10, description='$b_w$')
        self.ak_slider = ipw.FloatSlider(min=0,max=10, description='$a_k$')
        self.bk_slider = ipw.FloatSlider(min=0,max=10, description='$b_k$')

        self.vsf_plot = HandDrawFigure(
            self, 'vsf_vals',
            xdim=self.iops.grid.theta,
            labels={
                'title': 'Volume Scattering Function',
                'xlabel': 'theta',
                'ylabel': 'VSF'
            },
            #ylim=[0,1]
        )

    def init_logic(self):
        tr.link((self.iops, 'a_water'), (self.aw_slider, 'value'))
        tr.link((self.iops, 'b_water'), (self.bw_slider, 'value'))
        tr.link((self.iops, 'a_kelp'), (self.ak_slider, 'value'))
        tr.link((self.iops, 'b_kelp'), (self.bk_slider, 'value'))
        self.observe(self.set_iops_vsf, names='vsf_vals')
        #self.iops.observe(self.set_widget_vsf, names='vsf')

    def init_vals(self):
        # This should be a copy, just to set the inital value
        self.vsf_vals = self.iops.vsf_vals

    def set_iops_vsf(self, *args):
        self.iops.set_vsf(self.vsf_vals)

    def set_widget_vsf(self, *args):
        self.vsf_vals = self.iops.vsf

    def init_layout(self):
        self.children = [
            self.title,
            ipw.HBox([
                ipw.VBox([
                    self.aw_slider,
                    self.bw_slider,
                    self.ak_slider,
                    self.bk_slider
                ]),
                self.vsf_plot
            ])
        ]

class FrondWidget(ipw.VBox):
    frond = tr.Any()

    def __init__(self, frond):
        super().__init__()
        self.frond = frond
        self.init_elements()
        self.init_vals()
        self.init_layout()
        self.init_logic()

    def init_elements(self):
        self.title = ipw.HTML('<h3>Frond Parameters</h3>')
        self.fr_slider = ipw.FloatSlider(
            min=0,
            max=5,
            description='$f_r$'
        )
        self.fs_slider = ipw.FloatSlider(
            min=0,
            max=5,
            description='$f_s$'
        )

    def init_vals(self):
        self.fr_slider.value = self.frond.fr
        self.fs_slider.value = self.frond.fs

    def init_layout(self):
        self.children = [
            self.title,
            self.fr_slider,
            self.fs_slider
        ]

    def init_logic(self):
        tr.link(
            (self.fr_slider, 'value'),
            (self.frond, 'fr')
        )
        tr.link(
            (self.fs_slider, 'value'),
            (self.frond, 'fs')
        )


class VolumePlotWidget(ipw.VBox):
    kelp = tr.Any()
    kelp_figure = tr.Any()
    kelp_controls = tr.Any()
    light = tr.Any()
    light_figure = tr.Any()
    light_controls = tr.Any()

    def __init__(self, kelp, light):
        super().__init__()

        self.kelp = kelp
        self.light = light

        self.init_vals()
        self.init_elements()
        self.init_layout()
        self.init_logic()

    def init_vals(self):
        # if self.kelp.p_kelp is None:
        #     self.kelp.gen_kelp()

        # if self.light.irradiance is None:
        #     self.light.calculate_light_field()
        pass

    def init_elements(self):
        self.log_area = ipw.Output()
        self.title = ipw.HTML("<h3>Volume Plots</h3>")
        self.init_kelp_vis()
        self.init_light_vis()
        self.init_control_panel()

    def init_kelp_vis(self):
        self.kelp_header = ipw.HTML("<b>Kelp</b>")
        with self.log_area:
            self.kelp_controls, self.kelp_figure = ipv.quickvolshow(
                np.zeros([3, 3, 3])
            ).children

            self.kelp_figure.xlabel = 'x'
            self.kelp_figure.xlim = [
                self.kelp.grid.x.minval,
                self.kelp.grid.x.maxval
            ]

            self.kelp_figure.ylabel = 'z'
            self.kelp_figure.ylim = [
                self.kelp.grid.z.minval,
                self.kelp.grid.z.maxval
            ]

            self.kelp_figure.zlabel = 'y'
            self.kelp_figure.zlim = [
                self.kelp.grid.y.minval,
                self.kelp.grid.y.maxval
            ]


        self.calculate_kelp_button = ipw.Button(
            description="Calculate Kelp"
        )

    def init_light_vis(self):
        self.light_header = ipw.HTML("<b>Irradiance</b>")
        with self.log_area:
            self.light_controls, self.light_figure = ipv.quickvolshow(
                np.zeros([3, 3, 3])
            ).children

            self.light_figure.xlabel = 'x'
            self.light_figure.xlim = [
                self.light.grid.x.minval,
                self.light.grid.x.maxval
            ]

            # z-axis is flipped
            self.light_figure.ylabel = 'z'
            self.light_figure.ylim = [
                self.light.grid.z.maxval,
                self.light.grid.z.minval
            ]

            self.light_figure.zlabel = 'y'
            self.light_figure.zlim = [
                self.light.grid.y.minval,
                self.light.grid.y.maxval
            ]

        self.calculate_light_button = ipw.Button(
            description='Calculate Light'
        )

    def init_control_panel(self):
        self.control_panel = ipw.VBox([
            self.kelp_header,
            self.kelp_controls,
            self.light_header,
            self.light_controls
        ])

    def init_layout(self):
        self.kelp_vis = ipw.VBox([
            self.kelp_header,
            self.kelp_figure,
            self.calculate_kelp_button
        ])

        self.light_vis = ipw.VBox([
            self.light_header,
            self.light_figure,
            self.calculate_light_button
        ])

        self.children = [
            self.title,
            ipw.HBox([
                self.kelp_vis,
                self.light_vis
            ])
        ]

    def init_logic(self):
        # self.light.observe(
        #     self.load_irradiance,
        #     names='irradiance'
        # )
        # self.kelp.observe(
        #     self.load_kelp,
        #     names='p_kelp'
        # )

        self.calculate_kelp_button.on_click(
            self.calculate_kelp
        )
        self.calculate_light_button.on_click(
            self.calculate_light
        )

        self.link_figures(
            self.kelp_figure,
            self.light_figure
        )

    def link_figures(self, fig1, fig2):
        traits = [
            'camera_control',
            'camera_center',
            'anglex',
            'angley',
            'anglez',
            'camera_fov',
        ]
        return [
            self.link_trait(fig1, fig2, trait)
            for trait in traits
        ]

    def link_trait(self, t1, t2, trait_name):
        return tr.link(
            (t1, trait_name),
            (t2, trait_name)
        )

    def calculate_kelp(self, *args):
        self.kelp.gen_kelp()
        self.load_kelp()

    def calculate_light(self, *args):
        self.light.calculate_light_field()
        self.load_irradiance()

    def load_kelp(self, *args):
        self.update_vol_plot(
            self.kelp_figure,
            self.kelp.p_kelp
        )

    def load_irradiance(self, *args):
        self.update_vol_plot(
            self.light_figure,
            self.light.irradiance
        )

    def update_vol_plot(self, fig, data):
        with self.log_area:
            new = self.vol_transform(data)
            fig.data_min = np.min(new)
            fig.data_max = np.max(new)
            fig.volume_data = new


    def vol_transform(self, vol_data):
        "Transform 3D array for volume plot"
        # Move x axis to end
        return np.transpose(vol_data[:,:,::-1], (1, 2, 0))

