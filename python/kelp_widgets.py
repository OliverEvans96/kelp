import numpy as np
import bqplot as bq
import bqplot.interacts as bqi
import ipywidgets as ipw
import traitlets as tr
# # Widgets

class HandDrawFigure(bq.Figure):
    def __init__(self, traitful, trait_name, xvals=None, ylim=None, labels=None, color='blue'):
        """
        traitful: inherits from traitlets.HasTraits
        trait_name: string - name of trait of traitful, 1d numpy array of y values
        xvals: x coordinates of trait values, 1d numpy array
        ylim: length-2 list or tuple of y bounds: [ymin, ymax]
        labels: dict with keys: title, xlabel, ylabel
        """
        
        self.traitful = traitful
        yvals = getattr(traitful, trait_name)
        
        if xvals is None:
            xvals = np.arange(len(getattr(traitful, trait_name)), dtype=float)
        if ylim is None:
            ylim = (yvals.min(), yvals.max())
        if labels is None:
            labels={}
        if 'ylabel' not in labels.keys():
            labels['ylabel'] = trait_name
        if 'xlabel' not in labels.keys():
            labels['xlabel'] = 'x'
        if 'title' not in labels.keys():
            labels['title'] = ''
        
        self.title = labels['title']

        xscale = bq.LinearScale(min=xvals.min(), max=xvals.max())
        yscale = bq.LinearScale(min=ylim[0], max=ylim[1])
        xax = bq.Axis(scale=xscale, label=labels['xlabel'], grid_lines='none')
        yax = bq.Axis(scale=yscale, label=labels['ylabel'], orientation='vertical', grid_lines='none')

        line = bq.Lines(x=xvals, y=yvals, scales={'x': xscale, 'y': yscale}, 
                                colors=[color], interpolation='cardinal')
        handdraw = bqi.HandDraw(lines=line)
        
        def update_vals(change):
            with out_area:
                values[quant] = change['new']

        link = tr.link((line, 'y'), (traitful, trait_name))
        
        super().__init__(marks=[line], axes=[xax, yax], interaction=handdraw)

# Define variables over depth
class RopeWidget(ipw.VBox):

    rope = tr.Any()

    def __init__(self, rope):

        super().__init__()

        title = ipw.HTML("<h3>Rope Parameters</h3>")

        self.rope = rope
        grid = rope.grid

        z_quants = ['water_speeds', 'water_angles', 'frond_lengths', 'frond_stds']
        z_scale = bq.LinearScale(min=grid.zmin, max=grid.zmax)
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


class IOPWidget(ipw.VBox):
    iops = tr.Any()

    def __init__(self):
        aw_slider = ipw.FloatSlider(min=0,max=10, description='$a_w$')
        bw_slider = ipw.FloatSlider(min=0,max=10, description='$a_w$')
        ak_slider = ipw.FloatSlider(min=0,max=10, description='$a_k$')
        bk_slider = ipw.FloatSlider(min=0,max=10, description='$b_k$')

        self.children = [
            aw_slider,
            bw_slider,
            ak_slider,
            bk_slider
        ]

        links = [
            tr.link((aw_slider, 'value'), (iops, 'a_water')),
            tr.link((bw_slider, 'value'), (iops, 'a_water')),
            tr.link((ak_slider, 'value'), (iops, 'a_kelp')),
            tr.link((bk_slider, 'value'), (iops, 'b_kelp'))
        ]


