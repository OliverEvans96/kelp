import numpy as np
import bqplot as bq
import bqplot.interacts as bqi
import ipywidgets as ipw
import traitlets as tr
# # Widgets


# Define variables over depth
class RopeWidget(ipw.VBox):

    rope = tr.Any()

    def __init__(self, rope):

        super().__init__()

        title = ipw.HTML("<h3>Rope Parameters</h3>")

        self.rope = rope
        grid = rope.grid

        z_quants = ['vw', 'theta_w_rad', 'L_mean', 'L_std']
        z_scale = bq.LinearScale(min=grid.zmin, max=grid.zmax)
        z_ax = bq.Axis(scale=z_scale, label='Depth (z)', grid_lines='none')

        mins = {
            'vw': 0,
            'theta_w_rad': 0,
            'L_mean': 0,
            'L_std': 0
        }

        maxs = {
            'vw': 10,
            'theta_w_rad': 2*np.pi,
            'L_mean': 1,
            'L_std': 1
        }

        colors = {
            'vw': 'red',
            'theta_w_rad': 'green',
            'L_mean': 'blue',
            'L_std': 'yellow'
        }

        labels = {
            'vw': 'Water current velocity',
            'theta_w_rad': 'Water current angle',
            'L_mean': 'Frond length mean',
            'L_std': 'Frond length std. dev.'
        }


        values = {}

        values['L_mean'] = rope.frond_lengths
        values['L_std'] = rope.frond_stds
        values['vw'] = rope.water_speeds
        values['theta_w_rad'] = rope.water_angles

        ys = {}
        lines = {}
        handdraws = {}
        yax = {}
        figs = {}

        out_area = ipw.Output()

        for quant in z_quants:
            ys[quant] = bq.LinearScale(min=mins[quant], max=maxs[quant])
            lines[quant] = bq.Lines(x=grid.z, y=values[quant], scales={'x': z_scale, 'y': ys[quant]}, 
                                    colors=[colors[quant]], interpolation='cardinal')
            handdraws[quant] = bqi.HandDraw(lines=lines[quant])
            yax[quant] = bq.Axis(scale=ys[quant], label=labels[quant], orientation='vertical', grid_lines='none')
            figs[quant] = bq.Figure(marks=[lines[quant]], axes=[z_ax, yax[quant]], interaction=handdraws[quant])

            # Update values on handdraw
            # Define the function like this with default argument so that `quant` takes its
            # value at definition time, not evaluation time
            def update_vals(change, quant=quant):
                with out_area:
                    print()
                    print(quant)
                    values[quant] = change['new']
                    print('Updated!')
            lines[quant].observe(update_vals, names='y')

        links = [
            tr.link((lines['L_mean'], 'y'), (rope, 'frond_lengths')),
            tr.link((lines['L_std'], 'y'), (rope, 'frond_stds')),
            tr.link((lines['vw'], 'y'), (rope, 'water_speeds')),
            tr.link((lines['theta_w_rad'], 'y'), (rope, 'water_angles'))
        ]



        self.children = [
            title,
            ipw.HBox([
                figs['vw'],
                figs['theta_w_rad'],
            ]),
            ipw.HBox([
                figs['L_mean'],
                figs['L_std']
            ])
        ]
