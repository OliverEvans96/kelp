import datetime
import warnings
import itertools as it
import time
import concurrent.futures as cf

import dill
from arctic.exceptions import NoDataFoundException
from arctic import Arctic
from pymongo import MongoClient

import numpy as np
import pandas as pd
from sympy import divisors

import ipyparallel as ipp
import ipywidgets as ipw
import ipyvolume as ipv
import qgrid


## View

def param_grid(**kwargs):
    """Generate Cartesian product of keyword arguments"""
    prod = it.product(kwargs.items())
    name_list, values_list = zip(*kwargs.items())
    value_combinations = list(it.product(*values_list))
    return name_list, value_combinations

def param_df(**kwargs):
    """Generate pandas DataFrame from Cartesian product of keyword arguments"""
    names, values = param_grid(**kwargs)
    return pd.DataFrame(values, columns=names)

def param_qgrid(qgrid_layout=None, **kwargs):
    """Generate Qgrid table from Cartesian product of keyword arguments"""
    if not qgrid_layout:
        qgrid_layout=ipw.Layout()
    return qgrid.QGridWidget(df=param_df(**kwargs), layout=qgrid_layout)

def dict_list_from_df(df):
    """Turn each row into a dictionary indexed by column names"""
    return [{col: val for col, val in zip(df.columns, df.loc[ind,:])} for ind in df.index]

class ParamSpanRemoteConfig(object):
    def __init__(self):
        #self.init_db()
        self.init_engines()

    # def init_db(self):
    #     # Collect database info
    #     db_host = 'localhost'
    #     db_port = 31017
    #     lib_name = 'kale-param-spans'
    #     self.store = Arctic('mongodb://{}:{}'.format(db_host, db_port))
    #     self.store.initialize_library(lib_name)
    #     self.library = self.store[lib_name]
#
    #     # Package database info to send to engines
    #     self.db_info = [
    #         db_host,
    #         db_port,
    #         lib_name
    #     ]

    def init_engines(self):
        # Connect to controller
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.ipp_client = ipp.Client()

        # Establish views
        self.dview = self.ipp_client.direct_view()
        self.lview = self.ipp_client.load_balanced_view()

        # Tell engines their ids
        for i, engine in enumerate(self.ipp_client):
            engine.push({'engine_id': i})

        # Distribute database information
        # self.dview.push({'db_info': self.db_info})

        # def load_library():
        #     global store, library
        #     from arctic import Arctic
        #     db_host, db_port, lib_name = db_info
        #     store = Arctic('mongodb://{}:{}'.format(db_host, db_port))
        #     library = store[lib_name]
        # self.dview.apply(load_library)

class ParamSpanWidget(ipw.VBox):
    def __init__(self, param_span_name, remote_config, compute_func, vis_func, output_layout=None, qgrid_layout=None):
        self.compute_func = compute_func
        self.name = param_span_name
        self.remote_config = remote_config
        self.vis_func = vis_func
        self.output_layout = output_layout
        self.qgrid_layout = qgrid_layout

        super().__init__()

        self.init_executor()
        self.load_remote_config()
        self.init_widgets()
        self.init_layout()

    def init_executor(self):
        self.executor = cf.ThreadPoolExecutor()

    def load_remote_config(self):
        self.ipp_client = self.remote_config.ipp_client
        self.dview = self.remote_config.dview
        self.lview = self.remote_config.lview
        #self.db_info = self.remote_config.db_info
        #self.store = self.remote_config.store
        #self.library = self.remote_config.library
        self.lview = self.remote_config.lview

    def init_widgets(self):
        if not self.output_layout:
            self.output_layout = ipw.Layout(height='500px', border='1px solid', overflow_x='scroll', overflow_y='scroll')
        if not self.qgrid_layout:
            qgrid_layout = ipw.Layout()

        self.output = ipw.Output(layout=self.output_layout)
        # param_table is empty until set_params is called
        self.param_table = param_qgrid(self.qgrid_layout, **{'':[]})

    def init_logic(self):
        self.param_table.observe(self.visualize_wrapper, names='_selected_rows')

    def init_layout(self):
        self.children = [
            self.output,
            self.param_table
        ]

    def set_params(self, **all_params):
        """Provide parameter set to search over
        all_params = {
            'param1': [val1, val2, ...],
            'param2': [val3, val4, ...],
            ...
        }
        """
        self.param_table = param_qgrid(self.qgrid_layout, **all_params)
        self.init_logic()
        self.init_layout()

    def submit_and_store(self, paramset_id, paramset):
        """Submit one job and store the results, indexed by id."""
        fut = self.lview.apply(self.compute_func, **paramset)
        self.compute_futures[paramset_id] = fut
        self.results[paramset_id] = fut.result()

    def submit_computations(self, *change):
        # def compute_wrapper(compute_func, name, paramset_id, params, library):
        #     """Perform computation and send results to MongoDB for one set of params"""
        #     results = compute_func(**params)
        #
        #
        #     Index collection by paramset_id and paramspan name
        #     record_label = '{}-{}'.format(name, paramset_id)
        #     library.write(
        #         record_label,
        #         results
        #     )

        # Loop over all sets of parameters
        # paramset_id is the row index,
        # paramset is the dictionary of params
        paramitems = self.param_table.df.T.to_dict().items()
        self.results = [None]*len(paramitems)
        self.compute_futures = [None]*len(paramitems)
        self.save_futures = [None]*len(paramitems)
        for paramset_id, paramset in paramitems:
            # Submit task to IPyParallel
            # fut =  self.lview.apply(compute_wrapper, self.compute_func, self.name, paramset_id, paramset, self.library)
            # One thread per task so that results can be saved asynchronously when they return
            # Keep the outer future to make sure all results have been saved
            self.save_futures[paramset_id] = self.executor.submit(self.submit_and_store, paramset_id, paramset)

    def visualize_wrapper(self, *change):
        """Call visualization function and capture output"""
        # Do nothing if selection is empty:
        # Empty list evaluates to False
        if not self.param_table.get_selected_rows():
            return

        # Get params from selected row (take first row if >1 selected)
        # The ordering is not necessarily the same as in widget,
        # so weird things might happen if multiple rows are selected.
        paramset_id = self.param_table.get_selected_df().index[0]
        paramset = self.param_table.df.loc[paramset_id, :]

        # Search collection by parameters
        #record_label = '{}-{}'.format(self.name, paramset_id)
        # Clear screen if the results of this computation
        # are not available. Assume empty dict by default
        compute_results = []
        # try:
        #     compute_results = self.library.read(record_label).data
        # except NoDataFoundException:
        #     self.output.clear_output()
        future = self.compute_futures[paramset_id]
        # future is initially None
        if future and future.done():
            compute_results = self.results[paramset_id]

            # Avoid using output context to ensure that
            # only this function's output is included.
            @self.output.capture(clear_output=True, wait=True)
            def wrapper():
                self.vis_func(**compute_results)
            wrapper()

        else:
            @self.output.capture(clear_output=True, wait=True)
            def wrapper():
                print("Task {} not done: {}".format(paramset_id, self.compute_futures[paramset_id]))
            wrapper()

    def get_entries(self):
        """Get all results stored in database from parameter spans with this name"""
        return [entry for entry in self.library.list_symbols() if entry[:len(self.name)+1] == self.name+'-']

#     def delete_results(self):
#         """This will delete all results from parameter spans with this name!"""
#         for entry in self.get_entries():
#             self.library.delete(entry)

# User functions

def exp_compute(N, mean, std, color):
    import numpy as np
    from datetime import datetime
    x = np.random.normal(loc=mean, scale=std, size=N)
    realmean = np.mean(x)
    realstd = np.std(x)
    return {
        'engine_id': engine_id,
        'date': datetime.now().ctime(),
        'x': x,
        'N': N,
        'realmean': realmean,
        'realstd': realstd,
        'color': color
    }

def exp_viz(engine_id, date, x, N, realmean, realstd, color):
    print("Computed on engine {} at {}".format(engine_id, date))
    plt.figure(figsize=[8,5])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sns.distplot(x, color=color)
    plt.title(r'$N={}$, $\mu={:.2f}$, $\sigma={:.2f}$'.format(N, realmean, realstd))
    plt.show()
    print("Data: {}".format(x))

def kelp_visualize(duration, date, radiance, irradiance, p_kelp):
    print("Computation finished {}".format(date))
    print("Took {:.2f} seconds".format(duration))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kelpfig = ipv.figure()
        ipv.volshow(p_kelp.T, controls=False)
        # set initial angle
        kelpfig.anglex = 0.85 * np.pi
        kelpfig.angley = -0.30 * np.pi
        kelpfig.anglez = -0.67 * np.pi

        irradfig = ipv.figure()
        ipv.volshow(irradiance.T, controls=False)
        # set initial angle
        irradfig.anglex = 0.85 * np.pi
        irradfig.angley = -0.30 * np.pi
        irradfig.anglez = -0.67 * np.pi


    rad_text = "Rad Stats<br>---------<br>min: {:.3e}<br>max: {:.3e}<br>mean: {:.3e}<br>std: {:.3e}".format(
        radiance.min(),
        radiance.max(),
        radiance.mean(),
        radiance.std()
    )

    irrad_text = "Irrad Stats<br>-----------<br>min: {:.3e}<br>max: {:.3e}<br>mean: {:.3e}<br>std: {:.3e}".format(
        irradiance.min(),
        irradiance.max(),
        irradiance.mean(),
        irradiance.std()
    )

    display(ipw.HBox([
        kelpfig,
        irradfig,
        ipw.Box(layout=ipw.Layout(width='10px')),
        ipw.HTML("<tt>{}<br><br>{}</tt>".format(rad_text, irrad_text))
    ]))


def block_mean(large_arr, small_shape):
    """Calculate an array of block means of `large_arr` which has the shape `small_shape`"""

    if all(n_large % n_small == 0 for n_small, n_large in zip(small_shape, large_arr.shape)):
        # Try to abstract over number of dimensions
        avg = np.zeros(small_shape)
        block_sizes = [n_large // n_small for n_small, n_large in zip(small_shape, large_arr.shape)]
        #print("block_sizes = {}".format(block_sizes))

        # Loop over all combinations of indices
        for inds in it.product(*(range(n) for n in small_shape)):
            #print("inds = {}".format(inds))
            startinds = [ind * block_size for ind, block_size in zip(inds, block_sizes)]
            stopinds = [(ind+1) * block_size for ind, block_size in zip(inds, block_sizes)]
            slices = tuple(slice(startind, stopind) for startind, stopind in zip(startinds, stopinds))
            #print("startinds = {}".format(startinds))
            #print("stopinds = {}".format(stopinds))
            avg[inds] = large_arr[slices].mean()

        return avg

    else:
        raise IndexError("`small_shape` must divide `large_arr.shape` elementwise.")
