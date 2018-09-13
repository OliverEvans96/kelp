"""
Utilities for runs and studies.
Functions related to SQLite and NetCDF.
"""

import functools
import inspect
import itertools as it
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time

import ipyparallel as ipp
import netCDF4 as nc
import numpy as np

## Misc. utils ##

def get_random_unused_filename(dir='.', prefix='', suffix=''):
    with tempfile.NamedTemporaryFile(dir=dir, prefix=prefix, suffix=suffix, delete=False) as fh:
        filename = fh.name

    return filename

def get_git_commit_hash():
    """Get hash of current git branch."""
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()

def dict_from_args(func, args, kwargs={}):
    # Get function signature
    sig = inspect.signature(func)
    params = list(sig.parameters)

    arg_dict = dict(zip(params, args))
    return {**arg_dict, **kwargs}


## SQLite utils (non-project-specific) ##

def create_db_generic_table_from_tuple(conn, table_name, columns, prefix='.'):
    """
    Create sqlite db file, create table, and return connection.
    Assume table_name is unique.
    """

    create_table_template = ('''
    CREATE TABLE {table_name} (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        '''
    + (',\n'+8*' ').join([
        ' '.join(col)
        for col in columns
    ]) + '\n'+4*' '+');')

    create_table_command = create_table_template.format(table_name=table_name)

    print("CREATING TABLE")
    print(create_table_command)

    sys.stderr.write("Failed @ create_table with command: '{}'\n".format(create_table_command))

    conn.execute(create_table_command)
    conn.commit()

def create_db_table_from_dict(conn, table_name, data_dict, prefix='.'):
    """
    data_dict should be a `collections.OrderedDict`
    since the order of the columns is very important!
    """
    sql_type_dict = {
        float: 'REAL',
        int: 'INTEGER',
        str: 'CHAR(1024)'
    }

    # Convert numpy types to Python types
    for name, value in data_dict.items():
        if isinstance(value, np.number):
            data_dict[name] = np.asscalar(value)

    columns = (
        (name, sql_type_dict[type(value)])
        for name, value in data_dict.items()
    )

    return create_db_generic_table_from_tuple(conn, table_name, columns)

def update_col(conn, table_name, col_name, re_from, re_to):
    cur = conn.execute('select * from {}'.format(table_name))
    for row in cur:
        row_id = row[0]
        columns = [c[0] for c in cur.description]
        col_num = columns.index(col_name)
        col_in = row[col_num]
        col_out = re.sub(re_from, re_to, col_in)
        print("{}: '{}' -> '{}'".format(row_id, col_in, col_out))
        conn.execute(
            'UPDATE {} SET {}=? WHERE id=?'.format(table_name, col_name),
            (col_out, row_id)
        )
    conn.commit()
    print("Done!")

def get_table_names(conn):
    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
    # This is a list of 1-tuples which should be flattened.
    tables_preflatten = cursor.fetchall()
    tables = [table for tup in tables_preflatten for table in tup]

    return tables

def get_col_names(conn, table_name):
    cur = conn.execute('PRAGMA table_info({});'.format(table_name))
    # Each result from this query contains:
    # ('cid', 'name', 'type', 'notnull', 'dflt_value', 'pk')
    # We just want the names.
    names = [r[1] for r in cur.fetchall()]
    return names

def insert_generic_row(conn, data_dict, table_name):
    """Insert row, and create table if not present.
    NOTE: If the table does not yet exist, it's important
    that `data_dict` be a `collections.OrderedDict` so that
    the columns are created in the correct order.
    """

    # Create table if not present
    if table_name not in get_table_names(conn):
        print("Creating table {}".format(table_name))
        create_db_table_from_dict(conn, table_name, data_dict)
    else:
        print("Not creating table {}".format(table_name))

    print("col_names: {}".format(get_col_names(conn, table_name)))

    insert_template = (
        'INSERT INTO {table_name} VALUES ('
        + ', '.join([
            ':{}'.format(col_name)
            for col_name in get_col_names(conn, table_name)
        ])
        + ');'
    )

    insert_command = insert_template.format(table_name=table_name)

    print("Before try.")
    print(insert_command)

    try:
        print("Executing command:")
        print(insert_command)
        conn.execute(insert_command, {**data_dict, 'id': None})
    except sqlite3.OperationalError as e:
        print('Failed to insert row using command:')
        print("'{}'".format(insert_command))
        raise e


## NetCDF/SQLite (project-specific) ##

def create_dirs(study_dir):
    try:
        os.mkdir(study_dir)
    except FileExistsError:
        print("Appending existing study directory.")
    except OSError as err:
        raise err
    else:
        print("Creating new study directory.")

    os.makedirs(os.path.join(study_dir, 'data'), exist_ok=True)

def create_db_run_table(conn, table_name, prefix='.'):
    columns = (
        ('absorptance_kelp', 'REAL'),
        ('a_water', 'REAL'),
        ('b', 'REAL'),
        ('ns', 'REAL'),
        ('nz', 'REAL'),
        ('na', 'REAL'),
        ('num_dens', 'REAL'),
        ('kelp_dist', 'CHAR(32)'),
        ('fs', 'REAL'),
        ('fr', 'REAL'),
        ('ft', 'REAL'),
        ('max_length', 'REAL'),
        ('length_std', 'REAL'),
        ('zmax', 'REAL'),
        ('rope_spacing', 'REAL'),
        ('I0', 'REAL'),
        ('phi_s', 'REAL'),
        ('theta_s', 'REAL'),
        ('decay', 'REAL'),
        ('num_cores', 'INTEGER'),
        ('num_scatters', 'INTEGER'),
        ('fd_flag', 'INTEGER'),
        ('lis_opts', 'CHAR(256)'),
        ('date', 'CHAR(64)'),
        ('git_commit', 'CHAR(40)'),
        ('compute_time', 'REAL'),
        ('lis_iter', 'INTEGER'),
        ('lis_time', 'REAL'),
        ('lis_resid', 'REAL'),
        # Array data (NETCDF)
        ('data_path', 'CHAR(256')
    )

    return create_db_generic_table_from_tuple(conn, table_name, columns, prefix='.')

def run_not_present(completed_run_list, run_func, run_args, run_kwargs):
    """
    Check whether run has been completed already by checking
    against completed_run_list generated by get_completed_run_list
    (instead of checking every .db file for every run)
    """
    # Combine args and kwargs
    run_dict = dict_from_args(run_func, run_args, run_kwargs)

    # Assume this is the first time
    # until proven otherwise
    run_present = False

    # Check each run already recorded
    # Should only be one per file, but just in case
    for row_dict in completed_run_list:
        # Assume row matches until
        # a difference is found
        row_found = True

        # Check each function parameter
        # (assume they have the same names as db columns),
        # only need to check params in signature since
        # therer are probably more db columns than parameters.
        for param, run_val in run_dict.items():
            # Not a match if db value doesn't equal run value
            if row_dict[param] != run_val:
                row_found = False
                break

        # Stop looking through rows if we found a match
        if row_found:
            run_present = True
            break

    return not run_present

def get_completed_run_list(study_dir):
    """
    List runs which are already completed for the sake of
    determining whether to perform a run (don't want to redo work.)
    """
    data_dir = os.path.join(study_dir, 'data')

    completed_run_list = []

    # Look in each db file in data dir
    for filename in [fn for fn in os.listdir(data_dir) if re.match(r'.*\.db', fn)]:
        db_path = os.path.join(data_dir, filename)

        # Open database
        conn = sqlite3.connect(db_path)
        # If db was created but table was not, this will fail
        try:
            table_name = get_table_names(conn)[0]
        except IndexError:
            # (in which case we should skip to the next file)
            continue
        # Same action if DB is corrupt/not fully written
        except sqlite3.DatabaseError:
            continue

        table_cursor = conn.execute('select * from {}'.format(table_name))
        cols = [d[0] for d in table_cursor.description]

        completed_run_list += [
            dict(zip(cols, row_list))
            for row_list in table_cursor
        ]

    return completed_run_list

def insert_run(conn, table_name=None, **data_dict):
    # If table_name is not provided,
    # just use the first one we find
    if not table_name:
        table_name = get_table_names(conn)[0]

    insert_generic_row(conn, data_dict, table_name)

def combine_dbs(study_dir, table_name):
    data_dir = os.path.join(study_dir, 'data')
    dbs = [
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir) if re.match(r'.*\.db$', f)
    ]

    combined_db = os.path.join(study_dir, '{}.db'.format(table_name))
    print("Opening combined db: {}".format(combined_db))
    combined_conn = sqlite3.connect(combined_db)

    print("Connected.")
    for db in dbs:
        # Read from individual tables
        conn = sqlite3.connect(db)
        print("Combining {} (tables: {})".format(db, get_table_names(conn)))
        cursor = conn.execute('SELECT * FROM {}'.format(table_name))
        print("read.")
        columns = [tup[0] for tup in cursor.description]
        # Don't want to duplicate id column, so exclude it here.
        id_col = columns.index('id')
        columns_without_id = [
            col
            for i, col in enumerate(columns)
            if i != id_col
        ]

        # Write to combined table
        for row_tuple in cursor:
            row_without_id = [
                entry
                for i, entry in enumerate(row_tuple)
                if i != id_col
            ]
            row_dict = dict(zip(columns_without_id, row_without_id))
            insert_run(combined_conn, table_name, **row_dict)
        conn.close()

    combined_conn.commit()
    combined_conn.close()

def create_nc(data_path, **results):
    """Create netCDF file to store results"""

    # Create file
    rootgrp = nc.Dataset(data_path, 'w', format='NETCDF4_CLASSIC')

    # Get dimension sizes
    nx = results['nx']
    ny = results['ny']
    nz = results['nz']
    ntheta = results['ntheta']
    nphi = results['nphi']
    nomega = results['nomega']

    # Create Dimensions
    x = rootgrp.createDimension('x', nx)
    y = rootgrp.createDimension('y', ny)
    z = rootgrp.createDimension('z', nz)
    theta = rootgrp.createDimension('theta', ntheta)
    phi = rootgrp.createDimension('phi', nphi)
    omega = rootgrp.createDimension('omega', nomega)

    # Create dimension variables
    # xvals = rootgrp.createVariable('x','f8', ('x',))
    # yvals = rootgrp.createVariable('y','f8', ('y',))
    # zvals = rootgrp.createVariable('z','f8', ('z',))
    # thetavals = rootgrp.createVariable('theta','f8', ('theta',))
    # phivals = rootgrp.createVariable('phi','f8', ('phi',))

    # Assign dimension variables
    # xvals[:] = results.pop('x')
    # yvals[:] = results.pop('y')
    # zvals[:] = results.pop('z')
    # thetavals[:] = results.pop('theta')
    # phivals[:] = results.pop('phi')

    # Decide which dimensions to used based on the
    # number of dimensions in the array
    ndim_to_dim_dict = {
        1: ('z', ),
        3: ('x', 'y', 'z'),
        4: ('x', 'y', 'z', 'omega'),
    }

    # Datatype to use for array variables
    array_dtype = 'f8'

    # Loop over results
    for var_name, val in results.items():

        # Store array data
        ndim = len(np.shape(val))
        if ndim > 0:
            nc_arr_var = rootgrp.createVariable(
                var_name,
                array_dtype,
                ndim_to_dim_dict[ndim]
            )

            nc_arr_var[:] = results[var_name]

        # Scalar metadata
        else:
            # Convert numpy types to Python types
            if isinstance(val, np.number):
                val = np.asscalar(val)

            var_type = type(val)
            print("{}: {} ({})".format(var_name, val, var_type))
            # String handling for interoperability with fortran:
            # https://stackoverflow.com/a/37091930
            if var_type == str:
                # Strings treated as character arrays
                dim_name = '{}_dim'.format(var_name)
                char_dim = rootgrp.createDimension(dim_name, len(val))

                type_str = 'S{}'.format(len(val))
                char_val = nc.stringtochar(np.array([val], type_str))
                print("CREATE STR VAR: ('{}', '{}', '{}')".format(var_name, type_str, dim_name))
                var = rootgrp.createVariable(var_name, 'S1', (dim_name,))
                var[...] = char_val
            else:
                # Scalar int or float
                if var_type == float:
                    type_str = 'f4'
                elif var_type == int:
                    type_str = 'i4'

                var = rootgrp.createVariable(var_name, type_str)
                var[...] = val

    # Close file
    rootgrp.close()

    return data_path

def print_call(run_func, run_args, run_kwargs):
    func_name = run_func.__name__
    args_str = ', '.join([arg.__repr__() for arg in run_args])
    kwargs_str = ', '.join(['{}={}'.format(k, v.__repr__()) for k, v in run_kwargs.items()])
    sig_str = ', '.join([args_str, kwargs_str])
    print("Calling {}({})".format(func_name, sig_str))


## Decorators ##

def study_decorator(study_func):
    """Create directories before execution and
    merge dbs afterwards.

    Should be applied to functions like grid_study, etc.

    Function should return func_list, args_list, kwargs_list.
    Each of these will be applied via executor.apply().

    executor should be something like `ipp.Client().load_balanced_view()`
    with an `apply` method.

    Store results in sqlite database
    `study_name` will be database name
    Database file located at {base_dir}/{study_name}/{study_name.db}
    Other data files located at {base_dir}/{study_name}/data/*.nc
    """
    @functools.wraps(study_func)
    def wrapper(study_name, *study_args, base_dir=os.curdir, executor=None, **study_kwargs):

        study_dir = os.path.join(base_dir, study_name)
        if not executor:
            executor = ipp.Client().load_balanced_view()

        create_dirs(study_dir)
        study_calls = study_func(*study_args, **study_kwargs)

        run_futures = []

        # Read all .dbs first, then check each run against the list
        # (keeping all in memory is way better than reading every .db
        # every time, especially when there are lots of .dbs)
        print("Reading existing runs.")
        completed_run_list = get_completed_run_list(study_dir)
        print("Finished reading existing runs.")

        # Execute function calls from study function
        for run_func, run_args, run_kwargs in it.zip_longest(*study_calls):
            # Make args, kwargs optional
            if not run_args:
                run_args = ()
            if not run_kwargs:
                run_kwargs = {}

            # Add these keyword arguments for `run_wrapper`
            # Don't overwrite run_kwargs so as not to
            # confuse `run_not_present` with kwargs which
            # are not db columns
            appended_run_kwargs = {
                'study_dir': study_dir,
                'study_name': study_name,
                **run_kwargs
            }

            # called with these args
            if run_not_present(completed_run_list, run_func, run_args, run_kwargs):
                # Print function call before executing
                print_call(run_func, run_args, appended_run_kwargs)

                # Execute the function (pass to executor)
                run_futures.append(executor.apply(run_func, *run_args, **appended_run_kwargs))
            else:
                print("NOT ", end='')
                print_call(run_func, run_args, appended_run_kwargs)

            print()

        # TODO: I don't think this acually works

        # Once all functions have run, combine the results
        def wait_and_combine():
            for future in run_futures:
                future.wait()
                #print("{} futures done.".format(sum([f.done() for f in run_futures])))
            time.sleep(2)
            combine_dbs(study_dir, study_name)

        combine_thread = threading.Thread(target=wait_and_combine)
        combine_thread.start()

        # This will be returned immediately.
        # It can be probed with combined_thread.isAlive()
        return combine_thread, run_futures

    return wrapper

def run_decorator(run_func):
    """Run function and save results to database.
    The function is expected to return a dictionary
    of all values to be stored in the database.

    The .nc file will be written to data/`filename`.nc,
    and the temporary .db file will be written to `filename`.db,
    where `filename` is randomly generated.
    All .db files should later be merged into one db per directory.
    """
    @functools.wraps(run_func)
    def wrapper(*args, study_dir=None, study_name=None, **kwargs):
        scalar_params, results = run_func(*args, **kwargs)

        if not study_dir:
            raise ValueError("kwarg `study_dir` required for functions wrapped by `run_decorator`")
        if not study_name:
            raise ValueError("kwarg `study_name` required for functions wrapped by `run_decorator`")

        # Generate random name in correct directory
        # TODO: I'm a bit confused about the scope here.
        # Doesn't work without `ru.`
        data_path = ru.get_random_unused_filename(
            dir=os.path.join(study_dir, 'data'),
            suffix='.nc'
        )
        # Change suffix from .nc to .db for db path.
        db_path = re.sub(r'\.nc$', '.db', data_path)

        # SQL has no bool type
        for var_name, val in scalar_params.items():
            if isinstance(val, bool):
                scalar_params[var_name] = int(val)

        # Create data file containing results
        nc_dict = {
            **scalar_params,
            **results
        }
        # TODO: Scope again
        nc_data_path = ru.create_nc(data_path, **nc_dict)

        # Save to DB
        db_dict = {
            'data_path': nc_data_path,
            **scalar_params
        }

        conn = sqlite3.connect(db_path)
        # TODO: Scope again again
        ru.insert_run(conn, study_name, **db_dict)
        conn.commit()
        conn.close()

    return wrapper
