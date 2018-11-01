import os
import sqlite3
import netCDF4 as nc
import re
import subprocess
import time
import concurrent.futures as cf
from math import ceil

import kelp_compute
import run_utils as ru

def chunks(l, n):
    """Yield successive n-sized chunks from l.
    https://stackoverflow.com/a/312464/4228052
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]

def rename_table(conn, old_study_name, new_study_name):
    rename_query = 'ALTER TABLE {} RENAME TO {}'.format(
        old_study_name,
        new_study_name
    )
    print(
        "Executing '{}'".format(rename_query),
    )
    conn.execute(rename_query)
    conn.commit()

def rename_combined_db(study_dir, old_study_name, new_study_name):
    """
    Rename combined .db file (don't rename table or move to another directory).
    """
    old_combined_db_path = os.path.join(
        study_dir,
        '{}.db'.format(old_study_name)
    )

    new_combined_db_path = os.path.join(
        study_dir,
        '{}.db'.format(new_study_name)
    )

    os.rename(
        old_combined_db_path,
        new_combined_db_path
    )

def rewrite_data_path(conn, old_study_name, new_study_name):
    """
    Replace data_path in each row in table.
    /path/to/`old_study_name`/data/*.nc -> /path/to/`new_study_name`/data/*.nc

    ** Must be run BEFORE renaming table **
    (this function looks for a table called `old_study_name`)
    """

    select_query = 'SELECT id, data_path FROM {}'.format(old_study_name)
    update_query = 'UPDATE {} SET data_path = ? WHERE id = ?'.format(
        old_study_name
    )

    # Presumably there is only one row per db,
    # but loop just in case
    for row_id, old_data_path in conn.execute(select_query).fetchall():
        pattern = r'{}/data/(.*.nc)'.format(old_study_name)
        repl = r'{}/data/\1'.format(new_study_name)
        new_data_path = re.sub(pattern, repl, old_data_path)

        print("DP {} -> {}".format(old_data_path, new_data_path))

        # Write correct value to .db
        # print(
        #     "Executing '{}' with".format(update_query),
        #     (new_data_path, row_id)
        # )
        conn.execute(update_query, (new_data_path, row_id))

    conn.commit()

def data_path_to_data_filename(conn):
    """
    Replace data_path in each row in table.
    """

    study_name = ru.get_table_names(conn)[0]

    select_query = 'SELECT id, data_path FROM {}'.format(study_name)
    update_query = 'UPDATE {} SET data_path = ? WHERE id = ?'.format(
        study_name
    )

    # Presumably there is only one row per db,
    # but loop just in case
    for row_id, old_data_path in conn.execute(select_query).fetchall():
        new_data_path = os.path.basename(old_data_path)

        # Write correct value to .db
        #print(
        #    "Executing '{}' with".format(update_query),
        #    (new_data_path, row_id)
        #)
        conn.execute(update_query, (new_data_path, row_id))

    conn.commit()

def retry_dptdf(db_path, retries=5):
    print()
    print(db_path)
    conn = sqlite3.connect(db_path)
    for i in range(retries):
        try:
            data_path_to_data_filename(conn)
            conn.close()
            break
        except IndexError:
            print("ERROR: no tables.")
            break
        except sqlite3.OperationalError:
            time.sleep(1)
            print("Op Err. Retry {}".format(i+1))

def multi_dptdf(db_path_list):
    for db_path in db_path_list:
        retry_dptdf(db_path)

if __name__ == '__main__':
    base_dir = os.path.join(
        os.environ['SCRATCH'],
        'kelp-results',
    )

    print("Looking.")
    nthreads = 1000
    db_paths = subprocess.check_output(r'find {} | grep \.db$'.format(base_dir), shell=True).decode().split('\n')
    chunk_size = ceil(len(db_paths)/nthreads)
    path_chunks = chunks(db_paths, chunk_size)
    print("Paths: {}".format(db_paths))

    fut_list = []
    with cf.ThreadPoolExecutor() as ex:
        for i, chunk in enumerate(path_chunks):
            #print("CHUNK {}: {}".format(i, chunk))
            fut = ex.submit(multi_dptdf, chunk)
            fut_list.append(fut)
    #print([f.exception() for f in fut_list[:1]])
