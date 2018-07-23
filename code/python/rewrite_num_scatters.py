# Accidentaly wrote num_scatters to DB as np.int64 instead of int.
# Rewrite using values from .nc files.

import os
import sqlite3
import netCDF4 as nc
import kelp_compute

def rewrite_num_scatters(base_dir, study_name):
    study_dir = os.path.join(base_dir, study_name)
    data_dir = os.path.join(study_dir, 'data')

    db_path_list = [
        os.path.join(
            data_dir,
            basename
        )
        for basename in os.listdir(data_dir)
        if basename[-3:] == '.db'
    ]

    select_query_template = 'SELECT id, data_path FROM {}'
    update_query_template = 'UPDATE {} SET num_scatters = ? WHERE id = ?'
    select_query, update_query = map(
        lambda template: template.format(study_name),
        (select_query_template, update_query_template)
    )

    for db_path in db_path_list:
        print(db_path)
        conn = sqlite3.connect(db_path)
        # Presumably there is only one row per db,
        # but loop just in case
        for row_id, data_path in conn.execute(select_query).fetchall():
            # Get correct value from .nc
            ds = nc.Dataset(data_path)
            num_scatters = ds['num_scatters'][:]

            # Write correct value to .db
            conn.execute(update_query, (int(num_scatters), row_id))

        # Save and close
        conn.commit()
        conn.close()
        ds.close()

if __name__ == '__main__':
    base_dir = os.path.join(
        os.environ['SCRATCH'],
        'kelp-results',
    )
    study_name = 'as_single_tol12'
    study_dir = os.path.join(base_dir, study_name)
    # Rewrite individual dbs.
    rewrite_num_scatters(base_dir, study_name)
    # Rewrite combined db.
    print("Recombining")
    kelp_compute.combine_dbs(study_dir, study_name)
    print("Done")
    
