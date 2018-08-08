# Rename study.
# - study_dir name
# - combined .db name
# - table names in data_dir
# - combined table name, if it exists

import os
import sqlite3
import netCDF4 as nc
import kelp_compute
import re

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

        # Write correct value to .db
        print(
            "Executing '{}' with".format(update_query),
            (new_data_path, row_id)
        )
        conn.execute(update_query, (new_data_path, row_id))

    conn.commit()


def rename_study(base_dir, old_study_name, new_study_name):
    old_study_dir = os.path.join(base_dir, old_study_name)
    new_study_dir = os.path.join(base_dir, new_study_name)
    data_dir = os.path.join(old_study_dir, 'data')

    # Make sure this study does not already exist
    if os.path.exists(new_study_dir):
        raise OSError(
            "Rename failed. Study '{}' already exist.".format(
                new_study_dir
            )
        )

    db_path_list = [
        os.path.join(
            data_dir,
            basename
        )
        for basename in os.listdir(data_dir)
        if basename[-3:] == '.db'
    ]

    # Modify individual .dbs
    for db_path in db_path_list:
        print(db_path)
        conn = sqlite3.connect(db_path)
        rewrite_data_path(conn, old_study_name, new_study_name)
        rename_table(conn, old_study_name, new_study_name)
        conn.close()

    # Modify combined .db
    old_combined_db_path = os.path.join(
        old_study_dir,
        '{}.db'.format(old_study_name)
    )
    new_combined_db_path = os.path.join(
        old_study_dir,
        '{}.db'.format(new_study_name)
    )
    conn = sqlite3.connect(old_combined_db_path)
    rewrite_data_path(conn, old_study_name, new_study_name)
    rename_table(conn, old_study_name, new_study_name)
    os.rename(
        old_combined_db_path,
        new_combined_db_path
    )
    conn.close()

    # Rename directory
    os.rename(old_study_dir, new_study_dir)

if __name__ == '__main__':
    base_dir = os.path.join(
        os.environ['SCRATCH'],
        'kelp-results',
    )
    old_study_name = 'gs64_a01_b0_small'
    new_study_name = 'gs23_a01_b0'
    rename_study(base_dir, old_study_name, new_study_name)
