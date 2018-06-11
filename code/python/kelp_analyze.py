# 3rd party
import netCDF4 as nc
import pandas as pd

def table_to_df(conn, table_name):
    select_cmd = 'SELECT * FROM {table_name}'.format(table_name=table_name)
    cursor = conn.execute(select_cmd)
    columns = [col[0] for col in cursor.description]
    data = [x for x in cursor]

    df = pd.DataFrame(data, columns=columns)
    return df

def query_results(conn, table_name, **kwargs):
    # Sanitize string values
    for key, val in kwargs.items():
        if isinstance(val, str):
            kwargs[key] = '"{}"'.format(val)

    # Form SQL WHERE clause
    where_clause = ' AND '.join([
        '{}={}'.format(key, value)
        for key, value in kwargs.items()
    ])

    # Form entire query
    query = '''SELECT data_path FROM {table_name}
    WHERE {where_clause}
    '''.format(
        table_name=table_name,
        where_clause=where_clause
    )

    # Execute query
    cursor = conn.execute(query)
    results_list = [x for x in cursor]

    # Extract matching datasets
    datasets = [nc.Dataset(results[-1]) for results in results_list]
    return datasets
