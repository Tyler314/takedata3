import sqlite3


_DB_CONN = {}

def _get_conn(db_path):
    global _DB_CONN
    if db_path in _DB_CONN:
        return _DB_CONN[db_path]
    else:
        _DB_CONN[db_path] = sqlite3.Connection(
            db_path, check_same_thread=False)
        return _DB_CONN[db_path]


def fetchall(db_path, query, query_args, row_factory=None):
    """
    Makes a query lookup in the database, where the location is specified
    by 'db_path'. 'query' is the query itself, and it will have '?' characters
    inside of it, which serve as place holders for the elements within the
    tuple 'query_args'. It behaves the exact same way as if instead of '?', it
    were '{}', and you ran 'query.format(query_args)'.
    """
    conn = _get_conn(db_path)
    conn.row_factory=row_factory
    cursor = conn.cursor()
    cursor.execute(query, query_args)
    results = cursor.fetchall()
    cursor.close()
    return results


def run_commit(db_path, query, query_args):
    conn = _get_conn(db_path)
    cursor = conn.cursor()
    cursor.execute(query, query_args)
    conn.commit()
    cursor.close()

