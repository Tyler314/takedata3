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

