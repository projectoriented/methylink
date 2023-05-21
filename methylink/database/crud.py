import os
import sqlite3
import logging

LOG = logging.getLogger(__name__)

DB_NAME = "meth_tags.db"


class DataBaseAdapter:
    def __init__(self, db_name) -> None:
        self.db_name = db_name
        # Start database connection
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()


def create_database(db_name):
    if os.path.exists(db_name):
        LOG.info("Database already exists, removing.")
        os.remove(db_name)

    db = DataBaseAdapter(db_name=db_name)
    db.cur.execute("CREATE TABLE meth_tags (qname TEXT PRIMARY KEY, tag BLOB)")

    # commit changes and close connection
    db.conn.commit()
    db.conn.close()
    LOG.info(f"{db_name} created.")


def insert_one(qname, tag, db):
    db.cur.execute("INSERT INTO meth_tags VALUES (?, ?)", (qname, tag))

    LOG.debug(f"Inserted ({qname}, {tag}) into {DB_NAME}")

    db.conn.commit()
    db.conn.close()


def select_one(qname, db):
    result = db.cur.execute(
        "SELECT tag FROM meth_tags WHERE qname = ?", (qname,)
    ).fetchone()
    if result:
        LOG.debug(f"Found tag for {qname}")
        return result
    else:
        LOG.debug(f"Nothing returned for {qname}")
