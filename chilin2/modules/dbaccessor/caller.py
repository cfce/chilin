"""
Currently, use old sqlite3 database, replace later
"""

import sqlite3

def dbread():  ## read in ChiLin background data
    """

    :return:
    """


def dbwrite():
    """

    :return:
    """

def db_fastqc(db):
    qc_db = sqlite3.connect(db).cursor()
    qc_db.execute("SELECT median_quality FROM fastqc_info")
    history_data = [float(i[0]) for i in qc_db.fetchall()]
    return history_data

def db_bwa(db):
    db = sqlite3.connect(db).cursor()
    db.execute("select map_ratio from mapping")
    historyData = [str(i[0]) for i in (db.fetchall())]
    return historyData
