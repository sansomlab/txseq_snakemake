"""database.py - database upload for ruffus pipelines
=========================================================

Reference
---------

"""
import re
import os
import sqlite3
import sqlalchemy
from pathlib import Path
import txseq_snakemake.functions_from_cgat.database as database
import txseq_snakemake.functions_from_cgat.experiment as E
from txseq_snakemake.functions_from_cgat.iotools import snip, touch_file
#from cgatcore.pipeline.files import get_temp_file
#from cgatcore.pipeline.execution import run
from txseq_snakemake.functions_from_cgat.pipeline.parameters import get_params
#import txseq_snakemake.functions_from_cgat.tables as c_tables
#import txseq_snakemake.functions_from_cgat.csv2db as c_csv2db

def tablequote(track):
    '''quote a track name such that is suitable as a table name.'''
    return re.sub(r"[-(),\[\].]", "_", track)


def to_table(outfile):
    '''convert a filename from a load statement into a table name.

    This method checks if the filename ends with ".load". The suffix
    is then removed and the filename quoted so that it is suitable
    as a table name.

    Arguments
    ---------
    outfile : string
        A filename ending in ".load".

    Returns
    -------
    tablename : string

    '''
    assert outfile.endswith(".load")
    name = os.path.basename(outfile[:-len(".load")])
    return tablequote(name)


def build_load_statement(tablename, retry=True, options=""):
    """build a command line statement to upload data.

    Upload is performed via the :doc:`csv2db` script.

    The returned statement is suitable to use in pipe expression.
    This method is aware of the configuration values for database
    access and the chosen database backend.

    For example::

        load_statement = P.build_load_statement("data")
        statement = "cat data.txt | %(load_statement)s"
        P.run(statement)

    Arguments
    ---------
    tablename : string
        Tablename for upload
    retry : bool
        Add the ``--retry`` option to `csv2db.py`
    options : string
        Command line options to be passed on to `csv2db.py`

    Returns
    -------
    string

    """

    opts = []

    if retry:
        opts.append(" --retry ")

    params = get_params()
    opts.append("--database-url={}".format(params["database"]["url"]))

    db_options = " ".join(opts)
    load_statement = (
        "python -m txseq_snakemake.functions_from_cgat.csv2db {db_options} {options} --table={tablename}".format(**locals()))

    return load_statement


def load(infile,
         outfile=None,
         options="",
         tablename=None,
         retry=True,
         job_memory=None):
    """import data from a tab-separated file into database.

    The table name is given by outfile without the
    ".load" suffix.

    A typical load task in ruffus would look like this::

        @transform("*.tsv.gz", suffix(".tsv.gz"), ".load")
        def loadData(infile, outfile):
            P.load(infile, outfile)

    Upload is performed via the :doc:`csv2db` script.

    Arguments
    ---------
    infile : string
        Filename of the input data
    outfile : string
        Output filename. This will contain the logging information. The
        table name is derived from `outfile` if `tablename` is not set.
    options : string
        Command line options for the `csv2db.py` script.
    retry : bool
        If True, multiple attempts will be made if the data can
        not be loaded at the first try, for example if a table is locked.
    job_memory : string
        Amount of memory to allocate for job. If unset, uses the global
        default.
    """

    if job_memory is None:
        job_memory = "8G"

    if not tablename:
        tablename = to_table(outfile)

    statement = []

    if infile.endswith(".gz"):
        statement.append("zcat {infile}")
    else:
        statement.append("cat {infile}")

    statement.append(build_load_statement(tablename,
                                          options=options,
                                          retry=retry))

    statement = " | ".join(statement) + " > {outfile}"


    statement = statement.format(
        infile=infile,
        outfile=outfile,
    )

    return(statement)


def concatenate_and_load(infiles,
                         outfile,
                         regex_filename=None,
                         header=None,
                         cat="track",
                         has_titles=True,
                         missing_value="na",
                         retry=True,
                         tablename=None,
                         options="",
                         job_memory=None):
    """concatenate multiple tab-separated files and upload into database.

    The table name is given by outfile without the
    ".load" suffix.

    A typical concatenate and load task in ruffus would look like this::

        @merge("*.tsv.gz", ".load")
        def loadData(infile, outfile):
            P.concatenateAndLoad(infiles, outfile)

    Upload is performed via the :doc:`csv2db` script.

    Arguments
    ---------
    infiles : list
        Filenames of the input data
    outfile : string
        Output filename. This will contain the logging information. The
        table name is derived from `outfile`.
    regex_filename : string
        If given, *regex_filename* is applied to the filename to extract
        the track name. If the pattern contains multiple groups, they are
        added as additional columns. For example, if `cat` is set to
        ``track,method`` and `regex_filename` is ``(.*)_(.*).tsv.gz``
        it will add the columns ``track`` and method to the table.
    header : string
        Comma-separated list of values for header.
    cat : string
        Column title for column containing the track name. The track name
        is derived from the filename, see `regex_filename`.
    has_titles : bool
        If True, files are expected to have column titles in their first row.
    missing_value : string
        String to use for missing values.
    retry : bool
        If True, multiple attempts will be made if the data can
        not be loaded at the first try, for example if a table is locked.
    tablename: string
        Name to use for table. If unset derive from outfile.
    options : string
        Command line options for the `csv2db.py` script.
    job_memory : string
        Amount of memory to allocate for job. If unset, uses the global
        default. 
    """
    if job_memory is None:
        job_memory = "8G"

    if tablename is None:
        tablename = to_table(outfile)

   

    passed_options = options
    load_options, cat_options = ["--add-index=track"], []

    if regex_filename:
        cat_options.append("--regex-filename='%s'" % regex_filename)

    if header:
        load_options.append("--header-names=%s" % header)

    if not has_titles:
        cat_options.append("--no-titles")

    cat_options = " ".join(cat_options)
    load_options = " ".join(load_options) + " " + passed_options

    load_statement = build_load_statement(tablename,
                                          options=load_options,
                                          retry=retry)

    statement = ('''python -m txseq_snakemake.functions_from_cgat.tables \
    --cat={cat} \
    --missing-value={missing_value} \
    {cat_options} \
    {infiles} \
    | {load_statement} \
    > {outfile}''').format(
    cat=cat,
    missing_value=missing_value,
    cat_options=cat_options,
    infiles=" ".join(infiles),
    load_statement=load_statement,
    outfile=outfile,
    )

    return(statement)

