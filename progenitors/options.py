"""CLI argument parsing for the progenitor pipeline.

Settings: default --metadata and --yse-sql-query in progenitors/settings/pipeline.py
"""
import argparse
import logging
import sys

from .formatting_utils import is_number, parse_coord as _parse_coord
from .settings.pipeline import DEFAULT_METADATA, DEFAULT_YSE_SQL_QUERY


def message(msg):
    """Log completion banner (single line when possible)."""
    log = logging.getLogger("progenitors.cli")
    one_line = " ".join(line.strip() for line in msg.splitlines() if line.strip())
    log.info("done: %s", one_line)


def parse_two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError(None, "expected two space-separated floats")
    values = list(map(float, values))
    return tuple(values)


def parse_coord(ra, dec):
    """Parse coordinates; print errors to stderr-style console (CLI)."""
    ra_s = str(ra) if ra is not None else ''
    dec_s = str(dec) if dec is not None else ''
    log = logging.getLogger("progenitors.options")
    if (not (is_number(ra) and is_number(dec)) and
            (':' not in ra_s and ':' not in dec_s)):
        log.warning("coordinates not interpretable: %s %s", ra, dec)
        return None
    coord = _parse_coord(ra, dec)
    if coord is None:
        log.warning("coordinate parse failed: %s %s", ra, dec)
    return coord

def parse_arguments(usage=''):

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument(
        '--log-level',
        default=None,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO or PROGENITORS_LOG_LEVEL).',
    )
    parser.add_argument('--redo', default=False, action='store_true',
                        help='Redo metadata for all objects in database.')
    parser.add_argument(
        '--alert',
        default=False,
        action='store_true',
        help='With --update-classification: email when a transient newly appears on '
             'a core science tab (Type Ia, Ib/c, II-P/II-L, IIn, IIb) and passes the '
             'progenitor-interest mask.',
    )
    parser.add_argument('--always-update', default=False, action='store_true',
                        help='Send an email alert if there are new YSE objs.')
    parser.add_argument('--update-classification', default=False,
                        action='store_true',
                        help='Update the classification of objects using the sheet '
                             'Classification column after name-only crossmatch to the '
                             'TNS public objects catalog (requires TNS_* env vars; '
                             'catalog is downloaded once at pipeline start).')
    parser.add_argument('--update-tns-class', default=False, action='store_true',
                        help='Ensure TNS metadata (object type, classification refs) is '
                             'fetched: adds "tns" to metadata types if not already set.')
    parser.add_argument("--yse-sql-query", default=DEFAULT_YSE_SQL_QUERY, type=str,
                        help="YSE SQL query to use for grabbing YSE data.")
    parser.add_argument('--trim-metadata', default=False, action='store_true',
                        help='Trim keys from metadata that are not in table.')
    parser.add_argument('--redo-obj', default=None, type=str,
                        help='Redo metadata for these objects.  Input should '+\
                        'be a comma-separated list of objs (1987A,1993J).')
    parser.add_argument('--metadata', default=None, type=str,
                        help='Only do metadata for these data types.')
    parser.add_argument(
        '--profile',
        default=False,
        action='store_true',
        help='Log per-stage wall-clock timings for pipeline profiling.',
    )

    args = parser.parse_args()

    return(args)
