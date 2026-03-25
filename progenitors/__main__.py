"""Entry point for python -m progenitors and the progenitors console script.

Settings: default metadata and features in progenitors/settings/pipeline.py
"""
import logging
import sys
import time

from ._version import __version__
from .pipeline import main as pipeline_main
from . import options
from .pipeline_logging import setup_pipeline_logging
from .sheets import sheetproc
from .env_config import validate_env
from .settings.pipeline import DEFAULT_METADATA, DEFAULT_FEATURES


def main():
    """Parse CLI and run the progenitor pipeline."""
    start = time.time()
    args = options.parse_arguments()
    setup_pipeline_logging(level=args.log_level)
    log = logging.getLogger("progenitors.cli")
    log.info("progenitors %s", __version__)

    metadata = args.metadata.split(",") if args.metadata else list(DEFAULT_METADATA)
    metadata = [m.strip().lower() for m in metadata]
    if args.update_tns_class and 'tns' not in metadata:
        metadata.append('tns')

    # Require auth env vars up front; exit with clear message if missing
    features = list(DEFAULT_FEATURES)
    if args.alert:
        features.append('email')
    if 'tns' in metadata:
        features.append('tns')
    if args.update_classification and 'tns' not in features:
        features.append('tns')
    if 'ads' in metadata:
        features.append('ads')
    validate_env(features, credentials_path=sheetproc.params['credentials'])

    command = ' '.join(sys.argv)
    pipeline_main(
        redo=args.redo,
        update_classification=args.update_classification,
        alert=args.alert,
        yse_sql_query=args.yse_sql_query,
        always_update=args.always_update,
        trim_metadata=args.trim_metadata,
        redo_obj=args.redo_obj,
        metadata=metadata,
        profile=args.profile,
    )
    total_time = time.time() - start
    msg = f'Finished with: {command}\nIt took {total_time} seconds to complete this script.'
    options.message(msg)


if __name__ == '__main__':
    main()
