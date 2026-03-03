"""Entry point for python -m progenitors and the progenitors console script.

Settings: default metadata and features in progenitors/settings/pipeline.py
"""
import sys
import time
from .pipeline import main as pipeline_main
from . import options
from .sheets import sheetproc
from .env_config import validate_env
from .settings.pipeline import DEFAULT_METADATA, DEFAULT_FEATURES


def main():
    """Parse CLI and run the progenitor pipeline."""
    start = time.time()
    args = options.parse_arguments()
    metadata = args.metadata.split(",") if args.metadata else list(DEFAULT_METADATA)
    metadata = [m.strip().lower() for m in metadata]

    # Require auth env vars up front; exit with clear message if missing
    features = list(DEFAULT_FEATURES)
    if args.alert:
        features.append('email')
    if 'tns' in metadata:
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
    )
    total_time = time.time() - start
    msg = f'Finished with: {command}\nIt took {total_time} seconds to complete this script.'
    options.message(msg)


if __name__ == '__main__':
    main()
