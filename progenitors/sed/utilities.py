"""
SED utilities: dust maps, coordinate parsing, number formatting.

Core helpers live in ``progenitors.formatting_utils`` and ``progenitors.sfd_dust``.
"""
from ..formatting_utils import is_number, parse_coord, round_to_significant_figures
from ..sfd_dust import import_dustmap, get_sfd

# Legacy name: significant-figure rounding (not fixed decimal places).
round_to_n = round_to_significant_figures
