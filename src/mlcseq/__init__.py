"""mlcseq - a cleaned-up, tested re-implementation of the MLC-Seq algorithms
(Yuan et al., JACS 2024) for direct RNA sequencing by mass spectrometry."""
from .basecalling import BaseCalls, base_calling, homology_search
from .complementation import Ladder, ladder_complementation
from .composition import gap_compositions, has_composition
from .gap_fill import gap_fill
from .io import load_sample
from .mass_sum import mass_sum
from .matching import annotate_matches, drop_dups, match_masses

__version__ = "0.1.0"

__all__ = [
    "BaseCalls", "base_calling", "homology_search", "Ladder", "ladder_complementation",
    "gap_compositions", "has_composition", "gap_fill", "load_sample", "mass_sum",
    "annotate_matches", "drop_dups", "match_masses",
]
