"""Unmodified fold-only preprocessing, rebound to the new write boundary."""
from common import reference_module
_base=reference_module('background')
load_full_raw=_base.load_full_raw
raw_for_indices=_base.raw_for_indices
fit_background=_base.fit_background
apply_background=_base.apply_background
cached_background=_base.cached_background
