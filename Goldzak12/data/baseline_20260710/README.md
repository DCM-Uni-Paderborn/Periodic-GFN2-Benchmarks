# Previous LC12 dataset

This directory freezes the LC12/Goldzak12 tables reported on 2026-07-10 before
the benchmark was rerun with the synchronized CP2K/tblite build.

The CP2K outputs identify development revision `faf9aae912`. The exact linked
tblite source revision was not embedded in those outputs, and the original
binary path was subsequently rebuilt. Consequently, these files are suitable
as a numerical baseline but not as a fully reproducible build record. The new
run records source revisions, working-tree state, executable hashes, and the
complete protocol in `../build_provenance.json`.
