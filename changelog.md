
<a id='changelog-0.7.9'></a>
# Changes since release "0.7.9"

This is mainly a bug fix release.

## Contributors

- Bug report by @KatherineCaley, thanks Kath!
- Fixes by @GavinHuttley

## ENH

- The MAF parser now reads a file a chunk at a time, using
  `scinexus.io_util.iter_splitlines()`, instead of loading the whole file and
  splitting it into a list of lines. Parsing one 79 MB compressed alignment
  file drops from 1382 MB peak to 142 MB, and installing an alignment no
  longer holds a decompressed copy of each file in memory. Block ids and the
  sequences parsed out of them are unchanged.

## BUG

- `AlignDb.get_records_matching()` now returns every alignment block
  overlapping the query interval. It previously only matched blocks
  containing the query start or the query stop, so blocks lying wholly
  inside the interval were dropped. Queries bounded on one side only were
  worse affected, returning just the single block containing that bound.
- `get_alignment()` now picks the reference sequence segment that overlaps
  the query. An alignment block can contain several segments of the
  reference sequence, for example a duplicated region, and the segment was
  previously taken from an unordered set. That produced results which varied
  between runs, alignments with a stop before their start, and an
  `IndexError` from `IndelMap.get_align_index`. Where more than one segment
  of the reference overlaps the query, one alignment is now returned per
  segment.
- Installing an alignment no longer writes duplicate rows. Ensembl emits an
  alignment block once for each segment of the reference species it holds, so
  the same record reaches `add_records()` several times in one batch. The
  existing check only skipped block ids written by an earlier call, so those
  copies were all stored: 5,490 of 60,017 rows (9.1%) in a three species
  primate install. Query results were already correct, because the read path
  collapses identical records into a set, so an existing installed store does
  not need rebuilding.

- `AlignRecord.__eq__()` raised `ValueError` when comparing two records with
  the same coordinates but gap arrays of different shapes. Such records hash
  alike, so they met in the set built by `AlignDb.get_records_matching()`.
  They now compare with `numpy.array_equal()`, which tolerates a shape
  mismatch.

<a id='changelog-0.7.8'></a>
# Changes since release "0.7.8"

This is a major bug fix release.

## Contributors

- Ulises Hernandez, for reporting the gap placement bug

## BUG

- Gaps are placed correctly again when an alignment is read back. The gap
  lengths written by `eti install` are cumulative, matching how `IndelMap`
  stores them, but were being passed to it as per-gap lengths. Every gap after
  the first in a sequence came out too long and trailing bases were dropped,
  with no error raised. Alignments installed with 0.7.6 or 0.7.7 are read
  correctly now and do not need reinstalling. Any installed with 0.4.3 or
  earlier hold per-gap lengths and must be reinstalled with `eti install`.
- Extracting gaps during install no longer holds on to a scratch buffer the
  size of the whole gapped sequence for every record.
