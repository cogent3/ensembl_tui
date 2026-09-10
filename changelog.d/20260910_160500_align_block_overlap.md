### Contributors

- Katherine Caley, for reporting and diagnosing the alignment block overlap bug

### BUG

- `AlignDb.get_records_matching()` now returns every alignment block that
  overlaps the requested interval. It was testing whether the query start or
  stop fell inside a block, which matched only the two blocks containing the
  interval's endpoints and silently discarded every block lying wholly within
  it. Any region spanning more than one block was affected: for a 35.9 kb
  *D. melanogaster* gene tiled by 363 blocks, 2 were returned and 98.8% of the
  alignment was lost. Coverage scaled inversely with region length, so short
  regions looked correct. Queries passing only one of start or stop are
  unchanged. Callers of `eti alignments` should re-export.
