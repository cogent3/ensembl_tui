### BUG

- `AlignDb.get_records_matching()` now returns every alignment block that
  overlaps the requested interval. It was testing whether the query start or
  stop fell inside a block, which matched only the two blocks containing the
  interval's endpoints and silently discarded every block lying wholly within
  it.