import pickle

import cogent3
import cogent3_h5seqs as c3h5
import duckdb
import numpy
import pytest
from cogent3.core.location import IndelMap

from ensembl_tui import _align as eti_align
from ensembl_tui import _annotation as eti_annots
from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _ingest_align as eti_ingest_align
from ensembl_tui import _species as eti_species


def make_gene_attr(records: list[dict]) -> eti_annots.GeneView:
    schema = (
        "stable_id TEXT",
        "biotype TEXT",
        "seqid TEXT",
        "coord_system_name TEXT",
        "start INTEGER",
        "stop INTEGER",
        "strand TINYINT",
        "canonical_transcript_id INTEGER",
        "symbol TEXT",
        "gene_id INTEGER",
        "description TEXT",
    )
    columns = [c.split()[0] for c in schema]
    sql = f"""CREATE TABLE IF NOT EXISTS gene_attr ({",".join(schema)})"""
    conn = duckdb.connect(":memory:")
    conn.sql(sql)
    sql = "INSERT INTO gene_attr VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    rows = [[r.get(c) for c in columns] for r in records]
    conn.executemany(sql, parameters=rows)
    return eti_annots.GeneView(source=":memory:", db=conn)


def get_annotation_db() -> dict[str, eti_annots.Annotations]:
    gene_attr = {
        "s1": make_gene_attr(
            [
                {
                    "seqid": "s1",
                    "coord_system_name": "chromosome",
                    "biotype": "protein_coding",
                    "stable_id": "not-on-s2",
                    "start": 4,
                    "stop": 7,
                },
            ],
        ),
        "s2": make_gene_attr(
            [
                {
                    "seqid": "s2",
                    "coord_system_name": "chromosome",
                    "biotype": "protein_coding",
                    "stable_id": "includes-s2-gap",
                    "start": 2,
                    "stop": 6,
                },
            ],
        ),
        "s3": make_gene_attr(
            [
                {
                    "seqid": "s3",
                    "coord_system_name": "chromosome",
                    "biotype": "protein_coding",
                    "stable_id": "includes-s3-gap",
                    "start": 22,
                    "stop": 27,
                },
            ],
        ),
    }
    for sp, gv in gene_attr.items():
        bt = eti_annots.BiotypeView(source=gv.source, db=gv.conn)
        gene_attr[sp] = eti_annots.Annotations(source=":memory:", genes=gv, biotypes=bt)
    return gene_attr


def small_seqs():
    from cogent3 import make_aligned_seqs

    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    return make_aligned_seqs(
        seqs,
        moltype="dna",
        info={"species": {"s1": "human", "s2": "mouse", "s3": "dog"}},
    )


def make_records(start, end, block_id):
    aln = small_seqs()[start:end]
    records = []
    species = aln.info.species
    for seq in aln.seqs:
        seqid, seq_start, seq_end, seq_strand = seq.parent_coordinates(seq_coords=True)
        gs = seq.gapped_seq
        imap, s = gs.parse_out_gaps()
        if imap.num_gaps:
            gap_spans = numpy.array(
                [imap.gap_pos, imap.cum_gap_lengths],
                dtype=numpy.int32,
            ).T
        else:
            gap_spans = numpy.array([], dtype=numpy.int32)
        record = eti_align.AlignRecord(
            source="blah",
            species=species[seq.name],
            block_id=block_id,
            seqid=seqid,
            start=seq_start,
            stop=seq_end,
            strand=seq_strand,
            gap_spans=gap_spans,
        )
        records.append(record)
    return records


@pytest.fixture
def small_records():
    return make_records(1, 5, 0)


def empty_align_agg_gap_store():
    return eti_ingest_align.make_alignment_aggregator_db()


def test_aligndb_records_match_input(small_records):
    import copy

    orig_records = copy.deepcopy(small_records)
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(records=small_records, conn=agg)
    db = eti_align.AlignDb(db=agg, source=":memory:")
    got = next(iter(db.get_records_matching(species="human", seqid="s1")))
    assert got == set(orig_records)


def test_aligndb_records_skip_duplicated_block_ids(small_records):
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(conn=agg, records=small_records)
    sql = "SELECT COUNT(*) FROM align_blocks"
    count = agg.sql(sql).fetchone()[0]
    eti_ingest_align.add_records(conn=agg, records=small_records)
    assert agg.sql(sql).fetchone()[0] == count


# fixture to make synthetic genome and alignment db
# based on a given alignment
@pytest.fixture
def genomedbs_aligndb(small_records):
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(records=small_records, conn=agg)
    align_db = eti_align.AlignDb(source=":memory:", db=agg)
    seqs = small_seqs().degap()
    species = seqs.info.species
    data = seqs.to_dict()
    genomes = {}
    for name, seq in data.items():
        genome = c3h5.make_unaligned(
            "memory",
            in_memory=True,
            mode="w",
            alphabet=eti_genome.alphabet,
        )
        genome.add_seqs(seqs={name: seq})
        genomes[species[name]] = cogent3.make_unaligned_seqs(
            genome,
            annotation_db=None,
            info={"species": species[name]},
            moltype="dna",
        )

    return genomes, align_db


def test_building_alignment(genomedbs_aligndb, namer):
    genomes, align_db = genomedbs_aligndb
    got = next(
        iter(
            eti_align.get_alignment(
                align_db,
                genomes,
                ref_species="mouse",
                seqid="s2",
                namer=namer,
            ),
        ),
    )
    orig = small_seqs()[1:5]
    assert got.to_dict() == orig.to_dict()
    got.annotation_db.close()


@pytest.mark.parametrize(
    "kwargs",
    ({"ref_species": "dodo", "seqid": "s2"},),
)
def test_building_alignment_invalid_details(genomedbs_aligndb, kwargs):
    genomes, align_db = genomedbs_aligndb
    with pytest.raises(ValueError):
        list(eti_align.get_alignment(align_db, genomes, **kwargs))


def make_sample(two_aligns=False):
    aln = small_seqs()
    species = aln.info.species
    # make annotation db's
    annot_dbs = get_annotation_db()

    # we will reverse complement the s2 genome compared to the original
    # this means our coordinates for alignment records from that genome
    # also need to be rc'ed
    genomes = {}
    for seq in aln.seqs:
        name = seq.name
        seq = seq.seq
        if seq.name == "s2":
            seq = seq.rc()
            s2_genome = str(seq)
        genome = c3h5.make_unaligned(
            "memory",
            in_memory=True,
            mode="w",
            alphabet=eti_genome.alphabet,
        )
        genome.add_seqs(seqs={name: str(seq)})

        genomes[species[name]] = cogent3.make_unaligned_seqs(
            genome,
            annotation_db=annot_dbs[name],
            info={"species": species[name]},
            moltype="dna",
        )

    # define two alignment blocks that incorporate features
    align_records = _update_records(s2_genome, aln, "0", 1, 12)
    if two_aligns:
        align_records += _update_records(s2_genome, aln, "1", 22, 30)

    maker = eti_ingest_align.make_alignment_aggregator_db()
    eti_ingest_align.add_records(conn=maker, records=align_records)
    align_db = eti_align.AlignDb(db=maker, source=":memory:")
    return genomes, align_db


def _update_records(s2_genome, aln, block_id, start, end):
    # start, stop are the coordinates used to slice the alignment
    align_records = make_records(start, end, block_id)
    # in the alignment, s2 is in reverse complement relative to its genome
    # In order to be sure what "genome" coordinates are for s2, we first slice
    # the alignment
    aln = aln[start:end]
    # then get the ungapped sequence
    s2 = aln.get_seq("s2")
    # and reverse complement it ...
    selected = s2.rc()
    # so we can get the genome coordinates for this segment on the s2 genome
    start = s2_genome.find(str(selected))
    end = start + len(selected)
    for record in align_records:
        if record.seqid == "s2":
            record.start = start
            record.stop = end
            record.strand = -1
            break
    return align_records


@pytest.mark.parametrize(
    "start_end",
    [
        (None, None),
        (None, 11),
        (3, None),
        (3, 13),
    ],
)
@pytest.mark.parametrize(
    "species_coord",
    [
        ("human", "s1"),
        ("dog", "s3"),
    ],
)
def test_select_alignment_plus_strand(species_coord, start_end, namer):
    species, seqid = species_coord
    start, end = start_end
    aln = small_seqs()
    # the sample alignment db has an alignment block that
    # starts at 1 and ends at 12. The following slice is to
    # get the expected answer
    expect = aln[max(1, start or 1) : min(end or 12, 12)]
    # one sequence is stored in reverse complement
    genomes, align_db = make_sample()
    got = list(
        eti_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    assert len(got) == 1
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "start_end",
    [
        (None, None),
        (None, 5),
        (2, None),
        (2, 7),
    ],
)
def test_select_alignment_minus_strand(start_end, namer):
    species, seqid = "mouse", "s2"
    start, end = start_end
    aln = small_seqs()
    ft = aln.add_feature(
        biotype="custom",
        name="selected",
        seqid="s2",
        spans=[(max(1, start or 0), min(end or 12, 12))],
    )
    expect = aln[ft.map.start : min(ft.map.end, 12)]
    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query
    s2 = aln.get_seq("s2")
    s2_ft = next(iter(s2.get_features(name="selected")))
    if not any([start is None, end is None]):
        start = len(s2) - s2_ft.map.end
        end = len(s2) - s2_ft.map.start
    elif start == None != end:  # noqa E711
        start = len(s2) - s2_ft.map.end
        end = None
    elif start != None == end:  # noqa E711
        end = len(s2) - s2_ft.map.start
        start = None

    aln.annotation_db.close()

    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query

    genomes, align_db = make_sample(two_aligns=False)
    got = list(
        eti_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    # drop the strand info
    assert len(got) == 1, f"{s2_ft=}"
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_get_alignment_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    genomes, align_db = make_sample(two_aligns=False)
    got = next(
        iter(eti_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs)),
    )
    assert len(got.annotation_db) == 3


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_get_alignment_masked_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    kwargs["mask_features"] = ["gene"]
    genomes, align_db = make_sample(two_aligns=False)
    got = next(
        iter(eti_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs)),
    )
    assert len(got.annotation_db) == 3


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_align_db_get_records(coord):
    kwargs = dict(zip(("species", "seqid", "start", "stop"), coord, strict=False))
    # records are, we should get a single hit from each query
    # [('blah', 0, 'human', 's1', 1, 12, '+', array([], dtype=int32)),
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 1


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1"),
        ("mouse", "s2"),
        ("dog", "s3"),
    ],
)
def test_align_db_get_records_required_only(coord):
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
    # two hits for each species
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 2


def test_align_db_get_records_block_within_query():
    """blocks lying wholly inside the query interval must be returned

    Previously only blocks containing the query start or the query stop
    were matched, so all interior blocks were silently dropped.
    """
    agg = empty_align_agg_gap_store()
    # three consecutive, non-overlapping blocks tiling s1
    num_blocks = 3
    for i, (start, end) in enumerate(((1, 3), (3, 5), (5, 7))):
        eti_ingest_align.add_records(
            records=make_records(start, end, i),
            conn=agg,
        )
    db = eti_align.AlignDb(db=agg, source=":memory:")
    human = next(r for r in make_records(1, 7, 99) if r.species == "human")
    # a query spanning all three blocks, so the middle one is interior
    got = list(
        db.get_records_matching(
            species="human",
            seqid=human.seqid,
            start=human.start,
            stop=human.stop,
        ),
    )
    assert len(got) == num_blocks


def test_align_db_get_records_excludes_adjacent_blocks():
    """half-open bounds, so a block merely abutting the query is not a match"""
    agg = empty_align_agg_gap_store()
    for i, (start, end) in enumerate(((1, 3), (3, 5))):
        eti_ingest_align.add_records(records=make_records(start, end, i), conn=agg)
    db = eti_align.AlignDb(db=agg, source=":memory:")
    first = next(r for r in make_records(1, 3, 0) if r.species == "human")
    got = list(
        db.get_records_matching(
            species="human",
            seqid=first.seqid,
            start=first.start,
            stop=first.stop,
        ),
    )
    assert len(got) == 1


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s2"),
        ("mouse", "xx"),
        ("blah", "s3"),
    ],
)
def test_align_db_get_records_no_matches(coord):
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
    # no hits at all
    _, align_db = make_sample()
    got = list(align_db.get_records_matching(**kwargs))
    assert not got


def test_get_species():
    _, align_db = make_sample()
    assert set(align_db.get_species_names()) == {"dog", "human", "mouse"}


def test_write_alignments():
    genomes, align_db = make_sample(two_aligns=True)
    locations = eti_genome.get_gene_segments(
        annot_db=genomes["human"].annotation_db,
        species="human",
        stableids=["not-on-s2"],
    )
    app = eti_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
    )
    aln = app(locations[0])  # pylint: disable=not-callable
    assert len(aln[0]) == 3


@pytest.fixture
def db_align(DATA_DIR, tmp_dir):
    staging_path = tmp_dir / "staging"
    staging_path.mkdir(parents=True, exist_ok=True)
    align_name = "align_name"
    install_path = tmp_dir / "install"

    cfg = eti_config.Config(
        domain="main",
        staging_path=staging_path,
        install_path=install_path,
        species_dbs={},
        release="113",
        align_names=[align_name],
        tree_names=[],
        homologies=True,
        species_map=eti_species.make_species_map(None),
    )
    align_dir = cfg.staging_aligns / align_name
    align_dir.mkdir(parents=True, exist_ok=True)
    maf = align_dir / f"{align_name}.maf"
    maf.write_text((DATA_DIR / "tiny.maf").read_text())

    parquet_path = eti_ingest_align.install_alignment(cfg, align_name)
    return eti_align.AlignDb(source=parquet_path.parent)


def test_db_align(db_align):
    orig = len(db_align)
    source = db_align.source
    db_align.close()
    got = eti_align.AlignDb(source=source)
    assert len(got) == orig


@pytest.mark.parametrize("func", [str, repr])
def test_db_align_repr(db_align, func):
    got = func(db_align)
    assert "AlignDb" in got


def test_pickling_db(db_align):
    # should not fail
    pkl = pickle.dumps(db_align)  # nosec B301
    upkl = pickle.loads(pkl)  # nosec B301  # noqa: S301
    assert db_align.source == upkl.source


def test_aligndb_post_init_failure(tmp_path):
    with pytest.raises(FileNotFoundError):
        eti_align.AlignDb(source="/non/existent/directory")

    outfile = tmp_path / "textfile.txt"
    outfile.write_text("blah")
    with pytest.raises(OSError):
        eti_align.AlignDb(source="/non/existent/directory")


def test_aligndb_close(db_align):
    db_align.close()
    with pytest.raises(duckdb.ConnectionException):
        db_align.num_records()


def test_load_align_records():
    maf_record = {
        "species": "chlorocebus_sabaeus",
        "seqid": "28",
        "start": 2617243,
        "stop": 2645355,
        "strand": "-",
        "block_id": 20060000081647,
        "source": "10_primates.epo.other_6.maf.gz",
        # mixed case seq because this must be upper cased in seq2gap for cogent3
        # to deem it a valid seq
        "seq": "-acCAGGAAG",
    }
    got = eti_ingest_align.seq2gaps(maf_record.copy())
    assert (got.gap_spans == numpy.array([[0, 1]], dtype=numpy.int32)).all()


def _maf_record(seq: str, **kwargs) -> dict:
    record = {
        "species": "chlorocebus_sabaeus",
        "seqid": "28",
        "start": 0,
        "stop": len(seq) - seq.count("-"),
        "strand": "+",
        "block_id": 20060000081647,
        "source": "10_primates.epo.other_6.maf.gz",
        "seq": seq,
    }
    return record | kwargs


def test_seq2gaps_multiple_gaps():
    # gap_spans holds cumulative gap lengths, matching what IndelMap keeps
    # internally, so the second gap of 3 is stored as 2 + 3
    got = eti_ingest_align.seq2gaps(_maf_record("AC--GTA---CC"))
    assert (got.gap_spans == numpy.array([[2, 2], [5, 5]], dtype=numpy.int32)).all()


@pytest.mark.parametrize(
    "seq",
    (
        "AC--GTA---CC",
        "--ACGTA---CC",
        "AC--GTACC---",
        "A-C-G-T-A-C-",
        "-A-CG--TAC-C--G",
    ),
)
def test_seq2gaps_indelmap_roundtrip(seq):
    # get_alignment rebuilds an IndelMap from gap_spans this way, so the
    # rebuilt map must have the number of columns we started with
    record = eti_ingest_align.seq2gaps(_maf_record(seq))
    gap_pos, cum_gap_lengths = record.cum_gap_data
    imap = IndelMap(
        gap_pos=gap_pos,
        cum_gap_lengths=cum_gap_lengths,
        parent_length=record.stop - record.start,
    )
    assert len(imap) == len(seq)


def multi_gap_seqs():
    # s2 and s3 have two gaps each, which is what tells a cumulative gap
    # length apart from a per-gap one
    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGT--CAAATAATGAA",
        "s3": "GCTGAAGT--TGGAAGTTGCAAAT---GAA",
    }
    return cogent3.make_aligned_seqs(
        seqs,
        moltype="dna",
        info={"species": {"s1": "human", "s2": "mouse", "s3": "dog"}},
    )


@pytest.fixture
def ingested_multi_gap_aligndb():
    """genomes plus an AlignDb built through the real ingest path"""
    aln = multi_gap_seqs()
    species = aln.info.species
    records = [
        eti_ingest_align.seq2gaps(
            _maf_record(
                str(seq.gapped_seq),
                species=species[seq.name],
                seqid=seq.name,
                block_id=0,
                source="blah",
            ),
        )
        for seq in aln.seqs
    ]
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(records=records, conn=agg)
    align_db = eti_align.AlignDb(source=":memory:", db=agg)

    genomes = {}
    for name, seq in aln.degap().to_dict().items():
        genome = c3h5.make_unaligned(
            "memory",
            in_memory=True,
            mode="w",
            alphabet=eti_genome.alphabet,
        )
        genome.add_seqs(seqs={name: seq})
        genomes[species[name]] = cogent3.make_unaligned_seqs(
            genome,
            annotation_db=None,
            info={"species": species[name]},
            moltype="dna",
        )

    return genomes, align_db


@pytest.mark.parametrize(
    ("ref_species", "seqid"),
    (("human", "s1"), ("mouse", "s2"), ("dog", "s3")),
)
def test_building_alignment_multiple_gaps(
    ingested_multi_gap_aligndb,
    ref_species,
    seqid,
):
    # round trip: ingest writes the gaps, get_alignment reads them back, and
    # the alignment we get out is the one we put in. mouse and dog are the
    # references with two gaps of their own
    genomes, align_db = ingested_multi_gap_aligndb
    got = next(
        iter(
            eti_align.get_alignment(
                align_db,
                genomes,
                ref_species=ref_species,
                seqid=seqid,
            ),
        ),
    )
    orig = multi_gap_seqs()
    assert len(got) == len(orig)
    assert sorted(got.to_dict().values()) == sorted(orig.to_dict().values())
    got.annotation_db.close()


# tiny.maf has a single gap per record, so it cannot tell a per-gap length
# apart from a running total. these records have three, two and none.
# the trailing blank line is what real Ensembl maf files have
MULTI_GAP_MAF = """\
##maf version=1
# NOTE: sample written for testing purposes
# id: 20060000040557
a
s homo_sapiens.1        100    7 +  1000 AC--GTA---CC
s mouse.2               200   12 +  2000 ACTGGTAGGTCC
s rat.3                 300    8 +  3000 -CTGG--AGTC-

"""

ALIGN_WIDTH = 12


def test_install_path_multiple_gaps(tmp_path):
    # the whole ingest path, maf text through to AlignRecord. every record in
    # a block has to rebuild to the same number of alignment columns
    path = tmp_path / "multi_gap.maf"
    path.write_text(MULTI_GAP_MAF)
    records = eti_ingest_align.load_align_records()(path)
    assert len(records) == 3
    for record in records:
        gap_pos, cum_gap_lengths = record.cum_gap_data
        imap = IndelMap(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_gap_lengths,
            parent_length=record.stop - record.start,
        )
        assert len(imap) == ALIGN_WIDTH, record.species


def test_aln_seq_matches_genome(apes, apes_aligndb):
    genomes, align_db = apes, apes_aligndb
    hsap = genomes["homo_sapiens"]
    # now get an alignment for just 40 nucleotides, around the gene_name
    start = 37823381
    stop = start + 40
    gene_name = "ENSG00000128310"
    raw_seq_rc = str(hsap.seqs["22"][start:stop].rc())
    locus = eti_genome.genome_segment(
        species="homo_sapiens",
        seqid="22",
        start=start,
        stop=stop,
        strand=1,
        unique_id=gene_name,
    )
    maker = eti_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
        mask_features=None,
        shadow=None,
        mask_ref=True,
    )
    alns = maker.main(locus)
    assert len(alns) == 1
    aln = alns[0]
    expect_name = f"homo_sapiens:22:{start}-{stop}:-1"
    got = aln.get_seq(expect_name)
    assert str(got) == raw_seq_rc


def _get_annotated_align(
    genomes,
    align_db,
    species: str,
    seqid: str,
    strand: int,
    start: int,
    stop: int,
    unique_id: str,
    biotype: str,
    rc: bool = False,
):
    genome = genomes[species]
    features = list(
        genome.get_features(seqid="22", biotype=biotype, start=start, stop=stop),
    )
    feature = features[0].get_slice()
    if rc:
        feature = feature.rc()
    ft_seq = str(feature)
    ft_len = len(ft_seq)
    locus = eti_genome.genome_segment(
        species=species,
        seqid=seqid,
        start=start,  # repeat is 10bp from here
        stop=stop,  # ends 5 bp from here
        unique_id=unique_id,
        strand=strand,
    )
    build_align = eti_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
        mask_features=None,
        shadow=None,
        mask_ref=True,
    )
    alns = build_align.main(locus)
    assert len(alns) == 1
    aln = alns[0]

    exp_name = f"{species}:22:{start}-{stop}:{-1 if rc else 1}"
    unmasked = str(aln.seqs[exp_name].seq)
    # feature present in the unmasked sequence
    assert ft_seq in unmasked
    ft_start = unmasked.find(ft_seq)
    ft_end = ft_start + len(ft_seq)

    build_align = eti_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
        mask_features=biotype,
        shadow=None,
        mask_ref=True,
    )
    alns = build_align.main(locus)
    assert len(alns) == 1
    aln = alns[0]
    masked = str(aln.seqs[exp_name].seq)
    assert ft_seq not in masked
    assert masked.count("?") == ft_len
    # check masked region correct
    assert masked[:ft_start] == unmasked[:ft_start]
    assert masked[ft_end:] == unmasked[ft_end:]


def test_aln_cds_masked_repeat_plus_strand(apes, apes_aligndb):
    genomes, align_db = apes, apes_aligndb
    species = "homo_sapiens"
    # ENSG00000100346 has one simple repeat
    start = 39601558 - 10  # repeat is 10bp from here
    stop = 39601622 + 5  # ends 5 bp from here
    _get_annotated_align(
        genomes=genomes,
        align_db=align_db,
        species=species,
        seqid="22",
        start=start,
        stop=stop,
        strand=1,
        unique_id="blah",
        biotype="Simple_repeat",
        rc=False,
    )


def test_aln_cds_masked_repeat_reverse_complemented(apes, apes_aligndb):
    genomes, align_db = apes, apes_aligndb
    species = "homo_sapiens"
    # ENSG00000100412 has one simple repeat, alignment minus strand
    start = 41472602 - 10  # repeat is 10bp from here
    stop = 41472627 + 5  # ends 5 bp from here
    _get_annotated_align(
        genomes=genomes,
        align_db=align_db,
        species=species,
        seqid="22",
        start=start,
        stop=stop,
        strand=1,
        unique_id="blah",
        biotype="Simple_repeat",
        rc=True,
    )
