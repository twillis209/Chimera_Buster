import pytest

from Chimera_buster.chimera_buster import chimera_buster, readFastq

A = "TTTAACCGGTTAACCGGTTAACCGGTTTT"
B = "AAACCGGTTAACCGGTTAACCGGTTAAA"
A2 = "TTTGGCCAATTGGCCAATTGGCCAATTTT"
B2 = "AAAGGTTCCAAGGTTCCAAGGTTCCAAA"


def final_record(name, fwd, rev, size, cluster, seq="ACGT"):
    return (
        f">centroid={name};strand=+;umi_fwd=1;umi_rev=0;umi_fwd_seq={fwd};"
        f"umi_rev_seq={rev};seq=N;seqs={size};clusterid={cluster}\n{seq}\n"
    )


def write_inputs(tmp_path, records, prelim_ids):
    final = tmp_path / "final.fa"
    final.write_text("".join(records))
    prelim = tmp_path / "prelim.fa"
    prelim.write_text(
        "".join(
            f">centroid=r{i};strand=+;seqs={10 * i};clusterid={i}\nACGT\n"
            for i in prelim_ids
        )
    )
    return str(final), str(prelim)


def read_list(prefix, suffix):
    with open(f"{prefix}_{suffix}") as fh:
        return sorted(line.strip() for line in fh if line.strip())


# 1_0 and 2_0 are one molecule split into two clusters; 3_0 shares one UMI
# with 1_0 (a chimera); 4_0 is an unrelated molecule.
RECORDS = [
    final_record("1_0", A, B, 5, 0),
    final_record("2_0", A, B, 3, 1),
    final_record("3_0", A, B2, 2, 2),
    final_record("4_0", A2, B2, 4, 3),
]


@pytest.mark.parametrize("check", ["True", "False"])
def test_split_molecule_is_not_kept_twice(tmp_path, check):
    final, prelim = write_inputs(tmp_path, RECORDS, [1, 2, 3, 4])
    out = str(tmp_path / "CB")
    chimera_buster(final, prelim, out, 1, check)
    assert read_list(out, "nonchimera_list.txt") == ["1_0", "4_0"]


def test_check_clusters_reports_duplicates_separately(tmp_path):
    final, prelim = write_inputs(tmp_path, RECORDS, [1, 2, 3, 4])
    out = str(tmp_path / "CB")
    chimera_buster(final, prelim, out, 1, "True")
    assert read_list(out, "duplicate_list.txt") == ["2_0"]
    assert read_list(out, "chimera_list.txt") == ["3_0"]


def test_without_check_clusters_duplicates_count_as_chimeras(tmp_path):
    final, prelim = write_inputs(tmp_path, RECORDS, [1, 2, 3, 4])
    out = str(tmp_path / "CB")
    chimera_buster(final, prelim, out, 1, "False")
    assert read_list(out, "duplicate_list.txt") == []
    assert read_list(out, "chimera_list.txt") == ["2_0", "3_0"]


def test_check_clusters_with_zero_mismatch(tmp_path):
    final, prelim = write_inputs(tmp_path, RECORDS, [1, 2, 3, 4])
    out = str(tmp_path / "CB")
    chimera_buster(final, prelim, out, 0, "True")
    assert read_list(out, "nonchimera_list.txt") == ["1_0", "4_0"]
    assert read_list(out, "duplicate_list.txt") == ["2_0"]
    assert read_list(out, "chimera_list.txt") == ["3_0"]


def test_empty_input_writes_empty_outputs(tmp_path):
    final, prelim = write_inputs(tmp_path, [], [])
    out = str(tmp_path / "CB")
    chimera_buster(final, prelim, out, 1, "True")
    for suffix in ("nonchimera_list.txt", "chimera_list.txt", "duplicate_list.txt"):
        assert read_list(out, suffix) == []


def test_wrapped_fasta_is_read_whole(tmp_path):
    path = tmp_path / "wrapped.fa"
    path.write_text(">one\nACGT\nTTGG\n>two\nCC\n\n")
    assert readFastq(str(path)) == [[">one", "ACGTTTGG"], [">two", "CC"]]


def test_missing_prelim_cluster_names_the_cluster(tmp_path):
    final, prelim = write_inputs(tmp_path, RECORDS, [1, 2, 4])
    with pytest.raises(ValueError, match="3_0"):
        chimera_buster(final, prelim, str(tmp_path / "CB"), 1, "False")
