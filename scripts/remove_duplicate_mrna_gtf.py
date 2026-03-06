#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_gtf_attributes(attr_str):
    """
    Parse standard GTF attributes, e.g.
    transcript_id "g34.t2"; gene_id "g34";
    """
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attrs[key] = value.strip().strip('"')
    return attrs


def format_gtf_attributes(attrs):
    """
    Return a GTF-compliant attribute string.
    """
    parts = []
    for key, value in attrs.items():
        parts.append(f'{key} "{value}"')
    return "; ".join(parts) + ";"


def guess_transcript_id(feature, attr_str):
    """
    Extract transcript_id from column 9.

    Handles:
    1. Standard GTF:
       transcript_id "g34.t2"; gene_id "g34";
    2. Non-standard transcript/mRNA lines where column 9 is just:
       g34.t2
    """
    attrs = parse_gtf_attributes(attr_str)

    if "transcript_id" in attrs:
        return attrs["transcript_id"]

    if feature in {"transcript", "mRNA"}:
        raw = attr_str.strip()
        if raw and " " not in raw and ";" not in raw and '"' not in raw:
            return raw

    return None


def guess_gene_id(attr_str):
    attrs = parse_gtf_attributes(attr_str)
    return attrs.get("gene_id")


def first_pass_collect_duplicates(gtf_file):
    """
    First pass:
    collect transcript and mRNA coordinates by transcript_id.
    """
    feats = defaultdict(lambda: {"transcript": set(), "mRNA": set()})

    with open_maybe_gzip(gtf_file, "rt") as fh:
        for line_num, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                print(
                    f"WARNING: malformed line {line_num}, expected 9 columns",
                    file=sys.stderr,
                )
                continue

            seqname, source, feature, start, end, score, strand, frame, attr_str = cols

            if feature not in {"transcript", "mRNA"}:
                continue

            transcript_id = guess_transcript_id(feature, attr_str)
            if transcript_id is None:
                print(
                    f"WARNING: could not determine transcript_id at line {line_num}",
                    file=sys.stderr,
                )
                continue

            coord = (seqname, start, end, strand)
            feats[transcript_id][feature].add(coord)

    duplicates = set()
    duplicate_info = defaultdict(list)

    for transcript_id, d in feats.items():
        overlap = d["transcript"] & d["mRNA"]
        if overlap:
            duplicates.add(transcript_id)
            duplicate_info[transcript_id].extend(sorted(overlap))

    return duplicates, duplicate_info


def second_pass_write_clean_gtf(gtf_file, output_file, report_file, duplicates, duplicate_info):
    """
    Second pass:
    - remove duplicate mRNA lines
    - normalize transcript/mRNA attribute column to GTF-compliant format when possible
    - preserve all other lines
    """
    removed_count = 0
    written_count = 0

    with open_maybe_gzip(gtf_file, "rt") as in_fh, open_maybe_gzip(output_file, "wt") as out_fh, open(
        report_file, "w"
    ) as rep_fh:

        rep_fh.write("transcript_id\tgene_id\tseqname\tstart\tend\tstrand\tremoved_feature\tline_number\n")

        for line_num, line in enumerate(in_fh, start=1):
            raw_line = line.rstrip("\n")

            if not raw_line:
                out_fh.write("\n")
                continue

            if raw_line.startswith("#"):
                out_fh.write(raw_line + "\n")
                continue

            cols = raw_line.split("\t")
            if len(cols) != 9:
                out_fh.write(raw_line + "\n")
                continue

            seqname, source, feature, start, end, score, strand, frame, attr_str = cols

            transcript_id = None
            gene_id = guess_gene_id(attr_str)

            if feature in {"transcript", "mRNA"}:
                transcript_id = guess_transcript_id(feature, attr_str)

            # remove duplicate mRNA line
            if feature == "mRNA" and transcript_id in duplicates:
                coord = (seqname, start, end, strand)
                if coord in duplicate_info[transcript_id]:
                    rep_fh.write(
                        f"{transcript_id}\t{gene_id or '.'}\t{seqname}\t{start}\t{end}\t{strand}\tmRNA\t{line_num}\n"
                    )
                    removed_count += 1
                    continue

            # make transcript/mRNA lines GTF compliant when possible
            if feature in {"transcript", "mRNA"}:
                attrs = parse_gtf_attributes(attr_str)
                raw_id = guess_transcript_id(feature, attr_str)

                if "transcript_id" not in attrs and raw_id is not None:
                    attrs["transcript_id"] = raw_id

                # try to infer gene_id from transcript_id if absent and transcript_id looks like gene.t1
                if "gene_id" not in attrs and raw_id is not None:
                    if ".t" in raw_id:
                        attrs["gene_id"] = raw_id.rsplit(".t", 1)[0]

                if attrs:
                    cols[8] = format_gtf_attributes(attrs)

            out_fh.write("\t".join(cols) + "\n")
            written_count += 1

    return removed_count, written_count


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remove duplicate mRNA lines when the same transcript_id has both "
            "'transcript' and 'mRNA' at identical coordinates, producing a cleaned GTF."
        )
    )
    parser.add_argument("input_gtf", help="Input GTF file (.gtf or .gtf.gz)")
    parser.add_argument("output_gtf", help="Output cleaned GTF file")
    parser.add_argument("report_tsv", help="Output report TSV with removed duplicates")
    args = parser.parse_args()

    duplicates, duplicate_info = first_pass_collect_duplicates(args.input_gtf)

    removed_count, written_count = second_pass_write_clean_gtf(
        args.input_gtf,
        args.output_gtf,
        args.report_tsv,
        duplicates,
        duplicate_info,
    )

    print(f"Transcripts with transcript+mRNA duplicates: {len(duplicates)}", file=sys.stderr)
    print(f"Removed mRNA lines: {removed_count}", file=sys.stderr)
    print(f"Written lines to cleaned GTF: {written_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
