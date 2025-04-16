#!/usr/bin/env python3

import argparse
from collections import defaultdict

# === Argument Parsing ===
parser = argparse.ArgumentParser(description='Assign contig names based on PAF alignment to reference chromosomes.')
parser.add_argument('-i', '--input', required=True, help='Input PAF file')
parser.add_argument('-o', '--output', required=True, help='Output file for renamed contigs')
parser.add_argument('-m', '--min_length', type=int, default=1000, help='Minimum alignment length to consider (default: 1000)')
parser.add_argument('-c', '--min_coverage', type=float, default=0.2, help='Minimum proportion of contig aligned to a reference (default: 0.2)')
parser.add_argument('-g', '--gap_threshold', type=int, default=10000, help='Max gap between alignments to merge (default: 10kb)')

args = parser.parse_args()

# === Parameters ===
min_length = args.min_length
min_coverage = args.min_coverage
gap_threshold = args.gap_threshold
paf_file = args.input
output_file = args.output

# === Parse PAF ===
def parse_paf(paf_path):
    """
    Parse PAF file into a dictionary: {contig: [alignment_segments]}
    """
    contig_alignments = defaultdict(list)
    with open(paf_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            target_name = fields[5]
            alignment_length = int(fields[10])

            if alignment_length >= min_length:
                contig_alignments[query_name].append({
                    'target_name': target_name,
                    'query_start': query_start,
                    'query_end': query_end,
                    'query_length': query_length,
                    'alignment_length': alignment_length
                })
    return contig_alignments

# === Merge Intervals ===
def merge_intervals(intervals, gap_threshold=10000):
    """
    Merge overlapping or nearby intervals based on gap_threshold.
    """
    if not intervals:
        return []

    sorted_intervals = sorted(intervals, key=lambda x: x['query_start'])
    merged = []
    current = sorted_intervals[0]

    for next_seg in sorted_intervals[1:]:
        if next_seg['query_start'] <= current['query_end'] + gap_threshold:
            current['query_end'] = max(current['query_end'], next_seg['query_end'])
            current['alignment_length'] += next_seg['alignment_length']
        else:
            merged.append(current)
            current = next_seg
    merged.append(current)
    return merged

# === Assign Names and Write Outputs ===
def name_contigs(contig_alignments, output_path, min_coverage, gap_threshold):
    """
    Write renamed contigs to output, including part info and original contig name.
    Also writes a second file listing only target_name and query_name pairs (no duplicates).
    """
    mapping_output = output_path.rsplit('.', 1)[0] + "_mapping.tsv"
    seen_pairs = set()

    with open(output_path, 'w') as out, open(mapping_output, 'w') as map_out:
        out.write("region_label\ttarget_name\tquery_start\tquery_end\toriginal_contig\n")
        map_out.write("target_name\tquery_name\n")

        for contig, segments in contig_alignments.items():
            if not segments:
                continue

            query_length = segments[0]['query_length']
            target_groups = defaultdict(list)

            for seg in segments:
                target_groups[seg['target_name']].append(seg)

            for target, group in target_groups.items():
                merged_regions = merge_intervals(group, gap_threshold)
                total_aligned = sum(region['alignment_length'] for region in merged_regions)
                coverage = total_aligned / query_length

                if coverage >= min_coverage:
                    # Write mapping file (no duplicates)
                    if (target, contig) not in seen_pairs:
                        map_out.write(f"{target}\t{contig}\n")
                        seen_pairs.add((target, contig))

                    # Write detailed regions
                    for i, region in enumerate(merged_regions):
                        region_label = f"{contig}_part{i+1}_{target}" if len(merged_regions) > 1 else f"{contig}_{target}"
                        out.write(f"{region_label}\t{target}\t{region['query_start']}\t{region['query_end']}\t{contig}\n")

# === Run ===
contig_alignments = parse_paf(paf_file)
name_contigs(contig_alignments, output_file, min_coverage, gap_threshold)
