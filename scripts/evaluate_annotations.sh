#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./evaluate_annotations.sh inputs.tsv [out_dir]
#
# inputs.tsv: 2 columns, tab-separated (or space-separated):
#   col1 = GFF/GFF3 path
#   col2 = genome FASTA path (can be .fa/.fasta/... and can be .gz)
#
# Example:
#   /data/a.gff3    /data/genome.fa.gz

INPUT_TSV="${1:?ERROR: provide inputs.tsv with 3 columns: gff<TAB>fasta<TAB>tool}"
OUT_DIR="${2:-annotation_qc}"

THREADS="${THREADS:-16}"
LINEAGE="${LINEAGE:-eukaryota}"
LIBRARY_PATH="${LIBRARY_PATH:-}"

GFFREAD_BIN="${GFFREAD_BIN:-gffread}"
COMPLEASM_BIN="${COMPLEASM_BIN:-compleasm}"
SINGULARITY_BIN="${SINGULARITY_BIN:-singularity}"
ISOFORMTABLE_BIN="${ISOFORMTABLE_BIN:-/home/dmpachon/SugarcaneGenome_annotation/eval_annotations/software/get_isoforms_table_from_gff.py}"
OMARK_DB="${OMARK_DB:-software/omark/LUCA.h5}"
AGAT_SIF="${AGAT_SIF:-images/agat_1.6.1--pl5321hdfd78af_1.sif}"
TAXID="${TAXID:-4546}"

command -v "$GFFREAD_BIN" >/dev/null 2>&1 || { echo "ERROR: gffread not found ($GFFREAD_BIN)"; exit 2; }
command -v "$COMPLEASM_BIN" >/dev/null 2>&1 || { echo "ERROR: compleasm not found ($COMPLEASM_BIN)"; exit 2; }
command -v "$SINGULARITY_BIN" >/dev/null 2>&1 || { echo "ERROR: singularity not found ($SINGULARITY_BIN)"; exit 2; }
command -v "$ISOFORMTABLE_BIN" >/dev/null 2>&1 || { echo "ERROR: get_isoforms_table_from_gff.py not found ($ISOFORMTABLE_BIN)"; exit 2; }

mkdir -p "$OUT_DIR"

echo "Input table: $INPUT_TSV"
echo "Output dir : $OUT_DIR"
echo "Threads    : $THREADS"
echo "Lineage    : $LINEAGE"
echo "TaxID4omark: $TAXID"

[[ -n "$LIBRARY_PATH" ]] && echo "Library    : $LIBRARY_PATH"

# Robust line reader:
# - skips blank lines
# - allows comments starting with '#'
# - supports TAB or spaces as separator
line_no=0
while IFS=$'\t' read -r gff genome tool rest || [[ -n "${gff:-}" ]]; do
  ((line_no++)) || true

  # If the file is space-separated, the read above will put whole line in $gff.
  # Try splitting on whitespace in that case.
  if [[ -n "${gff:-}" && -z "${genome:-}" ]]; then
    # shellcheck disable=SC2206
    parts=( $gff )
    gff="${parts[0]:-}"
    genome="${parts[1]:-}"
    tool="${parts[2]:-}"
  fi

  # skip blanks / comments
  [[ -z "${gff:-}" ]] && continue
  [[ "${gff:0:1}" == "#" ]] && continue

  if [[ -z "${genome:-}" ]]; then
    echo "WARNING: line $line_no has no genome FASTA (need 2 columns). Skipping." >&2
    continue
  fi

  if [[ -z "${tool:-}" ]]; then
     echo "WARNING: line $line_no has no tool info (need 3 columns). Skipping." >&2
     continue
  fi

  if [[ ! -s "$gff" ]]; then
    echo "WARNING: GFF not found/empty: $gff (line $line_no). Skipping." >&2
    continue
  fi
  if [[ ! -s "$genome" ]]; then
    echo "WARNING: genome FASTA not found/empty: $genome (line $line_no). Skipping." >&2
    continue
  fi

  # Sample name for output folders: basename of gff without extension
  gff_bn="$(basename "$gff")"
  base="${gff_bn%.*}"
  annot_file_type="${gff_bn##*.}"

  echo
  echo "======================================================"
  echo "Line $line_no"
  echo "GFF   : $gff"
  echo "Genome: $genome"
  echo "Sample: $base"
  echo "Tool  : $tool"

  sample_out="$OUT_DIR/$tool/$base"
  seq_out="$sample_out/sequences"
  comp_out="$sample_out/compleasm"
  cleangff_out="$sample_out/cleanGFF/"
  statsgff_out="$sample_out/statsGFF/"
  omark_out="$sample_out/omark/"

  mkdir -p "$seq_out" "$comp_out" "$cleangff_out" "$statsgff_out" "$omark_out"

  transcripts_fa="$seq_out/${base}.transcripts.fa"
  cds_fa="$seq_out/${base}.cds.fa"
  proteins_fa="$seq_out/${base}.proteins.fa"
  cleangff_file="$cleangff_out/$base.AGAT.clean.gff3"
  cleangff_log="$cleangff_out/$base.AGAT.clean.log"
  statsgff_file="$statsgff_out/$base.AGAT.stats.txt"
  statsgff_log="$statsgff_out/$base.AGAT.stats.log"
  isoform_table="$omark_out/isoform.table"
  omamer_file="$omark_out/$base.omamer"
  omamer_log="$omark_out/$base.omamer.log"
  omark_log="$omark_out/$base.omark.log"

  # --------------------------------------------------
  # Prepare genome FASTA for gffread (must be a real file)
  # --------------------------------------------------
  tmp_genome=""
  genome_to_use="$genome"

  if [[ "$genome" == *.gz ]]; then
    # Cache decompressed genome under the sample output so reruns don't re-decompress
    tmp_genome="$seq_out/${base}.genome.decompressed.fa"

    if [[ -s "$tmp_genome" ]]; then
      echo "✓ Decompressed genome cache exists → $tmp_genome"
    else
      echo "→ Decompressing genome (.gz) to: $tmp_genome"
      gunzip -c "$genome" > "$tmp_genome"
    fi

    genome_to_use="$tmp_genome"
  fi

  # --------------------------------------------------
  # Step 0: normalize GFF and compute stats using AGAT
  # --------------------------------------------------
  if [[ -s "$cleangff_file" ]]; then
    echo "✓ Normalized GFF already exist → skipping agat_convert_sp_gxf2gxf.pl"
  else
    if [[ "$annot_file_type" == "gff" || "$annot_file_type" == "gtf" ]]; then
      echo "→ Running agat_convert_sp_gxf2gxf.pl"
      "$SINGULARITY_BIN" exec "$AGAT_SIF" agat_convert_sp_gxf2gxf.pl \
        --"$annot_file_type" "$gff" \
        -o "$cleangff_file" \
        --cpu "$THREADS" > "$cleangff_log" 2>&1
    else
       echo "Extension must be gff or gtt → skipping file."
       continue
    fi
  fi

  if [[ -s "$statsgff_file" ]]; then
    echo "✓ Stats for GFF already exist → skipping agat_sp_statistics.pl"
  else
    echo "→ Running agat_sp_statistics.pl"
    "$SINGULARITY_BIN" exec "$AGAT_SIF" agat_sp_statistics.pl \
      --gff "$cleangff_file" \
      -g "$genome_to_use" -d \
      -o "$statsgff_file" \
      --cpu "$THREADS" > "$statsgff_log" 2>&1
  fi

  # --------------------------------------------------
  # Step 1: gffread (only if proteins missing)
  # --------------------------------------------------
  if [[ -s "$proteins_fa" ]]; then
    echo "✓ Proteins already exist → skipping gffread"
  else
    echo "→ Running gffread"
    "$GFFREAD_BIN" \
      -g "$genome_to_use" \
      -w "$transcripts_fa" \
      -x "$cds_fa" \
      -y "$proteins_fa" \
      "$cleangff_file"

    if [[ ! -s "$proteins_fa" ]]; then
      echo "WARNING: Protein file empty after gffread (line $line_no). Skipping compleasm." >&2
      continue
    fi
  fi

  # --------------------------------------------------
  # Step 2: compleasm (skip if summary exists)
  # --------------------------------------------------
  summary_file="$comp_out/summary.txt"
  if [[ -s "$summary_file" ]]; then
    echo "✓ Compleasm output exists ($summary_file) → skipping"
  else
   echo "→ Running compleasm (protein mode)"
   if [[ -n "$LIBRARY_PATH" ]]; then
    "$COMPLEASM_BIN" protein \
      -p "$proteins_fa" \
      -l "$LINEAGE" \
      -o "$comp_out" \
      -t "$THREADS" \
      -L "$LIBRARY_PATH"
   else
    "$COMPLEASM_BIN" protein \
      -p "$proteins_fa" \
      -l "$LINEAGE" \
      -o "$comp_out" \
      -t "$THREADS"
   fi
  fi


  # --------------------------------------------------
  # Step 3: omark
  # --------------------------------------------------
  source ~/miniconda3/etc/profile.d/conda.sh
  conda activate omark
  if [[ -s "$isoform_table" ]]; then
    echo "✓ isoform_table exists ($summary_file) → skipping"
  else
    echo "→ Running get_isoforms_table_from_gff.py"
    python3 "$ISOFORMTABLE_BIN" "$cleangff_file" "$isoform_table"
  fi
  
  if  [[ -s "$omamer_file" ]]; then
    echo "✓ omamer_file exists ($omamer_file) → skipping"
  else
    echo "→ Running omamer"
    omamer search --db "$OMARK_DB" --query "$proteins_fa" --out "$omamer_file" --nthreads "$THREADS" > "$omamer_log" 2>&1
  fi

  if [[ -s "$omark_log" ]]; then
    echo "✓ omark_log exists ($omark_log) → skipping"
  else
    echo "→ Running omark"
    omark -f "$omamer_file" -d "$OMARK_DB" -o "$omark_out" --isoform_file "$isoform_table" --taxid "$TAXID" > "$omark_log" 2>&1
  fi  
  conda deactivate

  echo "✓ Done: $base"

done < "$INPUT_TSV"

echo
echo "All done."
