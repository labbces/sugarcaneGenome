#!/bin/bash

#SBATCH --job-name=SP803280_evalannot
#SBATCH --partition=long
#SBATCH --ntasks-per-node=20
#SBATCH --mem=150gb
#SBATCH --error=SP803280_evalannot-%j.err
#SBATCH --output=SP803280_evalannot-%j.out

export THREADS=20
export LINEAGE=liliopsida
export LIBRARY_PATH=/home/dmpachon/SugarcaneGenome_annotation/eval_annotations/mb_downloads/
export GFFREAD_BIN=/home/dmpachon/SugarcaneGenome_annotation/eval_annotations/software/gffread-0.12.7.Linux_x86_64/gffread
export COMPLEASM_BIN=/home/dmpachon/SugarcaneGenome_annotation/eval_annotations/software/compleasm_kit/compleasm.py
export AGAT_SIF=/home/dmpachon/SugarcaneGenome_annotation/eval_annotations/software/images/agat_1.6.1--pl5321hdfd78af_1.sif

./evaluate_annotations.sh evaluate_annotations.tsv annotations_qc
