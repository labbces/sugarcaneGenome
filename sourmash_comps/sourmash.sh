#!/bin/bash

#$ -q all.q
#$ -cwd
module load sourmash
cd /data/dmrp/
#sourmash sketch dna -p scaled=1000,k=31 --outdir /data/dmrp/sourmarsh_sigs/ /data/dmrp/SOFF_LAPpurple_vs_self/*.fasta
#sourmash sketch dna -p scaled=1000,k=31 --outdir /data/dmrp/sourmarsh_sigs/ /data/dmrp/SSPO_NPX_GCA_vs_self/*.fasta
#sourmash sketch dna -p scaled=1000,k=31 --outdir /data/dmrp/sourmarsh_sigs/ /data/dmrp/SSPO_AP85441_vs_self/*.fasta
#sourmash index sugarcane_db /data/dmrp/sourmarsh_sigs/
#sourmash sketch dna -p scaled=1000,k=31 --outdir /data/dmrp/sourmarsh_sigs/ /data/dmrp/HFHC.p_utg.fna.gz 
mkdir -p /data/dmrp/sourmarsh_comps/
for sig in  /data/dmrp/sourmarsh_sigs/*.sig; do
 outfile=`basename $sig`
 outfile=${outfile/.sig/.sourmash.search.txt}
 sourmash search $sig /data/dmrp/sourmarsh_sigs/ -o /data/dmrp/sourmarsh_comps/$outfile -n 0
done
