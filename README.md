# Sugarcane SP80-3280 Genome Assembly

## Sequencing data

### PacBio HiFi

### HiC

## Chromosome comparisons of parental species

We are comparing _S. officinarum_ LA-Purple with a basic number of chromosomes of 10,  each with 8 copies, so x=10, 2n=80, _S. spontaneum_ AP85-441 x=8, 2n=32, and _S. spontaneum_ NpX with x=10, 2n=40.

### kmer comparison with sourmash

We did a fast comparison of similarity among pairs of chromosomes of thes varieties, using [sourmarch](https://sourmash.readthedocs.io/en/latest/) with this [script](sourmash_comps/sourmash.sh). Results of the comparisons are in the [sourmarch_comps folder](sourmash_comps/). The sourmash results were processed with this [script](sourmash_comps/createChromosomeSimilarityGraph.py), and can be visualized in the following figure ([SVG](chrgraphs/chromosomesGraphSaccharum.svg), [Cytoscape](chrgraphs/chromosomesGraphSaccharum.sys)):

![Chromosome similarity graph](chrgraphs/chromosomesGraphSaccharum.png)

### Chr alignments with minimap2
