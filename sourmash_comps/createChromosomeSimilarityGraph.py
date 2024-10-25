#This script reads the sourmash files and 
# creates a graph, that help to identify 
# same chromossomes in the different species
import gzip

# File of filenames
fofn = "sourmashfiles.fofn"
chrnames="chrnames.txt"	# File with the names of the chromossomes
id2name={}    # Dictionary to store the names of the chromossomes

with open(chrnames, 'r') as cn:
    for line in cn:
        line = line.strip()
        fields = line.split('\t')
        id2name[fields[0]] = {
            'chrname': fields[1],
            'chr': fields[2]
        }

with open(fofn, 'r') as f, open('chrgraph.txt', 'w') as cg:
    cg.write(f"Node1\tNode2\tSim\tNode1Chrname\tNode1Chr\tNode2Chrname\tNode2Chr\n")
    for line in f:
        sourmashfile = line.strip()
        # Read the sourmash file
        with gzip.open(sourmashfile, 'rt') as s:
            for line in s:
                if line.startswith('similarity,md5,'):
                    continue
                else:
                    line = line.strip().split(',')
                    chromossome = line[2].split('/')[-1].split('.')[0]
                    chromossome2= line[4].split('/')[-1].split('.')[0]
                    chr1=f'{chromossome}.1'
                    chr2=f'{chromossome2}.1'
                    cg.write(f"{chr1}\t{chr2}\t{line[0]}\t{id2name[chr1]['chrname']}\t{id2name[chr1]['chr']}\t{id2name[chr2]['chrname']}\t{id2name[chr2]['chr']}\n")
                    
