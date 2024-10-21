import sys
import gzip

# File of filenames
fofn = "sourmashfiles.fofn"
# Dictionary to store groups
groups = {}
groupcount = 1
# Threshold for overlap percentage (e.g., 50%)
overlap_threshold = 0.5

with open(fofn, 'r') as f, open('chrgraph.txt', 'w') as cg:
    for line in f:
        sourmashfile = line.strip()
        listchrs = []
        
        # Read the sourmash file
        with gzip.open(sourmashfile, 'rt') as s:
            for line in s:
                if line.startswith('similarity,md5,'):
                    continue
                else:
                    line = line.strip().split(',')
                    chromossome = line[2].split('/')[-1].split('.')[0]
                    chromossome2= line[4].split('/')[-1].split('.')[0]
                    cg.write(f'{chromossome}.1\t{chromossome2}.1\t{line[0]}\n')
                    listchrs.append(chromossome)
        
        # Check if the group already exists and calculate overlap
        group_found = False
        for k, v in groups.items():
            set_v = set(v)
            set_listchrs = set(listchrs)
            overlap = set_v.intersection(set_listchrs)
            overlap_percentage = len(overlap) / min(len(set_v), len(set_listchrs))
            
            # If overlap percentage is higher than the threshold
            if overlap_percentage >= overlap_threshold:
                # print(f'Overlap found between {sourmashfile} and group {k}')
                group_found = True
                
                # Update the group to keep all unique elements (not just the larger list)
                combined_list = list(set_v.union(set_listchrs))
                groups[k] = combined_list
                # print(f'Group {k} updated with merged data from {sourmashfile}')
                break
        
        # If no significant overlap, create a new group
        if not group_found:
            print(f'New group created for {sourmashfile}')
            groups[groupcount] = listchrs
            groupcount += 1

# Print the groups to a file
with open('groups.txt', 'w') as g:
    for k, v in groups.items():
        for chr in v:
            sp="_".join(chr.split('_')[0:2])
            g.write(f'{k}\t{sp}\t{chr.upper()}.1\n')
