import sys
import gzip

# File of filenames
fofn = "sourmashfiles.fofn"
# Dictionary to store groups
groups = {}
groupcount = 1
# Threshold for overlap percentage (e.g., 50%)
overlap_threshold = 0.5

with open(fofn, 'r') as f:
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
                print(f'Overlap found between {sourmashfile} and group {k}')
                group_found = True
                
                # Keep the larger list
                if len(set_listchrs) > len(set_v):
                    groups[k] = listchrs  # Replace the group with the larger list
                    print(f'Group {k} updated with {sourmashfile}')
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
            g.write(f'{k}\t{chr}\n')
