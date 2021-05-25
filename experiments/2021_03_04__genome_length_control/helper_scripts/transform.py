# External expectations: 
#   - A list of all files in *one* of the size directories (e.g., mcsize__16)
#       - Just us an 'ls mcsize__16 > old_list.txt'
#       - Expected to be saved as "old_list.txt", otherwise change first with statement below

# Length of genome in orginal data
old_length = 400 
new_length = 100 # Change me for each genome size!
# Organism sizes for which the script should prep timing data
size_list = [8, 16, 32, 64, 128, 256, 512]

# Automatically calcluated variables
old_threshold = int(old_length) * 0.6
new_threshold = int(new_length) * 0.6
min_val = old_threshold - new_threshold
max_val = old_threshold + (new_length - new_threshold)

# All removes must come before all moves
rm_list = []
mv_list = []

with open('old_list.txt', 'r') as in_fp:
    with open('transfer.sh', 'w') as out_fp:
        # Each line should be in the form X.dat where X in [0, old_length]
        for line in in_fp:
            line = line.strip()
            if line == '':
                continue
            # Extract the number and calculate what it should move to in order to 
                # preserve the orginial restraint buffer value
                # e.g., 270 in a 400-bit genome has a restraint buffer value of 10 (270 - 260 = 10)
                    # Thus, for a 200-bit genome it would move to 120 + 10 = 130 
            val = int(line.split('.')[0])
            old_diff = val - old_threshold
            new_val = old_diff + new_threshold
            # If it's within the range of the new genome length, add it
            if new_val >= 0 and new_val <= new_length:
                mv_list.append([val, new_val])
            else: # Else, remove it
                rm_list.append(val)

        # Create a bash script to move and delete all files as needed
        out_fp.write('#!/bin/bash\n')
        for size in size_list:
            out_fp.write('cd mcsize__' + str(size) + '\n')
            out_fp.write('echo "' + str(size) + '"\n')
            for rm_val in rm_list:
                out_fp.write('rm ' + str(rm_val) + '.dat \n')
            for mv_pair in mv_list:
                out_fp.write('mv ' + str(mv_pair[0]) + '.dat ' + str(mv_pair[1]) + '.dat\n')
            out_fp.write('cd ..\n')
        

