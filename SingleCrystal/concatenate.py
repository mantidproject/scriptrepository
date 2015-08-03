# Combine integrate files.

import os

expname = 'sapphire'
filename = expname + '.integrate'
output = open(filename, 'w')

# Create a list of the run numbers
run_list = ['1001', '1002', '1003', '1006']

# Read and write the first integrate file with header.
run = run_list[0]
filename = expname + '_' + run + '.integrate'
input = open(filename, 'r')
file_all_lines = input.read()
output.write(file_all_lines)
input.close()
os.remove(filename)

# Read and write the rest of the integrate files without the header.
for run in run_list[1:]:
    filename = expname + '_' + run + '.integrate'
    input = open(filename, 'r')
    for line in input:
        if line[0] == '0': break
    output.write(line)
    for line in input:
        output.write(line)
    input.close()
    os.remove(filename)
