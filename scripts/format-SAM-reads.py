# Takes paths of Mason read simulator output files are formats them for easier reading.

import sys, os

if len(sys.argv) < 2:
    print "must supply path(s) of file(s) to be processed"
    os.Exit(1)

for filepath in sys.argv[1:]:
    with open(filepath, "r") as infile, open(filepath + ".proc", "w") as outfile:
        # drop header
        [infile.readline() for _ in range(3)]
        for line in infile:
            seq = line.split()[9]
            seq = seq.lower()
            outfile.write(seq + '\n')
