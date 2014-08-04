# Takes paths of Mason read simulator output files are formats them for easier reading.

import sys, os

if len(sys.argv) < 2:
    print "must supply path(s) of file(s) to be processed"
    os.Exit(1)

for filepath in sys.argv[1:]:
    with open(filepath, "r") as infile, open(filepath + ".proc", "w") as outfile:

        for line in infile:
            if line[0] != '>':
                outfile.write(line.strip().lower() + '\n')
