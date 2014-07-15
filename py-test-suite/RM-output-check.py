import sys

if sys.argv < 2:
    print "must supply RepeatMasker *.fa.out file(s) to be processed"
    sys.exit(1)

s = set()
numLines = 0
numShortLines = 0

for filename in sys.argv[1:]:
    with open(filename, "r") as f:
        _ = [f.readline() for i in range(3)]
        for line in f:
            numLines += 1
            fields = line.split()
            if len(fields) > 14:
                s.add(int(fields[14]))
            else:
                numShortLines += 1
                print "SHORT LINE", ("(" + str(numLines) + "):"), line

print "numLines:", numLines
print "numShortLines:", numShortLines
print "len(s):", len(s)
print "max(s):", max(s)
