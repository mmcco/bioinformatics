import sys

if sys.argv < 2:
    print "must supply RepeatMasker *.fa.out file(s) to be processed"
    sys.exit(1)

repeatMap = dict()

for filename in sys.argv[1:]:
    with open(filename, "r") as f:
        _ = [f.readline() for i in range(3)]
        for line in f:
            fields = line.split()
            repName, repClass, repID = fields[9], fields[10], int(fields[14])
            if repID not in repeatMap:
                repeatMap[repID] = [(repName, repClass)]
            elif (repName, repClass) not in repeatMap[repID]:
                repeatMap[repID].append((repName, repClass))

if max(repeatMap.values()) > 1:
    for repID, names in repeatMap.iteritems():
        if len(names) > 1:
            print repID, names
