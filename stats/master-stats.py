import sys
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
from collections import defaultdict

dir_name = sys.argv[1] + "-stats/"

matches = list()

with open(dir_name + "matches.txt", "r") as match_file:
    [match_file.readline() for _ in range(3)]

    for line in match_file:
        fields = line.split()
        if len(fields) == 5:
            fields[0], fields[1], fields[2], fields[4] = int(fields[0]), int(fields[1]), int(fields[2]), int(fields[4])
            matches.append(fields)
        else:
            print "malformed matches line"
            quit()

lens = defaultdict(int)

for match in matches:
    lens[match[2]] += 1

with open("match-lens.txt", "w") as outfile:
    for m_len, cnt in lens.iteritems():
        outfile.write(str(m_len) + " " + str(cnt) + '\n')

#plt.hist(lens_list)
#plt.savefig("match-lens.svg")
'''

depths = defaultdict(int)

with open(dir_name + "class-nodes.txt", "r") as cn_file:
    [cn_file.readline() for _ in range(2)]

    for line in cn_file:
        fields = line.split()
        depths[int(fields[2])] += 1

with open("cn-depths.txt", "w") as outfile:
    for depth, cnt in depths.iteritems():
        outfile.write(str(depth) + " " + str(cnt) + "\n")
'''
