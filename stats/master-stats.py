import sys
import matplotlib.pyplot as plt

dir_name = "../" + sys.argv[1] + "-reads/"

repeats = list()

with open(dir_name + "repeats.txt", "r") as repeat_file:
    repeat_file.readline()    # drop header

    for line in repeat_file:
        fields = line.split()
        if len(fields) == 4:
            fields[0], fields[2], fields[3] = int(fields[0]), int(fields[2]), int(fields[3])
            repeats.append(fields)
        else:
            print "malformed repeat line"
            quit()


