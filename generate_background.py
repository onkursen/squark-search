from itertools import count
from os import system
from time import time
from sys import argv

corrupted = \
{
    0: [323, 326, 329, 399, 438, 460, 463, 472, 475, 495, 514, 517, 520, 521, 527, 531],
    1: [828, 831],
    2: [25, 50, 70, 96, 128, 169, 177, 196, 245, 261, 340, 398, 442, 474, 497, 512, 588],
    3: [9, 29, 37, 48, 50, 55, 60, 61, 84, 124, 132, 137, 223, 224, 276, 316, 335, 340, 364, 367, 374, 382],
    4: [17, 36, 56, 58, 67, 74, 76, 78, 82, 84, 97, 106, 109, 113, 123, 125, 135, 143, 150],
    5: [15, 34, 35, 48, 60, 64, 75, 76, 93]
}

limits = \
{
    0: 616,
    1: 794,
    2: 555,
    3: 357,
    4: 228,
    5: 77
}

allowed6 = [289, 105, 380, 184, 132, 296, 128,  85, 297, 279, 345, 358, 278, 373, 353,  66, 399, 308, 225, 151, 339, 237, 354, 387, 302, 338,  20, 166, 314, 171, 396,  58, 290,  21, 390,  13, 352, 340, 369, 341, 426, 381, 365, 718, 398,  78, 318, 337,  50, 291, 377, 382, 452,  36, 465, 563,  63, 180, 335, 323, 388, 776, 731, 375, 414, 348, 137, 389, 416,  72,  45, 590, 505, 359, 572, 294, 356, 165, 555, 758, 714, 133, 313,  82, 117,  99,   9, 778, 363, 432, 245, 529, 795, 506, 611, 524, 315, 234, 331, 347, 334, 121, 409, 336, 325, 317, 209, 497, 312, 601, 343, 753, 322, 516,  70, 800, 779, 304, 767, 618, 503, 675, 298, 642, 593]

def recofy(filename, f):
    t = time()
    all_files = []; curr_file = []
    read_flag = False
    # Read file and parse Reco sections only

    for line in open(filename):
      line = line.strip()
      if line == 'EndReco':
        read_flag = False
        curr_file.append(line + '\n')
        f.write('\n'.join(curr_file))
        curr_file = []
      if read_flag:
        curr_file.append(line)
      if line == 'BeginReco':
        read_flag = True
        curr_file.append(line)
    print 'Done. Took %.3f secs.' % (time()-t)

directory =  "/x1/data/nkolev/project10/output/ttjPP/"
fraction = 0.5

output = open('ttj20.txt', 'w')
for num_jets in range(0,6):
    t0 = time()
    print 'Constructing for %d jets' % num_jets
    print 'Limit is %d' % limits[num_jets]
    num_files = 0
    locs = []
    for f in count(1):
        if f in corrupted[num_jets]: continue
        filename = '%sttj00%df%d.dat' % (directory, num_jets, f)
        recofy(filename, output)
        num_files += 1
        print num_files
        if num_files == limits[num_jets]*fraction:
            break
    print 'Total time for %d jets was %.3f secs.' % (num_jets, time()-t0)
output.close()