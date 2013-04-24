from time import time
from sys import argv
t = time()
filename = argv[1]
reco = '/home/onkursen/'+ filename.split('/')[-1][:-4] + '-reco.txt'
f = open(reco,'w')
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
f.close()
print 'Done. Took %.3f secs.' % (time()-t)
