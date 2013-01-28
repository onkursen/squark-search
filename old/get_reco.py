import sys

flag = False

g = open('recos.txt', 'w')

for line in open(sys.argv[1]):
  if line.strip() == 'EndReco':
    flag = False
    print 'end reco'
    break
  if flag:
    g.write(line)
  if line.strip() == 'BeginReco':
    flag = True
    print 'begin reco'