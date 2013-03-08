import itertools, math, sys, os

# splits event into individual entries to be parsed
def splitline(line):
  return [x for x in line.split(' ') if x != '']

# Write the elements of an array to a file, one element per line
def write_file(filename, ary):
  f = open(filename, 'w')
  for m in ary:
    f.write(str(m))
    f.write('\n')
  f.close()

# --------------------------------------------------

all_files = []; curr_file = []
read_flag = False
N = 0

M3A = []; M2A = []; M3B = []; M2B = []
b_angles = []; j1_angles = []; j2_angles = []
MISSING_E = []
enough_bottoms = 0

# ------------------------------
# READ FILE AND PARSE RECOS ONLY
# ------------------------------
for line in open(sys.argv[1]):
  if line.strip() == 'EndReco':
    read_flag = False
    all_files.append(curr_file)
    curr_file = []
    N += 1
  if read_flag:
    curr_file.append(line)
  if line.strip() == 'BeginReco':
    read_flag = True

for n in range(len(all_files)):
  # -----------------------------------
  # READ RECOS AND GET 1 EVENT PER LINE
  # -----------------------------------
  lines = [l[1:-1] for l in all_files[n]]
  if lines == []: continue

  events = []; curr_event = []
  for line in lines:
    # ignore BeginReco/EndReco or line w/ number of events
    if len(line) <= 10: continue
    curr_event.append(line)
    if line[-1] in ['T', 'F']: # marks the end of an event
      events.append(''.join(curr_event))
      curr_event = []

  # -----------------------------------
  # GET JETS AND BOTTOM QUARK EVENTS
  # -----------------------------------
  energies = []
  for event in events:
    entries = splitline(event)
    if entries[2] == '4': # jet
      energies.append(entries[6])


print '%d energies out of %d events' % (len(energies), len(all_files))

write_file('energies.txt', energies)