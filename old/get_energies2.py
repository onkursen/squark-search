import sys 
all_files = []; curr_file = []
read_flag = False
N = 0

def splitline(line):
  return [x for x in line.split(' ') if x != '']

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
    if entries[2] == '4':
      energies.append(float(entries[6]))

g = open('energies.txt', 'w')
for e in energies:
  g.write(str(e))
  g.write('\n')

print energies
print max(energies)