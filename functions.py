from math import sqrt, acos, log, tan
import os

# Splits event into individual entries to be parsed
def splitline(line):
  return [x for x in line.split(' ') if x != '']

# Adds sublists up pointwise to return a "net" list.
# Used most often for getting net energy-momentum 4-vector
def sumzip(list):
  return map(sum, zip(*list))

# Get only the energy-momentum components from an event group 
def get_pe(event):
  return map(float, splitline(event)[3:7])

def get_pT(event):
  return norm(get_pe(event)[:2])

# Calculate invariant mass given an energy-momentum 4-vector
# m^2 = E^2 - p^2 = E^2 - (px^2 + py^2 + pz^2)
def invariant_mass(pe):
  return sqrt(pe[3]**2 - (pe[0]**2 + pe[1]**2 + pe[2]**2))

# Calculate pseudorapidity
def pseudorapidity(event):
  pe = get_pe(event)
  return -1.0*log(tan(0.5 * abs(acos(pe[2]/pe[3]))))

# For now, calculate RMS error from mass of top quark and W boson
def top_W_error(M3, M2):
  m_top = 172.9; m_W = 80.385
  return sqrt((M3 - m_top)**2 + (M2 - m_W)**2)

# calculates the dot product of two vectors
def dot(x,y):
  return None if len(x) != len(y) else sum(s*t for s,t in zip(x,y))

# Calculates 2-norm of a vector
def norm(x):
  return sqrt(sum(x_i**2 for x_i in x))

def angle(x,y):
  return abs(acos( dot(x,y) / (norm(x) * norm(y)) ))

def get_reco(filename):
  reco = '/Users/mukul/Desktop/'+ filename.split('/')[-1][:-4] + '-reco.txt'
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

def get_jets_from_file(filename):
  reco = filename[:-4] + '-reco.txt'
  ifnotreco = '/Users/mukul/Desktop/'+ filename.split('/')[-1][:-4] + '-reco.txt'

  if not os.path.isfile(reco):
    raw_input("Reco only file doesn't exist; creating it now in:\n" + ifnotreco + '\nPress any key to continue.')
    print 'Creating reco file...'
    get_reco(filename)
    print 'Done. Move file to same directory as input file and run program again.'
    exit(0)
  else:
    raw_input('Reco file exists. Press any key to proceed.')

  all_files = []; curr_file = []
  read_flag = False

  asdf = open('output2.txt', 'w')

  # Read file and parse Reco sections only
  for line in open(reco):
    if line.strip() == 'EndReco':
      read_flag = False
      all_files.append(curr_file)
      curr_file = []
    if read_flag:
      curr_file.append(line)
    if line.strip() == 'BeginReco':
      read_flag = True

  events_by_file = []

  asdf = open('output2.txt', 'w')
  for f in all_files:
    # READ RECOS AND GET 1 EVENT PER LINE
    events = []; curr_event = []
    for line in f:
      # ignore BeginReco/EndReco or line w/ number of events
      if len(line) <= 12: continue
      curr_event.append(line[:-1])
      if line[-2] in ['T', 'F']: # marks the end of an event
        events.append(' '.join(curr_event))
        curr_event = []

    # GET JETS AND BOTTOM QUARK EVENTS
    bottoms = []; jets = []
    for event in events:
      entries = splitline(event)
      if entries[2] == '4': # Jet
        # Separate bottom quarks from non-tagged jets
        bottoms.append(event) if entries[-6] == '5.' else jets.append(event)

    asdf.write(str(map(len, [events, bottoms, jets])) + '\n')

    events_by_file.append((events, bottoms, jets))
  asdf.close()

  return events_by_file