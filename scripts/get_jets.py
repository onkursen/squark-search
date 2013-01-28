import itertools, math, sys, os

# splits event into individual entries to be parsed
def splitline(line):
  return [x for x in line.split(' ') if x != '']

# Adds sublists up pointwise to return a "net" list.
# Used most often for getting net energy-momentum 4-vector
def sumzip(list):
  return map(sum, zip(*list))

# Get only the energy-momentum components from an event group 
def get_pe(event):
  return map(float, splitline(event)[3:7])

# Calculate invariant mass given an energy-momentum 4-vector
def invariant_mass(pe):
  return math.sqrt(pe[3]**2 - (pe[0]**2 + pe[1]**2 + pe[2]**2))

# For now, calculate RMS error from mass of top quark and W boson
def top_W_error(M3, M2):
  m_top = 172.9; m_W = 80.385
  return math.sqrt((M3 - m_top)**2 + (M2 - m_W)**2)

# Write the elements of an array to a file, one element per line
def write_file(filename, ary):
  f = open(filename, 'w')
  for m in ary:
    f.write(str(m))
    f.write('\n')
  f.close()

# calculates the dot product of two vectors
def dot(x,y):
  if len(x) != len(y):
    return 
  return sum(x_i*y_i for x_i,y_i in zip(x,y))

# calculates 2-norm of a vector
def norm(x):
  return math.sqrt(sum(x_i**2 for x_i in x))

def angle(x,y):
  return abs(math.acos( dot(x,y) / (norm(x) * norm(y)) ))

def get_best_bjj(bottoms, jets, top_system):
  # ------------------------------------------
  # SELECT TWO bjj GROUPS WITH HIGHEST P_T
  # ------------------------------------------
  # keep track of bjj groups with top two p_t
  first = [0, '', '', '']
  second = [0, '', '', '']

  jet_combos = set(itertools.combinations(jets, 2)) # all pairs of jets

  for b in bottoms:
    bsplit = [x for x in b.split(' ') if x != '']
    for (j1, j2) in jet_combos:

      # p_x and p_y for both events
      p1, p2 = map(float, splitline(j1)[3:5]), map(float, splitline(j2)[3:5])

      p = map(sum, zip(p1, p2)) # net x,y-momentum
      p_t = math.sqrt(p[0]**2 + p[1]**2) # transverse momentum p_t

      # want to keep groups with two highest p_t
      if p_t > first[0]:
        second = first
        first = [p_t, b, j1, j2]
      elif p_t > second[0]:
        second = [p_t, b, j1, j2]

  # limiting case: if only 2 jets, then second will never be assigned
  if second[0] == 0: second = first

  # ----------------------------
  # CALCULATE M3 FOR BOTH GROUPS
  # ----------------------------
  pe1, pe2 = map(get_pe, first[1:]), map(get_pe, second[1:]) # individual energy-momentum 4-vectors
  pe_net1, pe_net2 = sumzip(pe1), sumzip(pe2) # net energy-momentum of both groups

  if pe_net1 == []: return
  M3_1, M3_2 = invariant_mass(pe_net1), invariant_mass(pe_net2) # invariant mass

  # ------------------------------------
  # CALCULATE M2 FOR JETS IN BOTH GROUPS
  # ------------------------------------
  pe_jets1, pe_jets2 = sumzip(pe1[1:]), sumzip(pe2[1:]) # net energy-momentum for jets
  M2_1, M2_2 = invariant_mass(pe_jets1), invariant_mass(pe_jets2) # invariant mass

  # -------------------------------------------------------
  # CALCULATE ERROR FROM b, W MASSES AND ASSIGN BEST TOP QUARK
  # COLLECT M3 AND M2 CORRESPONDING TO TOP QUARK A ACROSS MULTIPLE EVENTS
  # -------------------------------------------------------
  err1, err2 = top_W_error(M3_1, M2_1), top_W_error(M3_2, M2_2)
  current_system = [M3A,M2A] if top_system == 'A' else [M3B,M2B]
  if err1 <= err2:
    best_top = first 
    current_system[0].append(M3_1)
    current_system[1].append(M2_1)
  else:
    best_top = second
    current_system[0].append(M3_2)
    current_system[1].append(M2_2)

  return best_top

all_files = []; curr_file = []
read_flag = False
N = 0

M3A = []; M2A = []; M3B = []; M2B = []
b_angles = []; j1_angles = []; j2_angles = []
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
  bottoms = []; jets = []
  for event in events:
    entries = splitline(event)
    if entries[2] == '4': # jet
      # separate bottom quarks from non-tagged jets
      bottoms.append(event) if entries[13] == '5.' else jets.append(event)

  if len(bottoms) < 2:
    print 'Event %d has too few bottom quarks' % n
    continue
  else: enough_bottoms += 1

  # -----------------------------------
  # BEST BJJ COMBO = TOP QUARK A
  # BEST BJJ COMBO FROM REMAINING JETS = TOP QUARK B
  # -----------------------------------
  topA = get_best_bjj(bottoms, jets, 'A');
  if not topA: continue

  # Remove bjj of A from the set of bottoms and jets to consider for system B
  bottoms2 = [x for x in bottoms if x != topA[1]]
  jets2 = [x for x in jets if x not in topA[2:]]
  topB = get_best_bjj(bottoms2, jets2, 'B')
  if not topB: continue

  # -----------------------------------
  # AZIMUTHAL ANGLE CUTS
  # -----------------------------------

  # Sum transverse momentum components of ALL events in collision
  # Theoretically should be 0, but it won't be 
  pT_file_total = sumzip([get_pe(event) for event in events])[:2]

  # "missing" transverse momentum = negative of sum
  pT_missing = [-1*x for x in pT_file_total]

  # transverse momentum vector of each in bjj of system B
  pT_bjjB = [get_pe(topB[k])[:2] for k in range(1, 4)]

  # azimuthal angle in between each jet and missing transverse momentum
  angles = [angle(j, pT_missing) for j in pT_bjjB]

  b_angles.append(angles[0])
  j1_angles.append(angles[1])
  j2_angles.append(angles[2])
print '%d of %d events had enough bottom quarks' % (enough_bottoms, len(all_files))

# ------------------------------------------------------
# OUTPUT VALUES TO FILES FOR ROOT PLOTTING:
# 1. M3 OF TOP QUARK A
# 2. M2 OF TOP QUARK A
# 3. M3 OF TOP QUARK B
# 4. M2 OF TOP QUARK B
# 5. AZIMUTHAL ANGLES FOR BOTTOM QUARK B
# 6. AZIMUTHAL ANGLES FOR JET 1 IN SYSTEM B
# 7. AZIMUTHAL ANGLES FOR JET 2 IN SYSTEM B
# ------------------------------------------------------
write_file('output/m3a.txt', M3A)
write_file('output/m3b.txt', M3B)
write_file('output/m2a.txt', M2A)
write_file('output/m2b.txt', M2B)
write_file('output/angles_b.txt', b_angles)
write_file('output/angles_j1.txt', j1_angles)
write_file('output/angles_j2.txt', j2_angles)

# combine plots into one tar file
# os.system("tar cf plots.tar *.txt && rm *.txt")