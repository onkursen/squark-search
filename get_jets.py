import itertools, math, sys, os

# ------------------------------
# READ FILE AND PARSE RECOS ONLY
# ------------------------------

flag = False
N = 0
# filename = 'ttj006f393.dat'

for line in open(sys.argv[1]):
  if line.strip() == 'EndReco':
    flag = False
    # print 'end reco', N
    g.close()
    N += 1
  if flag:
    g.write(line)
  if line.strip() == 'BeginReco':
    flag = True
    g = open("recos%02d.txt" % N, 'w')
    # print 'begin reco'


M3A = []
M2A = []
M3B = []

enough = 0

for n in range(N):
  # -----------------------------------
  # READ RECOS AND GET 1 EVENT PER LINE
  # -----------------------------------
  lines = [l[1:-1] for l in open('recos%02d.txt' % n).readlines()]

  if lines == []: continue

  curr_event = []
  events = []
  # want to ignore BeginReco/EndReco or line w/ number of events
  for line in lines:
    if len(line) > 10:
      curr_event.append(line)

      if line[-1] in ['T', 'F']: # marks the end of an event
        events.append(''.join(curr_event))
        curr_event = []

  # -----------------------------------
  # GET JETS AND BOTTOM QUARK EVENTS
  # -----------------------------------

  # splits event into individual entries to be parsed
  def splitline(line):
    return [x for x in line.split(' ') if x != '']

  bottoms = []
  jets = []

  for event in events:
    entries = splitline(event)
    if entries[2] == '4': # jet
      # keep bottom quarks separate from other jets
      if entries[13] == '5.':
        bottoms.append(event)
      else:
        jets.append(event)

  if len(bottoms) < 2:
    print 'Event %d has too few bottom quarks' % n
    continue
  else:
    enough += 1

  # ------------------------------------------
  # SELECT TWO bjj GROUPS WITH HIGHEST P_T
  # ------------------------------------------

  # Adds sublists up pointwise to return a "net" list.
  # Used most often for getting net energy-momentum 4-vector
  def sumzip(list):
    return map(sum, zip(*list))

  # keep track of bjj groups with top two p_t
  first = [0, '', '', '']
  second = [0, '', '', '']

  jet_combos = set(itertools.combinations(jets, 2)) # all pairs of jets

  for b in bottoms:
    bsplit = [x for x in b.split(' ') if x != '']
    for (j1, j2) in jet_combos:

      # p_x and p_y for both events
      p1 = map(float, splitline(j1)[3:5])
      p2 = map(float, splitline(j2)[3:5])

      p = map(sum, zip(p1, p2)) # net x,y-momentum
      p_t = math.sqrt(p[0]**2 + p[1]**2) # transverse momentum p_t

      # want to keep groups with two highest p_t
      if p_t > first[0]:
        second = first
        first = [p_t, b, j1, j2]
      elif p_t > second[0]:
        second = [p_t, b, j1, j2]

  # ----------------------------
  # CALCULATE M3 FOR BOTH GROUPS
  # ----------------------------

  # Get only the energy-momentum components from an event group 
  def get_pe(event):
    return map(float, splitline(event)[3:7])

  # Calculate invariant mass given an energy-momentum 4-vector
  def invariant_mass(pe):
    return math.sqrt(pe[3]**2 - (pe[0]**2 + pe[1]**2 + pe[2]**2))

  # individual energy-momentum 4-vectors
  pe1 = map(get_pe, first[1:])
  pe2 = map(get_pe, second[1:])

  pe_net1, pe_net2 = sumzip(pe1), sumzip(pe2) # net energy-momentum of both groups

  if pe_net1 == []: continue

  M3_1, M3_2 = invariant_mass(pe_net1), invariant_mass(pe_net2) # invariant mass

  # print 'M3', M3_1, M3_2

  # ------------------------------------
  # CALCULATE M2 FOR JETS IN BOTH GROUPS
  # ------------------------------------

  pe_jets1, pe_jets2 = sumzip(pe1[1:]), sumzip(pe2[1:]) # net energy-momentum for jets
  M2_1, M2_2 = invariant_mass(pe_jets1), invariant_mass(pe_jets2) # invariant mass

  # print 'M2', M2_1, M2_2

  # -------------------------------------------------------
  # CALCULATE ERROR FROM b, W MASSES AND ASSIGN TOP QUARK A
  # COLLECT M3 AND M2 CORRESPONDING TO TOP QUARK A ACROSS MULTIPLE EVENTS
  # -------------------------------------------------------

  m_top = 172.9
  m_W = 80.385

  # For now, calculate RMS error
  err1 = math.sqrt((M3_1 - m_top)**2 + (M2_1 - m_W)**2)
  err2 = math.sqrt((M3_2 - m_top)**2 + (M2_2 - m_W)**2)

  # Assign top system to minimize error
  # topA = first if err1 <= err2 else second
  topA = first 
  if err1 > err2:
    topA = second
  # print 'Top quark A:'
  # for j in range(1,4): print topA[j]

  if err1 <= err2:
    M3A.append(M3_1)
    M2A.append(M2_1)
  else:
    M3A.append(M3_2)
    M2A.append(M2_2)

  topB = [x for x in events if x not in topA[1:]]
  peB = map(get_pe, topB)
  pe_netB = sumzip(peB)
  M3B.append(invariant_mass(pe_netB))

# ------------------------------------------------------
# OUTPUT ALL M3 AND M2 VALUES TO FILES FOR ROOT PLOTTING
# ------------------------------------------------------

f = open('m3a.txt', 'w')
for m in M3A:
  f.write(str(m))
  f.write('\n')
f.close()

f = open('m3b.txt', 'w')
for m in M3B:
  f.write(str(m))
  f.write('\n')
f.close()

f = open('m2a.txt', 'w')
for m in M2A:
  f.write(str(m))
  f.write('\n')
f.close()

print enough, ' events had enough bottom quarks'

# clear temporary files
os.system('rm recos*.txt')