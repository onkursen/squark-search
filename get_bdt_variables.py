# ------------------------------------------
# Onkur Sen
# 
# get_bdt_variables.py
#
# Usage: python get_bdt_variables.py [filename]
# 
# Reads all event runs from a file,
# separates each run into two bjj systems,
# and keeps track of important variables
# to feed into a boosted decision tree (BDT).
# ------------------------------------------

from functions import *
from itertools import combinations
from math import sqrt
from sys import argv
from time import time
import os

# Selects two bjj combinations with highest vectorially summed transverse momentum. Returns the combination with a smaller least-squares error from the top quark and W boson.
def get_best_bjj(bottoms, jets):
  
  if len(jets) < 2: return None, None, None 
  # ------------------------------------------
  # SELECT TWO bjj GROUPS WITH HIGHEST P_T
  # ------------------------------------------
  
  # keep track of bjj groups with top two p_t
  first = [0, '', '', '']; second = [0, '', '', ''] 
  jet_combos = set(combinations(jets, 2)) # all pairs of jets

  for b in bottoms:
    bsplit = [x for x in b.split(' ') if x != '']
    for (j1, j2) in jet_combos:

      # p_x and p_y for both events
      p1, p2 = map(float, splitline(j1)[3:5]), map(float, splitline(j2)[3:5])

      p = map(sum, zip(p1, p2))           # net x,y-momentum
      p_t = sqrt(p[0]**2 + p[1]**2)  # transverse momentum p_t

      # keep groups with two highest p_t
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
  if err1 <= err2:
    return first, M3_1, M2_1
  else:
    return second, M3_2, M2_2

def main():
  M3A = []; M2A = []; M3B = []; M2B = []
  ANGLES_B_A = []; ANGLES_J1_A = []; ANGLES_J2_A = [];
  ANGLES_B_B = []; ANGLES_J1_B = []; ANGLES_J2_B = [];
  MISSING_E = []
  enough_bottoms = 0; has_topA = 0; has_topB = 0; meets_cuts = 0
  events_with_params = []

  print 'Reading events from source file %s' % argv[1]
  t = time()
  events_by_file = get_jets_from_file(argv[1])
  num_events = len(events_by_file)
  print 'Done. Took %f secs.\n' % (time()-t)
  
  print 'Iterating through %d events.' % num_events
  t = time()
  
  variables = [
    'm3a',              # M3 OF TOP QUARK A
    'm2a',              # M2 OF TOP QUARK A
    'angles_b_a',       # AZIMUTHAL ANGLE FOR BOTTOM QUARK A
    'angles_j1_a',      # AZIMUTHAL ANGLE FOR JET 1 IN SYSTEM A
    'angles_j2_a',      # AZIMUTHAL ANGLE FOR JET 2 IN SYSTEM A
    'm3b',              # M3 OF TOP QUARK B
    'm2b',              # M2 OF TOP QUARK B
    'angles_b_b',       # AZIMUTHAL ANGLE FOR BOTTOM QUARK B
    'angles_j1_b',      # AZIMUTHAL ANGLE FOR JET 1 IN SYSTEM B
    'angles_j2_b',      # AZIMUTHAL ANGLE FOR JET 2 IN SYSTEM B
    'missing_eT'        # MISSING TRANVERSE ENERGY E_T
  ]

  output = open('bdt_variables-%s' % argv[1].split('/')[-1], 'w')
  
  print 'Variables used in BDT processing:'
  for (i, var) in zip(range(len(variables)), variables):
    print '%d. %s' % (i+1, var)
  output.write('\t'.join(variables) + '\n')


  for i in range(num_events):
    (events, bottoms, jets) = events_by_file[i]

    if len(bottoms) < 2: continue   # need one bottom quark for each top quark system
    else: enough_bottoms += 1

    # -----------------------------------
    # BEST BJJ COMBO = TOP QUARK A
    # BEST BJJ COMBO FROM REMAINING JETS = TOP QUARK B
    # -----------------------------------
    topA, m31, m21 = get_best_bjj(bottoms, jets);
    if not topA: continue
    else: has_topA += 1

    # Remove bjj of A from the set of bottoms and jets to consider for system B
    bottoms2 = [x for x in bottoms if x != topA[1]]
    jets2 = [x for x in jets if x not in topA[2:]]
    topB, m32, m22 = get_best_bjj(bottoms2, jets2)
    if not topB: continue
    else: has_topB += 1

    # ----------------
    # AZIMUTHAL ANGLES
    # ----------------

    # Sum transverse momentum components of ALL events in collision
    # Theoretically should be 0, but it won't be 
    pT_file_total = sumzip([get_pe(event) for event in events])[:2]

    # "missing" transverse momentum = negative of sum
    pT_missing = [-1*x for x in pT_file_total]
    missing_e = norm(pT_missing)

    # transverse momentum vector of each in bjj of system B
    pT_bjjA = [get_pe(topA[k])[:2] for k in range(1, 4)]
    pT_bjjB = [get_pe(topB[k])[:2] for k in range(1, 4)]

    # azimuthal angle in between each jet and missing transverse momentum
    anglesA = [angle(j, pT_missing) for j in pT_bjjA]
    anglesB = [angle(j, pT_missing) for j in pT_bjjB]

    output.write(
      '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' %
      (m31, m21, anglesA[0], anglesA[1], anglesA[2],
        m32, m22, anglesB[0], anglesB[1], anglesB[2],
        missing_e)
    )

    # if missing_e > 195 and \
    if missing_e > 145 and \
    40 <= m21 <= 120 and \
    120 <= m31 <= 220 and \
    40 <= m22 <= 120 and \
    110 <= m32 <= 230 and \
    anglesB[0] > 1.2 and \
    min(anglesB[1], anglesB[2]) > 0.7:
      meets_cuts += 1

  output.close()
  print 'Done. Took %f secs.\n' % (time()-t)

  print 'Total number of events:', num_events
  print 'Number of events with enough bottom quarks:', enough_bottoms
  print 'Number of events with top quark system A:', has_topA
  print 'Number of events with top quark system B:', has_topB
  print 'Number of events that passed through cuts:', meets_cuts
  eff = meets_cuts * 100.0 / num_events
  print 'Efficiency = %2.2f%%, Rejection = %2.2f%%\n' % (eff, 100-eff)

if __name__ == "__main__":
  main()