from ROOT import *

variables = ['MVA_BDT_S', 'MVA_BDT_B', 'MVA_BDT_S_high', 'MVA_BDT_B_high', 'MVA_BDT_effS', 'MVA_BDT_effB', 'MVA_BDT_effBvsS', 'MVA_BDT_rejBvsS', 'MVA_BDT_invBeffvsSeff', 'MVA_BDT_Train_S', 'MVA_BDT_Train_B', 'MVA_BDT_trainingEffS', 'MVA_BDT_trainingEffB', 'MVA_BDT_trainingEffBvsS', 'MVA_BDT_trainingRejBvsS', 'm3a__Signal', 'm3a__Background', 'm2a__Signal', 'm2a__Background', 'm3b__Signal', 'm3b__Background', 'm2b__Signal', 'm2b__Background', 'angles_b__Signal', 'angles_b__Background', 'angles_j1__Signal', 'angles_j1__Background', 'angles_j2__Signal', 'angles_j2__Background', 'missing_pT__Signal', 'missing_pT__Background', 'BoostWeight', 'BoostWeightVsTree', 'ErrFractHist', 'NodesBeforePruning', 'NodesAfterPruning']

g = TFile.Open('test.root').GetDirectory('Method_BDT').GetDirectory('BDT')
g.cd()

print 'Available variables are:'
print variables
var = raw_input("Choose which variable to plot: ")
while var != 'q':
  c1 = TCanvas("c1","c1",800,800)
  try:
    eval(var).Draw()
    raw_input('Plotting %s. Press any key to close.' % var)
  except NameError:
    print 'No file named %s' % var
  c1.Close()
  var = raw_input("Choose which variable to plot: ")

exit(0)

for v in variables:
  c1 = TCanvas("c1","c1",800,800)
  try:
    eval(v).Draw()
    raw_input('Plotting %s. Press any key to continue.' % v)
  except NameError:
    print 'No file named %s' % v
  c1.Close()