# ------------------------------------------
# Onkur Sen
# 
# process_bdt.py
#
# Usage: python process_bdt.py
# 
# Reads files with information on parameters
# that is processed from get_bdt_variables.py
# and feeds input into a boosted decision
# tree (BDT). Signal-background separation is
# plotted, and classification of events in
# another input file is attempted.
# ------------------------------------------

import ROOT
from time import strftime, localtime
from sys import argv

# ------------------------------------
# PREPARE DATA TO BE FED INTO BDT
# ------------------------------------

susy_source = argv[1]
ttj_source = argv[2]

susy = open('bdt_variables-%s.txt' % susy_source).read().split('\n')
ttj = open('bdt_variables-%s.txt' % ttj_source).read().split('\n')

# Remove variables and blank line
variables = susy.pop(0).split('\t'); susy.pop()
ttj.pop(0); ttj.pop()

print 'Number of SUSY data points from %s: %d' % (susy_source, len(susy))
print 'Number of TTJ data points from %s: %d' % (ttj_source, len(ttj))

# Fill ROOT nTuple with signal and background variables
ntuple = ROOT.TNtuple("ntuple","ntuple","%s:signal" % ':'.join(variables))

for event in susy:
  curr = map(float, event.split('\t')) + [1]
  ntuple.Fill(*curr)

for event in ttj:
  curr = map(float, event.split('\t')) + [0]
  ntuple.Fill(*curr)

raw_input('NTuple populated. Press any key to continue.')

# -------------------------------------------------------------------
# CREATE AND TRAIN BDT USING ROOT TMVA
# Code taken from: 
# http://aholzner.wordpress.com/2011/08/27/a-tmva-example-in-pyroot/
# -------------------------------------------------------------------
print 'NTuple prepared to be fed into BDT'

fout = ROOT.TFile("test.root","RECREATE")

factory = ROOT.TMVA.Factory(
  "TMVAClassification", fout,
  ":".join([
    "!V",
    "!Silent",
    "Color",
    "DrawProgressBar",
    "Transformations=I;D;P;G,D",
    "AnalysisType=Classification"]
  ))

for v in variables: factory.AddVariable(v,"F")

factory.AddSignalTree(ntuple)
factory.AddBackgroundTree(ntuple)

# cuts defining the signal and background sample
sigCut = ROOT.TCut("signal > 0.5")
bgCut = ROOT.TCut("signal <= 0.5")

factory.PrepareTrainingAndTestTree(
  sigCut,   # signal events
  bgCut,    # background events
  ":".join([
    "nTrain_Signal=0",
    "nTrain_Background=0",
    "nTrain_Signal=0",
    "nTrain_Background=0",
    "SplitMode=Alternate",
    "NormMode=NumEvents",
    "!V"
  ]))

method = factory.BookMethod(
  ROOT.TMVA.Types.kBDT,
  "BDT",
  ":".join([
   "!H",
   "!V",
   "NTrees=30", # approximately 1000 events fed to BDT; number of trees should be sqrt(num_events)
   "nEventsMin=150",
   "MaxDepth=3",
   "BoostType=AdaBoost",
   "AdaBoostBeta=0.5",
   "SeparationType=GiniIndex",
   "nCuts=20",
   "PruneMethod=NoPruning"
  ]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

raw_input('BDTs trained using ROOT TMVA in test.root. Press any key to open TBrowser.')
tb = ROOT.TBrowser()
raw_input('TBrowser opened. Press any key to close.')
exit(0)

# ------------------------------------
# PLOT HISTOGRAM OF SIGNAL VS. BACKGROUND 
# TO SHOW FEASIBILITY OF SEPARATION
# ------------------------------------

c1 = ROOT.TCanvas("c1","c1",800,800);

# fill histograms for signal and background from the test sample tree
ROOT.TestTree.Draw("BDT>>hSig(200,-0.5,0.5)","classID == 1","goff")  # signal
ROOT.TestTree.Draw("BDT>>hBg(200,-0.5,0.5)","classID == 0", "goff")  # background

ROOT.hSig.SetLineColor(ROOT.kBlue); ROOT.hSig.SetLineWidth(2)  # signal histogram
ROOT.hBg.SetLineColor(ROOT.kRed); ROOT.hBg.SetLineWidth(2)   # background histogram

# use a THStack to show both histograms
hs = ROOT.THStack("hs","SUSY Signal (Blue) vs. ttj Background (Red)")
hs.Add(ROOT.hSig)
hs.Add(ROOT.hBg)
hs.Draw()

c1.SaveAs("plots/bdt-separation.png")
raw_input("Signal-background separation plotted. Press any key to exit.")