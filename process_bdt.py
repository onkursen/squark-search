# Onkur Sen
# 
# Usage: python process_bdt.py [susy_filename] [ttj_filename]
# 
# Reads files with information on parameters
# that is processed from get_bdt_variables.py
# and feeds input into a boosted decision
# tree (BDT). Signal-background separation is
# plotted, and classification of events in
# another input file is attempted.

import ROOT
from time import strftime, localtime
from sys import argv, exit

# ### Fill Ntuple with signal and background data to be fed into BDT

try:
  susy_source = argv[1]
  ttj_source = argv[2]
except IndexError:
  exit('ERROR: Program requires two arguments; one for susy path and one for background path, e.g., "python process_bdt.py 100t400_8TeV-52.txt ttj006f3-6789.txt"')

susy = open(susy_source).read().split('\n')
ttj = open(ttj_source).read().split('\n')

# Remove variables at beginning and blank line at end
variables = susy.pop(0).split('\t'); susy.pop()
ttj.pop(0); ttj.pop()

print 'Number of SUSY data points from %s: %d' % (susy_source, len(susy))
print 'Number of TTJ data points from %s: %d' % (ttj_source, len(ttj))

# Number of decision trees should be on the order of sqrt(number of events)
numTrees = int(min(len(susy), len(ttj))**0.5)
print 'Number of decision trees used in BDT:', numTrees

# Fill ROOT nTuple with signal and background variables
ntuple = ROOT.TNtuple("ntuple","ntuple","%s:signal" % ':'.join(variables))

# Note: '*' needed for pointer
for event in susy:
  ntuple.Fill(*(map(float, event.split('\t')) + [1]))

for event in ttj:
  ntuple.Fill(*(map(float, event.split('\t')) + [0]))

raw_input('NTuple populated with signal and background events. Press any key to continue.')

# ### Train and test BDT using ROOT TMVA
# http://tmva.sourceforge.net/
# Code modified from: 
# http://aholzner.wordpress.com/2011/08/27/a-tmva-example-in-pyroot/

print 'Running BDT testing and training methods.'

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

# Cuts defining the signal and background sample
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
   "NTrees=%d" % numTrees,
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

raw_input('BDT results are in in test.root. Press any key to open TBrowser.')
tb = ROOT.TBrowser()
raw_input('TBrowser opened. Press any key to close.')