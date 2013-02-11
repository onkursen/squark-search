import ROOT
from time import strftime, localtime

# ------------------------------------
# PREPARE DATA TO BE FED INTO BDT
# ------------------------------------

def to_float_array(filename):
  return map(float, 
    open(filename).read().split('\n')[:-1])

# read files contained in `variables` array and append them
# to signal and background lists
variables = ['m3a', 'm2a', 'm3b', 'm2b', 'angles_b', 'angles_j1', 'angles_j2']
susy = []; ttj = [];
for v in variables:
  susy.append(to_float_array('output/output_susy/%s.txt' % v))
  ttj.append(to_float_array('output/output_ttj/%s.txt' % v))

# fill ROOT nTuple with signal and background variables
ntuple = ROOT.TNtuple("ntuple","ntuple","%s:signal" % ':'.join(variables))
for i in range(min([len(var) for var in susy])):
  curr = [var[i] for var in susy] + [1] # susy is signal
  ntuple.Fill(*curr)

for i in range(min([len(var) for var in ttj])):
  curr = [var[i] for var in ttj] + [0] # ttbar is background
  ntuple.Fill(*curr)

print 'NTuple prepared to be fed into BDT'

# -------------------------------------------------------------------
# CREATE AND TRAIN BDT USING ROOT TMVA
# Code taken from: 
# http://aholzner.wordpress.com/2011/08/27/a-tmva-example-in-pyroot/
# -------------------------------------------------------------------

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

for v in variables:
  factory.AddVariable(v,"F")

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
      "SplitMode=Random",
      "NormMode=NumEvents",
      "!V"
     ]))

method = factory.BookMethod(
  ROOT.TMVA.Types.kBDT,
  "BDT",
  ":".join([
     "!H",
     "!V",
     "NTrees=850",
     "nEventsMin=150",
     "MaxDepth=3",
     "BoostType=AdaBoost",
     "AdaBoostBeta=0.5",
     "SeparationType=GiniIndex",
     "nCuts=20",
     "PruneMethod=NoPruning",
     ]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

print 'BDTs trained using ROOT TMVA'

# ------------------------------------
# PLOT HISTOGRAM OF SIGNAL VS. BACKGROUND 
# TO SHOW FEASIBILITY OF SEPARATION
# ------------------------------------

c1 = ROOT.TCanvas("c1","c1",800,800);

# fill histograms for signal and background from the test sample tree
ROOT.TestTree.Draw("BDT>>hSig(22,-1.1,1.1)","classID == 0","goff")  # signal
ROOT.TestTree.Draw("BDT>>hBg(22,-1.1,1.1)","classID == 1", "goff")  # background

ROOT.hSig.SetLineColor(ROOT.kRed); ROOT.hSig.SetLineWidth(2)  # signal histogram
ROOT.hBg.SetLineColor(ROOT.kBlue); ROOT.hBg.SetLineWidth(2)   # background histogram

# use a THStack to show both histograms
hs = ROOT.THStack("hs","SUSY Signal (Red) vs. ttj Background (Blue)")
hs.Add(ROOT.hSig)
hs.Add(ROOT.hBg)
hs.Draw()

# legend = ROOT.TLegend(0.4,0.6,0.89,0.89); 
# legend.AddEntry(ROOT.hSig,"Signal")
# legend.AddEntry(ROOT.hBg,"Background")
# legend.Draw()

c1.SaveAs("plots/bdt %s.png" % strftime("%Y-%m-%d %H:%M:%S", localtime()))

raw_input("Press any key to close")
