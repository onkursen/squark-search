import ROOT

# create a TNtuple
ntuple = ROOT.TNtuple("ntuple","ntuple","x:y:signal")

# generate 'signal' and 'background' distributions
for i in range(10000):
    # throw a signal event centered at (1,1)
    ntuple.Fill(ROOT.gRandom.Gaus(1,1), # x
                ROOT.gRandom.Gaus(1,1), # y
                1)                      # signal
    
    # throw a background event centered at (-1,-1)
    ntuple.Fill(ROOT.gRandom.Gaus(-1,1), # x
                ROOT.gRandom.Gaus(-1,1), # y
                0)                       # background

# keeps objects otherwise removed by garbage collected in a list
gcSaver = []

# create a new TCanvas
gcSaver.append(ROOT.TCanvas())

# draw an empty 2D histogram for the axes
histo = ROOT.TH2F("histo","",1,-5,5,1,-5,5)
histo.Draw()

# draw the signal events in red
ntuple.SetMarkerColor(ROOT.kRed)
ntuple.Draw("y:x","signal > 0.5","same")

# draw the background events in blue
ntuple.SetMarkerColor(ROOT.kBlue)
ntuple.Draw("y:x","signal <= 0.5","same")

raw_input("1Press Enter to continue");

ROOT.TMVA.Tools.Instance()

# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.
fout = ROOT.TFile("test.root","RECREATE")

factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))

factory.AddVariable("x","F")
factory.AddVariable("y","F") 

factory.AddSignalTree(ntuple)
factory.AddBackgroundTree(ntuple)

# cuts defining the signal and background sample
sigCut = ROOT.TCut("signal > 0.5")
bgCut = ROOT.TCut("signal <= 0.5")

factory.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ":".join([
                                        "nTrain_Signal=0",
                                        "nTrain_Background=0",
                                        "SplitMode=Random",
                                        "NormMode=NumEvents",
                                        "!V"
                                       ]))

method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
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


raw_input("2Press Enter to continue");


reader = ROOT.TMVA.Reader()
import array
varx = array.array('f',[0]) ; reader.AddVariable("x",varx)
vary = array.array('f',[0]) ; reader.AddVariable("y",vary)

reader.BookMVA("BDT","weights/TMVAClassification_BDT.weights.xml")

# create a new 2D histogram with fine binning
histo2 = ROOT.TH2F("histo2","",200,-5,5,200,-5,5)

# loop over the bins of a 2D histogram
for i in range(1,histo2.GetNbinsX() + 1):
    for j in range(1,histo2.GetNbinsY() + 1):
        
        # find the bin center coordinates
        varx[0] = histo2.GetXaxis().GetBinCenter(i)
        vary[0] = histo2.GetYaxis().GetBinCenter(j)
        
        # calculate the value of the classifier
        # function at the given coordinate
        bdtOutput = reader.EvaluateMVA("BDT")
        
        # set the bin content equal to the classifier output
        histo2.SetBinContent(i,j,bdtOutput)

gcSaver.append(ROOT.TCanvas())
histo2.Draw("colz")

# draw sigma contours around means
for mean, color in (
    ((1,1), ROOT.kRed), # signal
    ((-1,-1), ROOT.kBlue), # background
    ):
    
    # draw contours at 1 and 2 sigmas
    for numSigmas in (1,2):
        circle = ROOT.TEllipse(mean[0], mean[1], numSigmas)
        circle.SetFillStyle(0)
        circle.SetLineColor(color)
        circle.SetLineWidth(2)
        circle.Draw()
        gcSaver.append(circle)

ROOT.gPad.Modified()

raw_input("3Press Enter to continue");

# fill histograms for signal and background from the test sample tree
ROOT.TestTree.Draw("BDT>>hSig(22,-1.1,1.1)","classID == 0","goff")  # signal
ROOT.TestTree.Draw("BDT>>hBg(22,-1.1,1.1)","classID == 1", "goff")  # background

ROOT.hSig.SetLineColor(ROOT.kRed); ROOT.hSig.SetLineWidth(2)  # signal histogram
ROOT.hBg.SetLineColor(ROOT.kBlue); ROOT.hBg.SetLineWidth(2)   # background histogram

# use a THStack to show both histograms
hs = ROOT.THStack("hs","")
hs.Add(ROOT.hSig)
hs.Add(ROOT.hBg)

# show the histograms
gcSaver.append(ROOT.TCanvas())
hs.Draw()

raw_input("4Press Enter to continue");