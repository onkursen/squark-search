# > Onkur Sen
# 
# > Usage: python process_bdt.py [susy_filename] [ttj_filename]
# 
# > Reads files with information on parameters
# > that is processed from get_bdt_variables.py
# > and feeds input into a boosted decision
# > tree (BDT). Signal-background separation is
# > plotted, and classification of events in
# > another input file is attempted.

import ROOT
from time import strftime, localtime
from sys import argv, exit
import os

# **Read input data and fill Ntuple with signal and background events**

try:
  susy_source = argv[1]
  ttj_source = argv[2]
except IndexError:
  exit('ERROR: Program requires two arguments; one for susy path and one for background path, e.g.:\npython process_bdt.py 100t400_8TeV-52.txt ttj006f3-6789.txt')

susy = open(susy_source).read().split('\n')
ttj = open(ttj_source).read().split('\n')

# Remove variables at beginning and blank line at end
variables = susy.pop(0).split('\t'); susy.pop()
ttj.pop(0); ttj.pop()

print 'Number of SUSY data points from %s: %d' % (susy_source, len(susy))
print 'Number of TTJ data points from %s: %d' % (ttj_source, len(ttj))

# Make trees the same length
if len(ttj) > len(susy):
  ttj = ttj[:len(susy)]
  print 'SUSY is limiting'
else:
  susy = susy[:len(ttj)]
  print 'TTJ is limiting'

# Number of decision trees should be on the order of sqrt(number of events)
numTrees = int(len(susy)**0.5)
print 'Number of decision trees used in BDT: %d\n' % numTrees

# Fill ROOT nTuple
ntuple = ROOT.TNtuple('ntuple','ntuple','%s:signal' % ':'.join(variables))

for event in susy:
  ntuple.Fill(*(map(float, event.split()) + [1]))

for event in ttj:
  # s = event.split()
  # print s
  # for e in s:
  #   raw_input(float(e))
  ntuple.Fill(*(map(float, event.split()) + [0]))

raw_input('NTuple populated with signal and background events. Press any key to continue.')

# **Train and test BDT using ROOT TMVA**
# http://tmva.sourceforge.net/
# Code modified from: 
# http://aholzner.wordpress.com/2011/08/27/a-tmva-example-in-pyroot/

print 'Running BDT testing and training methods.'

# Define output location and open buffers
if not os.path.isdir('root_files'): os.system('mkdir root_files')
output_location = 'root_files/%s-%s.root' % \
  (susy_source.split('/')[-1][:-4], ttj_source.split('/')[-1][:-4])
fout = ROOT.TFile(output_location, 'RECREATE')

factory = ROOT.TMVA.Factory(
  'TMVAClassification', fout,
  ':'.join([
    '!V',
    '!Silent',
    'Color',
    'DrawProgressBar',
    'Transformations=I;D;P;G,D',
    'AnalysisType=Classification']
  ))

for v in variables:
  factory.AddVariable(v,'F')

factory.AddSignalTree(ntuple)
factory.AddBackgroundTree(ntuple)

# Cuts defining the signal and background sample
sigCut = ROOT.TCut('signal > 0.5')
bgCut = ROOT.TCut('signal <= 0.5')

# options = ['testonSignalOnly', 'testonBackgroundOnly', 'testBoth']
# option = 2

train_test = 'nTrain_Signal=0:nTrain_Background=0'
# fraction_train = len(susy)/3
# train_test = 'nTrain_Signal=%d:nTrain_Background=%d' % (fraction_train, fraction_train)
# if option == 0:
#   train_test = 'nTrain_Signal=0:nTrain_Background=%d' % len(ttj)
# elif option == 1:
#   train_test = 'nTrain_Signal=%d:nTrain_Background=0' % len(susy)

factory.PrepareTrainingAndTestTree(
  sigCut,   # signal events
  bgCut,    # background events
  ':'.join([
    train_test,
    'SplitMode=Alternate',
    'NormMode=NumEvents',
    '!V'
  ]))

method = factory.BookMethod(
  ROOT.TMVA.Types.kBDT,
  'BDT',
  ':'.join([
   'H',
   '!V',
   # 'NTrees=%d' % numTrees,
   'NTrees=%d' % 70,
   'AdaBoostBeta=0.5',
   # 'PruneMethod=NoPruning',
   'nEventsMin=100',
   'PruneStrength=0.5'
  ]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

raw_input('BDT results are in %s. Press any key to open TBrowser.' % output_location)
tb = ROOT.TBrowser()
raw_input('TBrowser opened. Press any key to close.')