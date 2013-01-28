import ROOT, os, sys

filename = sys.argv[1]

histo = ROOT.TH1F(
  filename.upper(),
  "Distribution of %s for Top Quark %s" % 
  (filename[:-1].upper(), filename[-1].upper()),
  45, 0, 900);
for line in open("%s.txt" % filename):
  histo.Fill(float(line))

c1 = ROOT.TCanvas("c1","c1",800,800);
histo.Draw()
c1.SaveAs("plot_%s.png" % filename)
# os.system('open plot_%s.png' % filename)