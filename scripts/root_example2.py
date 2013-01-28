import ROOT
import numpy as n

x = n.linspace(0, 4*n.pi,101)
y = n.cos(x)

g = ROOT.TGraph(len(x), x,y)
g.SetTitle("cosine in x=[%.1f, %.1f]" % (x[0], x[-1]))
g.GetXaxis().SetTitle("x")
g.GetYaxis().SetTitle("y")
g.Draw("AL")
raw_input("Press Enter to quit");
