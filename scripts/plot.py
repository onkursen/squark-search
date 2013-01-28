from ROOT import TH1F, TCanvas, TLegend

def plot_var(var, directory):
  histo = TH1F(
    var.upper(),
    "Distribution of %s for Top Quark %s" % 
    (var[:-1].upper(), var[-1].upper()),
    45, 0, 900)
  for line in open("output/output_%s/%s.txt" % (directory, var)):
    histo.Fill(float(line))

  c1 = TCanvas("c1","c1",800,800)
  histo.Draw()
  c1.SaveAs("plots/plot_%s_%s.png" % (directory, var))

def plot_angle(var):
  c1 = TCanvas("c1","c1",800,800)

  susy = TH1F("Phi_SUSY", "Azimuthal Angles", 35, 0, 3.5)
  for line in open("output/output_susy/angles_%s.txt" % var):
    susy.Fill(float(line))
  susy.SetLineColor(4) # blue
  if susy.Integral()!=0: susy.Scale(1/susy.Integral())
  susy.SetMaximum(0.06)
  susy.GetYaxis().SetTitleSize(0.05)
  susy.SetLabelSize(0.05,"x")
  susy.GetXaxis().SetTitleSize(0.05)
  susy.SetXTitle("#Delta#phi")
  susy.SetYTitle("Fraction / 0.1")
  susy.Draw()

  ttj = TH1F("Phi_ttj", "Azimuthal Angles of Bottom Quark", 35, 0, 3.5)
  for line in open("output/output_ttj/angles_%s.txt" % var):
    ttj.Fill(float(line))
  ttj.SetLineColor(8) # green
  if (ttj.Integral()!=0): ttj.Scale(1/ttj.Integral())
  ttj.Draw("same")

  leg_hist = TLegend(0.7,0.8,0.9,0.95)
  leg_hist.AddEntry(susy,"SUSY","l")
  leg_hist.AddEntry(ttj,"ttj","l")
  leg_hist.Draw()

  c1.SaveAs("plots/plot_angles_%s.png" % var)

plot_var('m3a', 'susy')
plot_var('m3a', 'ttj')
print 'plotted m3a'
plot_var('m2a', 'susy')
plot_var('m2a', 'ttj')
print 'plotted m2a'
plot_var('m3b', 'susy')
plot_var('m3b', 'ttj')
print 'plotted m3b'
plot_var('m2b', 'susy')
plot_var('m2b', 'ttj')
print 'plotted m2b'
plot_angle('b')
print 'plotted b angle'
plot_angle('j1')
print 'plotted j1 angle'
plot_angle('j2')
print 'plotted j2 angle'