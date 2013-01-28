{
  #include <iostream.h>
  ifstream inp; double x;

  TH1F* h3 = new TH1F("Phi_SUSY", "Azimuthal Angles", 35, 0, 3.5);
  inp.open("susy_angles_j2.txt");
  while(!(inp >> x) == 0) {h3.Fill(x);}
  inp.close();

  h3.SetLineColor(4); // blue
  gStyle->SetOptStat(0);
  if (h3->Integral()!=0) {h3->Scale(1/h3->Integral());}
  h3->SetMaximum(0.06);
  
  h3->GetYaxis()->SetTitleSize(0.05);
  h3->SetLabelSize(0.05,"x");
  h3->GetXaxis()->SetTitleSize(0.05);
  h3->SetXTitle("#Delta#phi");
  h3->SetYTitle("Fraction / 0.1");

  h3.Draw();

  TH1F* h2 = new TH1F("Phi_ttj", "Azimuthal Angles of Bottom Quark", 35, 0, 3.5);
  inp.open("ttj_angles_j2.txt");
  while(!(inp >> x) == 0) {h2.Fill(x);}
  inp.close();
  h2.SetLineColor(8); // green
  if (h2->Integral()!=0) {h2->Scale(1/h2->Integral());}
  h2.Draw("same");

  leg_hist = new TLegend(0.7,0.8,0.9,0.95);
  // leg_hist->SetHeader("Some histograms");
  leg_hist->AddEntry(h3,"SUSY","l");
  leg_hist->AddEntry(h2,"ttj","l");
  leg_hist->Draw();

  c1->SaveAs("plot_angles_j2.png");

  return 0;
}