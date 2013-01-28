{
  #include <iostream.h>
  #include <iomanip.h>
  ifstream inp; double x;

  TH1F* h = new TH1F("M3A", "Distribution of M3 for Top Quark A", 45, 0, 900);
  inp.open("m3a.txt");
  while(!(inp >> x) == 0) {h.Fill(x);}
  inp.close();
  
  gStyle->SetOptStat();
  h.Draw();

  c1->SaveAs("plot_m3a.png");
  return 0;
}
