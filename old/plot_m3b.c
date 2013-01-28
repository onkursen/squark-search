{
  #include <iostream.h>
  #include <iomanip.h>
  ifstream inp; double x;

  TH1F* h = new TH1F("M3B", "Distribution of M3 for Top Quark B", 100, 0, 4000);
  inp.open("m3b.txt");
  while(!(inp >> x) == 0) {h.Fill(x);}
  inp.close();
  
  gStyle->SetOptStat();
  h.Draw();

  c1->SaveAs("plot_m3b.png");
  return 0;
}