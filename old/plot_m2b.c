{
  #include <iostream.h>
  ifstream inp; double x;

  TH1F* h2 = new TH1F("M2B", "Distribution of M2 for Top Quark B", 45, 0, 900);
  inp.open("m2b.txt");
  while(!(inp >> x) == 0) {h2.Fill(x);}
  inp.close();
  
  gStyle->SetOptStat();
  h2.Draw();

  c1->SaveAs("plot_m2b.png");
  return 0;
}