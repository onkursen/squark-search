{
  #include <iostream.h>
  #include <iomanip.h>
  ifstream inp; double x;

  TH1F* h = new TH1F("M3B", "Distribution of M3 for Top Quark B", 100, 0, 4000);
  inp.open("m3b.txt");
  while(!(inp >> x) == 0) {h.Fill(x);}
  h.Draw();
  inp.close();

  return 0;
}