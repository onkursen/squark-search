{
  #include <iostream.h>
  #include <iomanip.h>
  ifstream inp; double x;

  TH1F* h = new TH1F("M3A", "Distribution of M3 for Top Quark A", 45, 0, 900);
  inp.open("m3a.txt");
  while(!(inp >> x) == 0) {h.Fill(x);}
  h.Draw();
  inp.close();

  return 0;
}
