#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Runge_Kutta.h"

using namespace std;

Eigen::VectorXd runge(double (*fun)(double, double), double y0, double t0,
                      double end, long N) {
  double h = (end - t0) / N;
  Eigen::VectorXd ergebnis(N + 1);
  double k1, k2, k3, k4;
  ergebnis[0] = y0;
  for (long x = 0; x < N; x++) {
    k1 = fun(t0, y0);
    k2 = fun(t0 + h / 2, y0 + (h / 2) * k1);
    k3 = fun(t0 + h / 2, y0 + (h / 2) * k2);
    k4 = fun(t0 + h, y0 + h * k3);
    y0 = y0 + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    t0 = t0 + h;
    ergebnis[x + 1] = y0;
  }
  return ergebnis;
}

void funktion() {
  fstream myFile;
  myFile.open("nin_in.deva", ios::in);
  if (!myFile) {
    cout << "Nicht geoeffnet!\n";
  } else {
    cout << "geoeffnet" << endl;
    string s, b;
    int a, zaehler;
    stringstream ss;
    while (getline(myFile, s)) {

      while (s.find("Pointlist(") == string::npos)
        getline(myFile, s);
      stringstream str_strm;
      str_strm << s; // convert the string s into stringstream
      string temp_str;
      int temp_int;
      zaehler = 0;
      getline(str_strm, temp_str, '"');
      getline(str_strm, temp_str, '"');
      cout << temp_str << endl;
      while (!str_strm.eof()) {
        getline(str_strm, temp_str, ',');
        if (stringstream(temp_str) >> temp_int) { // try to convert string to
                                                  // int
          if (zaehler == 1)
            a = temp_int;
          zaehler++;
        }
        temp_str = ""; // clear temp string
      }
      break;
      //
    }
    cout << a;
  }
}
