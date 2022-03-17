#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "Diffusion1D.h"
#include "Diffusion_file.h"
#include "OneD_Poisson.h"
#include "Runge_Kutta.h"
#include "ThreeD_Poisson.h"
#include "TwoD_Poisson.h"

#define N 80
#define T 100
#define N2 80
#define N3 80

using namespace Eigen;
using namespace std;
typedef Eigen::SparseMatrix<double> SpMat;

double maximumNorm(Eigen::VectorXd vec1, Eigen::VectorXd vec2, long n) {
  double ergebnis = 0;
  for (long i = 0; i < vec1.size(); i++) {
    if (vec1[i] > vec2[i]) {
      if (vec1[i] - vec2[i] > ergebnis)
        ergebnis = vec1[i] - vec2[i];
    } else {
      if (vec2[i] - vec1[i] > ergebnis)
        ergebnis = vec2[i] - vec1[i];
    }
  }
  return ergebnis;
}

double Ladungsverteilung2D(double x, double y) { return 2 * sin(x) * cos(y); }
double Dirichlet(double x, double y) { return sin(x) * cos(y); }
double Neumannx(double x, double y) { return cos(x) * cos(y); }
double Neumanny(double x, double y) { return (-1) * sin(x) * sin(y); }

double Dirichlet3D(double x, double y, double z) {
  return sin(x) * sin(y) * sin(z);
}
double Ladunsverteilung3D(double x, double y, double z) {
  return 3 * sin(x) * sin(y) * sin(z);
}

double Ladungverteilung(double x) {
  if (0 <= x && x <= 0.2)
    return 0;
  else if (x <= 0.5)
    return 10e-11 / epsilon;
  else if (x <= 0.8)
    return -10e-11 / epsilon;
  else
    return 0;
}

double right(double t, double y) { return y + exp(t) * cos(t); }

double richtig(double t) { return exp(t) * sin(t); }

int main() {
  /*
      Eigen::VectorXd soll(N+1);
      double ende = 10, t0 = 0;
      Eigen::VectorXd ergebnis = runge(right,0, t0, ende, N);
      double h = (ende - t0)/N;
      for(long k = 0; k < N + 1; k++) {
          soll[k] = richtig(t0 + k*h);
      }
      cout << "Fehler :" << maximumNorm(soll, ergebnis, N);
  */

  double a = 0, b = M_PI, c = 0, d = M_PI, e = 0, f = M_PI, h1, h2, h3;

  h1 = (b - a) / N;
  h2 = (d - c) / N2;
  h3 = (f - e) / N3;

  // Diverse Testf�llte
  // 1-dimensional:

  /*
      ONED test(N,a,b,sin);
      Eigen::VectorXd richtig(N+1);
      for(long k = 0;k < N+1; k++) {
          richtig[k] = sin(a + k*h1);
      }
      test.setDirichlet_Rand1(0);
      test.setNeumann_Rand2(-1);
      Eigen::VectorXd ergebnis = test.solve();
      cout << "Maximumnorm : " << maximumNorm(richtig, ergebnis,N);
      test.write(ergebnis, true);
  */

  // 2-dimensional
  /*

      TWOD test(N, N2, a, b, c, d, Ladungsverteilung2D, Neumannx, Neumanny);
      //TWOD test(N,N2,a, b,c,d,Ladungsverteilung2D, Dirichlet);
      Eigen::VectorXd richtig((N+1)*(N2+1));
      //test.solvewitheFV();
      for(long i = 0; i < N+1; i++) {
         for(long j = 0; j < N2+1; j++) {
              richtig[i*(N2+1)+j] = Dirichlet(a+i*h1, c+j*h2);
         }
      }
      Eigen::VectorXd ergebnis = test.solve();
      cout << "Maximumnorm : " << maximumNorm(richtig, ergebnis,N);
     test.write(ergebnis, true);
  */

  // 3-dimensional

  /*
      ThreeD test(N,N2,N3, a,b,c,d,e,f,Ladunsverteilung3D, Dirichlet3D);
      Eigen::VectorXd richtig((N+1)*(N2+1)*(N3+1));
      for(long k = 0; k < N3 +1; k++) {
          for(long i = 0; i < N2+1; i++) {     //
              for(long j = 0; j < N+1; j++) {
                  richtig[k*(N+1)*(N2+1)+i*(N+1)+j] = Dirichlet3D(a+j*h1, c + i
     *h2, e + k*h3);
              }
          }
      }
      //test.solvewitheFV();
      Eigen::VectorXd ergebnis = test.solve();
      cout << "Maximumnorm : " << maximumNorm(richtig, ergebnis,N);

  */
  // Diffusion

  // Testfall 1
  /*

      Diffusion1D test(10,N,T,a,b,1, cos);
      Eigen::VectorXd ergebnis = test.solve();
      test.write(ergebnis, true);

  */
  // Testfall 2
  /*
      Diffusion1D test(150,N,T,a, b,0.01, cos);
      test.setoutFlux();
      Eigen::VectorXd ergebnis = test.solve();
      test.write(ergebnis, true);
  */

  Diff_gen test("../Input//nin_in.deva");
  //test.solve3();  //Errechnen der Ladunsverteilung und danach Potential berechnen
  test.test1();
  //cout << "Test" << endl;
  return 0;
}
