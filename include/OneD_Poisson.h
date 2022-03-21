#ifndef ONED_POISSON_H_INCLUDED
#define ONED_POISSON_H_INCLUDED

#include "../eigen-3.4.0/Eigen/Dense"
#include "../eigen-3.4.0/Eigen/Sparse"
#define epsilon 8.854e-12

class ONED {
private:
  // defaultm�ssig werden alle Randbedingungen auf 0 gesetzt
  double (*fun)(double);
  double d1, d2; // Dirichlet Randwerte;
  double n1, n2; // Neumann Bedingungen;

  long N;
  double h, a, b;    // Gitterbreite, Intervall der Funktion
  long rand1, rand2; // 0...homogene Dirichet, 1...
                     // Dirichlet, 2...Neumann, 3....gemischte Randwertprobleme

public:
  void setDirichlet_Rand1(double value);
  void setDirichlet_Rand2(double value);
  void setDirichlet(double a, double b);
  void setNeumann_Rand1(double value);
  void setNeumann_Rand2(double value);
  void setNeumann(double a, double b);
  Eigen::VectorXd
  solve(); // L�sen des Gleichungssystems und als R�ckgabe den Ergebnisvektor
  void write(Eigen::VectorXd vec, bool print);

  ONED(long N, double a, double b, double (*fun)(double));
  ~ONED();
};

#endif // 1D_POISSON_H_INCLUDED
