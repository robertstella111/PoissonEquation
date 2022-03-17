#include "ONED_Poisson.h"
#include "../eigen-3.4.0/Eigen/Core"
#include "../eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

ONED::ONED(long N, double a, double b, double (*fun)(double)) {
  (*this).N = N;
  (*this).a = a;
  (*this).b = b;
  (*this).fun = fun;
  n1 = 0;
  n2 = 0;
  d1 = 0;
  d2 = 0;
  rand1 = 0;
  rand2 = 0;
  h = (b - a) / N;
}
ONED::~ONED() { cout << "Beenden" << endl; }

void ONED::setDirichlet_Rand1(double value) {
  d1 = value;
  if (rand1 == 2)
    rand1 = 3;
  else
    rand1 = 1;
}
void ONED::setDirichlet_Rand2(double value) {
  d2 = value;
  if (rand2 == 2)
    rand2 = 3;
  else
    rand2 = 1;
}
void ONED::setDirichlet(double a, double b) {
  d1 = a;
  d2 = b;
  if (rand1 == 2)
    rand1 = 3;
  else
    rand1 = 1;
  if (rand2 == 2)
    rand2 = 3;
  else
    rand2 = 1;
}
void ONED::setNeumann_Rand1(double value) {
  n1 = value;
  if (rand1 == 1)
    rand1 = 3;
  else
    rand1 = 2;
}
void ONED::setNeumann_Rand2(double value) {
  n2 = value;
  if (rand2 == 1)
    rand2 = 3;
  else
    rand2 = 2;
}
void ONED::setNeumann(double a, double b) {
  n1 = a;
  n2 = b;
  if (rand1 == 1)
    rand1 = 3;
  else
    rand1 = 2;
  if (rand2 == 1)
    rand2 = 3;
  else
    rand2 = 2;
}

Eigen::VectorXd ONED::solve() {
  // Derzeit noch keine gemischten Randwertprobleme
  Eigen::VectorXd richtig(N + 1);
  if ((rand1 == 0 || rand1 == 1) && (rand2 == 0 || rand2 == 1)) {
    // N+1 x N+1 Matrix,
    SparseMatrix<double> mat(N + 1, N + 1);
    Eigen::VectorXd vec(N + 1);
    Eigen::VectorXd ergebnis(N + 1);
    for (long i = 0; i < N + 1; i++) {
      if (i == 0)
        vec[i] = d1;
      else if (i == N)
        vec[i] = d2;
      else
        vec[i] = pow(h, 2) * fun(a + i * h);
    }
    for (int k = 0; k < N + 1; k++) {
      if (k == 0)
        mat.insert(k, k) = 1;
      else if (k == N) {
        mat.insert(k, k) = 1;
      } else {
        mat.insert(k, k) = 2;
        mat.insert(k, k - 1) = -1;
        mat.insert(k, k + 1) = -1;
      }
    }
    // mat.setFromTriplets(tr.begin(), tr.end());
    BiCGSTAB<SparseMatrix<double>> lscg;
    lscg.setMaxIterations(N);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    return ergebnis;
  } else if (rand1 == 2 && rand2 == 2) {
    // N+1 x N+1 Matrix (Anfangswert u0 wird auf 0 gesetzt)
    SparseMatrix<double> mat(N + 1, N + 1);
    Eigen::VectorXd vec(N + 1);
    Eigen::VectorXd ergebnis(N + 1);
    for (long i = 0; i < N + 1; i++) {
      if (i == 0)
        vec[i] = -fun(a + i * h) + 2 * n1 / h;
      else if (i == N)
        vec[i] = -fun(a + i * h) - 2 * n2 / h;
      else
        vec[i] = pow(h, 2) * fun(a + i * h);
    }
    for (int k = 0; k < N + 1; k++) {
      if (k == 0) {
        mat.insert(k, k + 1) = 2 / pow(h, 2);
      } else if (k == N) {
        mat.insert(k, k - 1) = 2 / pow(h, 2);
        mat.insert(k, k) = -2 / pow(h, 2);
      } else {
        mat.insert(k, k) = 2;
        mat.insert(k, k - 1) = -1;
        mat.insert(k, k + 1) = -1;
      }
    }
    BiCGSTAB<SparseMatrix<double>> lscg;
    lscg.setMaxIterations(N * 2);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    return ergebnis;
  } else {
    // NxN Matrix, Bei einem Rand handelt es sich um Dirichlet beim anderen um
    // Neumann Bedingungen
    SparseMatrix<double> mat(N + 1, N + 1);
    Eigen::VectorXd vec(N + 1);
    Eigen::VectorXd ergebnis(N + 1);
    vector<Triplet<double>> tr;
    if (rand1 == 2) {
      // Neumann bei a u0
      for (long i = 0; i < N + 1; i++) {
        if (i == 0)
          vec[i] = -fun(a + i * h) + 2 * n1 / h;
        else if (i == N)
          vec[i] = d2;
        else
          vec[i] = pow(h, 2) * fun(a + i * h);
      }
      for (int k = 0; k < N + 1; k++) {
        if (k == 0) {
          mat.insert(k, k + 1) = 2 / pow(h, 2);
          mat.insert(k, k) = -2 / pow(h, 2);
        } else if (k == N) {
          mat.insert(k, k) = 1;
        } else {
          mat.insert(k, k) = 2;
          mat.insert(k, k - 1) = -1;
          mat.insert(k, k + 1) = -1;
        }
      }
      BiCGSTAB<SparseMatrix<double>> lscg;
      lscg.setMaxIterations(N * 2);
      lscg.compute(mat);
      ergebnis = lscg.solve(vec);
      return ergebnis;
    } else {
      cout << d1 << " " << n2;
      for (long i = 0; i < N + 1; i++) {
        if (i == 0)
          vec[i] = d1;
        else if (i == N)
          vec[i] = vec[i] = -fun(a + i * h) - 2 * n2 / h;
        else
          vec[i] = pow(h, 2) * fun(a + i * h);
      }
      for (int k = 0; k < N + 1; k++) {
        if (k == 0) {
          mat.insert(k, k) = 1;
        } else if (k == N) {
          mat.insert(k, k - 1) = 2 / pow(h, 2);
          mat.insert(k, k) = -2 / pow(h, 2);
        } else {
          mat.insert(k, k) = 2;
          mat.insert(k, k - 1) = -1;
          mat.insert(k, k + 1) = -1;
        }
      }
      BiCGSTAB<SparseMatrix<double>> lscg;
      lscg.setMaxIterations(N * 2);
      // lscg.setTolerance(1e-);
      lscg.compute(mat);
      ergebnis = lscg.solve(vec);
      return ergebnis;
    }
  }
}

void ONED::write(Eigen::VectorXd vec, bool print) {
  fstream myFile;
  myFile.open("Daten.txt", ios::out);
  if (!myFile) {
    cout << "Nicht geoeffnet!\n";
  } else {
    for (long i = 0; i < vec.size(); i++) {
      myFile << a + i * h << "|" << vec[i] << endl;
    }
    myFile.close();
    /*
    if(print == true) {
        system("ONED.py");
        system("start OneDPoissongleichung.pdf");
    }
    */
  }
}
