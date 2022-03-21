#include "../eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include <iostream>

#include "ThreeD_Poisson.h"

using namespace std;
using namespace Eigen;

ThreeD::ThreeD(long N1, long N2, long N3, double a, double b, double c,
               double d, double e, double f,
               double (*fun)(double, double, double),
               double (*rand)(double, double, double)) {
  // F�r Dirichlet
  dirichlet = rand;
  (*this).N1 = N1;
  (*this).N2 = N2;
  (*this).N3 = N3;
  (*this).a = a;
  (*this).b = b;
  (*this).c = c;
  (*this).d = d;
  (*this).e = e;
  (*this).f = f;
  h2 = (d - c) / N2;
  h1 = (b - a) / N1;
  h3 = (f - e) / N3;
  (*this).fun = fun;
  solvertype = 1;
}

Eigen::VectorXd ThreeD::solve() {
  if (solvertype == 1)
    return solve2();
  else
    return solve1();
}

void ThreeD::solvewitheFV() {
  solvertype = 2;
  return;
}
void ThreeD::solvewithFD() {
  solvertype = 1;
  return;
}

Eigen::VectorXd ThreeD::solve1() {
  SparseMatrix<double> mat((N1 + 1) * (N2 + 1) * (N3 + 1),
                           (N1 + 1) * (N2 + 1) * (N3 + 1));
  Eigen::VectorXd vec((N1 + 1) * (N2 + 1) * (N3 + 1));
  Eigen::VectorXd ergebnis((N1 + 1) * (N2 + 1) * (N3 + 1));
  vector<Triplet<double>> tr;
  // Dirichlet Randwerte
  // F�r alle folgenden Schleifen, j ist f�r die x-Achse, i ist f�r die y-Achse,
  // k ist f�r die z-Achse
  for (long k = 0; k < N3 + 1; k++) {
    for (long i = 0; i < N1 + 1; i++) {
      for (long j = 0; j < N2 + 1; j++) {
        if (i == 0 || i == N1 || j == 0 || j == N2 || k == 0 || k == N3) {
          vec[k * (N1 + 1) * (N2 + 1) + i * (N1 + 1) + j] =
              dirichlet(a + j * h1, c + i * h2, e + k * h3);
        } else
          vec[k * (N1 + 1) * (N2 + 1) + i * (N1 + 1) + j] =
              (-1) * h1 * h2 * h3 * fun(a + j * h1, c + i * h2, e + k * h3);
      }
    }
  }
  for (long k = 0; k < N3 + 1; k++) {
    for (long i = 0; i < N1 + 1; i++) { //
      for (long j = 0; j < N2 + 1; j++) {
        if (i == 0 || i == N1 || j == 0 || j == N2 || k == 0 || k == N3) {
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) = 1;
        } else {
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              (-2) * (h2 * h3) / h1 - (2 * h1 * h3) / h2 - (2 * h1 * h2) / h3;
          // mat.insert(k*(N2+1)*(N1+1)+i*(N1+1)+j,k*(N2+1)*(N1+1)+i*(N1+1)+j) =
          // ((-2)*h1*h3)/h2;
          // mat.insert(k*(N2+1)*(N1+1)+i*(N1+1)+j,k*(N2+1)*(N1+1)+i*(N1+1)+j) =
          // ((-2)*h1*h2)/h3;

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     (k + 1) * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              (h1 * h2) / h3;
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     (k - 1) * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              (h1 * h2) / h3;

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + (i + 1) * (N1 + 1) + j) =
              (h1 * h3) / h2;
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + (i - 1) * (N1 + 1) + j) =
              (h1 * h3) / h2;

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j + 1) =
              (h2 * h3) / h1;
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j - 1) =
              (h2 * h3) / h1;
        }
      }
    }
  }
  BiCGSTAB<SparseMatrix<double>> lscg;
  lscg.setMaxIterations((N1 + 1) * (N2 + 1));
  lscg.compute(mat);
  ergebnis = lscg.solve(vec);
  return ergebnis;
}

Eigen::VectorXd ThreeD::solve2() {
  SparseMatrix<double> mat((N1 + 1) * (N2 + 1) * (N3 + 1),
                           (N1 + 1) * (N2 + 1) * (N3 + 1));
  Eigen::VectorXd vec((N1 + 1) * (N2 + 1) * (N3 + 1));
  Eigen::VectorXd ergebnis((N1 + 1) * (N2 + 1) * (N3 + 1));
  // Dirichlet Randwerte
  // F�r alle folgenden Schleifen, j ist f�r die x-Achse, i ist f�r die y-Achse,
  // k - z-Achse
  for (long k = 0; k < N3 + 1; k++) {
    for (long i = 0; i < N1 + 1; i++) {
      for (long j = 0; j < N2 + 1; j++) {
        if (i == 0 || i == N1 || j == 0 || j == N2 || k == 0 || k == N3) {
          vec[k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j] =
              dirichlet(a + j * h1, c + i * h2, e + k * h3);
        } else
          vec[k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j] =
              (-1) * fun(a + j * h1, c + i * h2, e + k * h3);
      }
    }
  }
  for (long k = 0; k < N3 + 1; k++) {
    for (long i = 0; i < N1 + 1; i++) { //
      for (long j = 0; j < N2 + 1; j++) {
        if (i == 0 || i == N1 || j == 0 || j == N2 || k == 0 || k == N3) {
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) = 1;
        } else {
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              (-2) * (1 / pow(h1, 2) + 1 / pow(h2, 2) + 1 / pow(h3, 2));

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + (j + 1)) =
              1 / pow(h1, 2);
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + (j - 1)) =
              1 / pow(h1, 2);

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + (i + 1) * (N1 + 1) + j) =
              1 / pow(h2, 2);
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     k * (N2 + 1) * (N1 + 1) + (i - 1) * (N1 + 1) + j) =
              1 / pow(h2, 2);

          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     (k + 1) * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              1 / pow(h3, 2);
          mat.insert(k * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j,
                     (k - 1) * (N2 + 1) * (N1 + 1) + i * (N1 + 1) + j) =
              1 / pow(h3, 2);
          /*

          mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (1/pow(h1,2));
          mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (1/pow(h1,2));
          mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (1/pow(h2,2));
          mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (1/pow(h2,2));
          mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 *
          (1/pow(h1,2));

          */
        }
      }
    }
  }
  BiCGSTAB<SparseMatrix<double>> lscg;
  lscg.setMaxIterations((N1 + 1) * (N2 + 1));
  lscg.compute(mat);
  ergebnis = lscg.solve(vec);
  cout << "FD";
  return ergebnis;
}

/*
Eigen::VectorXd ThreeD::solvetest() {
        SparseMatrix<double> mat((N1+1)*(N2+1), (N1+1)*(N2+1));
        Eigen::VectorXd vec((N1+1)*(N2+1));
        Eigen::VectorXd ergebnis((N1+1)*(N2+1));
        for(long i = 0; i < N1+1; i++) {     //Vector Initialisieren
            for(long j = 0; j < N2+1; j++) {
                for(long z = 0; z < N3+1; z++) {
                    if(i == 0 || i == N1 || j == 0 || j == N2) {
                        if(i == 0 && j != 0 && j != N2) {     //linker Rand
                            vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c + j*h2)*(h2) -
h1*h2*fun(a+i*h1, c+j*h2)/2;
                        }
                        else if(i == N1 && j != 0 && j != N2) {     //rechter
Rand vec[i*(N2+1)+j] = (-1)* neumannx(a+ i*h1,c + j*h2)*(h2) - h1*h2*fun(a+i*h1,
c+j*h2)/2;
                        }
                        else if(j == 0 && i != 0 && i != N1) { // unter Rand
                            vec[i*(N2+1)+j] =  neumanny(a+ i*h1,c + j*h2)*(h1) -
h1*h2*fun(a+i*h1, c+j*h2)/2;
                        }
                        else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                            vec[i*(N2+1)+j] = (-1)* neumanny(a+ i*h1,c +
j*h2)*(h1) - h1*h2*fun(a+i*h1, c+j*h2)/2;
                        }
                        else {
                           if(i == 0 && j == 0) {
                                vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c +
j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 + neumanny(a+ i*h1,c + j*h2)*(h1/2);
                           }
                            else if(i == 0 && j == N2) {
                                vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c +
j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 - neumanny(a+ i*h1,c + j*h2)*(h1/2);
                            }
                            else if(i == N1 && j == N2) {
                                vec[i*(N2+1)+j] = (-1)* neumannx(a+ i*h1,c +
j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 - neumanny(a+ i*h1,c + j*h2)*(h1/2);
                            }
                            else if(i == N1 && j == 0) {
                                    vec[i*(N2+1)+j] =  (-1)*neumannx(a+ i*h1,c +
j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 + neumanny(a+ i*h1,c + j*h2)*(h1/2);
                            }
                        }
                    }
                    else vec[i*(N2+1)+j] = (-1) * h1 * h2 * fun(a + i*h1, c +
j*h2);
                }
            }
        }
        for(long i = 0; i < N1+1; i++) {     //
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) {
                    //Hier befinde ich mich am Rand!\n
                    if(i == 0 && j != 0 && j != N2) {     //linker Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/h1;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h1/h2) -
h2/h1;
                    }
                    else if(i == N1 && j != 0 && j != N2) {     //rechter Rand
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/h1;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h1/h2) -
h2/h1;
                    }
                    else if(j == 0 && i != 0 && i != N1) { // unter Rand
                         mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/h2;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h2/h1) -
h1/h2;
                    }
                    else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/h2;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h2/h1) -
h1/h2;
                    }
                    else {
                       if(i == 0 && j == 0) {   //linkes unteres Eck

                            mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                           // mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)*
(h1/h2+h2/h1);
                       }
                        else if(i == 0 && j == N2) {  //linkes oberes Eck
                             mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)*
(h1/h2+h2/h1);
                        }
                        else if(i == N1 && j == N2) {  //rechtes oberes Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)*
(h1/h2+h2/h1);
                        }
                        else if(i == N1 && j == 0) {  //rechtes unteres Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)*
(h1/h2+h2/h1);
                        }
                    }
                }
                else {
                    mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/h1;
                    mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/h1;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/h2;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/h2;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)* (h1/h2+h2/h1);
                }
            }
        }
        //l�sen
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(N1*N2);
        lscg.setTolerance(1e-12);
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        return ergebnis;
}
*/
