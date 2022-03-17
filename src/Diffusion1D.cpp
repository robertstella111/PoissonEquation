#include <iostream>
#include "../eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include <math.h>
#include <fstream>
#include "Diffusion1D.h"

using namespace std;
using namespace Eigen;



Diffusion1D::Diffusion1D(double D, long N, long T, double a, double b, double t, double (*verteilung) (double)) {
    (*this).verteilung = verteilung;
    (*this).N = N;
    (*this).T = T;
    (*this).a = a;
    (*this).b = b;
    (*this).t = t;
    (*this).D = D;
    dt = t/T;
    h = (b-a)/N;
    rand = 1;
}
void Diffusion1D::setEdgeConst() {
    rand = 1;
}
void Diffusion1D::setoutFlux() {
    rand = 2;
}
Eigen::VectorXd Diffusion1D::solve() {
    double lambda = (D*dt)/pow(h,2);
    double k1, k2;
    k1 = 1+2*lambda;
    k2 = (-1)*lambda;
    SparseMatrix<double> mat(N+1, N+1);
    Eigen::VectorXd vec(N+1);
    Eigen::VectorXd zwischen(N+1);
    Eigen::VectorXd ergebnis((N+1)*(T+1));
    if(rand == 1) {
        for(long i = 0; i < N+1; i++) {
            vec[i] = verteilung(a + i * h);
        }
        for(long i = 0; i < N+1; i++) {
            if(i == 0) { //erste Reihe
                mat.insert(0,0) = 1;
            }
            else if(i == N) {   //letzte Reihe
                mat.insert(i, i-1) = -1;
                mat.insert(i, i) = 1;
            }
            else {
                mat.insert(i,i) = k1;
                mat.insert(i,i+1) = k2;
                mat.insert(i,i-1) = k2;
            }
        }
        for(long i = 0; i < N+1; i++) {
            ergebnis[i] = vec[i];
        }
        vec[N] = 0;
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(N);
        lscg.compute(mat);
        for(long x = 1; x < T+1; x++) {
            zwischen = lscg.solve(vec);
            for(long i = 0; i < N+1; i++) {
                ergebnis[x*(N+1)+i] = zwischen[i];
            }
            vec = zwischen;
            vec[N] = 0;
        }
        return ergebnis;
    }
    else {
        for(long i = 0; i < N+1; i++) {
            vec[i] = verteilung(a + i*h);
        }
        for(long i = 0; i < N+1; i++) {
            if(i == 0) { //erste Reihe
                mat.insert(i, i+1) = 1;
                mat.insert(i, i) = -1;
            }
            else if(i == N) {   //letzte Reihe
                mat.insert(i, i-1) = -1;
                mat.insert(i, i) = 1;
            }
            else {
                mat.insert(i,i) = k1;
                mat.insert(i,i+1) = k2;
                mat.insert(i,i-1) = k2;
            }
        }
        for(long i = 0; i < N+1; i++) {
            ergebnis[i] = vec[i];
        }
        vec[0] = 0;
        vec[N] = 0;
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(N);
        lscg.compute(mat);
        for(long x = 1; x < T+1; x++) {
            zwischen = lscg.solve(vec);
            for(long i = 0; i < N+1; i++) {
                ergebnis[x*(N+1)+i] = zwischen[i];
            }
            vec = zwischen;
            vec[N] = 0;
            vec[0] = 0;
        }
        return ergebnis;
    }
}


void Diffusion1D::write(Eigen::VectorXd vec,bool print) {
    fstream myFile;
    myFile.open("Diff.txt", ios::out);
    if(!myFile) {
        cout << "Nicht geoeffnet!\n";
    }
    else {
        myFile << T+1 << "|" << N + 1 << endl;
       for(long i = 0; i < T+1; i++) {
            for(long j = 0; j < N+1; j++) {
                myFile <<  a+j*h << "|" << vec[i*(N+1)+j] << endl;
            }
        }
        if(print == true) {
            system("Diff.py");
        }
        myFile.close();
    }
}






