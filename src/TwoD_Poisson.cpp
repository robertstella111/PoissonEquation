#include <iostream>
#include <fstream>
#include "../eigen-3.4.0/Eigen/IterativeLinearSolvers"

#include "TwoD_Poisson.h"

using namespace std;
using namespace Eigen;

 TWOD::TWOD(long N1, long N2, double a, double b, double c, double d, double (*fun)(double, double),
            double (*randx) (double, double), double (*randy)(double, double)) {
    //Für Neumann
    (*this).rand = 2;
    neumannx = randx;
    neumanny = randy;
    (*this).N1 = N1;
    (*this).N2 = N2;
    (*this).a = a;
    (*this).b = b;
    (*this).c = c;
    (*this).d = d;
    h2 = (d-c)/N2;
    h1 = (b-a)/N1;
    (*this).fun = fun;
    solvertype = 1;
}
 TWOD::TWOD(long N1, long N2, double a, double b, double c, double d, double (*fun)(double, double),
            double (*rand) (double, double)) {
    //Für Dirichlet
     (*this).rand = 1;
    dirichlet = rand;
    (*this).N1 = N1;
    (*this).N2 = N2;
    (*this).a = a;
    (*this).b = b;
    (*this).c = c;
    (*this).d = d;
    h2 = (d-c)/N2;
    h1 = (b-a)/N1;
    (*this).fun = fun;
    solvertype = 1;
}

Eigen::VectorXd TWOD::solve() {
    if(solvertype == 1) return solve2();
    else return solve1();
}

void TWOD::solvewitheFV() {
    solvertype = 2;
    return;
}
void TWOD::solvewithFD() {
    solvertype = 1;
    return;
}


Eigen::VectorXd TWOD::solve1() {
    SparseMatrix<double> mat((N1+1)*(N2+1), (N1+1)*(N2+1));
    Eigen::VectorXd vec((N1+1)*(N2+1));
    Eigen::VectorXd ergebnis((N1+1)*(N2+1));
    vector<Triplet<double> > tr;
    if(rand == 1) {
        //Dirichlet Randwerte
        //Für alle folgenden Schleifen, i ist für die x-Achse, j ist für die y-Achse
        for(long i = 0; i < N1+1; i++) {
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) vec[i*(N2+1)+j] = dirichlet(a+i*h1, c+ j*h2);
                else vec[i*(N2+1)+j] = (-1) * h1 * h2 * fun(a + i*h1, c + j*h2);
            }
        }
        for(long i = 0; i < N1+1; i++) {     //
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) mat.insert(i*(N2+1)+j,i*(N2+1)+j) = 1;
                else {
                    mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/h1;
                    mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/h1;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/h2;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/h2;
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)* (h1/h2+h2/h1);
                }
            }
        }
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations((N1+1)*(N2+1));
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        return ergebnis;
    }
    else {
        SparseMatrix<double> mat((N1+1)*(N2+1), (N1+1)*(N2+1));
        Eigen::VectorXd vec((N1+1)*(N2+1));
        Eigen::VectorXd ergebnis((N1+1)*(N2+1));
        for(long i = 0; i < N1+1; i++) {     //Vector Initialisieren
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) {
                    if(i == 0 && j != 0 && j != N2) {     //linker Rand
                        vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c + j*h2)*(h2) - h1*h2*fun(a+i*h1, c+j*h2)/2;
                    }
                    else if(i == N1 && j != 0 && j != N2) {     //rechter Rand
                        vec[i*(N2+1)+j] = (-1)* neumannx(a+ i*h1,c + j*h2)*(h2) - h1*h2*fun(a+i*h1, c+j*h2)/2;
                    }
                    else if(j == 0 && i != 0 && i != N1) { // unter Rand
                        vec[i*(N2+1)+j] =  neumanny(a+ i*h1,c + j*h2)*(h1) - h1*h2*fun(a+i*h1, c+j*h2)/2;
                    }
                    else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                        vec[i*(N2+1)+j] = (-1)* neumanny(a+ i*h1,c + j*h2)*(h1) - h1*h2*fun(a+i*h1, c+j*h2)/2;
                    }
                    else {
                       if(i == 0 && j == 0) {
                            vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c + j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 + neumanny(a+ i*h1,c + j*h2)*(h1/2);
                       }
                        else if(i == 0 && j == N2) {
                            vec[i*(N2+1)+j] =  neumannx(a+ i*h1,c + j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 - neumanny(a+ i*h1,c + j*h2)*(h1/2);
                        }
                        else if(i == N1 && j == N2) {
                            vec[i*(N2+1)+j] = (-1)* neumannx(a+ i*h1,c + j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 - neumanny(a+ i*h1,c + j*h2)*(h1/2);
                        }
                        else if(i == N1 && j == 0) {
                                vec[i*(N2+1)+j] =  (-1)*neumannx(a+ i*h1,c + j*h2)*(h2/2) - h1*h2*fun(a+i*h1, c+j*h2)/4 + neumanny(a+ i*h1,c + j*h2)*(h1/2);
                        }
                    }
                }
                else vec[i*(N2+1)+j] = (-1) * h1 * h2 * fun(a + i*h1, c + j*h2);
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
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h1/h2) - h2/h1;
                    }
                    else if(i == N1 && j != 0 && j != N2) {     //rechter Rand
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/h1;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h1/h2) - h2/h1;
                    }
                    else if(j == 0 && i != 0 && i != N1) { // unter Rand
                         mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/h2;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h2/h1) - h1/h2;
                    }
                    else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/h2;
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h2/h1) - h1/h2;
                    }
                    else {
                       if(i == 0 && j == 0) {   //linkes unteres Eck

                            mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                           // mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1)* (h1/h2+h2/h1);
                       }
                        else if(i == 0 && j == N2) {  //linkes oberes Eck
                             mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)* (h1/h2+h2/h1);
                        }
                        else if(i == N1 && j == N2) {  //rechtes oberes Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)* (h1/h2+h2/h1);
                        }
                        else if(i == N1 && j == 0) {  //rechtes unteres Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = h2/(h1*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = h1/(h2*2);
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-1/2)* (h1/h2+h2/h1);
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
        //lösen
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(N1*N2);
        lscg.setTolerance(1e-12);
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        return ergebnis;
    }
}

Eigen::VectorXd TWOD::solve2() {
    SparseMatrix<double> mat((N1+1)*(N2+1), (N1+1)*(N2+1));
    Eigen::VectorXd vec((N1+1)*(N2+1));
    Eigen::VectorXd ergebnis((N1+1)*(N2+1));
    if(rand == 1) {
        //Dirichlet Randwerte
        //Für alle folgenden Schleifen, i ist für die x-Achse, j ist für die y-Achse
        for(long i = 0; i < N1+1; i++) {
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) vec[i*(N2+1)+j] = dirichlet(a+i*h1, c+ j*h2);
                else vec[i*(N2+1)+j] = (-1) * fun(a + i*h1, c + j*h2);
            }
        }
        for(long i = 0; i < N1+1; i++) {     //
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) mat.insert(i*(N2+1)+j,i*(N2+1)+j) = 1;
                else {
                    mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (1/pow(h1,2));
                    mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (1/pow(h1,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (1/pow(h2,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (1/pow(h2,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                }
            }
        }
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations((N1+1)*(N2+1));
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        return ergebnis;
    }
    else {
        //Neumann Randbedingungen
        SparseMatrix<double> mat((N1+1)*(N2+1), (N1+1)*(N2+1));
        Eigen::VectorXd vec((N1+1)*(N2+1));
        Eigen::VectorXd ergebnis((N1+1)*(N2+1));
        for(long i = 0; i < N1+1; i++) {     //Vector Initialisieren
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) {
                    if(i == 0 && j != 0 && j != N2) {     //linker Rand
                        vec[i*(N2+1)+j] = 2*neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2);
                    }
                    else if(i == N1 && j != 0 && j != N2) {     //rechter Rand
                        vec[i*(N2+1)+j] = -2*neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2);
                    }
                    else if(j == 0 && i != 0 && i != N1) { // unter Rand
                        vec[i*(N2+1)+j] = 2*neumanny(a+ i*h1,c + j*h2)/h2 - fun(a + i*h1, c + j*h2);
                    }
                    else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                        vec[i*(N2+1)+j] = -2*neumanny(a+ i*h1,c + j*h2)/h2 - fun(a + i*h1, c + j*h2);
                    }
                    else {
                        if(i == 0 && j == 0) {
                            vec[i*(N2+1)+j] =  2*neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2) +2* neumanny(a+ i*h1,c + j*h2)/h2;
                        }
                        else if(i == 0 && j == N2) {
                            vec[i*(N2+1)+j] =  2* neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2) -2* neumanny(a+ i*h1,c + j*h2)/h2;
                        }
                        else if(i == N1 && j == N2) {
                            vec[i*(N2+1)+j] =  -2*neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2) - 2*neumanny(a+ i*h1,c + j*h2)/h2;
                        }
                        else if(i == N1 && j == 0) {
                            vec[i*(N2+1)+j] =   -2*neumannx(a+ i*h1,c + j*h2)/h1 - fun(a + i*h1, c + j*h2) + 2*neumanny(a+ i*h1,c + j*h2)/h2;
                        }
                    }
                }
                else vec[i*(N2+1)+j] = (-1) * fun(a + i*h1, c + j*h2);
            }
        }
        for(long i = 0; i < N1+1; i++) {     //
            for(long j = 0; j < N2+1; j++) {
                if(i == 0 || i == N1 || j == 0 || j == N2) {
                    //Hier befinde ich mich am Rand!\n
                    if(i == 0 && j != 0 && j != N2) {     //linker Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (2/pow(h1,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (1/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (1/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                    }
                    else if(i == N1 && j != 0 && j != N2) {     //rechter Rand
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (2/pow(h1,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (1/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (1/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                    }
                    else if(j == 0 && i != 0 && i != N1) { // unter Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (1/pow(h1,2));
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (1/pow(h1,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (2/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                    }
                    else if(j == N2 && i != 0 && i != N1) { //oberer Rand
                        mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (1/pow(h1,2));
                        mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (1/pow(h1,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (2/pow(h2,2));
                        mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                    }
                    else {
                        if(i == 0 && j == 0) {   //linkes unteres Eck
                            mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (2/pow(h1,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (2/pow(h2,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - (2/pow(h1,2));
                        }
                        else if(i == 0 && j == N2) {  //linkes oberes Eck
                            mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (2/pow(h1,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (2/pow(h2,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - (2/pow(h1,2));
                        }
                        else if(i == N1 && j == N2) {  //rechtes oberes Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (2/pow(h1,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (2/pow(h2,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - (2/pow(h1,2));
                        }
                        else if(i == N1 && j == 0) {  //rechtes unteres Eck
                            mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (2/pow(h1,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (2/pow(h2,2));
                            mat.insert(i*(N2+1)+j,i*(N2+1)+j) = -(2/pow(h2,2)) - (2/pow(h1,2));
                        }
                    }
                }
                else {
                    mat.insert(i*(N2+1)+j,(i+1)*(N2+1)+j) = (1/pow(h1,2));
                    mat.insert(i*(N2+1)+j,(i-1)*(N2+1)+j) = (1/pow(h1,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j+1) = (1/pow(h2,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j-1) = (1/pow(h2,2));
                    mat.insert(i*(N2+1)+j,i*(N2+1)+j) = (-2)*(1/pow(h2,2)) - 2 * (1/pow(h1,2));
                }
            }
        }
        //lösen
        BiCGSTAB<SparseMatrix<double> > lscg;
        lscg.setMaxIterations(N1*N2);
        lscg.setTolerance(1e-12);
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        return ergebnis;
    }
}


 void TWOD::write(Eigen::VectorXd vec,bool print) {
    fstream myFile;
    myFile.open("TwoD.txt", ios::out);
    if(!myFile) {
        cout << "Nicht geoeffnet!\n";
    }
    else {
        myFile << N1+1 << "|" << N2 + 1 << endl;
        for(long i = 0; i < N1+1; i++) {
            for(long j = 0; j < N2+1; j++) {
                myFile <<  a+i*h1 << "|" << c+j*h2 << "|" << vec[i*(N2+1)+j] << endl;
            }
        }
        myFile.close();
        /*
        if(print == true) {
            system("TwoD.py");
            system("start Potential2D.pdf");
        }
        */
    }
 }
