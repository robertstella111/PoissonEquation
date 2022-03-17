#ifndef TWOD_POISSON_H_INCLUDED
#define TWOD_POISSON_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <string>
#define epsilon 8.854e-12

class TWOD {
    private:
        //defaultmässig werden alle Randbedingungen auf 0 gesetzt
        double( *fun)(double, double);      //Ladunsverteilung
        double (*dirichlet) (double, double);   //Rückgabe sind die Dirichlet Randwertbedinungen
        double (*neumannx) (double, double);     //Ableitung des Potentials nach x, Paramter y
        double (*neumanny) (double, double);    //Ableitung des Potentials nach y, Paramter x
        long N1, N2;
        double h1, h2, a, b, c, d; //Gitterbreite, Intervall der Funktion
        long rand;  // 1... Dirichlet, 2...Neumann
         Eigen::VectorXd solve2(); //solver finite Differenzen
         Eigen::VectorXd solve1();  //solver Finite Volumen
         int solvertype; // 1...FD, 2....FV

    public:
        Eigen::VectorXd solve();    //Lösen des Gleichungssystems und als Rückgabe den Ergebnisvektor
        Eigen::VectorXd solveN();
         Eigen::VectorXd solveN1();
        TWOD(long N1, long N2, double a, double b, double c, double d, double (*fun)(double, double),
            double (*randx) (double, double), double (*randy)(double, double));
        TWOD(long N1, long N2, double a, double b, double c, double d, double (*fun)(double, double),
            double (*rand) (double, double));
        void solvewitheFV();
        void solvewithFD();
        void write(Eigen::VectorXd vec,bool print);

};

#endif // TWOD_POISSON_H_INCLUDED
