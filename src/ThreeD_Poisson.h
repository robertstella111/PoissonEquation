#ifndef THREED_POISSON_H_INCLUDED
#define THREED_POISSON_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <string>
#define epsilon 8.854e-12

class ThreeD{
    private:
        //defaultmässig werden alle Randbedingungen auf 0 gesetzt
        double( *fun)(double, double, double);      //Ladunsverteilung
        double (*dirichlet) (double, double, double);   //Rückgabe sind die Dirichlet Randwertbedinungen
        double (*neumannx) (double, double, double);     //Ableitung des Potentials nach x, Paramter y
        double (*neumanny) (double, double, double);    //Ableitung des Potentials nach y, Paramter x
        double (*neumannz) (double, double, double);
        long N1, N2, N3;
        double h1, h2, h3, a, b, c, d, e, f; //Gitterbreite, Intervall der Funktion
        Eigen::VectorXd solve2(); //solver finite Differenzen
        Eigen::VectorXd solve1();  //solver Finite Volumen
        int solvertype; // 1...FD, 2....FV

    public:
        Eigen::VectorXd solve();    //Lösen des Gleichungssystems und als Rückgabe den Ergebnisvektor
        Eigen::VectorXd solvetest();
        ThreeD(long N1, long N2, long N3, double a, double b, double c, double d, double e, double f,
                double (*fun)(double, double, double), double (*rand) (double, double, double));
        void solvewitheFV();
        void solvewithFD();

};


#endif // THREED_POISSON_H_INCLUDED
