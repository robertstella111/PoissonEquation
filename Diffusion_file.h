#ifndef DIFFUSION_FILE_H_INCLUDED
#define DIFFUSION_FILE_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include "functions.h"
#include <string>
#include "read_mesh.h"
#define epsilon 8.854e-12
#define elementar 1.602e-19
#include "vec.h"
#include <algorithm>

class Diff_gen{
    private:
        //defaultmässig werden alle Randbedingungen auf 0 gesetzt
    Polyhedronlist2 polyhedronlist;
    vector<long> pointlist;



    public:
        Eigen::VectorXd solve();    //Lösen des Gleichungssystems und als Rückgabe den Ergebnisvektor
        Eigen::VectorXd solve2();      //explizit euler
        Eigen::VectorXd solve3();
        Diff_gen(string filename);
        void test1();
        Eigen::VectorXd Poisson(Eigen::VectorXd Ladunsverteilung);

};


#endif // DIFFUSION_FILE_H_INCLUDED
