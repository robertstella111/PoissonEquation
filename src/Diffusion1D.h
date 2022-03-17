#ifndef DIFFUSION1D_H_INCLUDED
#define DIFFUSION1D_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"


class Diffusion1D {
private:
    double h, dt, a, b, t, D;
    long N, T;
    double (*verteilung) (double);
    long rand;  //1... Randwert ist konstant, 2.... Steigung am Rand linken Rand = 0

public:
    Eigen::VectorXd solve();
    Diffusion1D(double D, long N, long T, double a, double b, double t, double (*verteilung) (double));
    void setEdgeConst();
    void setoutFlux();
    void write(Eigen::VectorXd vec,bool print);
};




#endif // DIFFUSION1D_H_INCLUDED
