#ifndef RUNGE_KUTTA_H_INCLUDED
#define RUNGE_KUTTA_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <string>
#define epsilon 8.854e-12

Eigen::VectorXd runge(double (*fun) (double,double), double y0, double t0, double end, long N);
void funktion();

class Points{
private:
public:
    std::string id;
    Eigen::VectorXd points;
    long anzahl;
};

#endif // RUNGE_KUTTA_H_INCLUDED
