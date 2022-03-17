#ifndef SURFACE_H_INCLUDED
#define SURFACE_H_INCLUDED

#include <iostream>
#include "vec.h"

class Surface{
    private:

    vec *vertex;
    int numbers;

    public:

    void init(int length);
	void import_vertex(vec* dat);
	double get_surface(vec n);

};

void tetrahedron_circumcenter(
        // In:
        const double a[3],
        const double b[3],
        const double c[3],
        const double d[3],
        // Out:
        double circumcenter[3]);

#endif // SURFACE_H_INCLUDED
