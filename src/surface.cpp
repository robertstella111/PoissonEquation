#include "surface.h"


void Surface::init(int length) {
    numbers = length;
}
void Surface::import_vertex(vec* dat) {
    vertex = dat;
}
double Surface::get_surface(vec n) {
    double surface = 0;
    for(long j = 0; j < numbers; j++) {
        if(j == numbers-1)   surface += (vertex[j] ^ vertex[0])*n;
        else surface += (vertex[j] ^ vertex[j+1])*n;
    }
    if(surface < 0) surface *= -1;
    return surface/2;
}




void tetrahedron_circumcenter(
        // In:
        const double a[3],
        const double b[3],
        const double c[3],
        const double d[3],
        // Out:
        double circumcenter[3])
{
    double denominator;


    // Use coordinates relative to point `a' of the tetrahedron.

    // ba = b - a
    double ba_x = b[0] - a[0];
    double ba_y = b[1] - a[1];
    double ba_z = b[2] - a[2];

    // ca = c - a
    double ca_x = c[0] - a[0];
    double ca_y = c[1] - a[1];
    double ca_z = c[2] - a[2];

    // da = d - a
    double da_x = d[0] - a[0];
    double da_y = d[1] - a[1];
    double da_z = d[2] - a[2];

    // Squares of lengths of the edges incident to `a'.
    double len_ba = ba_x * ba_x + ba_y * ba_y + ba_z * ba_z;
    double len_ca = ca_x * ca_x + ca_y * ca_y + ca_z * ca_z;
    double len_da = da_x * da_x + da_y * da_y + da_z * da_z;

    // Cross products of these edges.

    // c cross d
    double cross_cd_x = ca_y * da_z - da_y * ca_z;
    double cross_cd_y = ca_z * da_x - da_z * ca_x;
    double cross_cd_z = ca_x * da_y - da_x * ca_y;
    //cout << "Ca_y" << ca_y << "Da_z "<<da_z << "Da_y "<< da_y << "Ca_z " << ca_z   << endl;

    // d cross b
    double cross_db_x = da_y * ba_z - ba_y * da_z;
    double cross_db_y = da_z * ba_x - ba_z * da_x;
    double cross_db_z = da_x * ba_y - ba_x * da_y;

    // b cross c
    double cross_bc_x = ba_y * ca_z - ca_y * ba_z;
    double cross_bc_y = ba_z * ca_x - ca_z * ba_x;
    double cross_bc_z = ba_x * ca_y - ca_x * ba_y;

    // Calculate the denominator of the formula.
    denominator = 0.5 / (ba_x * cross_cd_x + ba_y * cross_cd_y + ba_z * cross_cd_z);
   // cout << "denominator:" << denominator << endl;
    // Calculate offset (from `a') of circumcenter.
    double circ_x = (len_ba * cross_cd_x + len_ca * cross_db_x + len_da * cross_bc_x) * denominator;
    double circ_y = (len_ba * cross_cd_y + len_ca * cross_db_y + len_da * cross_bc_y) * denominator;
    double circ_z = (len_ba * cross_cd_z + len_ca * cross_db_z + len_da * cross_bc_z) * denominator;

    circumcenter[0] = circ_x;
    circumcenter[1] = circ_y;
    circumcenter[2] = circ_z;

        // To interpolate a linear function at the circumcenter, define a
        // coordinate system with a xi-axis directed from `a' to `b',
        // an eta-axis directed from `a' to `c', and a zeta-axis directed
        // from `a' to `d'.  The values for xi, eta, and zeta are computed
        // by Cramer's Rule for solving systems of linear equations.
        /*
        denominator *= 2.0;
        *xi   = (circ_x * cross_cd_x + circ_y * cross_cd_y + circ_z * cross_cd_z) * denominator;
        *eta  = (circ_x * cross_db_x + circ_y * cross_db_y + circ_z * cross_db_z) * denominator;
        *zeta = (circ_x * cross_bc_x + circ_y * cross_bc_y + circ_z * cross_bc_z) * denominator;
*/
}

