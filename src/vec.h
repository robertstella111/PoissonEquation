#include <string>
#include <fstream>
#include <iostream>
#ifndef MY_VEC
#define MY_VEC

class vec;
using namespace std;
class vec
{
private:
	double x[3];
friend class VoronoiDiagram;
friend class Matrix12;
friend class Diff_gen;
//friend cell;
public:
	void import(double *);
	void clean();					// reset value to zero
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec operator/(const double&);
	double operator*(const vec&);
	vec operator^(const vec&);		// for cross product
	double operator %(const vec&);
	vec & operator=(const vec&);
	vec & operator=(double*);		// can replace import
	double norm();
	bool operator==(const vec&);    //Zum Vergleichen von vektoren

	//debug
    friend ostream &operator<<(ostream& os, const vec& vec) {
        os << "x= " << vec.x[0] << " y= " << vec.x[1] << " z= " << vec.x[2] << endl;
        return os;
    }
	void print();


};

#endif
