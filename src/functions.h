#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <math.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "read_mesh.h"
#include "vec.h"
 using namespace std;
 class Data;

long index_vec(vector<long> feld, long wert);
bool compare_points(vec a, vec b);
bool compare_koord(double a, double b);


class Matrix12{
    private:
    double x[9];

    public:
    double getkoeff(long zeile, long spalte);

    void init(long zeile, long spalte, double value);
    vec solve(double* ergebnis);


};

class List{
private:
    Data *list = NULL;
    long length;
    long last;
    long first;
public:
    List();
   void insert_front(bool tetraeder, long value);
   void insert_back(bool tetraeder, long value);
   void loeschen();
   long getSize();
   long getLast();
   long getFirst();
    bool getBool(long index);
    bool getLastBool();
    long getIndex(long index);
    void anpassen(vector<vec> circumcenter);

};

class Data{
private:
    bool tetraeder;
    long index;
    Data *next;
public:
    Data();
friend class List;

};


#endif // FUNCTIONS_H_INCLUDED
