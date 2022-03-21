#include "../include/functions.h"




vector<double> sub(vector<double> a, vector<double> b) {
vector<double> diff(3);
    for(long k = 0; k < 3; k++) {
        diff[k] = a[k]-b[k];
    }
    return diff;
}



    double Matrix12::getkoeff(long zeile, long spalte) {
        return x[(zeile-1)*3+(spalte-1)];
    }


void Matrix12::init(long zeile, long spalte, double value) {
     x[(zeile-1)*3+(spalte-1)] = value;
}


    vec Matrix12::solve(double *ergebnis) {
        double determinante;
        vec result;
        double zwischen[3];
        determinante =  getkoeff(1,1)*getkoeff(2,2)*getkoeff(3,3) +
                        getkoeff(3,1)*getkoeff(1,2)*getkoeff(2,3) +
                        getkoeff(2,1)*getkoeff(3,2)*getkoeff(1,3) -
                        getkoeff(1,1)*getkoeff(3,2)*getkoeff(2,3) -
                        getkoeff(3,1)*getkoeff(2,2)*getkoeff(1,3) -
                        getkoeff(2,1)*getkoeff(1,2)*getkoeff(3,3);
        zwischen[0] =   (ergebnis[0]*getkoeff(2,2)*getkoeff(3,3) +
                        ergebnis[2]*getkoeff(1,2)*getkoeff(2,3) +
                        ergebnis[1]*getkoeff(3,2)*getkoeff(1,3) -
                        ergebnis[0]*getkoeff(3,2)*getkoeff(2,3) -
                        ergebnis[2]*getkoeff(2,2)*getkoeff(1,3) -
                        ergebnis[1]*getkoeff(1,2)*getkoeff(3,3))/determinante;



        zwischen[1] =  (getkoeff(1,1)*ergebnis[1]*getkoeff(3,3) +
                        getkoeff(3,1)*ergebnis[0]*getkoeff(2,3) +
                        getkoeff(2,1)*ergebnis[2]*getkoeff(1,3) -
                        getkoeff(1,1)*ergebnis[2]*getkoeff(2,3) -
                        getkoeff(3,1)*ergebnis[1]*getkoeff(1,3) -
                        getkoeff(2,1)*ergebnis[0]*getkoeff(3,3))/determinante;

        zwischen[2] =  (getkoeff(1,1)*getkoeff(2,2)*ergebnis[2] +
                        getkoeff(3,1)*getkoeff(1,2)*ergebnis[1] +
                        getkoeff(2,1)*getkoeff(3,2)*ergebnis[0] -
                        getkoeff(1,1)*getkoeff(3,2)*ergebnis[1] -
                        getkoeff(3,1)*getkoeff(2,2)*ergebnis[0] -
                        getkoeff(2,1)*getkoeff(1,2)*ergebnis[2])/determinante;
        result.import(zwischen);
       //  cout << "Feherx  " << ((getkoeff(1,1)*zwischen[0] + getkoeff(1,2)*zwischen[1]+ getkoeff(1,3)*zwischen[2]) - ergebnis[0])*100/ergebnis[0] << endl ;
         // cout << "Feher y " << ((getkoeff(2,1)*zwischen[0] + getkoeff(2,2)*zwischen[1]+ getkoeff(2,3)*zwischen[2]) - ergebnis[1])*100/ergebnis[1] << endl ;
           //cout << "Feher z " << ((getkoeff(3,1)*zwischen[0] + getkoeff(3,2)*zwischen[1]+ getkoeff(3,3)*zwischen[2]) - ergebnis[2])*100/ergebnis[2] << endl ;
        return result;
    }


/*
void punktabstand(vector<double> basis, vector<double> p1, vector<double> p2, vector<double> p3,) {
    //Ob die Punktabstaende stimme
}
*/

long index_vec(vector<long> feld, long wert) {
    for(unsigned long k = 0; k < feld.size(); k++) {
        if(feld[k] == wert) return k;
    }
    return -1;
}
List::List() {
    length = 0;
    last = 0;
    first = 0;
}
void List::insert_front(bool tetraeder, long value) {
    first = value;
    Data *zwischen;
    zwischen = new Data[1];
    (*zwischen).tetraeder = tetraeder;
    (*zwischen).index = value;
    (*zwischen).next = list;
    list = zwischen;
    length++;
}

long List::getSize() {
    return length;
}

 void List::loeschen() {
    last = 0;
    first = 0;
    Data *zwischen;
    length = 0;
    if(list != NULL) {
        zwischen = list;
        while(list != NULL) {
            list = (*list).next;
            delete zwischen;
            zwischen = list;
        }
    }
    list = NULL;
 }


void List::insert_back(bool tetraeder, long value) {
Data *zwischen;
 last = value;
 if(list == NULL) {
        list = new Data();
        (*list).tetraeder = tetraeder;
        first = value;
        (*list).index = value;
    }
    else {
        zwischen = list;
        while((*zwischen).next != NULL) zwischen = (*zwischen).next;
        (*zwischen).next = new Data();
        zwischen = (*zwischen).next;
        (*zwischen).tetraeder = tetraeder;
        (*zwischen).index = value;
        (*zwischen).next = NULL;
    }
    length++;
}

Data::Data() {
    next = NULL;
}

long List::getLast() {
    return last;
}
long List::getFirst() {
    return first;
}

 bool compare_points(vec a, vec b) {
     double eps = 1e-12;
    if((a-b).norm() <= eps) return true;
    else return false;
 }



bool compare_koord(double a, double b) {
    double eps = 1e-14;
    double c = a-b;
    if(c < 0)  {
        if(-c <= eps) return true;
        else return false;
    }
    else {
       if(c <= eps) return true;
        else return false;
    }
}


bool List::getBool(long index) {
    if(index < length) {
        Data *zwischen;
        zwischen = list;
        for(long u = 0; u < index; u++) {
            zwischen = (*zwischen).next;
        }
        return (*zwischen).tetraeder;
    }
    else {
        cout << "Fehler" << endl;
        return false;
    }
}
long List::getIndex(long index) {
        if(index < length) {
        Data *zwischen;
        zwischen = list;
        for(long u = 0; u < index; u++) {
            zwischen = (*zwischen).next;
        }
        return (*zwischen).index;
    }
    else {
        cout << "Fehler2" << endl;
        return false;
    }
}

void List::anpassen(vector<vec> circumcenter) {
    Data *zwischen, *zwischen1;
    vec a, b;
    zwischen = list;
    /*
    for(long k = 0; k < length; k++) {
        cout <<  circumcenter[(*zwischen).index].x[0]<< "  " <<
        circumcenter[(*zwischen).index].x[1]<< "  " <<
        circumcenter[(*zwischen).index].x[2]<< "  " <<  endl;
        zwischen = (*zwischen).next;
    }
    zwischen = list;
    */
    for(long u = 0; u+1 < length; u++) {
        if((*zwischen).tetraeder && (*((*zwischen).next)).tetraeder) {
            a = circumcenter[(*zwischen).index];
            b = circumcenter[(*((*zwischen).next)).index];
            if(compare_points(a,b)) {
                zwischen1 = (*zwischen).next;
                (*zwischen).next = (*zwischen1).next;
                delete zwischen1;
                length--;
                u--;
            }
            else{
                zwischen = (*zwischen).next;
            }
        }
        else break;

    }
}


bool List::getLastBool() {
    Data *zwischen;
    zwischen = list;
    for(long u = 0; u + 1 < length; u++) {
        zwischen = (*zwischen).next;
    }
    return (*zwischen).tetraeder;

}
