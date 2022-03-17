#include "read_mesh.h"
#include "functions.h"
#include "surface.h"




PointList2::PointList2() {
    number = 0;
}


void PointList2::read_pointlist(string filename) {
    //Zum Einlesn der Punkte. Werden für die Tetraeder benötigt. Punkte werden als Vekotr abgespeichert
    fstream myFile;
    myFile.open(filename.c_str(), ios::in);
    if(!myFile) {
        cout << "Nicht geoeffnet!\n";
    }
    else {
        string read, read_tmp;
        long a;
        stringstream str_stm;
        vec v1;
        double input[3];
        vector<vec> v2;
        while(getline(myFile, read)) {
            if(read.find("Pointlists") != string::npos) {
                str_stm << read;
                getline(str_stm, read_tmp, '(');
                getline(str_stm, read_tmp, ')');
                stringstream(read_tmp) >> number;
                break;
            }
        }
        str_stm.str(string());
        str_stm.clear();
        for(long x = 0; x < number; x++) {
            while(getline(myFile, read)) {
                if(read.find("Pointlist(") != string::npos) {   //gefunden

                    str_stm << read;
                    getline(str_stm, read_tmp, '"');
                    getline(str_stm, read_tmp, '"');
                    name.push_back(read_tmp);
                    getline(str_stm, read_tmp, ',');
                    getline(str_stm, read_tmp, ',');
                    getline(str_stm, read_tmp, ',');
                    getline(str_stm, read_tmp, ')');
                    stringstream(read_tmp) >> a;
                    str_stm.str(string());
                    str_stm.clear();
                    break;
                }
            }
            //cout << a << endl;
            pointnumbers.push_back(a);
            getline(myFile, read);
            for(long i = 0; i < a; i++) {
                getline(myFile,read);
                //cout << read << endl;

                str_stm << read;
                getline(str_stm, read_tmp, ',');
                stringstream(read_tmp) >> input[0];
                getline(str_stm, read_tmp, ',');
                stringstream(read_tmp) >> input[1];
                getline(str_stm, read_tmp, ',');
                stringstream(read_tmp) >> input[2];
                v1.import(input);
                v2.push_back(v1);
                str_stm.str(string());
                str_stm.clear();
            }
            points.push_back(v2);
            v2.clear();
        }
        myFile.close();

    }
}

long PointList2::getnumber() {
    //derzeit unwichtig
    return number;
}

long PointList2::getlength(string search) {
    //Als Übergabeparameter die Bezeichnung der Punktgruppe. Gibt zurück wie viele Elemente diese
    //Gruppe besitzt
    long a;
    for(a = 0; a < number; a++) if(search.compare(name[a]) == 0) break;
    if(a >= number)return 0;
    return pointnumbers[a];
}



 vec PointList2::getPoint(string ref, long i) {
     //Gibt einen einzelnen Punkt als Vektor zurück
    long a;
    for(a = 0; a < number; a++) if(ref.compare(name[a]) == 0) break;
    if(a >= number)return vec();
    return points[a][i];
}

Polyhedronlist2::Polyhedronlist2() {
    number = 0;
    pointlist.read_pointlist("../Input/nin_in.deva");
}

bool Polyhedronlist2::neighbors(long index1, long index2, long mittelpunkt, long nachbar) {
    //Untersucht ob zwei Punkt eine gemeinsame Fläche besitzen
    //Die Punkt sind mittelpunkt und anchbar
    //index1 und index2 sind die Tetraedernummer die untersucht werden
    vector<long> list, list1;
    list = getPoints("simulationGrid",index1);
    list1 = getPoints("simulationGrid",index2);
    list1.erase(find(list1.begin(), list1.end(),mittelpunkt));
    list.erase(find(list.begin(), list.end(),mittelpunkt));
    list1.erase(find(list1.begin(), list1.end(),nachbar));
    list.erase(find(list.begin(), list.end(),nachbar));
    if(find(list.begin(),list.end(),list1[0]) == list.end() ||
        find(list.begin(),list.end(),list1[1]) == list.end()) {
            //Ein Element gleich
            return true;
    }
    else return false;
}
void Polyhedronlist2::read(string filename) {
    //liest Polyhedronlisten ein. Wichtig hierfür: Aus welchen Punkten besteht ein Tetraeder
    fstream myFile;
    myFile.open(filename.c_str(), ios::in);
    if(!myFile) {
        cout << "Nicht geoeffnet!\n";
    }
    else {
        string read, read_tmp;
        long a, anzahl, zwischen;
        vector<long> v1;
        vector<vector<long> > v2;
        //double b, c, d;
        stringstream str_stm;
        while(getline(myFile, read)) {
            if(read.find("Polyhedronlists") != string::npos) {
                str_stm << read;
                getline(str_stm, read_tmp, '(');
                getline(str_stm, read_tmp, ')');
                stringstream(read_tmp) >> number;

                break;
            }
        }
        str_stm.str(string());
        str_stm.clear();
        for(long x = 0; x < number; x++) {
            while(getline(myFile, read)) {
                if(read.find(" Polyhedronlist(") != string::npos) {   //gefunden
                    str_stm << read;
                    getline(str_stm, read_tmp, '"');
                    getline(str_stm, read_tmp, '"');
                    name.push_back(read_tmp);
                    getline(str_stm, read_tmp, ',');
                    getline(str_stm, read_tmp, ')');
                    stringstream(read_tmp) >> a;
                    length.push_back(a);
                    str_stm.str(string());
                    str_stm.clear();
                    break;
                }
            }
            getline(myFile, read);
            str_stm << read;
            getline(str_stm, read_tmp, '"');
            getline(str_stm, read_tmp, '"');
            ref.push_back(read_tmp);
            getline(myFile, read);
            str_stm.str(string());
            str_stm.clear();
            for(long i = 0; i < a; i++) {
                getline(myFile,read);
                if(i == 0) {
                    str_stm << read;
                    getline(str_stm, read_tmp, ',');
                    anzahl = 0;
                    while(stringstream(read_tmp) >> zwischen ) {
                        v1.push_back(zwischen);
                        anzahl++;
                        getline(str_stm, read_tmp, ',');
                    }
                    pointnumbers.push_back(anzahl);
                }
                else {
                    str_stm << read;
                    getline(str_stm, read_tmp, ',');
                    for(long j = 0; j < anzahl; j++) {
                        stringstream(read_tmp) >> zwischen;
                        v1.push_back(zwischen);
                        getline(str_stm, read_tmp, ',');
                    }
                }
                v2.push_back(v1);
                str_stm.str(string());
                str_stm.clear();
                v1.clear();
            }
            points.push_back(v2);
            v1.clear();     //Speicher wird nicht freigegeben da möglicherweise noch benötigt
            v2.clear();
        }
        myFile.close();
    }
}


vector<long> Polyhedronlist2::getPoints(string search, long index) {
    //Gibt die Punkte als vector<long> zurück aus denen ein bestimmter Tetraeder besteht. Herbei wird davon ausgegangen das alle
    //im simulationGrid vorkommen
    vector<long> ergebnis;
    long a;
    for(a = 0; a < number; a++) {
        if(search.compare(name[a]) == 0) break;
    }
    if(a >= number) return ergebnis;
    for(long x = 1; x < pointnumbers[a]; x++) {
        ergebnis.push_back(points[a][index][x]);
    }
    return ergebnis;
}


long Polyhedronlist2::getlength(string search) {
    long a;
    for(a = 0; a < number; a++) {
        if(search.compare(name[a]) == 0) break;
    }
    if(a == number) return 0;
    return length[a];
}

bool Polyhedronlist2::exist(string search, long tetraedernr, long point) {
    //Es wird untersucht ob in einem bestimmten Tetraeder ein bestimmter Punkt vorkommt
    long a;
    for(a = 0; a < number; a++) {
        if(search.compare(name[a]) == 0) break;
    }
    if(a == number) return false;
    for(long x = 1; x < pointnumbers[a]; x++) {
        if(points[a][tetraedernr][x] == point) return true;
    }
    return false;
}

VoronoiDiagram::VoronoiDiagram() {
    pointlist.read_pointlist("../Input/nin_in.deva");
    polyhedronlist.read("../Input/nin_in.deva");
    circumcenter.resize(polyhedronlist.getlength("simulationGrid"));
    volume.resize(polyhedronlist.getlength("simulationGrid"));
    surface.resize(pointlist.getlength("simulationGrid"));
    for(unsigned long u = 0; u < surface.size(); u++) surface[u].resize(pointlist.getlength("simulationGrid"));
    normalvec.resize(pointlist.getlength("simulationGrid"));
    for(unsigned long u = 0; u < normalvec.size(); u++) normalvec[u].resize(pointlist.getlength("simulationGrid"));

}


void VoronoiDiagram::getVoronoiDiagram() {
    long semic1[] = {
              43,     44,     45,     46,     47,     48,     49,     50,     51,     52,     53,     // #0
              54,     55,     56,     57,     58,     59,     60,     61,     62,     63,     64,     // #1
              65,     66,     67,     68,     69,     70,     71,     72,     73,     74,     75,     // #2
              76,     77,     78,     79,     80,     81,     82,     83,     84,     85,     86,     // #3
              87,     88,     89,     90,     91,     92,     93,     94,     95,     96,     97,
              98,     99,    100,    101,    102,    103,    104,    105,    106,    107,     // #0
             108,    109,    110,    111,    112,    113,    114,    115,    116,    117,     // #1
             118,    119,    120,    121,    122,    123,    124,    125,    126,    127,     // #2
             128,    129,    130,    131,    132,    133,    134,    135,    136,    137,     // #3
             138,    139,    140,    141,    142,    143,    144,    145,    146,    147,     // #4
             148,    149,    150,    151,    152,    153,    154,    155,    156,    157,     // #5
             158,    159,    160,    161,    162,    163,    164,    165,    166,    167,     // #6
             168,    169,    170,    171,    172,    173,    174,    175,    176,    177,     // #7
             178,    179,    180,    181,    182,    183,    184,    185,    186,    187,     // #8
             188,    189,    190,    191,    192,    193,
             194,    195,    196,    197,    198,    199,    200,    201,    202,    203,     // #0
             204,    205,    206,    207,    208,    209,    210,    211,    212,    213,     // #1
             214,    215,    216,    217,    218,    219,    220,    221,    222,    223,     // #2
             224,    225,    226,    227,    228,    229,    230,    231,    232,    233,     // #3
             234,    235,    236,    237,    238,    239,    240,    241,    242,    243,     // #4
             244,    245
    };
    vec I,J,K,L, new_point, *vertex;
    double vec1[3];
    Matrix12 matrix;
    vector <long> tetraeder, list,list1, neighbors, copy, tetraeder2;
    List reihenfolge;
    Surface surf;
    double testvolume = 0, roh[3], gesamtvolume = 0, oberflaeche = 0;

    roh[0] = 0;
    roh[1] = 0;
    roh[2] = 0;
    for(unsigned long j = 0; j < sizeof(semic1)/sizeof(long); j++) {
        //Berechne für alle Tetraeder den circumcenter und speichere diesen ab. Index des Circumcenters ist der Index des Tetraeder zu welchem
        //dieser gehört
        I = pointlist.getPoint("simulationGrid", polyhedronlist.getPoints("simulationGrid",semic1[j])[0]);
        J = pointlist.getPoint("simulationGrid", polyhedronlist.getPoints("simulationGrid",semic1[j])[1]);
        K = pointlist.getPoint("simulationGrid", polyhedronlist.getPoints("simulationGrid",semic1[j])[2]);
        L = pointlist.getPoint("simulationGrid", polyhedronlist.getPoints("simulationGrid",semic1[j])[3]);


        vec1[0] = ((I+J)*0.5)*(J-I);//dotprodct(sum(2,I,J), sub(J,I));
        vec1[1] = ((I+K)*0.5)*(K-I);//dotprodct(sum(2,I,K), sub(K,I));
        vec1[2] = ((I+L)*0.5)*(L-I);//dotprodct(sum(2,I,L), sub(L,I));
        matrix.init(1,1,(J-I).x[0]);
        matrix.init(1,2,(J-I).x[1]);
        matrix.init(1,3,(J-I).x[2]);
        matrix.init(2,1,(K-I).x[0]);
        matrix.init(2,2,(K-I).x[1]);
        matrix.init(2,3,(K-I).x[2]);
        matrix.init(3,1,(L-I).x[0]);
        matrix.init(3,2,(L-I).x[1]);
        matrix.init(3,3,(L-I).x[2]);

        circumcenter[semic1[j]] = matrix.solve(vec1);
    }
    //Alle circumcenter initialisiert
    for(long i = 0; i < pointlist.getlength("simulationGrid"); i++) {
        //Alle Punkte die existieren durchgehen, für den dieser Punkte seine Nachbaren finden und in welchem Tetraeder dieser überall vorkommt
        //Die NAchbarn werden in neighbors gespeichert. Die Tetraedernummern in tetraeder
        for(unsigned long j = 0; j < sizeof(semic1)/sizeof(long); j++) {
            if(polyhedronlist.exist("simulationGrid",semic1[j],i) == true) {
                tetraeder.push_back(semic1[j]);
                list = polyhedronlist.getPoints("simulationGrid",semic1[j]);
                 if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != i) {
                    //nicht gefuden
                    neighbors.push_back(list[0]);
                }
                if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != i) {
                    //nicht gefuden
                    neighbors.push_back(list[1]);
                }
                if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != i) {
                    //nicht gefuden
                    neighbors.push_back(list[2]);
                }
                if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != i) {
                    //nicht gefuden
                    neighbors.push_back(list[3]);
                }
            }
        }
        if(neighbors.size() != 0) {
            //Dieser Mittelpunkt befindet sich in unserem Berecih, da er Nachbarn besitzt

            for(unsigned j = 0; j < neighbors.size(); j++) {
                //Für alle Nachbarpunkte vom Mittelpunkt die Fläche ermitteln
                //Mittelpunkt: i, neighbors: neighbor[j]
                //heruaslöschen der Tetraeder in welchen der nAchbarpunkt nicht vorkommt
                copy = tetraeder;
                for(unsigned long z = 0; z < copy.size(); z++) {
                    if(polyhedronlist.exist("simulationGrid",copy[z],neighbors[j]) == false) {
                        copy.erase(copy.begin() + z);
                        z--;
                    }
                }

                reihenfolge.insert_back(true, copy[0]);
                copy.erase(copy.begin());


                //Die ertse tetraedernr der reihenfolge zuweisen. Die klasse reihenfolge speichert ab, in welcher Reihenfolge die Punkte korrekkt
                //durchlaufen werden. Zuerst werden alle Punkt in der liste hinten dran gehängt. Falls dies nicht mehr möglich ist, so werden die Punkte
                //vorne hinzugefügt
                for(unsigned long t = 0; t < copy.size();t++) {
                    list = polyhedronlist.getPoints("simulationGrid",reihenfolge.getFirst());
                    list1 = polyhedronlist.getPoints("simulationGrid",copy[t]);
                    list1.erase(find(list1.begin(), list1.end(),i));
                    list.erase(find(list.begin(), list.end(),i));
                    list1.erase(find(list1.begin(), list1.end(),neighbors[j]));
                    list.erase(find(list.begin(), list.end(),neighbors[j]));
                    //Aus beiden Listen den Nachbarn und den Mittelpunkt gelöescht
                    if(find(list.begin(),list.end(),list1[0]) != list.end() ||
                        find(list.begin(),list.end(),list1[1]) != list.end()) {
                        //Ein Element gleich
                        reihenfolge.insert_front(true,copy[t]);
                        copy.erase(copy.begin() + t);
                        t=-1;

                    }
               }

               for(unsigned long t = 0; t < copy.size();t++) {
                    list = polyhedronlist.getPoints("simulationGrid",reihenfolge.getLast());
                    list1 = polyhedronlist.getPoints("simulationGrid",copy[t]);
                    list1.erase(find(list1.begin(), list1.end(),i));

                    list.erase(find(list.begin(), list.end(),i));
                    list1.erase(find(list1.begin(), list1.end(),neighbors[j]));
                    list.erase(find(list.begin(), list.end(),neighbors[j]));

                    //Aus beiden Listen den Nachbarn und den Mittelpunkt gelöescht
                    if(find(list.begin(),list.end(),list1[0]) != list.end() ||
                        find(list.begin(),list.end(),list1[1]) != list.end()) {
                        //Ein Element gleich
                        reihenfolge.insert_back(true,copy[t]);
                        copy.erase(copy.begin() + t);
                       t = -1;
                    }
               }
               //Nun werden die Flächen der jeweiligen Voronoi Zellen an den Rändern ermittelt. Für Randzellen entöang der x-Achse wird die Fläche
               // in roh[0], entlang der y-Achse in roh[1] und entlag der z-Achse in roh[3] abgespeichert
                vertex = new vec[3];
                surf.init(3);

               if((compare_koord(pointlist.getPoint("simulationGrid",i).x[0],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[0],0) == true) ||
                  (compare_koord(pointlist.getPoint("simulationGrid",i).x[0],6e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[0],6e-7) == true)) {
                    //beide Punkte sind am Rand x = 0
                    if(reihenfolge.getSize()==2) {
                        //if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],6e-7)) cout << "Hallo"<<endl;
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = -1;
                        I.x[1] = 0;
                        I.x[2] = 0;
                        roh[0] += surf.get_surface(I);
                    }
                    else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[1]) ) {
                            //entlang der y Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = -1;
                            I.x[1] = 0;
                            I.x[2] = 0;
                            roh[0] += surf.get_surface(I);
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[2])) {
                            //entlang der z Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = -1;
                            I.x[1] = 0;
                            I.x[2] = 0;
                            roh[0] += surf.get_surface(I);

                        }
                        else cout << "Fehler" << endl;
                    }
                    else {

                        long index1 = -1;
                        long index2 = -1;
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(x_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                if(index1 == -1) index1 = u;
                                else {
                                        index2 = u;
                                        break;
                                }
                            }
                        }
                        new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = -1;
                        I.x[1] = 0;
                        I.x[2] = 0;
                        roh[0] += surf.get_surface(I);
                       // cout <<"Surface  " << surf.get_surface(I) << endl;

                    }
                }


                if((compare_koord(pointlist.getPoint("simulationGrid",i).x[1],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[1],0) == true)||
                        (compare_koord(pointlist.getPoint("simulationGrid",i).x[1],2e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[1],2e-7) == true)) {
                    //beide Punkte sind am Rand y
                    if(reihenfolge.getSize()==2) {
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = 0;
                        I.x[1] = 1;
                        I.x[2] = 0;
                        roh[1] += surf.get_surface(I);
                      }
                      else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[0]) ) {
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = 0;
                            I.x[1] = 1;
                            I.x[2] = 0;
                            roh[1] += surf.get_surface(I);
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[2])) {


                             new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = 0;
                            I.x[1] = 1;
                            I.x[2] = 0;
                            roh[1] += surf.get_surface(I);

                        }
                        else cout << "Fehler" << endl;
                    }
                    else {
                        //cout << reihenfolge.getSize() << endl;
                        long index1 = -1;
                        long index2 = -1;
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(y_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                if(index1 == -1) index1 = u;
                                else {
                                    index2 = u;
                                    break;
                                }
                            }
                        }
                       // cout << "Mittelpunkt " << i << " Nachbar " << neighbors[j] << " " << reihenfolge.getSize() << endl;
                        new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = 0;
                        I.x[1] = 1;
                        I.x[2] = 0;
                        roh[1] += surf.get_surface(I);
                     //   cout << surf.get_surface(I) << endl;
                    }
                }

                if((compare_koord(pointlist.getPoint("simulationGrid",i).x[2],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[2],0) == true)||
                        (compare_koord(pointlist.getPoint("simulationGrid",i).x[2],2e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[2],2e-7) == true)) {
                    //beide Punkte sind am Rand x = 0
                    if(reihenfolge.getSize()==2) {
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = 0;
                        I.x[1] = 0;
                        I.x[2] = 1;
                        roh[2] += surf.get_surface(I);
                      }
                      else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[1]) ) {
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = 0;
                            I.x[1] = 0;
                            I.x[2] = 1;
                            roh[2] += surf.get_surface(I);
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[0])) {
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            vertex[0] = new_point;
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            vertex[1] = new_point;
                            vertex[2] = pointlist.getPoint("simulationGrid",i);
                            surf.import_vertex(vertex);
                            I.x[0] = 0;
                            I.x[1] = 0;
                            I.x[2] = 1;
                            roh[2] += surf.get_surface(I);
                        }
                        else cout << "Fehler" << endl;
                    }
                    else {
                            long index1 = -1;
                            long index2 = -1;
                            for(long u = 0; u < reihenfolge.getSize(); u++) {
                                if(z_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                    if(index1 == -1) index1 = u;
                                    else {
                                            index2 = u;
                                            break;
                                    }
                                }
                            }
                         new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        vertex[0] = new_point;
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        vertex[1] = new_point;
                        vertex[2] = pointlist.getPoint("simulationGrid",i);
                        surf.import_vertex(vertex);
                        I.x[0] = 0;
                        I.x[1] = 0;
                        I.x[2] = 1;
                        roh[2] += surf.get_surface(I);

                    }
                }

                //Es werden die Flächen zwischen  zwei Punkten ausgerechnet. Mit diesen Flächen kann schlussendlich auch das Volumen berechnet werden
               if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],0) == false &&
                  compare_koord(pointlist.getPoint("simulationGrid",i).x[0],6e-7) == false &&
                  compare_koord(pointlist.getPoint("simulationGrid",i).x[1],0) == false &&
                  compare_koord(pointlist.getPoint("simulationGrid",i).x[1],2e-7) == false &&
                  compare_koord(pointlist.getPoint("simulationGrid",i).x[2],0) == false &&
                  compare_koord(pointlist.getPoint("simulationGrid",i).x[2],2e-7) == false) {

                        vertex = new vec[reihenfolge.getSize()];
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            vertex[u] = circumcenter[reihenfolge.getIndex(u)];
                        }
                        surf.init(reihenfolge.getSize());
                        surf.import_vertex(vertex);
                        I   = (pointlist.getPoint("simulationGrid",neighbors[j]) -
                           pointlist.getPoint("simulationGrid",i));
                        I = I * (1/I.norm());
                        J = circumcenter[reihenfolge.getFirst()];
                        testvolume += (I*J)*surf.get_surface(I);
                        surface[i][neighbors[j]] = surf.get_surface(I);

                }
               else if((compare_koord(pointlist.getPoint("simulationGrid",i).x[0],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[0],0) == true)||
                        (compare_koord(pointlist.getPoint("simulationGrid",i).x[0],6e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[0],6e-7) == true)) {
                    //beide Punkte sind am Rand x = 0
                    if(reihenfolge.getSize()==2) {

                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                      }
                      else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[1]) ) {
                            //entlang der y Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = otherpoints[reihenfolge.getLast()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[2])) {
                            //entlang der z Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = otherpoints[reihenfolge.getLast()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else cout << "Fehler" << endl;
                    }
                    else {
                        long index1 = -1;
                        long index2 = -1;
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(z_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                if(index1 == -1) index1 = u;
                                else {
                                        index2 = u;
                                        break;
                                }
                            }
                        }
                        new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                    }
                    //Flächenberechnung
                    vertex = new vec[reihenfolge.getSize()];
                    for(long u = 0; u < reihenfolge.getSize(); u++) {
                        if(reihenfolge.getBool(u))vertex[u] = circumcenter[reihenfolge.getIndex(u)];
                        else vertex[u] = otherpoints[reihenfolge.getIndex(u)];
                    }
                    surf.init(reihenfolge.getSize());
                    surf.import_vertex(vertex);
                    I   = (pointlist.getPoint("simulationGrid",neighbors[j]) -
                        pointlist.getPoint("simulationGrid",i));
                    I = I * (1/I.norm());
                    J = circumcenter[reihenfolge.getFirst()];
                    testvolume += (I*J)*surf.get_surface(I);
                    surface[i][neighbors[j]] = surf.get_surface(I);

                }
               else if((compare_koord(pointlist.getPoint("simulationGrid",i).x[1],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[1],0) == true)||
                        (compare_koord(pointlist.getPoint("simulationGrid",i).x[1],2e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[1],2e-7) == true)) {
                    //beide Punkte sind am Rand y
                    if(reihenfolge.getSize()==2) {
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                      }
                      else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[0]) ) {
                        //entlang der y Kante
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = otherpoints[reihenfolge.getLast()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[2])) {
                            //entlang der z Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = otherpoints[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else cout << "Fehler" << endl;
                    }
                    else {
                        long index1 = -1;
                        long index2 = -1;
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(z_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                if(index1 == -1) index1 = u;
                                else {
                                        index2 = u;
                                        break;
                                }
                            }
                        }
                        new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                    }
                    //Flächenberechnung
                     vertex = new vec[reihenfolge.getSize()];
                    for(long u = 0; u < reihenfolge.getSize(); u++) {
                        if(reihenfolge.getBool(u))vertex[u] = circumcenter[reihenfolge.getIndex(u)];
                        else vertex[u] = otherpoints[reihenfolge.getIndex(u)];
                    }
                    surf.init(reihenfolge.getSize());
                    surf.import_vertex(vertex);
                    I   = (pointlist.getPoint("simulationGrid",neighbors[j]) -
                        pointlist.getPoint("simulationGrid",i));
                    I = I * (1/I.norm());
                    J = circumcenter[reihenfolge.getFirst()];
                    surface[i][neighbors[j]] = surf.get_surface(I);
                    testvolume += (I*J)*surf.get_surface(I);
                }
               else if((compare_koord(pointlist.getPoint("simulationGrid",i).x[2],0) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[2],0) == true)||
                        (compare_koord(pointlist.getPoint("simulationGrid",i).x[2],2e-7) == true &&
                        compare_koord(pointlist.getPoint("simulationGrid",neighbors[j]).x[2],2e-7) == true)) {
                    //beide Punkte sind am Rand x = 0
                    if(reihenfolge.getSize()==2) {
                        new_point = circumcenter[reihenfolge.getLast()];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getFirst()];
                        new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                      }
                      else if(reihenfolge.getSize() == 1) {
                        //entlang einer kante
                        if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[1]) ) {
                            //entlang der y Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = otherpoints[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],
                                         pointlist.getPoint("simulationGrid",neighbors[j]).x[0])) {
                            //entlang der z Kante
                            new_point = circumcenter[reihenfolge.getLast()];
                            new_point.x[0] = pointlist.getPoint("simulationGrid",i).x[0];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = otherpoints[reihenfolge.getLast()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                            new_point = circumcenter[reihenfolge.getFirst()];
                            new_point.x[2] = pointlist.getPoint("simulationGrid",i).x[2];
                            reihenfolge.insert_back(false,add_newPoint(new_point));
                        }
                        else cout << "Fehler" << endl;
                    }
                     else {

                         long index1 = -1;
                        long index2 = -1;
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(z_Ebene(reihenfolge.getIndex(u),i,neighbors[j])) {
                                if(index1 == -1) index1 = u;
                                else {
                                        index2 = u;
                                        break;
                                }
                            }
                        }
                        new_point = circumcenter[reihenfolge.getIndex(index1)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                        new_point = circumcenter[reihenfolge.getIndex(index2)];
                        new_point.x[1] = pointlist.getPoint("simulationGrid",i).x[1];
                        reihenfolge.insert_back(false,add_newPoint(new_point));
                    }
                    //Flächenberechnung
                     vertex = new vec[reihenfolge.getSize()];
                    for(long u = 0; u < reihenfolge.getSize(); u++) {
                        if(reihenfolge.getBool(u))vertex[u] = circumcenter[reihenfolge.getIndex(u)];
                        else vertex[u] = otherpoints[reihenfolge.getIndex(u)];
                    }
                    surf.init(reihenfolge.getSize());
                    surf.import_vertex(vertex);
                    I   = (pointlist.getPoint("simulationGrid",neighbors[j]) -
                        pointlist.getPoint("simulationGrid",i));
                    I = I * (1/I.norm());
                    J = circumcenter[reihenfolge.getFirst()];
                    surface[i][neighbors[j]] = surf.get_surface(I);
                    testvolume += (I*J)*surf.get_surface(I);
                }
                else {
                    if(reihenfolge.getSize() >= 3) {
                        vertex = new vec[reihenfolge.getSize()];
                        for(long u = 0; u < reihenfolge.getSize(); u++) {
                            if(reihenfolge.getBool(u))vertex[u] = circumcenter[reihenfolge.getIndex(u)];
                            else vertex[u] = otherpoints[reihenfolge.getIndex(u)];
                        }
                        surf.init(reihenfolge.getSize());
                        surf.import_vertex(vertex);
                        I   = (pointlist.getPoint("simulationGrid",neighbors[j]) -
                            pointlist.getPoint("simulationGrid",i));
                        I = I * (1/I.norm());
                        J = circumcenter[reihenfolge.getFirst()];
                        surface[i][neighbors[j]] = surf.get_surface(I);
                        testvolume += (I*J)*surf.get_surface(I);
                    }
                    else cout << "Fehler" << endl;
                }

                //abseichern der Fläche zwischen 2 Punkten
                reihenfolge.loeschen();


             }
        }
        //Nun muss noch die Fläche für den Rand dazu gerechnet werden zum Volumen
            if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],0)) {

                I.x[0] = -1;
                I.x[1] = 0;
                I.x[2] = 0;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[0];
            }
            else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[0],6e-7)){
                I.x[0] = 1;
                I.x[1] = 0;
                I.x[2] = 0;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[0];
                //cout << roh[0] << endl;
            }

            if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],0)&&neighbors.size() != 0) {
                I.x[0] = 0;
                I.x[1] = -1;
                I.x[2] = 0;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[1];

            }
            else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[1],2e-7) && neighbors.size() != 0 ){
                I.x[0] = 0;
                I.x[1] = 1;
                I.x[2] = 0;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[1];
              //  cout << roh[1] << endl;
            }

            if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],0)) {
                I.x[0] = 0;
                I.x[1] = 0;
                I.x[2] = -1;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[2];
            }
           else if(compare_koord(pointlist.getPoint("simulationGrid",i).x[2],2e-7)) {
                I.x[0] = 0;
                I.x[1] = 0;
                I.x[2] = 1;
                J = pointlist.getPoint("simulationGrid",i);
                testvolume += (I*J)*roh[2];
            }
        if(testvolume < 0) testvolume = testvolume *(-1);

        gesamtvolume += testvolume/3;
        volume[i] = testvolume/3;
        oberflaeche = oberflaeche + roh[1];
        roh[0] = 0;
        roh[1] = 0;
        roh[2] = 0;
        testvolume = 0;
        neighbors.clear();
        tetraeder.clear();
    }
    cout << gesamtvolume << endl;
}


 long VoronoiDiagram::add_newPoint(vec new_point) {
    unsigned u;
    for(u = 0; u < otherpoints.size(); u++) {
        if(compare_points(otherpoints[u], new_point)) break;
    }
    if(u == otherpoints.size()) otherpoints.push_back(new_point);
    return u;
 }



bool VoronoiDiagram::x_Ebene(long tetraedernr, long i, long neighbor) {
    vector<long> points;
    points = polyhedronlist.getPoints("simulationGrid",tetraedernr);
    points.erase(find(points.begin(), points.end(), i));
    points.erase(find(points.begin(), points.end(), neighbor));
    if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[0],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[0],6e-7)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[0],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[0],6e-7)) return true;
    else return false;
}

bool VoronoiDiagram::y_Ebene(long tetraedernr, long i, long neighbor) {
     vector<long> points;
    points = polyhedronlist.getPoints("simulationGrid",tetraedernr);
    points.erase(find(points.begin(), points.end(), i));
    points.erase(find(points.begin(), points.end(), neighbor));
    if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[1],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[1],2e-7)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[1],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[1],2e-7)) return true;
    else return false;
}
bool VoronoiDiagram::z_Ebene(long tetraedernr, long i, long neighbor) {
     vector<long> points;
    points = polyhedronlist.getPoints("simulationGrid",tetraedernr);
    points.erase(find(points.begin(), points.end(), i));
    points.erase(find(points.begin(), points.end(), neighbor));
    if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[2],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[0]).x[2],2e-7)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[2],0)) return true;
    else if(compare_koord(pointlist.getPoint("simulationGrid",points[1]).x[2],2e-7)) return true;
    else return false;
}

 double VoronoiDiagram::getvolume(long mittelpunkt) {
    return volume[mittelpunkt];
 }
double VoronoiDiagram::getSurface(long mittelpunkt, long neighbor) {
    return surface[mittelpunkt][neighbor];
}

