#include "Diffusion_file.h"

#define T 8000

Diff_gen::Diff_gen(string filename) {
    polyhedronlist.read("nin_in.deva");

}

Eigen::VectorXd Diff_gen::solve() {

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


    double surface, abstand, zwischen, zwischen2;
    long totalsize;
    unsigned long j;
    PointList2 p;
    VoronoiDiagram voronoi;
    voronoi.getVoronoiDiagram();
    p.read_pointlist("nin_in.deva");
    vector<long> list;
    vector<long> tetraeder, neighbors;
    totalsize = sizeof(semic1)/sizeof(long);

    bool found;

    for(long x = 0; x < totalsize; x++) {
        if(x == 0) {
            for(long i = 0; i < 4; i++) pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
        }
        else {
             for(long i = 0; i < 4; i++) {

                for(j = 0; j < pointlist.size(); j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[x])[i] == pointlist[j]) {
                            break;
                    }
                }
                if(j == pointlist.size()) {
                        pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
                }
             }
        }
    }
    Eigen::SparseMatrix<double> mat(pointlist.size(), pointlist.size());
    Eigen::VectorXd vec(pointlist.size());
    Eigen::VectorXd ergebnis(pointlist.size());

    for(unsigned long x = 0; x < pointlist.size(); x++) {
            zwischen = 0;
            zwischen2 = 0;
            if(p.getPoint("simulationGrid",pointlist[x]).x[0] >= 5.999e-7 ) {
               //Maximum
                for(long i = 0; i < totalsize; i++) {
                    found = false;
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            found = true;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }

                    }
                }
                //Alle Punkte in neighbors und alle tetraedernummern gespeichert
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x]) - p.getPoint("simulationGrid",neighbors[i])).norm();
                     zwischen2  =1e14;

                    zwischen = -1 + surface/abstand;

                    mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;
                    //mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;

                }
                tetraeder.clear();
                neighbors.clear();
            }
            else if(p.getPoint("simulationGrid",pointlist[x]).x[0] <= 1e-9 ) {
                //Minimum
                for(long i = 0; i < totalsize; i++) {
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                           // cout << "found" << semic1[i] << endl;
                            break;
                        }
                    }
                    if(j < 4) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }

                    }
                }
                //Alle Punkte in neighbors und alle tetraedernummern gespeichert
                //cout << "Punkt: " << pointlist[x] << "  Nachbarn: " << neighbors.size() << " Tetraederanzahl:  " << tetraeder.size()<< endl;
                for(unsigned long i = 0; i < neighbors.size(); i++) {

                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x]) - p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen2  =1e14;

                    zwischen = -1 + surface/abstand;

                    mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;
                    //cout << "Index P1: " << x << " Index P2: " << index_vec(pointlist,neighbors[i]) << endl;
                }
                tetraeder.clear();
                neighbors.clear();
            }
            else {
                for(long i = 0;  i < totalsize; i++) {
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                           // cout << "found" << semic1[i] << endl;
                            break;
                        }
                    }
                    if(j < 4) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }

                    }
                }
                //Alle Punkte in neighbors und alle tetraedernummern gespeichert
                for(unsigned long i = 0; i < neighbors.size(); i++) {

                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x]) - p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen2 = 0;
                    zwischen += surface/abstand;
                    mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;
                    vec[x] = 0;
                }
                tetraeder.clear();
                neighbors.clear();
            }
             mat.insert(x,x) = -zwischen;
             vec[x] = zwischen2;
    }
    //loesen
    cout << "loese";
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    //lscg.setMaxIterations(400);
    //lscg.setTolerance(1e-23);
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(unsigned long k = 0; k < pointlist.size(); k++) cout << ergebnis[k] << endl;
    return ergebnis;
}


Eigen::VectorXd Diff_gen::solve2() {
//Sehr instabil beim Lösen!!!!!!!
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


    double surface, abstand, volume, D, surfaceintegral, dt = 1e-2;
    long totalsize;
    unsigned long j;
    PointList2 p;
    p.read_pointlist("nin_in.deva");
    vector<long> list;
    vector<long> tetraeder, neighbors;
    totalsize = sizeof(semic1)/sizeof(long);
    VoronoiDiagram voronoi;
    voronoi.getVoronoiDiagram();
    bool found;
    vector<double> v1, volumespeicher;
    vector<vector<double> > faktor;


    for(long x = 0; x < totalsize; x++) {
        if(x == 0) {
            for(long i = 0; i < 4; i++) pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
        }
        else {
          //  list = polyhedronlist.getPoints("simulationGrid",semic1[x]);
             for(long i = 0; i < 4; i++) {

                for(j = 0; j < pointlist.size(); j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[x])[i] == pointlist[j]) {
                            break;
                    }
                }
                if(j == pointlist.size()) {
                        pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
                }
             }
        }
    }
    Eigen::SparseMatrix<double> mat(pointlist.size(), pointlist.size());
    Eigen::VectorXd vec(pointlist.size());
    Eigen::VectorXd ergebnis(pointlist.size());
    //Startwert festlegen
    for(unsigned long x = 0; x < pointlist.size(); x++) {
         if(p.getPoint("simulationGrid",pointlist[x]).x[0] <= 1e-7) {
            ergebnis[x] = 1e17;
         }
         else if(p.getPoint("simulationGrid",pointlist[x]).x[0] < 5e-7) {
             ergebnis[x] = 1.1e15;
         }
         else {
             ergebnis[x] = 1e17;
         }
    }
    D = 36;
    dt = 1e-15/(T);
    for(long t = 0; t < T; t++) {

        //Anzahl wie oft LGS gelöst werden soll
        for(unsigned long x = 0; x < pointlist.size(); x++) {
            for(long i = 0; i < totalsize; i++) {
                found = false;
                for(j = 0; j < 4; j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                        found = true;
                    }
                }
                if(found == true) {
                    //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                    tetraeder.push_back(semic1[i]);
                    //Nun alle Nachbarpunkte zu diesem Punkt
                    list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                    if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[0]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[1]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[2]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[3]);
                    }
                }
            }
            surfaceintegral = 0;
                //Alle Tetraedernr mit dem Mittelpunkt in tetraeder gespeichert. Alle Nachbarpunkte in neighbors gespeichert
                for(unsigned long i = 0; i < neighbors.size(); i++) {

                    surface = voronoi.getSurface(pointlist[x],neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-p.getPoint("simulationGrid",neighbors[i])).norm();
                    surfaceintegral += (ergebnis[index_vec(pointlist,neighbors[i])]-ergebnis[x])*
                                        (surface/abstand);
                }
                volume = voronoi.getvolume(pointlist[x]);
                vec[x] = surfaceintegral*((D*dt)/volume) + ergebnis[x];
                volumespeicher.push_back(volume);
                tetraeder.clear();
                neighbors.clear();
        }
        //Alle Mittelpunkte durchgegangen
        ergebnis = vec;
        cout << t << endl;
        for(unsigned long k = 0; k < pointlist.size(); k++) cout << vec[k] << endl;
    }
    for(unsigned long k = 0; k < pointlist.size(); k++) cout << vec[k] << endl;
    return Poisson(ergebnis);
}


Eigen::VectorXd Diff_gen::solve3() {
    //Lösen mit impliciten Eulerverfahren
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

    VoronoiDiagram voronoi;
    voronoi.getVoronoiDiagram();
    double  surface, abstand, zwischen, volume, D, dt;
    long totalsize;
    unsigned long j;
    PointList2 p;
    p.read_pointlist("nin_in.deva");
    vector<long> list;
    vector<long> neighbors;
    totalsize = sizeof(semic1)/sizeof(long);
    bool found;
    vector<double> v1, volumespeicher;
    vector<vector<double> > faktor;

    //Einlesen welche Punkte in dem Gitte vorkommen
    for(long x = 0; x < totalsize; x++) {
        if(x == 0) {
            for(long i = 0; i < 4; i++) pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
        }
        else {
             for(long i = 0; i < 4; i++) {

                for(j = 0; j < pointlist.size(); j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[x])[i] == pointlist[j]) {
                            break;
                    }
                }
                if(j == pointlist.size()) {
                        pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
                }
             }
        }
    }
    //Matrix anlgen
    Eigen::SparseMatrix<double> mat(pointlist.size(), pointlist.size());
    Eigen::VectorXd vec(pointlist.size());
    Eigen::VectorXd ergebnis(pointlist.size());
    //Startwert festlegen
    for(unsigned long x = 0; x < pointlist.size(); x++) {
         if(p.getPoint("simulationGrid",pointlist[x]).x[0] <= 1e-7) {
            ergebnis[x] = 1e17;
         }
         else if(p.getPoint("simulationGrid",pointlist[x]).x[0] < 5e-7) {
             ergebnis[x] = 1e15;
         }
         else {
             ergebnis[x] = 1e17;
         }
    }

    D = 36;
    dt = 1e-15/T;
    for(long t = 0; t < T; t++) {
        if(t==0) {
        //Anzahl wie oft LGS gelöst werden soll
        //Für jeden einzelnen Punkt die Nachbaren durchgehen. und Oberflächenintegral bilden
        for(unsigned long x = 0; x < pointlist.size(); x++) {
            //Im Vector neighbors werden alle Nachbaren von dem Punkt pointlist[x] gespeichert
            for(long i = 0; i < totalsize; i++) {
                found = false;
                for(j = 0; j < 4; j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                        found = true;
                    }
                }
                if(found == true) {
                    //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                    //Nun alle Nachbarpunkte zu diesem Punkt
                    list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                    if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[0]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[1]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[2]);
                    }
                    if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                        //nicht gefuden
                        neighbors.push_back(list[3]);
                    }
                }
            }
            zwischen = 0;
                //Alle Nachbarpunkte in neighbors gespeichert
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen += surface/abstand;
                    mat.insert(x, index_vec(pointlist,neighbors[i])) = -(D*dt*surface)/abstand;
                }
                volume = voronoi.getvolume(pointlist[x]);
                mat.insert(x,x) = (D*dt)*zwischen + volume;
                vec[x] = ergebnis[x]*volume;
                neighbors.clear();
            }
        }
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
        lscg.setMaxIterations(1000);
        lscg.compute(mat);
        ergebnis = lscg.solve(vec);
        for(unsigned long x = 0; x < pointlist.size(); x++) vec[x] = ergebnis[x]*voronoi.getvolume(pointlist[x]);
    }
    for(unsigned long k = 0; k < pointlist.size(); k++) cout << ergebnis[k] << endl;
    return Poisson(ergebnis);
}

Eigen::VectorXd Diff_gen::Poisson(Eigen::VectorXd Ladunsverteilung) {
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

    double surface, abstand, zwischen, volume = 0;
    long totalsize;
    VoronoiDiagram voronoi;
    voronoi.getVoronoiDiagram();
    unsigned long j;
    PointList2 p;
    p.read_pointlist("nin_in.deva");
    vector<long> list;
    vector<long> tetraeder, neighbors;
    totalsize = sizeof(semic1)/sizeof(long);
    bool found;
    Eigen::SparseMatrix<double> mat(pointlist.size(), pointlist.size());
    Eigen::VectorXd vec(pointlist.size());
    Eigen::VectorXd ergebnis(pointlist.size());
    double u1 = 1, u2 = 0;      //Rndwerte
    for(unsigned long x = 0; x < pointlist.size(); x++) {
            zwischen = 0;
            if(p.getPoint("simulationGrid",pointlist[x]).x[0] >= 5.999e-7 ) {
               //Maximum
                for(long i = 0; i < totalsize; i++) {
                    found = false;

                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            //cout << polyhedronlist.getPoints("simulationGrid",semic1[i])[j] << "=" <<pointlist[x] << endl;
                            found = true;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }
                    }
                }
                volume = 0;
                for(unsigned long i = 0; i < neighbors.size(); i++) {

                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen = 1;
                    vec[x] = u2;
                }
                tetraeder.clear();
                neighbors.clear();
            }
            else if(p.getPoint("simulationGrid",pointlist[x]).x[0] <= 1e-14 ) {
                //Minimum
                for(long i = 0; i < totalsize; i++) {
                    found = false;
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            found = true;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }
                    }
                }
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen = +1;
                    vec[x] = u1;
                }
                tetraeder.clear();
                neighbors.clear();
            }
            else {
                for(long i = 0;  i < totalsize; i++) {
                    found = false;
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            found = true;
                            break;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }

                    }
                }
                volume = 0;
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen -= surface/abstand;
                    mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;
                }
                volume = voronoi.getvolume(pointlist[x]);
                vec[x] = (elementar*Ladunsverteilung[x]*1e6*volume)/epsilon;
                tetraeder.clear();
                neighbors.clear();
            }
            mat.insert(x,x) = zwischen;
    }
    //loesen
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    cout << "Potential" << endl;
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(unsigned long k = 0; k < pointlist.size(); k++) cout << ergebnis[k] << "V" << endl;
    return ergebnis;
}

double testfkt(double x, double y, double z) {
    double k = (2*M_PI)/(6e-7);
    return 19*pow(k,2)*cos(k*x)*cos(3*k*y)*cos(3*k*z);
}

double Randtest(double x, double y, double z) {
    double k = (2*M_PI)/(6e-7);
    return cos(k*x)*cos(3*k*y)*cos(3*k*z);
}

void Diff_gen::test1() {
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

    double surface, abstand, zwischen, volume = 0;
    long totalsize;
    VoronoiDiagram voronoi;
    voronoi.getVoronoiDiagram();
    unsigned long j;
    PointList2 p;
    p.read_pointlist("nin_in.deva");
    vector<long> list;
    vector<long> tetraeder, neighbors;
    totalsize = sizeof(semic1)/sizeof(long);
    bool found;

        //Einlesen welche Punkte in dem Gitte vorkommen
    for(long x = 0; x < totalsize; x++) {
        if(x == 0) {
            for(long i = 0; i < 4; i++) pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
        }
        else {
             for(long i = 0; i < 4; i++) {

                for(j = 0; j < pointlist.size(); j++) {
                    if(polyhedronlist.getPoints("simulationGrid",semic1[x])[i] == pointlist[j]) {
                            break;
                    }
                }
                if(j == pointlist.size()) {
                        pointlist.push_back(polyhedronlist.getPoints("simulationGrid",semic1[x])[i]);
                }
             }
        }
    }


    Eigen::SparseMatrix<double> mat(pointlist.size(), pointlist.size());
    Eigen::VectorXd vec(pointlist.size());
    Eigen::VectorXd ergebnis(pointlist.size());
    double u1 = 1, u2 = 0;      //Rndwerte
    double x_koord, y_koord, z_koord;
    for(unsigned long x = 0; x < pointlist.size(); x++) {
            zwischen = 0;
            if(p.getPoint("simulationGrid",pointlist[x]).x[0] >= 5.999e-7 ) {
               //Maximum
                for(long i = 0; i < totalsize; i++) {
                    found = false;

                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            //cout << polyhedronlist.getPoints("simulationGrid",semic1[i])[j] << "=" <<pointlist[x] << endl;
                            found = true;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }
                    }
                }
                volume = 0;
                for(unsigned long i = 0; i < neighbors.size(); i++) {

                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    x_koord = p.getPoint("simulationGrid",pointlist[x]).x[0];
                    y_koord = p.getPoint("simulationGrid",pointlist[x]).x[1];
                    z_koord = p.getPoint("simulationGrid",pointlist[x]).x[2];
                    zwischen = 1;
                    vec[x] = Randtest(x_koord,y_koord,z_koord);
                }
                tetraeder.clear();
                neighbors.clear();
            }
            else if(p.getPoint("simulationGrid",pointlist[x]).x[0] <= 1e-14 ||
                    compare_koord(p.getPoint("simulationGrid",pointlist[x]).x[1],0)||
                    compare_koord(p.getPoint("simulationGrid",pointlist[x]).x[1],2e-7)||
                    compare_koord(p.getPoint("simulationGrid",pointlist[x]).x[2],0)||
                    compare_koord(p.getPoint("simulationGrid",pointlist[x]).x[2],2e-7)) {
                //Minimum
                for(long i = 0; i < totalsize; i++) {
                    found = false;
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            found = true;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }
                    }
                }
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    x_koord = p.getPoint("simulationGrid",pointlist[x]).x[0];
                    y_koord = p.getPoint("simulationGrid",pointlist[x]).x[1];
                    z_koord = p.getPoint("simulationGrid",pointlist[x]).x[2];
                    zwischen = 1;
                    vec[x] = Randtest(x_koord,y_koord,z_koord);
                }
                tetraeder.clear();
                neighbors.clear();
            }
            else {
                for(long i = 0;  i < totalsize; i++) {
                    found = false;
                    for(j = 0; j < 4; j++) {
                        if(polyhedronlist.getPoints("simulationGrid",semic1[i])[j] == pointlist[x]) {
                            found = true;
                            break;
                        }
                    }
                    if(found == true) {
                        //in dem Dreieck semic1[i] exisitert der Punkt pointlist[x],
                        tetraeder.push_back(semic1[i]);
                        //Nun alle Nachbarpunkte zu diesem Punkt
                        list = polyhedronlist.getPoints("simulationGrid", semic1[i]);
                        if(find(neighbors.begin(), neighbors.end(),list[0]) == neighbors.end() && list[0] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[0]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[1]) == neighbors.end() && list[1] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[1]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[2]) == neighbors.end() && list[2] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[2]);
                        }
                        if(find(neighbors.begin(), neighbors.end(),list[3]) == neighbors.end() && list[3] != pointlist[x]) {
                            //nicht gefuden
                            neighbors.push_back(list[3]);
                        }

                    }
                }
                volume = 0;
                for(unsigned long i = 0; i < neighbors.size(); i++) {
                    surface = voronoi.getSurface(pointlist[x], neighbors[i]);
                    abstand = (p.getPoint("simulationGrid",pointlist[x])-
                                p.getPoint("simulationGrid",neighbors[i])).norm();
                    zwischen -= surface/abstand;
                    mat.insert(x, index_vec(pointlist,neighbors[i])) = surface/abstand;
                }
                volume = voronoi.getvolume(pointlist[x]);
                x_koord = p.getPoint("simulationGrid",pointlist[x]).x[0];
                y_koord = p.getPoint("simulationGrid",pointlist[x]).x[1];
                z_koord = p.getPoint("simulationGrid",pointlist[x]).x[2];

                vec[x] = testfkt(x_koord,y_koord,z_koord)*volume;
                cout << x << "  " << vec[x]<< endl;
                tetraeder.clear();
                neighbors.clear();
            }
            mat.insert(x,x) = zwischen;
    }
    //loesen
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > lscg;
    lscg.compute(mat);
    ergebnis = lscg.solve(vec);
    for(unsigned long k = 0; k < pointlist.size(); k++) cout << ergebnis[k] << endl;
}





