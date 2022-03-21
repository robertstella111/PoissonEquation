#ifndef READ_MESH_H_INCLUDED
#define READ_MESH_H_INCLUDED

#include "vec.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class PointList2 {
private:
  vector<vector<vec>> points;
  long number;
  vector<string> name;
  vector<long> pointnumbers;

public:
  PointList2();
  void read_pointlist(string filename);
  long getnumber();
  long getlength(string search);
  vec getPoint(string ref, long i);
};

class Polyhedronlist2 {
private:
  long number;         // Anzahl der Polyhedronlisten
  vector<string> name; // der Name der einzelnen Polyhedronlisten
  vector<string> ref;  // Gibt an auf welche Punktliste sich bezogen wird
  vector<vector<vector<long>>> points; // gesamtspeicher der Punkte
  vector<long> length; // Gibt an wie lange jede einzelne Polyhedronlist ist
  PointList2 pointlist;
  vector<long> pointnumbers; // Anzahl der Punkte

public:
  void read(string filename);
  bool neighbors(long index1, long index2, long mittelpunkt, long nachbar);
  Polyhedronlist2();
  long getlength(string search);
  vector<long> getPoints(string search, long index);
  bool exist(string search, long tetraedernr, long point);
};

class VoronoiDiagram {
private:
  PointList2 pointlist;
  Polyhedronlist2 polyhedronlist;
  long number;              // Anzahl Mittelpunkte
  long maximum;             // h�chste Tetraedernr+1 die im Mesh vorkommt
  vector<vec> circumcenter; // speichert f�r alle Mittelpunkte von 0 -
                            // Tetraedernummer-1 die volumen ab
  vector<double> volume;
  vector<vector<double>> surface;
  vector<vector<vec>> normalvec;
  vector<vec> otherpoints;
  long hash_funktion_tetraeder(long t);
  long hash_funktion_points(long t);
  long add_newPoint(vec new_point);
  bool x_Ebene(long tetraedernr, long i, long neighbor);
  bool y_Ebene(long tetraedernr, long i, long neighbor);
  bool z_Ebene(long tetraedernr, long i, long neighbor);

public:
  VoronoiDiagram();
  void getVoronoiDiagram();
  double getvolume(long mittelpunkt);
  double getSurface(long mittelpunkt, long neighbor);
};

#endif // READ_MESH_H_INCLUDED
