// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#pragma once
/* this is the study area */
#include "Point.h"
#include <vector>
#include <cstdlib>
#include <string.h>
using namespace std;

class Window
{
private:
	Point** _vertex; // a series of vetex (points) of the study area, in counter clockwise orderd
	int _nVertex;
	//int _nSegment;
	//double _xmin, _xmax, _ymin, _ymax;
	double _area;
	//int _type;
public:
	int type;
	int nSegment;
	double xmin, xmax, ymin, ymax;
public:
	Window();
	Window(string json, int typee);
	Window(Point** vertex, int n, int typee);
	~Window();
public:
	int getNumberOfVertex();
	int getNumberOfSegment();
	Point** getVertex();
	double getXmin();
	double getXmax();
	double getYmin();
	double getYmax();
	double getArea();
	int getType();
private:
	void calcArea();
};

