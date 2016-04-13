// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#pragma once
#include <vector>
#include "HashTable.h"

using namespace std;

class Point
{
public:
	double x; // X coord.
	double y; // Y coord.
	int count; // occurrences in the spp
public:
	HashTable* hashTab;
public:
	Point(); //default constructor
	Point(double xx, double yy); //overloaded constructor
	~Point(); //destructor
public:
	double getX(); // return X coord.
	double getY(); // return Y coord.
	double calcDist(Point* pnt); // compute distance to another point
	double calcDistX(Point* pnt); // compute distance on X axis
	double calcDistY(Point* pnt); // compute distance on Y axis
	bool withinRadius(Point* pnt, double r); // whether two points are within distance r

////// MPI //////
public:
	void packHashTable(double* distances, double* weights);
	void unpackHashTable(int n, double* distances, double* weights);
};

