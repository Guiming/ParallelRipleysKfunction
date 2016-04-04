// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <cmath>
#include "Point.h"
#include "Utility.h"

/*default constructor*/
Point::Point(){
}

/*overloaded constructor*/
Point::Point(double xx, double yy){
	x = xx;
	y = yy;
	count = 1;
	//hashTab = new HashTable(100);
}

/*destructor*/
Point::~Point(){
	if (hashTab != NULL){
		delete hashTab;
	}
}

/*return X coord.*/
double Point::getX(){
	return x;
}

/*return Y coord.*/
double Point::getY(){
	return y;
}

/*compute distance to another point*/
double Point::calcDist(Point* pnt){
	double dx = (pnt->x - x);
	double dy = (pnt->y - y);
	return sqrt(dx * dx + dy * dy);
	//return sqrt(pow((pnt->x - _x), 2) + pow((pnt->y - _y), 2));
}

/*compute distance on X axis*/
double Point::calcDistX(Point* pnt){
	return fabs(pnt->x - x);
}

/* compute distance on Y axis */
double Point::calcDistY(Point* pnt){
	return fabs(pnt->y - y);
}

/* whether two points are within distance r */
bool Point::withinRadius(Point* pnt, double r){
	// we know for sure this is gonna be false
	if (!(calcDistY(pnt) < r) || !(calcDistX(pnt) < r)){
		return false;
	}
	// otherwise need more calcluation to determine
	double d = calcDist(pnt);
	if (d < r){
		return true;
	}
	else{
		return false;
	}
}
