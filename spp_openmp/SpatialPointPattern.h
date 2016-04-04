// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#pragma once
#include "Point.h"
#include "Window.h"
#include <vector>
#include <string>
using namespace std;

class SpatialPointPattern
{
private: Point** _pnts; // an array of points
		 vector< vector<Point*> > _vpnts;
		 int _npnts; // number of points
		 double _xmin, _xmax, _ymin, _ymax;
		 double _height; // height of the rectanglar study area
		 double _width; // width of the rectanglar study area
		 double _rmax; // max distance for calculating K(r)
		 Window* _window;
public:
	SpatialPointPattern(); // default constructor
	SpatialPointPattern(int n, Window* w); // generate n random points
	SpatialPointPattern(string json, Window* w, double maxDist); // read points from geojson file, specify window
	SpatialPointPattern(string json, double maxDist); // read points from geojson file
	~SpatialPointPattern(); // default destructor
public: 
	Point* getPoint(int i); // return ith point in the 1 dimensional array
	Point* getPoint_vec(int i); // return ith point in the 2 dimensional vector
	Point* getPoint_vec(int i, int& curCol, int& curRow); // return ith point in the 2 dimensional vector
																						 // and the index to current col and row

	int getNumberOfPoints(); // return number of points
	double computeKr(double r); // the value of Ripley's K function given distance r without edge correction
	double* computeKr(int nr, double* rs, int num_threads, bool hash = true); // the value of Ripley's K function given distances rs without edge correction
	double computeKr_vec(double r);
	double* computeKr_vec(int nr, double* rs, int num_threads, bool hash = true);
	double getWidth();
	double getHeight();
	double getRmax();
private:
	void resetCount(); // reset count to 0
	void resetCount_vec(); 
public:
	void resample(); // bootstrap (resample the same number of points with replacement)
	void resample_vec(); 
};

