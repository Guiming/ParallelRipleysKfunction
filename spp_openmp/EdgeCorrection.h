// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#pragma once
#include "Point.h"
#include "Window.h"

class EdgeCorrection
{
private:
	Window* _w;
public:
	EdgeCorrection();
	EdgeCorrection(Window* w);
	~EdgeCorrection();
private:
	double ripleyBox(Point* pnt_i, double d_ij);
	double ripleyPoly(Point* pnt_i, double d_ij);
public:
	double ripleyCorrection(Point* pnt_i, double d_ij);
};

