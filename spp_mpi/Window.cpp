// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include "Window.h"
#include "json/json.h"

using namespace std;
using namespace Json;

Window::Window(){

}

Window::Window(string json, int typee){	
	Value root;
	Reader jsonReader = Reader();

	fstream ifs(json);
	bool parsingSuccessful = jsonReader.parse(ifs, root, false);
	if (!parsingSuccessful){
		printf("failed to read file %s\n", json.c_str());
		printf("%s\n", jsonReader.getFormattedErrorMessages().c_str());
	}

	Value coordinates = root["features"][0]["geometry"]["coordinates"][0];
	_nVertex = coordinates.size() - 1;
	nSegment = _nVertex;
	type = typee;
	_vertex = new Point*[_nVertex];

	double x, y;
	xmin = coordinates[_nVertex - 1][0].asDouble();
	xmax = xmin;
	ymin = coordinates[_nVertex - 1][1].asDouble();
	ymax = ymin;

	//printf("%2.2f, %2.2f\n", x0, y0);
	for (int i = _nVertex - 1; i >= 0; i--){
		x = coordinates[i][0].asDouble();
		y = coordinates[i][1].asDouble();
		_vertex[_nVertex - 1 - i] = new Point(x, y);

		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
	}
	calcArea();
	ifs.close();
}

Window::Window(Point** vertex, int n, int typee){
	_vertex = vertex;
	_nVertex = n;
	nSegment = _nVertex;
	type = typee;

	xmin = _vertex[0]->x;
	xmax = xmin;
	ymin = _vertex[0]->y;
	ymax = ymin;

	double x, y;
	for (int i = 1; i < _nVertex; i++){
		x = _vertex[i]->x;
		y = _vertex[i]->y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
	}
	calcArea();
}


Window::~Window(){
	if (_vertex != NULL){
		delete[] _vertex;
		_vertex = NULL;
	}
}

int Window::getNumberOfVertex(){
	return _nVertex;
}

int Window::getNumberOfSegment(){
	return nSegment;
}

Point** Window::getVertex(){
	return _vertex;
}

double Window::getXmin(){
	return xmin;
}
double Window::getXmax(){
	return xmax;
}
double Window::getYmin(){
	return ymin;
}
double Window::getYmax(){
	return ymax;
}

double Window::getArea(){
	return _area;
}

int Window::getType(){
	return type;
}

void Window::calcArea(){
	if (type == 1) // calculate area for rectangle
		_area = (xmax - xmin)*(ymax - ymin);
	else{ // calculate area for simple polygon
		double area = 0.0;
		for (int i = 0; i < _nVertex - 1; ++i)
			area += _vertex[i]->x * _vertex[i + 1]->y - _vertex[i + 1]->x * _vertex[i]->y;
		area += _vertex[_nVertex - 1]->x * _vertex[0]->y - _vertex[0]->x * _vertex[_nVertex - 1]->y;
		_area = abs(area) / 2.0;
	}
}
