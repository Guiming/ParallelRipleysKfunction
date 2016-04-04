// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <omp.h>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>

#include "SpatialPointPattern.h"
#include "EdgeCorrection.h"
#include "Utility.h"

#include "json/json.h"
#include <fstream>
using namespace std;
using namespace Json;

//extern int IGNORABLEDIFF;

bool sortY(Point* a, Point* b){
	return a->y < b->y;
}

bool sortX(vector<Point*> vp1, vector<Point*>  vp2){
	return vp1[0]->x < vp2[0]->x;
}

/* default constructor */
SpatialPointPattern::SpatialPointPattern(){
}

/* generate n random points */
SpatialPointPattern::SpatialPointPattern(int n, Window* w){
	_npnts = n;
	_window = w;
	_height = _window->getYmax() - _window->getYmin();
	_width = _window->getXmax() - _window->getXmin();
	double ripley = min(_height, _width) / 4;
	double rlarge = sqrt(1000 / (PI * _npnts / _window->getArea()));
	_rmax = min(ripley, ripley);
	_pnts = new Point*[_npnts];

	for (int i = 0; i < _npnts; i++){
		//tseed = clock();
		//srand((int)tseed % rand());
		double x0 = (double)rand() / (RAND_MAX)* 1000;
		double y0 = (double)rand() / (RAND_MAX)* 1000;
		//cout << x << " " << y << endl;

		// rotate by pi/2
		double x = x0 * cos(QUARTERPI) - y0 * sin(QUARTERPI);
		double y = y0 * cos(QUARTERPI) + x0 * sin(QUARTERPI);
		
		int g = 10;

		Point* pt = new Point(x, y);
		pt->hashTab = new HashTable(int(_rmax / g), g);
		_pnts[i] = pt;

		Point* pnt = new Point(x, y);

		pnt->hashTab = new HashTable(int(_rmax / g), g);
		//pnt->hashTab = new HashTable(100, 10);		

		bool found = false;
		for (int j = 0; j < _vpnts.size(); j++){
			int x0 = (int)(_vpnts[j][0]->x);
			if ((int)x == x0){
				_vpnts[j].push_back(pnt);
				found = true;
				break;
			}
		}
		if (!found){
			vector<Point*> pnts;
			pnts.push_back(pnt);
			_vpnts.push_back(pnts);
		}
	}

	// sort by Y
	for (int i = 0; i < _vpnts.size(); i++){
		sort(_vpnts[i].begin(), _vpnts[i].end(), 
			sortY
			//[](Point* a, Point* b) {return a->y < b->y;}
			);
	}
	// then sort by X
	sort(_vpnts.begin(), _vpnts.end(),
		sortX
		//[](vector<Point*> vp1, vector<Point*>  vp2) {return vp1[0]->x < vp2[0]->x;}
		);

	/* see whether sorting actually works */
	/*
	for (int i = 0; i < _vpnts.size(); i++){
		//printf("%d\n", _vpnts[i].size());
		for (int j = 0; j < _vpnts[i].size(); j++){
			Point* pnt = _vpnts[i][j];
			printf("%4.2f %4.2f\n", pnt->x, pnt->y);
		}
		//printf("\n");
	}
	*/
}

/* read points from geojson file, specify window*/
SpatialPointPattern::SpatialPointPattern(string json, Window* w, double maxDist){
	_window = w;
	_height = _window->getYmax() - _window->getYmin();
	_width = _window->getXmax() - _window->getXmin();
	double ripley = min(_height, _width) / 4;
	double rlarge = sqrt(1000 / (PI * _npnts / _window->getArea()));
	_rmax = min(ripley, ripley);

	Value root;
	Reader jsonReader = Reader();

	fstream ifs(json);
	bool parsingSuccessful = jsonReader.parse(ifs, root, false);
	if (!parsingSuccessful){
		printf("failed to read file %s\n", json.c_str());
		printf("%s\n", jsonReader.getFormattedErrorMessages().c_str());
	}

	Value points = root["features"];
	_npnts = points.size();
	_pnts = new Point*[_npnts];

	printf("DeltaD = %d\n", IGNORABLEDIFF);		

	for (int i = 0; i < _npnts; i++){

		double x, y;
		x = points[i]["geometry"]["coordinates"][0].asDouble();
		y = points[i]["geometry"]["coordinates"][1].asDouble();

		double factor = COEFFICIENT * IGNORABLEDIFF;
		//printf("Delta D = %ld\n", IGNORABLEDIFF);

		Point* pt = new Point(x, y);
		pt->hashTab = new HashTable(int(min(maxDist, _rmax) / factor), factor);
		_pnts[i] = pt;

		Point* pnt = new Point(x, y);		
		pnt->hashTab = new HashTable(int(min(maxDist, _rmax) / factor), factor);

		bool found = false;
		for (int j = 0; j < _vpnts.size(); j++){
			int x0 = TRANSFORM(_vpnts[j][0]->x);
			if (TRANSFORM(x) == x0){
				_vpnts[j].push_back(pnt);
				found = true;
				break;
			}
		}
		if (!found){
			vector<Point*> pnts;
			pnts.push_back(pnt);
			_vpnts.push_back(pnts);
		}
	}

	// sort by Y
	for (int i = 0; i < _vpnts.size(); i++){
		sort(_vpnts[i].begin(), _vpnts[i].end(),
			sortY
			//[](Point* a, Point* b) {return a->y < b->y;}
			);
	}
	// then sort by X
	sort(_vpnts.begin(), _vpnts.end(),
		sortX
		//[](vector<Point*> vp1, vector<Point*>  vp2) {return vp1[0]->x < vp2[0]->x;}
		);

	/* see whether sorting actually works */
	/*
	for (int i = 0; i < _vpnts.size(); i++){
		printf("%d\n", _vpnts[i].size());
		for (int j = 0; j < _vpnts[i].size(); j++){
			Point* pnt = _vpnts[i][j];
			printf("%4.2f %4.2f\n", pnt->x, pnt->y);
		}
		printf("\n");
	}
	*/
}

/* read points from geojson file, use minimum bounding box as default window*/
SpatialPointPattern::SpatialPointPattern(string json, double maxDist){

	Value root;
	Reader jsonReader = Reader();

	fstream ifs(json);
	bool parsingSuccessful = jsonReader.parse(ifs, root, false);
	if (!parsingSuccessful){
		printf("failed to read file %s\n", json.c_str());
		printf("%s\n", jsonReader.getFormattedErrorMessages().c_str());
	}

	Value points = root["features"];
	_npnts = points.size();
	_pnts = new Point*[_npnts];

	_xmin = _xmax = points[0]["geometry"]["coordinates"][0].asDouble();
	_ymin = _ymax = points[0]["geometry"]["coordinates"][1].asDouble();

	for (int i = 0; i < _npnts; i++){

		double x, y;
		x = points[i]["geometry"]["coordinates"][0].asDouble();
		y = points[i]["geometry"]["coordinates"][1].asDouble();

		if (x < _xmin)
			_xmin = x;
		if (x > _xmax)
			_xmax = x;
		if (y < _ymin)
			_ymin = y;
		if (y > _ymax)
			_ymax = y;

		double factor = COEFFICIENT * IGNORABLEDIFF;

		Point* pt = new Point(x, y);
		pt->hashTab = new HashTable(int(min(maxDist, _rmax) / factor), factor);
		_pnts[i] = pt;

		Point* pnt = new Point(x, y);
		pnt->hashTab = new HashTable(int(min(maxDist, _rmax) / factor), factor);

		bool found = false;
		for (int j = 0; j < _vpnts.size(); j++){
			int x0 = TRANSFORM(_vpnts[j][0]->x);
			if (TRANSFORM(x) == x0){
				_vpnts[j].push_back(pnt);
				found = true;
				break;
			}
		}
		if (!found){
			vector<Point*> pnts;
			pnts.push_back(pnt);
			_vpnts.push_back(pnts);
		}
	}

	// sort by Y
	for (int i = 0; i < _vpnts.size(); i++){
		sort(_vpnts[i].begin(), _vpnts[i].end(),
			sortY
			//[](Point* a, Point* b) {return a->y < b->y;}
			);
	}
	// then sort by X
	sort(_vpnts.begin(), _vpnts.end(),
		sortX
		//[](vector<Point*> vp1, vector<Point*>  vp2) {return vp1[0]->x < vp2[0]->x;}
		);

	// construct the minimum bounding box window
	Point** vertex = new Point*[4];
	vertex[0] = new Point(_xmin, _ymin);
	vertex[1] = new Point(_xmax, _ymin);
	vertex[2] = new Point(_xmax, _ymax);
	vertex[3] = new Point(_xmin, _ymax);
	_window = new Window(vertex, 4, 1); // "rectangle"

	_height = _window->getYmax() - _window->getYmin();
	_width = _window->getXmax() - _window->getXmin();
	double ripley = min(_height, _width) / 4;
	double rlarge = sqrt(1000 / (PI * _npnts / _window->getArea()));
	_rmax = min(ripley, ripley);

	/* see whether sorting actually works */
	/*
	for (int i = 0; i < _vpnts.size(); i++){
		printf("%d\n", _vpnts[i].size());
		for (int j = 0; j < _vpnts[i].size(); j++){
			Point* pnt = _vpnts[i][j];
			printf("%4.2f %4.2f\n", pnt->x, pnt->y);
		}
		printf("\n");
	}
	*/
}

/* default destructor */
SpatialPointPattern::~SpatialPointPattern(){

	// free space on points
	for (int i = 0; i < _npnts; i++){
		delete _pnts[i];
		_pnts[i] = NULL;
	}
	delete [] _pnts;
	_pnts = NULL;

	for (int i = 0; i < _vpnts.size(); i++){
		for (int j = 0; j < _vpnts[i].size(); j++){
			delete _vpnts[i][j];
			_vpnts[i][j] = NULL;
		}
		_vpnts[i].clear();
		vector<Point*>().swap(_vpnts[i]);
	}
	vector< vector<Point*> >().swap(_vpnts);

	// free space on window
	delete _window;
}

/* return ith point in the 1 dimensional array */
Point* SpatialPointPattern::getPoint(int i){
	if (i < _npnts){
		return _pnts[i];
	}
	else{
		return NULL;
	}
}

/* return ith point in the 2 dimensional vector */
Point* SpatialPointPattern::getPoint_vec(int i){
	if (i < _npnts){

		int idx = 0;
		int ncol = _vpnts.size();
		int nrow = 0;
		int col = 0;
		int row = 0;

		while(col < ncol){

			nrow = _vpnts[col].size();
			idx += nrow;

			if (i < idx) return _vpnts[col][nrow - (idx - i)];

			col++;
		}
	}
	else{
		return NULL;
	}
}

/* return ith point in the 2 dimensional vector and the index to current col and row */
Point* SpatialPointPattern::getPoint_vec(int i, int& curCol, int& curRow){
	if (i < _npnts){

		int idx = 0;
		int ncol = _vpnts.size();
		int nrow = 0;
		int col = 0;
		int row = 0;

		while(col < ncol){

			nrow = _vpnts[col].size();
			idx += nrow;

			if (i < idx) {
				curCol = col;
				curRow = nrow - (idx - i);
				return _vpnts[col][nrow - (idx - i)];
			}
			col++;
		}
	}
	else{
		return NULL;
	}
}

/*return number of points*/
int SpatialPointPattern::getNumberOfPoints(){
	return _npnts;
}

/* compute Ripley's K given distance r WITH edge correction */
double SpatialPointPattern::computeKr(double r){
	double sum = 0;
	int n = _npnts;

	EdgeCorrection* edge = new EdgeCorrection(_window);

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < n; i++){
		Point *pnt_i = getPoint(i);
		for (int j = i + 1; j < n; j++){
			Point *pnt_j = getPoint(j);
			double d = pnt_i->calcDist(pnt_j);
			if (d < r){
				double wij = edge->ripleyCorrection(pnt_i, d);
				double wji = edge->ripleyCorrection(pnt_j, d);
				double wtot = wij + wji;
				sum += wtot;
			}	
		}
	}
	return sqrt(_window->getArea() * sum / (PI * n * n));
}

/* compute Ripley's K given distances rs WITH edge correction */
double* SpatialPointPattern::computeKr(int nr, double* rs, int num_threads, bool hash){
	
	// one array for each thread
	double **counts = new double*[num_threads];
	for (int i = 0; i < num_threads; i++){
		double *count_i = new double[nr];
		for (int j = 0; j < nr; j++){
			count_i[j] = 0;
		}
		counts[i] = count_i;
	}

	// to store final result
	double* counts_total = new double[nr];
	for (int k = 0; k < nr; k++){
		counts_total[k] = 0;
	}
	
	int n = _npnts;
	EdgeCorrection* edge = new EdgeCorrection(_window);
	//Point *pnt_i, *pnt_j;

	int hits = 0; int nohits = 0;

#pragma omp parallel for schedule(dynamic) reduction(+:hits) reduction(+:nohits)
	for (int i = 0; i < n; i++){
		int thread_num = omp_get_thread_num();
		Point *pnt_i = getPoint(i);
		if (pnt_i->count > 0){
//#pragma omp parallel for schedule(dynamic)
			for (int j = i + 1; j < n; j++){
				//int thread_num = omp_get_thread_num();
				Point *pnt_j = getPoint(j);
				if (pnt_j->count > 0){
					double d = pnt_i->calcDist(pnt_j);
					if (d < rs[nr - 1]){
						double wi, wj, wtot;
						if (hash){
							// look it up in the hashTable
							Item* it_i = pnt_i->hashTab->getItemByKey(d);
							if (it_i != NULL){
								wi = it_i->weight;
								hits++;
							}
							else{ // if does not exist, insert entry
								wi = edge->ripleyCorrection(pnt_i, d);
								pnt_i->hashTab->insertItem(new Item (d, wi, NULL));
								nohits++;
							}

							Item* it_j = pnt_j->hashTab->getItemByKey(d);
							if (it_j != NULL){
								wj = it_j->weight;
								hits++;
							}
							else{ // if does not exist, insert entry
								wj = edge->ripleyCorrection(pnt_j, d);
								pnt_j->hashTab->insertItem(new Item (d, wj, NULL));
								nohits++;
							}
						}
						else{
							wi = edge->ripleyCorrection(pnt_i, d);
							wj = edge->ripleyCorrection(pnt_j, d);
						}

						wtot = pnt_i->count * pnt_j->count * (wi +  wj);
						for (int k = nr - 1; k >= 0; k--){
							if (d >= rs[k]) break;
							counts[thread_num][k] += wtot;
							//if (d < rs[k]) counts[thread_num][k] += wtot;
						}
					}
				}
			}
		}
		//pnt_i->hashTab->printTable();
	}

	//_pnts[0]->hashTab->printTable();
	//_pnts[2000]->hashTab->printTable();

	/*
	if (hash){
		printf("..average # hits: %d, ", hits / _npnts);
		printf("# not hits: %d\n", nohits / _npnts);
	}
	*/

	// aggregate to get reuslt
	double factor = _window->getArea() / (PI * n * n);
	for (int j = 0; j < nr; j++){
		for (int i = 0; i < num_threads; i++){
			counts_total[j] += counts[i][j];
		}	
		counts_total[j] = sqrt(factor * counts_total[j]);
	}

	// clean up
	for (int i = 0; i < num_threads; i++){
		delete counts[i];
		//free(counts[i]);
	}
	delete[] counts;

	return counts_total;
}

/* compute Ripley's K given distance r WITH edge correction */
double SpatialPointPattern::computeKr_vec(double r){
	double sum = 0;
	int n = _npnts;

	EdgeCorrection* edge = new EdgeCorrection(_window);

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < _vpnts.size(); i++){
		for (int j = 0; j < _vpnts[i].size(); j++){
			Point* pnt_ij = _vpnts[i][j]; // foci point
			double x_ij = pnt_ij->x;
			double y_ij = pnt_ij->y;

			// then we need to find the points within distance r with respect to this foci point 
			
			// examine the upwards rest in current column
			int size = _vpnts[i].size();
			for (int n = j + 1; n < size; n++){ // larger y
				if (_vpnts[i][n]->y - y_ij > r)
					break;
				double d = pnt_ij->calcDist(_vpnts[i][n]);
				if (d < r){
					double wij = edge->ripleyCorrection(pnt_ij, d);
					double wji = edge->ripleyCorrection(_vpnts[i][n], d);
					double wtot = wij + wji;
					sum += wtot;
				}			
			}

			// move to right to examine columns with larger x
			for (int m = i + 1; m < _vpnts.size(); m++){ 
				// columns with too large x need NO consideration
				if ((_vpnts[m][0]->x - F) - x_ij > r)
					break;

				//otherwise we DO need to look at this column
				// size of current column
				int size = _vpnts[m].size();
				double y_min = _vpnts[m][0]->y;
				double y_max = _vpnts[m][size - 1]->y;

				// look into this column, but avoid loop through the wholw column
				int start = min(j, size - 1);

				for (int n = start; n < size; n++){ // larger y
					if (_vpnts[m][n]->y - y_ij > r)
						break;
					double d = pnt_ij->calcDist(_vpnts[m][n]);
					if (d < r){
						double wij = edge->ripleyCorrection(pnt_ij, d);
						double wji = edge->ripleyCorrection(_vpnts[m][n], d);
						double wtot = wij + wji;
						sum += wtot;
					}
				}
				for (int n = start - 1; n >= 0; n--){ // smaller y
					if (y_ij - _vpnts[m][n]->y > r)
						break;
					double d = pnt_ij->calcDist(_vpnts[m][n]);
					if (d < r){
						double wij = edge->ripleyCorrection(pnt_ij, d);
						double wji = edge->ripleyCorrection(_vpnts[m][n], d);
						double wtot = wij + wji;
						sum += wtot;
					}
				}
			}

		}
	}
	return sqrt(_window->getArea() * sum / (PI * n * n));
}

/* compute Ripley's K given distances rs WITH edge correction */
double* SpatialPointPattern::computeKr_vec(int nr, double* rs, int num_threads, bool hash){

	// one array for each thread
	double **counts = new double*[num_threads];
	for (int i = 0; i < num_threads; i++){
		double *count_i = new double[nr];
		for (int j = 0; j < nr; j++){
			count_i[j] = 0;
		}
		counts[i] = count_i;
	}

	// to store final result
	double* counts_total = new double[nr];
	for (int k = 0; k < nr; k++){
		counts_total[k] = 0;
	}

	// total number of points
	int n = _npnts;
	// this is for edge correction
	EdgeCorrection* edge = new EdgeCorrection(_window);

	int hits = 0; int nohits = 0; // statistics on the hash table (if it is used)

///* in order to balance the workload in presence of skew
#pragma omp parallel for schedule (dynamic) reduction(+:hits) reduction(+:nohits)
	for(int idx = 0; idx < n; idx++){
		int thread_num = omp_get_thread_num();
		int i, j;
		Point* pnt_ij = getPoint_vec(idx, i, j); // foci point
		if (pnt_ij->count > 0){
//*/

/* unbalanced workload when data is skewed
#pragma omp parallel for schedule (dynamic) reduction(+:hits) reduction(+:nohits)
	for (int i = 0; i < _vpnts.size(); i++){
		int thread_num = omp_get_thread_num();
		for (int j = 0; j < _vpnts[i].size(); j++){
			Point* pnt_ij = _vpnts[i][j]; // foci point

			// skip those points with count = 0 (not part of the spp)
			while (pnt_ij->count == 0 && j < _vpnts[i].size() - 1){
				j++;
				pnt_ij = _vpnts[i][j];
			}
*/

			double x_ij = pnt_ij->x;
			double y_ij = pnt_ij->y;			

			// then we need to find the points within distance r with respect to this foci point 
			// examine the upwards rest in current column
			int size = _vpnts[i].size();
			for (int n = j + 1; n < size; n++){ // larger y
				//int thread_num = omp_get_thread_num();
				Point* pnt_ji = _vpnts[i][n];

				// skip those points with count = 0 (not part of the spp)
				while (pnt_ji->count == 0 && n < size - 1){
					n++;
					pnt_ji = _vpnts[i][n];
				}
				//if (pnt_ji->count > 0){

				if (pnt_ji->y - y_ij > rs[nr - 1])
					break;
				double d = pnt_ij->calcDist(pnt_ji);
				if (d < rs[nr - 1]){
					double wij, wji, wtot;
					if (hash){ // use hash table to avoid duplicated computation on weights
						// probe the hashTable
						Item* it_ij = pnt_ij->hashTab->getItemByKey(d);
						if (it_ij != NULL){
							wij = it_ij->weight;
							hits++;
						}
						else{ // if does not exist, insert entry
							wij = edge->ripleyCorrection(pnt_ij, d);
							pnt_ij->hashTab->insertItem(new Item(d, wij, NULL));
							nohits++;
						}
						// probe the hashTable
						Item* it_ji = pnt_ji->hashTab->getItemByKey(d);
						if (it_ji != NULL){
							wji = it_ji->weight;
							hits++;
						}
						else{ // if does not exist, insert entry
							wji = edge->ripleyCorrection(pnt_ji, d);
							pnt_ji->hashTab->insertItem(new Item(d, wji, NULL));
							nohits++;
						}
					}
					else{
						wij = edge->ripleyCorrection(pnt_ij, d);
						wji = edge->ripleyCorrection(pnt_ji, d);
					}

					wtot = pnt_ij->count * pnt_ji->count * (wij + wji);
					for (int k = nr - 1; k >= 0; k--){
						if (d >= rs[k]) break;
						counts[thread_num][k] += wtot;
						//if (d < rs[k]) counts[thread_num][k] += wtot;
					}
				}
			//}
			}

			int upperbd = i + 1;
			for (int m = i + 1; m < _vpnts.size(); m++){
				if ((_vpnts[m][0]->x - 2 * F) - x_ij > rs[nr - 1])
					break;
				upperbd++;
			}
			upperbd = min(upperbd, (int)_vpnts.size());

			// move to right to examine columns with larger x
			for (int m = i + 1; m < upperbd; m++){
				//int thread_num = omp_get_thread_num();
				// columns with too large x need NO consideration
				if ((_vpnts[m][0]->x - 2 * F) - x_ij > rs[nr - 1])
					break;
				//otherwise we DO need to look at this column
				// size of current column
				int size = _vpnts[m].size();
				double y_min = _vpnts[m][0]->y;
				double y_max = _vpnts[m][size - 1]->y;

				// cases where we CANNOT SKIP this column (not working properly)
				if (y_min < y_ij + rs[nr - 1] && y_max > y_ij - rs[nr - 1]){
					// look into this column, but avoid loop through the wholw column
					int start = min(j, size - 1); // simple way instead of doing the binary search

					for (int n = start; n < size; n++){ // larger y
						Point* pnt_ji = _vpnts[m][n];

						// skip those points with count = 0 (not part of the spp)
						while (pnt_ji->count == 0 && n < size - 1){
							n++;
							pnt_ji = _vpnts[m][n];
						}
						//if (pnt_ji->count > 0){

							if (pnt_ji->y - y_ij > rs[nr - 1])
								break;
							double d = pnt_ij->calcDist(pnt_ji);
							if (d < rs[nr - 1]){
								double wij, wji, wtot;
								if (hash){ // use hash table to avoid duplicated computation on weights
									// probe the hashTable
									Item* it_ij = pnt_ij->hashTab->getItemByKey(d);
									if (it_ij != NULL){
										wij = it_ij->weight;
										hits++;
									}
									else{ // if does not exist, insert entry
										wij = edge->ripleyCorrection(pnt_ij, d);
										pnt_ij->hashTab->insertItem(new Item(d, wij, NULL));
										nohits++;
									}
									// probe the hashTable
									Item* it_ji = pnt_ji->hashTab->getItemByKey(d);
									if (it_ji != NULL){
										wji = it_ji->weight;
										hits++;
									}
									else{ // if does not exist, insert entry
										wji = edge->ripleyCorrection(pnt_ji, d);
										pnt_ji->hashTab->insertItem(new Item(d, wji, NULL));
										nohits++;
									}
								}
								else{
									wij = edge->ripleyCorrection(pnt_ij, d);
									wji = edge->ripleyCorrection(pnt_ji, d);
								}
								wtot = pnt_ij->count * pnt_ji->count * (wij + wji);

								for (int k = nr - 1; k >= 0; k--){
									if (d >= rs[k]) break;
									counts[thread_num][k] += wtot;
									//if (d < rs[k]) counts[thread_num][k] += wtot;
								}
							}
						//}
					}
					for (int n = start - 1; n >= 0; n--){ // smaller y
						Point* pnt_ji = _vpnts[m][n];

						// skip those points with count = 0 (not part of the spp)
						while (pnt_ji->count == 0 && n > 0){
							n--;
							pnt_ji = _vpnts[m][n];
						}
						//if (pnt_ji->count > 0){

						if (y_ij - pnt_ji->y > rs[nr - 1])
							break;
						double d = pnt_ij->calcDist(pnt_ji);
						if (d < rs[nr - 1]){

							double wij, wji, wtot;

							if (hash){ // use hash table to avoid duplicated computation on weights
								// probe the hashTable
								Item* it_ij = pnt_ij->hashTab->getItemByKey(d);
								if (it_ij != NULL){
									wij = it_ij->weight;
									hits++;
								}
								else{ // if does not exist, insert entry
									wij = edge->ripleyCorrection(pnt_ij, d);
									pnt_ij->hashTab->insertItem(new Item(d, wij, NULL));
									nohits++;
								}
								// probe the hashTable
								Item* it_ji = pnt_ji->hashTab->getItemByKey(d);
								if (it_ji != NULL){
									wji = it_ji->weight;
									hits++;
								}
								else{ // if does not exist, insert entry
									wji = edge->ripleyCorrection(pnt_ji, d);
									pnt_ji->hashTab->insertItem(new Item(d, wji, NULL));
									nohits++;
								}
							}
							else{
								wij = edge->ripleyCorrection(pnt_ij, d);
								wji = edge->ripleyCorrection(pnt_ji, d);
							}

							wtot = pnt_ij->count * pnt_ji->count * (wij + wji);
							for (int k = nr - 1; k >= 0; k--){
								if (d >= rs[k]) break;
								counts[thread_num][k] += wtot;
								//if (d < rs[k]) counts[thread_num][k] += wtot;
							}
						}
					//}
					}
				}
			}
			//pnt_ij->hashTab->printTable();
		}
	}

	//int c = _vpnts.size();
	//int r = _vpnts[c / 2].size();
	//_vpnts[c / 2][r / 2]->hashTab->printTable();

	/*
	if (hash){
		printf("..average # hits: %d, ", hits / _npnts);
		printf("# not hits: %d\n", nohits / _npnts);
	}
	*/

	// aggregate to get reuslt
	double factor = _window->getArea()/ (PI * n * n);
	for (int j = 0; j < nr; j++){
		for (int i = 0; i < num_threads; i++){
			counts_total[j] += counts[i][j];
		}
		counts_total[j] = sqrt(factor * counts_total[j]);
	}

	// clean up
	for (int i = 0; i < num_threads; i++){
		delete counts[i];
		//free(counts[i]);
	}
	delete[] counts;

	return counts_total;
}

double SpatialPointPattern::getWidth(){ 
	return _width; 
}

double SpatialPointPattern::getHeight(){
	return _height;
}

double SpatialPointPattern::getRmax(){
	return _rmax;
}

// reset count to 0 (for resample)
void SpatialPointPattern::resetCount(){
	for (int i = 0; i < _npnts; i++){
		_pnts[i]->count = 0;
	}
}

// reset count to 0 (for resample)
void SpatialPointPattern::resetCount_vec(){
	int ncols = _vpnts.size();
	int nrows = 0;
	for (int col = 0; col < ncols; col++){
		nrows = _vpnts[col].size();
		for (int row = 0; row < nrows; row++)  {
			_vpnts[col][row]->count = 0;
		}
	}
}

// bootstrap (resample the same number of points with replacement)
void SpatialPointPattern::resample(){
	resetCount();
	int rnum = 0;
	
	srand(clock());
	int seed = rand();
	//ofstream f("D:/Dropbox/spp/output/random_numbers.csv");
	srand(seed); // seed the random number generator
	for (int i = 0; i < _npnts; i++){
		rnum = (int)(((double)rand() / RAND_MAX ) * (_npnts - 1));
		//f << rnum << "\n";
		_pnts[rnum]->count++;
	}
	//f.close();
}

void SpatialPointPattern::resample_vec(){
	resetCount_vec();
	int rnum = 0;
	srand(clock());
	int seed = rand();
	srand(seed); // seed the random number generator
	//ofstream f("D:/Dropbox/spp/output/random_numbers.csv");
	for (int i = 0; i < _npnts; i++){
		rnum = (int)(((double)rand() / RAND_MAX) * (_npnts - 1));
		//f << rnum << "\n";
		getPoint_vec(rnum)->count++;
		//getPoint_vec(i)->count = 1;
	}
	//f.close();
}
