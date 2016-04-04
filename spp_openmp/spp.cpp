// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "Point.h"
#include "Window.h"
#include "SpatialPointPattern.h"
#include "Utility.h"
#include "json/json.h"

using namespace std;
using namespace Json;

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)

int main(int argc, char* argv[]) {

	// timer
	auto start = chrono::high_resolution_clock::now();
	ofstream csvfs, perfs; // write result to csv file

	// initialization
	int max_num_threads = omp_get_max_threads();
	int num_threads = max_num_threads; // by default use all processors
	printf("this machine has %d cores\n", max_num_threads);

	// default working directory and filenames
	string winjson, pntjson, csvFn, perfFn;
	winjson = "data/us/us_sim50k.geojson";
	pntjson = "data/pnt/random_5000_sim50k.geojson";
	csvFn = "output/random_5000_sim50k.csv";
	
	// default distance intervals
	double maxDistance = 500000;
	int nIntervals = 2000;
	double startDist = 0;

	bool useHash = false;
	bool exactCorrection = false;
	bool coeff = 1;
	int f = 1000;
	int num_sims = 1;

	bool calc1D = true;
	bool sim1D = false;
	bool calc2D = true;
	bool sim2D = false;

	// keep track of Mean Absolute Error
	double errSum = 0.0;

	// deal with user input arguments
	if (argc < 16) {
		// Tell the user how to run the program
		printf("\nUsage: %s points.geojson study_area.geojson result.csv max_distance num_intervals exact_correction use_hash coefficient f num_sims calc_1d sim_1d calc_2d sim_2d num_threads\n", argv[0]);
		return 1;
	}
	else{
		// point .geojson file
		pntjson = argv[1];

		// study area .geojson file
		winjson = argv[2];

		// output .csv file
		csvFn = argv[3];

		// maximum distance
		maxDistance = (double)atoi(argv[4]);

		// number of distance intervals
		nIntervals = atoi(argv[5]);

		// exact correction?
		exactCorrection = atoi(argv[6]);

		// use hash?
		useHash = atoi(argv[7]);

		// coefficient
		coeff = atoi(argv[8]);

		// f
		f = atoi(argv[9]);

		// number of simulations
		num_sims = atoi(argv[10]);

		// calc on 1d list?
		calc1D = atoi(argv[11]);

		// sim_1d?
		sim1D = atoi(argv[12]);

		// calc on 2d list?
		calc2D = atoi(argv[13]);

		// sim_2d?
		sim2D = atoi(argv[14]);


		// number of processors for OpenMP
		num_threads = atoi(argv[15]);
		if (num_threads > max_num_threads)
			num_threads = max_num_threads;
	}

#undef COEFFICIENT
#define COEFFICIENT coeff

	printf("\nprogram started, initializing...\n");

	SpatialPointPattern* spp;
	if (exactCorrection){		
		Window* w = new Window(winjson, 2); // 1: rectangle, 2: polygon
		spp = new SpatialPointPattern(pntjson, w, maxDistance);
		printf("..edge correction will use EXACT BOUNDARY \n");
	}
	else{		
		spp = new SpatialPointPattern(pntjson, maxDistance);
		printf("..edge correction will use BOUNDING BOX \n");
	}

	if (useHash) printf("..hash table ENABLED \n");
	else printf("..hash table NOT ENABLED \n");

	if (calc1D) printf("..will run on 1D ARRAY \n");
	else printf("..will NOT run on 1D ARRAY \n");

	if (sim1D) printf("..will run simulations on 1D ARRAY \n");
	else printf("..will NOT run simulations on 1D ARRAY \n");

	if (calc2D) printf("..will run on 2D VECTOR \n");
	else printf("..will NOT run on 2D VECTOR \n");

	if (sim2D) printf("..will run simulations on 2D VECTOR \n");
	else printf("..will NOT run simulations on 2D VECTOR \n");

	int nPtns = spp->getNumberOfPoints();

	if (maxDistance > spp->getRmax()){
		maxDistance = spp->getRmax();
	}
	double binDist = maxDistance / nIntervals;

	double* rs = new double[nIntervals];
	for (int k = 0; k < nIntervals; k++){
		rs[k] = startDist + k * binDist;
	}

	printf("initialization took %.4f sec.\n", TIMEDIFF(chrono::high_resolution_clock::now(), start));

	// set openmp
	omp_set_num_threads(num_threads);

	perfFn = csvFn;
	perfFn = perfFn + "_" + static_cast<ostringstream*>(&(ostringstream() << exactCorrection))->str().c_str();
	perfFn = perfFn + "_" + static_cast<ostringstream*>(&(ostringstream() << useHash))->str().c_str();
	perfFn = perfFn + "_" + static_cast<ostringstream*>(&(ostringstream() << num_threads))->str().c_str();
	perfFn = perfFn + ".csv";
	perfs.open(perfFn);

	perfs << "num_points, max_distance, num_intervals, exact_correction, use_hash, coefficient, f, num_sims, t_1d, mae_1d, at_sims1d, ae_sims1d, t_2d, mae_2d, at_sims2d, ae_sims2d, num_threads\n";

	double t_1d = 0;
	double mae_1d = 0;
	double at_sims1d = 0; 
	double ae_sims1d = 0; 
	double t_2d = 0;
	double mae_2d = 0;
	double at_sims2d = 0; 
	double ae_sims2d = 0;

	printf("\ncomputation began...\n");
	printf("%d points over %d distance intervals with %d sim. on %d processors\n\n", nPtns, nIntervals, num_sims, num_threads);
	
	//============ test cases on 1d array =========================
	if (calc1D){
		// parallel on 1d array
		errSum = 0.0;
		start = chrono::high_resolution_clock::now();
		double* counts = spp->computeKr(nIntervals, rs, num_threads, useHash);
		t_1d = TIMEDIFF(chrono::high_resolution_clock::now(), start);
		
		printf("parallel computation on 1d array took %.4f sec. to run with MAE ", t_1d);
		csvfs.open(csvFn + "_1d");
		for (int k = 0; k < nIntervals; k++){
			csvfs << rs[k] << "," << counts[k] << "\n";
			errSum += abs(counts[k] - rs[k]);
		}
		csvfs.close();
		delete counts;
		mae_1d = errSum / nIntervals;
		printf("%2.4f\n", mae_1d);
	}

	if (sim1D){
		// bootstrapping
		for (int i = 0; i < num_sims; i++){

			// random resamping
			spp->resample();

			errSum = 0.0;
			start = chrono::high_resolution_clock::now();
			double* counts = spp->computeKr(nIntervals, rs, num_threads, useHash);
			double tmp = TIMEDIFF(chrono::high_resolution_clock::now(), start);
			
			at_sims1d += tmp;
			string str_i = static_cast<ostringstream*>(&(ostringstream() << i))->str();
			printf("--bootstrapping %s took %.4f sec. to run with MAE ", str_i.c_str(), tmp);

			csvfs.open(csvFn + "_1d_bs" + str_i.c_str());
			for (int k = 0; k < nIntervals; k++){
				csvfs << rs[k] << "," << counts[k] << "\n";
				errSum += abs(counts[k] - rs[k]);
			}
			csvfs.close();
			delete counts;

			double tmp_e = errSum / nIntervals;
			ae_sims1d += tmp_e;

			printf("%2.4f\n", tmp_e);
		}

		if (num_sims > 0) { 
			at_sims1d = at_sims1d / num_sims;
			ae_sims1d = ae_sims1d / num_sims; 
		}
		else ae_sims1d = 0;
	}


	//============ test cases on 2d vector =========================
	// parallel on 2d vector
	if (calc2D){
		errSum = 0.0;
		start = chrono::high_resolution_clock::now();
		double* counts2 = spp->computeKr_vec(nIntervals, rs, num_threads, useHash);
		t_2d = TIMEDIFF(chrono::high_resolution_clock::now(), start);
		
		printf("\nparallel computation on 2d vector took %.4f sec. to run with MAE ", t_2d);
		csvfs.open(csvFn + "_2d");
		for (int k = 0; k < nIntervals; k++){
			csvfs << rs[k] << "," << counts2[k] << "\n";
			errSum += abs(counts2[k] - rs[k]);
		}
		csvfs.close();
		delete counts2;
		mae_2d = errSum / nIntervals;
		printf("%2.4f\n", mae_2d);
	}

	if (sim2D){
		// bootstrapping
		for (int i = 0; i < num_sims; i++){

			// random resamping
			spp->resample_vec();

			errSum = 0.0;
			start = chrono::high_resolution_clock::now();
			double* counts2 = spp->computeKr_vec(nIntervals, rs, num_threads, useHash);
			double tmp = TIMEDIFF(chrono::high_resolution_clock::now(), start);
			
			at_sims2d += tmp;
			string str_i = static_cast<ostringstream*>(&(ostringstream() << i))->str();
			printf("--bootstrapping %s took %.4f sec. to run with MAE ", str_i.c_str(), tmp);

			csvfs.open(csvFn + "_2d_bs" + str_i.c_str());
			for (int k = 0; k < nIntervals; k++){
				csvfs << rs[k] << "," << counts2[k] << "\n";
				errSum += abs(counts2[k] - rs[k]);
			}
			csvfs.close();
			delete counts2;

			double tmp_e = errSum / nIntervals;
			ae_sims2d += tmp_e;

			printf("%2.4f\n", tmp_e);
		}

		if (num_sims > 0) { 
			at_sims2d = at_sims2d / num_sims;
			ae_sims2d = ae_sims2d / num_sims;
		}
		else ae_sims2d = 0;
	}
	perfs << nPtns << ", " << maxDistance << ", " << nIntervals << ", "
		<< exactCorrection << ", " << useHash << ", " << coeff << ", "
		<< f << ", " << num_sims << ", "
		<< fixed << setprecision(4) << t_1d << ", " << mae_1d << ", " << at_sims1d << ", " << ae_sims1d << ", "
		<< t_2d << ", " << mae_2d << ", " << at_sims2d << ", " << ae_sims2d << ", "
		<< num_threads << "\n";

	perfs.close();

	delete rs;
	delete spp;
}
