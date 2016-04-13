// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <omp.h>
#include <mpi.h>

#include "Point.h"
#include "Window.h"
#include "SpatialPointPattern.h"
#include "Utility.h"
#include "json/json.h"

using namespace std;
using namespace Json;


int main(int argc, char* argv[]) {

	/////////////// GLOBAL VARIABLES ///////////////////////////////////

	// MPI paras
	int rank, nproc, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	const int root = 0;
	MPI_Status status;
	
	// timers
	//time_t start, end;
	auto start = chrono::high_resolution_clock::now();
	auto end = chrono::high_resolution_clock::now();
	double mpi_start, mpi_end;
	double mpi_start_rev, mpi_end_rev;
	double t_inirun, t_sim;
	double t_inirun_rev, t_sim_rev;

	// OpenMP setting paras
	int max_num_threads, num_threads;

	// program paras
	double startDist;
	double maxDistance;
	int nIntervals;
	bool correction;
	bool useHash;
	bool exchangeHash;
	bool doSim;
	int nsim;
	int nsim_local;
	double errSum;
	const int num_paras = 8;
	double paras[num_paras];

	// data paras
	const string dataDir = "data/";
	const string outDir = "output/";
	Window* w;
	SpatialPointPattern* spp;
	int np, type, nw;
	double *xcoords, *ycoords, *wxcoords, *wycoords;
	double* rs, *counts;
	double *finalResult;

	// schedule. e.x. how much computation each node is gonna handle
	int *schedule; // composite schedule info for both initial run and subsequent simulations
	int job_index_range[2];
	double *ds, *ws;

	const bool check = false;
	const bool check_hash = false;
	const int check_point = 0;
	const bool bcast_rstore = false;
	
	// I/O
	const bool saveFiles = true;
	ofstream csvfs;
	/////////////// END OF GLOBAL VARIABLES ///////////////////////////////////

	/////////////// INITIALIZING MPI ///////////////////////////////////
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Get_processor_name(processor_name, &namelen);

    // starting time (need the earliest one)
    mpi_start = MPI_Wtime();
    MPI_Reduce(&mpi_start, &mpi_start_rev, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
	/////////////// END OF INITIALIZING MPI ///////////////////////////////////

    /////////////// SETTING UP OPENMP ///////////////////////////////////
	max_num_threads = omp_get_max_threads();
	num_threads = max_num_threads;
	omp_set_dynamic(0);
	omp_set_num_threads(num_threads);
	printf("\nprocess %d has %d CPU cores and %d of them will be used on %s\n", rank, max_num_threads, num_threads, processor_name);
	/////////////// END OF SETTING UP OPENMP ///////////////////////////////////

	/////////////////////// PARAMETERS AND SCHEDULE INFOR ///////////////////////////////
	// default parameters
	startDist = 0;
	maxDistance = 500000;
	nIntervals = 2000;
	correction = true;
	useHash = true;
	exchangeHash = true;
	doSim = true;
	nsim = 2;
	
	// default data
	string winjson, pntjson;
	winjson = dataDir + "us/us_sim50k.geojson";
	pntjson = dataDir + "pnt/random_5000_sim50k.geojson";
	
	// user input data and parameters
	if (argc < 10) {
		// Tell the user how to run the program
		printf("\nUsage: %s points.geojson study_area.geojson max_distance num_intervals exact_correction use_hash exchangeHash do_sim num_sims\n", argv[0]);
		return 1;
	}

	// point .geojson file
	pntjson = argv[1];
	// study area .geojson file
	winjson = argv[2];
	// maximum distance
	maxDistance = (double)atoi(argv[3]);
	// number of distance intervals
	nIntervals = atoi(argv[4]);
	// exact edge correction
	correction = atoi(argv[5]);
	// use hash table (i.e. reuse edge effect correction weights) or not
	useHash = atoi(argv[6]);
	// exchange hash among computing nodes or not
	exchangeHash = atoi(argv[7]);
	// do simulation?
	doSim = atoi(argv[8]);
	// number of simulation runs
	nsim = atoi(argv[9]);

	printf("\nprocess %d reads in data on %s...\n", rank, processor_name);

	paras[0] = startDist;
	paras[1] = maxDistance;
	paras[2] = nIntervals;
	paras[3] = correction;
	paras[4] = useHash;
	paras[5] = exchangeHash;
	paras[6] = doSim;
	paras[7] = nsim;

	if (correction){
		w = new Window(winjson, 2); // 1: rectangle, 2: polygon		
		spp = new SpatialPointPattern(pntjson, w, maxDistance);
	}
	else{
		spp = new SpatialPointPattern(pntjson, maxDistance);
		w = spp->getWindow();
	}
	
	// MPI schedule information	
	schedule = new int[nproc * 3];
	if (nproc > 1){	

		np = spp->getNumberOfPoints();	
		// initialize schedule info (for the initial run)
		int idx_accum = 0;
		// deploy less task for master node
		schedule[0] = idx_accum;			
		
		idx_accum = schedule[0] + (int)(np / (1.1 * nproc));
		
		//int quota = np * (1 - sqrt(0.5));
		//int quota = (2 * np - 1 - sqrt(2.0 * np * np - 2 * np + 1)) / 2;
		//idx_accum = schedule[0] + quota;

		schedule[1] = idx_accum - 1;
		// task granularity for remaining processes
		//int granularity = (np - idx_accum) / (nproc - 1);
		for (int i = 2; i < 2 * (nproc - 1); i += 2){
			int granularity = (np - idx_accum) / (2 * (nproc - i/2));
			schedule[i] = idx_accum;
			idx_accum = schedule[i] + granularity;
			schedule[i + 1] = idx_accum - 1;
			if(check) printf("--check: process %d from %d to %d on %s\n", i / 2, schedule[i], schedule[i + 1], processor_name);
		}
		schedule[2 * nproc - 2] = idx_accum;
		schedule[2 * nproc - 1] = np - 1;
		if(check) printf("--check: process %d from %d to %d on %s\n", nproc - 1, schedule[2 * nproc - 2], schedule[2 * nproc - 1], processor_name);

		// initialize schedule info (for simulations)
		int masterPieceSize = nsim / (1 * nproc);
		schedule[2 * nproc] = masterPieceSize;
		int pieceSize = (nsim - masterPieceSize) / (nproc - 1);
		for (int i = 1; i < nproc - 1; i++){
			schedule[2 * nproc + i] = pieceSize;
		}
		schedule[nproc * 3 - 1] = nsim - (masterPieceSize + (nproc - 2) * pieceSize);
	}
	else{
		// initialize schedule info
		schedule[0] = 0;
		schedule[1] = spp->getNumberOfPoints() - 1;
		schedule[2] = nsim;
	}
	/////////////////////// END OF PARAMETERS AND SCHEDULE INFOR ///////////////////////////////
	

	////////////////////////////// THE INITIAL RUN ////////////////////////////////////////
	printf("--process %d started computation on %s...\n", rank, processor_name);

	// first figure out how much computation to do based on schedule info
	job_index_range[0] = schedule[rank * 2];
	job_index_range[1] = schedule[rank * 2 + 1];
	printf("--process %d works on points %d to %d in the initial run on %s\n", rank, job_index_range[0], job_index_range[1], processor_name);
	
	nsim_local = schedule[2 * nproc + rank];
	if (doSim) printf("--process %d works on %d simulations on %s\n", rank, nsim_local, processor_name);

	if (maxDistance > spp->getRmax()){
		maxDistance = spp->getRmax();
	}
	double binDist = maxDistance / nIntervals;

	rs = new double[nIntervals];
	counts = new double[nIntervals];
	for (int k = 0; k < nIntervals; k++){
		rs[k] = startDist + k * binDist;
		counts[k] = 0;
	}	

	// fire the initial run and timing
	//start = time(0);	
	start = chrono::high_resolution_clock::now();
	
	spp->computeKr(nIntervals, rs, num_threads, useHash, job_index_range, counts);
	
	t_inirun = TIMEDIFF(chrono::high_resolution_clock::now(), start);
	MPI_Reduce(&t_inirun, &t_inirun_rev, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
	printf("--process %d partial computation on initial run took %.4f sec. on %s\n", rank, t_inirun, processor_name);
	
	if (check  && (check_point >= 0) && (check_point < np))
		printf("--check: hash table primitive with %d entries on %s\n", spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);

	// Broadcast the result to master node
	finalResult = new double[nIntervals];
	MPI_Reduce(counts, finalResult, nIntervals, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

	// Master node is responsible for aggregating resutls
	if (rank == root){
		int np = spp->getNumberOfPoints();
		double area = w->getArea();
		double factor = area / (PI * np * np);
		if(check) printf("--check: np, area, factor on Master: %d, %.2f, %.2f\n", np, area, factor);
		errSum = 0.0;		
		for (int i = 0; i < nIntervals; i++){
			finalResult[i] = sqrt(factor *finalResult[i]);
			errSum += abs(finalResult[i] - rs[i]);
		}
		printf("************************************************************************************\n");
		printf("Initial run result on process %d with MAE %.2f on %s\n", rank, errSum / nIntervals, processor_name);
		printf("************************************************************************************\n");

		if (saveFiles){
			csvfs.open(outDir + "initialRun.csv");
			for (int i = 0; i < nIntervals; i++){
				csvfs << rs[i] << "," << finalResult[i] << "\n";
			}
			csvfs.close();
		}
	}
	////////////////////////////// END OF THE INITIAL RUN ////////////////////////////////////////

	MPI_Barrier(MPI_COMM_WORLD);	
	////////////////////////////////// DO THE SIMULATIONS ///////////////////////////////////////
	if (doSim){
		//start = time(0);
		start = chrono::high_resolution_clock::now();
		if (check) printf("--check: started broadcasting hash tables at %.4f on %s\n", (double)time(0), processor_name);
		
		// first step: Broadcast hash tables [NOT RECOMMENDED. AKA exchangeHash = false is preferred]
		if (exchangeHash && nproc > 1 && useHash){ //// START OF BROADCASTING HASH TALBES
 	
			int np = spp->getNumberOfPoints();

			//// STAGE ONE OF EXCHANGE HASH TABLES
			for (int pnt = 0; pnt < np; pnt++){

				// if i'm NOT responsible for this point, MPI_SEND its hash table to THE responsible process
				if (pnt > job_index_range[1]){

					int n = spp->getPoint(pnt)->hashTab->getNumberOfItems(); // should be > 0

					// By Guiming @ 2015-12-20
					//n = n + 1; // + 1 for the keyLowerBound // should be > 1

					if (check && n == 0) printf("--check: broadcast empty hash table for point %d on %s\n", pnt, processor_name);
					//if (check && n == 1) printf("--check: broadcast empty hash table for point %d on %s\n", pnt, processor_name);

					ds = new double[n];
					ws = new double[n];
					for(int i = 0; i < n; i++){
						ds[i] = 0;
						ws[i] = 0;
					}

					// check hash table
					if(check_hash && pnt == check_point && rank == 1){
						printf("\n--check: point %d on process %d before pack/broadcast has %d entries on %s\n", check_point, rank,
							spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
						spp->getPoint(check_point)->hashTab->printTable();
					}

					spp->getPoint(pnt)->packHashTable(ds, ws);

					// check
					if (check && (pnt == check_point)) {
						printf("--check: hash table BEFORE exchange with %d entries on %s\n", spp->getPoint(pnt)->hashTab->getNumberOfItems(), processor_name);
						printf("--check: %d entries sent out on %s\n", n, processor_name);
						for (int x = 0; x < n; x++){
							printf("--check: %.2f %.2f on %s\n", ds[x], ws[x], processor_name);
						}
					}

					int p_respn = 0;
					for (int i = 0; i < nproc; i++){
						if (pnt >= schedule[2 * i] && pnt <= schedule[2 * i + 1]) p_respn = i;
					}
					MPI_Send(&n, 1, MPI_INT, p_respn, 0, MPI_COMM_WORLD);
					MPI_Send(ds, n, MPI_INT, p_respn, 1, MPI_COMM_WORLD);
					MPI_Send(ws, n, MPI_INT, p_respn, 2, MPI_COMM_WORLD);
					
					// clean up
					delete ds;
					delete ws;

					if (check && bcast_rstore)
						printf("--check: point %d hash table broadcasted on process %d with %d entries on %s\n", pnt, rank, n, processor_name);

					// check hash table
					if(check_hash && pnt == check_point && rank == 1){
						printf("\n--check: point %d on process %d after pack/broadcast has %d entries on %s\n", check_point, rank,
							spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
						spp->getPoint(check_point)->hashTab->printTable();
					}

					if(check) printf("--broadcast hash table for point %d pnt on process %d\n", pnt, rank);

				}

				// if i AM responsible for this point, MPI_RECEIVE other parts to form a complete hash table
				if (pnt >= job_index_range[0] && pnt <= job_index_range[1]){ //

					for(int p = 0; p < rank; p++){
						int n = 0;
						MPI_Recv(&n, 1, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
						ds = new double[n];
						ws = new double[n];
						for(int i = 0; i < n; i++){
							ds[i] = 0;
							ws[i] = 0;
						}

						MPI_Recv(ds, n, MPI_INT, p, 1, MPI_COMM_WORLD, &status);
						MPI_Recv(ws, n, MPI_INT, p, 2, MPI_COMM_WORLD, &status);

						spp->getPoint(pnt)->unpackHashTable(n, ds, ws);

						// check		
						if (check && (pnt == check_point)) {
							printf("--check: hash table AFTER exchange with %d entries on %s\n", spp->getPoint(pnt)->hashTab->getNumberOfItems(), processor_name);
							printf("--check: %d entries received on %s\n", n, processor_name);
							for (int x = 0; x < n; x++){
								printf("--check: %.2f %.2f on %s\n", ds[x], ws[x], processor_name);
							}
						}

						delete ds;
						delete ws;

						if (check && bcast_rstore)
							printf("--check: point %d hash table restored on process %d with %d entries on %s\n", pnt, rank, n, processor_name);
						
						// check hash table
						if(check_hash && pnt == check_point && rank == 0){
							printf("\n--check: point %d on process %d after receive has %d entries on %s\n", check_point, rank,
							spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
							spp->getPoint(check_point)->hashTab->printTable();
						}

						if(check) printf("--receive hash table for point %d from process %d on process %d\n", pnt, p, rank);
					}
				}
			}	
			MPI_Barrier(MPI_COMM_WORLD);

			//// STAGE TWO OF EXCHANGE HASH TABLES
			for (int pnt = 0; pnt < np; pnt++){
				// if i AM responsible for this point, MPI_BROADCAST its hash table to every other processes
				if (pnt >= job_index_range[0] && pnt <= job_index_range[1]){
	
					int n = spp->getPoint(pnt)->hashTab->getNumberOfItems(); // should be > 0

					// By Guiming @ 2015-12-20
					//n = n + 1; // + 1 for the keyLowerBound // should be > 1 

					if (check && n == 0) printf("--check: broadcast empty hash table for point %d on %s\n", pnt, processor_name);
					//if (check && n == 1) printf("--check: broadcast empty hash table for point %d on %s\n", pnt, processor_name);

					ds = new double[n];
					ws = new double[n];
					for(int i = 0; i < n; i++){
						ds[i] = 0;
						ws[i] = 0;
					}

					// check hash table
					if(check_hash && pnt == check_point && rank == 1){
						printf("\n--check: point %d on process %d before pack/broadcast has %d entries on %s\n", check_point, rank,
							spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
						spp->getPoint(check_point)->hashTab->printTable();
					}

					spp->getPoint(pnt)->packHashTable(ds, ws);

					// check
					if (check && (pnt == check_point)) {
						printf("--check: hash table BEFORE exchange with %d entries on %s\n", spp->getPoint(pnt)->hashTab->getNumberOfItems(), processor_name);
						printf("--check: %d entries sent out on %s\n", n, processor_name);
						for (int x = 0; x < n; x++){
							printf("--check: %.2f %.2f on %s\n", ds[x], ws[x], processor_name);
						}
					}

					// Broadcast to process other than itself
					MPI_Bcast(&n, 1, MPI_INT, rank, MPI_COMM_WORLD);
					MPI_Bcast(ds, n, MPI_DOUBLE, rank, MPI_COMM_WORLD);
					MPI_Bcast(ws, n, MPI_DOUBLE, rank, MPI_COMM_WORLD);
					
					// clean up
					delete ds;
					delete ws;

					if (check && bcast_rstore)
						printf("--check: point %d hash table broadcasted on process %d with %d entries on %s\n", pnt, rank, n, processor_name);


					// check hash table
					if(check_hash && pnt == check_point && rank == 1){
						printf("\n--check: point %d on process %d after pack/broadcast has %d entries on %s\n", check_point, rank,
							spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
						spp->getPoint(check_point)->hashTab->printTable();
					}
				}

				else{ // every other process MPI_BROADCAST to receive a complete hash table for pnt
					
					// which process is responsible for this point?
					int p_respn = 0;
					for (int i = 0; i < nproc; i++){
						if (pnt >= schedule[2 * i] && pnt <= schedule[2 * i + 1]) p_respn = i;
					}

					int n = 0;
					MPI_Bcast(&n, 1, MPI_INT, p_respn, MPI_COMM_WORLD);
					ds = new double[n];
					ws = new double[n];
					for(int i = 0; i < n; i++){
						ds[i] = 0;
						ws[i] = 0;
					}
					// Broadcast to process other than itself
					MPI_Bcast(ds, n, MPI_DOUBLE, p_respn, MPI_COMM_WORLD);
					MPI_Bcast(ws, n, MPI_DOUBLE, p_respn, MPI_COMM_WORLD);
					spp->getPoint(pnt)->unpackHashTable(n, ds, ws);

					// check		
					if (check && (pnt == check_point)) {
						printf("--check: hash table AFTER exchange with %d entries on %s\n", spp->getPoint(pnt)->hashTab->getNumberOfItems(), processor_name);
						printf("--check: %d entries received on %s\n", n, processor_name);
						for (int x = 0; x < n; x++){
							printf("--check: %.2f %.2f on %s\n", ds[x], ws[x], processor_name);
						}
					}

					delete ds;
					delete ws;

					if (check && bcast_rstore)
						printf("--check: point %d hash table restored on process %d with %d entries on %s\n", pnt, rank, n, processor_name);
					
					// check hash table
					if(check_hash && pnt == check_point && rank == 0){
						printf("\n--check: point %d on process %d after receive has %d entries on %s\n", check_point, rank,
						spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
						spp->getPoint(check_point)->hashTab->printTable();
					}
				}
			}	

		} //// END OF BROADCASTING HASH TALBES	
		
		// check hash table
		if(check_hash && rank == 1){
			printf("\n--check: point %d on process %d after broadcast all points has %d entries on %s\n", check_point, rank,
				spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
			spp->getPoint(check_point)->hashTab->printTable();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//end = time(0);
		end = chrono::high_resolution_clock::now();
		if (check) printf("--check: finished broadcasting hash tables at %.4f on %s\n", (double)time(0), processor_name);
		if (exchangeHash) printf("--exchanging hash tables took %.4f sec. on %s\n", TIMEDIFF(end, start), processor_name);

		// fire a test running for hashtable
		int np = spp->getNumberOfPoints();		
		int index[2] = {0, np - 1};
		if (check){
			printf("--check: number of points on process %d: %d on %s\n", rank, np, processor_name);
			printf("--check: points index in simulation: %d to %d on %s\n", index[0], index[1], processor_name);
		}

		////////////// simulation begins ////////////////////////////////
		//start = time(0);
		start = chrono::high_resolution_clock::now();
		int s_start = 0; 
		int s_end;
		if (rank == root){
			s_start = 0;
			s_end = s_start + schedule[2 * nproc + rank];
		}
		else{
			int tmp = 0;
			for (int i = 0; i < rank; i++){
				tmp = tmp + schedule[2 * nproc + i];
			}
			s_start = tmp;
			s_end = s_start + schedule[2 * nproc + rank];
		}

		printf("--process %d is responsible for simulation %d to %d on %s\n", rank, s_start, s_end - 1, processor_name);
		
		// check hashtable
		if(check_hash && rank == 1){
				printf("\n--check: point %d on process %d has %d entries before simulation on %s\n", check_point, rank,
					spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
				spp->getPoint(check_point)->hashTab->printTable();
		}

		for (int s = 0; s < nsim; s++){
			// do the simulation this process is responsible			
			//if (rank == 1 && s_start <= s && s < s_end){	// run simulaitons on one process, for debug
			if (s_start <= s && s < s_end){	
				double* counts_sim = new double[nIntervals];				
				for (int k = 0; k < nIntervals; k++){
					counts_sim[k] = 0;
				}
				// check
				if (check && check_point >= 0 && check_point < np)
					printf("--check: hash table BEFORE simulation with %d entries on %s\n", spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
				
				spp->resample();
				// start the timer
				//time_t start_inner = time(0);
				auto start_inner = chrono::high_resolution_clock::now();	
				
				// THIS PART IS CRUCIAL TO LOCATE THE BUG THAT IS IN POINT::PACKHASHTABLE()
				/* 
				int a = spp->getPoint(0)->hashTab->getNumberOfItems();
				int c = spp->getPoint(124)->hashTab->getNumberOfItems();
				printf("point 0: %d point 124: %d\n", a, c);
				if(false && rank == 0){
					printf("point 0 on process %d line 441\n", rank);	
					spp->getPoint(0)->hashTab->printTable();
				}*/

				if(check_hash && rank == 1){
					printf("--check: point %d on process %d has %d entries right before simulation on %s\n", check_point, rank,
						spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
					spp->getPoint(check_point)->hashTab->printTable();
				}		

				spp->computeKr(nIntervals, rs, num_threads, useHash, index, counts_sim);

				if(check_hash) printf("--check: *****I ran to this line on process %d right after simulation on %s********\n", rank, processor_name);

				//time_t end_inner = time(0);
				auto end_inner = chrono::high_resolution_clock::now();
				
				// check hash table
				if(check_hash && rank == 1){
					printf("\n--check: point %d on process %d has %d entries after simulation on %s\n",check_point, rank,
						spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);
					spp->getPoint(check_point)->hashTab->printTable();
				}

				// check
				if (check && check_point >= 0 && check_point < np)
					printf("--check: hash table AFTER simulation with %d entries on %s\n", spp->getPoint(check_point)->hashTab->getNumberOfItems(), processor_name);

				double area = spp->getWindow()->getArea();
				double factor = area / (PI * np * np);
				if (check){
					printf("--check: area on process %d: %.2f on %s\n", rank, area, processor_name);
					printf("--check: factor on process %d: %.2f on %s\n", rank, factor, processor_name);
				}
				errSum = 0.0;
				for (int i = 0; i < nIntervals; i++){
					counts_sim[i] = sqrt(factor *counts_sim[i]);
					errSum += abs(counts_sim[i] - rs[i]);
				}
				printf("--simulation on process %d took %.4f sec. to with MAE %.2f on %s\n", rank, TIMEDIFF(end_inner, start_inner), errSum / nIntervals, processor_name);

				//if (rank == root){ // COMMENTED OUT, SO each process does it's own IO					
					if (saveFiles){
						if(check) printf("--check: process %d writes simulation %d on %s\n", rank, s, processor_name);
						string str_s = static_cast<ostringstream*>(&(ostringstream() << s))->str();
						csvfs.open(outDir + "simulation_" + str_s + ".csv");
						for (int i = 0; i < nIntervals; i++){
							csvfs << rs[i] << "," << counts_sim[i] << "\n";
						}
						csvfs.close();
					}
					delete counts_sim;
				//}
				/*else{ // COMMENTED OUT, SO each process does it's own IO	
					MPI_Send(counts_sim, nIntervals, MPI_DOUBLE, root, s, MPI_COMM_WORLD);
					delete counts_sim;
				}*/
			}

			// Transfer or I/O
			/* // COMMENTED OUT, SO each process does it's own IO	
			if (rank == root){
				// which process is responsible for this simulation?
				int p_respn = 0;
				int sim_accu = 0;
				for (int i = 0; i < nproc; i++){
					sim_accu += schedule[2 * nproc + i];
					if (s < sim_accu){
						p_respn = i;
						break;
					}
				}
				if (p_respn != root){ 
					/// only run simulaitons on root, for debug
					double* counts_rec = new double[nIntervals];
					MPI_Recv(counts_rec, nIntervals, MPI_DOUBLE, p_respn, s, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					
					if (saveFiles){
						printf("--check: process %d writes simulation %d on %s\n", rank, s, processor_name);
						string str_s = static_cast<ostringstream*>(&(ostringstream() << s))->str();
						csvfs.open(outDir + "simulation_" + str_s + ".csv");
						for (int i = 0; i < nIntervals; i++){
							csvfs << rs[i] << "," << counts_rec[i] << "\n";
						}
						csvfs.close();
					}
					
					delete counts_rec;
				}
			}*/
		}
		//end = time(0);
		end = chrono::high_resolution_clock::now();
		//t_sim = double(end - start) / nsim_local;
		t_sim = TIMEDIFF(end, start) / nsim_local;
		MPI_Reduce(&t_sim, &t_sim_rev, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
		printf("--one simulation on process %d took %.4f sec. on average on %s\n", rank, t_sim, processor_name);
		
	}
	//////////////////////////////// END OF SIMULATIONS /////////////////////////////////

	//////////////////////////////////// CLEANING UP ////////////////////////////////////
	//MPI_Barrier(MPI_COMM_WORLD);
	delete schedule;
	delete counts;
	delete rs;
	delete finalResult;
	//delete spp; // MUST NOT have this line HERE!!!
	mpi_end= MPI_Wtime();
	MPI_Reduce(&mpi_end, &mpi_end_rev, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
	if(rank == root) {
		printf("==================================================================================\n");		
		printf("MPI program took %.4f sec. to run.\n", mpi_end_rev - mpi_start_rev);
		printf("initial run took %.4f sec. to run.\n", t_inirun_rev);
		printf("one simulation took %.4f sec. to run on average.\n", t_sim_rev);
		printf("==================================================================================\n");
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	//delete spp; // THIS CAUSE THE PROBLEM IN EXCHANGING HASH TABLES

	return 0;
}
