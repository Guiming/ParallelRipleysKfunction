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

//////////////////////// MPI //////////////////////
void Point::packHashTable(double* distances, double* weights){
	
	/*
	printf("in packHashTable: before packing ... \n");
	hashTab->printTable();
	*/
	
	int idx = 0;
	Item * head; //// TROUBLES IF PUT THESE TWO LINES IN THE LOOP BODY. NOT SURE WHY THOUGH
	LinkedList * list; ///// IT IS GOD DAMN HARD TO DEBUG AND FINALLY FIGURE, OR TRY, THIS OUT !!!!
	                   ///// MUST USE "LinkedList * list" INSTEAD OF "LinkedList list"
	
	// By Guiming @ 2015-12-20
	int n = hashTab->getLength();

	for (int i = 0; i < n; i++){
		list = &(hashTab->array)[i];
		if (list->getLength() != 0)
		{	
			//list.printList();
			head = list->getHead();
			Item * p = head;
			Item * q = head;
			while (q)
			{
				p = q;
				if (p != head){
					double d = p->distance;
					double w = p->weight;
					distances[idx] = d;
					weights[idx] = w;
					idx ++;
				}
				q = p->next;
			}
			//list.printList();
		}		
	}

	// By Guiming @ 2015-12-20
	// the last item is <keyLowerBound, 1.0>
	//distances[n] = hashTab->keyLowerBound;
	//weights[n] = 1.0;

	/*
	printf("in packHashTable: after packing ... \n");
	hashTab->printTable();
	*/
}

void Point::unpackHashTable(int n, double* distances, double* weights){
	//if (hashTab->getNumberOfItems() == 0){ // NOT !=0 
		for (int i = 0; i < n; i++){
			double d = distances[i];
			double w = weights[i];
			hashTab->insertItem(new Item(d, w, NULL));
		}
		//printf("--check: hash table restored \n");
	//}	
}
