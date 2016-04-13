// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

// some useful constants
#define EPSILON 0.00001
#define PI 3.141592653
#define TWOPI 6.283185307
#define QUARTERPI PI/4.0
#define IGNORABLEDIFF 500
#define COEFFICIENT 10

// yes or no
#ifndef YES
#define YES (0 == 0)
#define NO (!YES)
#endif

// timing
#define TIMEDIFF(T2,T1) chrono::duration_cast<chrono::milliseconds>(T2 - T1).count()/1000.0

// some useful operations
#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < EPSILON) ? 1 : 0)
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define BETWEEN(X,X0,X1) (((X) - (X0)) * ((X) - (X1)) <= 0)
#define UNDER(X,Y,X0,Y0,X1,Y1) (((Y1) - (Y0)) * ((X) - (X0)) >= ((Y) - (Y0)) * ((X1)- (X0)))
#define UNDERNEATH(X,Y,X0,Y0,X1,Y1) (((X0) < (X1)) ? UNDER(X,Y,X0,Y0,X1,Y1) : UNDER(X,Y,X1,Y1,X0,Y0))
#define TESTINSIDE(X,Y,X0,Y0,X1,Y1) (BETWEEN(X,X0,X1) && UNDERNEATH(X, Y, X0, Y0, X1, Y1))

// when sort the points on x, how many bins do i want?
#define F 100000
#define TRANSFORM(X) ((int)(X / F))

