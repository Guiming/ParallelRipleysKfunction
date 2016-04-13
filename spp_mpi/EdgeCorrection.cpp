// Copyright 2015 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

// Revised based on https://github.com/spatstat/spatstat/blob/master/src/corrections.c

#include <cmath>
#include <ctime>
#include "EdgeCorrection.h"
#include "Utility.h"

EdgeCorrection::EdgeCorrection(){

}

EdgeCorrection::EdgeCorrection(Window* w){
	_w = w;
}

EdgeCorrection::~EdgeCorrection(){

}

/* for edge correction */
double EdgeCorrection::ripleyCorrection(Point* pnt_i, double d_ij){
	if (_w->type == 1){
		return ripleyBox(pnt_i, d_ij);
	}
	else
		return ripleyPoly(pnt_i, d_ij);
}

/* for polygonal study area */
double EdgeCorrection::ripleyPoly(Point* pnt_i, double d_ij){
	//return 1.0;

	//time_t tt = time(0);

	int m, k, l, ncut, nchanges;
	double xcentre, ycentre, xx0, yy0, xx1, yy1, xx01, yy01;
	double x, y, radius, radius2, dx0, dx1, dy0;
	double a, b, c, t, det, sqrtdet, tmp;
	double theta[6], delta[7], tmid[7];
	double xtest, ytest, contrib, total;

	m = _w->nSegment;
	xcentre = pnt_i->x;
	ycentre = pnt_i->y;

	radius = d_ij;
	radius2 = radius * radius;
	total = 0.0;
	Point** vertex = _w->getVertex();
	for (k = 0; k < m; k++) {
		ncut = 0;
		xx0 = vertex[k]->x;
		yy0 = vertex[k]->y;
		if (k < m - 1){
			xx1 = vertex[k + 1]->x;
			yy1 = vertex[k + 1]->y;
		}
		else{
			xx1 = vertex[0]->x;
			yy1 = vertex[0]->y;
		}

		/* intersection with left edge */
		dx0 = xx0 - xcentre;
		det = radius2 - dx0 * dx0;
		if (det > 0) {
			sqrtdet = sqrt(det);
			y = ycentre + sqrtdet;
			if (y < yy0) {
				theta[ncut] = atan2(y - ycentre, dx0);
				ncut++;
			}
			y = ycentre - sqrtdet;
			if (y < yy0) {
				theta[ncut] = atan2(y - ycentre, dx0);
				ncut++;
			}
		}
		else if (det == 0) {
			if (ycentre < yy0) {
				theta[ncut] = atan2(0.0, dx0);
				ncut++;
			}
		}
		/* intersection with right edge */
		dx1 = xx1 - xcentre;
		det = radius2 - dx1 * dx1;
		if (det > 0) {
			sqrtdet = sqrt(det);
			y = ycentre + sqrtdet;
			if (y < yy1) {
				theta[ncut] = atan2(y - ycentre, dx1);
				ncut++;
			}
			y = ycentre - sqrtdet;
			if (y < yy1) {
				theta[ncut] = atan2(y - ycentre, dx1);
				ncut++;
			}
		}
		else if (det == 0) {
			if (ycentre < yy1) {
				theta[ncut] = atan2(0.0, dx1);
				ncut++;
			}
		}
		/* intersection with top segment */
		xx01 = xx1 - xx0;
		yy01 = yy1 - yy0;
		dy0 = yy0 - ycentre;
		a = xx01 * xx01 + yy01 * yy01;
		b = 2 * (xx01 * dx0 + yy01 * dy0);
		c = dx0 * dx0 + dy0 * dy0 - radius2;
		det = b * b - 4 * a * c;
		if (det > 0) {
			sqrtdet = sqrt(det);
			t = (sqrtdet - b) / (2 * a);
			if (t >= 0 && t <= 1) {
				x = xx0 + t * xx01;
				y = yy0 + t * yy01;
				theta[ncut] = atan2(y - ycentre, x - xcentre);
				++ncut;
			}
			t = (-sqrtdet - b) / (2 * a);
			if (t >= 0 && t <= 1) {
				x = xx0 + t * xx01;
				y = yy0 + t * yy01;
				theta[ncut] = atan2(y - ycentre, x - xcentre);
				++ncut;
			}
		}
		else if (det == 0) {
			t = -b / (2 * a);
			if (t >= 0 && t <= 1) {
				x = xx0 + t * xx01;
				y = yy0 + t * yy01;
				theta[ncut] = atan2(y - ycentre, x - xcentre);
				++ncut;
			}
		}
		/* for safety, force all angles to be in range [0, 2 * pi] */
		if (ncut > 0)
			for (l = 0; l < ncut; l++)
				if (theta[l] < 0)
					theta[l] += TWOPI;

		/* sort angles */
		if (ncut > 1) {
			do {
				nchanges = 0;
				for (l = 0; l < ncut - 1; l++) {
					if (theta[l] > theta[l + 1]) {
						/* swap */
						++nchanges;
						tmp = theta[l];
						theta[l] = theta[l + 1];
						theta[l + 1] = tmp;
					}
				}
			} while (nchanges > 0);
		}
		/* compute length of circumference inside polygon */
		if (ncut == 0) {
			/* entire circle is either in or out */
			xtest = xcentre + radius;
			ytest = ycentre;
			if (TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1))
				contrib = TWOPI;
			else
				contrib = 0.0;
		}
		else {
			/* find midpoints and lengths of pieces (adding theta = ) */
			delta[0] = theta[0];
			tmid[0] = theta[0] / 2;
			if (ncut > 1) {
				for (l = 1; l < ncut; l++) {
					delta[l] = theta[l] - theta[l - 1];
					tmid[l] = (theta[l] + theta[l - 1]) / 2;
				}
			}
			delta[ncut] = TWOPI - theta[ncut - 1];
			tmid[ncut] = (TWOPI + theta[ncut - 1]) / 2;
			contrib = 0.0;
			for (l = 0; l <= ncut; l++) {
				xtest = xcentre + radius * cos(tmid[l]);
				ytest = ycentre + radius * sin(tmid[l]);
				if (TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) {
					contrib += delta[l];
				}
			}
		}
		/* multiply by sign of trapezium */
		if (xx0  < xx1)
			contrib *= -1;

		total += contrib;
	}
	//printf("edge correction takes %2.3f ns.\n", (double)(time(0) - tt) * 1000000);
	return TWOPI / total;
}

/* for rectangluar study area */
double EdgeCorrection::ripleyBox(Point* pnt_i, double d_ij){
	
	int ncor;
	double xx, yy, x0, y0, x1, y1, dL, dR, dU, dD, aL, aU, aD, aR;
	double cL, cU, cD, cR, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
	double corner, extang;

	// rectangle extent
	x0 = _w->getXmin();
	y0 = _w->getYmin();
	x1 = _w->getXmax();
	y1 = _w->getYmax();

	// point i
	xx = pnt_i->x;
	yy = pnt_i->y;

	/*
	perpendicular distance from point to each edge of rectangle
	L = left, R = right, D = down, U = up
	*/
	dL = xx - x0;
	dR = x1 - xx;
	dD = yy - y0;
	dU = y1 - yy;

	/*
	test for corner of the rectangle
	*/
	ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
	corner = (ncor >= 2) ? YES : NO;

	/*
	angle between
	- perpendicular to edge of rectangle
	and
	- line from point to corner of rectangle
	*/
	bLU = atan2(dU, dL);
	bLD = atan2(dD, dL);
	bRU = atan2(dU, dR);
	bRD = atan2(dD, dR);
	bUL = atan2(dL, dU);
	bUR = atan2(dR, dU);
	bDL = atan2(dL, dD);
	bDR = atan2(dR, dD);

	/*
	half the angle subtended by the intersection between
	the circle of radius r[i,j] centred on point i
	and each edge of the rectangle (prolonged to an infinite line)
	*/
	aL = (dL < d_ij) ? acos(dL / d_ij) : 0.0;
	aR = (dR < d_ij) ? acos(dR / d_ij) : 0.0;
	aD = (dD < d_ij) ? acos(dD / d_ij) : 0.0;
	aU = (dU < d_ij) ? acos(dU / d_ij) : 0.0;

	/* apply maxima */
	cL = MIN(aL, bLU) + MIN(aL, bLD);
	cR = MIN(aR, bRU) + MIN(aR, bRD);
	cU = MIN(aU, bUL) + MIN(aU, bUR);
	cD = MIN(aD, bDL) + MIN(aD, bDR);

	/* total exterior angle over 2 pi */
	extang = (cL + cR + cU + cD) / TWOPI;

	/* add pi/2 for corners */
	if (corner)
		extang += 1 / 4;

	/* OK, now compute weight */
	return 1 / (1 - extang);
}
