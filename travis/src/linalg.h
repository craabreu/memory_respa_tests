/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2019 Martin Brehm
                  2012-2019 Martin Thomas
                  2016-2019 Sascha Gehrke

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/


#ifndef LINALG_H
#define LINALG_H


// This must always be the first include directive
#include "config.h"

#include "xdmatrixmn.h"
#include "xdvectorn.h"


int ComputeSVD(double **a, int m, int n, double *w, double **v);

int ComputeSVD_Flat(double *a, int m, int n, double *w, double *v);

void TestSVD();

void ComputePseudoInverse(int m, double *m_in, double *m_out);

CxDMatrixMN ComputePseudoInverse(CxDMatrixMN input);

void TestPseudoInverse();

/************************/

void Solve_LeastSquares_QR(double *a, double *b, double *x, int m, int n);

void TestLeastSquaresQR();

/************************/

void SolveLinearLU(CxDMatrixMN *a, CxDVectorN *b, CxDVectorN *x);

void TestLinearLU();


#endif


