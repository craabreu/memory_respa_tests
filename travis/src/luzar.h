/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2019 Martin Brehm
                  2012-2019 Martin Thomas
                  2016-2019 Sascha Gehrke

    This file was written by Sascha Gehrke.

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


#ifndef LUZAR_H
#define LUZAR_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"


class CLuzarCorrelation : public CxObject
{
public:
	CLuzarCorrelation();
	~CLuzarCorrelation();

	void Init(int, int);
	void Correlate(CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*);
	int CalcSize(int);

	double GetRMS(double, double, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, int);
	void Fit(CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, int);

	int m_iInput;
	int m_iDepth;
	int m_iFFTSize;
	CFFT *m_pFFT;
	CFFT *m_pFFT2;
	CFFT *m_pFFTback;
};


#endif

