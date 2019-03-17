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


#ifndef SPECTRUM_H
#define SPECTRUM_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "tools.h"
#include "backtrace.h"


class CSpectrum : public CxObject  
{
public:
	void SetIntegral(double f);
	void SetMaxValue(double f);
	void Multiply(double f);
	void MakeDB();
	void SetMaxRWL(double f);
	void Write(const char *pre, const char *s, const char *post);
	void FromComplex(double *f);
	void Create(int i);
	double *m_pComplex;
	double *m_pData;
	double m_fMaxRWL;
	long m_iSize;
	double m_fWaveNumber;
	CSpectrum();
	virtual ~CSpectrum();

};

#endif 
