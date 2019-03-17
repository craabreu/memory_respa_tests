/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2019 Martin Brehm
                  2012-2019 Martin Thomas
                  2016-2019 Sascha Gehrke

    This file was written by Martin Thomas.

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


#ifndef RAMAN_H
#define RAMAN_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"
#include "reordyn.h"


class CTimeStep;


// class CRamanDyn: public CReorDyn
// {
// public:
// 	CRamanDyn(int showMol, bool global = false);
// 	~CRamanDyn();
// 	
// 	void initialize();
// 	
// 	void getDipole0();
// 	void calcPolarizability(int fieldDirection);
// 	
// 	void finalize();
// 	
// private:
// 	CxObArray _dipole0;
// 	CxObArray _polarizability[3][3];
// 	CACF *_isoACF;
// 	CACF *_anisoACF;
// 	
// 	int _derivativeOrder;
// 	
// 	bool _global;
// };


class CRamanObservation: public CObservation
{
public:
	CRamanObservation(bool global = false);
	~CRamanObservation();
	
	void initialize();
	void getDipoleZero();
	void calcPolarizability(int fieldDirection);
	void process();
	void finalize();
	
private:
// 	CRamanDyn *_ramanDyn;
	
	char *_name;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	double _specResolution;
	bool _finiteDifferenceCorrection;
	bool _saveACF;
	bool _includeCross;
	int _quantumCorrection;
	double _quantumCorrectionTemperature;
	double _laserFreq;
	double _temperature;
	int m_fieldMode;

	
	CxDVec3Array _dipoleZero;
	CxObArray _polarizabilityCache;
};


bool gatherRaman();
bool initializeRaman();
void processRaman(CTimeStep *ts);
void finalizeRaman();

void parsePolarizability();

bool gatherPolarizabilityCalc();
bool initializePolarizabilityCalc();
void processPolarizabilityCalc(CTimeStep *ts);
void finalizePolarizabilityCalc();


#endif
