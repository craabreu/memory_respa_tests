/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2019 Martin Brehm
                  2012-2019 Martin Thomas
                  2016-2019 Sascha Gehrke

    This file was written by Martin Brehm and Martin Thomas.

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


#ifndef MAINTOOLS_H
#define MAINTOOLS_H


// This must always be the first include directive
#include "config.h"

#include "travis.h"
#include "tools.h"
#include "database.h"
#include "statistics.h"
#include "element.h"
#include "base64.h"


bool OpenInputTrajectory();
bool CloseInputTrajectory();
bool InputTrajectoryEOF();
bool SeekInputTrajectory(int step);


void WriteRevisionInfo();

void CheckSourceVersion();



int bqbtool_main(int argc, const char *argv[]);

	
class CAutoCorrelation : public CxObject
{
public:
	CAutoCorrelation();
	~CAutoCorrelation();

	void Init(int input, int depth, bool fft);
	void AutoCorrelate(CxDoubleArray *inp, CxDoubleArray *outp);
	void AutoCorrelate(std::vector<double> &inp, std::vector<double> &outp);
	void AutoCorrelateSqrt(CxDoubleArray *inp, CxDoubleArray *outp);

	int m_iInput;
	int m_iDepth;
	int m_iFFTSize;
	bool m_bFFT;
	CFFT *m_pFFT;
	CFFT *m_pFFT2;
};


class CCrossCorrelation : public CxObject                                              
{                                                                                      
public:                                                                                
	CCrossCorrelation();                                                            
	~CCrossCorrelation();                                                           
                                                                                       
	void Init(int input, int depth, bool fft);                                      
	void CrossCorrelate(CxDoubleArray *inp1, CxDoubleArray *inp2, CxDoubleArray *outp);
	void CrossCorrelate(std::vector<double> &inp1, std::vector<double> &inp2, std::vector<double> &outp);
	void CrossCorrelateSymmetric(CxDoubleArray *inp1, CxDoubleArray *inp2, CxDoubleArray *outp);
                                                                                       
	int m_iInput;                                                                   
	int m_iDepth;                                                                   
	int m_iFFTSize;                                                                 
	bool m_bFFT;                                                                    
	CFFT *m_pFFT;                                                                   
	CFFT *m_pFFT2;                                                                  
	CFFT *m_pFFTback;                                                               
};    
                                                                                 

/*void AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp, int depth);
void AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp, int depth, CFFT *fft, CFFT *fft2);*/
CAnalysisGroup* AddAnalysisGroup(const char *name);
void AddAnalysis(CAnalysisGroup* g, const char *name, const char *abbrev);
void InitAnalyses();
void DumpAnalyses();
//void UniteNb();
bool ParseAtom(const char *s, int refmol, int &ty, int &rty, int &atom);
bool ParseRefSystem(int refmol, const char *s, int points);
CTimeStep* GetTimeStep(int i);
CTimeStep** GetTimeStepAddress(int i);
void CalcVelocities();
void CalcVolumetricDataTimeDev();
void CalcCurrentDensity();
void CalcForces();
/*double AtomMass(char *s);
int AtomOrd(char *s);
double AtomRadius(char *s);*/
CElement* FindElement(const char *s, bool quiet);
double GuessBoxSize();
void strtolower(char *s);
void SortAtoms();
void SortElementsLabel();
void SortElementsMass();
bool SetAnalysis(const char *s);
bool ParseFunctions(const char *s);
bool ParsePeriodic(const char *s);
void WriteHeader();
void CommandLineHelp();
bool ParseArgs(int argc, const char *argv[]);
void ParsePassiveArgs(int argc, const char *argv[]);
//void VariablesToDatabase();
//void DatabaseToVariables();
//void WriteDefaultSettings(const char *s);
void CreateDatabaseDefaults();
void LoadSettings();
void InitDatabase();
void RECURSION_BuildCDF(CObservation *o, int channel, int om, CxDoubleArray **data, double *result);
CVirtualAtom* AddVirtualAtom(int mol);
void RemoveAllElements();
void RemoveAllAtoms();
void RemoveAllAnalyses();
void RemoveAllMolecules();
void RemoveAllObservations();
void GetTravisPath();
void ReorderAtoms(int molecule);
void ReorderLikeInput();
void DoubleBoxHelper(unsigned char tpx, unsigned char tpy, unsigned char tpz);
unsigned long GraceColor(int z, double bleach);
void parseCoreCharges();
bool setupWannier();
void ParseDipole();
void parseMagneticDipole();
void DipolGrimme(const char *s);


inline CxDVector3 FoldVector(CxDVector3 v)
{
	if (g_bBoxNonOrtho)
	{
		CxDVector3 w;
		w = g_mBoxToOrtho * v;
		while (w[0] > 0.5)
			w[0] -= 1.0;
		while (w[0] <= -0.5)
			w[0] += 1.0;
		while (w[1] > 0.5)
			w[1] -= 1.0;
		while (w[1] <= -0.5)
			w[1] += 1.0;
		while (w[2] > 0.5)
			w[2] -= 1.0;
		while (w[2] <= -0.5)
			w[2] += 1.0;
		return g_mBoxFromOrtho * w;
	} else
	{
		if (g_bPeriodicX)
		{
			while (v[0] > g_fBoxX/2)
				v[0] -= g_fBoxX;
			while (v[0] <= -g_fBoxX/2)
				v[0] += g_fBoxX;
		}
		if (g_bPeriodicY)
		{
			while (v[1] > g_fBoxY/2)
				v[1] -= g_fBoxY;
			while (v[1] <= -g_fBoxY/2)
				v[1] += g_fBoxY;
		}
		if (g_bPeriodicZ)
		{
			while (v[2] > g_fBoxZ/2)
				v[2] -= g_fBoxZ;
			while (v[2] <= -g_fBoxZ/2)
				v[2] += g_fBoxZ;
		}
	}
	return v;
}


inline CxDVector3 FoldVectorPositive(CxDVector3 v)
{
	if (g_bBoxNonOrtho)
	{
		CxDVector3 w;
		w = g_mBoxToOrtho * v;
		while (w[0] >= 1.0)
			w[0] -= 1.0;
		while (w[0] < 0)
			w[0] += 1.0;
		while (w[1] >= 1.0)
			w[1] -= 1.0;
		while (w[1] < 0)
			w[1] += 1.0;
		while (w[2] >= 1.0)
			w[2] -= 1.0;
		while (w[2] < 0)
			w[2] += 1.0;
		return g_mBoxFromOrtho * w;
	} else
	{
		if (g_bPeriodicX)
		{
			while (v[0] >= g_fBoxX)
				v[0] -= g_fBoxX;
			while (v[0] < 0)
				v[0] += g_fBoxX;
		}
		if (g_bPeriodicY)
		{
			while (v[1] >= g_fBoxY)
				v[1] -= g_fBoxY;
			while (v[1] < 0)
				v[1] += g_fBoxY;
		}
		if (g_bPeriodicZ)
		{
			while (v[2] >= g_fBoxZ)
				v[2] -= g_fBoxZ;
			while (v[2] < 0)
				v[2] += g_fBoxZ;
		}
	}
	return v;
}


inline double FoldedLength(CxDVector3 v)
{
	if (g_bBoxNonOrtho)
	{
		CxDVector3 w;
		w = g_mBoxToOrtho * v;
		while (w[0] > 0.5)
			w[0] -= 1.0;
		while (w[0] <= -0.5)
			w[0] += 1.0;
		while (w[1] > 0.5)
			w[1] -= 1.0;
		while (w[1] <= -0.5)
			w[1] += 1.0;
		while (w[2] > 0.5)
			w[2] -= 1.0;
		while (w[2] <= -0.5)
			w[2] += 1.0;
		return (g_mBoxFromOrtho * w).GetLength();
	} else
	{
		if (g_bPeriodicX)
		{
			while (v[0] > g_fBoxX/2)
				v[0] -= g_fBoxX;
			while (v[0] <= -g_fBoxX/2)
				v[0] += g_fBoxX;
		}
		if (g_bPeriodicY)
		{
			while (v[1] > g_fBoxY/2)
				v[1] -= g_fBoxY;
			while (v[1] <= -g_fBoxY/2)
				v[1] += g_fBoxY;
		}
		if (g_bPeriodicZ)
		{
			while (v[2] > g_fBoxZ/2)
				v[2] -= g_fBoxZ;
			while (v[2] <= -g_fBoxZ/2)
				v[2] += g_fBoxZ;
		}
		return v.GetLength();
	}
}


/*inline CxDVector3 FoldVector1(CxDVector3 v)
{
	while (v[0] >= 0.5)
		v[0] -= 1.0;
	while (v[0] < -0.5)
		v[0] += 1.0;

	while (v[1] >= 0.5)
		v[1] -= 1.0;
	while (v[1] < -0.5)
		v[1] += 1.0;

	while (v[2] >= 0.5)
		v[2] -= 1.0;
	while (v[2] < -0.5)
		v[2] += 1.0;

	return v;
}*/


inline double MeshRandom() {

	if (g_bProcAddMeshJitter)
		return ((rand()%20001)-10000)/2000000.0*g_fProcAddMeshJitter;
	else
		return 0;
}


void BuildAtomIndices();
bool DetermineTrajFormat();
void PrintSMode();
void PrintBMode();
void WriteCredits();
void WriteCredits_Long();
unsigned long CalcFFTSize(unsigned long i, bool silent);

//void FormatTime(unsigned long eta, char *buf);
void FormatTime(unsigned long eta, CxString *buf);

void RenderStructFormulas(int tries);
void RenderFormula(const char *s, int tries);
void InitGlobalVars();

void ParseVoronoiRadii();

void DumpNonOrthoCellData();
void ExtractXYZCellGeometry(const char *s);
void ExtractXYZCellGeometry3(const char *s);



const char* GetFileExtension(const char *s);

bool IsElementMetal(const char *s);
bool IsElementNobleGas(const char *s);


#endif
