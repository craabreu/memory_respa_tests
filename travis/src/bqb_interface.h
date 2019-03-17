/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2019.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



#ifndef BQB_INTERFACE_H
#define BQB_INTERFACE_H


// This must always be the first include directive
#include "bqb_config.h"



#define BQB_PL_SILENT     0  // Do not print anything on the screen
#define BQB_PL_QUIET      1
#define BQB_PL_STANDARD   2  // As in TRAVIS / bqbtool
#define BQB_PL_VERBOSE    3  // Print all information from Front-Ends
#define BQB_PL_DEBUG      4  // Also print all information from Compressor


class CBQBDriver;
class CBQBInterface;
class CBQBEngine;



CBQBInterface* BQBCreateInterface(unsigned long flags);

bool BQBDestroyInterface(CBQBInterface *i);


class CBQBInterface {
public:
	CBQBDriver* CreateDriver(unsigned long flags);

	bool DestroyDriver(CBQBDriver *d);

	CBQBEngine* CreateEngine(unsigned long flags);

	bool DestroyEngine(CBQBEngine *e);

	void SetPrintCallback( void (*fp)(const char *) );
	void SetPrintCallback( void (*fp)(const char *, ...) );

	void ResetPrintCallback();

	void SetEPrintCallback( void (*fp)(const char *) );
	void SetEPrintCallback( void (*fp)(const char *, ...) );

	void ResetEPrintCallback();

	void SetBPrintCallback( void (*fp)(const char *) );
	void SetBPrintCallback( void (*fp)(const char *, ...) );

	void ResetBPrintCallback();

	#ifdef __GNUG__  // Variadic Argument Type Checking of GCC

		// Note the implicit first "this" argument!

		void printf(const char *s, ...)  const __attribute__ ((format (printf, 2, 3)));

		void bprintf(const char *s, ...) const __attribute__ ((format (printf, 2, 3)));

		void eprintf(const char *s, ...) const __attribute__ ((format (printf, 2, 3)));

	#else

		void printf(const char *s, ...)  const;

		void bprintf(const char *s, ...) const;

		void eprintf(const char *s, ...) const;

	#endif

	void FlushLog() const;

	void SetPrintLevel(int i);

	int GetPrintLevel() const;

	bool IsPL(int i) const { return (m_iPrintLevel >= i); }


private:

	CBQBInterface();

	~CBQBInterface();

	int m_iPrintLevel;

	void (*m_pPrintCallback)(const char *);
	void (*m_pEPrintCallback)(const char *);
	void (*m_pBPrintCallback)(const char *);

	void (*m_pPrintCallbackVar)(const char *, ...);
	void (*m_pEPrintCallbackVar)(const char *, ...);
	void (*m_pBPrintCallbackVar)(const char *, ...);


	friend CBQBInterface* BQBCreateInterface(unsigned long flags);
	friend bool BQBDestroyInterface(CBQBInterface *i);
};



#endif


