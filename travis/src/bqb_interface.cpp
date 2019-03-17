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



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_interface.h"
#include "bqb_driver.h"
#include <stdarg.h>


const char *GetRevisionInfo_bqb_interface(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_interface() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



#define BQB_BUF_LEN 16384



CBQBInterface* BQBCreateInterface(unsigned long flags) {

	UNUSED(flags);

	CBQBInterface *p;

	p = new CBQBInterface();

	return p;
}


bool BQBDestroyInterface(CBQBInterface *i) {

	delete i;

	return true;
}


CBQBInterface::CBQBInterface() {

	m_pPrintCallback = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;

	m_iPrintLevel = BQB_PL_STANDARD;
}


CBQBInterface::~CBQBInterface() {

}


CBQBDriver* CBQBInterface::CreateDriver(unsigned long flags) {

	UNUSED(flags);

	CBQBDriver *p;

	p = new CBQBDriver(*this);

	return p;
}


bool CBQBInterface::DestroyDriver(CBQBDriver *d) {

	delete d;

	return true;
}


CBQBEngine* CBQBInterface::CreateEngine(unsigned long flags) {

	UNUSED(flags);

	CBQBEngine *p;

	p = new CBQBEngine(*this);

	return p;
}


bool CBQBInterface::DestroyEngine(CBQBEngine *e) {

	delete e;

	return true;
}


void CBQBInterface::printf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::printf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::printf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pPrintCallbackVar != NULL)
		m_pPrintCallbackVar("%s",buffer);
	else if (m_pPrintCallback != NULL)
		m_pPrintCallback(buffer);
	else
		::printf("%s",buffer);
}


void CBQBInterface::bprintf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::bprintf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::bprintf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pBPrintCallbackVar != NULL)
		m_pBPrintCallbackVar("%s",buffer);
	else if (m_pBPrintCallback != NULL)
		m_pBPrintCallback(buffer);
	else
		::printf("%s",buffer);
}


void CBQBInterface::eprintf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::eprintf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::eprintf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pEPrintCallbackVar != NULL)
		m_pEPrintCallbackVar("%s",buffer);
	else if (m_pEPrintCallback != NULL)
		m_pEPrintCallback(buffer);
	else
		::printf("%s",buffer);
}


void CBQBInterface::FlushLog() const {

}

	
void CBQBInterface::SetPrintCallback( void (*fp)(const char *) ) {

	m_pPrintCallback = fp;
	m_pEPrintCallback = fp;
	m_pBPrintCallback = fp;
}


void CBQBInterface::SetPrintCallback( void (*fp)(const char *, ...) ) {

	m_pPrintCallbackVar = fp;
	m_pEPrintCallbackVar = fp;
	m_pBPrintCallbackVar = fp;
}


void CBQBInterface::ResetPrintCallback() {

	m_pPrintCallback = NULL;
	m_pEPrintCallback = NULL;
	m_pBPrintCallback = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallbackVar = NULL;
}


void CBQBInterface::SetEPrintCallback( void (*fp)(const char *) ) {

	m_pEPrintCallback = fp;
}


void CBQBInterface::SetEPrintCallback( void (*fp)(const char *, ...) ) {

	m_pEPrintCallbackVar = fp;
}


void CBQBInterface::ResetEPrintCallback() {

	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
}


void CBQBInterface::SetBPrintCallback( void (*fp)(const char *) ) {

	m_pBPrintCallback = fp;
}


void CBQBInterface::SetBPrintCallback( void (*fp)(const char *, ...) ) {

	m_pBPrintCallbackVar = fp;
}


void CBQBInterface::ResetBPrintCallback() {

	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;
}


void CBQBInterface::SetPrintLevel(int i) {

	m_iPrintLevel = i;
}


int CBQBInterface::GetPrintLevel() const {

	return m_iPrintLevel;
}










