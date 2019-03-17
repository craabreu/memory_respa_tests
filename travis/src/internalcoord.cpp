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


// This must always be the first include directive
#include "config.h"

#include "internalcoord.h"


const char *GetRevisionInfo_internalcoord(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_internalcoord() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


CMolAngle::CMolAngle()
{
	m_faData.SetName("CMolAngle::m_faData");
}


CMolAngle::~CMolAngle()
{
}


CMolAngleGroup::CMolAngleGroup()
{
	m_oaAngles.SetName("CMolAngleGroup::m_oaAngles");
}


CMolAngleGroup::~CMolAngleGroup()
{
}


CMolBond::CMolBond()
{
	m_faData.SetName("CMolBond::m_faData");
}


CMolBond::~CMolBond()
{
}


CMolBondGroup::CMolBondGroup()
{
	m_oaBonds.SetName("CMolBondGroup::m_oaBonds");
}


CMolBondGroup::~CMolBondGroup()
{
}

