/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/rectwu.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                          RectangularWorkUnit                         */
/* ==================================================================== */

void RectangularWorkUnit::set(const WorkUnit *wu) {
	const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(wu);
	m_offset = rect->m_offset;
	m_size = rect->m_size;
    m_pipeOffset = rect->m_pipeOffset;
}

void RectangularWorkUnit::load(Stream *stream) {
	int data[4];
	stream->readIntArray(data, 4);
	m_offset.x = data[0];
	m_offset.y = data[1];
	m_size.x   = data[2];
	m_size.y   = data[3];
    m_pipeOffset = stream->readULong();
}

void RectangularWorkUnit::save(Stream *stream) const {
	int data[4];
	data[0] = m_offset.x;
	data[1] = m_offset.y;
	data[2] = m_size.x;
	data[3] = m_size.y;
	stream->writeIntArray(data, 4);
    stream->writeUInt(m_pipeOffset);
}

std::string RectangularWorkUnit::toString() const {
	std::ostringstream oss;
	oss << "RectangularWorkUnit[offset=" << m_offset.toString()
		<< ", size=" << m_size.toString() << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(RectangularWorkUnit, false, WorkUnit)



/* ==================================================================== */
/*                          ContiguousWorkUnit                          */
/* ==================================================================== */

void ContiguousWorkUnit::set(const WorkUnit *wu) {
    const ContiguousWorkUnit *rect = static_cast<const ContiguousWorkUnit *>(wu);
    m_offset = rect->m_offset;
    m_numSamples = rect->m_numSamples;
}

void ContiguousWorkUnit::load(Stream *stream) {
    m_offset = stream->readULong();
    m_numSamples = stream->readInt();
}

void ContiguousWorkUnit::save(Stream *stream) const {
    stream->writeULong(m_offset);
    stream->writeInt(m_numSamples);
}

std::string ContiguousWorkUnit::toString() const {
    std::ostringstream oss;
    oss << "ContiguousWorkUnit[offset=" << m_offset
        << ", numSamples=" << m_numSamples << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS(ContiguousWorkUnit, false, WorkUnit)
MTS_NAMESPACE_END
