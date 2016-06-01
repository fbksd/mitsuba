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

#pragma once
#if !defined(__MITSUBA_RENDER_RECTWU_H_)
#define __MITSUBA_RENDER_RECTWU_H_

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Work unit that specifies a rectangular region in an image.
 *
 * Used for instance in \ref BlockedImageProcess
 * \ingroup librender
 */
class MTS_EXPORT_RENDER RectangularWorkUnit : public WorkUnit {
public:
	inline RectangularWorkUnit() { }

	/* WorkUnit implementation */
	void set(const WorkUnit *wu);
	void load(Stream *stream);
	void save(Stream *stream) const;

	inline const Point2i &getOffset() const { return m_offset; }
	inline const Vector2i &getSize() const { return m_size; }
    inline size_t getPipeOffset() const { return m_pipeOffset; }

	inline void setOffset(const Point2i &offset) { m_offset = offset; }
	inline void setSize(const Vector2i &size) { m_size = size; }
    inline void setPipeOffset(size_t pipeOffset) { m_pipeOffset = pipeOffset; }

	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RectangularWorkUnit() { }
private:
	Point2i m_offset;
	Vector2i m_size;
    size_t m_pipeOffset;
};



/**
 * \brief Work unit that specifies a contiguous region of memory.
 *
 * Used in the rendering server for the benchmark system.
 */
class MTS_EXPORT_RENDER ContiguousWorkUnit : public WorkUnit {
public:
    inline ContiguousWorkUnit() { }

    /* WorkUnit implementation */
    void set(const WorkUnit *wu);
    void load(Stream *stream);
    void save(Stream *stream) const;

    inline size_t getOffset() const { return m_offset; }
    inline int getNumSamples() const { return m_numSamples; }

    inline void setOffset(size_t offset) { m_offset = offset; }
    inline void setNumSamples(int numSamples) { m_numSamples = numSamples; }

    std::string toString() const;

    MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~ContiguousWorkUnit() { }
private:
    size_t m_offset;
    int m_numSamples;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_RECTWU_H_ */
