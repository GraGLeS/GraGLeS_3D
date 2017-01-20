/*
	GraGLeS 2D A grain growth simulation utilizing level set approaches
    Copyright (C) 2015  Christian Miessen, Nikola Velinov

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
*/
#ifndef		__DIMENSIONAL_BUFFER__
#define		__DIMENSIONAL_BUFFER__

#include <vector>
#include <iostream>
#include <cmath>
#include "ExpandingVector.h"

/*!
 * \class DimensionalBuffer
 * \brief Class that encapsulates a dimensional buffer
 *
 * Dimensional buffer for ease in management of dimensional data. The data exists in space, where
 * the minimum x and y coordinates are 0 and the maximal coordinate value is bounded by the size
 * of signed 32 bit integer.
 */
template<class T>
class	DimensionalBuffer
{
public:
	/*!
	 * \brief Basic constructor. Just initializes the object with some default values.
	 */
	DimensionalBuffer() :
			m_xMin(0), m_xMax(1), m_yMin(0), m_yMax(1), m_zMin(0), m_zMax(1)
	{}
	/*!
	* \brief Constructor, which receives the coordinates of the managed region.
	*
	* This constructor takes the boundaries of the managed region as input, and properly
	* initializes the internal memory.
    */
	DimensionalBuffer(unsigned int min_x, unsigned int min_y, unsigned int min_z,
					  unsigned int max_x, unsigned int max_y, unsigned int max_z) :
		m_xMin(min_x), m_xMax(max_x), m_yMin(min_y),
		m_yMax(max_y), m_zMin(min_z), m_zMax(max_z)
	{
		resize(m_xMin, m_yMin, m_zMin, m_xMax, m_yMax, m_zMax);
	}
	/*!
	* \brief Default destructor.
    */
	~DimensionalBuffer()
	{}
	/*!
	* \brief Method that returns the value at the given coordinates.
	*
	* This method retrieves the value at the specified coordinates. If such a value does not exist
	* this method throws an out of bound exception.
	* \param row the y coordinate of the element.
	* \param column the x coordinate of the element.
	* \param depth the z coordinate of the element.
	*/
	T& getValueAt(unsigned int row, unsigned int column, unsigned int depth)
	{
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		return m_values.at((depth-m_zMin)*(m_xMax-m_xMin)*(m_yMax-m_yMin) + (row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin));
	}
	/*!
	* \brief Method that sets the value at the given coordinates.
	*
	* This method sets the value at the specified coordinates. If such a value does not exist
	* this method throws an out of bound exception.
	* \param row the y coordinate of the element.
	* \param column the x coordinate of the element.
	* \param depth the z coordinate of the element.
	* \param value the value to be set.
	*/
	void setValueAt(unsigned int row, unsigned int column, unsigned int depth, T value)
	{
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		m_values.at((depth-m_zMin)*(m_xMax-m_xMin)*(m_yMax-m_yMin) + (row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)) = value;
	}
	/*!
	*\brief Method that checks whether a point lies within the dimensional buffer.
	*
	* \param row the y coordinate of the point to test.
	* \param column the x coordinate of the point to test.
	*/
	bool isPointInside(unsigned int row, unsigned int column, unsigned int depth) const
	{
		return row >= m_yMin && row < m_yMax &&
			   column >= m_xMin && column < m_xMax &&
			   depth >= m_zMin && depth < m_zMax;
	}

	void clearValues(){
		m_values.clear();
	}
	/*!
	* \brief This method resizes the dimensions.
	*
	* This method resizes the dimensions and properly manages the internal data.
	* \param min_x the desired new minimal x coordinate.
	* \param min_y the desired new minimal y coordinate.
	* \param min_z the desired new minimal z coordinate.
	* \param max_x the desired new maximal x coordinate.
	* \param max_y the desired new maximal y coordinate.
	* \param max_z the desired new maximal z coordinate.
	*/
	void resize(unsigned int min_x, unsigned int min_y, unsigned int min_z,
			  	unsigned int max_x, unsigned int max_y, unsigned int max_z)
	{
		m_xMin = min_x;
		m_xMax = max_x;
		m_yMin = min_y;
		m_yMax = max_y;
		m_zMin = min_z;
		m_zMax = max_z;

		m_values.resize((m_xMax - m_xMin) * (m_yMax - m_yMin) * (m_zMax - m_zMin));
	}

	/*!
	* \brief This method resizes the dimensions.
	*
	* Resizes the current area to a square area and manages the internal memory.
	* \param maximumLength The maximal value for the x or y coordinates.
	*/
	void resizeToCube(unsigned int maximumLength)
	{
		int grid_size = maximumLength-1;
		unsigned int height = m_yMax - m_yMin;
		unsigned int width  = m_xMax - m_xMin;
		unsigned int depth  = m_zMax - m_zMin;
		//First resize the rectangle to a square
		if(width == height && width == depth)
			return;

		else{
			//		find largest dimension:
			int ref;
			if (width > height)	ref = width;
			else ref = height;

			if(ref < depth) ref =depth;


			if (ref > width)
			{
				int diff = ref - width;
				m_xMin -= diff/2;
				m_xMax += diff/2 + diff%2;
			}
			if (ref > height)
			{
				int diff = ref - height;
				m_yMin -= diff/2;
				m_yMax += diff/2 + diff%2;
			}
			if (ref > depth)
			{
				int diff = ref - depth;
				m_zMin -= diff/2;
				m_zMax += diff/2 + diff%2;
			}
		}


		//Now move the cube in the bounds if it has left them
		if( m_xMin < 0 )
		{
			int delta = -(int)m_xMin;
			m_xMax += delta;
			m_xMin += delta;
		}
		else if (m_xMax > grid_size)
		{
			int delta = grid_size - m_xMax;
			m_xMax += delta;
			m_xMin += delta;
		}
		if( m_yMin < 0 )
		{
			int delta = -(int)m_yMin;
			m_yMax += delta;
			m_yMin += delta;
		}
		else if (m_yMax > grid_size)
		{
			int delta = grid_size - m_yMax;
			m_yMax += delta;
			m_yMin += delta;
		}
		if( m_zMin < 0 )
		{
			int delta = -(int)m_zMin;
			m_zMax += delta;
			m_zMin += delta;
		}
		else if (m_zMax > grid_size)
		{
			int delta = grid_size - m_zMax;
			m_zMax += delta;
			m_zMin += delta;
		}


		resize(m_xMin, m_yMin, m_zMin, m_xMax, m_yMax, m_zMax);
	}
	/*!
	* \brief This method fills the area with the provided value.
	*/
	void clearValues(T value)
	{
		std::fill(m_values.begin(), m_values.end(), value);
	}
	/*!
	* \brief This method returns the total memory allocated by this buffer.
	*/
	int getTotalMemoryUsed() const
	{
		return m_values.capacity() * sizeof(T);
	}
	/*!
	* \brief This method returns the left boundary of the managed region.
	*/
	inline int getMinX() const
	{return m_xMin;}
	/*!
	* \brief This method returns the right boundary of the managed region.
	*/
	inline int getMaxX() const
	{return m_xMax;}
	/*!
	* \brief This method returns the bottom boundary of the managed region.
	*/
	inline int getMinY() const
	{return m_yMin;}
	/*!
	* \brief This method returns the top boundary of the managed region.
	*/
	inline int getMaxY() const
	{return m_yMax;}
	/*!
	* \brief This method returns the front boundary of the managed region.
	*/
	inline int getMinZ() const
	{return m_zMin;}
	/*!
	* \brief This method returns the back boundary of the managed region.
	*/
	inline int getMaxZ() const
	{return m_zMax;}
	/*!
	* \brief This method returns a pointer to the actual data stored in the buffer.
	*/
	inline T* getRawData()
	{return &m_values[0];}

private:
	int 	m_xMin;
	int 	m_xMax;
	int 	m_yMin;
	int 	m_yMax;
	int		m_zMin;
	int		m_zMax;
protected:
	ExpandingVector<T>	m_values;
};
#endif		//__DISTANCE_BUFFER__
