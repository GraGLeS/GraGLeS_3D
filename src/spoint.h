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
#ifndef __SPOINT__
#define	 __SPOINT__
#include <cmath>
/*!
 * \struct SPoint
 * \brief Structure used to represent a two dimensional point
 *
 * The point represented by this structure has coordinates of type double. No operators
 * are overloaded for this structure.
 */
struct SPoint {
	SPoint() :
		x(-1), y(-1), z(-1) {
	}
	SPoint(double _x, double _y, double _z) :
		x(_x), y(_y), z(_z) {
	}
	SPoint(const SPoint& other) :
		x(other.x), y(other.y), z(other.z) {
	}
	double squaredDistanceTo(const SPoint& other) const {
		return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y)
				+ (y - other.z) * (y - other.z);
	}
	double DistanceTo(const SPoint& other) {
		double dist = sqrt(
				(x - other.x) * (x - other.x) + (y - other.y) * (y - other.y)
						+ (z - other.z) * (z - other.z));

		return dist;
	}
	bool operator==(const SPoint &other) const {
		return (x == other.x) && (y == other.y) && (z == other.z);
	}
	bool operator<(const SPoint& other) const {
		return len() < other.len();
	}
	SPoint operator+(const SPoint& other) const {
		SPoint result(0, 0, 0);
		result.x = this->x + other.x;
		result.y = this->y + other.y;
		result.z = this->z + other.z;
		return result;
	}
	SPoint operator-(const SPoint& other) const {
		SPoint result(0, 0, 0);
		result.x = this->x - other.x;
		result.y = this->y - other.y;
		result.z = this->z - other.z;
		return result;
	}
	SPoint operator*(const double other) {
		SPoint result;
		result.x = this->x * other;
		result.y = this->y * other;
		result.z = this->z * other;
		return result;
	}
	SPoint projectPointToPlane(SPoint ST, SPoint RV1, SPoint RV2) {
		SPoint n;
		SPoint result(0, 0, 0);
		double d;
		double t;
		n = RV1.cross(RV2);
		d = n.dot(ST);
		t = (d - (n.dot(*this))) / n.lenSqr();
		result = (*this) + n * t;
		return result;
	}
	double dot(const SPoint& other) const {
		return x * other.x + y * other.y + z * other.z;
	}
	SPoint cross(const SPoint& other) const {
		return SPoint(y * other.z - z * other.y, z * other.x - x * other.z,
				x * other.y - y * other.x);
	}
	double len() const {
		return sqrt(x * x + y * y + z * z);
	}
	double lenSqr() const {
		return (x * x + y * y + z * z);
	}
	double x;
	double y;
	double z;
};

#endif	//__SPOINT__
