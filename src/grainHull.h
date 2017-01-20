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
#ifndef __GRAIN_HULL__
#define __GRAIN_HULL__

#include "triangle.h"
#include "marchingCubes.h"
#include "Eigen/Dense"
#include <vector>
#include "InterfacialElement.h"

using namespace std;

class LSbox;
class grainhdl;
struct Face;

class GrainHull {
private:
	vector<Triangle> m_actualHull;
	vector<Vector3d> m_normalVectors;
	vector<NeighborList> m_triangleNeighborLists;
	vector<unsigned int> m_neighbors;
	vector<GrainBoundary*> m_Grainboundary;
	vector<TripleLine*> m_TripleLines;
	vector<QuadrupleJunction*> m_QuadruplePoints;
	vector<HighOrderJunction*> m_HighOrderJunctions;
	double m_LD;
	double m_TripleLineLength;
public:
	friend class GrainBoundary;
	friend class TripleLine;
	friend class QuadrupleJunction;
	friend class HighOrderJunction;
	LSbox* m_owner;
	GrainHull(LSbox* owner);
	~GrainHull();
	bool generateHull();
	const NeighborList& getNeighborList(const Triangle& triangle);
	double computeGrainVolume();
	double computeSurfaceArea();
	const Triangle& projectPointToSurface(Vector3d& point);
	const vector<unsigned int>& getAllNeighbors();
	inline unsigned int getAllNeighborsCount() {
		return m_neighbors.size();
	}
	void plotNeighboringGrains(bool absoluteCoordinates,int timestep,int order);
	void plotContour(bool absoluteCoordinates, int timestep);
	HighOrderJunction* findHighOrderJunction(int key);
	QuadrupleJunction* findQuadrupleJunction(int key);
	TripleLine* findTripleLine(int key);
	GrainBoundary* findGrainBoundary(int key);
	//new:
	void clearInterfacialElements();
	void computeJunctionPosition();
	void correctJunctionPositionWithNeighborInformation();
	void computeGrainBoundaryElements();
	void subDivideTrianglesToInterfacialElements();
	Vector3d findClosestJunctionTo(Vector3d myposition);
	void mergeJunction();
	void switchBufferPositions();
	void computeInterfacialElementMesh();
	void meanWidth();
	void computeTriplelineLength();
	GBInfo projectPointToGrainBoundary(Vector3d& point, int id);
	void plotInterfacialElements(bool absoluteCoordinates, int timestep);
	bool IsNeighbor(int grain) {
		for (auto it : m_neighbors)
			if (it == grain)
				return true;
		return false;
	}

	vector<Face>* get_Faces();

	inline int get_numberOfQuadruplePoints() {
		return m_QuadruplePoints.size();
	}
	inline int get_numberOfHighOrderJunctions() {
		return m_HighOrderJunctions.size();
	}

	inline double getMeanWidth() {
		return m_LD;
	}
	inline double getTripleLineLength() {
		return m_TripleLineLength;
	}
	inline double getNumberOfTripleLines() {
		return m_TripleLines.size();
	}
	inline double getNumberOfQuadruplePoints() {
		return m_QuadruplePoints.size();
	}
	inline double getNumberOfHighOrderJunctions() {
		return m_HighOrderJunctions.size();
	}
	inline vector<Triangle> get_actualHull(){
		return m_actualHull;
	}
};

#endif //__GRAIN_HULL__
