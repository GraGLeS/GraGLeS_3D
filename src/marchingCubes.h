#ifndef		__MARCHING_CUBES_ALGORITHM__
#define		__MARCHING_CUBES_ALGORITHM__

#include "dimensionalBufferReal.h"
#include "triangle.h"
#include <vector>
using namespace std;

class LSbox;

struct GridCell {
	Vector3d points[8];
	double values[8];
};

#define NEIGHBOR_LIST_SIZE 8

struct NeighborList {
public:
	unsigned int neighbors[NEIGHBOR_LIST_SIZE];
	NeighborList(const NeighborList& other) {
		for (int i = 0; i < NEIGHBOR_LIST_SIZE; i++) {
			neighbors[i] = other.neighbors[i];
		}
	}
	NeighborList() {
		for (int i = 0; i < NEIGHBOR_LIST_SIZE; i++) {
			neighbors[i] = 0xFFFFFFFF;
		}
	}
	void PrintNeighborList(){
		for(int i=0; i<8; i++){
			cout << neighbors[i] << "	";
		}
		cout << endl;
	}
	unsigned int getNeighborsListCount() {
		int interactingGrains = 0;
		for (int j = 0; j < NEIGHBOR_LIST_SIZE; j++)
			interactingGrains += (neighbors[j] == 0xFFFFFFFF ? 0 : 1);
		return interactingGrains;
	}
};

class MarchingCubesAlgorithm {
public:
	MarchingCubesAlgorithm(DimensionalBufferReal& distance_buffer,
			LSbox* current_grain);
	~MarchingCubesAlgorithm() {
	}
	bool generateHull(vector<Triangle>& triangles,
			vector<NeighborList>& neighborLists);
	bool isInside(int row, int column, int depth);
	const vector<unsigned int>& getIdentifiedNeighbors() {
		return m_distinctGrains;
	}
private:

	bool polygonoizeCube(GridCell& inputData, vector<Triangle>& triangles);
	Vector3d VertexInterp(Vector3d p1, Vector3d p2, double valp1, double valp2);
	int generateAdditionalInformation(vector<NeighborList>& neighborLists);
	inline void insertDistinctInteger(unsigned int intger,
			NeighborList& nieghbors, int& elem_count) const;
	void recordNewGrains(NeighborList& grainIDs);
	int insertNeighborList(NeighborList& neighborList,
			vector<NeighborList>& neighborLists);
	DimensionalBufferReal& m_DistanceBuffer;
	LSbox* m_currentGrain;

	Vector3i m_leftBottomFront;
	vector<unsigned int> m_distinctGrains;
};

#endif		//__MARCHING_CUBES_ALGORITHM__
