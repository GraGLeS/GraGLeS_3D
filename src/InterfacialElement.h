/*
 * InterfacialElement.h
 *
 *  Created on: Nov 24, 2015
 *      Author: cm654063
 */

#ifndef INTERFACIALELEMENTS_H
#define INTERFACIALELEMENTS_H

#include <vector>
#include "triangle.h"
struct NeighborList;

class LSbox;
class grainhdl;
class Settings;
class MarchingCubesAlgorithm;
class myQuaternion;
class GrainHull;

using namespace std;

struct GBInfo {
	double energy;
	double mobility;
	GBInfo() {
	}
	GBInfo(double _mobility, double _energy) :
			energy(_energy), mobility(_mobility) {
	}
	GBInfo(const GBInfo& other) :
			mobility(other.mobility), energy(other.energy) {
	}
};

class InterfacialElement {
protected:
	vector<int> m_neighborIDs;
	vector<Triangle> m_Triangles;
	vector<Vector3d> m_barycenterTriangles;
	vector<Vector3d> m_UnitNormalTriangles;
	int m_Key_NeighborList;
	double m_mobility;
	double m_energy;
	GrainHull* m_owner;
public:
	InterfacialElement();
	InterfacialElement(int key, GrainHull* owner);
	virtual ~InterfacialElement();
	double computeMobilityMisori(double misori);
	double computeReadShockleyEnergy(double misori);

	virtual void computeEnergy()= 0;
	virtual void computeMobility()= 0;
	virtual Vector3d get_Position()=0;
	virtual void set_Position(Vector3d newposition)=0;
	virtual void set_BufferPosition(Vector3d newposition)=0;
	virtual void switch_BufferPosition()=0;

	void addBaryCenter(Vector3d current) {
		m_barycenterTriangles.push_back(current);
	}
	void addTriangle(Triangle current) {
		m_Triangles.push_back(current);
	}
	inline int get_m_Key_NeighborList() {
		return m_Key_NeighborList;
	}
	inline GBInfo get_GBInfo() {
		return GBInfo(m_energy, m_mobility);
	}
	inline double get_energy() {
		return m_energy;
	}
	inline double get_mobility() {
		return m_mobility;
	}
	inline vector<Triangle>* get_Triangles() {
		return &m_Triangles;
	}
	inline Triangle get_Triangles(int i) {
		return m_Triangles[i];
	}
	inline vector<Vector3d>* get_barycenterTriangles() {
		return &m_barycenterTriangles;
	}
	inline vector<Vector3d>* get_UnitNormalTriangles() {
		return &m_UnitNormalTriangles;
	}
	inline vector<int> get_NeighborIDs() {
		return m_neighborIDs;
	}

};

class QuadrupleJunction;

class HighOrderJunction: public InterfacialElement {
private:
	Vector3d m_position;
	Vector3d m_bufferPosition;
	//TODO:
public:
	friend class GrainHull;
	HighOrderJunction(int key, GrainHull *owner);
	HighOrderJunction(QuadrupleJunction* A, QuadrupleJunction* B,
			GrainHull *owner);
	~HighOrderJunction();
	void computeEnergy();
	void computeMobility();
	void computePosition();
	void mergeWith(InterfacialElement* B);
	inline Vector3d get_Position() {
		return m_position;
	}
	inline void set_Position(Vector3d newposition) {
		m_position = newposition;
	}
	inline void set_BufferPosition(Vector3d newposition) {
		m_bufferPosition = newposition;
	}
	inline void switch_BufferPosition() {
		m_position = m_bufferPosition;
	}
};

class QuadrupleJunction: public InterfacialElement {
private:
	Vector3d m_position;
	Vector3d m_bufferPosition;
public:
	friend class GrainHull;
	QuadrupleJunction(int key, GrainHull *owner);
	~QuadrupleJunction();
	void computeEnergy();
	void computeMobility();
	void computePosition();
	inline int get_FirstNeighbor() {
		return m_neighborIDs[0];
	}
	;
	inline int get_SecondNeighbor() {
		return m_neighborIDs[1];
	}
	;
	inline int get_ThirdNeighbor() {
		return m_neighborIDs[2];
	}
	;
	inline Vector3d get_Position() {
		return m_position;
	}
	;
	inline void set_Position(Vector3d newposition) {
		m_position = newposition;
	}
	inline void set_BufferPosition(Vector3d newposition) {
		m_bufferPosition = newposition;
	}
	inline void switch_BufferPosition() {
		m_position = m_bufferPosition;
	}
};

class TripleLine: public InterfacialElement {
private:
	vector<InterfacialElement*> m_vertices;
public:
	friend class GrainHull;
	TripleLine(int key, GrainHull *owner);
	TripleLine(int neighbor1, int neighbor2, GrainHull* owner);
	~TripleLine();
	void computeEnergy();
	void computeMobility();
	void findAdjacentJunctions(vector<QuadrupleJunction*>,
			vector<HighOrderJunction*>);
	inline int get_FirstNeighbor() {
		return m_neighborIDs[0];
	}
	inline int get_SecondNeighbor() {
		return m_neighborIDs[1];
	}
	inline vector<InterfacialElement*> get_vertices() {
		return m_vertices;
	}
	inline Vector3d get_Position() {
		return Vector3d(-1, -1, -1);;
	}
	inline void set_Position(Vector3d newposition) {

	}
	inline void set_BufferPosition(Vector3d newposition) {

	}
	inline void switch_BufferPosition() {

	}
};

class GrainBoundary: public InterfacialElement {
private:
	vector<TripleLine*> m_edges; // saves the indexes of edges in clockwise order
	Vector3d inclination;
	double area;
public:
	friend class GrainHull;
	GrainBoundary(int key, GrainHull *owner);
	~GrainBoundary();
	void computeEnergy();
	void computeMobility();
	void findAdjacentTripleLines(vector<TripleLine*>);
	void computeArea();
	inline vector<TripleLine*> get_edges() {
		return m_edges;
	}
	inline double get_area() {
		computeArea();
		return area;
	}
	inline Vector3d get_Position() {
		return Vector3d(-1, -1, -1);
	}
	inline void set_Position(Vector3d newposition) {

	}
	inline void set_BufferPosition(Vector3d newposition) {

	}
	inline void switch_BufferPosition() {

	}
};

#endif
