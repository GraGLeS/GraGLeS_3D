/*
 * InterfacialElement.cpp
 *
 *  Created on: Nov 24, 2015
 *      Author: cm654063
 */

#include "InterfacialElement.h"
#include "box.h"
#include "grainhdl.h"
#include "marchingCubes.h"
#include "Settings.h"
#include "myQuaternion.h"
#include "grainHull.h"
#include <fstream>

InterfacialElement::InterfacialElement(int key, GrainHull* owner) :
		m_Key_NeighborList(key), m_owner(owner) {
}

InterfacialElement::~InterfacialElement() {
}
double InterfacialElement::computeMobilityMisori(double misori) {
	if (Settings::UseMobilityModel == 0)
		return 1.0;
	if (Settings::UseMobilityModel == 1) {
		double mob = 1.0
				- (0.99 * exp(-5. * (pow(misori / (15 * PI / 180), 9.))));
		if (mob < 0.1)
			return 0.1;
		else
			return mob;
	} else
		return 1.0;
}

double InterfacialElement::computeReadShockleyEnergy(double misori) {
	if (Settings::UniqueGBEnergies == 1)
		return 1.0;
	double gamma_hagb = 1.0;
	double theta_ref = 15 * PI / 180;
	double gamma;
	if (misori > theta_ref)
		gamma = gamma_hagb;
	else
		gamma = gamma_hagb * (misori / theta_ref)
				* (1.0 - log(misori / theta_ref));
	return gamma;

}

GrainBoundary::GrainBoundary(int key, GrainHull* owner) :
		InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	if (currentID != m_owner->m_triangleNeighborLists[key].neighbors[0]) {
		m_neighborIDs.push_back(
				m_owner->m_triangleNeighborLists[key].neighbors[0]);
	} else
		m_neighborIDs.push_back(
				m_owner->m_triangleNeighborLists[key].neighbors[1]);
	computeMobility();
	computeEnergy();
}
GrainBoundary::~GrainBoundary() {
}

void GrainBoundary::computeEnergy() {
	LSbox* neighbor =
			m_owner->m_owner->get_grainHandler()->grains[m_neighborIDs[0]];
	double misori = m_owner->m_owner->computeMisorientation(neighbor);
	m_energy = computeReadShockleyEnergy(misori);
}

void GrainBoundary::computeMobility() {
	LSbox* neighbor =
			m_owner->m_owner->get_grainHandler()->grains[m_neighborIDs[0]];
	double misori = m_owner->m_owner->computeMisorientation(neighbor);
	m_mobility = computeMobilityMisori(misori);
}

void GrainBoundary::computeArea() {
	double surface = 0;
	for (unsigned int i = 0; i < m_Triangles.size(); i++) {
		Triangle& tri = m_Triangles[i];
		Vector3d AB = tri.points[0] - tri.points[1];
		Vector3d BC = tri.points[0] - tri.points[2];
		surface += AB.cross(BC).norm() / 2.0;
	}
	double h = m_owner->m_owner->get_grainHandler()->get_h();

	area = surface * h * h;
}

TripleLine::TripleLine(int key, GrainHull *owner) :
		InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	for (int i = 0; i < 3; i++) {
		if (currentID != m_owner->m_triangleNeighborLists[key].neighbors[i]) {
			m_neighborIDs.push_back(
					m_owner->m_triangleNeighborLists[key].neighbors[i]);
		}
	}
	computeMobility();
	computeEnergy();
}
TripleLine::TripleLine(int neighbor1, int neighbor2, GrainHull* owner) :
		InterfacialElement(-1, owner) {
	m_neighborIDs.push_back(neighbor1);
	m_neighborIDs.push_back(neighbor2);
	computeMobility();
	computeEnergy();
}

TripleLine::~TripleLine() {
}

void TripleLine::computeEnergy() {
	double sigma;
	double gamma[3] = { 0.0, 0.0, 0.0 };
	double theta_mis;
	LSbox* me = m_owner->m_owner;
	grainhdl* handler = m_owner->m_owner->get_grainHandler();
	LSbox* neighborGrains[2] = { handler->getGrainByID(m_neighborIDs[0]),
			handler->getGrainByID(m_neighborIDs[1]) };
	theta_mis = me->computeMisorientation(neighborGrains[0]);
	gamma[0] = computeReadShockleyEnergy(theta_mis);

	theta_mis = neighborGrains[0]->computeMisorientation(neighborGrains[1]);
	gamma[1] = computeReadShockleyEnergy(theta_mis);

	theta_mis = me->computeMisorientation(neighborGrains[1]);
	gamma[2] = computeReadShockleyEnergy(theta_mis);

	sigma = gamma[0] - gamma[1] + gamma[2];

	if (sigma < 0.01) {
//		cout << "negative sigma " << endl;
		sigma = 0.01;
	}
	if (sigma > 1.0)
		sigma = 1.0;
	m_energy = sigma;
}

void TripleLine::computeMobility() {
	double averageMobility = 0;
	double theta_mis;
	LSbox* me = m_owner->m_owner;
	grainhdl* handler = m_owner->m_owner->get_grainHandler();
	LSbox* neighborGrains[2] = { handler->getGrainByID(m_neighborIDs[0]),
			handler->getGrainByID(m_neighborIDs[1]) };
	theta_mis = me->computeMisorientation(neighborGrains[0]);
	averageMobility += computeMobilityMisori(theta_mis);
	theta_mis = neighborGrains[0]->computeMisorientation(neighborGrains[1]);
	averageMobility += computeMobilityMisori(theta_mis);
	theta_mis = me->computeMisorientation(neighborGrains[1]);
	averageMobility += computeMobilityMisori(theta_mis);
	averageMobility /= 3;
//TODO:
	double ds = 3 * sqrt(3 * handler->get_ds() * handler->get_ds());
// ds is the extension of the Tripleline - maximum 3 times the diagonal of a grid cell
	m_mobility = 1.
			/ ((1. / (ds * Settings::TripleLineDrag)) + 1. / averageMobility);

}

void GrainBoundary::findAdjacentTripleLines(vector<TripleLine*> Junctions) {
	int i = 0;
	for (const auto it : Junctions) {
		if (it->get_FirstNeighbor() == m_neighborIDs[0]
				|| it->get_SecondNeighbor() == m_neighborIDs[0]) {
			m_edges.push_back(&(*it));
			i++;
			if (i == 1)
				return;
		}
	}
}
void TripleLine::findAdjacentJunctions(vector<QuadrupleJunction*> JunctionsQ,
		vector<HighOrderJunction*> JunctionsH) {
	int i = 0;
	for (const auto it : JunctionsQ) {
		if (it->get_FirstNeighbor() == m_neighborIDs[0]
				|| it->get_SecondNeighbor() == m_neighborIDs[0]
				|| it->get_ThirdNeighbor() == m_neighborIDs[0])
			if (it->get_FirstNeighbor() == m_neighborIDs[1]
					|| it->get_SecondNeighbor() == m_neighborIDs[1]
					|| it->get_ThirdNeighbor() == m_neighborIDs[1]) {
				m_vertices.push_back(&(*it));
				i++;
				if (i == 2)
					return;
			}
	}
	for (const auto it : JunctionsH) {
		vector<int> neighborIDs = it->get_NeighborIDs();
		bool found_0 = false;
		bool found_1 = false;
		for (int k = 0; k < neighborIDs.size(); k++) {
			if (neighborIDs[k] == m_neighborIDs[0]) {
				found_0 = true;
			}
		}
		for (int k = 0; k < neighborIDs.size(); k++) {
			if (neighborIDs[k] == m_neighborIDs[1]) {
				found_1 = true;
			}
		}
		if (found_0 && found_1) {
			m_vertices.push_back(&(*it));
			i++;
			if (i == 2)
				return;
		}
	}
//if(m_owner->getID() == 1){
	if (i != 2) {
		m_vertices.push_back(NULL);
		//cout << "Tripleline has not found two adjacent junctions. Grain ID:"
		//		<< m_owner->m_owner->getID() << endl;
		//		string filename = string("Grain_") + to_string(
		//				(unsigned long long) m_owner->m_owner->getID()) + string(
		//				"Defect_Tripleline_key_") + to_string(
		//				(unsigned long long) m_Key_NeighborList) + string(".vtk");
		//		ofstream str(filename.c_str());
		//		str << "# vtk DataFile Version 3.0" << endl;
		//		str << "The normal vectors of the surface" << endl;
		//		str << "ASCII" << endl;
		//		str << "DATASET POLYDATA" << endl;
		//		str << endl;
		//		str << "POINTS " << m_Triangles.size() * 3 << " float" << endl;
		//		for (int j = 0; j < m_Triangles.size(); j++) {
		//			str << m_Triangles[j].points[0].x() << " "
		//					<< m_Triangles[j].points[0].y() << " "
		//					<< m_Triangles[j].points[0].z() << endl;
		//			str << m_Triangles[j].points[1].x() << " "
		//					<< m_Triangles[j].points[1].y() << " "
		//					<< m_Triangles[j].points[1].z() << endl;
		//			str << m_Triangles[j].points[2].x() << " "
		//					<< m_Triangles[j].points[2].y() << " "
		//					<< m_Triangles[j].points[2].z() << endl;
		//		}
		//		str << "POLYGONS " << m_Triangles.size() << " " << m_Triangles.size()
		//				* 4 << endl;
		//
		//		for (int j = 0; j < m_Triangles.size(); j++) {
		//			str << "3 " << 3 * j << " " << 3 * j + 1 << " " << 3 * j + 2
		//					<< endl;
		//		}
		//		str.close();
		//		//int time = m_owner->m_owner->get_grainHandler()->get_loop();
		//		//		m_owner->plotContour(true,time);
		//		cout << m_owner->m_TripleLines.size() << endl;
		//for	(auto it : m_owner->m_TripleLines) {
		//		if(it->get_m_Key_NeighborList() ==0)
		//		cout << it->get_m_Key_NeighborList() << " || " << it->get_NeighborIDs()[0] << " || " << it->get_NeighborIDs()[1] << " || " << it->get_NeighborIDs()[2] << " || " << endl;
		//		string filename = string("Grain_") + to_string(
		//						(unsigned long long) m_owner->m_owner->getID()) + string(
		//						"Tripleline_key_") + to_string(
		//						(unsigned long long) m_Key_NeighborList) + string(".vtk");
		//		ofstream str(filename.c_str());
		//		str << "# vtk DataFile Version 3.0" << endl;
		//		str << "The normal vectors of the surface" << endl;
		//		str << "ASCII" << endl;
		//		str << "DATASET POLYDATA" << endl;
		//		str << endl;
		//		str << "POINTS " << it->get_Triangles()->size() * 3 << " float" << endl;
		//		for (int j = 0; j < it->get_Triangles()->size(); j++) {
		//			str << it->get_Triangles(j).points[0].x() << " "
		//			<< it->get_Triangles(j).points[0].y() << " "
		//			<< it->get_Triangles(j).points[0].z() << endl;
		//			str << it->get_Triangles(j).points[1].x() << " "
		//			<< it->get_Triangles(j).points[1].y() << " "
		//			<< it->get_Triangles(j).points[1].z() << endl;
		//			str << it->get_Triangles(j).points[2].x() << " "
		//			<< it->get_Triangles(j).points[2].y() << " "
		//			<< it->get_Triangles(j).points[2].z() << endl;
		//		}
		//		str << "POLYGONS " << it->get_Triangles()->size() << " " << it->get_Triangles()->size()
		//		* 4 << endl;
		//
		//		for (int j = 0; j < it->get_Triangles()->size(); j++) {
		//			str << "3 " << 3 * j << " " << 3 * j + 1 << " " << 3 * j + 2
		//			<< endl;
		//		}
		//		str.close();
		//	}
	}

	return;
}

HighOrderJunction::HighOrderJunction(QuadrupleJunction* A, QuadrupleJunction* B,
		GrainHull *owner) :
		InterfacialElement(A->get_m_Key_NeighborList(), owner) {
	int currentID = owner->m_owner->getID();
	m_neighborIDs.push_back(A->get_FirstNeighbor());
	m_neighborIDs.push_back(A->get_SecondNeighbor());
	m_neighborIDs.push_back(A->get_ThirdNeighbor());
	m_Triangles.reserve(
			A->get_Triangles()->size() + B->get_Triangles()->size());
	m_Triangles.insert(m_Triangles.end(), A->get_Triangles()->begin(),
			A->get_Triangles()->end());
	m_barycenterTriangles.reserve(
			A->get_barycenterTriangles()->size()
					+ B->get_barycenterTriangles()->size());
	m_barycenterTriangles.insert(m_barycenterTriangles.end(),
			A->get_barycenterTriangles()->begin(),
			A->get_barycenterTriangles()->end());

	m_UnitNormalTriangles.reserve(
			A->get_UnitNormalTriangles()->size()
					+ B->get_UnitNormalTriangles()->size());
	m_UnitNormalTriangles.insert(m_UnitNormalTriangles.end(),
			A->get_UnitNormalTriangles()->begin(),
			A->get_UnitNormalTriangles()->end());
	m_Triangles.insert(m_Triangles.end(), B->get_Triangles()->begin(),
			B->get_Triangles()->end());
	m_barycenterTriangles.insert(m_barycenterTriangles.end(),
			B->get_barycenterTriangles()->begin(),
			B->get_barycenterTriangles()->end());
	m_UnitNormalTriangles.insert(m_UnitNormalTriangles.end(),
			B->get_UnitNormalTriangles()->begin(),
			B->get_UnitNormalTriangles()->end());

	bool neighbor = false;
	for (int i = 0; i < 3; i++) {
		if (m_neighborIDs[i] == B->get_FirstNeighbor())
			neighbor = true;
	}
	if (neighbor == true)
		m_neighborIDs.push_back(B->get_FirstNeighbor());
	neighbor = false;
	for (int i = 0; i < 3; i++) {
		if (m_neighborIDs[i] == B->get_SecondNeighbor())
			neighbor = true;
	}
	if (neighbor == true)
		m_neighborIDs.push_back(B->get_SecondNeighbor());
	neighbor = false;
	for (int i = 0; i < 3; i++) {
		if (m_neighborIDs[i] == B->get_ThirdNeighbor())
			neighbor = true;
	}
	if (neighbor == true)
		m_neighborIDs.push_back(B->get_ThirdNeighbor());

	computePosition();
	computeMobility();
	computeEnergy();
}

HighOrderJunction::HighOrderJunction(int key, GrainHull* owner) :
		InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	int j = 0;
	j = m_owner->m_triangleNeighborLists[key].neighbors[0];
	for (int i = 0;
			i < m_owner->m_triangleNeighborLists[key].getNeighborsListCount();
			i++) {
		if (currentID != m_owner->m_triangleNeighborLists[key].neighbors[i]) {
			m_neighborIDs.push_back(
					m_owner->m_triangleNeighborLists[key].neighbors[i]);
			j++;
		}
	}
	computeMobility();
	computeEnergy();
}

HighOrderJunction::~HighOrderJunction() {
}

void HighOrderJunction::computeEnergy() {
	m_energy = 1.0;

}
void HighOrderJunction::computeMobility() {
	m_mobility = 1.0;
}

void HighOrderJunction::mergeWith(InterfacialElement* B) {
	m_Triangles.insert(m_Triangles.end(), B->get_Triangles()->begin(),
			B->get_Triangles()->end());
	m_barycenterTriangles.insert(m_barycenterTriangles.end(),
			B->get_barycenterTriangles()->begin(),
			B->get_barycenterTriangles()->end());
	m_UnitNormalTriangles.insert(m_UnitNormalTriangles.end(),
			B->get_UnitNormalTriangles()->begin(),
			B->get_UnitNormalTriangles()->end());

//merge Neighbor IDs
	for (auto it : B->get_NeighborIDs()) {
		bool found = false;
		for (auto it2 : m_neighborIDs) {
			if (it2 == it)
				found = true;
		}
		m_neighborIDs.push_back(it);
	}
	computePosition();
	computeMobility();
	computeEnergy();
}

void HighOrderJunction::computePosition() {
	m_position = Vector3d(0, 0, 0);
	if (m_barycenterTriangles.size() <= 0) {
		cout << "high order junction has no triangles" << endl;
		char c;
		cin >> c;
	}
	for (int i = 0; i < m_barycenterTriangles.size(); i++) {
		m_position += m_barycenterTriangles[i];
	}
	m_position /= (double) m_barycenterTriangles.size();
}
QuadrupleJunction::QuadrupleJunction(int key, GrainHull* owner) :
		InterfacialElement(key, owner) {
	int currentID = owner->m_owner->getID();
	int j = 0;
	for (int i = 0; i < 4; i++) {
		if (currentID != m_owner->m_triangleNeighborLists[key].neighbors[i]) {
			m_neighborIDs.push_back(
					m_owner->m_triangleNeighborLists[key].neighbors[i]);
			j++;
		}
	}
	computeMobility();
	computeEnergy();
}

QuadrupleJunction::~QuadrupleJunction() {
}

void QuadrupleJunction::computeEnergy() {
//compute average weight for all possible Triplets
	TripleLine T1(m_neighborIDs[0], m_neighborIDs[1], m_owner);
	TripleLine T2(m_neighborIDs[1], m_neighborIDs[2], m_owner);
	TripleLine T3(m_neighborIDs[0], m_neighborIDs[2], m_owner);
	m_energy = (T1.get_energy() + T2.get_energy() + T3.get_energy()) / 3.0;
//	if (m_energy < 1.0)
//		m_energy = 1.0;
}

void QuadrupleJunction::computeMobility() {

	m_mobility = 1.0;
	TripleLine T1(m_neighborIDs[0], m_neighborIDs[1], m_owner);
	TripleLine T2(m_neighborIDs[1], m_neighborIDs[2], m_owner);
	TripleLine T3(m_neighborIDs[0], m_neighborIDs[2], m_owner);
	m_mobility = (T1.get_mobility() + T2.get_mobility() + T3.get_mobility())
			/ 3.0;
//	if (m_mobility < 1.0)
//		m_mobility = 1.0;
}

void QuadrupleJunction::computePosition() {
	m_position = Vector3d(0, 0, 0);
	if (m_barycenterTriangles.size() <= 0) {
		cout << "quadruple junction has no triangles" << endl;
		char c;
		cin >> c;
	}
	for (int i = 0; i < m_barycenterTriangles.size(); i++) {
		m_position += m_barycenterTriangles[i];
	}
	m_position /= (double) m_barycenterTriangles.size();
}

