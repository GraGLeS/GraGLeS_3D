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
#ifndef BOX_H
#define BOX_H

#include "dimensionalBufferIDLocal.h"
#include "dimensionalBufferReal.h"
#include "junction.h"
#include "dimensionalBuffer.h"
#include "pooledDimensionalBufferDouble.h"
#include "spoint.h"
#include "grainHull.h"
#include "charasteristic.h"
#include "triangle.h"
#include "RTree.h"
#include "fftw3.h"
#include "Eigen/Dense"
#include "myQuaternion.h"

#ifdef USE_MKL
#include "mkl_dfti.h"
#endif

using namespace std;
using namespace Eigen;

struct TextureData;
class grainhdl;
class DimensionalBufferReal;
class MarchingSquaresAlgorithm;
class myQuaternion;

enum E_BUFFER_SELECTION {
	E_INPUT_DISTANCE, E_OUTPUT_DISTANCE, E_IDLOCAL, E_INVALID_BUFFER
};

/*!
 * \class LSbox
 * \brief Class encapsulating a Level Set Box.
 *
 * LSbox class contains the coordinates of the box in the actual grid. <br>
 * For each point that the LSbox covers, it stores: <br>
 * - Distances to the actual grain boundary. <br>
 * - The ID of the closest grain.
 *
 * The class also stores the coordinates of the LSbox, Euler angles that represent
 * the orientation, the volume of the grain, <br>
 * the energy of the grain and a pointer to the \b grainhdl object.
 *
 */
class LSbox {
private:
	unsigned int m_ID;
	bool m_exists;
	grainhdl* m_grainHandler;
	GrainHull m_explicitHull;
	bool m_isMotionRegular;
	bool m_intersectsBoundaryGrain;
	DimensionalBufferIDLocal m_IDLocal;
	unsigned int m_neighborCount;
	myQuaternion* m_orientationQuat;
	double m_volume;
	double m_energy;
	double m_surface;
	int m_newXMin;
	int m_newXMax;
	int m_newYMin;
	int m_newYMax;
	int m_newZMin;
	int m_newZMax;
	double m_StoredElasticEnergy;
	SPoint m_centroid;
	DimensionalBufferReal* m_inputDistance;
	DimensionalBufferReal* m_outputDistance;
	vector<unsigned int> m_comparisonList;
	vector<unsigned int> m_secondOrderNeighbours;
	double m_magneticEnergy;
#ifdef USE_FFTW
	fftwp_plan m_backwardsPlan;
	fftwp_plan m_forwardPlan;
#elif defined USE_MKL
	DFTI_DESCRIPTOR_HANDLE m_handle;
	DFTI_DESCRIPTOR_HANDLE m_b_handle;
	MKL_LONG m_dimensions;
#endif
public:
	//Constructors to document
	LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner);
	LSbox(int id, vector<Vector3d>& hull, grainhdl* owner);
	LSbox(int id, const vector<Vector3d>& vertices,
			DimensionalBuffer<int>& IDField, grainhdl* owner);
	LSbox(int id, const vector<Vector3d>& vertices, myQuaternion ori,
			grainhdl* owner, double StoredEnergy);
	//Dtors
	virtual ~LSbox();
	void calculateDistanceFunction(DimensionalBuffer<int>& IDField);
	void executeDERedistancing();
	void executeSurfaceRedistancing();
	void executePreRedistancing();
	void executeFinalRedistancing();
	void executeNewRedistancing();
	void executeRedistancing();
	void extractContour();
	void calculateMagneticEnergy();
	inline bool isMotionRegular() const {
		return m_isMotionRegular;
	}
	void executeComparison();
	double getDistanceFromInputBuff(int i, int j, int k);
	DimensionalBufferReal& getInputDistanceBuffer() {
		return *m_inputDistance;
	}
	void executeSetComparison();
	void computeSecondOrderNeighbours();
	void computeDirectNeighbours(
			const RTree<unsigned int, int, 3, float>& rtree);
	void computeGrainVolume();
	void computeSurfaceArea();
	void computeSurfaceElements();
	void computeVolumeAndEnergy();
	double getGBEnergyTimesGBMobility(int i, int j);
	double getGBEnergyTimesGBMobility(LSbox* neighbour);
	double getGBEnergy(LSbox* neighbour);
	double getWeigthFromHandler(int i, int j);
	void constructBoundarySectors();
	double getWeight(int i, int j, bool minimal = false);
	int getDirectNeighbourCount() {
		return m_neighborCount;
	}
	bool IsNeighbor(int candidate) {
		return m_explicitHull.IsNeighbor(candidate);
	}
	vector<int> getDirectNeighbourIDs();
	vector<double> getGBLengths();
	//map<int, double>& getDiscreteEnergyDistribution() { }
	bool checkIntersection(LSbox* box2);
	void preallocateMemory(ExpandingVector<char>& memory_dump);
	Vector3d findClosestJunctionTo(Vector3d myposition);

	void resizeIDLocalToDistanceBuffer();
	void recalculateIDLocal();
	void setIDLocal(int ID);
	//Debug printing functions
	void plotBoxInterfacialElements(bool absoluteCoordinates = false);
	void plotBoxContour(bool absoluteCoordinates = false);
	void plotNeighboringGrains(bool absoluteCoordinates);
	void plotBoxVolumetric(string identifier,
			E_BUFFER_SELECTION bufferSelection, double h);
	void plotBoxIDLocal();

	double computeMisorientation(LSbox* grain_2);
	double computeMisorientation(unsigned int grainID);
	void resizeGrid(int newSize, double h_old);

	void initConvoMemory(ExpandingVector<char>& memory_dump);
	void createConvolutionPlans(ExpandingVector<char>& memory_dump);
	void executeConvolution(ExpandingVector<char>& mem_pool);

	void correctJunctionPositionWithNeighborInformation();
	void computeInterfacialElementMesh();
	void switchBufferPositions();

	void copyDataToContainer(DimensionalBuffer<unsigned int> * container, int threadID);

#ifdef USE_FFTW
	void makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);
	void makeFFTPlans(float *in, float* out,fftwf_complex *fftTemp, fftwf_plan *fftplan1, fftwf_plan *fftplan2);
	void convolutionGeneratorFFTW(fftwp_complex *fftTemp, fftwp_plan fftplan1, fftwp_plan fftplan2);
	void executeFFTW(fftw_plan fftplan);
	void executeFFTW(fftwf_plan fftplan);
	void cleanupConvolution();
#elif defined USE_MKL
	void convolutionGeneratorMKL(MKL_Complex16* fftTemp);
#endif
	void switchInNOut();
	//todo: refactor with a proper name
	void boundaryCondition();
	inline bool intersectsBoundaryGrain() const {
		return m_intersectsBoundaryGrain;
	}
	void updateFirstOrderNeigbors();
	double GBmobilityModel(double thetaMis);
	bool isNeighbour(LSbox* candidate);
	bool BoundaryIntersection();

	//todo: Analyze if function is required
	vector<double> linearRegression(vector<SPoint>& points2D);
	void calculateTriangleCentroid(vector<SPoint>& triangleCentroid,
			vector<SPoint> triangle);
	void calculateCentroid(SPoint& centroid, vector<GrainJunction> junctions);

	double get_h();

	void outputMemoryUsage(ofstream& output);

	inline bool grainExists() const {
		return m_exists;
	}
	inline double getVolume() const {
		return m_volume;
	}
	inline double getEnergy() const {
		return m_energy;
	}
	inline double getSurface() const {
		return m_surface;
	}
	inline double getMeanWidth() {
		return m_explicitHull.getMeanWidth();
	}
	inline double getTripleLineLength() {
		return m_explicitHull.getTripleLineLength();
	}
	inline unsigned int getID() const {
		return m_ID;
	}
	inline int getMinX() const {
		return m_outputDistance->getMinX();
	}
	inline int getMaxX() const {
		return m_outputDistance->getMaxX();
	}
	inline int getMinY() const {
		return m_outputDistance->getMinY();
	}
	inline int getMaxY() const {
		return m_outputDistance->getMaxY();
	}
	inline int getMinZ() const {
		return m_outputDistance->getMinZ();
	}
	inline int getMaxZ() const {
		return m_outputDistance->getMaxZ();
	}
	inline SPoint getCentroid() const {
		return m_centroid;
	}
	inline const myQuaternion* getOrientationQuat() {
		return m_orientationQuat;
	}
	inline double* getOrientationBunge() {
		return m_orientationQuat->Quaternion2Euler();
	}
	inline grainhdl* get_grainHandler() {
		return m_grainHandler;
	}
	int getNeighbourAt(int i, int j, int k);
	inline double get_magneticEnergy() {
		return m_magneticEnergy;
	}
	inline double get_SEE() {
		return m_StoredElasticEnergy;
	}
	inline double get_NumberOfTriplelines() {
		return m_explicitHull.getNumberOfTripleLines();
	}
	inline double get_NumberOfQuadruplePoints() {
		return m_explicitHull.getNumberOfQuadruplePoints();
	}
	inline double get_NumberOfHighOrderJunctions() {
		return m_explicitHull.getNumberOfHighOrderJunctions();
	}
	inline vector<Face>* get_Faces() {
		return m_explicitHull.get_Faces();
	}
	inline vector<unsigned int> get_Neighbors(){
		return m_secondOrderNeighbours;
	}
	inline vector<Triangle> get_actualHull(){
		return m_explicitHull.get_actualHull();
	}

	inline void setRandomStoredElasticEnergy(){
		m_StoredElasticEnergy = (double)rand()/(double)RAND_MAX * 10e-14;
	}

	inline const vector<unsigned int>& getAllNeighbors(){
		return m_explicitHull.getAllNeighbors();
	}

	TextureData collectTextureData();
};
#endif
