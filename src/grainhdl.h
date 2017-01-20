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

#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "ExpandingVector.h"
#include "dimensionalBuffer.h"
#include "spoint.h"
#include "Settings.h"
#include "Eigen/Dense"
#include <omp.h>
#include <numa.h>



#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
// #define min(x,y) (x<y?x:y)
// #define max(x,y) (x>y?x:y)

using namespace std;
using namespace Eigen;

class LSbox;
class mathMethods;
class myQuaternion;
class IGrainScheduler;
/*!
 * \class grainhdl
 * \brief Class that manages the grain growth simulation.
 */
class grainhdl {
protected:
	int ngrains;
	double dt;
	double h;
	double TimeSlope;

	int realDomainSize;
	int ngridpoints;
	int grid_blowup;

	int Mode;
	IGrainScheduler* m_grainScheduler;
public:
	int currentNrGrains;
	mathMethods* mymath;
	unsigned int loop;
	//! control variable for research mode
	bool loadCurvature;
	unsigned int loadCurvatureLoop;
	bool convolutionCorrection;
	bool calcCentroid;
	bool calcRegression;
	E_RESEARCH_PROJECT project;
	bool constantE;
	vector<double> time;
	vector<double> totalenergy;
	vector<int> nr_grains;
	vector<double> discreteEnergyDistribution;
	DimensionalBuffer<int>* IDField;

	//! A 2D vector which stores weights.
	vector<vector<double> > weightsMatrix;

	double *ST;
	double *part_pos;
	double delta;
	double *bunge;
	double deviation;
	double BoundaryGrainTube;
	double Realtime;
	vector<myQuaternion> myOrientationSpace;
	vector<double> myOrientationSpaceVolumeFracs;
	vector<LSbox*> grains;
	LSbox* boundary;

	grainhdl();
	~grainhdl();

	void initializeSimulation();
	void read_HeaderCPG();
	void readOriFile();

	void VOROMicrostructure();
	void readMicrostructureFromVertex();
	void readMicrostructure();
	void read_voxelized_microstructure();
	void saveMicrostructure();

	void createParamsForSim(const char* param_filename,
			const char* vertex_dum_filename = NULL);

	void find_neighbors();

	void distanceInitialisation();
	void convolution(double& planOverhead);
	void createConvolutionPlans();
	void destroyConvolutionPlans();
	void save_conv_step();
	void comparison_box();
	void countGrains();

	void updateSecondOrderNeighbors();
	void level_set();
	void redistancing();

	void run_sim();
	void save_NrGrainsStats();
	void clear_mem();
	void save_Texture();
	void save_id();
	void plot_contour();
	void gridCoarsement();
	void updateGridAndTimeVariables(double newGridSize) ;
	void find_correctTimestepSize();
	void tweakIDLocal();

	void saveAllSurfaces();
	void saveNetworkAsVoxelContainer();
	void switchDistancebuffer();

	void saveSpecialContourEnergies(int id);
	void removeGrain(int id);
	// 	wrapper functions:

	void set_h(double hn);
	void set_realDomainSize(int realDomainSizen);
	//! Used if points are set manually
	int read_ScenarioPoints();

	inline LSbox* getGrainByID(unsigned int ID) {
		if (ID == 0)
			return boundary;
		else if (ID > 0 && ID < grains.size())
			return grains[ID];
		else
			return NULL;
	}

	inline long get_ngrains() {
		return ngrains;
	}
	inline int get_realDomainSize() {
		return realDomainSize;
	}
	inline int get_ngridpoints() {
		return ngridpoints;
	}
	inline double get_h() {
		return h;
	}
	inline int get_grid_blowup() {
		return grid_blowup;
	}
	inline int get_loop() {
		return loop;
	}
	inline double get_dt() {
		return dt;
	}
	inline double getBoundaryGrainTube() {
		return BoundaryGrainTube;
	}
	inline double get_ds() {
		return 1 / (double)realDomainSize;
	}
	inline double get_TimeSlope() {
		return TimeSlope;
	}

protected:
	void initEnvironment();
	void initNUMABindings();
	void buildBoxVectors(vector<vector<Vector3d>*>& hulls);
	void buildBoxVectors(int* ID, vector<vector<Vector3d>*>& hulls,
			myQuaternion* Quaternionen, double* StoredEnergy);
	int m_ThreadPoolCount;
	vector<ExpandingVector<char> > m_ThreadMemPool;
};
#endif
