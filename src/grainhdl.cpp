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
#include "grainhdl.h"
#include "Settings.h"
#include "spoint.h"
#include "RTree.h"
#include "utilities.h"
#include "box.h"
#include "mymath.h"
#include "voro++/include/voro++/voro++.hh"
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"
#include "Eigen/Dense"
#include "IGrainScheduler.h"
#include "myQuaternion.h"
#include "IterativeGrainScheduler.h"
#include "SquaresGrainScheduler.h"
#include <string>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include "Structs.h"
#ifdef USE_MKL
#include "mkl.h"
#endif

using namespace std;

using namespace rapidxml;
using namespace Eigen;
using namespace voro;

grainhdl::grainhdl() :
		m_ThreadPoolCount(0) {
}
grainhdl::~grainhdl() {
	delete mymath;
}

void grainhdl::initializeSimulation() {
	initEnvironment();
	mymath = new mathMethods();
	// 	readInit();

	//! The NumberOfParticles passed via parameters.xml is altered
	//! for MicrostructureGenMode 4
	if (Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE) {
		read_HeaderCPG();
	} else {
		ngrains = Settings::NumberOfParticles;
		currentNrGrains = ngrains;
		realDomainSize =
				(pow((double) Settings::NumberOfParticles, 1. / 3.0)
						* pow(PI * 4 / 3, 1. / 3) / 2
						* Settings::NumberOfPointsPerGrain) + 1;
	}

	discreteEnergyDistribution.resize(Settings::DiscreteSamplingRate);
	fill(discreteEnergyDistribution.begin(), discreteEnergyDistribution.end(),
			0);

	updateGridAndTimeVariables(realDomainSize);

	//Recalculate the setting parameters for the sector radiuses
	//	Settings::ConstantSectorRadius *= h;
	//	Settings::InterpolatingSectorRadius *= h;

	boundary = new LSbox(0, 0, 35 * (PI / 180), 0, this);
	// 	(*boundary).plot_box(false,2,"no.gnu");
	//!grains.resize(Settings::NumberOfParticles + 1);
	grains.resize(Settings::NumberOfParticles + 1);
	grains[0] = boundary;
	switch (Settings::MicrostructureGenMode) {
	case E_GENERATE_WITH_VORONOY: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = 15 * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;
		if (Settings::UseMagneticField)
			readOriFile();
		VOROMicrostructure();
		break;
	}
	case E_READ_VOXELIZED_MICROSTRUCTURE: {
		cout << "Starting to read microstructure input files" << endl;
		read_voxelized_microstructure();
		break;
	}
	default: {
		throw std::runtime_error("Unknown microstructure generation mode!");
	}

	}
	// 	construct_boundary();
	//program options:
	cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
	cout << "Number of Grains: " << ngrains << endl;
	cout << "simulated Timesteps: " << Settings::NumberOfTimesteps << endl;
	cout << "DELTA TUBE: " << delta << endl;
	cout << "Timestepwidth " << dt << endl;
	cout << "Number of Gridpoints: " << realDomainSize << " in [0,1]" << endl
			<< endl;

	cout << endl << "******* start simulation: *******" << endl << endl;
}

void grainhdl::readOriFile() {
	FILE * OriFromFile;
	OriFromFile = fopen(Settings::AdditionalFilename.c_str(), "r");
	int id, N = 0;
	char c;
	// count number of orientations
	do {
		c = fgetc(OriFromFile);
		if (c == '\n')
			N++;
	} while (c != EOF);
	N--;
	rewind(OriFromFile);
	// read over header
	do {
		c = fgetc(OriFromFile);

	} while (c != '\n');

	double vol, euler[3];
	myOrientationSpace.resize(N);
	myOrientationSpaceVolumeFracs.resize(N);
	for (int i = 0; i < N; i++) {
		fscanf(OriFromFile, "%lf \t %lf \t %lf \t %lf\n", &euler[0], &euler[1],
				&euler[2], &vol);
		myOrientationSpace[i].euler2Quaternion(euler);
		myOrientationSpaceVolumeFracs[i] = vol;
	}
}

void grainhdl::read_HeaderCPG() {
	FILE * compressedGrainInfo;
	compressedGrainInfo = fopen(Settings::AdditionalFilename.c_str(), "r");
	if (compressedGrainInfo == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}
	// Read header
	char c;
	string name;
	char buffer[100];
	int I_D;
	int DX, DY, DZ;
	int NX, NY, NZ;

	for (int i = 1; i <= 12; i++) {
		if (i == 1 || i == 12) {
			do {
				c = fgetc(compressedGrainInfo);
			} while (c != '\n');
			continue;
		}
		if (i == 2) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &I_D);
			//cout << buffer << "\t" << I_D << endl;
		}
		if (i == 3) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DX);
			//cout << buffer << "\t" << DX << endl;
		}
		if (i == 4) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DY);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 5) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DZ);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 6) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NX);
			//cout << buffer << "\t" << NX << endl;
		}
		if (i == 7) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NY);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 8) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NZ);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 9) {
			int grains = 0;
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer,
					&Settings::NumberOfParticles);
			//cout << buffer << "\t" << grains << endl;
		}
		if (i == 10 || i == 11)
			fscanf(compressedGrainInfo, "\n");
	}
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	realDomainSize = NX;
	Settings::NumberOfPointsPerGrain = (int) (NX / pow(ngrains, 1 / 3.0) + 0.5);
	// half open container of VORO++
	fclose(compressedGrainInfo);
}

void grainhdl::VOROMicrostructure() {

	timeval time;
	double timer, voro_time;
	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;

	stringstream filename, plotfiles;

	grains.resize(Settings::NumberOfParticles + 1);

	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	if (randbedingung == false)
		realDomainSize -= 1;

	voronoicell_neighbor c;
	int blocks = (int) (pow((Settings::NumberOfParticles / 8), (1 / 3.)) + 1);
	cout << "blocks: " << blocks << endl;
	if (blocks < 1)
		blocks = 1;
	container con(0, 1, 0, 1, 0, 1, blocks, blocks, blocks, randbedingung,
			randbedingung, randbedingung, 10);
	c_loop_all vl(con);

	//Randomly add particles into the container
	double x, y, z;
	for (int i = 0; i < ngrains; i++) {
		x = rnd();
		y = rnd();
		z = rnd();
		/**********************************************************/
		/*
		 * Double pyramid
		 */
//
//			if (i == 0) {
//				x = 0.5;
//				y = 0.5;
//				z = 0.5;
//			} else if (i == 1) {
//				x = 0.25;
//				y = 0.25;
//				z = 0.25;
//			} else if (i == 2) {
//				x = 0.75;
//				y = 0.25;
//				z = 0.25;
//			} else if (i == 3) {
//				x = 0.25;
//				y = 0.75;
//				z = 0.25;
//			} else if (i == 4) {
//				x = 0.75;
//				y = 0.75;
//				z = 0.25;
//			} else if (i == 5) {
//				x = 0.25;
//				y = 0.25;
//				z = 0.75;
//			} else if (i == 6) {
//				x = 0.25;
//				y = 0.75;
//				z = 0.75;
//			} else if (i == 7) {
//				x = 0.75;
//				y = 0.25;
//				z = 0.75;
//			} else if (i == 8) {
//				x = 0.75;
//				y = 0.75;
//				z = 0.75;
//			}
		/*
		 * Tetraeder
		 */
//			if (i == 0) {
//				x = 0.5;
//				y = 0.5;
//				z = 0.5;
//			} else if (i == 1) {
//				x = 0.5;
//				y = 0.5;
//				z = 0.0;
//			} else if (i == 2) {
//				x = 0.5;
//				y = 0.9713;
//				z = 0.6669;
//			} else if (i == 3) {
//				x = 0.319;
//				y = 0.3427;
//				z = 0.6669;
//			} else if (i == 4) {
//				x = 0.6928;
//				y = 0.1337;
//				z = 0.6669;
//			}
		con.put(i, x, y, z);

	}
	vector<vector<Vector3d>*> initialHulls;
	vector<double> cellCoordinates;

	if (vl.start()) {

		initialHulls.resize(ngrains + 1);
		do {
			double cur_x, cur_y, cur_z;
			con.compute_cell(c, vl);
			//new: get the grain_id
			int box_id = vl.pid() + 1;
			vl.pos(cur_x, cur_y, cur_z);
			c.vertices(cur_x, cur_y, cur_z, cellCoordinates);
			initialHulls.at(box_id) = new vector<Vector3d>;
			initialHulls.at(box_id)->reserve(1000);
			for (unsigned int i = 0; i < cellCoordinates.size() / 3; i++) {
				initialHulls.at(box_id)->push_back(
						Vector3d(cellCoordinates.at(3 * i + 1),
								cellCoordinates.at(3 * i),
								cellCoordinates.at(3 * i + 2)));
			}
		} while (vl.inc());

		try {
			IDField = new DimensionalBuffer<int>(0, 0, 0, ngridpoints,
					ngridpoints, ngridpoints);
		} catch (exception& e) {
			cout << "Unable to intialize ID Local! Reason:\n";
			cout << e.what() << endl;
			cout << "Simulation will now halt." << endl;
			exit(2);
			;
		}
		double x, y, z, rx, ry, rz;
		int cell_id;
		for (int k = 0; k < ngridpoints; k++) {
			for (int i = 0; i < ngridpoints; i++) {
				for (int j = 0; j < ngridpoints; j++) {
					x = double((j - grid_blowup) * h);
					y = double((i - grid_blowup) * h);
					z = double((k - grid_blowup) * h);
					if (i < grid_blowup || j < grid_blowup || k < grid_blowup
							|| i >= ngridpoints - 1 - grid_blowup
							|| j >= ngridpoints - 1 - grid_blowup
							|| k >= ngridpoints - 1 - grid_blowup) {
						IDField->setValueAt(i, j, k, 0);
					} else if (con.find_voronoi_cell(x, y, z, rx, ry, rz,
							cell_id)) {
						int box_id = cell_id + 1;
						IDField->setValueAt(i, j, k, box_id);
					} else {
						IDField->setValueAt(i, j, k, 0);
					}
				}
			}
		}

	} else {
		throw runtime_error("Voronoy container error at start() method!");
	}

	gettimeofday(&time, NULL);
	voro_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;
	cout << "Voronoi construction time: " << voro_time << endl;

	buildBoxVectors(initialHulls);
	for (int i = 0; i < initialHulls.size(); i++) {
		delete initialHulls[i];
	}
//	con.draw_particles("VoronoyP.gnu");
//	con.draw_cells_gnuplot("VoronoyC.gnu");
}

void grainhdl::readMicrostructure() {
}

void grainhdl::readMicrostructureFromVertex() {
}

void grainhdl::read_voxelized_microstructure() {
	FILE * compressedGrainInfo;
	compressedGrainInfo = fopen(Settings::AdditionalFilename.c_str(), "rt");
	if (compressedGrainInfo == NULL) {
		cout << "Could not read from file in Settings::AdditionalFilename !";
		exit(2);
	}
	rewind(compressedGrainInfo);
	// Read header

	char c;
	string name;
	char buffer[100];
	int ID_offset;
	int DX, DY, DZ;
	int NX, NY, NZ;

	for (int i = 1; i <= 12; i++) {
		if (i == 1 || i == 12) {
			do {
				c = fgetc(compressedGrainInfo);
			} while (c != '\n');
			continue;
		}
		if (i == 2) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &ID_offset);
			//cout << buffer << "\t" << I_D << endl;
		}
		if (i == 3) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DX);
			//cout << buffer << "\t" << DX << endl;
		}
		if (i == 4) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DY);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 5) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DZ);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 6) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NX);
			//cout << buffer << "\t" << NX << endl;
		}
		if (i == 7) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NY);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 8) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NZ);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 9) {
			int grains = 0;
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &grains);
			//cout << buffer << "\t" << grains << endl;
		}
		if (i == 10 || i == 11)
			fscanf(compressedGrainInfo, "\n");
	}

	double bunge[3];

	int *x, *y, *z;
	int* ID;
	int nvertices = 8;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	int *counts;
	vector<vector<Vector3d>*> vertices;
	vertices.resize(ngrains + 1);
	myQuaternion* Quaternionen = new myQuaternion[ngrains + 1];
	ID = new int[ngrains + 1];
	counts = new int[ngrains + 1];
	x = new int[ngrains + 1];
	y = new int[ngrains + 1];
	z = new int[ngrains + 1];
	double* storedEnergy = new double[ngrains + 1];
	int id = 0;
	for (int nn = 1; nn <= ngrains; nn++) {
		vertices[nn] = new vector<Vector3d>;
		//ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, zmin, ymin, ymax, zmax
		fscanf(
				compressedGrainInfo,
				"%d\t %d\t %d\t %d\t %lf\t %lf\t %lf\t %d\t %d\t %d\t %d\t %d\t %d\t %d \t %lf\n",
				&id, &x[nn], &y[nn], &z[nn], &bunge[0], &bunge[1], &bunge[2],
				&xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &counts[nn],
				&storedEnergy[nn]);
//		printf(
//				"%d\t %d\t %d\t %d\t %lf\t %lf\t %lf\t %d\t %d\t %d\t %d\t %d\t %d\t %d \t %lf\n",
//				id, x[nn], y[nn], z[nn], bunge[0], bunge[1], bunge[2], xmin,
//				xmax, ymin, ymax, zmin, zmax, counts[nn], storedEnergy[nn]);

		if (counts[nn] == 0) {
			ID[nn] = -1;
//			char buf;
//			cin >> buf;
			continue;
		}
		ID[nn] = id - (ID_offset - 1);

		Quaternionen[nn].euler2Quaternion(bunge);
		vertices[nn]->push_back(Vector3d(xmin, ymin, zmin));
		vertices[nn]->push_back(Vector3d(xmin, ymax, zmin));
		vertices[nn]->push_back(Vector3d(xmax, ymin, zmin));
		vertices[nn]->push_back(Vector3d(xmax, ymax, zmin));
		vertices[nn]->push_back(Vector3d(xmin, ymin, zmax));
		vertices[nn]->push_back(Vector3d(xmin, ymax, zmax));
		vertices[nn]->push_back(Vector3d(xmax, ymin, zmax));
		vertices[nn]->push_back(Vector3d(xmax, ymax, zmax));
	}

	fclose(compressedGrainInfo);

	//BINARY read-in

	FILE * voxelized_data;
	long Size = NX * NY * NZ;

	voxelized_data = fopen(Settings::ReadFromFilename.c_str(), "rb");
	if (voxelized_data == NULL) {
		cout
				<< "Could not read from specified file: Settings::ReadFromFilename !";
		exit(2);
	}

	IDField = new DimensionalBuffer<int>(0, 0, 0, ngridpoints, ngridpoints,
			ngridpoints);
	cout << "Allocated an IDfield of " << ngridpoints
			<< "³ grid points consuming "
			<< IDField->getTotalMemoryUsed() / 1.e6 << " MB." << endl;
	for (int k = 0; k < ngridpoints; k++) {
		for (int j = 0; j < ngridpoints; j++) {
			for (int i = 0; i < ngridpoints; i++) {
				//<<<<<<< HEAD
				//				if (i < grid_blowup || j < grid_blowup || k < grid_blowup
				//						|| i >= ngridpoints - 1 - grid_blowup
				//						|| j >= ngridpoints - 1 - grid_blowup
				//						|| k >= ngridpoints - 1 - grid_blowup)
				//=======
				if (i < grid_blowup || j < grid_blowup || k < grid_blowup
						|| i >= ngridpoints - grid_blowup
						|| j >= ngridpoints - grid_blowup
						|| k >= ngridpoints - grid_blowup)
					IDField->setValueAt(j, i, k, 0);
				else {
					unsigned int box_id;
					fread(&box_id, sizeof(unsigned int), 1, voxelized_data);
					box_id = box_id - (ID_offset - 1);
					IDField->setValueAt(j, i, k, box_id);
				}

			}
		}
	}
	//	fstream datei("Voxel.dat", ios::out);
	//	for (int i = 0; i < ngridpoints; i++) {
	//		for (int j = 0; j < ngridpoints; j++) {
	//			datei << IDField->getValueAt(i, j, k) << " | ";
	//		}
	//		datei << endl;
	//	}
	//	fclose(voxelized_data);

	buildBoxVectors(ID, vertices, Quaternionen, storedEnergy);
	for (auto it : vertices) {
		delete it;
	}
	delete[] ID;
	delete[] Quaternionen;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] counts;
	delete[] storedEnergy;
}

void grainhdl::distanceInitialisation() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
#pragma omp critical
		{
			cout << "Thread " << omp_get_thread_num() << " owns "
					<< workload.size() << " grains." << endl;
		}
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] != NULL)
					if (grains[id]->grainExists() != false)
						grains[id]->calculateDistanceFunction(*IDField);
		}
	}
}
void grainhdl::createConvolutionPlans() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] != NULL)
					grains[id]->preallocateMemory(
							m_ThreadMemPool[omp_get_thread_num()]);
		}
	}
#ifdef USE_FFTW
	for(unsigned int i=0; i<Settings::MaximumNumberOfThreads; i++)
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
		for(auto & id : workload)
		if(grains[id] != NULL)
		grains[id]->createConvolutionPlans(m_ThreadMemPool[i]);
	}
#endif
}

void grainhdl::convolution(double& planOverhead) {
	double timer = 0;
	timeval time;

	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;
	createConvolutionPlans();
	gettimeofday(&time, NULL);
	planOverhead += time.tv_sec + time.tv_usec / 1000000.0 - timer;

#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] != NULL)
					grains[id]->executeConvolution(
							m_ThreadMemPool[omp_get_thread_num()]);
		}
	}
#ifdef USE_FFTW
	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;
	for (unsigned int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL)
		grains[i]->cleanupConvolution();
	}
	gettimeofday(&time, NULL);
#endif
	planOverhead += time.tv_sec + time.tv_usec / 1000000.0 - timer;
}

void grainhdl::comparison_box() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] != NULL) {
					grains[id]->executeComparison();
					grains[id]->executeSetComparison();
				}
		}
	}
}

void grainhdl::level_set() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			if (grains[id]->grainExists() == false) {
				delete grains[id];
				grains[id] = NULL;
			} else
				grains[id]->extractContour();
		}
	}
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
//			grains[id]->correctJunctionPositionWithNeighborInformation();
		}
	}
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
//			grains[id]->switchBufferPositions();
		}
	}
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			grains[id]->computeInterfacialElementMesh();
		}
	}
}

void grainhdl::redistancing() {
	currentNrGrains = 0;
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
//		grains[id]->plotBoxVolumetric("Pre",
//				E_INPUT_DISTANCE, get_h());
//		cout << "Surface Redistancing" << endl;
//		grains[id]->executeSurfaceRedistancing();
//		grains[id]->executeDERedistancing();
//		cout << "Pre Redistancing" << endl;
//		grains[id]->executePreRedistancing();
//		cout << "Final Redistancing" << endl;
//		grains[id]->executeFinalRedistancing();
			grains[id]->executeRedistancing();
//		grains[id]->switchInNOut();
//		grains[id]->executeNewRedistancing();
//		grains[id]->executeCombinedRedistancing();
//		grains[id]->switchInNOut();
//		grains[id]->plotBoxVolumetric("Post",
//				E_INPUT_DISTANCE, get_h());
		}
	}
}

void grainhdl::save_Texture() {
	double totalSurface = 0;
	double total_energy = 0.0;
	if (Settings::NeighborTracking && loop > 0) {
		stringstream filename;
		filename << "Texture" << "_" << loop << ".bin";
		FILE* binaryTexture = fopen(filename.str().c_str(), "wb");
		filename.str(std::string());
		filename.clear();
		filename << "Faces" << "_" << loop << ".bin";
		FILE* binaryFaces = fopen(filename.str().c_str(), "wb");

		for (int i = 1; i < grains.size(); i++) {
			if (grains[i] != NULL && grains[i]->grainExists()) {
				totalSurface += grains[i]->getSurface() * 0.5;
				total_energy += grains[i]->getEnergy() * 0.5;
				TextureData to_write = grains[i]->collectTextureData();
				fwrite(&to_write, sizeof(TextureData), 1, binaryTexture);
				vector<Face>* myfaces = grains[i]->get_Faces();
				for (auto it : *myfaces) {
					fwrite(&(it), sizeof(Face), 1, binaryFaces);
				}
				delete myfaces;
			}
		}
		fclose(binaryFaces);
		fclose(binaryTexture);
	} else {
		string filename = string("Texture_")
				+ to_string((unsigned long long) loop) + string(".ori");
		ofstream file;
		file.open(filename.c_str());
		double *euler;
		for (int i = 1; i < grains.size(); i++) {
			if (grains[i] != NULL && grains[i]->grainExists()) {
				totalSurface += grains[i]->getSurface() * 0.5;
				total_energy += grains[i]->getEnergy() * 0.5;
				euler =
						grains[i]->getOrientationQuat()->Quaternion2EulerConst();
				file << grains[i]->getID() << "\t"
						<< grains[i]->getDirectNeighbourCount() << "\t"
						<< grains[i]->intersectsBoundaryGrain() << "\t"
						<< grains[i]->getVolume() << "\t"
						<< grains[i]->getSurface() << "\t"
						<< grains[i]->getEnergy() << "\t"
						<< grains[i]->getMeanWidth() << "\t"
						<< grains[i]->getTripleLineLength() << "\t"
						<< grains[i]->get_NumberOfTriplelines() << "\t"
						<< grains[i]->get_NumberOfQuadruplePoints() << "\t"
						<< grains[i]->get_NumberOfHighOrderJunctions() << "\t";
				if (Settings::UseMagneticField)
					file << grains[i]->get_magneticEnergy();
				else
					file << grains[i]->get_SEE();
				file << "\t" << euler[0] << "\t" << euler[1] << "\t" << euler[2]
						<< "\n";
				delete[] euler;
			}
		}
	}
	nr_grains.push_back(currentNrGrains);
	time.push_back(Realtime);
	totalenergy.push_back(total_energy);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << currentNrGrains
			<< endl;
	cout << "Amount of free Energy in the Network :" << total_energy << endl;
	cout << "Total GB Surface in the Network :" << totalSurface << endl << endl
			<< endl;
}

void grainhdl::run_sim() {
	double parallelRest = 0;
	double convo_time = 0;
	double comparison_time = 0;
	double levelset_time = 0;
	double redistancing_time = 0;
	double plan_overhead = 0;
	double output = 0;

	timeval time;
	double timer;
	distanceInitialisation();
	Realtime = 0;
	find_neighbors();

	for (loop = Settings::StartTime;
			loop <= Settings::StartTime + Settings::NumberOfTimesteps; loop++) {
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		gridCoarsement();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		convolution(plan_overhead);
		gettimeofday(&time, NULL);
		convo_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		switchDistancebuffer();
		updateSecondOrderNeighbors();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		if (Settings::DecoupleGrains != 1) {
			gettimeofday(&time, NULL);
			timer = time.tv_sec + time.tv_usec / 1000000.0;
			comparison_box();
			gettimeofday(&time, NULL);
			comparison_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

			gettimeofday(&time, NULL);
			timer = time.tv_sec + time.tv_usec / 1000000.0;
			switchDistancebuffer();
			gettimeofday(&time, NULL);
			parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;
		} else {
			tweakIDLocal();
		}

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		level_set();
		gettimeofday(&time, NULL);
		levelset_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		redistancing();
		gettimeofday(&time, NULL);
		redistancing_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		countGrains();
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		if (((loop - Settings::StartTime) % int(Settings::AnalysisTimestep))
				== 0 || loop == Settings::NumberOfTimesteps) {
//			saveNetworkState();
			save_Texture();
			save_NrGrainsStats();
			if (loop > 0
					&& ((loop - Settings::StartTime)
							% int(
									Settings::AnalysisTimestep
											* Settings::PlotInterval)) == 0) {
			saveNetworkAsVoxelContainer();
//				plot_contour();
			}
		}
//		plot_contour();
		gettimeofday(&time, NULL);
		output += time.tv_sec + time.tv_usec / 1000000.0 - timer;
		Realtime += (dt
				* (Settings::Physical_Domain_Size
						* Settings::Physical_Domain_Size) / (
//								TimeSlope *
				Settings::HAGB_Energy * Settings::HAGB_Mobility));

		if (currentNrGrains < Settings::BreakupNumber) {
			cout
					<< "Network has coarsed to less than specified by Settings::BreakupNumber. "
					<< "Remaining Grains: " << currentNrGrains
					<< ". Break and save." << endl;
			break;
		}
		// 	utils::CreateMakeGif();
	}
	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << Realtime << endl;
	cout << "Detailed timings: " << endl;
	cout << "Convolution time: " << convo_time << endl;
	cout << "     Of which plan overhead is: " << plan_overhead << endl;
	cout << "Comparison time: " << comparison_time << endl;
	cout << "Redistancing time: " << redistancing_time << endl;
	cout << "Levelset time: " << levelset_time << endl;
	cout << "GridCoarse/SwitchBuffer/UpNeigh: " << parallelRest << endl;
	cout
			<< "Sum parallel regions: "
			<< convo_time + comparison_time + levelset_time + parallelRest
					+ redistancing_time << endl;
	cout << "I/O time: " << output << endl;
}

void grainhdl::countGrains() {
	currentNrGrains = 0;
	for (auto it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		currentNrGrains++;
	}
}

void grainhdl::plot_contour() {
//#pragma omp parallel
//	{
//		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
//				omp_get_thread_num());
//		for (auto id : workload) {
//			if (id <= Settings::NumberOfParticles)
//				if (grains[id] == NULL)
//					continue;
//			grains[id]->plotBoxContour(true);
//		}
//	}

	int id = Settings::NeighbourhoodGrain;
	grains[id]->plotBoxContour(true);
	grains[id]->plotNeighboringGrains(true);
}

void grainhdl::saveMicrostructure() {
}
void grainhdl::createParamsForSim(const char* param_filename,
		const char* vertex_dump_filename) {
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(
			doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(
			doc_tree.allocate_attribute("encoding", "utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(
			Settings::generateXMLParametersNode(&doc_tree, vertex_dump_filename,
					loop, currentNrGrains));
	ofstream output;
	output.open(param_filename);
	output << doc_tree;
	output.close();

}

void grainhdl::save_NrGrainsStats() {
//	ofstream myfile;
//	myfile.open("NrGrains&EnergyStatistics.txt");
//	for (unsigned int i = 0; i < nr_grains.size(); i++) {
//		myfile << time[i] << "\t";
//		myfile << nr_grains[i] << "\t";
//		//myfile << totalenergy[i] << "\t";
//		myfile << realDomainSize << endl;
//	}
//	myfile.close();
//
//}
// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);
	ofstream myfile;
	if (loop == 0)
		myfile.open("NrGrains&EnergyStatistics.txt");
	else
		myfile.open("NrGrains&EnergyStatistics.txt", ios::app);
	myfile << time.back() << "\t";
	myfile << nr_grains.back() << "\t";
	myfile << totalenergy.back() << "\t";
	myfile << realDomainSize << endl;
	myfile.close();

// 	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
//cout << "number of distanzmatrices: "<< domains.size() << endl;

}

void grainhdl::updateSecondOrderNeighbors() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			grains[id]->computeSecondOrderNeighbours();
		}
	}
}

void grainhdl::find_neighbors() {
	RTree<unsigned int, int, 3, float> tree;
	int min[3], max[3];
	for (unsigned int i = 1; i <= Settings::NumberOfParticles; i++) {
		if (grains[i] == NULL || grains[i]->grainExists() == false)
			continue;
		min[0] = grains[i]->getMinX();
		min[1] = grains[i]->getMinY();
		min[2] = grains[i]->getMinZ();
		max[0] = grains[i]->getMaxX();
		max[1] = grains[i]->getMaxY();
		max[2] = grains[i]->getMaxZ();
		tree.Insert(min, max, i);
	}
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			grains[id]->computeDirectNeighbours(tree);
		}
	}
}

void grainhdl::saveSpecialContourEnergies(int id) {
}

void grainhdl::tweakIDLocal() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles) {
				if (grains[id] == NULL)
					continue;
				grains[id]->setIDLocal(boundary->getID());
			}
		}
	}
}

void grainhdl::saveNetworkAsVoxelContainer() {

	DimensionalBuffer<unsigned int> * container = new DimensionalBuffer<
			unsigned int>(0, 0, 0, realDomainSize, realDomainSize,
			realDomainSize);
	string filename = string("Microstructure_")
			+ to_string((unsigned long long) loop) + string(".uds");
	ofstream file;
	file.open(filename.c_str());
	for (int thread = 0; thread < omp_get_max_threads(); thread++) {
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				thread);
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			grains[id]->copyDataToContainer(container, thread);
		}
	}
	for (auto it : grains) {
		if (it == NULL)
			continue;
		if (it->grainExists() == false || it->getID() == 0)
			continue;
		// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
		double *bunge = it->getOrientationBunge();
		file << it->getID() << "\t"
				<< ((it->getMaxX() - it->getMinX()) / 2 + it->getMinX()) * h
				<< "\t"
				<< ((it->getMaxY() - it->getMinY()) / 2 + it->getMinY()) * h
				<< "\t"
				<< ((it->getMaxZ() - it->getMinZ()) / 2 + it->getMinZ()) * h
				<< "\t" << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2]
				<< "\t" << it->getMinX() << "\t" << it->getMaxX() << "\t"
				<< it->getMinY() << "\t" << it->getMaxY() << "\t"
				<< it->getMinZ() << "\t" << it->getMaxZ() << "\t"
				<< it->getVolume() << "\t" << it->get_SEE() << endl;
		delete[] bunge;
	}
	file.close();

	ofstream OutFile;

	stringstream filename2;
	FILE* binaryFile, *File;
	filename2.clear();
	filename2 << "Container_size_" << realDomainSize << "_t_" << loop << ".raw";
	binaryFile = fopen(filename2.str().c_str(), "wb");
	fwrite(container->getRawData(), sizeof(unsigned int),
			(realDomainSize * realDomainSize * realDomainSize), binaryFile);
	fclose(binaryFile);
	filename2.str( std::string() );
	filename2.clear();
	filename2 << "ContainerVolumeEnergy_size_" << realDomainSize << "_t_"
			<< loop << ".raw";
	File = fopen(filename2.str().c_str(), "wb");
	for (int k = 0; k < realDomainSize; k++) {
		for (int i = 0; i < realDomainSize; i++) {
			for (int j = 0; j < realDomainSize; j++) {
				double energy;
				if (Settings::UseStoredElasticEnergy) {
					LSbox* access = getGrainByID(
							container->getValueAt(i, j, k));
					if (access == NULL)
						continue;
					energy = access->get_SEE();
					fwrite(&energy, sizeof(double), 1, binaryFile);
				} else if (Settings::UseMagneticField) {
					LSbox* access = getGrainByID(
							container->getValueAt(i, j, k));
					if (access == NULL)
						continue;
					energy = access->get_magneticEnergy();
					fwrite(&energy, sizeof(double), 1, binaryFile);
				}
			}
		}
	}
	fclose(File);
	delete container;
}

void grainhdl::save_id() {
}

void grainhdl::saveAllSurfaces() {

}
void grainhdl::removeGrain(int id) {
	grains[id] = NULL;
}

void grainhdl::switchDistancebuffer() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] == NULL)
					continue;
			grains[id]->switchInNOut();
		}
	}
}

void grainhdl::clear_mem() {
	if (ST != NULL) {
		delete[] ST;
	}
}

void grainhdl::buildBoxVectors(vector<vector<Vector3d>*>& hulls) {
	m_grainScheduler->buildGrainWorkloads(hulls, ngridpoints);
	bool exceptionHappened = false;
	string error_message;

#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				try {
					LSbox* grain = new LSbox(id, *hulls[id], *IDField, this);
					grains[id] = grain;
				} catch (exception& e) {
#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ")
								+ to_string((unsigned long long) id)
								+ string(" failed at timestep ")
								+ to_string((unsigned long long) loop)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}
		}

		if (exceptionHappened) {
			throw runtime_error(error_message);
		}
	}
}

void grainhdl::buildBoxVectors(int* ID, vector<vector<Vector3d>*>& hulls,
		myQuaternion* Quaternionen, double* storedEnergy) {
	m_grainScheduler->buildGrainWorkloads(hulls, ngridpoints);
	bool exceptionHappened = false;
	string error_message;
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfParticles && ID[id] > 0)
				try {
					LSbox* grain = new LSbox(ID[id], *hulls[id],
							Quaternionen[id], this, storedEnergy[id]);
					grains[id] = grain;
				} catch (exception& e) {
#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ")
								+ to_string((unsigned long long) id)
								+ string(" failed at timestep ")
								+ to_string((unsigned long long) loop)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}
		}

		if (exceptionHappened) {
			throw runtime_error(error_message);
			exit(2);
		}
	}
}

void grainhdl::set_h(double hn) {
	h = hn;
}
void grainhdl::set_realDomainSize(int realDomainSizen) {
	realDomainSize = realDomainSizen;
	ngridpoints = realDomainSize + 2 * grid_blowup;
}
/**
 * This function analyzes the input file for MicrostructureGenMode 4.
 * The amount of lines of the input file is determined. This number
 * indicates the number of points specified in the file, i.e the number of
 * grains. This means one pair of x-y-coordinates each in every line. A
 * tabulator is used as the separator.
 *
 * @return the amount of points in the input file
 */
int grainhdl::read_ScenarioPoints() {

	int counter = 0;
	string line;
	ifstream reader(Settings::ReadFromFilename.c_str());
	while (std::getline(reader, line)) {
		counter++;
	}
	return counter;
}

void grainhdl::initEnvironment() {
//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

	m_ThreadPoolCount = Settings::MaximumNumberOfThreads;
	m_ThreadMemPool.resize(m_ThreadPoolCount);

//These lines might need to be moved if spatial distribution of grains is utilized
//At best the grain scheduler should be configurable through the parameters file

//m_grainScheduler = new IterativeGrainScheduler(Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);

//choose grain scheduler:
	if (Settings::GrainScheduler == E_ITERATIVE) {
		m_grainScheduler = new IterativeGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	} else if (Settings::GrainScheduler == E_SQUARES) {
		m_grainScheduler = new SquaresGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	} else if (Settings::GrainScheduler == E_DEFAULT_SCHEDULER) {
		m_grainScheduler = new IterativeGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	}
	initNUMABindings();
#pragma omp parallel
	{
		double max_size = Settings::NumberOfPointsPerGrain
				* Settings::NumberOfPointsPerGrain * 50;
		int power_of_two = 1 << (int) (ceil(log2(max_size)) + 0.5);
		//!int power_of_two = 1 << (int) (ceil(log2(2<<20)) + 0.5); //!27
		m_ThreadMemPool[omp_get_thread_num()].resize(power_of_two);
	}
}

struct NUMANode {
	int num_cpus;
	int numa_cpus[64];
};

unsigned int my_numa_bitmask_weight(const struct bitmask *mask) {
	unsigned int weight = 0;
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			weight++;
		}
	}
	return weight;
}

void grainhdl::initNUMABindings() {
	vector<NUMANode> nodes;
	nodes.reserve(16);
	numa_available();
// returns a mask of CPUs on which the current task is allowed to run.
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			printf("We are allowed to used node %d\n", j);
			NUMANode node;
			memset(&node, 0xFF, sizeof(node));
			//converts a node number to a bitmask of CPUs.
			//The user must pass a bitmask structure with a mask buffer long enough to represent all possible cpu's
			numa_node_to_cpus(j, cpus);
			node.num_cpus = my_numa_bitmask_weight(cpus);
			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {

				if (numa_bitmask_isbitset(cpus, i)
						&& numa_bitmask_isbitset(numa_all_cpus_ptr, i)) {
					node.numa_cpus[cpuCounter] = i;
					cpuCounter++;
				}
			}
			nodes.push_back(node);
		}
	}
	numa_free_cpumask(cpus);
#pragma omp parallel
	{
		int threadID = omp_get_thread_num();
		for (unsigned int i = 0; i < nodes.size(); i++) {
			if (threadID < nodes.at(i).num_cpus) {
#pragma omp critical
				{
					printf("Will bind thread %d to cpu %d\t",
							omp_get_thread_num(),
							nodes.at(i).numa_cpus[threadID]);
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					printf(res == 0 ? "Done\n" : "Failed\n");
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}

void grainhdl::find_correctTimestepSize() {
	double my_max = 0;
	double my_min = 1000000;
	if (Settings::UseMagneticField) {
		for (int i = 1; i < ngrains; i++) {
			if (grains[i]->get_magneticEnergy() < my_min)
				my_min = grains[i]->get_magneticEnergy();
			if (grains[i]->get_magneticEnergy() > my_max)
				my_max = grains[i]->get_magneticEnergy();
		}
	} else if (Settings::UseStoredElasticEnergy) {
		for (int i = 1; i < ngrains; i++) {
			if (grains[i]->get_SEE() < my_min)
				my_min = grains[i]->get_SEE();
			if (grains[i]->get_SEE() > my_max)
				my_max = grains[i]->get_SEE();
		}
	} else
		return;
	if (ngrains == 1)
		my_min = 0.0;
	double Energy_deltaMAX = (my_max - my_min);
	double m_dt_Correction = 0.5 / (double) realDomainSize / Energy_deltaMAX
			/ dt;
	if (m_dt_Correction > 1.0)
		m_dt_Correction = 1.0;
	dt *= m_dt_Correction;
}

void grainhdl::updateGridAndTimeVariables(double newGridSize) {
	realDomainSize = newGridSize;
	switch (Settings::ConvolutionMode) {
	case E_LAPLACE: {
		dt = 0.8 / double(realDomainSize * realDomainSize);
		break;
	}
	case E_LAPLACE_RITCHARDSON: {
		dt = 0.8 / double(realDomainSize * realDomainSize);
		break;
	}
	case E_GAUSSIAN: {
		dt = Settings::TimeStepSize / double(realDomainSize * realDomainSize);
		TimeSlope = 1 / 1.21;
		break;
	}
	default: {
		throw std::runtime_error("Unknown convolution mode!");
	}
	}
	grid_blowup = Settings::DomainBorderSize;
	delta = Settings::DomainBorderSize / double(realDomainSize);
	ngridpoints = realDomainSize + 2 * grid_blowup;
	h = 1.0 / realDomainSize;
//find_correctTimestepSize();
}

void grainhdl::gridCoarsement() {
	int newSize = pow(currentNrGrains, (1 / 3.0)) * pow(PI * 4 / 3, 1. / 3) / 2
			* Settings::NumberOfPointsPerGrain;
	if ((double) currentNrGrains / (double) ngrains
			< Settings::GridCoarsementGradient && loop != 0
			&& Settings::GridCoarsement && newSize < realDomainSize) {
		double h_old = h;
		updateGridAndTimeVariables(newSize);
		cout << "coarsing the current grid in Timestep: " << loop << endl;
		cout << "newSize :" << newSize << endl << endl;
#pragma omp parallel
		{
			vector<unsigned int>& workload =
					m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				if (id <= Settings::NumberOfParticles)
					if (grains[id] == NULL)
						continue;

				grains[id]->resizeGrid(newSize, h_old);

			}

		}
		//! DISCREPANCY: Compare to the application of dt in the convolution, time decreasing factor 0.8

		ngrains = currentNrGrains;
#pragma omp parallel
		{
			vector<unsigned int>& workload =
					m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				if (id <= Settings::NumberOfParticles)
					if (grains[id] == NULL)
						continue;
				grains[id]->recalculateIDLocal();
			}
		}
#pragma omp parallel
		{
			vector<unsigned int>& workload =
					m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				if (id <= Settings::NumberOfParticles)
					if (grains[id] == NULL)
						continue;
				grains[id]->extractContour();
			}
		}
	} else {
		switchDistancebuffer();
	}
}
