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
#include "box.h"
#include "Settings.h"
#include "dimensionalBufferReal.h"
#include "pooledDimensionalBufferReal.h"
#include "utilities.h"
#include "grainhdl.h"
#include "mymath.h"
#include "marchingCubes.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include "Structs.h"
#include <sys/time.h>

#define PERIODIC(x, f) (((x)+f)%f)

LSbox::LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_surface(0), m_explicitHull(this) {
	m_orientationQuat = new myQuaternion();
	double euler[3] = { phi1, PHI, phi2 };
	m_orientationQuat->euler2Quaternion(euler);
	m_inputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);
	m_outputDistance = new DimensionalBufferReal(0, 0, 0, 0, 0, 0);
	m_magneticEnergy = 0;

	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
}

LSbox::LSbox(int id, vector<Vector3d>& hull, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_surface(0), m_explicitHull(this) {
	int grid_blowup = owner->get_grid_blowup();
	m_magneticEnergy = 0;
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	double h = owner->get_h();
	// determine size of grain
	m_orientationQuat = new myQuaternion();
#pragma omp critical
	{
		if (Settings::UseTexture) {
			double newOri[3];
			(*(m_grainHandler->mymath)).newOrientationFromReference(
					m_grainHandler->bunge, m_grainHandler->deviation, newOri);
			m_orientationQuat->euler2Quaternion(newOri);
		} else
			m_orientationQuat->randomOriShoemakeQuat(m_grainHandler->mymath);

	}
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;
	int zmin = ymin;
	int zmax = 0;

	double x, y, z;
	for (unsigned int k = 0; k < hull.size(); k++) {
		x = hull.at(k)[0];
		y = hull.at(k)[1];
		z = hull.at(k)[2];

		if (y / h < ymin)
			ymin = grid_blowup + y / h;
		if (y / h > ymax)
			ymax = grid_blowup + y / h + 1;
		if (x / h < xmin)
			xmin = grid_blowup + x / h;
		if (x / h > xmax)
			xmax = grid_blowup + x / h + 1;
		if (z / h < zmin)
			zmin = grid_blowup + z / h;
		if (z / h > zmax)
			zmax = grid_blowup + z / h + 1;
	}

	if (xmin == m_grainHandler->get_ngridpoints()) {
		cout << "no bounding box could be found for grain: " << m_ID << endl;
		m_exists = false;
		return;
	}
	xmax += grid_blowup;
	xmin -= grid_blowup;
	ymax += grid_blowup;
	ymin -= grid_blowup;
	zmax += grid_blowup;
	zmin -= grid_blowup;

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);

	m_inputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());

	resizeIDLocalToDistanceBuffer();

}

LSbox::LSbox(int id, const vector<Vector3d>& vertices,
		DimensionalBuffer<int>& IDField, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_surface(0), m_explicitHull(this) {
	int grid_blowup = owner->get_grid_blowup();
	m_magneticEnergy = 0;
	double h = owner->get_h();
	// determine size of grain
	m_orientationQuat = new myQuaternion();
#pragma omp critical
	{
		if (Settings::UseMagneticField) {
			double number;
			int i = 0;
			do {
				number = rnd();
				m_orientationQuat = new myQuaternion();
				*m_orientationQuat = m_grainHandler->myOrientationSpace[i];
				i++;
			} while (number
					>= ((m_grainHandler->myOrientationSpaceVolumeFracs[i])));
		} else {
			if (Settings::UseTexture) {
				double newOri[3];
				(*(m_grainHandler->mymath)).newOrientationFromReference(
						m_grainHandler->bunge, m_grainHandler->deviation,
						newOri);
				m_orientationQuat->euler2Quaternion(newOri);
			} else
				m_orientationQuat->randomOriShoemakeQuat(
						m_grainHandler->mymath);
		}
	}
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;
	int zmin = ymin;
	int zmax = 0;
	//	for (int i = IDField.getMinY(); i < IDField.getMaxY(); i++)
	//		for (int j = IDField.getMinX(); j < IDField.getMaxX(); j++)
	//			for (int k = IDField.getMinZ(); k < IDField.getMaxZ(); k++) {
	//				if (m_ID == IDField.getValueAt(i, j, k)) {
	//					xmax = max(j, xmax);
	//					xmin = min(j, xmin);
	//					ymax = max(i, ymax);
	//					ymin = min(i, ymin);
	//					zmax = max(k, zmax);
	//					zmin = min(k, zmin);
	//				}
	//			}
	//
	//	xmax += grid_blowup;
	//	xmin -= grid_blowup;
	//	ymax += grid_blowup;
	//	ymin -= grid_blowup;
	//	zmax += grid_blowup;
	//	zmin -= grid_blowup;

	int z, y, x;
	for (int k = 0; k < vertices.size(); k++) {
		y = vertices[k][0] / h + 0.5;
		x = vertices[k][1] / h + 0.5;
		z = vertices[k][2] / h + 0.5;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;
	zmax += 2 * grid_blowup;

	if (xmin == m_grainHandler->get_ngridpoints()) {
		cout << "no bounding box could be found for grain: " << m_ID << endl;
		m_exists = false;
		return;
	}

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);

	m_inputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());

	resizeIDLocalToDistanceBuffer();

	m_StoredElasticEnergy = (double)rand()/(double)RAND_MAX * 10e-14;

}

LSbox::LSbox(int id, const vector<Vector3d>& vertices, myQuaternion ori,
		grainhdl* owner, double StoredEnergy) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_explicitHull(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_surface(0) {
	m_StoredElasticEnergy = Settings::DislocEnPerM * StoredEnergy
			/ Settings::HAGB_Energy * Settings::Physical_Domain_Size;
	m_orientationQuat = new myQuaternion(ori.get_q0(), ori.get_q1(),
			ori.get_q2(), ori.get_q3());
	//m_grainBoundary.getRawBoundary() = vertices;
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	int grid_blowup = m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	// determine size of grain
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;
	int zmax = 0;
	int zmin = xmin;

	double z, y, x;
	for (int k = 0; k < vertices.size(); k++) {
//		cout << "vertice xyz: " << vertices[k][0] << "  " << vertices[k][1] << "  "
//				<< vertices[k][2] << "  " << endl;
		x = vertices[k][0];
		y = vertices[k][1];
		z = vertices[k][2];
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (z < zmin)
			zmin = z;
		if (z > zmax)
			zmax = z;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;
	zmax += 2 * grid_blowup;

	if (ymax > m_grainHandler->get_ngridpoints())
		ymax = m_grainHandler->get_ngridpoints();
	if (xmax > m_grainHandler->get_ngridpoints())
		xmax = m_grainHandler->get_ngridpoints();
	if (zmax > m_grainHandler->get_ngridpoints())
		zmax = m_grainHandler->get_ngridpoints();
	if (ymin < 0)
		ymin = 0;
	if (xmin < 0)
		xmin = 0;
	if (zmin < 0)
		zmin = 0;
//	cout << m_ID << "constructed a box with size: " << xmin << "  " << xmax << "  "
//			<< ymin << "  " << ymax << "  "<< zmin << "  " << zmax << "  " << endl;
	if (xmin == m_grainHandler->get_ngridpoints()) {
		cout << m_ID << "no bounding box could be found for grain: " << m_ID
				<< endl;
		m_exists = false;
		return;
	}
	m_inputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, zmin, xmax, ymax,
			zmax);

	m_inputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	resizeIDLocalToDistanceBuffer();
}

LSbox::~LSbox() {
	if (m_orientationQuat != NULL)
		delete m_orientationQuat;
	delete m_inputDistance;
	delete m_outputDistance;
}

double LSbox::get_h() {
	return m_grainHandler->get_h();
}

void LSbox::calculateMagneticEnergy() {
	int i, j;
	double ND[3] = { 0.0, 0.0, 1.0 };
	double *euler = new double[3];
	euler = m_orientationQuat->Quaternion2Euler();
	double p1 = euler[0];
	double t = euler[1];
	double p2 = euler[2];
	double cosine = 0;

	double cAxis[3] = { //Rotated ND in order to represent the unit vector with c orientation
			0.0, 0.0, 0.0 };

	double rotMatrix[9] = { cos(p1) * cos(p2) - sin(p1) * sin(p2) * cos(t), sin(
			p1) * cos(p2) + cos(p1) * sin(p2) * cos(t), sin(p2) * sin(t), -cos(
			p1) * sin(p2) - sin(p1) * cos(p2) * cos(t), -sin(p1) * sin(p2)
			+ cos(p1) * cos(p2) * cos(t), cos(p2) * sin(t), sin(p1) * sin(t),
			-cos(p1) * sin(t), cos(t) };

	double trans[3][3];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			trans[j][i] = rotMatrix[i * 3 + j];

	/*    for(i=0;i<3;i++)
	 printf("%f\t%f\t%f\n",rotMatrix[i][0],rotMatrix[i][1],rotMatrix[i][2]);
	 printf("\n");

	 for(i=0;i<3;i++)
	 printf("%f\t%f\t%f\n",trans[i][0],trans[i][1],trans[i][2]);
	 printf("\n");   */

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			cAxis[i] += trans[i][j] * ND[j];

	cosine = (Settings::MagneticVector_x * cAxis[0])
			+ (Settings::MagneticVector_y * cAxis[1])
			+ (Settings::MagneticVector_z * cAxis[2]);
	if (cosine > 1)
		cosine = 1.0;
	//Both vectors are unit vectors no need to normalize
	//QUICKASSERT(fabs(cosine) <= 1);
	//if different from cosine should be RANDOMCOLOR;
	//See the setColor function for more information
	//printf("cosine=%lf\n",cosine);
	//getchar();
	m_magneticEnergy = 0.5 * Settings::VacuumPermeability
			* Settings::deltaMagSys * Settings::MagneticForceField
			* Settings::MagneticForceField * cosine * cosine;
	m_magneticEnergy = m_magneticEnergy / Settings::HAGB_Energy
			* Settings::Physical_Domain_Size;
	delete[] euler;
}

//TODO: Distance function must be properly implemented using point to mesh distance
void LSbox::calculateDistanceFunction(DimensionalBuffer<int>& IDField) {
	int min = m_grainHandler->get_grid_blowup();
	int max = m_grainHandler->get_ngridpoints() - min - 1;
	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ();
			k++) {
		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
				i++) {
			for (int j = m_inputDistance->getMinX();
					j < m_inputDistance->getMaxX(); j++) {
				if (i < min || i > max || j < min || j > max || k < min
						|| k > max)
					m_inputDistance->setValueAt(i, j, k,
							-m_grainHandler->get_h());
				else if (m_ID == m_grainHandler->IDField->getValueAt(i, j, k))
					m_inputDistance->setValueAt(i, j, k,
							m_grainHandler->get_h());
				else
					m_inputDistance->setValueAt(i, j, k,
							-m_grainHandler->get_h());
			}
		}
	}
	executeRedistancing();
}

// Convolution und Helperfunctions
/**************************************/
/**************************************/

void LSbox::preallocateMemory(ExpandingVector<char>& memory_dump) {
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int desired_size = 0;
#ifdef USE_FFTW
	desired_size = n * n * (floor(n / 2) + 1) * sizeof(fftwp_complex);
#elif defined USE_MKL
	desired_size = n * n * (floor(n / 2) + 1) * sizeof(MKL_Complex16);
#endif
	memory_dump.expand(desired_size);
}

void LSbox::executeConvolution(ExpandingVector<char>& mem_pool) {

	if (grainExists() != true)
		return;
	//  set references for the convolution step
#ifdef USE_FFTW
	fftwp_complex *fftTemp = (fftwp_complex*) &mem_pool[0];
	convolutionGeneratorFFTW(fftTemp, m_forwardPlan, m_backwardsPlan);
#elif defined USE_MKL
	MKL_Complex16* fftTemp = (MKL_Complex16*) &mem_pool[0];
	convolutionGeneratorMKL(fftTemp);
#endif
	resizeIDLocalToDistanceBuffer();
	m_IDLocal.clear();
	if (!Settings::DisableConvolutionCorrection && m_grainHandler->loop != 0
			&& m_isMotionRegular == true) {
		double actualGrainRadius = pow(getVolume() / PI, 1 / 3)
				/ m_grainHandler->get_h();
		double rLimit = 20.0;
		//		if (Settings::ConvolutionMode == E_LAPLACE)
		//			rLimit = 20.0;
		//		else if (Settings::ConvolutionMode == E_LAPLACE_RITCHARDSON)
		//			rLimit = 25.0;
		//		else
		//		if (Settings::ConvolutionMode == E_GAUSSIAN)
		//			rLimit = 20.0;

		int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax,
				intersec_zmin, intersec_zmax;
		if (m_IDLocal.getMinX() < m_outputDistance->getMinX())
			intersec_xmin = m_outputDistance->getMinX();
		else
			intersec_xmin = m_IDLocal.getMinX();

		if (m_IDLocal.getMinY() < m_outputDistance->getMinY())
			intersec_ymin = m_outputDistance->getMinY();
		else
			intersec_ymin = m_IDLocal.getMinY();

		if (m_IDLocal.getMinZ() < m_outputDistance->getMinZ())
			intersec_zmin = m_outputDistance->getMinZ();
		else
			intersec_zmin = m_IDLocal.getMinZ();

		if (m_IDLocal.getMaxX() > m_outputDistance->getMaxX())
			intersec_xmax = m_outputDistance->getMaxX();
		else
			intersec_xmax = m_IDLocal.getMaxX();

		if (m_IDLocal.getMaxY() > m_outputDistance->getMaxY())
			intersec_ymax = m_outputDistance->getMaxY();
		else
			intersec_ymax = m_IDLocal.getMaxY();

		if (m_IDLocal.getMaxZ() > m_outputDistance->getMaxZ())
			intersec_zmax = m_outputDistance->getMaxZ();
		else
			intersec_zmax = m_IDLocal.getMaxZ();

		for (int k = intersec_zmin; k < intersec_zmax; k++)
			for (int i = intersec_ymin; i < intersec_ymax; i++)
				for (int j = intersec_xmin; j < intersec_xmax; j++) {
					double val = m_inputDistance->getValueAt(i, j, k);
					if (abs(val) < m_grainHandler->delta) {
						Vector3d point(i, j, k); //its x y z => j i k
						IDChunkMinimal grain = m_IDLocal.getValueAt(i, j, k);
						double radiuscorrection = 1.0;
						//						if (actualGrainRadius < rLimit) {
						//							const double yInterceptBottom = 0.77;
						//							double radiuscorrection = 1.0;
						//							double const a = 4.5;
						//							radiuscorrection = 1
						//									- exp(-actualGrainRadius / a)
						//											* (1 - yInterceptBottom);
						//						}

						GBInfo localGB(1, 1);
						localGB = m_explicitHull.projectPointToGrainBoundary(
								point, grain.grainID);
						if(localGB.energy* localGB.mobility > 1.0)
						{
							cout << "attention:" << localGB.energy* localGB.mobility  << endl;

						}
						m_outputDistance->setValueAt(i, j, k,
								val
										+ (m_outputDistance->getValueAt(i, j, k)
												- val) * localGB.energy
												* localGB.mobility
												* radiuscorrection);
						if (Settings::UseMagneticField) {
							LSbox* neighbor = m_grainHandler->getGrainByID(
									grain.grainID);
							double f_magneticEnergy = 0;
							if (neighbor != NULL)
								f_magneticEnergy = (m_magneticEnergy
										- neighbor->get_magneticEnergy())
										* localGB.mobility
										* m_grainHandler->get_dt()/ (m_grainHandler->get_TimeSlope()
												*m_grainHandler->get_TimeSlope());

							m_outputDistance->setValueAt(i, j, k,
									(m_outputDistance->getValueAt(i, j, k)
											- f_magneticEnergy));
						}
						if (Settings::UseStoredElasticEnergy) {
							LSbox* neighbor = m_grainHandler->getGrainByID(
									grain.grainID);
							double f_StoredEnergy = 0;
							if (neighbor != NULL)
								f_StoredEnergy = (m_StoredElasticEnergy
										- neighbor->get_SEE())
										* localGB.mobility
										* m_grainHandler->get_dt() / (m_grainHandler->get_TimeSlope()
												*m_grainHandler->get_TimeSlope());

							m_outputDistance->setValueAt(i, j, k,
									(m_outputDistance->getValueAt(i, j, k)
											- f_StoredEnergy));
						}
					}
				}
	}
	resizeIDLocalToDistanceBuffer();
	m_IDLocal.clear();
}

void LSbox::setIDLocal(int ID) {
	IDChunkMinimal SetThisID;
	SetThisID.grainID = ID;
	m_IDLocal.clearValues(SetThisID);
}

double LSbox::getGBEnergyTimesGBMobility(int i, int j) {
	return -1;
}

double LSbox::getGBEnergyTimesGBMobility(LSbox* neighbour) {
	return -1;
}

double LSbox::getGBEnergy(LSbox* neighbour) {
	return -1;
}

void LSbox::resizeIDLocalToDistanceBuffer() {
	int xmaxId = m_outputDistance->getMaxX();
	int xminId = m_outputDistance->getMinX();
	int ymaxId = m_outputDistance->getMaxY();
	int yminId = m_outputDistance->getMinY();
	int zminId = m_outputDistance->getMinZ();
	int zmaxId = m_outputDistance->getMaxZ();
	m_IDLocal.resize(xminId, yminId, zminId, xmaxId, ymaxId, zmaxId);
}
#ifdef USE_FFTW
void LSbox::makeFFTPlans(double *in, double* out, fftw_complex *fftTemp,
		fftw_plan *fftplan1, fftw_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int m = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
	int l = m_outputDistance->getMaxZ() - m_outputDistance->getMinZ();

	*fftplan1 = fftw_plan_dft_r2c_3d(n, m, l, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_3d(n, m, l, fftTemp, out, FFTW_ESTIMATE);

	/*
	 The flags argument is usually either FFTW_MEASURE or FFTW_ESTIMATE. FFTW_MEASURE
	 instructs FFTW to run and measure the execution time of several FFTs in order to find the
	 best way to compute the transform of size n. This process takes some time (usually a few
	 seconds), depending on your machine and on the size of the transform. FFTW_ESTIMATE,
	 on the contrary, does not run any computation and just builds a reasonable plan that is
	 probably sub-optimal. In short, if your program performs many transforms of the same size
	 and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate. */
}

void LSbox::makeFFTPlans(float *in, float* out, fftwf_complex *fftTemp,
		fftwf_plan *fftplan1, fftwf_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int m = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
	int l = m_outputDistance->getMaxZ() - m_outputDistance->getMinZ();

	*fftplan1 = fftwf_plan_dft_r2c_3d(n, m, l, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftwf_plan_dft_c2r_3d(n, m, l, fftTemp, out, FFTW_ESTIMATE);

}

void LSbox::cleanupConvolution() {
	fftw_destroy_planp(m_forwardPlan);
	fftw_destroy_planp(m_backwardsPlan);
}

void LSbox::createConvolutionPlans(ExpandingVector<char>& memory_dump) {
	fftwp_complex *fftTemp = (fftwp_complex*) &memory_dump[0];

	makeFFTPlans(m_inputDistance->getRawData(), m_outputDistance->getRawData(),
			fftTemp, &m_forwardPlan, &m_backwardsPlan);
}

void LSbox::convolutionGeneratorFFTW(fftwp_complex *fftTemp,
		fftwp_plan fftplan1, fftwp_plan fftplan2) {

	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	double dt = m_grainHandler->get_dt();
	int n2 = floor(n / 2) + 1;
	int nn = m_grainHandler->get_realDomainSize();
	double nsq = nn * nn;
	double G;
	int j2;
	int i2;

	double local_nsq =n * n;

	executeFFTW(fftplan1);
	//	Forward DFT

	switch (Settings::ConvolutionMode) {
	case E_LAPLACE: {
		double l = 2.0 * PI / n;
		for (int k = 0; k < n; k++) {
			double coslk = cos(k*l);
			for (int i = 0; i < n; i++) {
				double cosli = cos(l*i);
				for (int j = 0; j < n2; j++) {
//					G = 2.0 * (cosli + coslk + cos(l * j) - 3.0) * nsq;
					G = exp(2.0 * (cosli + coslk + cos(l * j) - 3.0) * nsq / local_nsq) /(local_nsq*n);
//					G = 1.0 / (1.0 + (dt * G)) / (n * n);
					fftTemp[j + n2 * (i + n * k)][0] = fftTemp[j + n2 * (i + n
							* k)][0] * G;
					fftTemp[j + n2 * (i + n * k)][1] = fftTemp[j + n2 * (i + n
							* k)][1] * G;
				}
			}
		}
		break;
	}
	case E_GAUSSIAN: {
		//			Convolution with Normaldistribution
		for (int k = 0; k < n; k++) {
			int k2 = min(k, n - k);
			for (int i = 0; i < n; i++) {
				i2 = min(i, n - i);
				for (int j = 0; j < n2; j++) {
					j2 = min(j, n - j);
					G = exp( -sqrt(2/PI)*
							(i2 * i2 + j2 * j2 + k2 * k2) * 4.0 * dt * nsq
							/ local_nsq * PI * PI) / (local_nsq*n);
					fftTemp[j + n2 * (i + n * k)][0] = fftTemp[j + n2 * (i + n
							* k)][0] * G;
					fftTemp[j + n2 * (i + n * k)][1] = fftTemp[j + n2 * (i + n
							* k)][1] * G;
				}
			}
		}
		break;
		}
		default:
		throw runtime_error("Unknown convolution mode!");
	}

	executeFFTW(fftplan2);
	//	Inverse DFT
}

void LSbox::executeFFTW(fftw_plan fftplan) {
	fftw_execute(fftplan);
}

void LSbox::executeFFTW(fftwf_plan fftplan) {
	fftwf_execute(fftplan);
}

#elif defined USE_MKL
void LSbox::convolutionGeneratorMKL(MKL_Complex16* fftTemp)
{
	if (m_grainHandler->get_loop()==0) {
		m_dimensions=0;
	}

	DFTI_CONFIG_VALUE precision;

#if PRECISION > 0
	precision = DFTI_SINGLE;
#else
	precision = DFTI_DOUBLE;
#endif
	bool update_backward_plan = false;
	if (m_dimensions != m_outputDistance->getMaxX() - m_outputDistance->getMinX()) {
		const MKL_LONG N = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
		m_dimensions = N;
		MKL_LONG dimensions[3] = {N,N,N};
		MKL_LONG input_strides[4] = {0,N,N*N,1};
		MKL_LONG output_strides[4] = {0,(N/2+1),(N/2+1)*N,1};
		const MKL_LONG PO = N*N*N;
		const MKL_LONG SO = N*N*(N/2+1);
		update_backward_plan = true;
		if(m_grainHandler->get_loop()!=0)
		DftiFreeDescriptor(&m_handle);
		DftiCreateDescriptor(&m_handle, precision, DFTI_REAL, 3, dimensions);
		DftiSetValue(m_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		DftiSetValue(m_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
		DftiSetValue(m_handle, DFTI_INPUT_STRIDES, input_strides);
		DftiSetValue(m_handle, DFTI_OUTPUT_STRIDES, output_strides);
		//DftiSetValue(m_handle, DFTI_INPUT_DISTANCE, PO) ;
		//DftiSetValue(m_handle, DFTI_OUTPUT_DISTANCE, SO) ;
		DftiCommitDescriptor(m_handle);
	}
	DftiComputeForward(m_handle, m_inputDistance->getRawData(), fftTemp);

	int n = m_dimensions;
	double dt = m_grainHandler->get_dt();
	int n2 = floor(n / 2) + 1;
	int nn = (*m_grainHandler).get_realDomainSize();
	double nsq = nn *nn;
	double G;
	int j2;
	int i2;

	double local_nsq =n * n;

	switch (Settings::ConvolutionMode) {
	case E_LAPLACE: {
		double l = 2.0 * PI / n;
		for (int k = 0; k < n; k++) {
			double coslk = cos(k*l);
			for (int i = 0; i < n; i++) {
				double cosli = cos(l*i);
				for (int j = 0; j < n2; j++) {
					//					G = 2.0 * (cosli + coslk + cos(l * j) - 3.0) * nsq;
					G = exp(2.0 * (cosli + coslk + cos(l * j) - 3.0) * nsq / local_nsq) /(local_nsq*n);
					//					G = 1.0 / (1.0 + (dt * G)) / (n * n);

					fftTemp[j + n2 * (i + n * k)].real = fftTemp[j + n2 * (i + n * k)].real * G;
					fftTemp[j + n2 * (i + n * k)].imag = fftTemp[j + n2 * (i + n * k)].imag * G;
				}
			}
		}
		break;
	}
	case E_GAUSSIAN: {
		double local_nsq = n * n;
		//			Convolution with Normaldistribution
		for (int k = 0; k < n; k++) {
			int k2 = min(k, n - k);
			for (int i = 0; i < n; i++) {
				i2 = min(i, n - i);
				for (int j = 0; j < n2; j++) {
					j2 = min(j, n - j);
					G = exp(-sqrt(2/PI)*(i2 * i2 + j2 * j2 + k2 * k2) * 4.0 * dt * nsq
							/ local_nsq * PI * PI) / (local_nsq * n);
					fftTemp[j + n2 * (i + n * k)].real = fftTemp[j + n2 * (i + n * k)].real * G;
					fftTemp[j + n2 * (i + n * k)].imag = fftTemp[j + n2 * (i + n * k)].imag * G;
				}
			}
		}
		break;
	}
	default: {
		throw runtime_error("Unknown convolution mode!");
		break;
	}

	}
	if (update_backward_plan) {
		if(m_grainHandler->get_loop()!=0)
		DftiFreeDescriptor(&m_b_handle);
		MKL_LONG N = m_dimensions;
		MKL_LONG dimensions[3] = {N,N,N};
		MKL_LONG input_strides[4] = {0,(N/2+1),(N/2+1)*N,1};
		MKL_LONG output_strides[4] = {0,N,N*N,1};
		const MKL_LONG PO = N*N*N;
		const MKL_LONG SO = N*N*(N/2+1);
		DftiCreateDescriptor(&m_b_handle, precision, DFTI_REAL, 3, dimensions);
		DftiSetValue(m_b_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		DftiSetValue(m_b_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
		DftiSetValue(m_b_handle, DFTI_INPUT_STRIDES, input_strides);
		DftiSetValue(m_b_handle, DFTI_OUTPUT_STRIDES, output_strides);
		//DftiSetValue(m_b_handle, DFTI_INPUT_DISTANCE, PO) ;
		//DftiSetValue(m_b_handle, DFTI_OUTPUT_DISTANCE, SO) ;
		DftiCommitDescriptor(m_b_handle);
	}
	DftiComputeBackward(m_b_handle, fftTemp, m_outputDistance->getRawData());

}
#endif
/**************************************/
/**************************************/

/**************************************/
/**************************************/

void LSbox::switchInNOut() {
	DimensionalBufferReal* temp;
	temp = m_inputDistance;
	m_inputDistance = m_outputDistance;
	m_outputDistance = temp;
}

// Comparison + Helperfunctions
/**************************************/
/**************************************/

void LSbox::executeSetComparison() {
	m_newXMin = m_outputDistance->getMaxX();
	m_newXMax = m_outputDistance->getMinX();
	m_newYMin = m_outputDistance->getMaxY();
	m_newYMax = m_outputDistance->getMinY();
	m_newZMin = m_outputDistance->getMaxZ();
	m_newZMax = m_outputDistance->getMinZ();

	for (int k = m_outputDistance->getMinZ(); k < m_outputDistance->getMaxZ();
			k++) {
		for (int i = m_outputDistance->getMinY();
				i < m_outputDistance->getMaxY(); i++) {
			for (int j = m_outputDistance->getMinX();
					j < m_outputDistance->getMaxX(); j++) {
				if (abs(m_inputDistance->getValueAt(i, j, k))
						< 0.7 * m_grainHandler->delta) {
					m_outputDistance->setValueAt(i, j, k,
							0.5
									* (m_inputDistance->getValueAt(i, j, k)
											- m_outputDistance->getValueAt(i, j,
													k)));
				} else
					m_outputDistance->setValueAt(i, j, k,
							m_inputDistance->getValueAt(i, j, k));

				if (m_outputDistance->getValueAt(i, j, k) >= 0) {
					if (i < m_newYMin)
						m_newYMin = i;
					if (i > m_newYMax)
						m_newYMax = i;
					if (j < m_newXMin)
						m_newXMin = j;
					if (j > m_newXMax)
						m_newXMax = j;
					if (k < m_newZMin)
						m_newZMin = k;
					if (k > m_newZMax)
						m_newZMax = k;
				}
			}
		}
	}

	m_newXMin -= m_grainHandler->get_grid_blowup();
	m_newXMax += m_grainHandler->get_grid_blowup();
	m_newYMin -= m_grainHandler->get_grid_blowup();
	m_newYMax += m_grainHandler->get_grid_blowup();
	m_newZMin -= m_grainHandler->get_grid_blowup();
	m_newZMax += m_grainHandler->get_grid_blowup();

	int ngridpoints = m_grainHandler->get_ngridpoints();

	if(m_newXMin<0)
		m_newXMin=0;
	if(m_newYMin<0)
		m_newYMin=0;
	if(m_newZMin<0)
		m_newZMin=0;
	if(m_newXMax>=ngridpoints)
		m_newXMax = ngridpoints-1;
	if(m_newYMax>=ngridpoints)
		m_newYMax = ngridpoints-1;
	if(m_newZMax>=ngridpoints)
		m_newZMax = ngridpoints-1;
}

bool LSbox::checkIntersection(LSbox* box2) {
	if (m_inputDistance->getMinX() > box2->m_inputDistance->getMaxX()
			|| m_inputDistance->getMaxX() < box2->m_inputDistance->getMinX()
			|| m_inputDistance->getMinY() > box2->m_inputDistance->getMaxY()
			|| m_inputDistance->getMaxY() < box2->m_inputDistance->getMinY()
			|| m_inputDistance->getMinZ() > box2->m_inputDistance->getMaxZ()
			|| m_inputDistance->getMaxZ() < box2->m_inputDistance->getMinZ())
		return false;
	return true;
}

void LSbox::executeComparison() {
	if (grainExists() != true)
		return;
	m_outputDistance->clearValues(-1.0);
	m_secondOrderNeighbours = m_comparisonList;
	int x_min_new, x_max_new, y_min_new, y_max_new, z_min_new, z_max_new;

	for (unsigned int neighs = 0; neighs < m_secondOrderNeighbours.size();
			neighs++) {
		LSbox* neighbor = m_grainHandler->getGrainByID(
				m_secondOrderNeighbours[neighs]);
//		int x_min_new, x_max_new, y_min_new, y_max_new, z_min_new, z_max_new;

		if (m_inputDistance->getMinX() < neighbor->m_inputDistance->getMinX())
			x_min_new = neighbor->m_inputDistance->getMinX();
		else
			x_min_new = m_inputDistance->getMinX();

		if (m_inputDistance->getMaxX() > neighbor->m_inputDistance->getMaxX())
			x_max_new = neighbor->m_inputDistance->getMaxX();
		else
			x_max_new = m_inputDistance->getMaxX();

		if (m_inputDistance->getMinY() < neighbor->m_inputDistance->getMinY())
			y_min_new = neighbor->m_inputDistance->getMinY();
		else
			y_min_new = m_inputDistance->getMinY();

		if (m_inputDistance->getMaxY() > neighbor->m_inputDistance->getMaxY())
			y_max_new = neighbor->m_inputDistance->getMaxY();
		else
			y_max_new = m_inputDistance->getMaxY();

		if (m_inputDistance->getMinZ() < neighbor->m_inputDistance->getMinZ())
			z_min_new = neighbor->m_inputDistance->getMinZ();
		else
			z_min_new = m_inputDistance->getMinZ();

		if (m_inputDistance->getMaxZ() > neighbor->m_inputDistance->getMaxZ())
			z_max_new = neighbor->m_inputDistance->getMaxZ();
		else
			z_max_new = m_inputDistance->getMaxZ();

		for (int k = z_min_new; k < z_max_new; k++) {
			for (int i = y_min_new; i < y_max_new; i++) {
				for (int j = x_min_new; j < x_max_new; j++) {
					double dist = neighbor->getDistanceFromInputBuff(i, j, k);
					if (dist > m_outputDistance->getValueAt(i, j, k)) {
						m_outputDistance->setValueAt(i, j, k, dist);
						m_IDLocal.getValueAt(i, j, k).grainID =
								neighbor->getID();
					}
				}
			}
		}
	}

	vector<int> Neighbors;

	if(m_ID == 1){
		for (int k = z_min_new; k < z_max_new; k++) {
			for (int i = y_min_new; i < y_max_new; i++) {
				for (int j = x_min_new; j < x_max_new; j++) {
					if(find(Neighbors.begin(),Neighbors.end(),m_IDLocal.getValueAt(i,j,k).grainID) == Neighbors.end())
						Neighbors.push_back(m_IDLocal.getValueAt(i,j,k).grainID);
				}
			}
		}
	}

	if (BoundaryIntersection()) {
		m_intersectsBoundaryGrain = true;
		boundaryCondition();
	} else
		m_intersectsBoundaryGrain = false;
}

bool LSbox::BoundaryIntersection() {
	int xMinBoundary = m_grainHandler->get_grid_blowup()
			+ m_grainHandler->getBoundaryGrainTube();
	int yMinBoundary = xMinBoundary;
	int zMinBoundary = xMinBoundary;

	int xMaxBoundary = m_grainHandler->get_ngridpoints()
			- m_grainHandler->get_grid_blowup()
			- m_grainHandler->getBoundaryGrainTube();
	int yMaxBoundary = xMaxBoundary;
	int zMaxBoundary = xMaxBoundary;

	if (m_outputDistance->getMinX() > xMinBoundary
			&& m_outputDistance->getMaxX() < xMaxBoundary
			&& m_outputDistance->getMinY() > yMinBoundary
			&& m_outputDistance->getMaxY() < yMaxBoundary
			&& m_outputDistance->getMinZ() > zMinBoundary
			&& m_outputDistance->getMaxZ() < zMaxBoundary)
		return false;
	else {
		return true;
	}
}

double LSbox::getDistanceFromInputBuff(int i, int j, int k) {
	return m_inputDistance->getValueAt(i, j, k);
}

//Differs largely from latest version of the 2D code! Needs to be checked
void LSbox::boundaryCondition() {
	int grid_blowup = m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	int m = m_grainHandler->get_ngridpoints();
	double distZMin, distZMax, distXMin, distXMax, distX;
	double distYMin, distYMax, distY, distZ;
	double dist = 0;
	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ();
			k++) {
		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
				i++) {
			for (int j = m_inputDistance->getMinX();
					j < m_inputDistance->getMaxX(); j++) {
				distXMin = -(j - grid_blowup);
				distYMin = -(i - grid_blowup);
				distZMin = -(k - grid_blowup);
				distXMax = (j - (m - 1 - grid_blowup));
				distYMax = (i - (m - 1 - grid_blowup));
				distZMax = (k - (m - 1 - grid_blowup));
				////=======
				//	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ(); k++) {
				//		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY(); i++) {
				//			for (int j = m_inputDistance->getMinX(); j
				//					< m_inputDistance->getMaxX(); j++) {
				distXMin = -(j - grid_blowup + 1);
				distYMin = -(i - grid_blowup + 1);
				distZMin = -(k - grid_blowup + 1);
				//				distXMax = (j - (m - grid_blowup));
				//				distYMax = (i - (m - grid_blowup));
				//				distZMax = (k - (m - grid_blowup));

				if (abs(distXMin) < abs(distXMax))
					distX = distXMin;
				else
					distX = distXMax;

				if (abs(distYMin) < abs(distYMax))
					distY = distYMin;
				else
					distY = distYMax;

				if (abs(distZMin) < abs(distZMax))
					distZ = distZMin;
				else
					distZ = distZMax;

				// the point is outside in one of the 8 corners:
				if (distX > 0 && distY > 0 && distZ > 0)
					dist = sqrt(
							(double) distX * distX + distY * distY
									+ distZ * distZ);

				// the point is inside the domain - > value is maximum of the negative distances:
				else if (distX < 0 && distY < 0 && distZ < 0) {
					dist = max(distX, distY);
					dist = max(dist, distZ);
				}

				// the point is outside in x direction
				else if (distX >= 0) {
					if (distY < 0 && distZ < 0) // next to one x-plane
						dist = distX;
					else if (distY < 0 && distZ >= 0)
						dist = sqrt((double) distX * distX + distZ * distZ);
					else if (distY >= 0 && distZ < 0)
						dist = sqrt((double) distX * distX + distY * distY);
				}

				else if (distY >= 0) {
					if (distZ < 0 && distX < 0) // next to one y-plane
						dist = distY;
					else if (distX < 0 && distZ >= 0)
						dist = sqrt((double) distY * distY + distZ * distZ);
				}

				else if (distZ >= 0) {
					if (distY < 0 && distX < 0) // next to one z-plane
						dist = distZ;
				}

				else if (distX == 0 || distY == 0 || distZ == 0)
					dist = 0; // one or more are zero the other lower

				if (dist > -grid_blowup)
					if (dist * h > m_outputDistance->getValueAt(i, j, k)) {
						m_outputDistance->setValueAt(i, j, k, dist * h);
						m_IDLocal.getValueAt(i, j, k).grainID = 0;
					}
			}
		}
	}
}

void LSbox::computeSecondOrderNeighbours() {
	vector<unsigned int> neighbourCandidates;
	if (grainExists() != true)
		return;

	bool just_in;
	m_comparisonList.clear();

	for (unsigned int i = 0; i < m_secondOrderNeighbours.size(); i++) {
		if (m_grainHandler->getGrainByID(m_secondOrderNeighbours[i])->grainExists()
				== true)
			m_comparisonList.push_back(m_secondOrderNeighbours[i]);

		LSbox* currentNeighbor = m_grainHandler->getGrainByID(
				m_secondOrderNeighbours[i]);

		for (unsigned int j = 0;
				j < currentNeighbor->m_secondOrderNeighbours.size(); j++) {
			LSbox* hisNeighbor = m_grainHandler->getGrainByID(
					currentNeighbor->m_secondOrderNeighbours[j]);
			if (hisNeighbor->grainExists() == true)
				if (checkIntersection(hisNeighbor)) {
					neighbourCandidates.push_back(hisNeighbor->getID());
				}
		}
	}
	for (unsigned int i = 0; i < neighbourCandidates.size(); i++) {
		just_in = false;
		if (neighbourCandidates[i] == getID())
			continue;
		if (m_grainHandler->getGrainByID(neighbourCandidates[i])
				== m_grainHandler->boundary)
			continue;

		for (unsigned int j = 0; j < m_comparisonList.size(); j++) {
			if (m_comparisonList[j] == neighbourCandidates[i]) {
				just_in = true;
				break;
			}
		}
		if ((!just_in))
			m_comparisonList.push_back(neighbourCandidates[i]);
	}
	neighbourCandidates.clear();
}

/**************************************/
// end of Comparison
/**************************************/

// Find Contour operates on inputDistance
/**************************************/
/**************************************/

void LSbox::extractContour() {
	m_exists = m_explicitHull.generateHull();
	if (!m_exists) {
		return;
	}
	if (Settings::DecoupleGrains != 1) {
		m_outputDistance->resize(m_newXMin, m_newYMin, m_newZMin, m_newXMax,
				m_newYMax, m_newZMax);
		m_outputDistance->resizeToCube(m_grainHandler->get_ngridpoints());
	}
	m_neighborCount = m_explicitHull.getAllNeighborsCount();
//	if (m_grainHandler->get_loop() > 28)
//		m_explicitHull.plotContour(true, m_grainHandler->get_loop());
	computeGrainVolume();
	computeSurfaceArea();
	computeSurfaceElements();
//	m_explicitHull.plotInterfacialElements(true, m_grainHandler->get_loop());
}

void LSbox::updateFirstOrderNeigbors() {
	if (grainExists() != true)
		return;
	return;
}
void LSbox::computeGrainVolume() {
	m_volume = m_explicitHull.computeGrainVolume();
}

void LSbox::computeSurfaceArea() {
	m_surface = m_explicitHull.computeSurfaceArea();
}
void LSbox::computeSurfaceElements() {

	m_explicitHull.computeGrainBoundaryElements();
	m_explicitHull.subDivideTrianglesToInterfacialElements();
	m_explicitHull.computeJunctionPosition();
}

void LSbox::correctJunctionPositionWithNeighborInformation() {
	m_explicitHull.correctJunctionPositionWithNeighborInformation();
}

void LSbox::computeInterfacialElementMesh() {
	m_explicitHull.computeInterfacialElementMesh();
//	m_explicitHull.plotInterfacialElements(true, m_grainHandler->get_loop());
}

void LSbox::switchBufferPositions() {
	m_explicitHull.switchBufferPositions();
}

Vector3d LSbox::findClosestJunctionTo(Vector3d myposition) {
	Vector3d neighborJunction;
	neighborJunction = m_explicitHull.findClosestJunctionTo(myposition);
	return neighborJunction;
}

void LSbox::computeVolumeAndEnergy() {
	if (grainExists() != true || m_isMotionRegular == false)
		return;
}

/*
 * new redistancing
 */

double pos_sq(double x){
	if(x>0) return x*x;
	else return 0;
}

double neg_sq(double x){
	if(x<0) return x*x;
	else return 0;
}

double max(double x, double y){
	if(x>y) return x;
	else return y;
}

int sign(double val){
	return (0 < val) - (val < 0);
}

double norm_grad(DimensionalBufferReal* sdf, int i, int j, int k, double h){
	double a,b,c,d,e,f;
	double grad_x_sq, grad_y_sq, grad_z_sq;

	int ip = i+1;
	int im = i-1;
	int jp = j+1;
	int jm = j-1;
	int kp = k+1;
	int km = k-1;

	a = sdf->getValueAt(j,i,k) - sdf->getValueAt(j,im,k);
	b = sdf->getValueAt(j,ip,k) - sdf->getValueAt(j,i,k);

	c = sdf->getValueAt(j,i,k) - sdf->getValueAt(jm,i,k);
	d = sdf->getValueAt(jp,i,k) - sdf->getValueAt(j,i,k);

	e = sdf->getValueAt(j,i,k) - sdf->getValueAt(j,i,km);
	f = sdf->getValueAt(j,i,kp) - sdf->getValueAt(j,i,k);

	if(sdf->getValueAt(j,i,k)>0){
		grad_x_sq = max(pos_sq(a), neg_sq(b));
		grad_y_sq = max(pos_sq(c), neg_sq(d));
		grad_z_sq = max(pos_sq(e), neg_sq(f));
	}else{
		grad_x_sq = max(pos_sq(b), neg_sq(a));
		grad_y_sq = max(pos_sq(d), neg_sq(c));
		grad_z_sq = max(pos_sq(f), neg_sq(e));
	}

	return sqrt(grad_x_sq + grad_y_sq + grad_z_sq)/h;
}

void LSbox::executeSurfaceRedistancing(){
	if (grainExists() != true)
			return;

	double tolerance = 0.001;
	double max_dt = tolerance+1;
	double h = m_grainHandler->get_h();
	double dt = 0.95*h;

	int x_Min = m_inputDistance->getMinX();
	int x_Max = m_inputDistance->getMaxX();
	int y_Min = m_inputDistance->getMinY();
	int y_Max = m_inputDistance->getMaxY();
	int z_Min = m_inputDistance->getMinZ();
	int z_Max = m_inputDistance->getMaxZ();

	m_outputDistance->resize(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	DimensionalBufferReal* d = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);
	DimensionalBufferReal* curvature = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	/*
	 * Calculate d for all points in the direct vicinity of the surface
	 */

	timeval time1;
	timeval time2;

	if(getID()==1){
		gettimeofday(&time1, NULL);
	}

	for(int i=x_Min+1; i<x_Max-1;i++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int k=z_Min+1; k<z_Max-1;k++){
				if(abs(m_inputDistance->getValueAt(j,i,k))<m_grainHandler->delta){
					double ip,im,jp,jm,kp,km,current;
					ip=m_inputDistance->getValueAt(j,i+1,k);
					im=m_inputDistance->getValueAt(j,i-1,k);
					jp=m_inputDistance->getValueAt(j+1,i,k);
					jm=m_inputDistance->getValueAt(j-1,i,k);
					kp=m_inputDistance->getValueAt(j,i,k+1);
					km=m_inputDistance->getValueAt(j,i,k-1);
					current=m_inputDistance->getValueAt(j,i,k);

					double ipjp, ipjm, imjp, imjm, jpkp, jpkm, jmkp, jmkm, kpip, kpim, kmip, kmim;

					ipjp=m_inputDistance->getValueAt(j+1,i+1,k);
					ipjm=m_inputDistance->getValueAt(j-1,i+1,k);
					imjp=m_inputDistance->getValueAt(j+1,i-1,k);
					imjm=m_inputDistance->getValueAt(j-1,i-1,k);
					jpkp=m_inputDistance->getValueAt(j+1,i,k+1);
					jpkm=m_inputDistance->getValueAt(j+1,i,k-1);
					jmkp=m_inputDistance->getValueAt(j-1,i,k+1);
					jmkm=m_inputDistance->getValueAt(j-1,i,k-1);
					kpip=m_inputDistance->getValueAt(j,i+1,k+1);
					kpim=m_inputDistance->getValueAt(j,i-1,k+1);
					kmip=m_inputDistance->getValueAt(j,i+1,k-1);
					kmim=m_inputDistance->getValueAt(j,i-1,k-1);

					if(
							current*ip<0 ||
							current*im<0 ||
							current*jp<0 ||
							current*jm<0 ||
							current*kp<0 ||
							current*km<0
					){

						/*
						 * den is the inverse of the norm of the gradient of the level set function
						 */
						double den;
						den = 1/sqrt(pow((ip-im)/(2*h),2)+
								pow((jp-jm)/(2*h),2)+
								pow((kp-km)/(2*h),2));

						/*
						 * in this step the curvature is calculated
						 */

						curvature->setValueAt(j,i,k,
								-(ip+im + jp+jm + kp+km-6*current)/(h*h) +
								pow(den,2) * ( (ip -im)* (4*(ip+im-2*current) +
												ipjp-ipjm-
												imjp+imjm+
												kpip-kmip-
												kpim+kmim)
											  +(jp-jm) * (4*(jp+jm-2*current) +
												ipjp-ipjm-
												imjp+imjm+
												jpkp-jmkp-
												jpkm+jmkm)
											  +(kp-km) * (4*(kp+km-2*current) +
												jpkp-jmkp-
												jpkm+jmkm+
												kpip-kmip-
												kpim+kmim)
										) / (8*pow(h,3))
						);
						d->setValueAt(j,i,k,
								current*den);
					}
				}
			}
		}
	}
//	if(getID()==1){
//			gettimeofday(&time2, NULL);
//			cout << "Time for calculating the curvature:" << endl;
//			cout << time2.tv_sec-time1.tv_sec << ":" << time2.tv_usec-time1.tv_usec << endl;
//	}

	/*
	 * Choose one side of the surface and correct d on this side by using the values of d
	 * on the other side of the surface
	 * This should ensure that the zero isosurface does not shift
	 */

	for(int i=x_Min+1; i<x_Max-1;i++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int k=z_Min+1; k<z_Max-1;k++){
				if(abs(m_inputDistance->getValueAt(j,i,k))<m_grainHandler->delta){
					double ip,im,jp,jm,kp,km,current;
					ip=m_inputDistance->getValueAt(j,i+1,k);
					im=m_inputDistance->getValueAt(j,i-1,k);
					jp=m_inputDistance->getValueAt(j+1,i,k);
					jm=m_inputDistance->getValueAt(j-1,i,k);
					kp=m_inputDistance->getValueAt(j,i,k+1);
					km=m_inputDistance->getValueAt(j,i,k-1);
					current=m_inputDistance->getValueAt(j,i,k);
					if(
							current*ip<0 ||
							current*im<0 ||
							current*jp<0 ||
							current*jm<0 ||
							current*kp<0 ||
							current*km<0
					){
						if(curvature->getValueAt(j,i,k)*current<0){
							int M=0;
							double d_tilde=0;
							if(current*ip<0){
								d_tilde+=d->getValueAt(j,i+1,k)*current/ip;
								M++;
							}
							if(current*im<0){
								d_tilde+=d->getValueAt(j,i-1,k)*current/im;
								M++;
							}
							if(current*jp<0){
								d_tilde+=d->getValueAt(j+1,i,k)*current/jp;
								M++;
							}
							if(current*jm<0){
								d_tilde+=d->getValueAt(j-1,i,k)*current/jm;
								M++;
							}
							if(current*kp<0){
								d_tilde+=d->getValueAt(j,i,k+1)*current/kp;
								M++;
							}
							if(current*km<0){
								d_tilde+=d->getValueAt(j,i,k-1)*current/km;
								M++;
							}
							d->setValueAt(j,i,k,d_tilde/(double)M);
						}
					}
				}
			}
		}
	}

	/*
	 * Set the initial conditions for the corrected level set function
	 */

	for(int i=x_Min; i<x_Max;i++){
		for(int j=y_Min; j<y_Max; j++){
			for(int k=z_Min; k<z_Max;k++){
				m_outputDistance->setValueAt(j,i,k,m_inputDistance->getValueAt(j,i,k));
			}
		}
	}

//	if(getID()==1){
//			gettimeofday(&time1, NULL);
//			cout << "Time for calculating d on one side of the surface:" << endl;
//			cout << time1.tv_sec-time2.tv_sec << ":" << time1.tv_usec-time2.tv_usec << endl;
//	}

	/*
	 * In the direct vicinity of the surface set the level set function to d.
	 * By doing this the gradient of the levelset function is set to one.
	 */

	for(int i=x_Min+1; i<x_Max-1;i++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int k=z_Min+1; k<z_Max-1;k++){
				if(
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i+1,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i-1,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j+1,i,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j-1,i,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k+1)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k-1)<0
				){
					m_outputDistance->setValueAt(j,i,k,d->getValueAt(j,i,k));
				}
			}
		}
	}

	d->clearValues();
	curvature->clearValues();
	delete d;
	delete curvature;
}

void LSbox::executePreRedistancing(){
	if (grainExists() != true)
		return;

	double h = m_grainHandler->get_h();
	double candidate, i_slope, distToZero;

	m_outputDistance->clearValues(-1.0);

	// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax,
	intersec_zmin, intersec_zmax;

	intersec_xmin = m_inputDistance->getMinX();
	intersec_ymin = m_inputDistance->getMinY();
	intersec_zmin = m_inputDistance->getMinZ();
	intersec_xmax = m_inputDistance->getMaxX();
	intersec_ymax = m_inputDistance->getMaxY();
	intersec_zmax = m_inputDistance->getMaxZ();


	// first to updates layer by layer to take advantage of the order of point in memory - there are aligned layer by layer.
	for (int k = intersec_zmin+1; k < intersec_zmax - 1; k++) {
		for (int i = intersec_ymin+1; i < m_outputDistance->getMaxY()-1; i++) {
			for (int j = intersec_xmin+1; j < m_outputDistance->getMaxX() - 1;
					j++) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){
					// x-direction forward
					if (j < intersec_xmax - 1 && i < intersec_ymax) {
						if (m_inputDistance->getValueAt(i, j, k)
								* m_inputDistance->getValueAt(i, j + 1, k) <= 0.0) {
							// interpolate
							i_slope = (m_inputDistance->getValueAt(i, j + 1, k)
									- m_inputDistance->getValueAt(i, j, k)) / h;
							distToZero = -m_inputDistance->getValueAt(i, j, k)
													/ i_slope;
							if (abs(m_outputDistance->getValueAt(i, j, k))
									> abs(distToZero))
								m_outputDistance->setValueAt(i, j, k,
										-distToZero * sgn(i_slope));
						}
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i,
												j + 1, k)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j + 1, k)))
							m_outputDistance->setValueAt(i, j + 1, k, candidate);
					} else {
						candidate = m_outputDistance->getValueAt(i, j, k)
												+ (sgn(m_outputDistance->getValueAt(i, j + 1, k))
														* h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j + 1, k)))
							m_outputDistance->setValueAt(i, j + 1, k, candidate);
					}
				}
			}
		}

		for (int i = intersec_ymin+1; i < m_outputDistance->getMaxY()-1; i++) {
			for (int j = intersec_xmax - 1; j > m_outputDistance->getMinX()+1;
					j--) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){

					// x-direction outputDistanceward
					//check for sign change
					if (j > intersec_xmin && i < intersec_ymax) {
						// calculate new distance candidate and assign if appropriate
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i,
												j - 1, k)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j - 1, k)))
							m_outputDistance->setValueAt(i, j - 1, k, candidate);
					} else {
						//<<<<<<< HEAD
						//					candidate = m_outputDistance->getValueAt(i, j, k)
						//							+ sgn(m_outputDistance->getValueAt(i, j - 1, k))
						//									* h;
						//					if (abs(candidate)
						//							< abs(m_outputDistance->getValueAt(i, j - 1, k)))
						//=======
						candidate = m_outputDistance->getValueAt(i, j, k)
												+ sgn(m_outputDistance->getValueAt(i, j - 1, k))
												* h;
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j - 1, k)))
							//>>>>>>> 66c3d6672001fb0db664a7cc036f13ecf8da0d05
							m_outputDistance->setValueAt(i, j - 1, k, candidate);
					}
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin+1; j < m_outputDistance->getMaxX()-1; j++) {
			for (int i = intersec_ymin+1; i < m_outputDistance->getMaxY() - 1;
					i++) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){
					if (j < intersec_xmax && i < intersec_ymax - 1) {
						if (m_inputDistance->getValueAt(i, j, k)
								* m_inputDistance->getValueAt(i + 1, j, k) <= 0.0) {
							// interpolate
							i_slope = (m_inputDistance->getValueAt(i + 1, j, k)
									- m_inputDistance->getValueAt(i, j, k)) / h;
							distToZero = -m_inputDistance->getValueAt(i, j, k)
													/ i_slope;
							if (abs(m_outputDistance->getValueAt(i, j, k))
									> abs(distToZero))
								m_outputDistance->setValueAt(i, j, k,
										-distToZero * sgn(i_slope));
						}
						// calculate new distance candidate and assign if appropriate
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i + 1,
												j, k)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i + 1, j, k)))
							m_outputDistance->setValueAt(i + 1, j, k, candidate);
					} else {
						candidate = m_outputDistance->getValueAt(i, j, k)
												+ (sgn(m_outputDistance->getValueAt(i + 1, j, k))
														* h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i + 1, j, k)))
							m_outputDistance->setValueAt(i + 1, j, k, candidate);
					}
				}
			}
		}

		for (int j = intersec_xmin+1; j < m_outputDistance->getMaxX()-1; j++) {
			for (int i = intersec_ymax - 1; i > m_outputDistance->getMinY()+1;
					i--) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){
					if (j < intersec_xmax && i > intersec_ymin) {
						// calculate new distance candidate and assign if appropriate
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i - 1,
												j, k)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i - 1, j, k)))
							m_outputDistance->setValueAt(i - 1, j, k, candidate);
					} else {
						candidate = m_outputDistance->getValueAt(i, j, k)
												+ (sgn(m_outputDistance->getValueAt(i - 1, j, k))
														* h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i - 1, j, k)))
							m_outputDistance->setValueAt(i - 1, j, k, candidate);
					}
				}
			}
		}
	}

	// update the point into the third dimensio. the strategy has to change to avoid unneccesary cache loads
	// the idea is to compare all points in one layer first to the next and go on:
	//TODO redist into the third direction
	// z forward:
	for (int k = intersec_zmin+1; k < m_outputDistance->getMaxZ() - 1; k++) {
		for (int i = intersec_ymin+1; i < m_outputDistance->getMaxY()-1; i++) {
			for (int j = intersec_xmin+1; j < m_outputDistance->getMaxX()-1; j++) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){
					// x-direction forward
					if (k < intersec_zmax - 1 && i < intersec_ymax
							&& j < intersec_xmax) {
						if (m_inputDistance->getValueAt(i, j, k)
								* m_inputDistance->getValueAt(i, j, k + 1) <= 0.0) {
							// interpolate
							i_slope = (m_inputDistance->getValueAt(i, j, k + 1)
									- m_inputDistance->getValueAt(i, j, k)) / h;
							distToZero = -m_inputDistance->getValueAt(i, j, k)
													/ i_slope;
							if (abs(m_outputDistance->getValueAt(i, j, k))
									> abs(distToZero))
								m_outputDistance->setValueAt(i, j, k,
										-distToZero * sgn(i_slope));
						}
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i, j,
												k + 1)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j, k + 1)))
							m_outputDistance->setValueAt(i, j, k + 1, candidate);
					} else {
						candidate = m_outputDistance->getValueAt(i, j, k)
												+ (sgn(m_outputDistance->getValueAt(i, j, k + 1))
														* h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j, k + 1)))
							m_outputDistance->setValueAt(i, j, k + 1, candidate);
					}
				}
			}
		}
	}


	// z backward:
	for (int k = intersec_zmax-2; k > m_outputDistance->getMinZ() + 1; k--) {
		for (int i = intersec_ymin+1; i < m_outputDistance->getMaxY()-1; i++) {
			for (int j = intersec_xmin+1; j < m_outputDistance->getMaxX()-1; j++) {
				if(
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j+1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j-1,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i+1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i-1,j,k)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k+1)>0 &&
						m_inputDistance->getValueAt(i,j,k)*m_inputDistance->getValueAt(i,j,k-1)>0
				){
					// x-direction forward
					if (k > intersec_zmin && i < intersec_ymax
							&& j < intersec_xmax) {
						if (m_inputDistance->getValueAt(i, j, k)
								* m_inputDistance->getValueAt(i, j, k - 1) <= 0.0) {
							// interpolate
							i_slope = (m_inputDistance->getValueAt(i, j, k - 1)
									- m_inputDistance->getValueAt(i, j, k)) / h;
							distToZero = -m_inputDistance->getValueAt(i, j, k)
															/ i_slope;
							if (abs(m_outputDistance->getValueAt(i, j, k))
									> abs(distToZero))
								m_outputDistance->setValueAt(i, j, k,
										-distToZero * sgn(i_slope));
						}
						candidate =
								m_outputDistance->getValueAt(i, j, k)
								+ (sgn(
										m_inputDistance->getValueAt(i, j,
												k - 1)) * h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j, k - 1)))
							m_outputDistance->setValueAt(i, j, k - 1, candidate);
					} else {
						candidate = m_outputDistance->getValueAt(i, j, k)
														+ (sgn(m_outputDistance->getValueAt(i, j, k - 1))
																* h);
						if (abs(candidate)
								< abs(m_outputDistance->getValueAt(i, j, k - 1)))
							m_outputDistance->setValueAt(i, j, k - 1, candidate);
					}
				}
			}
		}
	}

}

void LSbox::executeFinalRedistancing(){
	if (grainExists() != true)
		return;

	double tolerance = 0.001;
	double max_dt = tolerance+1;
	double h = m_grainHandler->get_h();
	double dt = 0.95*h;

	int x_Min = m_inputDistance->getMinX();
	int x_Max = m_inputDistance->getMaxX();
	int y_Min = m_inputDistance->getMinY();
	int y_Max = m_inputDistance->getMaxY();
	int z_Min = m_inputDistance->getMinZ();
	int z_Max = m_inputDistance->getMaxZ();

	DimensionalBufferReal* sdf_dt = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	int iter=0;

	while(max_dt > tolerance){

		max_dt = 0.0;
		for(int i=x_Min+1; i<x_Max-1; i++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int k=z_Min+1; k<z_Max-1; k++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					double ip,im,jp,jm,kp,km,current;
					ip=m_inputDistance->getValueAt(j,i+1,k);
					im=m_inputDistance->getValueAt(j,i-1,k);
					jp=m_inputDistance->getValueAt(j+1,i,k);
					jm=m_inputDistance->getValueAt(j-1,i,k);
					kp=m_inputDistance->getValueAt(j,i,k+1);
					km=m_inputDistance->getValueAt(j,i,k-1);
					current=m_inputDistance->getValueAt(j,i,k);

					if(
							current*ip>0 &&
							current*im>0 &&
							current*jp>0 &&
							current*jm>0 &&
							current*kp>0 &&
							current*km>0
					){

						double s = current/
								sqrt(current*
										current + h*h);
						sdf_dt->setValueAt(j,i,k,s *
								(h - norm_grad(m_outputDistance, i, j, k, h)));
						if(
								abs(jp)<m_grainHandler->delta &&
								abs(jm)<m_grainHandler->delta &&
								abs(ip)<m_grainHandler->delta &&
								abs(im)<m_grainHandler->delta &&
								abs(kp)<m_grainHandler->delta &&
								abs(km)<m_grainHandler->delta
						)
							if(abs(sdf_dt->getValueAt(j,i,k))>max_dt) max_dt = abs(sdf_dt->getValueAt(j,i,k));
					}
				}
			}
		}

		for(int i=x_Min+1; i<x_Max-1; i++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int k=z_Min+1; k<z_Max-1; k++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					if(
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i+1,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i-1,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j+1,i,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j-1,i,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k+1)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k-1)>0
					){
						m_outputDistance->setValueAt(j,i,k,
								m_outputDistance->getValueAt(j,i,k)+ dt * sdf_dt->getValueAt(j,i,k));
					}
				}
			}
		}
		iter++;
	}


//	if(getID()==1){
//			gettimeofday(&time1, NULL);
//			cout << "Time for the correction of the remaining level set function:" << endl;
//			cout << time1.tv_sec-time2.tv_sec << ":" << time1.tv_usec-time2.tv_usec << endl;
//
//			cout << "Number of iterations: "<< iter << endl;
//	}

	m_outputDistance->clampValues(-m_grainHandler->delta,
			m_grainHandler->delta);

	sdf_dt->clearValues();
	delete sdf_dt;
}

void LSbox::executeDERedistancing(){
	if (grainExists() != true)
		return;

	double h = m_grainHandler->get_h();
	double tolerance = 0.001*h;
	double max_dt = tolerance+1;
	double dt = 0.95*h;

	int x_Min = m_inputDistance->getMinX();
	int x_Max = m_inputDistance->getMaxX();
	int y_Min = m_inputDistance->getMinY();
	int y_Max = m_inputDistance->getMaxY();
	int z_Min = m_inputDistance->getMinZ();
	int z_Max = m_inputDistance->getMaxZ();

	m_outputDistance->resize(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	DimensionalBufferReal* sdf_dt = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	for(int i=x_Min; i<x_Max;i++){
		for(int j=y_Min; j<y_Max; j++){
			for(int k=z_Min; k<z_Max;k++){
				m_outputDistance->setValueAt(j,i,k,m_inputDistance->getValueAt(j,i,k));
			}
		}
	}

	int iter=0;

	while(max_dt > tolerance){

		max_dt = 0.0;
		for(int i=x_Min+1; i<x_Max-1; i++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int k=z_Min+1; k<z_Max-1; k++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					double ip,im,jp,jm,kp,km,current;
					ip=m_inputDistance->getValueAt(j,i+1,k);
					im=m_inputDistance->getValueAt(j,i-1,k);
					jp=m_inputDistance->getValueAt(j+1,i,k);
					jm=m_inputDistance->getValueAt(j-1,i,k);
					kp=m_inputDistance->getValueAt(j,i,k+1);
					km=m_inputDistance->getValueAt(j,i,k-1);
					current=m_inputDistance->getValueAt(j,i,k);


					double s = current/
							sqrt(current*
									current + h*h);
					sdf_dt->setValueAt(j,i,k,s *
							(h - norm_grad(m_outputDistance, i, j, k, h)));
					if(
							abs(jp)<m_grainHandler->delta &&
							abs(jm)<m_grainHandler->delta &&
							abs(ip)<m_grainHandler->delta &&
							abs(im)<m_grainHandler->delta &&
							abs(kp)<m_grainHandler->delta &&
							abs(km)<m_grainHandler->delta
					)
						if(abs(sdf_dt->getValueAt(j,i,k)/s)>max_dt) max_dt = abs(sdf_dt->getValueAt(j,i,k)/s);
				}
			}
		}
		cout << max_dt/h << endl;

		for(int i=x_Min+1; i<x_Max-1; i++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int k=z_Min+1; k<z_Max-1; k++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					m_outputDistance->setValueAt(j,i,k,
							m_outputDistance->getValueAt(j,i,k)+ dt * sdf_dt->getValueAt(j,i,k));
				}
			}
		}
		iter++;
	}


	if(getID()==1){
//			gettimeofday(&time1, NULL);
//			cout << "Time for the correction of the remaining level set function:" << endl;
//			cout << time1.tv_sec-time2.tv_sec << ":" << time1.tv_usec-time2.tv_usec << endl;

			cout << "Number of iterations: "<< iter << endl;
	}

	m_outputDistance->clampValues(-m_grainHandler->delta,
			m_grainHandler->delta);

	sdf_dt->clearValues();
	delete sdf_dt;
}

void LSbox::executeNewRedistancing(){
	if (grainExists() != true)
			return;

	double tolerance = 0.0001;
	double max_dt = tolerance+1;
	double h = m_grainHandler->get_h();
	double dt = 0.75*h;

	int x_Min = m_inputDistance->getMinX();
	int x_Max = m_inputDistance->getMaxX();
	int y_Min = m_inputDistance->getMinY();
	int y_Max = m_inputDistance->getMaxY();
	int z_Min = m_inputDistance->getMinZ();
	int z_Max = m_inputDistance->getMaxZ();

	m_outputDistance->resize(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	DimensionalBufferReal* d = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);
	DimensionalBufferReal* curvature = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	/*
	 * Calculate d for all points in the direct vicinity of the surface
	 */

	timeval time1;
	timeval time2;

	if(getID()==1){
		gettimeofday(&time1, NULL);
	}

	for(int k=z_Min+1; k<z_Max-1;k++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int i=x_Min+1; i<x_Max-1;i++){
				if(abs(m_inputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
					continue;

				double ip,im,jp,jm,kp,km,current;
				ip=m_inputDistance->getValueAt(j,i+1,k);
				im=m_inputDistance->getValueAt(j,i-1,k);
				jp=m_inputDistance->getValueAt(j+1,i,k);
				jm=m_inputDistance->getValueAt(j-1,i,k);
				kp=m_inputDistance->getValueAt(j,i,k+1);
				km=m_inputDistance->getValueAt(j,i,k-1);
				current=m_inputDistance->getValueAt(j,i,k);

				double ipjp, ipjm, imjp, imjm, jpkp, jpkm, jmkp, jmkm, kpip, kpim, kmip, kmim;

				ipjp=m_inputDistance->getValueAt(j+1,i+1,k);
				ipjm=m_inputDistance->getValueAt(j-1,i+1,k);
				imjp=m_inputDistance->getValueAt(j+1,i-1,k);
				imjm=m_inputDistance->getValueAt(j-1,i-1,k);
				jpkp=m_inputDistance->getValueAt(j+1,i,k+1);
				jpkm=m_inputDistance->getValueAt(j+1,i,k-1);
				jmkp=m_inputDistance->getValueAt(j-1,i,k+1);
				jmkm=m_inputDistance->getValueAt(j-1,i,k-1);
				kpip=m_inputDistance->getValueAt(j,i+1,k+1);
				kpim=m_inputDistance->getValueAt(j,i-1,k+1);
				kmip=m_inputDistance->getValueAt(j,i+1,k-1);
				kmim=m_inputDistance->getValueAt(j,i-1,k-1);

				if(
						current*ip>0 &&
						current*im>0 &&
						current*jp>0 &&
						current*jm>0 &&
						current*kp>0 &&
						current*km>0
				)
					continue;

				/*
				 * den is the inverse of the norm of the gradient of the level set function
				 */
				double den;
				den = 1/sqrt(pow((ip-im)/(2*h),2)+
						pow((jp-jm)/(2*h),2)+
						pow((kp-km)/(2*h),2));

				/*
				 * in this step the curvature is calculated
				 */

				curvature->setValueAt(j,i,k,
						-(ip+im + jp+jm + kp+km-6*current)/(h*h) +
						pow(den,2) * ( (ip -im)* (4*(ip+im-2*current) +
								ipjp-ipjm-
								imjp+imjm+
								kpip-kmip-
								kpim+kmim)
								+(jp-jm) * (4*(jp+jm-2*current) +
										ipjp-ipjm-
										imjp+imjm+
										jpkp-jmkp-
										jpkm+jmkm)
										+(kp-km) * (4*(kp+km-2*current) +
												jpkp-jmkp-
												jpkm+jmkm+
												kpip-kmip-
												kpim+kmim)
						) / (8*pow(h,3))
				);
				d->setValueAt(j,i,k,
						current*den);
			}
		}
	}
	if(getID()==1){
		gettimeofday(&time2, NULL);
		cout << "Time for calculating the curvature:" << endl;
		cout << time2.tv_sec-time1.tv_sec << ":" << time2.tv_usec-time1.tv_usec << endl;
	}

	/*
	 * Choose one side of the surface and correct d on this side by using the values of d
	 * on the other side of the surface
	 * This should ensure that the zero isosurface does not shift
	 */

	for(int k=z_Min+1; k<z_Max-1;k++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int i=x_Min+1; i<x_Max-1;i++){
				if(abs(m_inputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
					continue;
				double ip,im,jp,jm,kp,km,current;
				ip=m_inputDistance->getValueAt(j,i+1,k);
				im=m_inputDistance->getValueAt(j,i-1,k);
				jp=m_inputDistance->getValueAt(j+1,i,k);
				jm=m_inputDistance->getValueAt(j-1,i,k);
				kp=m_inputDistance->getValueAt(j,i,k+1);
				km=m_inputDistance->getValueAt(j,i,k-1);
				current=m_inputDistance->getValueAt(j,i,k);
				if(
						current*ip>0 &&
						current*im>0 &&
						current*jp>0 &&
						current*jm>0 &&
						current*kp>0 &&
						current*km>0
				)
					continue;
				if(curvature->getValueAt(j,i,k)*current<0){
					int M=0;
					double d_tilde=0;
					if(current*ip>0){
					}else{
						d_tilde+=d->getValueAt(j,i+1,k)*current/ip;
						M++;
					}
					if(current*im>0){
					}else{
						d_tilde+=d->getValueAt(j,i-1,k)*current/im;
						M++;
					}
					if(current*jp>0){
					}else{
						d_tilde+=d->getValueAt(j+1,i,k)*current/jp;
						M++;
					}
					if(current*jm>0){
					}else{
						d_tilde+=d->getValueAt(j-1,i,k)*current/jm;
						M++;
					}
					if(current*kp>0){
					}else{
						d_tilde+=d->getValueAt(j,i,k+1)*current/kp;
						M++;
					}
					if(current*km>0){
					}else{
						d_tilde+=d->getValueAt(j,i,k-1)*current/km;
						M++;
					}
					d->setValueAt(j,i,k,d_tilde/(double)M);
				}
			}
		}
	}

	/*
	 * Set the initial conditions for the corrected level set function
	 */

	for(int k=z_Min; k<z_Max;k++){
		for(int j=y_Min; j<y_Max; j++){
			for(int i=x_Min; i<x_Max;i++){
				m_outputDistance->setValueAt(j,i,k,m_inputDistance->getValueAt(j,i,k));
			}
		}
	}

	if(getID()==1){
			gettimeofday(&time1, NULL);
			cout << "Time for calculating d on one side of the surface:" << endl;
			cout << time1.tv_sec-time2.tv_sec << ":" << time1.tv_usec-time2.tv_usec << endl;
	}

	/*
	 * In the direct vicinity of the surface set the level set function to d.
	 * By doing this the gradient of the levelset function is set to one.
	 */

	for(int k=z_Min+1; k<z_Max-1;k++){
		for(int j=y_Min+1; j<y_Max-1; j++){
			for(int i=x_Min+1; i<x_Max-1;i++){
				if(
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i+1,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i-1,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j+1,i,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j-1,i,k)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k+1)<0 ||
						m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k-1)<0
				){
					m_outputDistance->setValueAt(j,i,k,d->getValueAt(j,i,k));
				}
			}
		}
	}



	if(getID()==1){
			gettimeofday(&time2, NULL);
			cout << "Time for correcting the levelset function in the vicinity of the surface:" << endl;
			cout << time2.tv_sec-time1.tv_sec << ":" << time2.tv_usec-time1.tv_usec << endl;
	}

	DimensionalBufferReal* sdf_dt = new DimensionalBufferReal(x_Min,y_Min,z_Min,x_Max,y_Max,z_Max);

	/*
	 * Set the gradient of the level set function everywhere else to 1.
	 */

	int iter=0;

	while(max_dt > tolerance){

//		string Gradient;
//		Gradient = "Gradient";
//		Gradient += to_string(iter);
//		ofstream GradientFile;

//		if(getID()==1 && m_grainHandler->get_loop() == 1){
//			GradientFile.open(Gradient.c_str());
//		}

		max_dt = 0.0;
		for(int k=z_Min+1; k<z_Max-1; k++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int i=x_Min+1; i<x_Max-1; i++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					double ip,im,jp,jm,kp,km,current;
					ip=m_inputDistance->getValueAt(j,i+1,k);
					im=m_inputDistance->getValueAt(j,i-1,k);
					jp=m_inputDistance->getValueAt(j+1,i,k);
					jm=m_inputDistance->getValueAt(j-1,i,k);
					kp=m_inputDistance->getValueAt(j,i,k+1);
					km=m_inputDistance->getValueAt(j,i,k-1);
					current=m_inputDistance->getValueAt(j,i,k);

					if(
							current*ip>0 &&
							current*im>0 &&
							current*jp>0 &&
							current*jm>0 &&
							current*kp>0 &&
							current*km>0
					){

//						double s = current/
//								sqrt(current*
//										current + h*h);

						int s = sign(current);
						sdf_dt->setValueAt(j,i,k,s *
								(1 - norm_grad(m_outputDistance, i, j, k, h)));
						if(
								abs(jp)<m_grainHandler->delta &&
								abs(jm)<m_grainHandler->delta &&
								abs(ip)<m_grainHandler->delta &&
								abs(im)<m_grainHandler->delta &&
								abs(kp)<m_grainHandler->delta &&
								abs(km)<m_grainHandler->delta
						)
							if(abs(sdf_dt->getValueAt(j,i,k))>max_dt) max_dt = abs(sdf_dt->getValueAt(j,i,k));
					}
//					if(getID()==1 && m_grainHandler->get_loop() == 1){
//								GradientFile << sdf_dt->getValueAt(j,i,k) << endl;
//					}
				}
			}
		}

		for(int k=z_Min+1; k<z_Max-1; k++){
			for(int j=y_Min+1; j<y_Max-1; j++){
				for(int i=x_Min+1; i<x_Max-1; i++){
					if(abs(m_outputDistance->getValueAt(j,i,k))>=m_grainHandler->delta)
						continue;
					if(
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i+1,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i-1,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j+1,i,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j-1,i,k)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k+1)>0 &&
							m_inputDistance->getValueAt(j,i,k)*m_inputDistance->getValueAt(j,i,k-1)>0
					){
						m_outputDistance->setValueAt(j,i,k,
								m_outputDistance->getValueAt(j,i,k)+ dt * sdf_dt->getValueAt(j,i,k));
					}
				}
			}
		}

//		if(getID()==1 && m_grainHandler->get_loop() == 1){
//			GradientFile.close();
//		}

		iter++;
	}


	if(getID()==1){
			gettimeofday(&time1, NULL);
			cout << "Time for the correction of the remaining level set function:" << endl;
			cout << time1.tv_sec-time2.tv_sec << ":" << time1.tv_usec-time2.tv_usec << endl;

			cout << "Number of iterations: "<< iter << endl;
	}

	m_outputDistance->clampValues(-m_grainHandler->delta,
			m_grainHandler->delta);

	d->clearValues();
	curvature->clearValues();
	sdf_dt->clearValues();
	delete d;
	delete curvature;
	delete sdf_dt;

}

/**************************************/
//  Redistancing
/**************************************/
void LSbox::executeRedistancing() {
	if (grainExists() != true)
		return;

	double h = m_grainHandler->get_h();
	double candidate, i_slope, distToZero;

	m_outputDistance->clearValues(-1.0);

	// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax,
			intersec_zmin, intersec_zmax;

	if (m_inputDistance->getMinX() < m_outputDistance->getMinX())
		intersec_xmin = m_outputDistance->getMinX();
	else
		intersec_xmin = m_inputDistance->getMinX();

	if (m_inputDistance->getMinY() < m_outputDistance->getMinY())
		intersec_ymin = m_outputDistance->getMinY();
	else
		intersec_ymin = m_inputDistance->getMinY();

	if (m_inputDistance->getMinZ() < m_outputDistance->getMinZ())
		intersec_zmin = m_outputDistance->getMinZ();
	else
		intersec_zmin = m_inputDistance->getMinZ();

	if (m_inputDistance->getMaxX() < m_outputDistance->getMaxX())
		intersec_xmax = m_inputDistance->getMaxX();
	else
		intersec_xmax = m_outputDistance->getMaxX();

	if (m_inputDistance->getMaxY() < m_outputDistance->getMaxY())
		intersec_ymax = m_inputDistance->getMaxY();
	else
		intersec_ymax = m_outputDistance->getMaxY();

	if (m_inputDistance->getMaxZ() < m_outputDistance->getMaxZ())
		intersec_zmax = m_inputDistance->getMaxZ();
	else
		intersec_zmax = m_outputDistance->getMaxZ();

	// first to updates layer by layer to take advantage of the order of point in memory - there are aligned layer by layer.
	for (int k = intersec_zmin; k < intersec_zmax - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX() - 1;
					j++) {
				// x-direction forward
				if (j < intersec_xmax - 1 && i < intersec_ymax) {
					if (m_inputDistance->getValueAt(i, j, k)
							* m_inputDistance->getValueAt(i, j + 1, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j + 1, k)
								- m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i,
													j + 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j + 1, k)))
						m_outputDistance->setValueAt(i, j + 1, k, candidate);
				} else {
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ (sgn(m_outputDistance->getValueAt(i, j + 1, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j + 1, k)))
						m_outputDistance->setValueAt(i, j + 1, k, candidate);
				}
			}
		}

		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
			for (int j = intersec_xmax - 1; j > m_outputDistance->getMinX();
					j--) {
				// x-direction outputDistanceward
				//check for sign change
				if (j > intersec_xmin && i < intersec_ymax) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i,
													j - 1, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j - 1, k)))
						m_outputDistance->setValueAt(i, j - 1, k, candidate);
				} else {
					//<<<<<<< HEAD
					//					candidate = m_outputDistance->getValueAt(i, j, k)
					//							+ sgn(m_outputDistance->getValueAt(i, j - 1, k))
					//									* h;
					//					if (abs(candidate)
					//							< abs(m_outputDistance->getValueAt(i, j - 1, k)))
					//=======
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ sgn(m_outputDistance->getValueAt(i, j - 1, k))
									* h;
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j - 1, k)))
						//>>>>>>> 66c3d6672001fb0db664a7cc036f13ecf8da0d05
						m_outputDistance->setValueAt(i, j - 1, k, candidate);
				}
			}
		}

		// y-direction forward
		for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
			for (int i = intersec_ymin; i < m_outputDistance->getMaxY() - 1;
					i++) {
				if (j < intersec_xmax && i < intersec_ymax - 1) {
					if (m_inputDistance->getValueAt(i, j, k)
							* m_inputDistance->getValueAt(i + 1, j, k) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i + 1, j, k)
								- m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i + 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i + 1, j, k)))
						m_outputDistance->setValueAt(i + 1, j, k, candidate);
				} else {
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ (sgn(m_outputDistance->getValueAt(i + 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i + 1, j, k)))
						m_outputDistance->setValueAt(i + 1, j, k, candidate);
				}
			}
		}

		for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
			for (int i = intersec_ymax - 1; i > m_outputDistance->getMinY();
					i--) {
				if (j < intersec_xmax && i > intersec_ymin) {
					// calculate new distance candidate and assign if appropriate
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i - 1,
													j, k)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i - 1, j, k)))
						m_outputDistance->setValueAt(i - 1, j, k, candidate);
				} else {
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ (sgn(m_outputDistance->getValueAt(i - 1, j, k))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i - 1, j, k)))
						m_outputDistance->setValueAt(i - 1, j, k, candidate);
				}
			}
		}
	}
	// update the point into the third dimensio. the strategy has to change to avoid unneccesary cache loads
	// the idea is to compare all points in one layer first to the next and go on:
	//TODO redist into the third direction
	// z forward:
	for (int k = intersec_zmin; k < m_outputDistance->getMaxZ() - 1; k++) {
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
				// x-direction forward
				if (k < intersec_zmax - 1 && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance->getValueAt(i, j, k)
							* m_inputDistance->getValueAt(i, j, k + 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j, k + 1)
								- m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i, j,
													k + 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j, k + 1)))
						m_outputDistance->setValueAt(i, j, k + 1, candidate);
				} else {
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ (sgn(m_outputDistance->getValueAt(i, j, k + 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j, k + 1)))
						m_outputDistance->setValueAt(i, j, k + 1, candidate);
				}
			}
		}
	}
	// z backward:
	for (int k = intersec_zmax - 1; k > m_outputDistance->getMinZ(); k--) {
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
			for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
				// x-direction forward
				if (k > intersec_zmin && i < intersec_ymax
						&& j < intersec_xmax) {
					if (m_inputDistance->getValueAt(i, j, k)
							* m_inputDistance->getValueAt(i, j, k - 1) <= 0.0) {
						// interpolate
						i_slope = (m_inputDistance->getValueAt(i, j, k - 1)
								- m_inputDistance->getValueAt(i, j, k)) / h;
						distToZero = -m_inputDistance->getValueAt(i, j, k)
								/ i_slope;
						if (abs(m_outputDistance->getValueAt(i, j, k))
								> abs(distToZero))
							m_outputDistance->setValueAt(i, j, k,
									-distToZero * sgn(i_slope));
					}
					candidate =
							m_outputDistance->getValueAt(i, j, k)
									+ (sgn(
											m_inputDistance->getValueAt(i, j,
													k - 1)) * h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j, k - 1)))
						m_outputDistance->setValueAt(i, j, k - 1, candidate);
				} else {
					candidate = m_outputDistance->getValueAt(i, j, k)
							+ (sgn(m_outputDistance->getValueAt(i, j, k - 1))
									* h);
					if (abs(candidate)
							< abs(m_outputDistance->getValueAt(i, j, k - 1)))
						m_outputDistance->setValueAt(i, j, k - 1, candidate);
				}
			}
		}
	}

	// for all layers the redist is done - compare into the depth to do
	m_outputDistance->clampValues(-m_grainHandler->delta,
			m_grainHandler->delta);

	m_inputDistance->resize(m_outputDistance->getMinX(),
			m_outputDistance->getMinY(), m_outputDistance->getMinZ(),
			m_outputDistance->getMaxX(), m_outputDistance->getMaxY(),
			m_outputDistance->getMaxZ());
}

/**************************************/
// end of redist
/**************************************/

/**************************************/
// plot the box and all its properties
/**************************************/
void LSbox::resizeGrid(int newSize, double h) {
	int ngridpointsNew = m_grainHandler->get_ngridpoints();
	double hn = m_grainHandler->get_h();
	int grid_blowup = m_grainHandler->get_grid_blowup();

	m_newXMin = (m_outputDistance->getMinX() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;
	m_newXMax = (m_outputDistance->getMaxX() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;
	m_newYMin = (m_outputDistance->getMinY() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;
	m_newYMax = (m_outputDistance->getMaxY() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;
	m_newZMin = (m_outputDistance->getMinZ() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;
	m_newZMax = (m_outputDistance->getMaxZ() - grid_blowup) * (h / hn) + 0.5
			+ grid_blowup;

	double xl, xr, yo, yu, zu, zo;

	double pointx, pointy, pointz;
	//	m_inputDistance->clearValues(-1.0);
	m_inputDistance->resize(m_newXMin, m_newYMin, m_newZMin, m_newXMax,
			m_newYMax, m_newZMax);
	m_inputDistance->resizeToCube(ngridpointsNew);

	for (int k = m_inputDistance->getMinZ(); k < m_inputDistance->getMaxZ();
			k++) {
		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
				i++) {
			for (int j = m_inputDistance->getMinX();
					j < m_inputDistance->getMaxX(); j++) {
				pointx = (j - grid_blowup) * (hn / h);
				pointy = (i - grid_blowup) * (hn / h);
				pointz = (k - grid_blowup) * (hn / h);

				xl = floor(pointx);
				xr = floor(pointx + 1);
				yu = floor(pointy);
				yo = floor(pointy + 1);
				zu = floor(pointz);
				zo = floor(pointz + 1);

				if (xr + grid_blowup > m_outputDistance->getMaxX() - 2
						|| yo + grid_blowup > m_outputDistance->getMaxY() - 2
						|| zo + grid_blowup > m_outputDistance->getMaxZ() - 2
						|| yu + grid_blowup < m_outputDistance->getMinY()
						|| xl + grid_blowup < m_outputDistance->getMinX()
						|| zu + grid_blowup < m_outputDistance->getMinZ()) {
					m_inputDistance->setValueAt(i, j, k,
							-m_grainHandler->delta);
					continue;
				}

				double ro, ru, newDistVal_zu, newDistVal_zo;
				double trilinearInterpolation;

				//bilinear interpolation in the plane zu:
				ro = (xr - pointx)
						* m_outputDistance->getValueAt(yo + grid_blowup,
								xl + grid_blowup, zu + grid_blowup)
						+ (pointx - xl)
								* m_outputDistance->getValueAt(yo + grid_blowup,
										xr + grid_blowup, zu + grid_blowup);
				ru = (xr - pointx)
						* m_outputDistance->getValueAt(yu + grid_blowup,
								xl + grid_blowup, zu + grid_blowup)
						+ (pointx - xl)
								* m_outputDistance->getValueAt(yu + grid_blowup,
										xr + grid_blowup, zu + grid_blowup);
				newDistVal_zu = (yo - pointy) * ru + (pointy - yu) * ro;

				//bilinear interpolation in the plane zo:
				ro = (xr - pointx)
						* m_outputDistance->getValueAt(yo + grid_blowup,
								xl + grid_blowup, zo + grid_blowup)
						+ (pointx - xl)
								* m_outputDistance->getValueAt(yo + grid_blowup,
										xr + grid_blowup, zo + grid_blowup);
				ru = (xr - pointx)
						* m_outputDistance->getValueAt(yu + grid_blowup,
								xl + grid_blowup, zo + grid_blowup)
						+ (pointx - xl)
								* m_outputDistance->getValueAt(yu + grid_blowup,
										xr + grid_blowup, zo + grid_blowup);
				newDistVal_zo = (yo - pointy) * ru + (pointy - yu) * ro;

				trilinearInterpolation = (zo - pointz) * newDistVal_zu
						+ (pointz - zu) * newDistVal_zo;

				m_inputDistance->setValueAt(i, j, k, trilinearInterpolation);

			}
		}
	}
	//		m_outputDistance->clearValues(-1.0);
	m_outputDistance->resize(m_inputDistance->getMinX(),
			m_inputDistance->getMinY(), m_inputDistance->getMinZ(),
			m_inputDistance->getMaxX(), m_inputDistance->getMaxY(),
			m_inputDistance->getMaxZ());
	//m_outputDistance->resizeToCube(ngridpointsNew);

}

void LSbox::recalculateIDLocal() {
	resizeIDLocalToDistanceBuffer();
	executeComparison();
}

double LSbox::computeMisorientation(LSbox* grain_2) {
	double result = 0.0;

	//	if (Settings::LatticeType == E_CUBIC) {
	result = m_orientationQuat->misorientationCubicQxQ(
			grain_2->m_orientationQuat);
	//	} else if (Settings::LatticeType == E_HEXAGONAL) {
	//		result
	//				= m_grainHandler->m_misOriHdl->calculateMisorientation_hexagonal(
	//						m_orientationQuat, grain_2->m_orientationQuat);
	//	}
	if (result > 1 * PI / 180.0)
		return result;
	else
		return 1 * PI / 180.0;
}
double LSbox::computeMisorientation(unsigned int grainID) {
	return computeMisorientation(m_grainHandler->getGrainByID(grainID));
}

double LSbox::GBmobilityModel(double thetaMis) {
	return 1.0;
}

bool LSbox::isNeighbour(LSbox* candidate) {
	return true;
}

double LSbox::getWeight(int i, int j, bool minimal) {
	return -1;
}

void LSbox::calculateCentroid(SPoint& centroid,
		vector<GrainJunction> junctions) {

	int nVertices = junctions.size();
	SPoint cent;
	cent.x = 0.0;
	cent.y = 0.0;
	double areaSigned = 0.0;
	double x0 = 0.0; // First vertex' x value
	double y0 = 0.0; // Second vertex' Y value
	double x1 = 0.0; // Next vertex' x value
	double y1 = 0.0; // Next vertex' y value
	double loopArea = 0.0; // Intermediate area

	int i, j = 0;
	for (i = 0, j = nVertices - 1; i < nVertices; j = i++) {

		x0 = junctions[i].coordinates.x;
		y0 = junctions[i].coordinates.y;
		x1 = junctions[j].coordinates.x;
		y1 = junctions[j].coordinates.y;
		loopArea = x0 * y1 - x1 * y0;
		areaSigned += loopArea;
		cent.x += (x0 + x1) * loopArea;
		cent.y += (y0 + y1) * loopArea;
	}

	areaSigned *= 0.5;
	cent.x /= (6.0 * areaSigned);
	cent.y /= (6.0 * areaSigned);

	centroid = cent;
}

double LSbox::getWeigthFromHandler(int i, int j) {
	return m_grainHandler->weightsMatrix[i][j];
}

vector<double> LSbox::linearRegression(vector<SPoint>& points2D) {

	//! the resulting vector contains two elements. The first element is
	//! the slope of the regression line and the second element is the
	//! y-intercept of this line.

	vector<double> linearRegression;
	linearRegression.reserve(2);

	int numberPoints = points2D.size();

	double meanX = 0.0;
	double meanY = 0.0;

	vector<SPoint>::iterator iter;
	for (iter = points2D.begin(); iter != points2D.end(); iter++) {

		meanX += (*iter).x;
		meanY += (*iter).y;
	}

	meanX /= double(numberPoints);
	meanY /= double(numberPoints);

	double covarianceXY = 0.0;
	double varianceX = 0.0;
	for (int i = 0; i < numberPoints; i++) {
		covarianceXY += (points2D[i].x - meanX) * (points2D[i].y - meanY);
		varianceX += (points2D[i].x - meanX) * (points2D[i].x - meanX);
	}

	linearRegression.push_back(covarianceXY / varianceX);
	linearRegression.push_back(meanY - ((covarianceXY / varianceX) * meanX));

	return linearRegression;
}

void LSbox::calculateTriangleCentroid(vector<SPoint>& triangleCentroid,
		vector<SPoint> triangle) {

	int nVertices = triangle.size();
	SPoint cent;
	cent.x = 0.0;
	cent.y = 0.0;
	double areaSigned = 0.0;
	double x0 = 0.0; // First vertex' x value
	double y0 = 0.0; // Second vertex' Y value
	double x1 = 0.0; // Next vertex' x value
	double y1 = 0.0; // Next vertex' y value
	double loopArea = 0.0; // Intermediate area

	int i, j = 0;
	for (i = 0, j = nVertices - 1; i < nVertices; j = i++) {

		x0 = triangle[i].x;
		y0 = triangle[i].y;
		x1 = triangle[j].x;
		y1 = triangle[j].y;
		//cout << "Junction type of i is " << junctions[i].junction_type << " Junction type of j is " << junctions[j].junction_type <<endl;
		loopArea = x0 * y1 - x1 * y0;
		areaSigned += loopArea;
		cent.x += (x0 + x1) * loopArea;
		cent.y += (y0 + y1) * loopArea;
	}

	areaSigned *= 0.5;
	cent.x /= (6.0 * areaSigned);
	cent.y /= (6.0 * areaSigned);

	triangleCentroid.push_back(cent);
}

void LSbox::outputMemoryUsage(ofstream& output) {
	output << m_inputDistance->getTotalMemoryUsed() << " "
			<< m_outputDistance->getTotalMemoryUsed() << " "
			<< m_IDLocal.getTotalMemoryUsed() << endl;

	output << m_outputDistance->getMaxX() - m_outputDistance->getMinX() << " "
			<< m_outputDistance->getMaxY() - m_outputDistance->getMinY() << " "
			<< m_outputDistance->getMaxZ() - m_outputDistance->getMinZ()
			<< endl;
	output << m_inputDistance->getMaxX() - m_inputDistance->getMinX() << " "
			<< m_inputDistance->getMaxY() - m_inputDistance->getMinY() << " "
			<< m_inputDistance->getMaxZ() - m_inputDistance->getMinZ() << endl;
	output << m_IDLocal.getMaxX() - m_IDLocal.getMinX() << " "
			<< m_IDLocal.getMaxY() - m_IDLocal.getMinY() << " "
			<< m_IDLocal.getMaxZ() - m_IDLocal.getMinZ() << endl;
}

vector<int> LSbox::getDirectNeighbourIDs() {
	return vector<int>();
}

vector<double> LSbox::getGBLengths() {
	return vector<double>();
}

void LSbox::computeDirectNeighbours(
		const RTree<unsigned int, int, 3, float>& tree) {
	int min[3], max[3];
	min[0] = getMinX();
	min[1] = getMinY();
	min[2] = getMinZ();
	max[0] = getMaxX();
	max[1] = getMaxY();
	max[2] = getMaxZ();
	vector<unsigned int> intersectingGrains;
	tree.Search(min, max, intersectingGrains);
	for (unsigned int k = 0; k < intersectingGrains.size(); k++) {
		if (m_ID != intersectingGrains[k]) {
			m_secondOrderNeighbours.push_back(intersectingGrains[k]);
		}
	}
}

void LSbox::plotBoxInterfacialElements(bool absoluteCoordinates) {
	m_explicitHull.plotInterfacialElements(false, m_grainHandler->loop);
}

void LSbox::plotBoxContour(bool absoluteCoordinates) {
	m_explicitHull.plotContour(false, m_grainHandler->loop);
}

void LSbox::plotNeighboringGrains(bool absoluteCoordinates) {
	m_explicitHull.plotNeighboringGrains(false, m_grainHandler->loop,Settings::NeighbourhoodOrder);
}

void LSbox::plotBoxVolumetric(string identifier,
		E_BUFFER_SELECTION bufferSelection, double h) {
	string filename = string("GrainVolume_")
			+ to_string((unsigned long long) m_ID) + string("Timestep_")
			+ to_string((unsigned long long) m_grainHandler->loop) + identifier
			+ string(".vtk");

	FILE* output = fopen(filename.c_str(), "wt");
	if (output == NULL) {
		throw runtime_error("Unable to save box hull!");
	}
	DimensionalBufferReal* distance;
	switch (bufferSelection) {
	case E_INPUT_DISTANCE:
		distance = m_inputDistance;
		break;
	case E_OUTPUT_DISTANCE:
		distance = m_outputDistance;
		break;
	default:
		throw runtime_error(
				string(
						"Invalid buffer selection in plotBoxVolumetric for grain")
						+ to_string((unsigned long long) m_ID)
						+ string(" at timestep ")
						+ to_string((unsigned long long) m_grainHandler->loop));
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
			"vtk output\n"
			"ASCII\n"
			"DATASET STRUCTURED_GRID");
	int totalPoints = (distance->getMaxX() - distance->getMinX())
			* (distance->getMaxY() - distance->getMinY())
			* (distance->getMaxZ() - distance->getMinZ());
	fprintf(output, "DIMENSIONS %d %d %d\n",
			distance->getMaxX() - distance->getMinX(),
			distance->getMaxY() - distance->getMinY(),
			distance->getMaxZ() - distance->getMinZ());
	fprintf(output, "POINTS %d float\n", totalPoints);
	double ii, jj, kk;
	for (int i = distance->getMinY(); i < distance->getMaxY(); i++)
		for (int j = distance->getMinX(); j < distance->getMaxX(); j++)
			for (int k = distance->getMinZ(); k < distance->getMaxZ(); k++) {
				//ii = ((double)(i-m_grainHandler->get_grid_blowup())-(double)(distance->getMaxY()+distance->getMinY())/2.)*h;
				//jj = ((double)(j-m_grainHandler->get_grid_blowup())-(double)(distance->getMaxX()+distance->getMinX())/2.)*h;
				//kk = ((double)(k-m_grainHandler->get_grid_blowup())-(double)(distance->getMaxZ()+distance->getMinZ())/2.)*h;
				ii = (i - m_grainHandler->get_grid_blowup()) * h;
				jj = (j - m_grainHandler->get_grid_blowup()) * h;
				kk = (k - m_grainHandler->get_grid_blowup()) * h;
				fprintf(output, "%f %f %f\n", jj, ii, kk);
			}

	fprintf(output, "\nPOINT_DATA %d\n", totalPoints);
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Distance 1 %d float\n", totalPoints);
	for (int i = distance->getMinY(); i < distance->getMaxY(); i++)
		for (int j = distance->getMinX(); j < distance->getMaxX(); j++)
			for (int k = distance->getMinZ(); k < distance->getMaxZ(); k++)
				fprintf(output, "%f\n", (distance->getValueAt(i, j, k)));
	fclose(output);
}

void LSbox::plotBoxIDLocal() {
	string filename = string("GrainIDLocal_")
			+ to_string((unsigned long long) m_ID) + string("Timestep_")
			+ to_string((unsigned long long) m_grainHandler->loop)
			+ string(".vtk");
	FILE* output = fopen(filename.c_str(), "wt");
	if (output == NULL) {
		throw runtime_error("Unable to save box hull!");
	}

	fprintf(output, "%s\n", "# vtk DataFile Version 3.0\n"
			"vtk output\n"
			"ASCII\n"
			"DATASET STRUCTURED_GRID");
	int totalPoints = (m_IDLocal.getMaxX() - m_IDLocal.getMinX())
			* (m_IDLocal.getMaxY() - m_IDLocal.getMinY())
			* (m_IDLocal.getMaxZ() - m_IDLocal.getMinZ());
	fprintf(output, "DIMENSIONS %d %d %d\n",
			m_IDLocal.getMaxX() - m_IDLocal.getMinX(),
			m_IDLocal.getMaxY() - m_IDLocal.getMinY(),
			m_IDLocal.getMaxZ() - m_IDLocal.getMinZ());
	fprintf(output, "POINTS %d int\n", totalPoints);
	for (int i = m_IDLocal.getMinY(); i < m_IDLocal.getMaxY(); i++)
		for (int j = m_IDLocal.getMinX(); j < m_IDLocal.getMaxX(); j++)
			for (int k = m_IDLocal.getMinZ(); k < m_IDLocal.getMaxZ(); k++)
				fprintf(output, "%d %d %d\n", j, i, k);

	fprintf(output, "\nPOINT_DATA %d\n", totalPoints);
	fprintf(output, "FIELD FieldData 1\n");
	fprintf(output, "Distance 1 %d float\n", totalPoints);
	for (int i = m_IDLocal.getMinY(); i < m_IDLocal.getMaxY(); i++)
		for (int j = m_IDLocal.getMinX(); j < m_IDLocal.getMaxX(); j++)
			for (int k = m_IDLocal.getMinZ(); k < m_IDLocal.getMaxZ(); k++)
				fprintf(output, "%d\n", m_IDLocal.getValueAt(i, j, k).grainID);
	fclose(output);
}
int LSbox::getNeighbourAt(int i, int j, int k) {
	if (m_IDLocal.isPointInside(i, j, k)) {
		return m_IDLocal.getValueAt(i, j, k).grainID;
	} else {
		return -1;
	}
}

void LSbox::copyDataToContainer(DimensionalBuffer<unsigned int> * container,
		int threadID) {
	int gridblowup = m_grainHandler->get_grid_blowup();
	for (int k = getMinZ(); k < getMaxZ(); k++) {
		for (int i = getMinY(); i < getMaxY(); i++) {
			for (int j = getMinX(); j < getMaxX(); j++) {
				if (0. <= m_outputDistance->getValueAt(i, j, k)) {
					if (i >= gridblowup && j >= gridblowup && k >= gridblowup
							&& i < container->getMaxY() + gridblowup
							&& j < container->getMaxX() + gridblowup
							&& k < container->getMaxZ() + gridblowup)
						if (Settings::StoreTaskDistribution)
							container->setValueAt(i - gridblowup,
									j - gridblowup, k - gridblowup,
									(unsigned int) threadID);
						else
							container->setValueAt(i - gridblowup,
									j - gridblowup, k - gridblowup,
									(unsigned int) getID());

				}
			}
		}
	}
}

TextureData LSbox::collectTextureData() {
	double x = (getMaxX() - getMinX()) / 2 + getMinX();
	x *= m_grainHandler->get_h();
	double y = (getMaxY() - getMinY()) / 2 + getMinY();
	y *= m_grainHandler->get_h();
	double z = (getMaxZ() - getMinZ()) / 2 + getMinZ();
	z *= m_grainHandler->get_h();
	TextureData newdata(m_ID, m_explicitHull.getAllNeighborsCount(),
			m_explicitHull.IsNeighbor(0), m_volume, m_surface, m_energy,
			m_StoredElasticEnergy, m_orientationQuat, x, y, z);
	return newdata;
}
