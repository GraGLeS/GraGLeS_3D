#ifndef		__TRINAGLE_H__
#define 	__TRINAGLE_H__

#include "Eigen/Dense"
#include <iostream>
using namespace Eigen;
using namespace std;

struct Triangle {
	Vector3d points[3];
	int additionalData;
	Vector3d computeBarycenter() {
		Vector3d barycenter;
		barycenter[0] = (points[0][0] + points[1][0] + points[2][0]) / 3;
		barycenter[1] = (points[0][1] + points[1][1] + points[2][1]) / 3;
		barycenter[2] = (points[0][2] + points[1][2] + points[2][2]) / 3;
		return barycenter;
	}
	Vector3d computeOuterUnitNormal() {
		Vector3d outerUnitNormal;
		Vector3d R;
		Vector3d S;
		R[1] = points[0][0] - points[1][0];
		R[1] = points[0][1] - points[1][1];
		R[3] = points[0][2] - points[1][2];
		S[1] = points[0][0] - points[2][0];
		S[1] = points[0][1] - points[2][1];
		S[3] = points[0][2] - points[2][2];
		outerUnitNormal[0] = R[2] * S[3] - R[3] * S[2];
		outerUnitNormal[1] = R[3] * S[1] - R[1] * S[3];
		outerUnitNormal[2] = R[1] * S[2] - R[2] * S[1];
		return outerUnitNormal;
	}
	bool IsNeighbor(Triangle &T2) {
		int NNeighbors = 0;
		int t1[2];
		int t2[2];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (points[i] == T2.points[j]) {
					t1[NNeighbors] = i;
					t2[NNeighbors] = j;
					NNeighbors++;
				}
			}
		}
		if (NNeighbors == 2) {

			/*
			 * Check if the triangles have the correct orientation
			 */

			if ((t1[1] - t1[0] + 3) % 3 - (t2[1] - t2[0] + 3) % 3 == 0) {
				Vector3d temp;
				temp = T2.points[t2[1]];
				T2.points[t2[1]] = T2.points[t2[0]];
				T2.points[t2[0]] = temp;
			}
			return true;
		} else
			return false;
	}

	double calculateMeanWidthComponent(Triangle &T2, Vector3d normal1,
			Vector3d normal2) {
		double scalarProduct;
		double beta, epsilon;
		if (normal1.norm() == 0 || normal2.norm() == 0)
			return 0.0;
		normal1.normalize();
		normal2.normalize();
		scalarProduct = normal1.dot(normal2);

		/*
		 * It is possible that through an numerical error the scalar product of two parallel vectors
		 * becomes greater than one. In this case the value of the scalar product is set to one.
		 */

		if (scalarProduct > 1) {
			scalarProduct = 1.;
		}
		if (scalarProduct < -1) {
			scalarProduct = -1.;
		}
		beta = acos(scalarProduct);

		/*
		 * Find the points that are shared by both triangles
		 */

		int sharedPoints[2];
		int temp = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (points[i] == T2.points[j]) {
					sharedPoints[temp] = i;
					temp++;
				}
			}
		}

		/*
		 * Orient the shared line as it is oriented in triangle T1
		 */
		Vector3d sharedLine;
		if ((sharedPoints[1] - sharedPoints[0] + 3) % 3 == 1) {
			sharedLine = points[sharedPoints[1]] - points[sharedPoints[0]];
		} else {
			sharedLine = points[sharedPoints[0]] - points[sharedPoints[1]];
		}

		epsilon = sharedLine.norm();

		/*
		 * If the cross product of the two normal vectors is parallel to the shared line the angle is positive
		 * if it is antiparallel the angle is negative this is accounted for in the next if-statement
		 */
		Vector3d crossprod = normal1.cross(normal2);
		if (crossprod.dot(sharedLine)> 0) {
			beta *= -1;
		}
		if (beta != beta) {
			cout << "beta" << endl;
			cout << normal1.transpose() * normal2 << endl;
			cout << normal1.transpose() << endl;
			cout << normal2.transpose() << endl;
		}
		if (epsilon != epsilon)
			cout << "epsilon" << endl;
		return beta * epsilon;

	}
};

#endif		//__TRIANGLE_H__
