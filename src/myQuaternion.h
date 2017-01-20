/*
 * myQuaternion.h
 *
 *  Created on: 17.09.2015
 *      Author: miessen
 */

#include "mymath.h"
#include "ggLS.h"

#define SWAP(A, B) {struct tempStruct { char C[sizeof(A)];} swap_tmp;\
    swap_tmp = *( struct tempStruct*) &A;\
    *( struct tempStruct*) &A = *( struct tempStruct*) &B;\
    *( struct tempStruct*) &B = swap_tmp;}



#ifndef myQuaternion_H_
#define myQuaternion_H_

class mathMethods;

//myQuaternionen class
class myQuaternion{
private:
	double q0;
	double q1;
	double q2;
	double q3;
public:
	myQuaternion(void);
	myQuaternion(double q0, double q1, double q2, double q3);
	myQuaternion(double alpha,double x,double y,double z, bool radiant);
	~myQuaternion(void);
	myQuaternion& operator= ( const myQuaternion &rhs);
	//myQuaternion& myQuaternion::operator+ (Quarternion const& lhs, Quarternion const& rhs);
	//myQuaternion& myQuaternion::operator- (Quarternion const& lhs, Quarternion const& rhs);
	myQuaternion operator* (const myQuaternion &rhs);

	bool operator == (myQuaternion const& rhs);
	bool operator != (myQuaternion const& rhs);

	myQuaternion Inverse();
	double getNorm(void);
	void Normalize();
	void Sort(void);
	myQuaternion misorientationQuaternionCubic( myQuaternion* p );
	double rotationAngleBetweenQuaternions( myQuaternion *p );
	void euler2Quaternion(double *euler);
	double* Quaternion2EulerConst(void) const;
	double* Quaternion2Euler(void);
	void randomOriShoemakeQuat(mathMethods* math);
	double misorientationCubicQxQ(myQuaternion* p);

	inline double get_q0() {return q0;};
	inline double get_q1() {return q1;};
	inline double get_q2() {return q2;};
	inline double get_q3() {return q3;};
};


#endif /* myQuaternion_H_ */
