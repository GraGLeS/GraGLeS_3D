#ifndef		__IGRAIN_SCHEDULER__
#define		__IGRAIN_SCHEDULER__

#include <vector>
#include "Eigen/Dense"

using namespace std;


class IGrainScheduler
{
public:
	IGrainScheduler(){}
	virtual ~IGrainScheduler(){}
	virtual void buildGrainWorkloads(vector<vector<Eigen::Vector3d>*>& hulls, int n_gridpoints) = 0;
	virtual std::vector<unsigned int>&	getThreadWorkload(int threadID) = 0;
};

#endif		//__IGRAIN_SCHEDULER__
