class myQuaternion;

struct VolEvolution {
	double dA;
	int nVertex;

	VolEvolution(double da, int nr) :
			dA(da), nVertex(nr) {
	}
};

struct Face {
	double area;
	unsigned int grainA;
	unsigned int grainB;
	Face(double _area, unsigned int _grainA, unsigned int _grainB) :
			area(_area), grainA(_grainA), grainB(_grainB) {
	}
};

struct TextureData {
	double volume;
	double surfaceArea;
	double GBEnergy;
	double BulkEnergy;
	double phi1;
	double PHI;
	double phi2;
	double x;
	double y;
	double z;
	unsigned int id;
	unsigned int NeighbourCount;
	unsigned int intersectsBoundaryGrain;
	TextureData(unsigned int _id, unsigned int _NeighbourCount,
			bool _intersectsBoundaryGrain, double _volume, double _surfaceArea,
			double _GBEnergy, double _BulkEnergy, myQuaternion *ori, double _x,
			double _y, double _z) :
			id(_id), NeighbourCount(_NeighbourCount), intersectsBoundaryGrain(
					_intersectsBoundaryGrain), volume(_volume), surfaceArea(
					_surfaceArea), GBEnergy(_GBEnergy), BulkEnergy(_BulkEnergy),
					x(_x), y(_y), z(_z) {
		double* newori = ori->Quaternion2Euler();
		phi1 = newori[0];
		PHI = newori[1];
		phi2 = newori[2];
		delete[] newori;

	}
	~TextureData() {
	}
};
