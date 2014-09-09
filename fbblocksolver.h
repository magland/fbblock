#ifndef fbblocksolver_H
#define fbblocksolver_H

#include "fbglobal.h"
#include "arrays.h"
#include "fberrorestimator.h"
#include "fbblock.h"

class FBBlockSolverPrivate;
class FBBlockSolver {
public:
	friend class FBBlockSolverPrivate;
	FBBlockSolver();
	virtual ~FBBlockSolver();
	void setEpsilon(fbreal epsilon);
	void setMaxIterations(int val);
	void setNumThreads(int val);
	void setUsePreconditioner(bool val);
	void setStiffnessMatrix(const FBArray2D<float> &stiffness_matrix);
	void setYoungsModulus(float val);
	void setVoxelVolume(float val);
	void setBVFMap(const FBArray3D<unsigned char> &bvf_map); // N1 x N2 x N3
	void setInitialDisplacementsOnFreeVariables(const FBSparseArray4D &displacements); //3x(N1+1)x(N2+1)x(N3+1)
	void setInitialDisplacementsOnFreeVariables(const FBArray4D<float> &displacements); //3x(N1+1)x(N2+1)x(N3+1)
	void setInitialDisplacements(FBMacroscopicStrain &strain);
	void setFixedVariables(const FBSparseArray4D &fixed_variables); //3x(N1+1)x(N2+1)x(N3+1)
	long setFixedVariables(FBMacroscopicStrain &macroscopic_strain);
	void setResolution(QList<fbreal> &res);
	void solve(); 
	void clear();
	
	void solveNonlinear(float step_size,int num_steps,int num_iterations_per_step);
	
	int getNumIterations();
	QList<double> getStress();	
	void getDisplacements(FBSparseArray4D &displacements); //3x(N1+1)x(N2+1)x(N3+1)
	void getDisplacements(FBArray4D<float> &displacements); //3x(N1+1)x(N2+1)x(N3+1)
	void getForces(FBSparseArray4D &forces); //3x(N1+1)x(N2+1)x(N3+1)
	void getEnergy(FBSparseArray4D &energy); //1 x N1 x N2 x N3
	FBErrorEstimator *errorEstimator();
private:
	FBBlockSolverPrivate *d;
};

#endif
