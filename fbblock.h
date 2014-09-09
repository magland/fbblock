#ifndef fbblock_H
#define fbblock_H

#include "arrays.h"
#include "nonlinearadjuster.h"

/*

*	*	*	*	*	*	*
*	#	#	#	#	#	*
*	#	&	&	&	#	*
*	#	&	&	&	#	*
*	#	&	&	&	#	*
*	#	#	#	#	#	*
*	*	*	*	*	*	*

In this case Nx=Ny=5
* = outer-interface vertices: (Nx+2)(Ny+2)(Nz+2) - NxNyNz    (inner-interfaces of neighbor block)
# = inner-interface vertices: NxNyNz - (Nx-2)(Ny-2)(Nz-2) (outer-interfaces of neighbor blocks)
& = internal vertices: (Nx-2)(Ny-2)(Nz-2)
& and # are the owned vertices: NxNyNz

number of internal cells: (Nx-1)(Ny-1)(Nz-1)
number of external cells: (Nx+1)(Ny+1)(Nz+1) - (Nx-1)(Ny-1)(Nz-1)
Total number of cells: (Nx+1)(Ny+1)(Nz+1)

*/

struct FBBlockSetupParameters {
	//input
	int Nx,Ny,Nz; //this block owns all vertices within a Nx x Ny x Nz grid
	FBArray3D<unsigned char> BVF; //(Nx+1)x(Ny+1)x(Nz+1) array of bvf values ranging between 0 and 100 (including outer interface)
	FBArray4D<float> X0; //(Nx+2)x(Ny+2)x(Nz+2)x3 array of initial displacements (including outer interface)
	FBArray4D<unsigned char> fixed; //(Nx+2)x(Ny+2)x(Nz+2)x3 array of fixed variables (boundary conditions) (including outer interface)
	FBArray2D<float> stiffness_matrix; //24x24
	float youngs_modulus;
	float voxel_volume;
	bool use_preconditioner;
	float resolution[3];
	int block_x_position;
	int block_y_position;
	int block_z_position;
	
	//output
	//FBSparseArray4D p_on_inner_interface; // (Nx+2)x(Ny+2)x(Nz+2)x3, with values only on the free variables of the inner interface
	FBArray3D<float> p_on_top_inner_interface; // 3x(Nx)x(Ny), with values only on the free variables of the top inner interface
	FBArray3D<float> p_on_bottom_inner_interface; // 3x(Nx)x(Ny), with values only on the free variables of the bottom inner interface
	double bnorm2;
	double rnorm2;
};

struct FBBlockIterateStepAParameters {
	//input
	FBArray3D<float> p_on_top_outer_interface; //3x(Nx)x(Ny), with values only on the free variables of the top outer interface
	FBArray3D<float> p_on_bottom_outer_interface; //3x(Nx)x(Ny), with values only on the free variables of the bottom outer interface
	//FBSparseArray4D p_on_outer_interface; // (Nx+2)x(Ny+2)x(Nz+2)x3, with values only on the free variables of the outer interface
	
	//output
	//the following inner products are computed on the free variables of the "owned" vertices
	double r_Ap;
	double p_Ap;
	double Ap_Ap;
	double r_z;
};

struct FBBlockIterateStepBParameters {
	//input
	double alpha;
	double beta;	
	int WN[3];
	
	//output
	FBArray3D<float> p_on_top_inner_interface; //3x(Nx)x(Ny), with values only on the free variables of the top inner interface
	FBArray3D<float> p_on_bottom_inner_interface; //3x(Nx)x(Ny), with values only on the free variables of the bottom inner interface
	//FBSparseArray4D p_on_inner_interface; // (Nx+2)x(Ny+2)x(Nz+2)x3, with values only on the free variables of the inner interface
	double r_r;
	double bb_bb;
	QList<double> stress;
		
	//Notes:
	//  In this step, the FBBlock will internally update r, p, x, and Ap based on alpha, beta, and p_on_outer_interface
};

class FBBlockPrivate;
class FBBlock {
public:
	friend class FBBlockPrivate;
	FBBlock(int block_num);
	virtual ~FBBlock();
	void setup(FBBlockSetupParameters &P);
	void iterate_step_A(FBBlockIterateStepAParameters &P);
	void iterate_step_B(FBBlockIterateStepBParameters &P);
	void setResolution(QList<float> &res);
	void setNonlinearAdjuster(NonlinearAdjuster *X);
	void setBvfThreshold(float thresh);
	
	int Nx() const;
	int Ny() const;
	int Nz() const;
	int xPosition() const;
	int yPosition() const;
	int zPosition() const;
	
	QList<double> getStress();
	float getDisplacement(int x,int y,int z,int d);
	float getForce(int x,int y,int z,int d);
	void computeEnergyMap(FBSparseArray4D &E);
	long variableCount();
	long ownedVariableCount();
	long ownedFreeVariableCount();
	void clearArrays(); //clears all arrays, except for displacements, residuals and variable indices
	void clearArrays2(); //clears displacements, residuals and variable indices
	
private:
	FBBlockPrivate *d;
};


#endif
