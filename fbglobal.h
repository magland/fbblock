#ifndef fbglobal_H
#define fbglobal_H

#include "textfile.h"
#include "arrays.h"

typedef float fbreal;

struct FBMacroscopicStrain {
	//components of the macroscopic strain tensor
	fbreal eps11,eps22,eps33;
	fbreal eps12,eps13,eps23;
	
	//First index is the face, second index is the direction
	//For example, boundaryRestrictions[1][0]=true means that the Y faces are restricted in the X direction.
	bool boundaryRestrictions[3][3]; 
};

void compute_stiffness_kernel(FBArray2D<float> &stiffness_matrix,const fbreal youngs_modulus,const fbreal poissons_ratio,const QList<fbreal> &resolution);
fbreal initial_displacement(int i1,int i2,int i3,int dir,fbreal resolution[3],FBMacroscopicStrain &macroscopic_strain);
bool is_vertex(FBArray3D<unsigned char> &bvfmap,long i1,long i2,long i3);
bool is_element(FBArray3D<unsigned char> &bvfmap,long i1,long i2,long i3);


#endif
