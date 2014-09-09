#include "fbblock.h"
#include <QList>
#include <QDebug>
#include <QTime>
#include "fbglobal.h"
#include "fbtimer.h"
#include <math.h>
#include "mda_io.h"

struct FBBlockElement {
	long ref_indices[4];
	unsigned char bvf;
	float strain;
};

struct FBVertexLocation {
	int x,y,z;
	long ref_index; 
};

class FBBlockPrivate {
public:
	FBBlock *q;
	FBArray2D<float> m_stiffness_matrix; //24x24
	float m_youngs_modulus; //only for reference when computing the element strains for nonlinear analysis
	float m_voxel_volume;
	int m_Nx,m_Ny,m_Nz;
	long m_num_variables;
	FBArray1D<float> m_x;
	FBArray1D<float> m_r;
	FBArray1D<float> m_p;
	FBArray1D<float> m_Ap;
	FBArray1D<unsigned char> m_free;
	FBArray1D<unsigned char> m_vertex_type; //1 = internal, 2 = inner-interface, 3=outer-interface
	FBArray1D<float> m_preconditioner;
	bool m_use_precondioner;
	QVector<FBBlockElement> m_elements;
	QVector<FBVertexLocation> m_outer_vertex_locations; 
	QVector<FBVertexLocation> m_inner_vertex_locations;
	FBArray3D<long> m_variable_indices;
	FBArray3D<unsigned char> m_bvf_map;
	float m_resolution[3];
	int m_block_x_position;
	int m_block_y_position;
	int m_block_z_position;
	QString m_block_id;
	NonlinearAdjuster *m_nonlinear_adjuster;
	
	double inner_product_on_owned_free_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2);
	double inner_product_on_owned_free_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2,const FBArray1D<float> &V3);
	double inner_product_on_owned_fixed_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2);
	void multiply_by_A(FBArray1D<float> &Y,const FBArray1D<float> &X); //Y=AX
	void compute_preconditioner(const FBArray1D<float> &X,FBArray1D<float> &C);
	QList<double> compute_stress();
};

FBBlock::FBBlock(int block_num) 
{
	d=new FBBlockPrivate;
	d->q=this;
	d->m_Nx=d->m_Ny=d->m_Nz=0;
	d->m_num_variables=0;
	d->m_use_precondioner=false;
	d->m_nonlinear_adjuster=0;
	d->m_youngs_modulus=1;
	d->m_voxel_volume=1;
	d->m_block_id=QString("block%1").arg(block_num);
	block_num++;
	for (int i=0; i<3; i++) d->m_resolution[i]=1;
}

FBBlock::~FBBlock()
{
	delete d;
}

void check_for_nan(QString str,FBArray1D<float> &X) {
	bool found=false;
	for (int j=0; j<X.length(); j++) {
		if (isnan(X.ptr[j])) {
			found=true;
		}
	}
	if (found) {
		qDebug()  << "nan found" << str;
	}
}

void FBBlock::setup(FBBlockSetupParameters &P) {
	d->m_bvf_map=P.BVF;
	
	//set the stiffness_matrix and Nx,Ny,Nz
	d->m_stiffness_matrix=P.stiffness_matrix;
	d->m_youngs_modulus=P.youngs_modulus;
	d->m_voxel_volume=P.voxel_volume;
	d->m_Nx=P.Nx;
	d->m_Ny=P.Ny;
	d->m_Nz=P.Nz;
	d->m_use_precondioner=P.use_preconditioner;
	for (int i=0; i<3; i++) d->m_resolution[i]=P.resolution[i];
	d->m_block_x_position=P.block_x_position;
	d->m_block_y_position=P.block_y_position;
	d->m_block_z_position=P.block_z_position;
	
	//determine which vertices are needed
	FBArray3D<unsigned char> vertex_occupancy; 
	vertex_occupancy.allocate(P.Nx+2,P.Ny+2,P.Nz+2);
	for (int zz=0; zz<P.Nz+1; zz++)
	for (int yy=0; yy<P.Ny+1; yy++)
	for (int xx=0; xx<P.Nx+1; xx++) {
		if (P.BVF.value(xx,yy,zz)) {
			for (int dzz=0; dzz<=1; dzz++)
			for (int dyy=0; dyy<=1; dyy++)
			for (int dxx=0; dxx<=1; dxx++) {
				vertex_occupancy.setValue(1,xx+dxx,yy+dyy,zz+dzz);
			}
		}
	}

	//assign the variable indices
	d->m_num_variables=0;
	d->m_variable_indices.allocate(P.Nx+2,P.Ny+2,P.Nz+2);
	d->m_variable_indices.setAll(-1);
	for (int zz=0; zz<P.Nz+2; zz++)
	for (int yy=0; yy<P.Ny+2; yy++)
	for (int xx=0; xx<P.Nx+2; xx++) {
		if (vertex_occupancy.value(xx,yy,zz)) {
			d->m_variable_indices.setValue(d->m_num_variables,xx,yy,zz);
			d->m_num_variables+=3; //each vertex contains 3 directions
		}
	}
	
	//allocate the vectors, and define m_free and m_vertex_type vectors, and initialize m_x.
	//Also, set up m_outer_vertex_locations and m_inner_vertex_locations
	
	if (!d->m_num_variables) {
		qWarning() << "block is empty.";
		return; //when there is no bvf on the whole block
	}
	
	d->m_x.allocate(d->m_num_variables);
	d->m_r.allocate(d->m_num_variables);
	d->m_p.allocate(d->m_num_variables);
	d->m_Ap.allocate(d->m_num_variables);
	d->m_free.allocate(d->m_num_variables);
	d->m_vertex_type.allocate(d->m_num_variables);	
	for (int zz=0; zz<P.Nz+2; zz++)
	for (int yy=0; yy<P.Ny+2; yy++)
	for (int xx=0; xx<P.Nx+2; xx++) {
		if (vertex_occupancy.value(xx,yy,zz)) {
			for (int dd=0; dd<3; dd++) {
				long varind=d->m_variable_indices.value(xx,yy,zz)+dd;
				if (varind>=0) {
					if (!P.fixed.value(xx,yy,zz,dd)) d->m_free.ptr[varind]=1;
					if ((xx>=2)&&(xx<=P.Nx-1)&&(yy>=2)&&(yy<=P.Ny-1)&&(zz>=2)&&(zz<=P.Nz-1)) {
						d->m_vertex_type.ptr[varind]=1; //internal vertex
					}
					else if ((xx>=1)&&(xx<=P.Nx)&&(yy>=1)&&(yy<=P.Ny)&&(zz>=1)&&(zz<=P.Nz)) {
						d->m_vertex_type.ptr[varind]=2; //inner interface vertex
						if (dd==0) {
							FBVertexLocation VL;
							VL.x=xx;
							VL.y=yy;
							VL.z=zz;
							VL.ref_index=varind;
							d->m_inner_vertex_locations << VL;
						}
					}
					else {
						d->m_vertex_type.ptr[varind]=3; //outer interface vertex
						if (dd==0) {
							FBVertexLocation VL;
							VL.x=xx;
							VL.y=yy;
							VL.z=zz;
							VL.ref_index=varind;
							d->m_outer_vertex_locations << VL;
						}
					}
					d->m_x.ptr[varind]=P.X0.value(xx,yy,zz,dd);
				}
			}
		}
	}
	
	//set up the FBBlockElement list
	for (int zz=0; zz<P.Nz+1; zz++)
	for (int yy=0; yy<P.Ny+1; yy++)
	for (int xx=0; xx<P.Nx+1; xx++) {
		if (P.BVF.value(xx,yy,zz)) {
			FBBlockElement E0;
			E0.bvf=P.BVF.value(xx,yy,zz);
			E0.strain=0;
			E0.ref_indices[0]=(long)d->m_variable_indices.value(xx,yy,zz);
			E0.ref_indices[1]=(long)d->m_variable_indices.value(xx,yy+1,zz);
			E0.ref_indices[2]=(long)d->m_variable_indices.value(xx,yy,zz+1);
			E0.ref_indices[3]=(long)d->m_variable_indices.value(xx,yy+1,zz+1);
			d->m_elements << E0;
		}
	}
	
	//initialize r = -Ax (note that x is defined even on the fixed variables, so we don't need b)
	d->multiply_by_A(d->m_r,d->m_x); //r=Ax
	for (long ii=0; ii<d->m_num_variables; ii++) {
		d->m_r.ptr[ii]*=-1; //multiply r by -1: r = -Ax
	}
	//r is not defined on the outer interface; zeros there
	P.rnorm2=d->inner_product_on_owned_free_variables(d->m_r,d->m_r); // to compare with bnorm2 as the reference norm
	
	if (d->m_use_precondioner) {
		d->m_preconditioner.allocate(d->m_num_variables);	
		d->compute_preconditioner(d->m_x,d->m_preconditioner);
	}
	
	//define p equal to r on the free variables only; zeros everywhere else
	for (long ii=0; ii<d->m_num_variables; ii++) {
		if (d->m_free.ptr[ii]) {
			if ((d->m_use_precondioner)&&(d->m_preconditioner.ptr[ii])) {
				d->m_p.ptr[ii]=d->m_r.ptr[ii]/d->m_preconditioner.ptr[ii];
			}
			else d->m_p.ptr[ii]=d->m_r.ptr[ii];
		}
	}
	//here, p is not defined on outer interface
	
	//set p on the inner interface (free variables only)
	/*P.p_on_inner_interface.allocate(DATA_TYPE_FLOAT,3,P.Nx+2,P.Ny+2,P.Nz+2);
	for (int pass=1; pass<=3; pass++)
	for (int ii=0; ii<d->m_inner_vertex_locations.count(); ii++) { 
		FBVertexLocation *VL=&d->m_inner_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				if (pass<=2) P.p_on_inner_interface.setupIndex(pass,dd,VL->x,VL->y,VL->z);
				else if (pass==3) {
					P.p_on_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x,VL->y,VL->z);
				}
			}
		}
	}*/
	P.p_on_top_inner_interface.allocate(3,d->m_Nx,d->m_Ny);
	P.p_on_bottom_inner_interface.allocate(3,d->m_Nx,d->m_Ny);
	for (int ii=0; ii<d->m_inner_vertex_locations.count(); ii++) { 
		FBVertexLocation *VL=&d->m_inner_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				if (VL->z==1) P.p_on_top_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x-1,VL->y-1);
				else if (VL->z==d->m_Nz) P.p_on_bottom_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x-1,VL->y-1);
			}
		}
	}
}
void FBBlock::iterate_step_A(FBBlockIterateStepAParameters &P) {
	//update p on the outer interface (only free variables)
	/*for (int ii=0; ii<d->m_outer_vertex_locations.count(); ii++) {
		FBVertexLocation *VL=&d->m_outer_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				d->m_p.ptr[varind]=P.p_on_outer_interface.value(dd,VL->x,VL->y,VL->z);
			}
		}
	}*/
	
	if (d->m_nonlinear_adjuster) {
		//if we are adjusting the A matrix, like in a nonlinear simulation, then we NEED to reinitialize the residual at the beginning of each iteration!
		//initialize r = -Ax (note that x is defined even on the fixed variables, so we don't need b)
		d->multiply_by_A(d->m_r,d->m_x); //r=Ax
		for (long ii=0; ii<d->m_num_variables; ii++) {
			d->m_r.ptr[ii]*=-1; //multiply r by -1: r = -Ax
		}
	}
	
	for (int ii=0; ii<d->m_outer_vertex_locations.count(); ii++) {
		FBVertexLocation *VL=&d->m_outer_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				if (VL->z==0) d->m_p.ptr[varind]=P.p_on_top_outer_interface.value(dd,VL->x-1,VL->y-1);
				else if (VL->z==d->m_Nz+1) d->m_p.ptr[varind]=P.p_on_bottom_outer_interface.value(dd,VL->x-1,VL->y-1);
			}
		}
	}
	//now p is defined everywhere
	
	//update x,r,p, and Ap according to alpha and beta
	//now, p is defined everywhere, so we can do the following multiplication
	FBTimer::startTimer(QString("step_A_multipy_by_A-thread-%1").arg(d->m_block_id));
	d->multiply_by_A(d->m_Ap,d->m_p); 
	FBTimer::stopTimer(QString("step_A_multipy_by_A-thread-%1").arg(d->m_block_id));
	//now Ap is defined on the owned vertices
	
	//here's the output
	FBTimer::startTimer(QString("step_A_inner_products-thread-%1").arg(d->m_block_id));
	if (d->m_use_precondioner) {
		P.r_z=d->inner_product_on_owned_free_variables(d->m_r,d->m_r,d->m_preconditioner);
		P.r_Ap=d->inner_product_on_owned_free_variables(d->m_r,d->m_Ap,d->m_preconditioner);
		P.Ap_Ap=d->inner_product_on_owned_free_variables(d->m_Ap,d->m_Ap,d->m_preconditioner);
	}
	else {
		P.r_z=d->inner_product_on_owned_free_variables(d->m_r,d->m_r);
		P.r_Ap=d->inner_product_on_owned_free_variables(d->m_r,d->m_Ap);
		P.Ap_Ap=d->inner_product_on_owned_free_variables(d->m_Ap,d->m_Ap);
	}
	P.p_Ap=d->inner_product_on_owned_free_variables(d->m_p,d->m_Ap);
	FBTimer::stopTimer(QString("step_A_inner_products-thread-%1").arg(d->m_block_id));
}
void FBBlock::iterate_step_B(FBBlockIterateStepBParameters &P) {
	FBTimer::startTimer(QString("step_B_update_p-thread-%1").arg(d->m_block_id));
	
	for (long ii=0; ii<d->m_num_variables; ii++) {
		d->m_r.ptr[ii]=d->m_r.ptr[ii]-d->m_Ap.ptr[ii]*P.alpha; //r is never valid on the outer interface
		if (d->m_free.ptr[ii]) {
			d->m_x.ptr[ii]=d->m_x.ptr[ii]+d->m_p.ptr[ii]*P.alpha; //x is valid everywhere			
			if ((d->m_use_precondioner)&&(d->m_preconditioner.ptr[ii])) {
				qDebug() << __FUNCTION__ << __LINE__ << "using preconditioner";
				d->m_p.ptr[ii]=d->m_p.ptr[ii]*P.beta+d->m_r.ptr[ii]/d->m_preconditioner.ptr[ii]; //p is now not valid on the outer interface
			}
			else {
				d->m_p.ptr[ii]=d->m_p.ptr[ii]*P.beta+d->m_r.ptr[ii]; //p is now not valid on the outer interface
			}
		}
	}
	
	FBTimer::stopTimer(QString("step_B_update_p-thread-%1").arg(d->m_block_id));

	FBTimer::startTimer(QString("step_B_p_inner_products-thread-%1").arg(d->m_block_id));
	P.r_r=d->inner_product_on_owned_free_variables(d->m_r,d->m_r);
	P.bb_bb=d->inner_product_on_owned_fixed_variables(d->m_r,d->m_r);
	FBTimer::stopTimer(QString("step_B_p_inner_products-thread-%1").arg(d->m_block_id));
	
	FBTimer::startTimer(QString("step_B_compute_stress-thread-%1").arg(d->m_block_id));
	P.stress=d->compute_stress();
	FBTimer::stopTimer(QString("step_B_compute_stress-thread-%1").arg(d->m_block_id));

	//here's the output
	//set p on the inner interface (free variables only)
	/*P.p_on_inner_interface.allocate(DATA_TYPE_FLOAT,3,d->m_Nx+2,d->m_Ny+2,d->m_Nz+2);
	for (int pass=1; pass<=3; pass++)
	for (int ii=0; ii<d->m_inner_vertex_locations.count(); ii++) {
		FBVertexLocation *VL=&d->m_inner_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				if (pass<=2) P.p_on_inner_interface.setupIndex(pass,dd,VL->x,VL->y,VL->z);
				else if (pass==3) P.p_on_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x,VL->y,VL->z);
			}
		}
	}*/
	
	FBTimer::startTimer(QString("step_B_p_on_inner_interface-thread-%1").arg(d->m_block_id));
	P.p_on_top_inner_interface.allocate(3,d->m_Nx,d->m_Ny);
	P.p_on_bottom_inner_interface.allocate(3,d->m_Nx,d->m_Ny);
	for (int ii=0; ii<d->m_inner_vertex_locations.count(); ii++) {
		FBVertexLocation *VL=&d->m_inner_vertex_locations[ii];
		for (int dd=0; dd<3; dd++) {
			long varind=VL->ref_index+dd;
			if (d->m_free.ptr[varind]) {
				if (VL->z==1) P.p_on_top_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x-1,VL->y-1);
				else if (VL->z==d->m_Nz) P.p_on_bottom_inner_interface.setValue(d->m_p.ptr[varind],dd,VL->x-1,VL->y-1);
			}
		}
	}
	FBTimer::stopTimer(QString("step_B_p_on_inner_interface-thread-%1").arg(d->m_block_id));
	
	if (d->m_nonlinear_adjuster) {
		FBTimer::startTimer(QString("step_B_compute_strains-%1").arg(d->m_block_id));
		FBSparseArray4D energy_map0;
		computeEnergyMap(energy_map0);
		//bool strain_is_bad=false;
		long ct=0;
		for (long zz=0; zz<d->m_Nz+1; zz++)
		for (long yy=0; yy<d->m_Ny+1; yy++)
		for (long xx=0; xx<d->m_Nx+1; xx++)  {
			if (d->m_bvf_map.value(xx,yy,zz)) {
				float energy0=energy_map0.value(0,xx,yy,zz);
				float bvf_factor=d->m_elements[ct].bvf*1.0/100;
				float YM=d->m_youngs_modulus;
				float voxel_volume=d->m_voxel_volume;
				d->m_elements[ct].strain=sqrt(2*qAbs(energy0)/(voxel_volume*YM*bvf_factor));
				//if ((d->m_elements[ct].strain)&&((d->m_elements[ct].strain>1.15)||(d->m_elements[ct].strain<0.85))) strain_is_bad=true;
				/*float len_x=-getDisplacement(xx+1,yy,zz,0)+getDisplacement(xx,yy,zz,0)+1; //this needs to be fixed, for resolution!
				float len_y=-getDisplacement(xx,yy+1,zz,1)+getDisplacement(xx,yy,zz,1)+1;
				float len_z=-getDisplacement(xx,yy,zz+1,2)+getDisplacement(xx,yy,zz,2)+1;
				float vol0=len_x*len_y*len_z;
				float volume_ratio=qAbs(vol0/d->m_voxel_volume);
				if (volume_ratio<1) d->m_elements[ct].strain=1-volume_ratio;
				else d->m_elements[ct].strain=0;
				if (ct==0) qDebug() << "@#@#@#@#" << len_x << len_y << len_z << vol0 << volume_ratio << d->m_elements[ct].strain;*/
				ct++;
			}
		}
		/*if (strain_is_bad) {
			write_mda(energy_map0,"c:/dev/projects/fbblock/test/energy.mda");
			FBArray4D<float> displacements;
			displacements.allocate(3,d->m_Nx+2,d->m_Ny+2,d->m_Nz+2);
			for (long zz=0; zz<d->m_Nz+2; zz++)
			for (long yy=0; yy<d->m_Ny+2; yy++)
			for (long xx=0; xx<d->m_Nx+2; xx++) 
			for (int dd=0; dd<3; dd++) {
				displacements.setValue(getDisplacement(xx,yy,zz,dd),dd,xx,yy,zz);
			}
			write_mda(displacements,"c:/dev/projects/fbblock/test/displacements.mda");
			
			exit(0);
		}*/
		FBTimer::stopTimer(QString("step_B_compute_strains-%1").arg(d->m_block_id));
		
	}
}
	
double FBBlockPrivate::inner_product_on_owned_free_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2) {
	double ret=0;
	for (long ii=0; ii<m_num_variables; ii++) {
		if ((m_vertex_type.ptr[ii]!=3)&&(m_free.ptr[ii])) {
			ret+=V1.ptr[ii]*V2.ptr[ii];
		}
	}
	return ret;
}
double FBBlockPrivate::inner_product_on_owned_free_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2,const FBArray1D<float> &V3) {
	double ret=0;
	for (long ii=0; ii<m_num_variables; ii++) {
		if ((m_vertex_type.ptr[ii]!=3)&&(m_free.ptr[ii])) {
      			if (V3.ptr[ii]) ret+=V1.ptr[ii]*V2.ptr[ii]/V3.ptr[ii];
      			else ret+=V1.ptr[ii]*V2.ptr[ii];
		}
	}
	return ret;
}
	
double FBBlockPrivate::inner_product_on_owned_fixed_variables(const FBArray1D<float> &V1,const FBArray1D<float> &V2) {
	double ret=0;
	for (long ii=0; ii<m_num_variables; ii++) {
		if ((m_vertex_type.ptr[ii]!=3)&&(!m_free.ptr[ii])) {
			ret+=V1.ptr[ii]*V2.ptr[ii];
		}
	}
	return ret;
}
void FBBlockPrivate::multiply_by_A(FBArray1D<float> &Y,const FBArray1D<float> &X) { //Y=AX
	Y.setAll(0);
	
	float *stiffness_matrix_data=(float *)malloc(sizeof(float)*24*24);
	int ct=0;
	for (int rr=0; rr<24; rr++)
	for (int cc=0; cc<24; cc++) {
		stiffness_matrix_data[ct]=m_stiffness_matrix.value(rr,cc);
		ct++;
	}
			
	for (long i=0; i<m_elements.count(); i++) {
		float X0[24];
		float Y0[24];
		long varinds[24];
		FBBlockElement *E0=&m_elements[i];
		{
			{
				for (int kk=0; kk<4; kk++) {
					for (int jj=0; jj<6; jj++)
						varinds[kk*6+jj]=E0->ref_indices[kk]+jj;
				}
				for (int kk=0; kk<24; kk++) {
					X0[kk]=X.ptr[varinds[kk]];
					Y0[kk]=0;
				}
			}
			//optimized matrix multiplication
			int ct=0;
			for (int rr=0; rr<24; rr++)
			for (int cc=0; cc<24; cc++) {
				Y0[rr]+=stiffness_matrix_data[ct]*X0[cc];
				ct++;
			}
			float bvf_factor=E0->bvf*1.0/100;
			if (m_nonlinear_adjuster) bvf_factor*=m_nonlinear_adjuster->computeAdjustment(E0->strain);
			//if ((i==0)&&(m_nonlinear_adjuster)) qDebug() << "The bvf factor is:" << bvf_factor << "for a strain of" << E0->strain;
			for (int kk=0; kk<24; kk++) {
				if (m_vertex_type.ptr[varinds[kk]]!=3) {
					Y.ptr[varinds[kk]]+=Y0[kk]*bvf_factor;
				}
			}
		}
	}
	
	free(stiffness_matrix_data);
}

void FBBlockPrivate::compute_preconditioner(const FBArray1D<float> &X,FBArray1D<float> &C) { 
	for (long i=0; i<m_elements.count(); i++) {
		FBBlockElement *E0=&m_elements[i];
		long varinds[24];
		for (int kk=0; kk<4; kk++) {
			for (int jj=0; jj<6; jj++)
				varinds[kk*6+jj]=E0->ref_indices[kk]+jj;
		}
		float X0[24];
		float Y0[24];
		float C0[24];
		for (int kk=0; kk<24; kk++) {
			if (!m_free.ptr[varinds[kk]])
				X0[kk]=X.ptr[varinds[kk]];
			else 
				X0[kk]=0; 
			Y0[kk]=0;
			C0[kk]=0;
		}
		float bvf_factor=E0->bvf*1.0/100;
		if (m_nonlinear_adjuster) bvf_factor*=m_nonlinear_adjuster->computeAdjustment(E0->strain);
		for (int rr=0; rr<24; rr++) {
			for (int cc=0; cc<24; cc++) {
				Y0[rr]+=m_stiffness_matrix.value(rr,cc)*bvf_factor*X0[cc];
				if (cc==rr) C0[rr]=m_stiffness_matrix.value(rr,cc)*bvf_factor;
			}
		}
		for (int kk=0; kk<24; kk++) {
			if ((m_vertex_type.ptr[varinds[kk]]!=3)&&(m_free.ptr[varinds[kk]])) {
				C.ptr[varinds[kk]]+=C0[kk];
			}
		}
	}	
}

float FBBlock::getDisplacement(int xx,int yy,int zz,int dd) {	
	long varind=d->m_variable_indices.value(xx,yy,zz); 
	if (varind>=0) varind+=dd;
	else return 0;
	return d->m_x.ptr[varind];	
}
float FBBlock::getForce(int xx,int yy,int zz,int dd) {	
	long varind=d->m_variable_indices.value(xx,yy,zz); 
	if (varind>=0) varind+=dd;
	else return 0;
	return d->m_r.ptr[varind];	
}

long FBBlock::ownedFreeVariableCount() {
	long ret=0;
	for (long i=0; i<d->m_free.length();i++)
		if ((d->m_free.ptr[i])&&(d->m_vertex_type.ptr[i]!=3)) ret++;
	return ret;
}
long FBBlock::variableCount() {
	return d->m_num_variables;
}
void FBBlock::clearArrays() {
	d->m_free.clear();
	d->m_Ap.clear();
	d->m_p.clear();
	d->m_vertex_type.clear();
	d->m_elements.clear();
	d->m_inner_vertex_locations.clear();
	d->m_outer_vertex_locations.clear();
}
void FBBlock::clearArrays2() {
	d->m_x.clear();
	d->m_r.clear();
	d->m_variable_indices.clear();
}
void FBBlock::setResolution(QList<float> &res) {
	for (int i=0; i<3; i++) d->m_resolution[i]=res[i];
}
QList<double> FBBlockPrivate::compute_stress() {
	QList<double> ret;
	for (int j=0; j<6; j++) ret << 0;
	//was there a bug here? used to go up to i3<m_Nz+1
	for (long i3=0; i3<m_Nz; i3++) 
	for (long i2=0; i2<m_Ny; i2++) 
	for (long i1=0; i1<m_Nx; i1++) {
		float fx=q->getForce(i1+1,i2+1,i3+1,0);
		float fy=q->getForce(i1+1,i2+1,i3+1,1);
		float fz=q->getForce(i1+1,i2+1,i3+1,2);
		if ((fx)||(fy)||(fz)) {
			ret[0]+=fx*(m_block_x_position+i1)*m_resolution[0]; //sigma_11
			ret[1]+=fy*(m_block_y_position+i2)*m_resolution[1]; //sigma_22
			ret[2]+=fz*(m_block_z_position+i3)*m_resolution[2]; //sigma_33
			ret[3]+=(fx*(m_block_y_position+i2)*m_resolution[1]+fy*(m_block_x_position+i1)*m_resolution[0])*0.5; //sigma_12
			ret[4]+=(fx*(m_block_z_position+i3)*m_resolution[2]+fz*(m_block_x_position+i1)*m_resolution[0])*0.5; //sigma_13
			ret[5]+=(fy*(m_block_z_position+i3)*m_resolution[2]+fz*(m_block_y_position+i2)*m_resolution[1])*0.5; //sigma_23
		}
	}
	return ret;
}
long FBBlock::ownedVariableCount() {
	long ret=0;
	for (long i=0; i<d->m_free.length();i++)
		if (d->m_vertex_type.ptr[i]!=3) ret++;
	return ret;
}
QList<double> FBBlock::getStress() {
	return d->compute_stress();
}

void FBBlock::computeEnergyMap(FBSparseArray4D &E) {
	float *stiffness_matrix_data=(float *)malloc(sizeof(float)*24*24);
	int ct=0;
	for (int rr=0; rr<24; rr++)
	for (int cc=0; cc<24; cc++) {
		stiffness_matrix_data[ct]=d->m_stiffness_matrix.value(rr,cc);
		ct++;
	}	
	
	E.allocate(DATA_TYPE_FLOAT,1,d->m_Nx+1,d->m_Ny+1,d->m_Nz+1);
	for (int pass=1; pass<=2; pass++) {
		for (long i3=0; i3<d->m_Nz+1; i3++)
		for (long i2=0; i2<d->m_Ny+1; i2++)
		for (long i1=0; i1<d->m_Nx+1; i1++)  {
			if (is_element(d->m_bvf_map,i1,i2,i3)) {
				E.setupIndex(pass,0,i1,i2,i3);
			}
		}
	}	
	
	for (long zz=0; zz<d->m_Nz+1; zz++)
	for (long yy=0; yy<d->m_Ny+1; yy++)
	for (long xx=0; xx<d->m_Nx+1; xx++)  {
		if (is_element(d->m_bvf_map,xx,yy,zz)) {
			fbreal bvf=d->m_bvf_map.value(xx,yy,zz);
			long ref_indices[4];
			ref_indices[0]=(long)d->m_variable_indices.value(xx,yy,zz);
			ref_indices[1]=(long)d->m_variable_indices.value(xx,yy+1,zz);
			ref_indices[2]=(long)d->m_variable_indices.value(xx,yy,zz+1);
			ref_indices[3]=(long)d->m_variable_indices.value(xx,yy+1,zz+1);
			
			float X0[24];
			long varinds[24];
			{
				for (int kk=0; kk<4; kk++) {
					for (int jj=0; jj<6; jj++)
						varinds[kk*6+jj]=ref_indices[kk]+jj;
				}
				for (int kk=0; kk<24; kk++) {
					X0[kk]=d->m_x.ptr[varinds[kk]];
				}
			}
			//optimized matrix multiplication
			double energy0=0;
			int ct=0;
			for (int rr=0; rr<24; rr++)
			for (int cc=0; cc<24; cc++) {
				energy0+=stiffness_matrix_data[ct]*X0[cc]*X0[rr];
				ct++;
			}
			energy0*=bvf*1.0/100*(-1)*0.5;
			E.setValue(energy0,0,xx,yy,zz);
			/*fbreal bvf=d->m_bvf_map.value(i1,i2,i3);
			for (int a3=0; a3<=1; a3++)
			for (int a2=0; a2<=1; a2++)
			for (int a1=0; a1<=1; a1++) 
			for (int adir=0; adir<3; adir++) {
				long avarind=d->m_variable_indices.value(i1+a1,i2+a2,i3+a3)+adir;
				for (int b3=0; b3<=1; b3++)
				for (int b2=0; b2<=1; b2++)
				for (int b1=0; b1<=1; b1++) 
				for (int bdir=0; bdir<3; bdir++) {
					long bvarind=d->m_variable_indices.value(i1+b1,i2+b2,i3+b3)+bdir;
					float val0=d->m_x.ptr[bvarind]
							*d->m_x.ptr[avarind]
							*d->m_stiffness_matrix.value(adir+a1*3+a2*6+a3*12,bdir+b1*3+b2*6+b3*12)
							*bvf/100*(-1);
					E.setValue(
							E.value(0,i1,i2,i3)+val0,
							0,i1,i2,i3);
				}
			}*/
		}
	}
	free(stiffness_matrix_data);
}
int FBBlock::Nx() const {
	return d->m_Nx;
}
int FBBlock::Ny() const {
	return d->m_Ny;
}
int FBBlock::Nz() const {
	return d->m_Nz;
}
int FBBlock::xPosition() const {
	return d->m_block_x_position;
}
int FBBlock::yPosition() const {
	return d->m_block_y_position;
}
int FBBlock::zPosition() const {
	return d->m_block_z_position;
}
void FBBlock::setNonlinearAdjuster(NonlinearAdjuster *X) {
	d->m_nonlinear_adjuster=X;
}

