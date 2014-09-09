#include "fbblocksolver.h"
#include "fbblock.h"
#include <QList>
#include <QQueue>
#include <stdio.h>
#include <QFile>
#include <QDebug>
#include <QTime>
#include <math.h>
#include <QThread>
#include <QTest>
#include "mda_io.h"
#include "fbtimer.h"
#include "nonlinearadjuster.h"

struct BlockInfo {
	int xmin,xmax;
	int ymin,ymax;
	int zmin,zmax;
	FBArray3D<float> p_on_top_inner_interface;
	FBArray3D<float> p_on_bottom_inner_interface;
};

class FBBlockSolverPrivate {
public:
	FBBlockSolver *q;
	//FBSparseArray4D m_p; //the p-vector, i.e. search direction, only defined on the interface vertices
	FBArray2D <float> m_stiffness_matrix;
	float m_youngs_modulus;
	float m_voxel_volume;
	FBArray3D <unsigned char> m_bvf_map; //N1 x N2 x N3
	FBSparseArray4D m_initial_displacements; //3x(N1+1)x(N2+1)x(N3+1)
	FBSparseArray4D m_fixed_variables; //3x(N1+1)x(N2+1)x(N3+1)
	QVector<FBBlock *> m_blocks;
	long m_num_iterations;
	FBErrorEstimator m_error_estimator;
	QList<FBBlockIterateStepAParameters> m_PPP_A;
	QList<FBBlockIterateStepBParameters> m_PPP_B;
	QList<BlockInfo> m_block_infos;
	
	fbreal m_epsilon;
	int m_max_iterations;
	int m_num_threads;	
	bool m_use_precondioner;
	float m_resolution[3];
	
	NonlinearAdjuster *m_nonlinear_adjuster;
	
	void do_iterations();
};

FBBlockSolver::FBBlockSolver() 
{
	d=new FBBlockSolverPrivate;
	d->q=this;
	d->m_youngs_modulus=1;
	d->m_voxel_volume=1;
	d->m_epsilon=0.001F;
	d->m_max_iterations=0;
	d->m_num_threads=1;
	d->m_use_precondioner=false;
	d->m_nonlinear_adjuster=0;
	for (int i=0; i<3; i++) d->m_resolution[i]=1;
	
	d->m_num_iterations=0;
}

FBBlockSolver::~FBBlockSolver()
{
	qDeleteAll(d->m_blocks);
	delete d;
}
void FBBlockSolver::setEpsilon(fbreal epsilon) {d->m_epsilon=epsilon;}
void FBBlockSolver::setMaxIterations(int val) {d->m_max_iterations=val;}
void FBBlockSolver::setNumThreads(int val) {d->m_num_threads=val;}
void FBBlockSolver::setUsePreconditioner(bool val) {d->m_use_precondioner=val;}
void FBBlockSolver::setStiffnessMatrix(const FBArray2D<float> &stiffness_matrix) {
	d->m_stiffness_matrix=stiffness_matrix;
}
void FBBlockSolver::setYoungsModulus(float val) {
	d->m_youngs_modulus=val;
}
void FBBlockSolver::setVoxelVolume(float val) {
	d->m_voxel_volume=val;
}
void FBBlockSolver::setBVFMap(const FBArray3D<unsigned char> &bvf_map) {
	d->m_bvf_map=bvf_map;
}
void FBBlockSolver::setInitialDisplacementsOnFreeVariables(const FBSparseArray4D &displacements) {
	d->m_initial_displacements.resetIteration();
	while (d->m_initial_displacements.advanceIteration()) {
		int i1=d->m_initial_displacements.currentIndex1();
		int i2=d->m_initial_displacements.currentIndex2();
		int i3=d->m_initial_displacements.currentIndex3();
		int i4=d->m_initial_displacements.currentIndex4();
		if (!d->m_fixed_variables.value(i1,i2,i3,i4)) {
			d->m_initial_displacements.setValue(displacements.value(i1,i2,i3,i4),i1,i2,i3,i4);
		}
	}
}
void FBBlockSolver::setInitialDisplacementsOnFreeVariables(const FBArray4D<float> &displacements) {
	d->m_initial_displacements.resetIteration();
	while (d->m_initial_displacements.advanceIteration()) {
		int i1=d->m_initial_displacements.currentIndex1();
		int i2=d->m_initial_displacements.currentIndex2();
		int i3=d->m_initial_displacements.currentIndex3();
		int i4=d->m_initial_displacements.currentIndex4();
		if (!d->m_fixed_variables.value(i1,i2,i3,i4)) {
			d->m_initial_displacements.setValue(displacements.value(i1,i2,i3,i4),i1,i2,i3,i4);
		}
	}
}

void FBBlockSolver::setFixedVariables(const FBSparseArray4D &fixed_variables) {
	d->m_fixed_variables=fixed_variables;
}

bool is_on_an_interface(long x,long y,long z,int block_size) {
	if (x<0) return 0;
	if (y<0) return 0;
	if (z<0) return 0;
	int x2=x%block_size;
	int y2=y%block_size;
	int z2=z%block_size;
	if ((x2==0)||(x2==block_size-1)) return true;
	if ((y2==0)||(y2==block_size-1)) return true;
	if ((z2==0)||(z2==block_size-1)) return true;
	return false;
}

class FBBlockSolverThread : public QThread {
public:
	QList<FBBlockIterateStepAParameters *> step_A_parameters;
	QList<FBBlockIterateStepBParameters *> step_B_parameters;
	QList<FBBlock *> blocks;
	bool do_step_A;
	bool do_step_B;
	
	FBBlockSolverThread() {
		do_step_A=false;
		do_step_B=false;
	}
	void run() {
		if (do_step_A) {
			for (int i=0; i<blocks.count(); i++) {
				blocks[i]->iterate_step_A(*(step_A_parameters[i]));
			}
		}
		else if (do_step_B) {
			for (int i=0; i<blocks.count(); i++) {
				blocks[i]->iterate_step_B(*(step_B_parameters[i]));
			}
		}
	}
};

bool is_on_an_interface(long x,long y,long z,const QList<BlockInfo> &infos) {
	for (int i=0; i<infos.count(); i++) {
		BlockInfo II=infos[i];
		if ((II.xmin==x)||(II.xmin-1==x)) return true;
		if ((II.xmax==x)||(II.xmax+1==x)) return true;
		if ((II.ymin==y)||(II.ymin-1==y)) return true;
		if ((II.ymax==y)||(II.ymax+1==y)) return true;
		if ((II.zmin==z)||(II.zmin-1==z)) return true;
		if ((II.zmax==z)||(II.zmax+1==z)) return true;
	}
	return false;
}

void FBBlockSolver::solve() {
	FBTimer::startTimer("solve");
	
	FBTimer::startTimer("setup");
	//number of blocks needed = (Nx) x (Ny) x (Nz)
	int N1=d->m_bvf_map.N1();
	int N2=d->m_bvf_map.N2();
	int N3=d->m_bvf_map.N3();
	
	QList<long> slice_vertex_count;
	long total_vertex_count=0;
	for (int z=0; z<N3+1; z++) {
		long vertex_count=0;
		for (int y=0; y<N2+1; y++)
		for (int x=0; x<N1+1; x++) {
			if (is_vertex(d->m_bvf_map,x,y,z)) {
				vertex_count++;
				total_vertex_count++;
			}
		}
		slice_vertex_count << vertex_count;
	}

	d->m_block_infos.clear();
	double num_vertices_per_thread=((double)total_vertex_count)*1.0/d->m_num_threads;
	int z0=-1;
	for (int ithread=0; ithread<d->m_num_threads; ithread++) {
		BlockInfo block_info0;
		block_info0.xmin=-1; block_info0.xmax=N1+1;
		block_info0.ymin=-1; block_info0.ymax=N2+1;
		block_info0.zmin=z0+1;
		if (ithread==d->m_num_threads-1) {
			block_info0.zmax=N3;
		}
		else {
			int zmax=z0+1;
			if (zmax<N3+1) {
				long block_vertex_count=slice_vertex_count[zmax];
				bool done=false;
				while (!done) {
					if (zmax>=N3) done=true;
					else {
						double diff0=qAbs(block_vertex_count-num_vertices_per_thread);
						zmax++;
						if (zmax<=N3) block_vertex_count+=slice_vertex_count[zmax];
						double diff1=qAbs(block_vertex_count-num_vertices_per_thread);
						if (diff1>=diff0) {
							zmax--;
							done=true;
						}
					}
				}
			}
			block_info0.zmax=zmax;
			z0=block_info0.zmax;
		}
		if (block_info0.zmin<=N3) d->m_block_infos << block_info0;
	}
	
	qDeleteAll(d->m_blocks);
	d->m_blocks.clear();
	
	//Set up the blocks and d->m_p
	/*d->m_p.allocate(DATA_TYPE_FLOAT,3,d->m_bvf_map.N1()+1,d->m_bvf_map.N2()+1,d->m_bvf_map.N3()+1);
	for (int pass=1; pass<=2; pass++) {
		for (long z=0; z<d->m_bvf_map.N3()+1; z++)
		for (long y=0; y<d->m_bvf_map.N2()+1; y++)
		for (long x=0; x<d->m_bvf_map.N1()+1; x++) {
			if (is_on_an_interface(x,y,z,d->m_block_infos)) {
				if (is_vertex(d->m_bvf_map,x,y,z)) {
					for (int dd=0; dd<3; dd++) {
						d->m_p.setupIndex(pass,dd,x,y,z);
					}
				}
			}
		} 
	}*/

	long num_variables=0;

	for (int iii=0; iii<d->m_block_infos.count(); iii++) {
		FBBlock *B=new FBBlock(iii);
		FBBlockSetupParameters PP;
		PP.use_preconditioner=d->m_use_precondioner;
		BlockInfo Info0=d->m_block_infos[iii];
		for (int i=0; i<3; i++) PP.resolution[i]=d->m_resolution[i];
		PP.Nx=Info0.xmax-Info0.xmin+1;
		PP.Ny=Info0.ymax-Info0.ymin+1;
		PP.Nz=Info0.zmax-Info0.zmin+1;
		PP.block_x_position=Info0.xmin;
		PP.block_y_position=Info0.ymin;
		PP.block_z_position=Info0.zmin; 
		//BVF
		PP.BVF.allocate(PP.Nx+1,PP.Ny+1,PP.Nz+1);
		for (int zz=Info0.zmin-1; zz<Info0.zmax+1; zz++)
		for (int yy=Info0.ymin-1; yy<Info0.ymax+1; yy++)
		for (int xx=Info0.xmin-1; xx<Info0.xmax+1; xx++) {
			int xx0=xx-(Info0.xmin-1); int yy0=yy-(Info0.ymin-1); int zz0=zz-(Info0.zmin-1);
			unsigned char val=d->m_bvf_map.value(xx,yy,zz);
			PP.BVF.setValue(val,xx0,yy0,zz0);
		}
		//fixed variables
		PP.fixed.allocate(PP.Nx+2,PP.Ny+2,PP.Nz+2,3);
		for (int zz=Info0.zmin-1; zz<=Info0.zmax+1; zz++)
		for (int yy=Info0.ymin-1; yy<=Info0.ymax+1; yy++)	
		for (int xx=Info0.xmin-1; xx<=Info0.xmax+1; xx++)
		for (int dd=0; dd<3; dd++) {
			int xx0=xx-(Info0.xmin-1); int yy0=yy-(Info0.ymin-1); int zz0=zz-(Info0.zmin-1);
			if (d->m_fixed_variables.value(dd,xx,yy,zz)) {
				PP.fixed.setValue(1,xx0,yy0,zz0,dd);
			}
		}
		//stiffness_matrix
		PP.stiffness_matrix=d->m_stiffness_matrix;
		PP.youngs_modulus=d->m_youngs_modulus;
		PP.voxel_volume=d->m_voxel_volume;

		//Initial displacements		
		PP.X0.allocate(PP.Nx+2,PP.Ny+2,PP.Nz+2,3);
		for (int zz=Info0.zmin-1; zz<=Info0.zmax+1; zz++)
		for (int yy=Info0.ymin-1; yy<=Info0.ymax+1; yy++)	
		for (int xx=Info0.xmin-1; xx<=Info0.xmax+1; xx++)
		for (int dd=0; dd<3; dd++) {
			int xx0=xx-(Info0.xmin-1); int yy0=yy-(Info0.ymin-1); int zz0=zz-(Info0.zmin-1);
			PP.X0.setValue(d->m_initial_displacements.value(dd,xx,yy,zz),xx0,yy0,zz0,dd);   
		}
		//setup
		B->setup(PP);
		/*//set m_p on inner interfaces
		PP.p_on_inner_interface.resetIteration();
		while (PP.p_on_inner_interface.advanceIteration()) {
			int dd0=PP.p_on_inner_interface.currentIndex1();
			int xx0=PP.p_on_inner_interface.currentIndex2();
			int yy0=PP.p_on_inner_interface.currentIndex3();
			int zz0=PP.p_on_inner_interface.currentIndex4();			
			float val0=PP.p_on_inner_interface.currentValue();
			d->m_p.setValue(val0,dd0,Info0.xmin-1+xx0,Info0.ymin-1+yy0,Info0.zmin-1+zz0);
		}*/
		//store block in list		
		d->m_block_infos[iii].p_on_top_inner_interface=PP.p_on_top_inner_interface;
		d->m_block_infos[iii].p_on_bottom_inner_interface=PP.p_on_bottom_inner_interface;
		d->m_blocks << B;

		num_variables+=B->ownedFreeVariableCount();
	}
	printf("Total number of variables: %ld\n",num_variables);
	printf("Using %d blocks.\n",d->m_blocks.count());
	
	printf("Setting up the Step A Parameters...\n");
	d->m_PPP_A.clear();
	for (long i=0; i<d->m_blocks.count(); i++) {
		//FBBlock *B0=d->m_blocks[i];
		FBBlockIterateStepAParameters PP;
		//PP.p_on_outer_interface.allocate(DATA_TYPE_FLOAT,3,B0->Nx()+2,B0->Ny()+2,B0->Nz()+2);
		/*for (int pass=1; pass<=2; pass++) {
			for (int zz=0; zz<PP.p_on_outer_interface.N4(); zz++) 
			for (int yy=0; yy<PP.p_on_outer_interface.N3(); yy++) 
			for (int xx=0; xx<PP.p_on_outer_interface.N2(); xx++)
			for (int dd=0; dd<3; dd++) {
				if ((xx==0)||(xx==PP.p_on_outer_interface.N2()-1)
				   ||(yy==0)||(yy==PP.p_on_outer_interface.N3()-1)
				   ||(zz==0)||(zz==PP.p_on_outer_interface.N4()-1)) {
				   	//perhaps optimize this later by checking to see if this vertex is occupied
					PP.p_on_outer_interface.setupIndex(pass,dd,xx,yy,zz);
				}
			}
		}*/
		d->m_PPP_A << PP;
	}
	
	printf("Setting up the Step B Parameters...\n");
	d->m_PPP_B.clear();
	for (int i=0; i<d->m_blocks.count(); i++) {
		FBBlockIterateStepBParameters PP;
		d->m_PPP_B << PP;
	}
	FBTimer::stopTimer("setup");

	d->do_iterations();
	
	FBTimer::stopTimer("solve");
}

class MyNonlinearAdjuster : public NonlinearAdjuster {
public:
	float eps_yield;
	float computeAdjustment(float eps) {
		if (eps<2*eps_yield) return 1-(eps/(2*eps_yield))*0.95;
		else return 0.05;
		/*float transition_length=eps_yield*1;
		if (eps<eps_yield) return 1;
		else if (eps<eps_yield+transition_length) return 1-(eps-eps_yield)/(transition_length)*0.95;
		else return 0.05;*/
	}
};

void FBBlockSolver::solveNonlinear(float step_size,int num_steps,int num_iterations_per_step) {
	solve(); //first do the linear simulation
	
	
	MyNonlinearAdjuster adjuster;
	
	qDebug() << getStress();
	
	for (int ii=1; ii<=num_steps; ii++) {
		float eps=step_size*ii;
		adjuster.eps_yield=0.01/eps;
		d->m_max_iterations=num_iterations_per_step;		
		d->m_nonlinear_adjuster=&adjuster;
		d->m_num_iterations=0;
		d->m_epsilon=0;
		d->do_iterations();
		printf("%g\n",eps);
		qDebug() << getStress();
	}
}

long FBBlockSolver::setFixedVariables(FBMacroscopicStrain &macroscopic_strain) {
	long N1=d->m_bvf_map.N1();
	long N2=d->m_bvf_map.N2();
	long N3=d->m_bvf_map.N3();

	d->m_fixed_variables.allocate(DATA_TYPE_BYTE,3,N1+1,N2+1,N3+1);

	for (int pass=1; pass<=3; pass++)
	for (long i3=0; i3<N3+1; i3++)
	for (long i2=0; i2<N2+1; i2++)
	for (long i1=0; i1<N1+1; i1++) {
		if (is_vertex(d->m_bvf_map,i1,i2,i3)) {
			bool fix1=false,fix2=false,fix3=false;
			if ((i1==0) ||(i1==N1)) {
				if (macroscopic_strain.boundaryRestrictions[0][0])
					fix1=true;
				if (macroscopic_strain.boundaryRestrictions[0][1])
					fix2=true;
				if (macroscopic_strain.boundaryRestrictions[0][2])
					fix3=true;
			}
			if ((i2==0) ||(i2==N2)) {
				if (macroscopic_strain.boundaryRestrictions[1][0])
					fix1=true;
				if (macroscopic_strain.boundaryRestrictions[1][1])
					fix2=true;
				if (macroscopic_strain.boundaryRestrictions[1][2])
					fix3=true;
			}
			if ((i3==0) ||(i3==N3)) {
				if (macroscopic_strain.boundaryRestrictions[2][0])
					fix1=true;
				if (macroscopic_strain.boundaryRestrictions[2][1])
					fix2=true;
				if (macroscopic_strain.boundaryRestrictions[2][2])
					fix3=true;
			}
			if (pass<=2) {
				if (fix1) {d->m_fixed_variables.setupIndex(pass,0,i1,i2,i3);}
				if (fix2) {d->m_fixed_variables.setupIndex(pass,1,i1,i2,i3);}
				if (fix3) {d->m_fixed_variables.setupIndex(pass,2,i1,i2,i3);}
			}
			else if (pass==3) {
				if (fix1) {d->m_fixed_variables.setValue(1,0,i1,i2,i3);}
				if (fix2) {d->m_fixed_variables.setValue(1,1,i1,i2,i3);}
				if (fix3) {d->m_fixed_variables.setValue(1,2,i1,i2,i3);}
			}
		}	
	}
	
	printf("Counting elements...\n");
	long num_elements=0;
	for (long i3=0; i3<N3+1; i3++)
	for (long i2=0; i2<N2+1; i2++)
	for (long i1=0; i1<N1+1; i1++) {
		if (is_element(d->m_bvf_map,i1,i2,i3))
			num_elements++;
	}
	return(num_elements);
}
void FBBlockSolver::setInitialDisplacements(FBMacroscopicStrain &strain) {
	long N1=d->m_bvf_map.N1();
	long N2=d->m_bvf_map.N2();
	long N3=d->m_bvf_map.N3();
	d->m_initial_displacements.allocate(DATA_TYPE_FLOAT,3,N1+1,N2+1,N3+1);	
	for (int pass=1; pass<=3; pass++)
	for (long i3=0; i3<N3+1; i3++)
	for (long i2=0; i2<N2+1; i2++)
	for (long i1=0; i1<N1+1; i1++) {
		if (is_vertex(d->m_bvf_map,i1,i2,i3)) {
			for (int dd=0; dd<3; dd++) {
				if (pass<=2) d->m_initial_displacements.setupIndex(pass,dd,i1,i2,i3);
				else if (pass==3) {
					d->m_initial_displacements.setValue(initial_displacement(i1,i2,i3,dd,d->m_resolution,strain),dd,i1,i2,i3);
				}
			}
		}	
	}
}
void FBBlockSolver::getDisplacements(FBSparseArray4D &displacements) {
	//here we must retrieve the displacements from the blocks.
	
	displacements.allocate(DATA_TYPE_FLOAT,3,d->m_bvf_map.N1()+1,d->m_bvf_map.N2()+1,d->m_bvf_map.N3()+1);
	for (int pass=1; pass<=2; pass++) {
		for (long z=0; z<d->m_bvf_map.N3()+1; z++)
		for (long y=0; y<d->m_bvf_map.N2()+1; y++)
		for (long x=0; x<d->m_bvf_map.N1()+1; x++)
		if (is_vertex(d->m_bvf_map,x,y,z)) {
			for (int dd=0; dd<3; dd++) {
				displacements.setupIndex(pass,dd,x,y,z);
			}
		}
	}
	
	for (long i=0; i<d->m_blocks.count(); i++) {
		long x0=d->m_blocks[i]->xPosition();
		long y0=d->m_blocks[i]->yPosition();
		long z0=d->m_blocks[i]->zPosition();
		for (int kk=0; kk<d->m_blocks[i]->Nz(); kk++)
		for (int jj=0; jj<d->m_blocks[i]->Ny(); jj++)
		for (int ii=0; ii<d->m_blocks[i]->Nx(); ii++) {
			if (is_vertex(d->m_bvf_map,x0+ii,y0+jj,z0+kk)) {
				for (int dd=0; dd<3; dd++) { 
					fbreal displacement0=d->m_blocks[i]->getDisplacement(ii+1,jj+1,kk+1,dd);
					displacements.setValue(displacement0,dd,x0+ii,y0+jj,z0+kk);
				}
			}
		}
	}
}
void FBBlockSolver::getDisplacements(FBArray4D<float> &displacements) {
	//here we must retrieve the displacements from the blocks.
	
	displacements.allocate(3,d->m_bvf_map.N1()+1,d->m_bvf_map.N2()+1,d->m_bvf_map.N3()+1);
	
	for (long i=0; i<d->m_blocks.count(); i++) {
		long x0=d->m_blocks[i]->xPosition();
		long y0=d->m_blocks[i]->yPosition();
		long z0=d->m_blocks[i]->zPosition();
		for (int kk=0; kk<d->m_blocks[i]->Nz(); kk++)
		for (int jj=0; jj<d->m_blocks[i]->Ny(); jj++)
		for (int ii=0; ii<d->m_blocks[i]->Nx(); ii++) {
			if (is_vertex(d->m_bvf_map,x0+ii,y0+jj,z0+kk)) {
				for (int dd=0; dd<3; dd++) { 
					fbreal displacement0=d->m_blocks[i]->getDisplacement(ii+1,jj+1,kk+1,dd);
					displacements.setValue(displacement0,dd,x0+ii,y0+jj,z0+kk);
				}
			}
		}
	}
}
void FBBlockSolver::getForces(FBSparseArray4D &forces) {
	//here we must retrieve the forces from the blocks.
	
	forces.allocate(DATA_TYPE_FLOAT,3,d->m_bvf_map.N1()+1,d->m_bvf_map.N2()+1,d->m_bvf_map.N3()+1);
	for (int pass=1; pass<=2; pass++) {
		for (long z=0; z<d->m_bvf_map.N3()+1; z++)
		for (long y=0; y<d->m_bvf_map.N2()+1; y++)
		for (long x=0; x<d->m_bvf_map.N1()+1; x++)
		if (is_vertex(d->m_bvf_map,x,y,z)) {
			for (int dd=0; dd<3; dd++) {
				forces.setupIndex(pass,dd,x,y,z);
			}
		}
	}
	
	for (long i=0; i<d->m_blocks.count(); i++) {
		long x0=d->m_blocks[i]->xPosition();
		long y0=d->m_blocks[i]->yPosition();
		long z0=d->m_blocks[i]->zPosition();
		for (int kk=0; kk<d->m_blocks[i]->Nz(); kk++)
		for (int jj=0; jj<d->m_blocks[i]->Ny(); jj++)
		for (int ii=0; ii<d->m_blocks[i]->Nx(); ii++) {
			if (is_vertex(d->m_bvf_map,x0+ii,y0+jj,z0+kk)) {
				for (int dd=0; dd<3; dd++) { 
					fbreal force0=d->m_blocks[i]->getForce(ii+1,jj+1,kk+1,dd);
					forces.setValue(force0,dd,x0+ii,y0+jj,z0+kk);
				}
			}
		}
	}
}
void FBBlockSolver::setResolution(QList<fbreal> &res) {
	for (int i=0; i<3; i++) d->m_resolution[i]=res[i];
}
/*class GetStressThread : public QThread {
public:
	//input
	FBBlock *block;
	//output
	QList<double> output;
	void run() {
		output=block->getStress();
	}
};*/
/*QList<double> FBBlockSolver::getStress() {
	//this is the volume of the entire bvf map
	double stress_denominator=(d->m_bvf_map.N1()*d->m_bvf_map.N2()*d->m_bvf_map.N3()*d->m_resolution[0]*d->m_resolution[1]*d->m_resolution[2]);
	
	QList<GetStressThread *> threads;
	QList<double> ret;	
	for (int i=0; i<6; i++) ret << 0;
	int num_blocks=d->m_blocks.count();
	for (int ii=0; ii<num_blocks; ii++) {
		GetStressThread *T0=new GetStressThread;
		T0->block=d->m_blocks[ii];
		T0->start();
		threads << T0;
	}
	bool done=false;
	while (!done) {
		done=true;
		for (int ii=0; ii<num_blocks; ii++) {
			if (!threads[ii]->isFinished()) done=false;
		}
		QTest::qWait(10);
	}
	for (int ii=0; ii<num_blocks; ii++) {
		for (int i=0; i<6; i++) ret[i]+=threads[ii]->output[i];
	}
	qDeleteAll(threads);
	for (int i=0; i<6; i++) ret[i]/=stress_denominator;
	return ret;
}*/
QList<double> FBBlockSolver::getStress() {
	//this is the volume of the entire bvf map
	double stress_denominator=(d->m_bvf_map.N1()*d->m_bvf_map.N2()*d->m_bvf_map.N3()*d->m_resolution[0]*d->m_resolution[1]*d->m_resolution[2]);

	QList<double> ret;	
	for (int jj=0; jj<6; jj++) ret << 0;
	for (int ii=0; ii<d->m_blocks.count(); ii++) {
		double val0=0;
		for (int jj=0; jj<6; jj++) ret[jj]+=d->m_PPP_B[ii].stress[jj];
	}
	for (int jj=0; jj<6; jj++) ret[jj]/=stress_denominator;
	return ret;
}
int FBBlockSolver::getNumIterations() {
	return d->m_num_iterations;
}
void FBBlockSolver::getEnergy(FBSparseArray4D &energy) {
	//here we must retrieve the energy maps from the blocks.
	
	energy.allocate(DATA_TYPE_FLOAT,1,d->m_bvf_map.N1(),d->m_bvf_map.N2(),d->m_bvf_map.N3());
	for (int pass=1; pass<=2; pass++) {
		for (long z=0; z<d->m_bvf_map.N3(); z++)
		for (long y=0; y<d->m_bvf_map.N2(); y++)
		for (long x=0; x<d->m_bvf_map.N1(); x++)
		if (is_element(d->m_bvf_map,x,y,z)) {
			energy.setupIndex(pass,0,x,y,z);
		}
	}
	
	for (long i=0; i<d->m_blocks.count(); i++) {
		FBSparseArray4D energy_map0;
		d->m_blocks[i]->computeEnergyMap(energy_map0);
		long x0=d->m_blocks[i]->xPosition();
		long y0=d->m_blocks[i]->yPosition();
		long z0=d->m_blocks[i]->zPosition();
		for (int kk=0; kk<d->m_blocks[i]->Nz()+1; kk++)
		for (int jj=0; jj<d->m_blocks[i]->Ny()+1; jj++)
		for (int ii=0; ii<d->m_blocks[i]->Nx()+1; ii++) {
			if (is_element(d->m_bvf_map,x0-1+ii,y0-1+jj,z0-1+kk)) {
				fbreal energy0=energy_map0.value(0,ii,jj,kk);
				energy.setValue(energy0,0,x0-1+ii,y0-1+jj,z0-1+kk);
			}
		}
	}
}
FBErrorEstimator *FBBlockSolver::errorEstimator() {
	return &d->m_error_estimator;
}
void FBBlockSolverPrivate::do_iterations() {
	for (int ii=0; ii<m_blocks.count(); ii++) {
		m_blocks[ii]->setNonlinearAdjuster(m_nonlinear_adjuster);
	}

	FBTimer::startTimer("iterations");	
	int num_times_below_epsilon=0;
	while (((m_num_iterations<m_max_iterations)||(m_max_iterations<=0))&&(num_times_below_epsilon<5)) {
		FBTimer::startTimer("setup_for_A");
		//Set parameters for step A, i.e., p_on_outer_interface
		/*{
			for (long iblock=0; iblock<m_blocks.count(); iblock++) {
				FBBlock *B0=m_blocks[iblock];
				for (long zz=B0->zPosition()-1; zz<=B0->zPosition()+B0->Nz(); zz+=B0->Nz()+1) {
					for (long yy=B0->yPosition()-1; yy<=B0->yPosition()+B0->Ny(); yy++)
					for (long xx=B0->xPosition()-1; xx<=B0->xPosition()+B0->Nx(); xx++) 
					for (int dd=0; dd<3; dd++) {
						float val0=m_p.value(dd,xx,yy,zz);
						if (val0) {
							m_PPP_A[iblock].p_on_outer_interface.setValue(val0,dd,xx-(B0->xPosition()-1),yy-(B0->yPosition()-1),zz-(B0->zPosition()-1));
						}
					}
				}
			}
		}*/	
		for (int iblock=0; iblock<m_blocks.count(); iblock++) {
			if (iblock-1>=0) m_PPP_A[iblock].p_on_top_outer_interface=m_block_infos[iblock-1].p_on_bottom_inner_interface;
			if (iblock+1<m_blocks.count()) m_PPP_A[iblock].p_on_bottom_outer_interface=m_block_infos[iblock+1].p_on_top_inner_interface;
		}
		//iterate_step_A
		double r_Ap=0;
		double p_Ap=0;
		double Ap_Ap=0;
		double r_z=0;
		QList<FBBlockSolverThread *> threads;
		//define the threads
		for (int i=0; i<m_num_threads; i++) {
			FBBlockSolverThread *T0=new FBBlockSolverThread;
			threads << T0;
		}
		//add the parameters and blocks to the threads
		for (int i=0; i<m_blocks.count(); i++) {
			int thread_number=0;
			if (m_num_threads>1) thread_number=i%m_num_threads;
			threads[thread_number]->step_A_parameters << &m_PPP_A[i];
			threads[thread_number]->step_B_parameters << &m_PPP_B[i];
			threads[thread_number]->blocks << m_blocks[i];
		}
		FBTimer::stopTimer("setup_for_A");
		FBTimer::startTimer("step_A");
		//start the threads
		for (int i=0; i<m_num_threads; i++) {
			threads[i]->do_step_A=true;
			threads[i]->do_step_B=false;
			threads[i]->start();
		}
		//wait until the threads are finished
		{
			bool done=false;
			while (!done) {
				QTest::qWait(10);
				done=true;
				for (int i=0; i<m_num_threads; i++) {
					if (!threads[i]->isFinished()) {
						done=false;
					}
				}
			}
		}
		FBTimer::stopTimer("step_A");
		FBTimer::startTimer("setup_for_B");
		//define the numbers
		for (long i=0; i<m_blocks.count(); i++) {
			r_z+=m_PPP_A[i].r_z;
			r_Ap+=m_PPP_A[i].r_Ap;
			p_Ap+=m_PPP_A[i].p_Ap;
			Ap_Ap+=m_PPP_A[i].Ap_Ap;
		}		
		//iterate_step_B
		for (long i=0; i<m_blocks.count(); i++) {
			m_PPP_B[i].alpha=r_z/p_Ap; 
			if (r_z!=0) m_PPP_B[i].beta=(r_z-2*m_PPP_B[i].alpha*r_Ap+m_PPP_B[i].alpha*m_PPP_B[i].alpha*Ap_Ap)/r_z;
			else m_PPP_B[i].beta=0;
			m_PPP_B[i].WN[0]=m_bvf_map.N1();
			m_PPP_B[i].WN[1]=m_bvf_map.N2();
			m_PPP_B[i].WN[2]=m_bvf_map.N3();
		}	
		FBTimer::stopTimer("setup_for_B");
		FBTimer::startTimer("step_B");
		//start the threads
		for (int i=0; i<m_num_threads; i++) {
			threads[i]->do_step_A=false;
			threads[i]->do_step_B=true;
			threads[i]->start();
		}
		//wait until the threads are finished
		{
			bool done=false;
			while (!done) {
				QTest::qWait(10);
				done=true;
				for (int i=0; i<m_num_threads; i++) {
					if (!threads[i]->isFinished()) {
						done=false;
					}
				}
			}	
		}
		FBTimer::stopTimer("step_B");
		FBTimer::startTimer("after_B");
		qDeleteAll(threads);
		threads.clear();
		
		/*FBTimer::startTimer("update_p_on_inner_interfaces");
		//update m_p on inner interfaces
		for (long i=0; i<m_blocks.count(); i++) {		
			int x0=m_blocks[i]->xPosition();
			int y0=m_blocks[i]->yPosition();
			int z0=m_blocks[i]->zPosition();
			m_PPP_B[i].p_on_inner_interface.resetIteration();
			while (m_PPP_B[i].p_on_inner_interface.advanceIteration()) {
				int dd0=m_PPP_B[i].p_on_inner_interface.currentIndex1();
				int xx0=m_PPP_B[i].p_on_inner_interface.currentIndex2();
				int yy0=m_PPP_B[i].p_on_inner_interface.currentIndex3();
				int zz0=m_PPP_B[i].p_on_inner_interface.currentIndex4();
				float val0=m_PPP_B[i].p_on_inner_interface.currentValue();
				m_p.setValue(val0,dd0,x0-1+xx0,y0-1+yy0,z0-1+zz0);
			}				
		}
		FBTimer::stopTimer("update_p_on_inner_interfaces");*/
		for (int i=0; i<m_blocks.count(); i++) {
			m_block_infos[i].p_on_top_inner_interface=m_PPP_B[i].p_on_top_inner_interface;
			m_block_infos[i].p_on_bottom_inner_interface=m_PPP_B[i].p_on_bottom_inner_interface;
		}
		
		m_num_iterations++;
		
		FBTimer::startTimer("get_stress");
		QList<double> stress0=q->getStress();
		FBTimer::stopTimer("get_stress");
		m_error_estimator.addStressData(stress0);
		/*qDebug()  << QString("Iteration %1, Stress: (%2,%3,%4,%5,%6,%7), Est. Rel. Err.: %8").arg(m_num_iterations)
						.arg(stress0[0],0,'g',4).arg(stress0[1],0,'g',4).arg(stress0[2],0,'g',4)
						.arg(stress0[3],0,'g',4).arg(stress0[4],0,'g',4).arg(stress0[5],0,'g',4)
						.arg(m_error_estimator.estimatedRelativeError(),0,'g',4);
		*/				
		if (m_error_estimator.estimatedRelativeError()<m_epsilon) 
			num_times_below_epsilon++;
		else
			num_times_below_epsilon=0;
		FBTimer::stopTimer("after_B");
		
	}
	FBTimer::stopTimer("iterations");
}
void FBBlockSolver::clear() {
	for (long i=0; i<d->m_blocks.count(); i++) {
		d->m_blocks[i]->clearArrays();
	}
}

