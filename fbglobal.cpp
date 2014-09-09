#include "fbglobal.h"

struct PolynomialTerm {
	fbreal coefficient;
	QList<int> exponents; // x, (1-x), y, (1-y), z, (1-z)
};

class Polynomial {
public:
	Polynomial() {}
	Polynomial(const Polynomial &P) {copy_from(P);}
	virtual ~Polynomial() {}
	int termCount() const {return m_terms.count();}
	PolynomialTerm term(int ind) const {return m_terms[ind];}
	void addTerm(PolynomialTerm T) {m_terms << T;}
	void addTerm(fbreal coeff,int e1,int e2,int e3,int e4,int e5,int e6);
	Polynomial operator*(const Polynomial &P);
	void operator=(const Polynomial &P) {copy_from(P);}
	fbreal integral(fbreal a,fbreal b,fbreal c);
private:
	QList<PolynomialTerm> m_terms;
	
	void copy_from(const Polynomial &P) {
		m_terms.clear();
		for (int j=0; j<P.termCount(); j++) {
			addTerm(P.term(j));
		}
	}
};


Polynomial Polynomial::operator*(const Polynomial &P) {
	Polynomial ret;
	for (int i=0; i<termCount(); i++)
		for (int j=0; j<P.termCount(); j++) {
			PolynomialTerm T;
			T.exponents.clear();
			for (int k=0; k<6; k++) {
				T.exponents << term(i).exponents[k]+P.term(i).exponents[k];
			}
			T.coefficient=term(i).coefficient*P.term(i).coefficient;
			ret.addTerm(T);
		}
	return ret;
}

fbreal get_factor(int i1,int i2,fbreal res) {
	if ((i1==0)&&(i2==0)) return 1.0F/res;
	if ((i1==0)&&(i2==1)) return 1.0F/2;
	if ((i1==1)&&(i2==0)) return 1.0F/2;
	if ((i1==1)&&(i2==1)) return 1.0F/6*res;
	if ((i1==2)&&(i2==0)) return 1.0F/3*res;
	if ((i1==0)&&(i2==2)) return 1.0F/3*res;
	return 0;
}

fbreal term_integral(PolynomialTerm T,fbreal a,fbreal b,fbreal c) {
	return T.coefficient
			*get_factor(T.exponents[0],T.exponents[1],a)
			*get_factor(T.exponents[2],T.exponents[3],b)
			*get_factor(T.exponents[4],T.exponents[5],c);
}

fbreal Polynomial::integral(fbreal a,fbreal b,fbreal c) {
	fbreal ret=0;
	for (int j=0; j<termCount(); j++) {
		ret+=term_integral(term(j),a,b,c);
	}
	return ret;
}

void Polynomial::addTerm(fbreal coeff,int e1,int e2,int e3,int e4,int e5,int e6) {
	PolynomialTerm T;
	T.coefficient=coeff;
	T.exponents << e1 << e2 << e3 << e4 << e5 << e6;
	addTerm(T);
}

Polynomial get_eps11_term(int a1,int a2,int a3,int dir) {
	Polynomial ret;
	if (dir!=0) return ret;
	PolynomialTerm T;
	if (a1==0) T.coefficient=-1; else T.coefficient=1;
	for (int j=0; j<6; j++) T.exponents << 0;
	if (a2==0) T.exponents[3]=1; else T.exponents[2]=1;
	if (a3==0) T.exponents[5]=1; else T.exponents[4]=1;
	ret.addTerm(T);
	return ret;
}

Polynomial get_eps22_term(int a1,int a2,int a3,int dir) {
	Polynomial ret;
	if (dir!=1) return ret;
	PolynomialTerm T;
	if (a2==0) T.coefficient=-1; else T.coefficient=1;
	for (int j=0; j<6; j++) T.exponents << 0;
	if (a1==0) T.exponents[1]=1; else T.exponents[0]=1;
	if (a3==0) T.exponents[5]=1; else T.exponents[4]=1;
	ret.addTerm(T);
	return ret;
}

Polynomial get_eps33_term(int a1,int a2,int a3,int dir) {
	Polynomial ret;
	if (dir!=2) return ret;
	PolynomialTerm T;
	if (a3==0) T.coefficient=-1; else T.coefficient=1;
	for (int j=0; j<6; j++) T.exponents << 0;
	if (a1==0) T.exponents[1]=1; else T.exponents[0]=1;
	if (a2==0) T.exponents[3]=1; else T.exponents[2]=1;
	ret.addTerm(T);
	return ret;
}

Polynomial get_eps12_term(int a1,int a2,int a3,int dir) {
	if (dir==0) return get_eps22_term(a1,a2,a3,1);
	else if (dir==1) return get_eps11_term(a1,a2,a3,0);
	Polynomial ret;
	return ret;
}

Polynomial get_eps13_term(int a1,int a2,int a3,int dir) {
	if (dir==0) return get_eps33_term(a1,a2,a3,2);
	else if (dir==2) return get_eps11_term(a1,a2,a3,0);
	Polynomial ret;
	return ret;
}

Polynomial get_eps23_term(int a1,int a2,int a3,int dir) {
	if (dir==1) return get_eps33_term(a1,a2,a3,2);
	else if (dir==2) return get_eps22_term(a1,a2,a3,1);
	Polynomial ret;
	return ret;
}

void compute_stiffness_kernel(FBArray2D<float> &stiffness_matrix,const fbreal youngs_modulus,const fbreal poissons_ratio,const QList<fbreal> &resolution) {
	fbreal factor=youngs_modulus/((1+poissons_ratio)*(1-2*poissons_ratio));
	fbreal stiffness_tensor[6][6];
	for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
		stiffness_tensor[i][j]=0;
	stiffness_tensor[0][0]=stiffness_tensor[1][1]=stiffness_tensor[2][2]=factor*(1-poissons_ratio);
	stiffness_tensor[0][1]=stiffness_tensor[0][2]=stiffness_tensor[1][0]=stiffness_tensor[1][2]=stiffness_tensor[2][0]=stiffness_tensor[2][1]=factor*poissons_ratio;
	stiffness_tensor[3][3]=stiffness_tensor[4][4]=stiffness_tensor[5][5]=factor*(1-2*poissons_ratio)/2;

	for (int i=0; i<24; i++)
		for (int j=0; j<24; j++) {
			stiffness_matrix.setValue(0,i,j);
		}
	for (int a1=0; a1<=1; a1++)
	for (int a2=0; a2<=1; a2++)
	for (int a3=0; a3<=1; a3++)
	for (int adir=0; adir<3; adir++) {
		int aind=adir+a1*3+a2*6+a3*12;
		for (int b1=0; b1<=1; b1++)
		for (int b2=0; b2<=1; b2++)
		for (int b3=0; b3<=1; b3++)
		for (int bdir=0; bdir<3; bdir++) {
			int bind=bdir+b1*3+b2*6+b3*12;	
			Polynomial deps_da[6]; //derivative of strain w.r.t. a
			Polynomial deps_db[6]; //derivative of strain w.r.t. b
			
			deps_da[0]=get_eps11_term(a1,a2,a3,adir);
			deps_db[0]=get_eps11_term(b1,b2,b3,bdir);
			deps_da[1]=get_eps22_term(a1,a2,a3,adir);
			deps_db[1]=get_eps22_term(b1,b2,b3,bdir);
			deps_da[2]=get_eps33_term(a1,a2,a3,adir);
			deps_db[2]=get_eps33_term(b1,b2,b3,bdir);
			deps_da[3]=get_eps12_term(a1,a2,a3,adir);
			deps_db[3]=get_eps12_term(b1,b2,b3,bdir);
			deps_da[4]=get_eps13_term(a1,a2,a3,adir);
			deps_db[4]=get_eps13_term(b1,b2,b3,bdir);
			deps_da[5]=get_eps23_term(a1,a2,a3,adir);
			deps_db[5]=get_eps23_term(b1,b2,b3,bdir);

			fbreal coeff=0;			
			for (int i=0; i<6; i++)
			for (int j=0; j<6; j++) {
				Polynomial P;
				P=deps_da[i]*deps_db[j];
				coeff += -P.integral(resolution[0],resolution[1],resolution[2])*stiffness_tensor[i][j];
			}
			
			stiffness_matrix.setValue(stiffness_matrix.value(aind,bind)+coeff,aind,bind);
		}
	}
}
bool is_element(FBArray3D<unsigned char> &bvfmap,long i1,long i2,long i3) {
	if ((i1<0)||(i1>=bvfmap.N1())) return false;
	if ((i2<0)||(i2>=bvfmap.N2())) return false;
	if ((i3<0)||(i3>=bvfmap.N3())) return false;
	return (bvfmap.value(i1,i2,i3)!=0);
}
bool is_vertex(FBArray3D<unsigned char> &bvfmap,long i1,long i2,long i3) {
	for (int d1=-1; d1<=0; d1++)
	for (int d2=-1; d2<=0; d2++)
	for (int d3=-1; d3<=0; d3++)
		if (is_element(bvfmap,i1+d1,i2+d2,i3+d3))
			return true;
	return false;
}

fbreal initial_displacement(int i1,int i2,int i3,int dir,fbreal resolution[3],FBMacroscopicStrain &macroscopic_strain) {
	fbreal x=i1*resolution[0];
	fbreal y=i2*resolution[1];
	fbreal z=i3*resolution[2];
	if (dir==0) return macroscopic_strain.eps11*x+macroscopic_strain.eps12/2*y+macroscopic_strain.eps13/2*z;
	if (dir==1) return macroscopic_strain.eps12/2*x+macroscopic_strain.eps22*y+macroscopic_strain.eps23/2*z;
	if (dir==2) return macroscopic_strain.eps13/2*x+macroscopic_strain.eps23/2*y+macroscopic_strain.eps33*z;
	return 0;
}
