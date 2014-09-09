#include "fbsparsearray4d.h"
#include <QList>
#include <QDebug>


class FBSparseArray4DPrivate {
public:
	FBSparseArray4D *q;
	long m_N1,m_N2,m_N3,m_N4;
	long m_N1N2,m_N1N2N3,m_N1N2N3N4;
	FBSparseArray1D m_data;
	long m_current_i1;
	long m_current_i2;
	long m_current_i3;
	long m_current_i4;
	//If you add more members, be sure to handle them in the copy constructor!
	
	void copy_from(const FBSparseArray4D &X);
};

void FBSparseArray4DPrivate::copy_from(const FBSparseArray4D &X) {
	m_N1=X.d->m_N1;
	m_N2=X.d->m_N2;
	m_N3=X.d->m_N3;
	m_N4=X.d->m_N4;
	m_N1N2=X.d->m_N1N2;
	m_N1N2N3=X.d->m_N1N2N3;
	m_N1N2N3N4=X.d->m_N1N2N3N4;
	m_data=X.d->m_data;
}

FBSparseArray4D::FBSparseArray4D() 
{
	d=new FBSparseArray4DPrivate;
	d->q=this;
	d->m_N1=d->m_N2=d->m_N3=d->m_N4=0;
	d->m_N1N2=d->m_N1N2N3=d->m_N1N2N3N4=0;
}
FBSparseArray4D::FBSparseArray4D(const FBSparseArray4D &X) {
	d=new FBSparseArray4DPrivate;
	d->q=this;
	d->m_N1=d->m_N2=d->m_N3=d->m_N4=0;
	d->m_N1N2=d->m_N1N2N3=d->m_N1N2N3N4=0;
	d->copy_from(X);	
}

void FBSparseArray4D::operator=(const FBSparseArray4D &X) {
	d->copy_from(X);
}

FBSparseArray4D::~FBSparseArray4D()
{
	delete d;
}
void FBSparseArray4D::allocate(int data_type,int N1,int N2,int N3,int N4) {
	d->m_data.allocate(data_type,N1*N2*N3*N4);
	d->m_N1=N1; d->m_N2=N2; d->m_N3=N3; d->m_N4=N4;
	d->m_N1N2=N1*N2; d->m_N1N2N3=N1*N2*N3; d->m_N1N2N3N4=N1*N2*N3*N4;
}
void FBSparseArray4D::setValue(float val,long i1,long i2,long i3,long i4) {
	if ((i1<0)||(i1>=d->m_N1)) return;
	if ((i2<0)||(i2>=d->m_N2)) return;
	if ((i3<0)||(i3>=d->m_N3)) return;
	if ((i4<0)||(i4>=d->m_N4)) return;
	d->m_data.setValue(val,i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4);
}
float FBSparseArray4D::value(long i1,long i2,long i3,long i4) const {
	if ((i1<0)||(i1>=d->m_N1)) return 0;
	if ((i2<0)||(i2>=d->m_N2)) return 0;
	if ((i3<0)||(i3>=d->m_N3)) return 0;
	if ((i4<0)||(i4>=d->m_N4)) return 0;
	return d->m_data.value(i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4);
}
float FBSparseArray4D::value1(long i) const {
	return d->m_data.value(i);
}
int FBSparseArray4D::N1() const {
	return d->m_N1;
}
int FBSparseArray4D::N2() const {
	return d->m_N2;
}
int FBSparseArray4D::N3() const {
	return d->m_N3;
}
int FBSparseArray4D::N4() const {
	return d->m_N4;
}
void FBSparseArray4D::resetIteration() const {
	d->m_current_i1=d->m_current_i2=d->m_current_i3=d->m_current_i4=-1;
	d->m_data.resetIteration();
}
bool FBSparseArray4D::advanceIteration() const {
	if (!d->m_data.entryCount()) return false;
	if (!d->m_N1N2N3N4) return false;
	long oldind=d->m_data.currentIndex();
	if (!d->m_data.advanceIteration()) return false;
	long newind=d->m_data.currentIndex();
	if (oldind<0) {
		d->m_current_i1=newind%d->m_N1;
		d->m_current_i2=(newind/d->m_N1)%d->m_N2;
		d->m_current_i3=(newind/d->m_N1N2)%d->m_N3;
		d->m_current_i4=(newind/d->m_N1N2N3)%d->m_N4;
		return true;
	}
	else {
		long diffind=newind-oldind;
		d->m_current_i1+=diffind;
		if (d->m_current_i1>=d->m_N1) {
			while (d->m_current_i1>=d->m_N1) {
				d->m_current_i1-=d->m_N1;
				d->m_current_i2++;
			}
		}
		else return true;
	
		if (d->m_current_i2>=d->m_N2) {
			while (d->m_current_i2>=d->m_N2) {
				d->m_current_i2-=d->m_N2;
				d->m_current_i3++;
			}
		}
		else return true;
	
		if (d->m_current_i3>=d->m_N3) {
			while (d->m_current_i3>=d->m_N3) {
				d->m_current_i3-=d->m_N3;
				d->m_current_i4++;
			}
		}
		else return true;
	
		if (d->m_current_i4>=d->m_N4) {
			qWarning() << "Unexpected problem in FBSparseArray4D::advanceIteration().";
			return false;
		}
		else return true;
	}
}
int FBSparseArray4D::currentIndex1() const {
	return d->m_current_i1;
}
int FBSparseArray4D::currentIndex2() const {
	return d->m_current_i2;
}
int FBSparseArray4D::currentIndex3() const {
	return d->m_current_i3;
}
int FBSparseArray4D::currentIndex4() const {
	return d->m_current_i4;
}
float FBSparseArray4D::currentValue() const {
	return d->m_data.currentValue();
}
bool FBSparseArray4D::setupIndex(int pass,long i1,long i2,long i3,long i4) {
	return d->m_data.setupIndex(pass,i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4);
}







