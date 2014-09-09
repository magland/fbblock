#ifndef arrays_H
#define arrays_H

#include "fbsparsearray4d.h"
#include <malloc.h>
#include <QVector>
#include <QDebug>

template <class TT> class FBArray1DPrivate;
template <class T>
class FBArray1D {
public:
	template <class TT> friend class FBArray1DPrivate;
	FBArray1D();
	FBArray1D(const FBArray1D &X);
	void operator=(const FBArray1D &X);
	virtual ~FBArray1D();
	void allocate(long N);
	void setAll(T val);	
	long length() const;
	void clear();
	
	T *ptr;
private:
	FBArray1DPrivate <T> *d;
};

template <class TT>
class FBArray1DPrivate {
public:
	FBArray1D<TT> *q;
	long m_N;
	void copy_from(const FBArray1D<TT> &X) {
		if (X.d->m_N!=m_N) q->allocate(X.d->m_N);
		for (long i=0; i<m_N; i++)
			q->ptr[i]=X.d->q->ptr[i];
	}
};

template <class TT> class FBArray4DPrivate;
template <class T>
class FBArray4D {
public:
	template <class TT> friend class FBArray4DPrivate;
	FBArray4D();
	FBArray4D(const FBArray4D &X);
	void operator=(const FBArray4D &X);
	virtual ~FBArray4D();
	void allocate(long N1,long N2,long N3,long N4);	
	T value(long i1,long i2,long i3,long i4) const ;
	T value1(long i1) const ;
	void setValue(T val,long i1,long i2,long i3,long i4);
	void setValue1(T val,long i);
	void incrementValue(T val,long i1,long i2,long i3,long i4);
	void setAll(T val);
	long N1() const;
	long N2() const;
	long N3() const;
	long N4() const;
	long debug_length_of_data() const;
	void clear();
private: 
	FBArray4DPrivate <T> *d;
};

template <class T>
class FBArray2D : public FBArray4D<T> {
public:
	FBArray2D() {}
	FBArray2D<T>(const FBArray2D &X) : FBArray4D<T>(X) {};
	void allocate(long N1,long N2) {FBArray4D<T>::allocate(N1,N2,1,1);}
	T value(long i1,long i2) const {return FBArray4D<T>::value(i1,i2,0,0);}
	void setValue(T val,long i1,long i2) {FBArray4D<T>::setValue(val,i1,i2,0,0);}
};
template <class T>
class FBArray3D : public FBArray4D<T> {
public:
	FBArray3D() {};
	FBArray3D<T>(const FBArray3D &X) : FBArray4D<T>(X) {};
	void allocate(long N1,long N2,long N3) {FBArray4D<T>::allocate(N1,N2,N3,1);}
	T value(long i1,long i2,long i3) const {return FBArray4D<T>::value(i1,i2,i3,0);}
	void setValue(T val,long i1,long i2,long i3) {FBArray4D<T>::setValue(val,i1,i2,i3,0);}
};

template <class T> FBArray1D<T>::FBArray1D() 
{
	d=new FBArray1DPrivate<T>;
	d->q=this;
	ptr=0;
	d->m_N=0;
}
template <class T> FBArray1D<T>::FBArray1D(const FBArray1D &X) {
	d=new FBArray1DPrivate<T>;
	d->q=this;
	ptr=0;
	d->m_N=0;
	d->copy_from(X);
}
template <class T> void FBArray1D<T>::operator=(const FBArray1D &X) {
	d->copy_from(X);
}

template <class T> FBArray1D<T>::~FBArray1D()
{
	if (ptr) free(ptr);
	delete d;
}
template <class T> void FBArray1D<T>::allocate(long N) {
	d->m_N=0;
	if (ptr) free(ptr);
	if (N==0) {
		qWarning() << "Allocating array of size zero";
		ptr=0; return;
	}	
	ptr=(T *)malloc(sizeof(T)*N);
	if (ptr) {
		d->m_N=N;
		setAll(0);
	}
	else {
		ptr=0;
		qWarning() << "Error allocating!";
	}
}
template <class T> void FBArray1D<T>::setAll(T val) {
	for (long i=0; i<d->m_N; i++) ptr[i]=val;
}
template <class T> long FBArray1D<T>::length() const {
	return d->m_N;
}
template <class T> void FBArray1D<T>::clear() {
	if (ptr) free(ptr); ptr=0;
}
template <class TT>
class FBArray4DPrivate {
public:
	FBArray4D<TT> *q;
	long m_N1,m_N2,m_N3,m_N4;
	long m_N1N2;
	long m_N1N2N3;
	FBArray1D<TT> m_data;
	void copy_from(const FBArray4D<TT> &X) {
		m_N1=X.d->m_N1;
		m_N2=X.d->m_N2;
		m_N3=X.d->m_N3;
		m_N4=X.d->m_N4;
		m_N1N2=X.d->m_N1N2;
		m_N1N2N3=X.d->m_N1N2N3;
		m_data=X.d->m_data;
	}
};

template <class T> FBArray4D<T>::FBArray4D() 
{
	d=new FBArray4DPrivate<T>;
	d->q=this;
	d->m_N1=d->m_N2=d->m_N3=d->m_N4=0;
	d->m_N1N2=0;
	d->m_N1N2N3=0;
}
template <class T> FBArray4D<T>::FBArray4D(const FBArray4D &X) {
	d=new FBArray4DPrivate<T>;
	d->q=this;
	d->m_N1=d->m_N2=d->m_N3=d->m_N4=0;
	d->m_N1N2=0;
	d->m_N1N2N3=0;
	d->copy_from(X);
}
template <class T> void FBArray4D<T>::operator=(const FBArray4D &X) {
	d->copy_from(X);
}

template <class T> FBArray4D<T>::~FBArray4D()
{
	delete d;
}
template <class T> void FBArray4D<T>::allocate(long N1,long N2,long N3,long N4) {
	d->m_N1=N1;
	d->m_N2=N2;
	d->m_N3=N3;
	d->m_N4=N4;
	d->m_N1N2=d->m_N1*d->m_N2;
	d->m_N1N2N3=d->m_N1*d->m_N2*d->m_N3;
	d->m_data.allocate(N1*N2*N3*N4);
}
template <class T> T FBArray4D<T>::value(long i1,long i2,long i3,long i4) const {
	if ((i1<0)||(i1>=d->m_N1)) return 0;
	if ((i2<0)||(i2>=d->m_N2)) return 0;
	if ((i3<0)||(i3>=d->m_N3)) return 0;
	if ((i4<0)||(i4>=d->m_N4)) return 0;
	T ret=d->m_data.ptr[i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4];
	return ret;
}
template <class T> T FBArray4D<T>::value1(long i) const {
	T ret=d->m_data.ptr[i];
	return ret;
}
template <class T> void FBArray4D<T>::setValue(T val,long i1,long i2,long i3,long i4) {
	if ((i1<0)||(i1>=d->m_N1)) return;
	if ((i2<0)||(i2>=d->m_N2)) return;
	if ((i3<0)||(i3>=d->m_N3)) return;
	if ((i4<0)||(i4>=d->m_N4)) return;
	d->m_data.ptr[i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4]=val;
}
template <class T> void FBArray4D<T>::incrementValue(T val,long i1,long i2,long i3,long i4) {
	if ((i1<0)||(i1>=d->m_N1)) return;
	if ((i2<0)||(i2>=d->m_N2)) return;
	if ((i3<0)||(i3>=d->m_N3)) return;
	if ((i4<0)||(i4>=d->m_N4)) return;
	d->m_data.ptr[i1+d->m_N1*i2+d->m_N1N2*i3+d->m_N1N2N3*i4]+=val;
}
template <class T> void FBArray4D<T>::setValue1(T val,long i) {
	d->m_data.ptr[i]=val;
}
template <class T> void FBArray4D<T>::setAll(T val) {
	d->m_data.setAll(val);
}
template <class T> long FBArray4D<T>::N1() const {
	return d->m_N1;
}
template <class T> long FBArray4D<T>::N2() const {
	return d->m_N2;
}
template <class T> long FBArray4D<T>::N3() const {
	return d->m_N3;
}
template <class T> long FBArray4D<T>::N4() const {
	return d->m_N4;
}
template <class T> void FBArray4D<T>::clear() 
{
	d->m_data.clear();
	d->m_N1=d->m_N2=d->m_N3=d->m_N4=0;
	d->m_N1N2=0;
	d->m_N1N2N3=0;
}
template <class T> long FBArray4D<T>::debug_length_of_data() const {
	return d->m_data.length();
}



#endif


