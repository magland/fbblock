#ifndef fbsparsearray4d_H
#define fbsparsearray4d_H

#include "fbsparsearray1d.h"

class FBSparseArray4DPrivate;
class FBSparseArray4D {
public:
	friend class FBSparseArray4DPrivate;
	FBSparseArray4D();
	FBSparseArray4D(const FBSparseArray4D &X);
	void operator=(const FBSparseArray4D &X);
	virtual ~FBSparseArray4D();
	void allocate(int data_type,int N1,int N2,int N3,int N4);
	bool setupIndex(int pass,long i1,long i2,long i3,long i4); //pass=1,2
	void setValue(float val,long i1,long i2,long i3,long i4);
	float value(long i1,long i2,long i3,long i4) const;
	float value1(long i) const;
	int N1() const;
	int N2() const;
	int N3() const;
	int N4() const;
	void resetIteration() const;
	bool advanceIteration() const;
	int currentIndex1() const;
	int currentIndex2() const;
	int currentIndex3() const;
	int currentIndex4() const;
	float currentValue() const;
private:
	FBSparseArray4DPrivate *d;
};

#endif
