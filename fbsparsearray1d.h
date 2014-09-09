#ifndef fbsparsearray1d_H
#define fbsparsearray1d_H

#define DATA_TYPE_FLOAT 1
#define DATA_TYPE_BYTE 2

class FBSparseArray1DPrivate;
class FBSparseArray1D {
public:
	friend class FBSparseArray1DPrivate;
	FBSparseArray1D();
	FBSparseArray1D(const FBSparseArray1D &X);
	void operator=(const FBSparseArray1D &X);
	virtual ~FBSparseArray1D();
	void allocate(int data_type,long N);
	bool setupIndex(int pass,long ind); //pass=1,2
	void clear();
	long length() const;
	long entryCount() const;
	void setValue(float val,long ind);
	float value(long ind) const;
	void resetIteration() const;
	bool advanceIteration() const;
	long currentIndex() const;
	float currentValue() const;
private:
	FBSparseArray1DPrivate *d;
};

#endif
