#ifndef mda_H
#define mda_H

#ifndef MDA_EXPORT
#define MDA_EXPORT
#endif

#include <stdlib.h>
#include <stdio.h>
#include <QHash>
#include <QString>

#define MAX_MDA_DIMS 50
#define MDA_COMMAND_MASK -3
#define MDA_TYPE_COMPLEX -1
#define MDA_TYPE_BYTE -2
#define MDA_TYPE_REAL -3
#define MDA_TYPE_SHORT -4
#define MDA_TYPE_INT32 -5
#define MDA_TYPE_UINT16 -6

typedef float real;
class MDA_EXPORT AbstractComplex {
public:
	virtual real re() const {return 0;}
	virtual real im() const {return 0;}	
};
class MDA_EXPORT MDAComplex : public AbstractComplex {
public:
	MDAComplex(real rp=0,real ip=0) {m_re=rp; m_im=ip;}
	MDAComplex(const AbstractComplex &X) {m_re=X.re(); m_im=X.im();}
	virtual ~MDAComplex() {};
	
	real re() const {return m_re;}
	real im() const {return m_im;}
	void setRe(real val) {m_re=val;}
	void setIm(real val) {m_im=val;}
	
	void operator=(const AbstractComplex &X) {setRe(X.re()); setIm(X.im());}
	void operator=(real X) {setRe(X); setIm(0);}
	
private:
	real m_re,m_im;
};
//operatations on complex numbers
/*! Norm squared*/
real MDA_EXPORT normsqr(const AbstractComplex &X);
/*! Complex product.*/
MDAComplex MDA_EXPORT operator*(const AbstractComplex &X1, const AbstractComplex &X2);
/*! Complex quotient.*/
MDAComplex MDA_EXPORT operator/(const AbstractComplex &X1, const AbstractComplex &X2);
/*! Complex sum.*/
MDAComplex MDA_EXPORT operator+(const AbstractComplex &X1, const AbstractComplex &X2);
/*! Complex difference.*/
MDAComplex MDA_EXPORT operator-(const AbstractComplex &X1, const AbstractComplex &X2);
/*! Product of complex and real.*/
MDAComplex MDA_EXPORT operator*(const AbstractComplex &X1, real X2);
/*! Quotient of complex and real.*/
MDAComplex MDA_EXPORT operator/(const AbstractComplex &X1, real X2);
/*! Sum of complex and real.*/
MDAComplex MDA_EXPORT operator+(const AbstractComplex &X1, real X2);
/*! Difference of complex and real.*/
MDAComplex MDA_EXPORT operator-(const AbstractComplex &X1, real X2);
/*! Product of real and complex.*/
MDAComplex MDA_EXPORT operator*(real X1, const AbstractComplex &X2);
/*! Quotient of real and complex.*/
MDAComplex MDA_EXPORT operator/(real X1, const AbstractComplex &X2);
/*! Sum of real and complex.*/
MDAComplex MDA_EXPORT operator+(real X1, const AbstractComplex &X2);
/*! Difference of real and complex.*/
MDAComplex MDA_EXPORT operator-(real X1, const AbstractComplex &X2);
/*! Complex comparison.*/
bool MDA_EXPORT operator==(const AbstractComplex &X1, const AbstractComplex &X2);
/*! Complex/real comparison.*/
bool MDA_EXPORT operator==(const AbstractComplex &X1, real X2);
/*! Real/complex comparison.*/
bool MDA_EXPORT operator==(real X1, const AbstractComplex &X2);
/*! Absolute value of a complex number*/
real MDA_EXPORT abs(const AbstractComplex &X);
/*! Complex conjugate*/
MDAComplex MDA_EXPORT conj(const AbstractComplex &X);
/*! Complex exponential function*/
MDAComplex MDA_EXPORT exp(const AbstractComplex &X);
/*! Complex logarithmic function*/
MDAComplex MDA_EXPORT log(const AbstractComplex &X);
/*! Complex square root function*/
MDAComplex MDA_EXPORT sqrt(const AbstractComplex &X);
/*! Complex power*/
MDAComplex MDA_EXPORT power(const AbstractComplex &X1, const AbstractComplex &X2);
MDAComplex MDA_EXPORT operator-(const AbstractComplex &X);




typedef int integer;

/*#ifdef __LP64__
	typedef int int32;
#else
	typedef qint32 int32;
#endif*/

class MDA_EXPORT MdaElement : public AbstractComplex {
public:
	void *m_data;
	int m_data_type;
	
	real re() const;
	real im() const;
	void setRe(real val);
	void setIm(real val);
	
	void operator=(const MdaElement &X) {setRe(X.re()); setIm(X.im());}
	void operator=(const MDAComplex &X) {setRe(X.re()); setIm(X.im());}
	void operator=(real val) {setRe(val); setIm(0);}

	
};

class Mda;
//these commented out on 12/14/12
//typedef Mda Vector;
//typedef Mda Matrix;
//typedef Mda RealMda;
//typedef real NonnegativeReal;

struct complex_struct {
	real re,im;
};


/*!
	\brief Multi-dimensional array - basic data type in ChainLink
*/
class MDA_EXPORT Mda {
	friend class MdaElement;
public:
	//true if this object represents a ":" vector, i.e. all indices
	bool all_indices;
	Mda();
	Mda(const Mda &X);
	virtual ~Mda();
	///set equal to a 1x1 array consisting of one real number
	void operator=(real val);
	///set equal to a 1x1 array consisting of one complex number
	void operator=(const AbstractComplex &val);
	///set equal to Mda X
	void operator=(const Mda &X);
	///allocate with datatype=MDA_TYPE_REAL,MDA_TYPE_COMPLEX,MDA_TYPE_BYTE,MDA_TYPE_INT32,MDA_TYPE_SHORT,MDA_TYPE_BYTE, of dimensions S1xS2x...
	void allocate(int datatype,qint32 S1,qint32 S2=-1,qint32 S3=-1,qint32 S4=-1,qint32 S5=-1,qint32 S6=-1,qint32 S7=-1,qint32 S8=-1);
	///allocate with datatype=MDA_TYPE_REAL,MDA_TYPE_COMPLEX,MDA_TYPE_BYTE,MDA_TYPE_INT32,MDA_TYPE_SHORT,MDA_TYPE_BYTE, of dimensions dims
	void allocate(int datatype,int num_dimensions,qint32 *dims);
	///clear and deallocate, if allocate_one_by_one=true then allocates a 1x1 array of 0
	void clear(bool allocate_one_by_one=true);
	///set all entries equal to val
	void setAll(const AbstractComplex &val);
	///set all entries equal to val
	void setAll(real val);
	///returns true if (ind1,ind2,...) is in the range of this array, i.e. is valid index
	bool inRange(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1);
	///access element at given 1d index... the element can either be set or used in an expression
	MdaElement operator[](qint64 index);
	///access element at given index... the element can either be set or used in an expression
	MdaElement operator()(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1);
	///access element at given index... the element can either be set or used in an expression
	MdaElement operator()(qint32 *index);
	///returns the complex number at given index
	MDAComplex get(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1) const;
	///returns the complex number at given index
	MDAComplex get(qint32 *index) const;
	///returns reference to the mask character at given index
	unsigned char &mask(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1);
	///returns reference to the mask character at given index
	unsigned char &mask(qint32 *index);
	///returns the mask character at given index
	unsigned char getmask(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1) const;
	///returns the mask character at given index
	unsigned char getmask(qint32 *index) const;
	///returns the number of dimensions in the array
	int dimCount() const;
	qint32 *dims() {return m_dims;}
	qint64 size() const;
	qint32 size(int dim) const;
	int dataType() const;
	void convertToComplex();
	void convertToReal();
	void convertToByte();
	
	bool read(const QString &filename);
	bool write(const QString &filename) const;
	bool readSlices(const QString &filename,long minslice,long maxslice);
	static bool readSize(const QString &filename,Mda &ret);
	
	///Shift with wrapping
	void shift(qint32 *shifts);
	void squeeze();
	void reshape(int dim_count,qint32 *dims);
	void reshape(qint32 S1,qint32 S2=-1,qint32 S3=-1,qint32 S4=-1,qint32 S5=-1,qint32 S6=-1,qint32 S7=-1,qint32 S8=-1);
public:
	/*! Pointer to the actual data if dataType()==MDA_TYPE_BYTE*/
	unsigned char *data_byte;
	/*! Pointer to the actual data if dataType()==MDA_TYPE_SHORT*/
	short int *data_short;
	/*! Pointer to the actual data if dataType()==MDA_TYPE_UINT16*/
	quint16 *data_uint16;
	/*! Pointer to the actual data if dataType()==MDA_TYPE_INT32*/
	qint32 *data_int32;
	/*! Pointer to the actual data if dataType()==MDA_TYPE_REAL*/
	real *data_real;
	/*! Pointer to the actual data if dataType()==MDA_TYPE_COMPLEX*/
	complex_struct *data_complex;
	/*! Pointer to mask data. The mask is simply an array of unsigned characters with the same dimensions as the mda.*/
	unsigned char *data_mask;
private:
	int m_dim_count;
	qint32 m_dims[MAX_MDA_DIMS];
	qint64 m_dimprod;
	int m_data_type;
	
	void initialize(int datatype);
	void allocate_mask();
	qint64 get_index(qint64 ind1,qint32 ind2=-1,qint32 ind3=-1,qint32 ind4=-1,qint32 ind5=-1,qint32 ind6=-1,qint32 ind7=-1,qint32 ind8=-1) const;
	qint64 get_index(qint32 *index) const;
	void copy_from(const Mda &X);
};

#endif
