#include "fbsparsearray1d.h"
#include <QDebug>

#define SPARSE_STAGE_EMPTY 0
#define SPARSE_STAGE_ALLOCATED 1
#define SPARSE_STAGE_STEP1 2
#define SPARSE_STAGE_STEP2 3
#define SPARSE_STAGE_LOCKED 4

//note that this has to be strictly less than 256 because of the m_block_entry_counts being between 0 and 255
#define SPARSE_BLOCK_SIZE 255

class FBSparseArray1DPrivate {
public:
	FBSparseArray1D *q;
	long m_N;
	long m_block_count;
	int m_data_type;
	int m_stage;
	long m_entry_count;
	float *m_data_float;
	unsigned char *m_data_byte;
	unsigned char *m_offsets;
	long *m_block_indices;
	unsigned char *m_block_entry_counts;
	long m_last_step2_index;
	long m_last_step2_index_in_block;
	bool m_has_error;
	long m_current_block_index;
	int m_current_index_in_block;
	
	void finalize_step_1();
	void finalize_step_2();
	int find_offset_index(long blocknum,unsigned char offset);
	void copy_from(const FBSparseArray1D &X);
};
void FBSparseArray1DPrivate::copy_from(const FBSparseArray1D &X) {
	q->allocate(X.d->m_data_type,X.d->m_N);
	
	for (int pass=1; pass<=2; pass++) {
		X.resetIteration();
		while (X.advanceIteration()) {
			q->setupIndex(pass,X.currentIndex());
		}
	}
	
	X.resetIteration();
	while (X.advanceIteration()) {
		q->setValue(X.currentValue(),X.currentIndex());
	}
}

FBSparseArray1D::FBSparseArray1D() 
{
	d=new FBSparseArray1DPrivate;
	d->q=this;
	d->m_N=0;
	d->m_data_type=DATA_TYPE_FLOAT;
	d->m_stage=SPARSE_STAGE_EMPTY;
	d->m_data_byte=0;
	d->m_data_float=0;
	d->m_offsets=0;
	d->m_block_indices=0;
	d->m_block_entry_counts=0;
	d->m_entry_count=0;
	d->m_block_count=0;
	d->m_has_error=false;
	d->m_current_block_index=0;
	d->m_current_index_in_block=0;
}
FBSparseArray1D::FBSparseArray1D(const FBSparseArray1D &X) {
	d=new FBSparseArray1DPrivate;
	d->q=this;
	d->m_N=0;
	d->m_data_type=DATA_TYPE_FLOAT;
	d->m_stage=SPARSE_STAGE_EMPTY;
	d->m_data_byte=0;
	d->m_data_float=0;
	d->m_offsets=0;
	d->m_block_indices=0;
	d->m_block_entry_counts=0;
	d->m_entry_count=0;
	d->m_block_count=0;
	d->m_has_error=false;
	d->m_current_block_index=0;
	d->m_current_index_in_block=0;
	d->copy_from(X);
}
FBSparseArray1D::~FBSparseArray1D()
{
	clear();
	delete d;
}
void FBSparseArray1D::clear() {
	if (d->m_data_byte) free(d->m_data_byte); d->m_data_byte=0;
	if (d->m_data_float) free(d->m_data_float); d->m_data_float=0;
	if (d->m_offsets) free(d->m_offsets); d->m_offsets=0;
	if (d->m_block_indices) free(d->m_block_indices); d->m_block_indices=0;
	if (d->m_block_entry_counts) free(d->m_block_entry_counts); d->m_block_entry_counts=0;
	d->m_stage=SPARSE_STAGE_EMPTY;
	d->m_data_type=DATA_TYPE_FLOAT;
	d->m_N=0;
	d->m_block_count=0;
	d->m_entry_count=0;
	d->m_has_error=false;
}
void FBSparseArray1D::allocate(int data_type,long N) {
	clear();
	d->m_N=N;
	d->m_block_count=N/SPARSE_BLOCK_SIZE+1;
	d->m_data_type=data_type;
	d->m_block_indices=(long *)malloc(sizeof(long)*d->m_block_count);
	d->m_block_entry_counts=(unsigned char *)malloc(sizeof(unsigned char)*d->m_block_count);
	for (long i=0; i<d->m_block_count; i++) {
		d->m_block_indices[i]=0;
		d->m_block_entry_counts[i]=0;
	}
	d->m_stage=SPARSE_STAGE_ALLOCATED;
	
}
bool FBSparseArray1D::setupIndex(int pass,long ind) {
	if (d->m_has_error) return false;
	if ((ind>=d->m_N)||(ind<0)) {
		qWarning() << "problem in setupIndex: index out of range" << ind << d->m_N;
		return false;
	}
	if (pass==1) {	
		if (d->m_stage!=SPARSE_STAGE_ALLOCATED) {
			qWarning() << "Incorrect call to FBSparseArray1D::setupIndex() (pass=1). This function should be called after allocate() and before setupIndex2().";
			exit(0);
			d->m_has_error=true;
			return false;
		}
		long block_num=ind/SPARSE_BLOCK_SIZE;
		d->m_block_entry_counts[block_num]++;
	}
	else if (pass==2) {
		if (d->m_stage==SPARSE_STAGE_ALLOCATED)
			d->finalize_step_1();
		if (d->m_stage!=SPARSE_STAGE_STEP1) {
			qWarning() << "Incorrect call to FBSparseArray1D::setupIndex() (pass=2). This function should be called after setupIndex() (pass=1) and before setValue() and value().";
			exit(0);
			d->m_has_error=true;
			return false;
		}
		if (ind<=d->m_last_step2_index) {
			qWarning() << "Incorrect call to FBSparseArray1D::setupIndex() (pass=2). Index out of order:" << d->m_last_step2_index << ind;
			exit(0);
			d->m_has_error=true;
			return false;
		}
		long block_index=ind/SPARSE_BLOCK_SIZE;
		long last_block_index=d->m_last_step2_index/SPARSE_BLOCK_SIZE;
		int index_in_block=0;
		if (d->m_last_step2_index<0) last_block_index=-1;
		if (block_index==last_block_index) index_in_block=d->m_last_step2_index_in_block+1;
		
		if (index_in_block>=d->m_block_entry_counts[block_index]) {
			qWarning() << "Incorrect call to FBSparseArray1D::setupIndex2() (pass=2). Index in block is too large." << index_in_block << d->m_block_entry_counts[block_index] << block_index;
			exit(0);
			d->m_has_error=true;
			return false;
		}
		d->m_offsets[d->m_block_indices[block_index]+index_in_block]=ind%SPARSE_BLOCK_SIZE;
		
		d->m_last_step2_index=ind;
		d->m_last_step2_index_in_block=index_in_block;
	}
	else {
		qWarning() << "Error in setupIndex(), invalid pass number:" << pass;
		exit(0);
		d->m_has_error=true;
		return false;
	}
	return true;
}

long FBSparseArray1D::length() const {
	if (d->m_has_error) return 0;
	return d->m_N;
}
void FBSparseArray1D::setValue(float val,long ind) {
	if (d->m_has_error) return;
	if (d->m_stage==SPARSE_STAGE_STEP1) d->finalize_step_2();
	if (d->m_stage!=SPARSE_STAGE_STEP2) {
		qWarning() << "Incorrect call to FBSparseArray1D::setValue(). This function should be called after setupIndex() (pass=2).";
		exit(0);
		d->m_has_error=true;
		return;
	}
	long blockind=ind/SPARSE_BLOCK_SIZE;
	int offset_index=d->find_offset_index(blockind,ind%SPARSE_BLOCK_SIZE);
	if (offset_index<0) {
		qWarning() << "Incorrect call to FBSparseArray1D::setValue(). Index does not exist in sparse array:" << ind;
		exit(0);
		d->m_has_error=true;
		return;
	}
	if (d->m_data_type==DATA_TYPE_FLOAT)
		d->m_data_float[d->m_block_indices[blockind]+offset_index]=val;
	else if (d->m_data_type==DATA_TYPE_BYTE)
		d->m_data_byte[d->m_block_indices[blockind]+offset_index]=(unsigned char)val;
}
float FBSparseArray1D::value(long ind) const {
	if (d->m_has_error) return 0;
	if ((ind<0)||(ind>=d->m_N)) return 0;
	if (d->m_stage==SPARSE_STAGE_STEP1) d->finalize_step_2();
	if (d->m_stage!=SPARSE_STAGE_STEP2) {
		qWarning() << "Incorrect call to FBSparseArray1D::value(). This function should be called after setupIndex2() (pass=2).";
		exit(0);
		d->m_has_error=true;
		return 0;
	}
	long blockind=ind/SPARSE_BLOCK_SIZE;
	int offset_index=d->find_offset_index(blockind,ind%SPARSE_BLOCK_SIZE);
	if (offset_index<0) return 0;
	if (d->m_data_type==DATA_TYPE_FLOAT)
		return d->m_data_float[d->m_block_indices[blockind]+offset_index];
	else if (d->m_data_type==DATA_TYPE_BYTE)
		return d->m_data_byte[d->m_block_indices[blockind]+offset_index];
	else 
		return 0;
}
void FBSparseArray1DPrivate::finalize_step_1() {
	m_entry_count=0;
	for (long blockind=0; blockind<m_block_count; blockind++) {
		m_block_indices[blockind]=m_entry_count;
		m_entry_count+=m_block_entry_counts[blockind];
	}
	if (m_entry_count) {
		if (m_data_type==DATA_TYPE_FLOAT) {
			m_data_float=(float *)malloc(sizeof(float)*m_entry_count);
			for (long i=0; i<m_entry_count; i++) m_data_float[i]=0;
		}
		else if (m_data_type==DATA_TYPE_BYTE) {
			m_data_byte=(unsigned char *)malloc(sizeof(unsigned char)*m_entry_count);
			for (long i=0; i<m_entry_count; i++) m_data_byte[i]=0;
		}
		m_offsets=(unsigned char *)malloc(sizeof(unsigned char)*m_entry_count);
	}
	
	m_last_step2_index=-1;
	m_last_step2_index_in_block=-1;
	m_stage=SPARSE_STAGE_STEP1;
}
void FBSparseArray1DPrivate::finalize_step_2() {
	m_stage=SPARSE_STAGE_STEP2;
}
int FBSparseArray1DPrivate::find_offset_index(long blocknum,unsigned char offset) {
	unsigned char *offsets=&m_offsets[m_block_indices[blocknum]];
	int N=m_block_entry_counts[blocknum];
	if (N==0) return -1;
	int minval=0;
	int maxval=N-1;
	int testval=N/2;
	bool done=false;
	while (!done) {
		int old_testval=testval;
		if (offsets[testval]==offset) return testval;
		else if (offsets[testval]<offset) {
			minval=testval+1;
			testval=(testval+1+maxval)/2;
		}
		else if (offsets[testval]>offset) {
			maxval=testval-1;
			testval=(testval-1+minval)/2;
		}
		if (testval==old_testval) return -1;
	}
	return -1; //shouldn't ever get here anyway
}
long FBSparseArray1D::entryCount() const {
	if (d->m_has_error) return 0;
	return d->m_entry_count;
}
void FBSparseArray1D::resetIteration() const {
	d->m_current_block_index=-1;
	d->m_current_index_in_block=-1;
}
bool FBSparseArray1D::advanceIteration() const {
	if (d->m_current_block_index<0) {
		d->m_current_block_index=0;
		d->m_current_index_in_block=0;
		if (d->m_current_block_index>=d->m_block_count) return false;
		while (d->m_current_index_in_block>=d->m_block_entry_counts[d->m_current_block_index]) {
			d->m_current_block_index++;
			if (d->m_current_block_index>=d->m_block_count) return false;
		}
	}
	else {
		if (d->m_current_block_index>=d->m_block_count) return false;
		d->m_current_index_in_block++;
		while (d->m_current_index_in_block>=d->m_block_entry_counts[d->m_current_block_index]) {
			d->m_current_block_index++;
			d->m_current_index_in_block=0;
			if (d->m_current_block_index>=d->m_block_count) return false;	
		}
	}
	return true;
}
long FBSparseArray1D::currentIndex() const {
	if (d->m_current_block_index<0) return -1;
	return d->m_current_block_index*SPARSE_BLOCK_SIZE+d->m_offsets[d->m_block_indices[d->m_current_block_index]+d->m_current_index_in_block];
}
float FBSparseArray1D::currentValue() const {
	if (d->m_data_type==DATA_TYPE_FLOAT)
		return d->m_data_float[d->m_block_indices[d->m_current_block_index]+d->m_current_index_in_block];
	else if (d->m_data_type==DATA_TYPE_BYTE)
		return d->m_data_byte[d->m_block_indices[d->m_current_block_index]+d->m_current_index_in_block];
	else
		return 0;
}
void FBSparseArray1D::operator=(const FBSparseArray1D &X) {
	d->copy_from(X);
}





