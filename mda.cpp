#include "mda.h"
#include <stdio.h>
#include <QDebug>
#include <QTime>
#include <QFileInfo>
//#include "jniftiimage2.h"

real MdaElement::re() const {
	switch (m_data_type) {
		case MDA_TYPE_BYTE: return *((unsigned char *)m_data); break;
		case MDA_TYPE_SHORT: return *((qint16 *)m_data); break;
		case MDA_TYPE_UINT16: return *((quint16 *)m_data); break;
		case MDA_TYPE_INT32: return *((qint32 *)m_data); break;
		case MDA_TYPE_REAL: return *((real *)m_data); break;
		case MDA_TYPE_COMPLEX: return ((complex_struct *)m_data)->re; break;
	};
	return 0;
}
real MdaElement::im() const {
	if (m_data_type==MDA_TYPE_COMPLEX) {
		return ((complex_struct *)m_data)->im;
	}
	else return 0;
}


void MdaElement::setRe(real val) {
	switch (m_data_type) {
		case MDA_TYPE_BYTE: *((unsigned char *)m_data)=(unsigned char)val; break;
		case MDA_TYPE_SHORT: *((qint16 *)m_data)=(qint16)val; break;
		case MDA_TYPE_UINT16: *((quint16 *)m_data)=(quint16)val; break;
		case MDA_TYPE_INT32: *((qint32 *)m_data)=(qint32)val; break;
		case MDA_TYPE_REAL: *((real *)m_data)=(real)val; break;
		case MDA_TYPE_COMPLEX: ((complex_struct *)m_data)->re=val; break;
	};
}

void MdaElement::setIm(real val) {
	if (m_data_type==MDA_TYPE_COMPLEX) {
		((complex_struct *)m_data)->im=val;
	};
}

void Mda::initialize(int datatype) {
	data_byte=0;
	data_real=0;
	data_short=0;
	data_uint16=0;
	data_int32=0;
	data_complex=0;
	data_mask=0;
	allocate(datatype,1,1);
}

Mda::Mda() 
{
	all_indices=false;
	initialize(MDA_TYPE_REAL);
}
Mda::Mda(const Mda &X) {
	all_indices=false;
	initialize(MDA_TYPE_REAL);
	copy_from(X);
}

void Mda::operator=(real val) {
	allocate(MDA_TYPE_REAL,1,1);
	data_real[0]=val;
}
void Mda::operator=(const AbstractComplex &val) {
	allocate(MDA_TYPE_COMPLEX,1,1);
	data_complex[0].re=val.re();
	data_complex[0].im=val.im();
}

void Mda::operator=(const Mda &X) {
	copy_from(X);
}

Mda::~Mda()
{
	clear(false);
}

void Mda::allocate(int datatype,int num_dimensions,qint32 *dims) {
	clear(false);
	m_data_type=datatype;
	m_dim_count=num_dimensions;
	m_dimprod=0;
	if (num_dimensions<=0) return;
	m_dimprod=1;
	for (int j=0; j<num_dimensions; j++) {
		m_dimprod*=dims[j];
		m_dims[j]=dims[j];
	}
	if (m_dimprod==0) m_dimprod=1;
	switch (datatype) {
		case MDA_TYPE_BYTE: 
			data_byte=(unsigned char *)malloc(sizeof(unsigned char)*m_dimprod);
			break;
		case MDA_TYPE_SHORT: 
			data_short=(qint16 *)malloc(sizeof(qint16)*m_dimprod);
			break;
		case MDA_TYPE_UINT16: 
			data_uint16=(quint16 *)malloc(sizeof(quint16)*m_dimprod);
			break;
		case MDA_TYPE_INT32: 
			data_int32=(qint32 *)malloc(sizeof(qint32)*m_dimprod);
			break;
		case MDA_TYPE_REAL: 
			data_real=(real *)malloc(sizeof(real)*m_dimprod);
			break;
		case MDA_TYPE_COMPLEX: 
			data_complex=(complex_struct *)malloc(sizeof(complex_struct)*m_dimprod);
			break;
	}
}

void Mda::allocate(int datatype,qint32 S1,qint32 S2,qint32 S3,qint32 S4,qint32 S5,qint32 S6,qint32 S7,qint32 S8) {
	qint32 dims[MAX_MDA_DIMS];
	int num_dimensions=0;
	if (S1>=0) {
		num_dimensions++;
		dims[0]=S1;
	if (S2>=0) {
		num_dimensions++;
		dims[1]=S2;
	if (S3>=0) {
		num_dimensions++;
		dims[2]=S3;
	if (S4>=0) {
		num_dimensions++;
		dims[3]=S4;
	if (S5>=0) {
		num_dimensions++;
		dims[4]=S5;
	if (S6>=0) {
		num_dimensions++;
		dims[5]=S6;
	if (S7>=0) {
		num_dimensions++;
		dims[6]=S7;
	if (S8>=0) {
		num_dimensions++;
		dims[7]=S8;
	}}}}}}}}
	allocate(datatype,num_dimensions,dims);
}



void Mda::setAll(const AbstractComplex &val) {
	qint32 N=size();
	for (int j=0; j<N; j++)
		(*this)[j]=val;
}
void Mda::setAll(real val) {
	qint32 N=size();
	for (int j=0; j<N; j++)
		(*this)[j]=val;
}

void Mda::clear(bool allocate_one_by_one) {
	if (data_byte) free(data_byte); data_byte=0;
	if (data_short) free(data_short); data_short=0;
	if (data_uint16) free(data_uint16); data_uint16=0;
	if (data_int32) free(data_int32); data_int32=0;
	if (data_real) free(data_real); data_real=0;
	if (data_complex) free(data_complex); data_complex=0;
	if (data_mask) free(data_mask); data_mask=0;
	m_data_type=0;
	m_dim_count=0;
	m_dimprod=0;
	if (allocate_one_by_one) allocate(MDA_TYPE_REAL,1,1);
}



MdaElement Mda::operator[](const qint64 index) {
	if (index<0) {
		qWarning() << "WARNING!!! Index is negative" << __FUNCTION__ << index;
	}
	MdaElement ret;
	ret.m_data_type=m_data_type;
	{
		switch (m_data_type) {
			case MDA_TYPE_BYTE: ret.m_data=&data_byte[index]; break;
			case MDA_TYPE_SHORT: ret.m_data=&data_short[index]; break;
			case MDA_TYPE_UINT16: ret.m_data=&data_uint16[index]; break;
			case MDA_TYPE_INT32: ret.m_data=&data_int32[index]; break;
			case MDA_TYPE_REAL: ret.m_data=&data_real[index]; break;
			case MDA_TYPE_COMPLEX: ret.m_data=&data_complex[index]; break;
			default: ret.m_data=0;
		};
	}
	return ret;
}

qint64 Mda::get_index(qint32 *index) const {
	qint64 ret=0;
	qint64 factor=1;
	for (int j=0; j<dimCount(); j++) {
		ret+=index[j]*factor;
		factor*=size(j);
	}
	return ret;
}

qint64 Mda::get_index (qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) const {
	qint64 index=ind1;
	qint64 factor=1;
	if ((ind2>=0)&&(m_dim_count>=2)) {
		factor*=m_dims[0];
		index+=ind2*factor;
	if ((ind3>=0)&&(m_dim_count>=3)) {
		factor*=m_dims[1];
		index+=ind3*factor;
	if ((ind4>=0)&&(m_dim_count>=4)) {
		factor*=m_dims[2];
		index+=ind4*factor;
	if ((ind5>=0)&&(m_dim_count>=5)) {
		factor*=m_dims[3];
		index+=ind5*factor;
	if ((ind6>=0)&&(m_dim_count>=6)) {
		factor*=m_dims[4];
		index+=ind6*factor;
	if ((ind7>=0)&&(m_dim_count>=7)) {
		factor*=m_dims[5];
		index+=ind7*factor;
	if ((ind8>=0)&&(m_dim_count>=8)) {
		factor*=m_dims[6];
		index+=ind8*factor;
	}}}}}}}
	return index;
}

bool Mda::inRange(qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) {
	if (ind2<0) return ((0<=ind1)&&(ind1<m_dimprod));
	if (m_dim_count>=1) {
		if ((ind1<0)||(ind1>=m_dims[0])) return false;
	}
	else return true;
	if (m_dim_count>=2) {
		if ((ind2<0)||(ind2>=m_dims[1])) return false;
	}
	else return true;
	if (m_dim_count>=3) {
		if ((ind3<0)||(ind3>=m_dims[2])) return false;
	}
	else return true;
	if (m_dim_count>=4) {
		if ((ind4<0)||(ind4>=m_dims[3])) return false;
	}
	else return true;
	if (m_dim_count>=5) {
		if ((ind5<0)||(ind5>=m_dims[4])) return false;
	}
	else return true;
	if (m_dim_count>=6) {
		if ((ind6<0)||(ind6>=m_dims[5])) return false;
	}
	else return true;
	if (m_dim_count>=7) {
		if ((ind7<0)||(ind7>=m_dims[6])) return false;
	}
	else return true;
	if (m_dim_count>=8) {
		if ((ind8<0)||(ind8>=m_dims[7])) return false;
	}
	else return true;
	return false;
}

MdaElement Mda::operator()(qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) {
	qint64 index=get_index(ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8);
	if (index<0) {
		qWarning() << "Index is negative" << __FUNCTION__ << ind1 << ind2 << ind3 << ind4;
	}
	return (*this)[index];
}
MdaElement Mda::operator()(qint32 *index_in) {
	qint64 index=get_index(index_in);
	if (index<0) {
		qWarning() << "Index is negative" << __FUNCTION__ << index_in[0] << index_in[1] << index_in[2] << index_in[3];
	}
	return (*this)[index];
}
MDAComplex Mda::get(qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) const {
	qint64 index;
	if (ind2<0) index=ind1;
	else index=get_index(ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8);
	if (index<0) return MDAComplex(0,0);
	if (index>=m_dimprod) return MDAComplex(0,0);
	{
		switch (m_data_type) {
			case MDA_TYPE_BYTE: return MDAComplex(data_byte[index],0); break;
			case MDA_TYPE_SHORT: return MDAComplex(data_short[index],0); break;
			case MDA_TYPE_UINT16: return MDAComplex(data_uint16[index],0); break;
			case MDA_TYPE_INT32: return MDAComplex(data_int32[index],0); break;
			case MDA_TYPE_REAL: return MDAComplex(data_real[index],0); break;
			case MDA_TYPE_COMPLEX: return MDAComplex(data_complex[index].re,data_complex[index].im); break;
			default: return MDAComplex(0,0);
		};
	}
	return MDAComplex(0,0);
}
MDAComplex Mda::get(qint32 *index) const {
	qint64 index2=get_index(index);
	return get(index2);
}

unsigned char &Mda::mask(qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) {
	qint64 index=get_index(ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8);
	{
		if (!data_mask) allocate_mask();
		return data_mask[index];
	}
}
unsigned char &Mda::mask(qint32 *index_in) {
	qint64 index=get_index(index_in);
	{
		if (!data_mask) allocate_mask();
		return data_mask[index];
	}
}

unsigned char Mda::getmask(qint64 ind1,qint32 ind2,qint32 ind3,qint32 ind4,qint32 ind5,qint32 ind6,qint32 ind7,qint32 ind8) const {
	qint64 index=get_index(ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8);
	{
		if (!data_mask) return 0;
		return data_mask[index];
	}
}
unsigned char Mda::getmask(qint32 *index_in) const {
	qint64 index=get_index(index_in);
	return getmask(index);
}


int Mda::dimCount() const {
	return m_dim_count;
}
qint64 Mda::size() const {
	return m_dimprod;
}
qint32 Mda::size(int dim) const {
	if (dim>=m_dim_count) return 1;
	return m_dims[dim];
}

int Mda::dataType() const {
	return m_data_type;
}

void Mda::allocate_mask() {
	if ((!data_mask)&&(m_dim_count>0)&&(m_dimprod>0)) {
		data_mask=(unsigned char *)malloc(m_dimprod*sizeof(unsigned char));
		for (int j=0; j<m_dimprod; j++)
			data_mask[j]=0;
	}
}

bool Mda::read(const QString &filename) {

	QString suf=QFileInfo(filename).suffix();
	/*
	if (suf=="nii") {
		JNiftiImage2 NN;
		NN.read(filename);
		allocate(MDA_TYPE_REAL,NN.N1(),NN.N2(),NN.N3(),NN.N4());
		for (int i4=0; i4<NN.N4(); i4++) {
			Array3D tmp=NN.dataXYZ(i4);
			for (int i3=0; i3<NN.N3(); i3++)
			for (int i2=0; i2<NN.N2(); i2++)
			for (int i1=0; i1<NN.N1(); i1++) {
				(*this)(i1,i2,i3,i4)=tmp.value(i1,i2,i3);
			}
		}
		return true;
	}
	*/
	
	
	QString fname=filename;
	#ifdef WIN32
	fname.replace("/","\\");
	#endif
	clear();
	FILE *inf=fopen(fname.toAscii().data(),"rb");
	if (!inf) return false;
	qint32 hold_num_dims;
	qint32 hold_dims[MAX_MDA_DIMS];
	if (fread(&hold_num_dims,sizeof(qint32),1,inf)<=0)
		return false;
	qint32 data_type;
	if (hold_num_dims<0) {
		data_type=hold_num_dims;
		qint32 num_bytes;
		fread(&num_bytes,sizeof(qint32),1,inf);
		fread(&hold_num_dims,sizeof(qint32),1,inf);
	}
	else
		data_type=MDA_TYPE_COMPLEX;
	if (hold_num_dims>MAX_MDA_DIMS) return false;
	if (hold_num_dims<=0) return false;
	qint32 j;
	for (j=0; j<hold_num_dims; j++) {
		qint32 holdval;
		fread(&holdval,sizeof(qint32),1,inf);
		hold_dims[j]=holdval;
	}
	{
		allocate(data_type,hold_num_dims,hold_dims);
		qint32 N=size();
		if (data_type==MDA_TYPE_COMPLEX) {
			long bufsize=10000; //must be even!
			real *real_buf=(real *)malloc(4*bufsize);
			for (long ii=0; ii<N*2; ii+=bufsize) {
				long len=qMin(bufsize,N*2-ii);
				fread(real_buf,4,len,inf);
				for (int kk=0; kk<len/2; kk++) {
					data_complex[ii/2+kk].re=real_buf[kk*2];
					data_complex[ii/2+kk].im=real_buf[kk*2+1];
				}
			}
			free(real_buf);
		}
		else if (data_type==MDA_TYPE_REAL) {
			fread(data_real,4,N,inf);
		}
		else if (data_type==MDA_TYPE_SHORT) {
			fread(data_short,2,N,inf);
		}
		else if (data_type==MDA_TYPE_UINT16) {
			fread(data_uint16,2,N,inf);
		}
		else if (data_type==MDA_TYPE_INT32) {
			fread(data_int32,4,N,inf);
		}
		else if (data_type==MDA_TYPE_BYTE) {
			fread(data_byte,1,N,inf);
		}
		/*float re,im;
		for (j=0; j<N; j++) {
			if ((data_type==MDA_TYPE_COMPLEX)&&(data_complex)) {
				fread(&re,sizeof(float),1,inf);
				fread(&im,sizeof(float),1,inf);
				data_complex[j].re=re;
				data_complex[j].im=im;
			}
		}*/
		qint32 com;
		if (fread(&com,sizeof(qint32),1,inf)) {
			if (com==MDA_COMMAND_MASK) {
				allocate_mask();
				fread(data_mask,1,N,inf);
				//for (j=0; j<N; j++)
				//	fread(&data_mask[j],sizeof(char),1,inf);
			}
		}
	}
	fclose(inf);
	return true;
}

bool Mda::readSize(const QString &filename,Mda &ret) {
	QString fname=filename;
	#ifdef WIN32
	fname.replace("/","\\");
	#endif
	FILE *inf=fopen(fname.toAscii().data(),"rb");
	if (!inf) return false;
	qint32 hold_num_dims;
	qint32 hold_dims[MAX_MDA_DIMS];
	if (fread(&hold_num_dims,sizeof(qint32),1,inf)<=0)
		return false;
	qint32 data_type;
	if (hold_num_dims<0) {
		data_type=hold_num_dims;
		qint32 num_bytes;
		fread(&num_bytes,sizeof(qint32),1,inf);
		fread(&hold_num_dims,sizeof(qint32),1,inf);
	}
	else
		data_type=MDA_TYPE_COMPLEX;
	if (hold_num_dims>MAX_MDA_DIMS) return false;
	if (hold_num_dims<=0) return false;
	qint32 j;
	for (j=0; j<hold_num_dims; j++) {
		qint32 holdval;
		fread(&holdval,sizeof(qint32),1,inf);
		hold_dims[j]=holdval;
	}

	ret.allocate(MDA_TYPE_REAL,hold_num_dims,1);
	for (int j=0; j<hold_num_dims; j++)
		ret[j]=hold_dims[j];	
	
	fclose(inf);
	return true;
}

bool Mda::readSlices(const QString &filename,long minslice,long maxslice) {
	QString fname=filename;
	#ifdef WIN32
	fname.replace("/","\\");
	#endif
	if (maxslice<minslice) return false;
	if (minslice<0) return false;
	clear();
	FILE *inf=fopen(fname.toAscii().data(),"rb");
	if (!inf) return false;
	qint32 hold_num_dims;
	qint32 hold_dims[MAX_MDA_DIMS];
	qint32 hold_dims2[MAX_MDA_DIMS];
	if (fread(&hold_num_dims,sizeof(qint32),1,inf)<=0) {
		fclose(inf);
		return false;
	}
	qint32 data_type;
	qint32 num_bytes=8;
	if (hold_num_dims<0) {
		data_type=hold_num_dims;
		fread(&num_bytes,sizeof(qint32),1,inf);
		fread(&hold_num_dims,sizeof(qint32),1,inf);
	}
	else
		data_type=MDA_TYPE_COMPLEX;
	if (hold_num_dims>MAX_MDA_DIMS) {
		fclose(inf);
		return false;
	}
	if (hold_num_dims<=0) {
		fclose(inf);
		return false;
	}
	qint32 j;
	for (j=0; j<hold_num_dims; j++) {
		qint32 holdval;
		fread(&holdval,sizeof(qint32),1,inf);
		hold_dims[j]=holdval;
	}
	if (hold_num_dims<3) {
		fclose(inf);
		return false;
	}
	if (maxslice>=hold_dims[2]) {
		fclose(inf);
		return false;
	}
	hold_dims2[0]=hold_dims[0];
	hold_dims2[1]=hold_dims[1];
	hold_dims2[2]=maxslice-minslice+1;
	long to_skip=minslice*hold_dims[0]*hold_dims[1];
	
	allocate(data_type,hold_num_dims,hold_dims2);
	qint32 N=size();
	fseek(inf,to_skip*num_bytes,SEEK_CUR);
	/*float re,im;
	for (j=0; j<N; j++) {
		if ((data_type==MDA_TYPE_COMPLEX)&&(data_complex)) {
			fread(&re,sizeof(float),1,inf);
			fread(&im,sizeof(float),1,inf);
			data_complex[j].re=re;
			data_complex[j].im=im;
		}
		else if ((data_type==MDA_TYPE_REAL)&&(data_real))
			fread(&data_real[j],4,1,inf);
		else if ((data_type==MDA_TYPE_SHORT)&&(data_short))
			fread(&data_short[j],2,1,inf);
		else if ((data_type==MDA_TYPE_INT32)&&(data_int32))
			fread(&data_int32[j],4,1,inf);
		else if ((data_type==MDA_TYPE_BYTE)&&(data_byte))
			fread(&data_byte[j],1,1,inf);
	}*/
	if (data_type==MDA_TYPE_COMPLEX) {
		long bufsize=10000; //must be even!
		real *real_buf=(real *)malloc(4*bufsize);
		for (long ii=0; ii<N*2; ii+=bufsize) {
			long len=qMin(bufsize,N*2-ii);
			fread(real_buf,4,len,inf);
			for (int kk=0; kk<len/2; kk++) {
				data_complex[ii/2+kk].re=real_buf[kk*2];
				data_complex[ii/2+kk].im=real_buf[kk*2+1];
			}
		}
		free(real_buf);
	}
	else if (data_type==MDA_TYPE_REAL) {
		fread(data_real,4,N,inf);
	}
	else if (data_type==MDA_TYPE_SHORT) {
		fread(data_short,2,N,inf);
	}
	else if (data_type==MDA_TYPE_UINT16) {
		fread(data_uint16,2,N,inf);
	}
	else if (data_type==MDA_TYPE_INT32) {
		fread(data_int32,4,N,inf);
	}
	else if (data_type==MDA_TYPE_BYTE) {
		fread(data_byte,1,N,inf);
	}
	/*qint32 com;
	if (fread(&com,sizeof(qint32),1,inf)) {
		if (com==MDA_COMMAND_MASK) {
			allocate_mask();
			fseek(inf,to_skip*num_bytes,SEEK_CUR);
			for (j=0; j<N; j++)
				fread(&data_mask[j],sizeof(char),1,inf);
		}
	}*/
	fclose(inf);
	return true;
}

bool Mda::write(const QString &filename) const {
	QString fname=filename;
	#ifdef WIN32
	fname.replace("/","\\");
	#endif
	FILE *outf=fopen(fname.toAscii().data(),"wb");
	if (!outf) return false;
	if (m_data_type!=MDA_TYPE_COMPLEX) {
		qint32 data_type=m_data_type;
		fwrite(&data_type,sizeof(qint32),1,outf);
		qint32 num_bytes=1;
		switch (data_type) {
			case MDA_TYPE_BYTE: num_bytes=1; break;
			case MDA_TYPE_SHORT: num_bytes=2; break;
			case MDA_TYPE_UINT16: num_bytes=2; break;
			case MDA_TYPE_INT32: num_bytes=4; break;
			case MDA_TYPE_REAL: num_bytes=4; break;
			case MDA_TYPE_COMPLEX: num_bytes=8; break;
		}
		fwrite(&num_bytes,sizeof(qint32),1,outf);
	}
	qint32 num_dims=m_dim_count;
	fwrite(&num_dims,sizeof(qint32),1,outf);
	qint32 j;
	for (j=0; j<num_dims; j++) {
		qint32 holdval=m_dims[j];
		fwrite(&holdval,sizeof(qint32),1,outf);
	}
	{
		qint32 N=size();
		if (m_data_type==MDA_TYPE_COMPLEX) {
			long bufsize=10000; //must be even!
			real *real_buf=(real *)malloc(4*bufsize);
			for (long ii=0; ii<N*2; ii+=bufsize) {
				long len=qMin(bufsize,N*2-ii);
				for (int kk=0; kk<len/2; kk++) {
					real_buf[kk*2]=data_complex[ii/2+kk].re;
					real_buf[kk*2+1]=data_complex[ii/2+kk].im;
				}
				fwrite(real_buf,4,len,outf);
			}
			free(real_buf);
		}
		else if (m_data_type==MDA_TYPE_REAL) {
			fwrite(data_real,4,N,outf);
		}
		else if (m_data_type==MDA_TYPE_SHORT) {
			fwrite(data_short,2,N,outf);
		}
		else if (m_data_type==MDA_TYPE_UINT16) {
			fwrite(data_uint16,2,N,outf);
		}
		else if (m_data_type==MDA_TYPE_INT32) {
			fwrite(data_int32,4,N,outf);
		}
		else if (m_data_type==MDA_TYPE_BYTE) {
			fwrite(data_byte,1,N,outf);
		}
		/*float re,im;
		for (j=0; j<N; j++) {
			if ((m_data_type==MDA_TYPE_COMPLEX)&&(data_complex)) {
				re=data_complex[j].re;
				im=data_complex[j].im;
				fwrite(&re,sizeof(float),1,outf);
				fwrite(&im,sizeof(float),1,outf);
			}
			else if ((m_data_type==MDA_TYPE_REAL)&&(data_real))
				fwrite(&(data_real)[j],4,1,outf);
			else if ((m_data_type==MDA_TYPE_SHORT)&&(data_short))
				fwrite(&(data_short)[j],2,1,outf);
			else if ((m_data_type==MDA_TYPE_INT32)&&(data_int32))
				fwrite(&(data_int32)[j],2,1,outf);
			else if ((m_data_type==MDA_TYPE_BYTE)&&(data_byte))
				fwrite(&(data_byte)[j],1,1,outf);
		}*/
		if (data_mask) {
			qint32 com=MDA_COMMAND_MASK;
			fwrite(&com,sizeof(qint32),1,outf);
			//for (j=0; j<N; j++)
			//	fwrite(&data_mask[j],sizeof(char),1,outf);
			fwrite(data_mask,1,N,outf);
		}
		else {
			qint32 com=0;
			fwrite(&com,sizeof(qint32),1,outf);
		}
	}
	fclose(outf);
	return true;
}

void Mda::convertToComplex() {
	if (dataType()==MDA_TYPE_COMPLEX) return;
	
	{
		data_complex=(complex_struct *)malloc(sizeof(complex_struct)*m_dimprod);	
		if (dataType()==MDA_TYPE_REAL) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_complex[j].re=data_real[j];
				data_complex[j].im=0;
			}
			free(data_real);
			data_real=0;
		}
		else if (dataType()==MDA_TYPE_BYTE) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_complex[j].re=(real)data_byte[j];
				data_complex[j].im=0;
			}
			free(data_byte);
			data_byte=0;
		}
		else if (dataType()==MDA_TYPE_SHORT) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_complex[j].re=(real)data_short[j];
				data_complex[j].im=0;
			}
			free(data_short);
			data_short=0;
		}
		else if (dataType()==MDA_TYPE_UINT16) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_complex[j].re=(real)data_uint16[j];
				data_complex[j].im=0;
			}
			free(data_uint16);
			data_uint16=0;
		}
		else if (dataType()==MDA_TYPE_INT32) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_complex[j].re=(real)data_int32[j];
				data_complex[j].im=0;
			}
			free(data_int32);
			data_int32=0;
		}
		else {
			free(data_complex);
			data_complex=0;
			qWarning() << "Unexpected problem:" << __FILE__ << __LINE__ << __FUNCTION__;
			return;
		}
	}
	m_data_type=MDA_TYPE_COMPLEX;
}

void Mda::convertToReal() {
	if (dataType()==MDA_TYPE_REAL) return;
	
	{
		data_real=(real *)malloc(sizeof(real)*m_dimprod);	
		if (dataType()==MDA_TYPE_COMPLEX) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_real[j]=data_complex[j].re;
			}
			free(data_complex);
			data_complex=0;
		}
		else if (dataType()==MDA_TYPE_BYTE) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_real[j]=(real)data_byte[j];
			}
			free(data_byte);
			data_byte=0;
		}
		else if (dataType()==MDA_TYPE_SHORT) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_real[j]=(real)data_short[j];
			}
			free(data_short);
			data_short=0;
		}
		else if (dataType()==MDA_TYPE_UINT16) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_real[j]=(real)data_uint16[j];
			}
			free(data_uint16);
			data_uint16=0;
		}
		else if (dataType()==MDA_TYPE_INT32) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_real[j]=(real)data_int32[j];
			}
			free(data_int32);
			data_int32=0;
		}
		else {
			free(data_real);
			data_real=0;
			qWarning() << "Unexpected problem:" << __FILE__ << __LINE__ << __FUNCTION__;
			return;
		}
	}
	m_data_type=MDA_TYPE_REAL;
}
void Mda::convertToByte() {
	if (dataType()==MDA_TYPE_BYTE) return;
	
	{
		data_byte=(unsigned char *)malloc(sizeof(unsigned char)*m_dimprod);	
		if (dataType()==MDA_TYPE_COMPLEX) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_byte[j]=(unsigned char)data_complex[j].re;
			}
			free(data_complex);
			data_complex=0;
		}
		else if (dataType()==MDA_TYPE_REAL) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_byte[j]=(unsigned char)data_real[j];
			}
			free(data_real);
			data_real=0;
		}
		else if (dataType()==MDA_TYPE_SHORT) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_byte[j]=(unsigned char)data_short[j];
			}
			free(data_short);
			data_short=0;
		}
		else if (dataType()==MDA_TYPE_UINT16) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_byte[j]=(unsigned char)data_uint16[j];
			}
			free(data_uint16);
			data_uint16=0;
		}
		else if (dataType()==MDA_TYPE_INT32) {
			for (qint32 j=0; j<m_dimprod; j++) {
				data_byte[j]=(unsigned char)data_int32[j];
			}
			free(data_int32);
			data_int32=0;
		}
		else {
			free(data_byte);
			data_byte=0;
			qWarning() << "Unexpected problem:" << __FILE__ << __LINE__ << __FUNCTION__;
			return;
		}
	}
	m_data_type=MDA_TYPE_BYTE;
}

void Mda::copy_from(const Mda &X) {
	clear(true);
	qint32 hold_dims[MAX_MDA_DIMS];
	for (int j=0; j<X.dimCount(); j++)
		hold_dims[j]=X.size(j);
	{
		allocate(X.dataType(),X.dimCount(),hold_dims);
		switch (X.dataType()) {
			case MDA_TYPE_BYTE: memcpy(data_byte,X.data_byte,sizeof(unsigned char)*X.size());
											break;
			case MDA_TYPE_SHORT: memcpy(data_short,X.data_short,sizeof(qint16)*X.size());
											break;
			case MDA_TYPE_UINT16: memcpy(data_uint16,X.data_uint16,sizeof(quint16)*X.size());
											break;
			case MDA_TYPE_INT32: memcpy(data_int32,X.data_int32,sizeof(qint32)*X.size());
											break;
			case MDA_TYPE_REAL: memcpy(data_real,X.data_real,sizeof(real)*X.size());
											break;
			case MDA_TYPE_COMPLEX: memcpy(data_complex,X.data_complex,sizeof(complex_struct)*X.size());
											break;
		}
		if (X.data_mask) {
			allocate_mask();
			memcpy(data_mask,X.data_mask,sizeof(unsigned char)*X.size());
		}
	}
}

void Mda::shift(qint32 *shifts) {
	Mda hold;
	hold=(*this);
	qint32 inds[MAX_MDA_DIMS],inds2[MAX_MDA_DIMS];
	qint32 j;
	for (j=0; j<m_dim_count; j++) {
		inds[j]=0;
		if (size(j)==0) return;
		while (shifts[j]<0) shifts[j]+=size(j);
		inds2[j]=(inds[j]+shifts[j])%size(j);
	}
	

	qint32 k;
	bool done=false;
	while (!done) {
		hold(inds2)=(*this)(inds);

		//increment
		k=-1;
		do {
			k++;
			if (k<m_dim_count) {
				inds[k]++;
				if (inds[k]>=m_dims[k]) {
					inds[k]=0;
				}
				inds2[k]=(inds[k]+shifts[k])%size(k);
			}
		} while ((k<hold.dimCount())&&(inds[k]==0));
		if (k>=hold.dimCount()) done=true;
	}
	

	copy_from(hold);
	
}

void Mda::squeeze() {
	qint32 new_dims[MAX_MDA_DIMS];
	qint32 old_dim_count=m_dim_count;
	qint32 j,ct=0;
	for (j=0; j<m_dim_count; j++) {
		if (m_dims[j]>1) {
			new_dims[ct]=m_dims[j];
			ct++;
		}
	}
	m_dim_count=ct;
	for (j=0; j<m_dim_count; j++) {
		m_dims[j]=new_dims[j];
	}
	for (j=m_dim_count; j<old_dim_count; j++)
		m_dims[j]=1;
	if (m_dim_count<1) {
		m_dim_count=1;
		m_dims[0]=1;
	}
}



void Mda::reshape(int dim_count,qint32 *dims) {
	long prod=1;
	for (int j=0; j<dim_count; j++) prod*=dims[j];
	if (prod!=m_dimprod) return;
	m_dim_count=dim_count;
	for (int j=0; j<dim_count; j++) m_dims[j]=dims[j];
}
void Mda::reshape(qint32 S1,qint32 S2,qint32 S3,qint32 S4,qint32 S5,qint32 S6,qint32 S7,qint32 S8) {
	qint32 dims[MAX_MDA_DIMS];
	int num_dimensions=0;
	if (S1>=0) {
		num_dimensions++;
		dims[0]=S1;
	if (S2>=0) {
		num_dimensions++;
		dims[1]=S2;
	if (S3>=0) {
		num_dimensions++;
		dims[2]=S3;
	if (S4>=0) {
		num_dimensions++;
		dims[3]=S4;
	if (S5>=0) {
		num_dimensions++;
		dims[4]=S5;
	if (S6>=0) {
		num_dimensions++;
		dims[5]=S6;
	if (S7>=0) {
		num_dimensions++;
		dims[6]=S7;
	if (S8>=0) {
		num_dimensions++;
		dims[7]=S8;
	}}}}}}}}
	reshape(num_dimensions,dims);
}

#include <math.h>
#include <QTime>

real normsqr(const AbstractComplex &X) {
	return X.re()*X.re()+X.im()*X.im();
}
MDAComplex operator*(const AbstractComplex &X1, const AbstractComplex &X2) {
	return MDAComplex(X1.re()*X2.re()-X1.im()*X2.im(),X1.re()*X2.im()+X1.im()*X2.re());
}
MDAComplex operator/(const AbstractComplex &X1, const AbstractComplex &X2) {
	real ns=normsqr(X2);
	if (ns==0) return MDAComplex(0,0);
	else return X1*conj(X2)/ns;
}
MDAComplex operator+(const AbstractComplex &X1, const AbstractComplex &X2) {
	return MDAComplex(X1.re()+X2.re(),X1.im()+X2.im());
}
MDAComplex operator-(const AbstractComplex &X1, const AbstractComplex &X2)  {
	return MDAComplex(X1.re()-X2.re(),X1.im()-X2.im());
}
MDAComplex operator*(const AbstractComplex &X1, real X2) {
	return MDAComplex(X1.re()*X2,X1.im()*X2);
}
MDAComplex operator/(const AbstractComplex &X1, real X2) {
	if (X2==0) return MDAComplex(0,0);
	return MDAComplex(X1.re()/X2,X1.im()/X2);
}
MDAComplex operator+(const AbstractComplex &X1, real X2) {
	return MDAComplex(X1.re()+X2,X1.im());
}
MDAComplex operator-(const AbstractComplex &X1, real X2) {
	return MDAComplex(X1.re()-X2,X1.im());
}
MDAComplex operator*(real X1, const AbstractComplex &X2) {
	return X2*X1;
}
MDAComplex operator/(real X1, const AbstractComplex &X2) {
	return MDAComplex(X1,0)/X2;
}
MDAComplex operator+(real X1, const AbstractComplex &X2) {
	return MDAComplex(X1+X2.re(),X2.im());
}
MDAComplex operator-(real X1, const AbstractComplex &X2) {
	return MDAComplex(X1-X2.re(),-X2.im());
}
bool operator==(const AbstractComplex &X1, const AbstractComplex &X2) {
	return ((X1.re()==X2.re())&&(X1.im()==X2.im()));
}
bool operator==(const AbstractComplex &X1, real X2) {
	return ((X1.re()==X2)&&(X1.im()==0));
}
bool operator==(real X1, const AbstractComplex &X2) {
	return (X2==X1);
}
real abs(const AbstractComplex &X) {
	return (real)sqrt(normsqr(X));
}
MDAComplex conj(const AbstractComplex &X) {
	return MDAComplex(X.re(),-X.im());
}
MDAComplex exp(const AbstractComplex &X) {
	real r=(real)exp(X.re());
	real theta=X.im();
	return MDAComplex(r*(real)cos(theta),r*(real)sin(theta));
}
MDAComplex log(const AbstractComplex &X) {
	real theta=(real)atan2(X.im(),X.re());
	real r=(real)log(normsqr(X))/2;
	return MDAComplex(r,theta);
}
MDAComplex sqrt(const AbstractComplex &X) {
	if (X==0) return MDAComplex(0,0);
	return exp(log(X)*0.5);
}
MDAComplex power(const AbstractComplex &X1, const AbstractComplex &X2) {
	if (X1==0) return MDAComplex(0,0);
	return exp(log(X1)*X2);
}

MDAComplex operator-(const AbstractComplex &X) {
	return X*(-1);
}



