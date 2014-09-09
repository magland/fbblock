#include "mda_io.h"
#include "mda.h"

bool read_mda(FBArray4D<unsigned char> &X,QString fname) {
	Mda tmp;
	if (!tmp.read(fname.toAscii().data())) return false;
	X.allocate(tmp.size(0),tmp.size(1),tmp.size(2),tmp.size(3));
	for (long ii=0; ii<tmp.size(); ii++) {
		X.setValue1(tmp.get(ii).re(),ii);
	}
	return true;
}
bool read_mda(FBArray3D<float> &X,QString fname) {
	Mda tmp;
	if (!tmp.read(fname.toAscii().data())) return false;
	X.allocate(tmp.size(0),tmp.size(1),tmp.size(2));
	for (long ii=0; ii<tmp.size(); ii++) {
		X.setValue1(tmp.get(ii).re(),ii);
	}
	return true;
}
bool write_mda(FBSparseArray4D &X,QString fname) {
	Mda tmp;
	tmp.allocate(MDA_TYPE_REAL,X.N2(),X.N3(),X.N4(),X.N1());
	for (long z=0; z<X.N4(); z++)
	for (long y=0; y<X.N3(); y++)
	for (long x=0; x<X.N2(); x++)
	for (long w=0; w<X.N1(); w++) {
		tmp(x,y,z,w)=X.value(w,x,y,z);
	}
	return tmp.write(fname.toAscii().data());
}
bool read_mda(FBSparseArray4D &X,QString fname) {
	Mda tmp;
	if (!tmp.read(fname)) return false;
	X.allocate(DATA_TYPE_FLOAT,tmp.size(3),tmp.size(0),tmp.size(1),tmp.size(2));
	for (int pass=1; pass<=3; pass++) {
		for (int i4=0; i4<X.N4(); i4++)
		for (int i3=0; i3<X.N3(); i3++)
		for (int i2=0; i2<X.N2(); i2++)
		for (int i1=0; i1<X.N1(); i1++) {
			if (tmp.get(i2,i3,i4,i1).re()) {
				if (pass<=2) X.setupIndex(pass,i1,i2,i3,i4);
				else X.setValue(tmp.get(i2,i3,i4,i1).re(),i1,i2,i3,i4);
			}
		}
	}
	return true;
}
bool write_mda(FBArray4D<float> &X,QString fname) {
	Mda tmp;
	tmp.allocate(MDA_TYPE_REAL,X.N1(),X.N2(),X.N3(),X.N4());
	long NN=tmp.size();
	for (long ii=0; ii<NN; ii++) {
		tmp(ii)=X.value1(ii);
	}
	return tmp.write(fname.toAscii().data());
}
bool write_mda(FBArray3D<unsigned char> &X,QString fname) {
	Mda tmp;
	tmp.allocate(MDA_TYPE_REAL,X.N1(),X.N2(),X.N3());
	long NN=tmp.size();
	for (long ii=0; ii<NN; ii++) {
		tmp(ii)=X.value1(ii);
	}
	return tmp.write(fname.toAscii().data());
}

