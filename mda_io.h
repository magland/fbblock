#ifndef mda_io_H
#define mda_io_H

#include "arrays.h"
#include "fbsparsearray4d.h"

bool write_mda(FBArray3D<unsigned char> &X,QString fname);
bool write_mda(FBArray4D<float> &X,QString fname);
bool read_mda(FBArray4D<unsigned char> &X,QString fname);
bool read_mda(FBArray3D<float> &X,QString fname);
bool write_mda(FBSparseArray4D &X,QString fname);
bool read_mda(FBSparseArray4D &X,QString fname);

#endif
