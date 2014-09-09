#include "nonlinearadjuster.h"

class NonlinearAdjusterPrivate {
public:
	NonlinearAdjuster *q;
};

NonlinearAdjuster::NonlinearAdjuster() 
{
	d=new NonlinearAdjusterPrivate;
	d->q=this;
}

NonlinearAdjuster::~NonlinearAdjuster()
{
	delete d;
}

