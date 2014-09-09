#ifndef nonlinearadjuster_H
#define nonlinearadjuster_H

class NonlinearAdjusterPrivate;
class NonlinearAdjuster {
public:
	friend class NonlinearAdjusterPrivate;
	NonlinearAdjuster();
	virtual ~NonlinearAdjuster();
	virtual float computeAdjustment(float eps)=0;
private:
	NonlinearAdjusterPrivate *d;
};

#endif
