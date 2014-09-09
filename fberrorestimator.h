#ifndef fberrorestimator_H
#define fberrorestimator_H

#include <QList>

class FBErrorEstimatorPrivate;
class FBErrorEstimator {
public:
	friend class FBErrorEstimatorPrivate;
	FBErrorEstimator();
	virtual ~FBErrorEstimator();
	void addStressData(const QList<double> &stress);
	double estimatedRelativeError(int iteration=-1) const;
	QList<double> stressData(int iteration);
private:
	FBErrorEstimatorPrivate *d;
};

#endif
