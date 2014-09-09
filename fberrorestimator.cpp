#include "fberrorestimator.h"
#include <math.h>
#include <QDebug>

typedef QList<double> DoubleList;

class FBErrorEstimatorPrivate {
public:
	FBErrorEstimator *q;
	QList<DoubleList> m_stress_records;
	QList<double> m_estimated_relative_errors;
};

FBErrorEstimator::FBErrorEstimator() 
{
	d=new FBErrorEstimatorPrivate;
	d->q=this;
}

FBErrorEstimator::~FBErrorEstimator()
{
	delete d;
}
void FBErrorEstimator::addStressData(const QList<double> &stress) {
	d->m_stress_records << stress;
}
void estimate_slope_intercept(double &slope,double &intercept,const QList<double> &data) {
	int N=data.count();
	double Sx=0;
	double Sy=0;
	double Sxx=0;
	double Sxy=0;
	for (int i=0; i<N; i++) {
		Sx+=i;
		Sy+=data[i];
		Sxx+=i*i;
		Sxy+=i*data[i];
	}
	double numer=(N*Sxy-Sx*Sy);
	double denom=(N*Sxx-Sx*Sx);
	if (!denom) {
		slope=intercept=0;
		return;
	}
	else {
		slope=numer/denom;
		intercept=Sy/N-slope*Sx/N;
	}
}

double compute_estimated_relative_error(const QList<DoubleList> &stress_records,int iteration) {	
	int min_iterations_to_use=30;
	int num_iterations_to_use=30;
	if (iteration+1<min_iterations_to_use) return 1;
	QList<double> last_stress=stress_records[iteration];
	int max_index=0;
	double max_val=0;
	for (int ii=0; ii<last_stress.count(); ii++) {
		if (qAbs(last_stress[ii])>max_val) {
			max_val=qAbs(last_stress[ii]);
			max_index=ii;
		}
	}
	
	//We assume that Y[j] = a + b*exp(-c*j)
	QList<double> Y; //vector of stresses
	for (int ii=qMax(iteration+1-num_iterations_to_use,0); ii<iteration+1; ii++) {
		Y << (stress_records[ii])[max_index];
	}
	//compute the derivative: Y_prime[j] = -b*c*exp(-c*j)
	QList<double> Y_prime;
	for (int ii=0; ii<Y.count()-1; ii++) {
		Y_prime << Y[ii+1]-Y[ii];
	}
	//take the log of the absolute value: log_Y_prime[j] = log(|b*c|) - c*j
	QList<double> log_Y_prime;
	for (int ii=0; ii<Y_prime.count(); ii++) {
		if (qAbs(Y_prime[ii])) log_Y_prime << log(qAbs(Y_prime[ii]));
		else log_Y_prime << -50;
	}
	//Now we fit log_Y_prime to a line
	double slope,intercept;
	estimate_slope_intercept(slope,intercept,log_Y_prime);
	//We solve c = -slope, and log(abs(b)) = intercept-log(abs(c))
	double c=-slope; if (c<=0) return 1;
	double log_abs_b=intercept-log(qAbs(c));
	//the estimated error is abs(b)*exp(-c*N)
	//so the estimated log error is log|b| - c*N
	double est_log_error=log_abs_b-c*(log_Y_prime.count()-1);
	//the log  relative error is this minus log(|Y[N]|)
	double log_rel_error=est_log_error-log(qAbs(Y[Y.count()-1]));
	if (log_rel_error<=-20) return exp(-20);
	if (log_rel_error>=0) return 1;
	return exp(log_rel_error); //we avoid computing exp until very end
}

double FBErrorEstimator::estimatedRelativeError(int iteration) const {
	if (iteration<0) iteration=d->m_stress_records.count()-1;
	if (iteration>=d->m_stress_records.count()) return 1;
	while (iteration>=d->m_estimated_relative_errors.count()) {
		float val0=compute_estimated_relative_error(d->m_stress_records,iteration);
		d->m_estimated_relative_errors << val0;
	}
	if ((iteration>=0)&&(iteration<d->m_estimated_relative_errors.count()))
		return d->m_estimated_relative_errors[iteration];
	else
		return 1;
}
QList<double> FBErrorEstimator::stressData(int iteration) {
	return d->m_stress_records.value(iteration);
}

