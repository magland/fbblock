#include "fbtimer.h"
#include <QTime>
#include <QHash>
#include <QStringList>

struct TimerData {
	long ms_elapsed;
	QTime time;	
};

class FBTimerPrivate {
public:
	FBTimer *q;
	QHash<QString,TimerData> m_timers;
};

FBTimer::FBTimer() 
{
	d=new FBTimerPrivate;
	d->q=this;
}

FBTimer::~FBTimer()
{
	delete d;
}
void FBTimer::startTimer(QString timer_name) {
	FBTimerPrivate *dd=instance()->d;
	if (!dd->m_timers.contains(timer_name)) {
		TimerData TD;
		TD.ms_elapsed=0;
		dd->m_timers[timer_name]=TD;
	}
	dd->m_timers[timer_name].time.start();
}
void FBTimer::stopTimer(QString timer_name) {
	FBTimerPrivate *dd=instance()->d;
	if (!dd->m_timers.contains(timer_name)) return;
	dd->m_timers[timer_name].ms_elapsed+=dd->m_timers[timer_name].time.restart();
}
FBTimer *FBTimer::instance() {
	static FBTimer *ret=0;
	if (!ret) ret=new FBTimer;
	return ret;
}
double FBTimer::elapsed(QString timer_name) {
	FBTimerPrivate *dd=instance()->d;
	if (!dd->m_timers.contains(timer_name)) return 0;
	return dd->m_timers[timer_name].ms_elapsed/1000.0;
}
QStringList FBTimer::timerNames() {
	FBTimerPrivate *dd=instance()->d;
	QStringList ret=dd->m_timers.keys();
	qSort(ret);
	return ret;
}

