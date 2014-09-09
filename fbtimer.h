#ifndef fbtimer_H
#define fbtimer_H

#include <QString>

class FBTimerPrivate;
class FBTimer {
public:
	friend class FBTimerPrivate;
private:
	FBTimer();
public:
	static FBTimer *instance();
	virtual ~FBTimer();
	static void startTimer(QString timer_name);
	static void stopTimer(QString timer_name);
	static double elapsed(QString timer_name);
	static QStringList timerNames();
private:
	FBTimerPrivate *d;
};

#endif
