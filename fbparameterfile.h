#ifndef fbparameterfile_H
#define fbparameterfile_H

#include <QHash>
#include <QString>
#include <QStringList>
#include "fbblocksolver.h"

class FBParameterFile {
public:
	FBParameterFile();
	virtual ~FBParameterFile();
	void read(QString fname);
	QString getString(QString key);
	int getInteger(QString key);
	fbreal getReal(QString key);
	QList<fbreal> getArray(QString key);
	QStringList keys() {return m_data.keys();}
private:
	QHash<QString,QString> m_data;
};

#endif
