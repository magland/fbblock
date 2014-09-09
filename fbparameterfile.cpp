#include "fbparameterfile.h"
#include <QStringList>
#include <QRegExp>
#include <QFile>
#include <QDebug>

FBParameterFile::FBParameterFile() 
{
}

FBParameterFile::~FBParameterFile()
{
}

void FBParameterFile::read(QString fname) {
	m_data.clear();
	QString txt=read_text_file(fname);
	if (txt.isEmpty()) return;
	
	QStringList lines=txt.split("\n");
	foreach (QString line,lines) {
		int ind1=line.indexOf("#");
		if (ind1>=0) line=line.mid(0,ind1);
		int ind=line.indexOf("=");
		if (ind>=0) {
			QString key=line.mid(0,ind).trimmed();
			QString val=line.mid(ind+1).trimmed();
			m_data[key]=val;
		}
	}
}

QString FBParameterFile::getString(QString key) {
	if (m_data.contains(key)) return m_data[key];
	else return "";
}
int FBParameterFile::getInteger(QString key) {
	if (m_data.contains(key)) return m_data[key].toInt();
	else return 0;
}
fbreal FBParameterFile::getReal(QString key) {
	if (m_data.contains(key)) return m_data[key].toDouble();
	else return 0;
}
QList<fbreal> FBParameterFile::getArray(QString key) {
	QList<fbreal> ret;
	QString str=getString(key);
	QStringList list=str.split(QRegExp("\\s+"));
	foreach (QString word,list) {
		ret << word.toDouble();
	}
	return ret;
}

