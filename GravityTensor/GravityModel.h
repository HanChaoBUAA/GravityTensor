#ifndef __GRAVITYMODEL_H__
#define __GRAVITYMODEL_H__

#include <malloc.h>
#include <string.h>

#define ID_GraData_None     0
#define ID_GraData_Normal   1

class CGravityData
{
public:
	int GrvDataKind;
	int m_NN,m_MM;
	double *CeD,*SeD;
	double Gm,RefDistance;

public:
	CGravityData();
	~CGravityData();

	void CopyData(CGravityData*pGrvData);
	void CloneData(CGravityData*pGrvData,int Status=0);

public:
	double *GreateData(int N,int M,double *);

	void InitData(int DegreeNum);

	int SetData(int n,int m,double ce,double se,int DataKind);
	static double BetaFac(int n,int m);

	virtual double *GetCe(){return CeD;}
	virtual double *GetSe(){return SeD;}

};


class CGravityNormalData:public CGravityData
{
public:

	CGravityNormalData();
	CGravityNormalData(CGravityData*pGrvData);

	int LoadData(char*FileName);

};

class CGrNormalData_WGS84_EGM96:public CGravityNormalData
{
public:
	CGrNormalData_WGS84_EGM96();

};

extern CGrNormalData_WGS84_EGM96 GrData_WGS84_EGM96;

#endif
