#ifndef __GRAVITYVW_H__
#define __GRAVITYVW_H__

#include "ObjectVector.h"
#include "GravityModel.h"

class CSphericalHNormalFactor
{
public:
	virtual double Beta(int n,int m){return 1.0;}
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q){return 1.0;}
};

class CGravNormalFactor:public CSphericalHNormalFactor
{
public:
	virtual double Beta(int n,int m);
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q);
};

class CBelikovFactor:public CSphericalHNormalFactor
{
public:
	virtual double Beta(int n,int m);
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q);
};

class CLinLin2Factor:public CSphericalHNormalFactor
{
public:
	virtual double Beta(int n,int m);
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q);
};

class CHanChao1Factor:public CSphericalHNormalFactor
{
public:
	virtual double Beta(int n,int m);
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q);
};

class CHanChao2Factor:public CSphericalHNormalFactor
{
public:
	double K,L;
	CHanChao2Factor(){K=1,L=1;}
	virtual double Beta(int n,int m);
	virtual double BetaB(int n,int m);
	virtual double PPi(int n,int m,int p,int q);
	virtual void SetFactor(double K_in,double L_in);
};

extern CSphericalHNormalFactor NoFactor;
extern CGravNormalFactor       FullNormalFactor;
extern CBelikovFactor          BelikovFactor;
extern CLinLin2Factor          LinLin2Factor;

class CSphericalHTensorData
{
private:
	CSphericalHNormalFactor *Fac;
public:
	int DegreeNumber,OrderNumber,TensorOrder;
	int DxN,DyN,DzN;
	double *CeD,*SeD;
	double Re,GMR;

public:

	CSphericalHTensorData();
	~CSphericalHTensorData();

	CSphericalHTensorData(CGravityData*pGrvData,int DegreeNum=-1,int OrderNum=-1,CSphericalHNormalFactor *pFac=NULL,int FacID=0);

	static CSphericalHTensorData * CreateCoefData(CGravityData*pGrvData, int DegreeNum=-1, int OrderNum=-1,CSphericalHTensorData *pSCxyz=NULL,CSphericalHNormalFactor *pFac=NULL,int FacID=0);

	static CSphericalHTensorData * CreateCoefData(CSphericalHTensorData*pGrvData, const char Axis='n', int DegreeNum=-1, int OrderNum=-1, CSphericalHTensorData *pSCxyz=NULL);

	void SetFactor(CSphericalHNormalFactor* pFac);

protected:

	static double * InitCeSeData(double *CSeD,int DegreeNum);
	static void SetTensorCoefData(CSphericalHTensorData*pGrvData,
		int DegreeNum,int OrderNum,
		CSphericalHTensorData *pSCxyz,
		const char Axis);

	void SetGravityCoefData(CGravityData*pGrvData,
		int Degree,int Order,
		CSphericalHNormalFactor *pFac=NULL,int FacID=0);

};


class CSphericalH
{
private:
	CSphericalHNormalFactor *Fac;

	static double ComputeTensorValue(CSphericalHTensorData *pTensor,double *VV,double *WW,
		int desiredDegree=-1,int desiredOrder=-1);

	static int ComputeVW(double *FixedPos,int desiredDegree,int desiredOrder,double R_ref,
		double *VV,double *WW,
		int CurrDegree,int CurrOrder,
		CSphericalHNormalFactor *pRef);
protected:
	int	CurrStatus;
	int VW_Degree,VW_Order,VMCurrNN,VMCurrMM;
	double *VV,*WW;

	double CurrPos[3];

	static double * SetVWSpace(double *VW,int CurrDegree,int desiredDegree);

	virtual void SetFactor(CSphericalHNormalFactor* pFac);

	int SetDegreeOrder(int GrvDegree,int GrvOrder);
	int InitSphHData(int Degree);

public:
	CSphericalH();
	~CSphericalH();

	int    SetPosition(double Pos[3]);
	double ComputeTensorValue(CSphericalHTensorData *pTensor,int desiredDegree=-1,int desiredOrder=-1);

};

typedef CObjVec<CSphericalHTensorData,20> CGravityTensorData;

class CGravityTensor:public CSphericalH
{
protected:
	int Status;
	int GrvFacID;
	CGravityData *pGrvData;
	CSphericalHNormalFactor *GrvFac;
	CGravityTensorData TensorDataSet;
	int TensorOrder,GrvDegree,GrvOrder;

	static int  TensorNumber(int N);
	static int  TensorIndex(int N,int a,int b);
	static int  DeFacIndex(int Index,int N,int *a,int *b,int *c);
	static char GetParent(int Index,int N,int *a,int *b,int *c);

public:
	CGravityTensor();
	~CGravityTensor();

	CGravityTensor(CGravityData *pData,CSphericalHNormalFactor *pFac=NULL,int FacID=0);
	virtual void SetFactor(CSphericalHNormalFactor* pFac);

	int SetGrData(CGravityData *pData,CSphericalHNormalFactor *pFac=NULL,int FacID=0);
	int SetDegreeOrder(int TensorOrder,int GrvDegree,int GrvOrder);

	int SetPosition(double *Pos);
	int GetTensorValue(int a,int b,int c, double *Value);

};

#endif
