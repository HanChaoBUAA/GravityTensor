#include <math.h>
#include "BaseDefine.h"
#include "GravityModel.h"

CGravityData::CGravityData()
{
	Gm = 3.986004418E+14;
	RefDistance = 6378137.0;
	m_NN = 0;
	m_MM = 0;
	CeD = NULL;
	SeD = NULL;
}

CGravityData::~CGravityData()
{
	M_FREE(CeD);
	M_FREE(SeD);
	m_NN = 0;
	m_MM = 0;
}

void CGravityData::InitData(int DegreeNum)
{

}

double *CGravityData::GreateData(int N,int M,double *scD)
{
	int i;
	if(N!=m_NN){
		M_FREE(scD);
		m_NN = N;
	}
	if(scD==NULL){
		scD = (double*)malloc((GrIdx(N,N)+1)*sizeof(double));
	}
	if(scD!=NULL){
		for(i=0;i<=GrIdx(N,N);i++){
			scD[i] = 0.0;
		}
	}
	return scD;
}

double CGravityData::BetaFac(int n,int m)
{
	double Value=0;
	if(m<0||n<2){
		Value = 1.0;
	}
	else if(m==0){
		Value = sqrt(2*n+1.0);
	}
	else if(m<=n){
		int i;
		Value = 2*(2*n+1.0);
		Value = 1.0/Value;
		for(i=n+m;i>n-m;i--){
			Value=i*Value;
		}
		Value = sqrt(1.0/Value);
	}
	return Value;
}

int CGravityData::SetData(int n,int m,double ce,double se,int DataKind)
{
	if(CeD!=NULL&&SeD!=NULL&&n<=m_NN && n<=m_MM){
		if(DataKind ==GrvDataKind){
			CeD[GrIdx(n,m)] = ce;
			SeD[GrIdx(n,m)] = se;
			return 1;
		}
		else {
			return 0;
		}
	}
	return 0;
}

void CGravityData::CopyData(CGravityData*pGrvData)
{
	M_FREE(CeD);
	M_FREE(SeD);
	if(pGrvData!=NULL){
		int n,m;
		Gm			= pGrvData->Gm;
		RefDistance = pGrvData->RefDistance;

		m_NN = pGrvData->m_NN;
		m_MM = pGrvData->m_NN;
		CeD = GreateData(m_NN,m_MM,CeD);
		SeD = GreateData(m_NN,m_MM,SeD);

		for(n=2;n<=m_NN;n++){
			for(m=0;m<=min(m_MM,n);m++){
				SetData(n,m,pGrvData->CeD[GrIdx(n,m)],pGrvData->SeD[GrIdx(n,m)],pGrvData->GrvDataKind);
			}
		}
	}
	else{
		m_NN = 0;
		m_MM = 0;
	}
}

void CGravityData::CloneData(CGravityData*pGrvData,int Status)
{
	if(Status==0){
		M_FREE(CeD);
		M_FREE(SeD);
		if(pGrvData!=NULL){
			Gm = pGrvData->Gm;
			RefDistance = pGrvData->RefDistance;
			m_NN = pGrvData->m_NN;
			m_MM = pGrvData->m_NN;
			CeD = GreateData(m_NN,m_MM,CeD);
			SeD = GreateData(m_NN,m_MM,SeD);
			memcpy(CeD,pGrvData->CeD,(GrIdx(m_NN,m_NN)+1)*sizeof(double));
			memcpy(SeD,pGrvData->SeD,(GrIdx(m_NN,m_NN)+1)*sizeof(double));
		}
		else{
			m_NN = 0;
			m_MM = 0;
		}
	}
	else if(pGrvData!=NULL){
		memcpy(CeD,pGrvData->CeD,(GrIdx(m_NN,m_NN)+1)*sizeof(double));
		memcpy(SeD,pGrvData->SeD,(GrIdx(m_NN,m_NN)+1)*sizeof(double));
	}
}

CGravityNormalData::CGravityNormalData()
{
	GrvDataKind = ID_GraData_Normal;
}

CGravityNormalData::CGravityNormalData(CGravityData*pGrvData)
{
	GrvDataKind = ID_GraData_Normal;
	CopyData(pGrvData);
}