#include <math.h>
#include <string.h>
#include "BaseDefine.h"
#include "GravityVW.h"

double * zero_v(double *v)
{
	v[1] = 0;
	v[2] = 0;
	v[3] = 0;
	return v+1;
}

int isequal_v(double a[3], double b[3])
{
	if(a==b){
		return 1;
	}
	else{
		if(memcmp(a,b,3*sizeof(double))==0){
			return 1;
		}
		else{
			return 0;
		}
	}
}

double dot(double *a, double *b)
{
	double ret_val;

	--b;
	--a;

	ret_val = a[1] * b[1] + a[2] * b[2] + a[3] * b[3];

	return ret_val;
}

CSphericalHNormalFactor NoFactor;
CGravNormalFactor       FullNormalFactor;
CBelikovFactor          BelikovFactor;
CLinLin2Factor          LinLin2Factor;

double CSphericalHNormalFactor::BetaB(int n,int m)
{
	double Value=1;
	if(m<0||n<0){
		Value = 1.0;
	}
	else if(m>n){
		Value = 1.0;
	}
	else if(m==0){
		Value = sqrt(2*n+1.0);
	}
	else if(m<=n){
		int i;
		Value = 2*(2*n+1.0);
		Value = 1.0/Value;
		for(i=n-m+1;i<=n+m;i++){
			Value=i*Value;
		}
		Value = 1.0/sqrt(Value);
	}
	return Value;
}

double CGravNormalFactor::Beta(int n,int m)
{
	double Value=1;
	if(m<0||n<0){
		Value = 1.0;
	}
	else if(m>n){
		Value = 1.0;
	}
	else if(m==0){
		Value = 1.0/sqrt(2*n+1.0);
	}
	else if(m<=n){
		int i;
		Value = 2*(2*n+1.0);
		Value = 1.0/Value;
		for(i=n-m+1;i<=n+m;i++){
			Value=i*Value;
		}
		Value = sqrt(Value);
	}
	return Value;
}

double CGravNormalFactor::BetaB(int n,int m)
{
	return 1;
}

double CGravNormalFactor::PPi(int n,int m,int p,int q)
{
	double Value=1;
	if(p<0||n<0||q<0||m<0){
		Value = 1.0;
	}
	else if(p<q||n<m){
		Value = 1.0;
	}
	else{
		int i;
		Value = (2*n+1.0)/(2*p+1.0);
		if(m==q){
		}
		else{
			if(m!=0){
				Value = Value * 2;
			}
			if(q!=0){
				Value = Value * 0.5;
			}
		}
		if(p+q<n+m){
			for(i=p+q+1;i<=n+m;i++){
				Value = Value / i;
			}
		}
		else if(p+q>n+m){
			for(i=n+m+1;i<=p+q;i++){
				Value = Value * i;
			}
		}
		if(p-q<n-m){
			for(i=p-q+1;i<=n-m;i++){
				Value = Value * i;
			}
		}
		else if(p-q>n-m){
			for(i=n-m+1;i<=p-q;i++){
				Value = Value / i;
			}
		}
		Value = sqrt(Value);
	}
	return Value;
}

double CBelikovFactor::Beta(int n,int m)
{
	double Value=1;
	if(m<0||n<1){
		Value = 1.0;
	}
	else{
		int i;
		Value = 1.0/pow(2.0,n);
		for(i=n+1;i<=n+m;i++){
			Value=i*Value;
		}
	}
	return Value;
}

double CBelikovFactor::BetaB(int n,int m)
{
	double Value=1;
	int i;

	if(n<0||m<0){
		Value = 1;
	}
	else if(n<m){
		Value = 1;
	}
	else if(m>0){
		Value = 1;
		for(i=n+1;i<=n+m;i++){
			Value = Value*i;
		}
		Value = 1.0/Value;
		for(i=n-m+1;i<=n;i++){
			Value = Value*i;
		}
		Value = 1.0/Value;
		Value = sqrt(2*(2*n+1)*Value);
	}
	else{
		Value = sqrt(2*n+1.0);
	}
	Value = Value/pow(2.0,n);

	//    Value = Beta(n,m)/FullNormalFactor.Beta(n,m);

	return Value;
}

double CBelikovFactor::PPi(int n,int m,int p,int q)
{
	double Value=1;
	int i;

	if(p<0||n<0||q<0||m<0){
		Value = 1.0;
	}
	else if(p<q||n<m){
		Value = 1.0;
	}
	else{
		if(n>p){
			Value = pow(2.0,n-p);
			for(i=p+1;i<=n;i++){
				Value = Value*i;
			}
		}
		else if(n<p){
			Value = pow(2.0,p-n);
			for(i=n+1;i<=p;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
		if(n+m>p+q){
			Value = 1.0/Value;
			for(i=p+q+1;i<=n+m;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
		else if(n+m<p+q){
			for(i=n+m+1;i<=p+q;i++){
				Value = Value*i;
			}
		}
	}
	return Value;
}

double CLinLin2Factor::Beta(int n,int m)
{
	double Value=1.0;
	if(m<0||n<1||n<m){
		Value = 1.0;
	}
	else{
		int i;
		Value = 1.0;
		for(i=1;i<=n-m;i++){
			Value=i*Value;
		}
		Value = 1.0/Value;
		for(i=1;i<=n;i++){
			Value=(2*i-1)*Value;
		}
	}
	return Value;
}

double CLinLin2Factor::BetaB(int n,int m)
{
	double Value=1;

	Value = Beta(n,m)/FullNormalFactor.Beta(n,m);

	return Value;
}

double CLinLin2Factor::PPi(int n,int m,int p,int q)
{
	double Value=1.0;
	int i;

	if(p<0||n<0||q<0||m<0){
		Value = 1.0;
	}
	else if(p<q||n<m){
		Value = 1.0;
	}
	else{
		Value = 1.0;
		if(n>p){
			for(i=p+1;i<=n;i++){
				Value = Value*(2*i-1);
			}
			Value = 1.0/Value;
		}
		else if(n<p){
			for(i=n+1;i<=p;i++){
				Value = Value*(2*i-1);
			}
		}
		if(n-m>p-q){
			for(i=p-q+1;i<=n-m;i++){
				Value = Value*i;
			}
		}
		else if(n-m<p-q){
			Value = 1.0/Value;
			for(i=n-m+1;i<=p-q;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
	}
	return Value;
}

double CHanChao1Factor::Beta(int n,int m)
{
	double Value=1;
	if(m<0||n<0){
		Value = 1.0;
	}
	else if(m>n){
		Value = 1.0;
	}
	else if(m==0){
		Value = 1.0/(2*n+1.0);
	}
	else if(m<=n){
		int i;
		Value = (2*n+1.0);
		Value = 1.0/Value;
		for(i=n-m+1;i<=n+m;i++){
			Value=i*Value;
		}
	}
	return Value;
}

double CHanChao1Factor::BetaB(int n,int m)
{
	double Value=1;

	Value = Beta(n,m)/FullNormalFactor.Beta(n,m);

	return Value;
}

double CHanChao1Factor::PPi(int n,int m,int p,int q)
{
	double Value=1.0;
	int i;

	if(p<0||n<0||q<0||m<0){
		Value = 1.0;
	}
	else if(p<q||n<m){
		Value = 1.0;
	}
	else{
		Value = 1.0;
		if(n!=p){
			Value = (2*n+1.0)/(2*p+1);
		}
		if(n+m>p+q){
			Value = 1.0/Value;
			for(i=p+q+1;i<=n+m;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
		else if(n+m<p+q){
			for(i=n+m+1;i<=p+q;i++){
				Value = Value*i;
			}
		}
		if(n-m>p-q){
			for(i=p-q+1;i<=n-m;i++){
				Value = Value*i;
			}
		}
		else if(n-m<p-q){
			Value = 1.0/Value;
			for(i=n-m+1;i<=p-q;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
	}
	return Value;
}

void CHanChao2Factor::SetFactor(double K_in,double L_in)
{
	if(K>0){
		K=K_in;
	}
	if(L>0){
		L=L_in;
	}
}

double CHanChao2Factor::Beta(int n,int m)
{
	double Value=1;
	if(m<0||n<0){
		Value = 1.0;
	}
	else if(m>n){
		Value = 1.0;
	}
	else if(m<=n){
		int i;
		Value = pow(K,n);
		if(m!=0){
			Value = Value * pow(L,m);
		}	
		for(i=n-m;i<n+m;i++){
			Value = Value * i;
		}	
	}
	return Value;
}

double CHanChao2Factor::BetaB(int n,int m)
{
	double Value=1;

	Value = Beta(n,m)/FullNormalFactor.Beta(n,m);

	return Value;
}

double CHanChao2Factor::PPi(int n,int m,int p,int q)
{
	double Value=1.0;
	int i=0;

	if(p<0||n<0||q<0||m<0){
		Value = 1.0;
	}
	else if(p<q||n<m){
		Value = 1.0;
	}
	else{
		int i;
		Value = 1.0;
		if(p>n){
			for(i=n+1;i<p;i++){
				Value = Value * K * i;
			}
		}	
		else if(p<n){
			for(i=p+1;i<n;i++){
				Value = Value * K * i;
			}
			Value = 1.0/Value;
		}	
		if(q-m!=0){
			Value = Value * pow(L,q-m);
		}	
		if(n-m<p-q){
			for(i=n-m+1;i<=p-q;i++){
				Value = Value*i;
			}
		}
		else if(n-m>p-q){
			Value = 1.0/Value;
			for(i=p-q+1;i<=n-m;i++){
				Value = Value*i;
			}
			Value = 1.0/Value;
		}
	}
	return Value;
}

CSphericalHTensorData::CSphericalHTensorData()
{
	Fac = &NoFactor;
	CeD = NULL;
	SeD = NULL;
	DegreeNumber = 0;
	OrderNumber  = 0;
	TensorOrder  = 0;
	DxN = 0;
	DyN = 0;
	DzN = 0;
	Re  = 1.0;
	GMR = 1.0;
}

CSphericalHTensorData::~CSphericalHTensorData()
{
	Fac = NULL;
	DegreeNumber = 0;
	TensorOrder  = 0;
	DxN = 0;
	DyN = 0;
	DzN = 0;
	M_FREE(CeD);
	M_FREE(SeD);
	Re  = 1.0;
	GMR = 1.0;
}

CSphericalHTensorData::CSphericalHTensorData(CGravityData*pData,int Degree,int Order,CSphericalHNormalFactor *pFac,int FacID)
{
	CeD = NULL;
	SeD = NULL;
	SetGravityCoefData(pData,Degree,Order,pFac,FacID);
}

void CSphericalHTensorData::SetFactor(CSphericalHNormalFactor* pFac)
{
	Fac = pFac;
}

double * CSphericalHTensorData::InitCeSeData(double *CSeD,int DegreeNum)
{
	int i;
	M_FREE(CSeD);
	if(DegreeNum>0){
		CSeD = (double*)malloc((GrIdx(DegreeNum,DegreeNum)+1)*sizeof(double));
		for(i=0;i<GrIdx(DegreeNum,DegreeNum)+1;i++){
			CSeD[i] = 0.0;
		}
	}
	return CSeD;
}

void CSphericalHTensorData::SetGravityCoefData(CGravityData*pData,int Degree,int Order,CSphericalHNormalFactor *pFac,int FacID)
{
	int n,m;

	if(pData==NULL){
		DegreeNumber = 0;
	}
	CeD = InitCeSeData(CeD,DegreeNumber);
	SeD = InitCeSeData(SeD,DegreeNumber);
	if(DegreeNumber == 0){
		return ;
	}

	if(Degree<0){
		Degree = pData->m_NN;
		Order  = min(Degree,pData->m_MM);
	}
	else{
		Degree = min(Degree,pData->m_NN);
		Degree = max(Degree,2);
		Order = min(min(Order,Degree),pData->m_MM);
		Order = max(Order,0);
	}

	DegreeNumber = Degree;
	OrderNumber  = Order;
	TensorOrder  = 0;

	DxN = 0;
	DyN = 0;
	DzN = 0;

	CeD = InitCeSeData(CeD,DegreeNumber);
	SeD = InitCeSeData(SeD,DegreeNumber);

	if(pFac!=NULL){
		Fac=pFac;
	}
	if(FacID){
		double Beta;
		for(n=2;n<=DegreeNumber;n++){
			for(m=0;m<=min(n,OrderNumber);m++){
				Beta = Fac->BetaB(n,m);
				CeD[GrIdx(n,m)] = pData->CeD[GrIdx(n,m)]*Beta;
				SeD[GrIdx(n,m)] = pData->SeD[GrIdx(n,m)]*Beta;
			}
		}
	}
	else{
		for(n=2;n<=DegreeNumber;n++){
			for(m=0;m<=min(n,OrderNumber);m++){
				CeD[GrIdx(n,m)] = pData->CeD[GrIdx(n,m)];
				SeD[GrIdx(n,m)] = pData->SeD[GrIdx(n,m)];
			}
		}
	}

	GMR = pData->Gm/pData->RefDistance;
	Re  = pData->RefDistance;

}

void CSphericalHTensorData::SetTensorCoefData(CSphericalHTensorData*pData,
											  int DegreeNum,int OrderNum,
											  CSphericalHTensorData *pSCxyz,
											  const char Axis)
{
	int n,m,s;
	double Fac1,Fac2;
	CSphericalHNormalFactor *pFac;

	if(pData==NULL||pSCxyz==NULL){
		return;
	}

	if(DegreeNum<1){
		DegreeNum = pData->DegreeNumber;
	}
	DegreeNum = min(pData->DegreeNumber,DegreeNum);
	DegreeNum = max(2,DegreeNum);

	if(OrderNum<1){
		OrderNum = pData->OrderNumber;
	}
	OrderNum = min(pData->OrderNumber,OrderNum);
	OrderNum = max(2,OrderNum);

	if(Axis=='n'){
		pSCxyz->CeD = InitCeSeData(pSCxyz->CeD,DegreeNum);
		pSCxyz->SeD = InitCeSeData(pSCxyz->SeD,DegreeNum);
		pSCxyz->DegreeNumber = DegreeNum;
		pSCxyz->OrderNumber  = OrderNum;
		pSCxyz->TensorOrder  = pData->TensorOrder;
		pSCxyz->GMR = pData->GMR;
		pSCxyz->Re  = pData->Re;
		pSCxyz->Fac = pData->Fac;
	}
	else if(Axis=='x'||Axis=='y'||Axis=='z'){
		DegreeNum = DegreeNum + 1;
		OrderNum  = OrderNum + 1;
		pSCxyz->CeD = InitCeSeData(pSCxyz->CeD,DegreeNum);
		pSCxyz->SeD = InitCeSeData(pSCxyz->SeD,DegreeNum);
		pSCxyz->DegreeNumber = DegreeNum;
		pSCxyz->OrderNumber  = OrderNum;
		pSCxyz->TensorOrder  = pData->TensorOrder + 1;
		pSCxyz->GMR = pData->GMR/pData->Re;
		pSCxyz->Re  = pData->Re;
		pSCxyz->Fac = pData->Fac;
	}
	else{
		return;
	}

	if(Axis=='n'){
		for(n=2+pData->TensorOrder;n<=pSCxyz->DegreeNumber;n++){
			for(m=0;m<=min(n,OrderNum);m++){
				pSCxyz->CeD[GrIdx(n,m)] =pData->CeD[GrIdx(n,m)];
				pSCxyz->SeD[GrIdx(n,m)] =pData->SeD[GrIdx(n,m)];
			}
		}
		pSCxyz->DxN = pData->DxN;
		pSCxyz->DyN = pData->DyN;
		pSCxyz->DzN = pData->DzN;
		return;
	}

	pFac = pData->Fac;

	if(Axis=='x'){
		for(n=2+pSCxyz->TensorOrder;n<=pSCxyz->DegreeNumber;n++){
			Fac1 = pFac->PPi(n-1,1,n,0);
			pSCxyz->CeD[GrIdx(n,0)]		= 0.5*        n   *(n-1)*pData->CeD[GrIdx(n-1,1)]*Fac1;
			pSCxyz->SeD[GrIdx(n,0)]		= 0.0;
			if(OrderNum>1){
				Fac1 = pFac->PPi(n-1,2,n,1);
				Fac2 = pFac->PPi(n-1,0,n,1);
				pSCxyz->CeD[GrIdx(n,1)]		= 0.5*(  (n-1)*(n-2)*pData->CeD[GrIdx(n-1,2)]*Fac1
					-2*pData->CeD[GrIdx(n-1,0)]*Fac2);
				pSCxyz->SeD[GrIdx(n,1)]		= 0.5*   (n-1)*(n-2)*pData->SeD[GrIdx(n-1,2)]*Fac1;
			}
			else{
				Fac2 = pFac->PPi(n-1,0,n,1);
				pSCxyz->CeD[GrIdx(n,1)]		=                  - pData->CeD[GrIdx(n-1,0)]*Fac2;
				pSCxyz->SeD[GrIdx(n,1)]		= 0.0;
			}
			s = min(n,OrderNum);
			for(m=2;m<=s;m++){
				if(m+1>n-1){
					Fac2 = pFac->PPi(n-1,m-1,n,m);
					pSCxyz->CeD[GrIdx(n,m)] =-0.5*               pData->CeD[GrIdx(n-1,m-1)]*Fac2;
					pSCxyz->SeD[GrIdx(n,m)] =-0.5*               pData->SeD[GrIdx(n-1,m-1)]*Fac2;
				}
				else{
					Fac1 = pFac->PPi(n-1,m+1,n,m);
					Fac2 = pFac->PPi(n-1,m-1,n,m);
					pSCxyz->CeD[GrIdx(n,m)] = 0.5*((n-m)*(n-m-1)*pData->CeD[GrIdx(n-1,m+1)]*Fac1
						-pData->CeD[GrIdx(n-1,m-1)]*Fac2);
					pSCxyz->SeD[GrIdx(n,m)] = 0.5*((n-m)*(n-m-1)*pData->SeD[GrIdx(n-1,m+1)]*Fac1
						-pData->SeD[GrIdx(n-1,m-1)]*Fac2);
				}
			}
		}
		pSCxyz->DxN = pData->DxN+1;
		pSCxyz->DyN = pData->DyN;
		pSCxyz->DzN = pData->DzN;
	}
	else if(Axis=='y'){
		for(n=2+pSCxyz->TensorOrder;n<=pSCxyz->DegreeNumber;n++){
			Fac1 = pFac->PPi(n-1,1,n,0);
			pSCxyz->CeD[GrIdx(n,0)]		= 0.5*        n   *(n-1)*pData->SeD[GrIdx(n-1,1)]*Fac1;
			pSCxyz->SeD[GrIdx(n,0)]		= 0.0;
			if(OrderNum>1){
				Fac1 = pFac->PPi(n-1,2,n,1);
				Fac2 = pFac->PPi(n-1,0,n,1);
				pSCxyz->CeD[GrIdx(n,1)]		= 0.5*   (n-1)*(n-2)*pData->SeD[GrIdx(n-1,2)]*Fac1;
				pSCxyz->SeD[GrIdx(n,1)]		=-0.5*(  (n-1)*(n-2)*pData->CeD[GrIdx(n-1,2)]*Fac1
					+2.0*pData->CeD[GrIdx(n-1,0)]*Fac2);
			}
			else{
				Fac2 = pFac->PPi(n-1,0,n,1);
				pSCxyz->CeD[GrIdx(n,1)]		= 0.0;
				pSCxyz->SeD[GrIdx(n,1)]		=               -    pData->CeD[GrIdx(n-1,0)]*Fac2;
			}
			s = min(n,OrderNum);
			for(m=2;m<=s;m++){
				if(m+1>n-1){
					Fac2 = pFac->PPi(n-1,m-1,n,m);
					pSCxyz->CeD[GrIdx(n,m)] = 0.5*               pData->SeD[GrIdx(n-1,m-1)]*Fac2;
					pSCxyz->SeD[GrIdx(n,m)] =-0.5*               pData->CeD[GrIdx(n-1,m-1)]*Fac2;
				}
				else{
					Fac1 = pFac->PPi(n-1,m+1,n,m);
					Fac2 = pFac->PPi(n-1,m-1,n,m);
					pSCxyz->CeD[GrIdx(n,m)] = 0.5*((n-m)*(n-m-1)*pData->SeD[GrIdx(n-1,m+1)]*Fac1
						+pData->SeD[GrIdx(n-1,m-1)]*Fac2);
					pSCxyz->SeD[GrIdx(n,m)] =-0.5*((n-m)*(n-m-1)*pData->CeD[GrIdx(n-1,m+1)]*Fac1
						+pData->CeD[GrIdx(n-1,m-1)]*Fac2);
				}
			}
		}
		pSCxyz->DxN = pData->DxN;
		pSCxyz->DyN = pData->DyN+1;
		pSCxyz->DzN = pData->DzN;
	}
	else if(Axis=='z'){
		for(n=2+pSCxyz->TensorOrder;n<=pSCxyz->DegreeNumber;n++){
			s = min(n,OrderNum);
			m = s;
			pSCxyz->CeD[GrIdx(n,m)]		= 0.0;
			pSCxyz->SeD[GrIdx(n,m)]		= 0.0;
			for(m=0;m<s;m++){
				Fac1 = pFac->PPi(n-1,m,n,m);
				pSCxyz->CeD[GrIdx(n,m)] =			 -(n-m) *pData->CeD[GrIdx(n-1,m)]*Fac1;
				pSCxyz->SeD[GrIdx(n,m)] =			 -(n-m) *pData->SeD[GrIdx(n-1,m)]*Fac1;
			}
		}
		pSCxyz->DxN = pData->DxN;
		pSCxyz->DyN = pData->DyN;
		pSCxyz->DzN = pData->DzN+1;
	}

}

CSphericalHTensorData * CSphericalHTensorData::CreateCoefData(CSphericalHTensorData*pData,
															  const char Axis,
															  int DegreeNum,int OrderNum,
															  CSphericalHTensorData *pSCxyz)
{
	if(pData==NULL){
		return NULL;
	}

	if(pSCxyz==NULL){
		pSCxyz = new CSphericalHTensorData;
	}

	SetTensorCoefData(pData,DegreeNum,OrderNum, pSCxyz, Axis);

	return pSCxyz;

}

CSphericalHTensorData * CSphericalHTensorData::CreateCoefData(CGravityData*pData,
															  int DegreeNum,int OrderNum,
															  CSphericalHTensorData *pSCxyz,
															  CSphericalHNormalFactor *pFac,int FacID)
{
	if(pData==NULL){
		return NULL;
	}

	if(DegreeNum<1){
		DegreeNum = pData->m_NN;
	}
	DegreeNum = min(pData->m_NN,DegreeNum);
	DegreeNum = max(2,DegreeNum);

	if(pSCxyz==NULL){
		pSCxyz = new CSphericalHTensorData;
	}

	pSCxyz->SetGravityCoefData(pData,DegreeNum,OrderNum,pFac,FacID);

	return pSCxyz;

}

CSphericalH::CSphericalH()
{
	Fac = &NoFactor;

	VW_Degree = 0;
	VW_Order  = 0;

	VMCurrNN = 0;
	VMCurrMM = 0;

	VV = NULL;
	WW = NULL;

	zero_v(CurrPos);
	CurrStatus = 0;

}

CSphericalH::~CSphericalH()
{
	M_FREE(VV);
	M_FREE(WW);
}

void CSphericalH::SetFactor(CSphericalHNormalFactor* pFac)
{
	Fac = pFac;
}

int CSphericalH::SetDegreeOrder(int GrvDegree_in,int GrvOrder_in)
{
	if(GrvDegree_in>VW_Degree){
		VV = SetVWSpace(VV,VW_Degree,GrvDegree_in);
		WW = SetVWSpace(WW,VW_Degree,GrvDegree_in);
		VW_Degree = GrvDegree_in;
		VW_Order  = GrvOrder_in;
	}
	return 1;
}

int CSphericalH::SetPosition(double Pos[3])
{
	if(Pos==NULL){
		return 0;
	}
	if(!isequal_v(Pos,CurrPos)){
		memcpy(CurrPos,Pos,3*sizeof(double));
		CurrStatus = 1;
	}
	return 1;
}

double CSphericalH::ComputeTensorValue(CSphericalHTensorData *pData,int desiredDegree,int desiredOrder)
{
	if(desiredDegree<0){
		desiredDegree = pData->DegreeNumber;
	}
	if(desiredOrder<0){
		desiredOrder = pData->OrderNumber;
	}
	//	if(InitSphHData(desiredDegree)){
	if(SetDegreeOrder(desiredDegree,desiredOrder)){
		if(CurrStatus==1){
			VMCurrNN = 0;
			VMCurrMM = 0;
		}
		if(desiredDegree>VMCurrNN||desiredOrder>VMCurrMM){
			ComputeVW(CurrPos,desiredDegree,desiredOrder,pData->Re,
				VV,WW,VMCurrNN,VMCurrMM,Fac);
			VMCurrNN = max(VMCurrNN,desiredDegree);
			VMCurrMM = max(VMCurrMM,desiredOrder);
			CurrStatus = 0;
		}
		//		}
		return ComputeTensorValue(pData,VV,WW,desiredDegree,desiredOrder);
	}
	return 0;
}

double * CSphericalH::SetVWSpace(double *VW,int CurrDegree,int desiredDegree)
{
	int i;
	double * newVW=NULL;

	if(CurrDegree==desiredDegree){
		return VW;
	}

	if(VW==NULL){
		CurrDegree = 0;
	}

	desiredDegree = max(2,desiredDegree);
	desiredDegree = max(CurrDegree,desiredDegree);

	newVW = (double*)malloc((LegIdx(desiredDegree,desiredDegree)+1)*sizeof(double));

	if(VW!=NULL&&CurrDegree>0){
		memcpy(newVW,VW,(LegIdx(CurrDegree,CurrDegree)+1)*sizeof(double));
	}

	if(VW!=NULL&&CurrDegree>0){
		for(i=LegIdx(CurrDegree,CurrDegree)+1;i<LegIdx(desiredDegree,desiredDegree)+1;i++){
			newVW[i] = 0.0;
		}
	}
	else{
		for(i=0;i<LegIdx(desiredDegree,desiredDegree)+1;i++){
			newVW[i] = 0.0;
		}
	}

	M_FREE(VW);

	return newVW;

}

int CSphericalH::InitSphHData(int GrvDegree)
{
	VV = SetVWSpace(VV,0,GrvDegree);
	WW = SetVWSpace(WW,0,GrvDegree);
	VW_Degree = GrvDegree;
	VW_Order  = GrvDegree;
	VMCurrNN = 0;
	VMCurrMM = 0;
	if(VV!=NULL && WW!=NULL){
		return 1;
	}
	else{
		return 0;
	}
}

int CSphericalH::ComputeVW(double *R_Fixed,int desiredDegree,int desiredOrder,double R_ref,
						   double *VV,double *WW,int CurrDegree,int CurrOrder,
						   CSphericalHNormalFactor *pFac)
{

	double* r_bf, r_sqr,rho, x0,y0,z0;
	double Fac1,Fac2;
	int n,m;

	if(desiredDegree<=CurrDegree && desiredOrder<=CurrOrder){
		return 1;
	}

	if(pFac==NULL){
		pFac = &NoFactor;
	}

	r_bf = R_Fixed;
	r_sqr =  dot(r_bf, r_bf);
	rho   =  R_ref * R_ref / r_sqr;

	x0 = R_ref * r_bf[0] / r_sqr;
	y0 = R_ref * r_bf[1] / r_sqr;
	z0 = R_ref * r_bf[2] / r_sqr;

	if(CurrDegree==0){
		// Calculate zonal terms V(n,0); set W(n,0)=0.0
		VV[LegIdx(0,0)] = R_ref / sqrt(r_sqr)/pFac->Beta(0,0);
		WW[LegIdx(0,0)] = 0.0;

		VV[LegIdx(1,0)] = z0 * VV[LegIdx(0,0)]/pFac->Beta(1,0);
		WW[LegIdx(1,0)] = 0.0;
	}

	for(n = max(2,CurrDegree+1); n <= desiredDegree; n++){
		Fac1 = pFac->PPi(n,0,n-1,0);
		Fac2 = pFac->PPi(n,0,n-2,0);
		VV[LegIdx(n,0)] = ((2*n-1) * z0 * Fac1*VV[LegIdx(n-1,0)] - (n - 1) * rho * Fac2*VV[LegIdx(n-2,0)]) /n;
		WW[LegIdx(n,0)] = 0.0;
	}

	for (m = max(1,CurrOrder+1); m <= desiredOrder; m++){
		Fac1 = pFac->PPi(m,m,m-1,m-1);
		VV[LegIdx(m,m)] = (2 * m - 1) * Fac1 * ( x0 * VV[LegIdx(m-1,m-1)] - y0 * WW[LegIdx(m-1,m-1)] );
		WW[LegIdx(m,m)] = (2 * m - 1) * Fac1 * ( x0 * WW[LegIdx(m-1,m-1)] + y0 * VV[LegIdx(m-1,m-1)] );

		if (m < desiredDegree){
			Fac1 = pFac->PPi(m+1,m,m,m);
			VV[LegIdx(m+1,m)] = (2 * m + 1) * z0 * Fac1 * VV[LegIdx(m,m)];
			WW[LegIdx(m+1,m)] = (2 * m + 1) * z0 * Fac1 * WW[LegIdx(m,m)];
		}

		for (n = (m+2); n <= desiredDegree; n++){
			Fac1 = pFac->PPi(n,m,n-1,m);
			Fac2 = pFac->PPi(n,m,n-2,m);
			VV[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*VV[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*VV[LegIdx(n-2,m)]) / (n-m);
			WW[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*WW[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*WW[LegIdx(n-2,m)]) / (n-m);
		}
	}

	if(CurrDegree<desiredDegree){
		for (m = 1; m < CurrOrder; m++){
			for (n = CurrDegree+1; n <= desiredDegree; n++){
				Fac1 = pFac->PPi(n,m,n-1,m);
				Fac2 = pFac->PPi(n,m,n-2,m);
				VV[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*VV[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*VV[LegIdx(n-2,m)]) / (n-m);
				WW[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*WW[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*WW[LegIdx(n-2,m)]) / (n-m);
			}
		}
		if(CurrOrder>=1){
			m = CurrOrder;
			if(CurrOrder == CurrDegree && CurrOrder>=1){
				Fac1 = pFac->PPi(m+1,m,m,m);
				VV[LegIdx(m+1,m)] = (2 * m + 1) * z0 * Fac1*VV[LegIdx(m,m)];
				WW[LegIdx(m+1,m)] = (2 * m + 1) * z0 * Fac1*WW[LegIdx(m,m)];
				for (n = CurrDegree+2; n <= desiredDegree; n++){
					Fac1 = pFac->PPi(n,m,n-1,m);
					Fac2 = pFac->PPi(n,m,n-2,m);
					VV[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*VV[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*VV[LegIdx(n-2,m)]) / (n-m);
					WW[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*WW[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*WW[LegIdx(n-2,m)]) / (n-m);
				}
			}
			else{
				for (n = CurrDegree+1; n <= desiredDegree; n++){
					Fac1 = pFac->PPi(n,m,n-1,m);
					Fac2 = pFac->PPi(n,m,n-2,m);
					VV[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*VV[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*VV[LegIdx(n-2,m)]) / (n-m);
					WW[LegIdx(n,m)] = ((2*n-1)*z0*Fac1*WW[LegIdx(n-1,m)] - (n+m-1)*rho*Fac2*WW[LegIdx(n-2,m)]) / (n-m);
				}
			}
		}
	}

	return 1;

}  // End of method 'SphericalHarmonicGravity::ComputeVW()'

double CSphericalH::ComputeTensorValue(CSphericalHTensorData *pData,double *VV,double *WW,int desiredDegree,int desiredOrder)
{
	int n,m,G_idx,L_idx;
	double Value;
	if(pData==NULL||VV==NULL||WW==NULL){
		return 0.0;
	}
	if(desiredDegree<0){
		desiredDegree = pData->DegreeNumber;
	}
	if(desiredOrder<0){
		desiredOrder = desiredDegree;
	}
	Value = 0.0;
	for(n=2+pData->TensorOrder;n<=desiredDegree;n++){
		G_idx = GrIdx(n,0);L_idx = LegIdx(n,0);
		for(m=0;m<=min(n,desiredOrder);m++){
			Value += pData->CeD[G_idx]*VV[L_idx]+pData->SeD[G_idx]*WW[L_idx];
			G_idx++;L_idx++;
		}
	}
	return Value*pData->GMR;
}

CGravityTensor::CGravityTensor()
{
	Status   = 0;
	pGrvData = NULL;
	GrvFac   = &NoFactor;
	GrvFacID = 0;

	TensorDataSet.RemoveAll();
	TensorOrder = -1;
	GrvDegree   = -1;
	GrvOrder    = -1;
}

CGravityTensor::~CGravityTensor()
{
	Status = 0;
	pGrvData = NULL;
	GrvFac = NULL;
	TensorDataSet.RemoveAll();
	TensorOrder = -1;
	GrvDegree   = -1;
	GrvOrder    = -1;
}

CGravityTensor::CGravityTensor(CGravityData *pData,CSphericalHNormalFactor *pFac,int FacID)
{
	SetGrData(pData,pFac,FacID);
}

void CGravityTensor::SetFactor(CSphericalHNormalFactor* pFac)
{
	GrvFac = pFac;
	CSphericalH::SetFactor(pFac);
}

int CGravityTensor::SetGrData(CGravityData *pData,CSphericalHNormalFactor *pFac,int FacID)
{
	Status = 0;
	pGrvData = pData;
	GrvFacID = FacID;
	TensorOrder = -1;
	GrvDegree   = -1;
	GrvOrder    = -1;
	if(pFac!=NULL){
		CGravityTensor::SetFactor(pFac);
	}
	TensorDataSet.RemoveAll();
	return 1;
}

int CGravityTensor::SetDegreeOrder(int TensorOrder_in,int GrvDegree_in,int GrvOrder_in)
{

	if(pGrvData!=NULL){
		GrvDegree = min(GrvDegree_in,pGrvData->m_NN);
		GrvDegree = max(GrvDegree,2);
		GrvOrder = min(min(GrvOrder_in,GrvDegree),pGrvData->m_MM);
		GrvOrder = max(GrvOrder,0);
		CSphericalHTensorData *pTensorData=NULL,*pParentTensorData;
		TensorDataSet.RemoveAll();
		if(TensorOrder_in>=0){
			int a,b,c,Index,PatrentIndex,Num=0;
			char ch;
			TensorOrder = 0;
			pTensorData = new CSphericalHTensorData(pGrvData,GrvDegree,GrvOrder,GrvFac,GrvFacID);
			if(pTensorData!=NULL){
				TensorDataSet.PutObject(pTensorData);
			}
			else{
				return 0;
			}
			TensorOrder++;
			while(TensorOrder<=TensorOrder_in){
				Num = TensorNumber(TensorOrder);
				for(Index=0;Index<Num;Index++){
					ch = GetParent(Index,TensorOrder,&a,&b,&c);
					if(ch!='n'){
						PatrentIndex = TensorIndex(TensorOrder-1,a,b);
						pParentTensorData = TensorDataSet.GetObject(PatrentIndex);
						if(pParentTensorData!=NULL){
							pTensorData = CSphericalHTensorData::CreateCoefData(pParentTensorData,ch);
							if(pTensorData!=NULL){
								TensorDataSet.PutObject(pTensorData);
							}
							else{
								return TensorOrder-1;
							}
						}
					}
					else{
						return TensorOrder-1;
					}
				}
				TensorOrder++;
			}

			CSphericalH::SetDegreeOrder(GrvDegree_in+TensorOrder_in,GrvOrder_in+TensorOrder_in);

			return TensorOrder_in;
		}
		else{
			TensorOrder = -1;
		}
	}

	return -1;

}

int CGravityTensor::TensorNumber(int N)
{
	if(N<0){
		return 0;
	}
	else{
		return (N+2)*(N+1)/2;
	}
}

int CGravityTensor::TensorIndex(int N,int a,int b)
{
	int Index;
	if(a<0||b<0){
		return -1;
	}
	else if(a+b>N){
		return -1;
	}
	else if(a>N){
		return -1;
	}

	Index = (N+1)*(N+2)*(N+3)/6;// 前 N 阶张量数据空间
	Index = Index - (N-a+1)*(N-a+2)/2;
	Index = Index + b;

	return Index;

}

int CGravityTensor::DeFacIndex(int Index,int N,int *a,int *b,int *c)
{
	int i;

	for(i=0;i<=N;i++){
		if(0<=Index && Index<=N-i){
			*a = i;
			*b = Index;
			*c = N - *a -*b;
			return 1;
		}
		Index = Index - (N-i+1);
	}
	return 0;

}

char CGravityTensor::GetParent(int Index,int N,int *a,int *b,int *c)
{
	DeFacIndex(Index,N,a,b,c);
	if(*c>0){
		*c=*c-1;
		return 'z';
	}
	else if(*b>0){
		*b=*b-1;
		return 'y';
	}
	else if(*a>0){
		*a=*a-1;
		return 'x';
	}
	return 'n';
}

int CGravityTensor::SetPosition(double *Pos_in)
{
	if(pGrvData!=NULL && Pos_in!=NULL ){
		return CSphericalH::SetPosition(Pos_in);
	}
	else{
		return 0;
	}
}

int CGravityTensor::GetTensorValue(int a,int b,int c, double *Value)
{
	int N,Index;
	CSphericalHTensorData *pTensorData=NULL;

	a=abs(a);b=abs(b);c=abs(c);N=a+b+c;

	if(N>TensorOrder){
		return -1;
	}
	Index = TensorIndex(N,a,b);
	pTensorData = TensorDataSet.GetObject(Index);
	if(pTensorData!=NULL){
		*Value = ComputeTensorValue(pTensorData);
		return Index;
	}
	else{
		return -1;
	}


}
