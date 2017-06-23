#ifndef __OBJECTVECTOR_H__
#define __OBJECTVECTOR_H__

#include <malloc.h>
#include <string.h>

#endif

#include "Iteration.h"

template <typename OBJECT>
class ObjVecData
{
public:
	int Status;
	OBJECT *pObj;
};

#ifndef dimabs
#define dimabs(x) ((x)>0?(x):-(x))
#endif

#define ObjPagePTR ObjVecPage<OBJECT,PtsPG>*

template <typename OBJECT,int PtsPG>
struct ObjVecPage
{
	int EmptyNum;
	ObjPagePTR PrevPage;
	ObjVecData<OBJECT> Pts[dimabs(PtsPG)];
	ObjPagePTR NextPage;
};

#ifndef NULL
#define NULL 0
#endif

template <typename OBJECT,int PtsPG>
class CObjVec
{
public:

	CObjVec();

	virtual ~CObjVec();

    void PutObject( OBJECT *pt,int status=0);

	OBJECT * ReplaceCurrObject(OBJECT *pt,CIteration *itr=NULL);

	void RemoveAll(int Flag=0);

	bool IsExist(OBJECT *pt,int Flag=0,CIteration *itr=NULL);

	bool MoveFirst(CIteration *itr=NULL) const;

	bool MoveForward(CIteration *itr=NULL)  const;

	bool GotoIndex(long Index,CIteration *itr=NULL)  const;

	OBJECT * CurrObject(CIteration *itr=NULL)  const;

	OBJECT * GetObject(long index,CIteration *itr=NULL)  const;

private:
	void ReConstruct(void);
	ObjPagePTR NewPage(void);
	void AppendNewPage(void);

	mutable CIteration itr0;
	long TotalPages, TotalObjects;

	ObjPagePTR FirstPage;
	ObjPagePTR LastPage;

	ObjVecData<OBJECT> * ReplaceCurrData(ObjVecData<OBJECT> *pt,CIteration *itr=NULL);

};

template <typename OBJECT,int PtsPG>
CObjVec<OBJECT,PtsPG>::CObjVec()
{
	LastPage = FirstPage = NULL;
	ReConstruct();
}

template <typename OBJECT,int PtsPG>
inline void CObjVec<OBJECT,PtsPG>::ReConstruct(void)
{
	if(FirstPage==NULL){
		LastPage = FirstPage = NewPage();
	}
	else{
		ObjPagePTR p;
		while(LastPage->PrevPage!=NULL){
			p = LastPage;
			LastPage = LastPage->PrevPage;
			LastPage->NextPage = NULL;
			free(p);
		}
		for(int i=0;i<dimabs(PtsPG);i++){
			FirstPage->Pts[i].pObj = NULL;
		}
		FirstPage->PrevPage = NULL;
		FirstPage->NextPage = NULL;
		FirstPage->EmptyNum = dimabs(PtsPG);
	}
	TotalPages    = 1;
	TotalObjects  = 0;
	itr0.Initialize();
}

template <typename OBJECT,int PtsPG>
CObjVec<OBJECT,PtsPG>::~CObjVec()
{
	RemoveAll();
	ObjPagePTR p;
	while(FirstPage!=NULL){
		p = FirstPage;
		FirstPage = FirstPage->NextPage;
		free(p);
	}
}

template <typename OBJECT,int PtsPG>
ObjPagePTR CObjVec<OBJECT,PtsPG>::NewPage(void)
{
	ObjPagePTR p;
	if((p=(ObjPagePTR)malloc(sizeof(ObjVecPage<OBJECT,PtsPG>)))==NULL){
		throw;
	}
	p->EmptyNum = dimabs(PtsPG);
	for(int i=0;i<dimabs(PtsPG);i++){
		p->Pts[i].Status = 0;
		p->Pts[i].pObj = NULL;
	}
	p->PrevPage = NULL;
	p->NextPage = NULL;
	return p;
}

template <typename OBJECT,int PtsPG>
void CObjVec<OBJECT,PtsPG>::AppendNewPage()
{
	LastPage->NextPage = NewPage();
	LastPage->NextPage->PrevPage = LastPage;
	TotalPages++;
	LastPage = LastPage->NextPage;
}

template <typename OBJECT,int PtsPG>
bool CObjVec<OBJECT,PtsPG>::IsExist(OBJECT *pt_i,int Flag,CIteration *itr)
{
	CIteration itr1;
	bool Status = false;
	if(Flag==0){
		itr1 = itr0;
		itr  = &itr1;
	}
	else{
		if(itr==NULL){
			itr = &itr0;
		}
	}
	if(MoveFirst(itr)){
		OBJECT*pt;
		while(true){
			pt = CurrObject(itr);
			if(pt==pt_i){
				Status = true;
			}
			if(!MoveForward(itr)){
				break;
			}
		}
	}
	return Status;
}

template <typename OBJECT,int PtsPG>
inline void CObjVec<OBJECT,PtsPG>::PutObject(OBJECT *pt,int status)
{
	CIteration *itr = &itr0;
	if(pt==NULL){
		return;
	}
	else if(pt!=NULL && IsExist(pt)){
		return;
	}
	else if(pt!=NULL && IsExist(NULL,1)){
		ObjVecData<OBJECT> Pt;
		Pt.pObj = pt;
		Pt.Status = status;
		ReplaceCurrData(&Pt,itr);
	}
	else{
		itr->Idx = TotalObjects%dimabs(PtsPG);
		if(TotalObjects!=0 && itr->Idx%dimabs(PtsPG) == 0){
			AppendNewPage();
		}
		itr->Page  = LastPage;
		((ObjPagePTR)(itr->Page))->Pts[itr->Idx].Status = status;
		((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj = pt;
		((ObjPagePTR)(itr->Page))->EmptyNum--;
		TotalObjects++;
	}
}

template <typename OBJECT,int PtsPG>
ObjVecData<OBJECT> * CObjVec<OBJECT,PtsPG>::ReplaceCurrData(ObjVecData<OBJECT> *pt,CIteration *itr)
{
	if(itr==NULL){
		itr = &itr0;
	}
	if(itr->Idx<0){
		return NULL;
	}
	else if(pt==NULL){
		ObjVecData<OBJECT> Pt1;
		Pt1.Status = ((ObjPagePTR)(itr->Page))->Pts[itr->Idx].Status;
		Pt1.pObj = ((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj;
		if(Pt1.pObj!=NULL){
			((ObjPagePTR)(itr->Page))->Pts[itr->Idx].Status = pt->Status;
			((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj = pt->pObj;
			((ObjPagePTR)(itr->Page))->EmptyNum++;
		}
		pt->pObj = Pt1.pObj;
		pt->Status = Pt1.Status; 
		return pt;
	}
	else if(!IsExist(pt->pObj)){
		ObjVecData<OBJECT> Pt1;
		Pt1.pObj = ((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj;
		Pt1.Status = ((ObjPagePTR)(itr->Page))->Pts[itr->Idx].Status;
		((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj = pt->pObj;
		((ObjPagePTR)(itr->Page))->Pts[itr->Idx].Status = pt->Status;
		if(Pt1.pObj==NULL){
			((ObjPagePTR)(itr->Page))->EmptyNum--;
		}
		pt->pObj = Pt1.pObj;
		pt->Status = Pt1.Status; 
		return pt;
	}
	else{
		return NULL;
	}
}

template <typename OBJECT,int PtsPG>
inline bool CObjVec<OBJECT,PtsPG>::MoveFirst(CIteration *itr) const
{
	if(itr==NULL){
		itr = &itr0;
	}
	if(TotalObjects > 0){
		itr->Page = FirstPage;
		itr->Idx  = 0;
		return true;
	}
	else{
		itr->Page = NULL;
		itr->Idx  = -1;
		return false;
	}
}

template <typename OBJECT,int PtsPG>
inline bool CObjVec<OBJECT,PtsPG>::MoveForward(CIteration *itr) const
{
	if(itr==NULL){
		itr = &itr0;
	}
	if(itr->Page==NULL){
		return false;
	}
	else if(((ObjPagePTR)(itr->Page))->NextPage!=NULL){
		if(itr->Idx<dimabs(PtsPG)-1){
			itr->Idx++;
			return true;
		}
		else{
			itr->Page  = ((ObjPagePTR)(itr->Page))->NextPage;
			itr->Idx = 0;
			return true;
		}
	}
	else{
		if((int)itr->Idx<(TotalObjects%dimabs(PtsPG))-1){
			itr->Idx++;
			return true;
		}
		else{
			return false;
		}
	}
}

template <typename OBJECT,int PtsPG>
inline bool CObjVec<OBJECT,PtsPG>::GotoIndex(long index,CIteration *itr) const
{
	if(itr==NULL){
		itr = &itr0;
	}
	index++;
	if(index > TotalObjects || index <=0 ){
		return false;
	}
	itr->Page = FirstPage;
	while(index > dimabs(PtsPG)){
		itr->Page = ((ObjPagePTR)(itr->Page))->NextPage;
		index -= dimabs(PtsPG);
	}
	itr->Idx = index-1;
	return true;
}

template <typename OBJECT,int PtsPG>
inline OBJECT * CObjVec<OBJECT,PtsPG>::CurrObject(CIteration *itr) const
{
	if(itr==NULL){
		itr = &itr0;
	}
	if(itr->Page==NULL){
		return NULL;
	}
	else{
		return ((ObjPagePTR)(itr->Page))->Pts[itr->Idx].pObj;
	}
}

template <typename OBJECT,int PtsPG>
inline OBJECT * CObjVec<OBJECT,PtsPG>::GetObject(long index,CIteration *itr) const
{
	if(itr==NULL){
		itr = &itr0;
	}
	if(GotoIndex(index)){
		return CurrObject(itr);
	}
	else{
		return NULL;
	}
}

template <typename OBJECT,int PtsPG>
void CObjVec<OBJECT,PtsPG>::RemoveAll(int Flag)
{
	CIteration *itr = &itr0;
	if(Flag==1){

	}
	MoveFirst(itr);
	ObjVecData<OBJECT> Pt;
	while(true){
		Pt.pObj = NULL;
		Pt.Status = 0;
		ReplaceCurrData(&Pt,itr);
		if(PtsPG>0 && Pt.pObj!=NULL && Pt.Status==0){
			delete Pt.pObj;
		}
		if(!MoveForward(itr)){
			break;
		}
	}
	ReConstruct();
	return;
}