#ifndef __ITERATION_H__
#define __ITERATION_H__

class CIteration
{
public:
    int Idx;
	void * Page;
    void Initialize(){Idx=-1;Page=NULL;}
	bool operator==(const CIteration & itr)
	{
		return (this==&itr) || ((Idx==itr.Idx)&&(Page==itr.Page));
	}
	bool operator!=(const CIteration & itr)
	{
		return ((Idx!=itr.Idx)||(Page!=itr.Page));
	}
};

#endif