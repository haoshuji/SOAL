#ifndef CModelPE_H
#define CModelPE_H

#include "Model.h"

class CModelPE:public CModel{
public:	
	void Learning(CResult *result, CData *data, CParameter *par);
	
	CModelPE();
	
	~CModelPE();
};

#endif
