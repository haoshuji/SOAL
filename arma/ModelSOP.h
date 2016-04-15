#ifndef CModelSOP_H
#define CModelSOP_H

#include "Model.h"

class CModelSOP:public CModel{
public:	
	void Learning(CResult *result, CData *data, CParameter *par);
	
	CModelSOP();
	
	~CModelSOP();
};

#endif
