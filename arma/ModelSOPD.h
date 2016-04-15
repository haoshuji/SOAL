#ifndef CModelSOPD_H
#define CModelSOPD_H

#include "Model.h"

class CModelSOPD:public CModel{
public:	
	void Learning(CResult *result, CData *data, CParameter *par);
	
	CModelSOPD();
	
	~CModelSOPD();
};

#endif
