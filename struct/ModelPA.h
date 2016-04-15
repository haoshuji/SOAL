#ifndef CModelPA_H
#define CModelPA_H

#include "Model.h"

class CModelPA:public CModel{
public:	
	void Learning(CResult *result, CData *data, CParameter *par);
	
	CModelPA();
	
	~CModelPA();
};

#endif
