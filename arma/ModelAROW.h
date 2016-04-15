#pragma once
#include "Model.h"
class CModelAROW :
	public CModel
{
public:
	CModelAROW();
	~CModelAROW();
	void Learning(CResult *result, CData *data, CParameter *par);
};

