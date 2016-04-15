#pragma once
#include "Model.h"
class CModelAROWD :
	public CModel
{
public:
	CModelAROWD();
	~CModelAROWD();
	void Learning(CResult *result, CData *data, CParameter *par);
};

