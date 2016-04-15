#pragma once
#include "Model.h"
class CModelAROWC :
	public CModel
{
public:
	CModelAROWC();
	~CModelAROWC();
	void Learning(CResult *result, CData *data, CParameter *par);
};

