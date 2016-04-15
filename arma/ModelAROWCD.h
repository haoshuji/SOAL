#pragma once
#include "Model.h"
class CModelAROWCD :
	public CModel
{
public:
	CModelAROWCD();
	~CModelAROWCD();
	void Learning(CResult *result, CData *data, CParameter *par);
};

