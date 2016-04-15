/*
 * =====================================================================================
 *
 *       Filename:  Model.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/03/2014 09:38:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuji Hao (), hao.shuji@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#include "Model.h"

string CModel::Get_alg_name()
{
    return alg_name;
}
void CModel::SetAlgName(string alg_name){
	this->alg_name = alg_name;
}