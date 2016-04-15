 #ifndef Model_H
 #define Model_H
/*
 * =====================================================================================
 *        Class:  CModel
 *  Description:  abstract class for specific online active learning algorithm
 * =====================================================================================
 */

#include "Result.h"
#include "Data.h"
#include "Parameter.h"
#include "armadillo"
#include <algorithm>
#include <iostream>
#include "Global.h"
#include "Vector.h"
using namespace arma;
using namespace std;

class CModel
{
    public:
        string alg_name;

		CModel(){}                            /* constructor */
        std::string Get_alg_name();
		void SetAlgName(string alg_name);
        virtual void Learning(CResult *result, CData *data, CParameter *par){};        
		~CModel(){}
        
}; /* -----  end of class CModel  ----- */

#endif
