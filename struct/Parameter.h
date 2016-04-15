#ifndef PARAMETER_H
#define PARAMETER_H
#include "Global.h"

class CParameter
{
	public:
	//initialize by constructor and set by argv
	int num_alg, num_fold, index_binary_class;
	int num_que, num_ticks, time_ticks;
	
	int *permutation;
	//0 \leq bbq_k \leq 1, 
	double  b_start, b_start2,b;
	
	double *w;
	//vec Sigma_vec;

	double PAI_C, PAII_C;
	double AROW_r, AROW_eta; //parameters proposed in our algorithms
	double AROWC_r; //parameter r in crammer AROW algorithm
	CParameter();

	void ImportParameters(std::string file_fullpath);

	void Generate_permutation(int n);

	void Initialize(int d, int n);
	
	void Reset(int d, int n);

	~CParameter();
	int Mul2Bin;
	int Norm2One;
	int find_AROW_r;	
	int find_AROW_eta;
	int find_PAI_C;
	int find_PAII_C;
	int find_AROWC_r;

	double que_increase_speed,que_increase_speed2;
	int Full_Matrix;
	bool F1_or_acc;
};

#endif
