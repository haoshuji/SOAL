#include "Parameter.h"
#include <cmath>
#include <iostream>

CParameter::CParameter()
	: Mul2Bin(0)
	, find_AROW_r(1)
	, find_AROW_eta(1)
	, find_PAI_C(1)
	, find_PAII_C(1)
	, find_AROWC_r(1)
	, que_increase_speed(2.0)
	, F1_or_acc(false)
{
    num_alg = 0;
	num_fold = 5;
	index_binary_class = 1;
	num_que = 10;
	num_ticks = 10;	
	b_start = -10.0;
	b = 1.0;	
	index_binary_class = 1;
	Norm2One = 1;
	Full_Matrix = 1;

	AROW_r = 0.1;	AROW_eta = 1;
	AROWC_r = 1;	PAI_C = 1;
	PAII_C = 1;

}

void CParameter::ImportParameters(std::string file_fullpath){
	
	std::string name_par;
	double value_par;

	std::ifstream infile(file_fullpath.c_str());

	while (infile >> name_par >> value_par)
	{
		if (!name_par.compare("num_alg"))
			this->num_alg = (int)value_par;		
		else if (!name_par.compare("num_fold"))
			this->num_fold = (int)value_par;
		else if (!name_par.compare("index_binary_class"))
			this->index_binary_class = (int)value_par;
		else if (!name_par.compare("num_que"))
			this->num_que = (int)value_par;
		else if (!name_par.compare("que_increase_speed"))
			this->que_increase_speed = (double)value_par;
		else if (!name_par.compare("b_start"))
			this->b_start = (double)value_par;
		else if (!name_par.compare("que_increase_speed2"))
			this->que_increase_speed2 = (double)value_par;
		else if (!name_par.compare("b_start2"))
			this->b_start2 = (double)value_par;
		else if (!name_par.compare("num_ticks"))
			this->num_ticks = (int)value_par;				
		else if (!name_par.compare("PAI_C"))
			this->PAI_C = (double)value_par;
		else if (!name_par.compare("PAII_C"))
			this->PAII_C = (double)value_par;
		else if (!name_par.compare("find_PAI_C"))
			this->find_PAI_C = (int)value_par;
		else if (!name_par.compare("find_PAII_C"))
			this->find_PAII_C = (int)value_par;
		else if (!name_par.compare("AROW_r"))
			this->AROW_r = (double)value_par;
		else if (!name_par.compare("AROWC_r"))
			this->AROWC_r = (double)value_par;
		else if (!name_par.compare("AROW_eta"))
			this->AROW_eta = (double)value_par;
		else if (!name_par.compare("find_AROW_r"))
			this->find_AROW_r = (int)value_par;
		else if (!name_par.compare("find_AROW_eta"))
			this->find_AROW_eta = (int)value_par;
		else if (!name_par.compare("find_AROWC_r"))
			this->find_AROWC_r = (int)value_par;
		else if (!name_par.compare("Norm2One"))
			this->Norm2One = (int)value_par;		
		else if (!name_par.compare("Mul2Bin"))
			this->Mul2Bin = (int)value_par;
		else if (!name_par.compare("Full_Matrix"))
			this->Full_Matrix= (int)value_par;
		else if (!name_par.compare("F1_or_acc"))
			this->F1_or_acc= (bool)value_par;
		else
		{
			cout << "Unknonw par in setting file" << file_fullpath << endl;
			exit(EXIT_FAILURE);
		}
	}
	infile.close();
	std::cout << "Parsed setting file" << std::endl;
}

void CParameter::Initialize(int d, int n)
{
	time_ticks = (int)floor((double)n/num_ticks - 0.4);
//	if (diag == 0)
		//Sigma = eye<mat>(d, d);	
        //
	//Sigma_vec = Reset<vec>(d);
	if (this->Full_Matrix == 1)
		Sigma = arma::eye<Matrix>(d, d);
	else		
		Sigma = arma::ones<Matrix>(d, 1);
	w = arma::zeros<Vector>(d);
	//w = Reset<vec>(d);
	permutation = new int[n];
	for (int i = 0; i < n; i++)
	{
		permutation[i] = i;
	}
}
void CParameter::Reset(int d, int n){
	if (this->Full_Matrix == 1)
		Sigma.ones();
	else
		Sigma.eye();

	this->w.zeros();	
}
void CParameter::Generate_permutation(int n)
{
    if(permutation != NULL)
    {
	    for(int t=0;t<n;t++)
	        permutation[t] = t;
	    
		for (int i = 0; i<n; i++)
		{
			int j = rand() % (n - i) + i;
			int t = permutation[j];
			permutation[j] = permutation[i];
			permutation[i] = t;
		}
	}
	else
		exit(1);
}

CParameter::~CParameter()
{
	if (permutation != NULL)
	{
		delete []permutation;
	}

}
