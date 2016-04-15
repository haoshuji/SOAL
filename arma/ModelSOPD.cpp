#include "ModelSOPD.h"

CModelSOPD::CModelSOPD(){    
}

CModelSOPD::~CModelSOPD(){}

void CModelSOPD::Learning(CResult *result, CData *data, CParameter *par)
{
	int t_p=0, t_n=0, f_p=0, f_n = 0, que_num = 0;
    double hat_y_t = 0.0, p_t = 0.0, y_t = 0.0;
    	
    Matrix x_t = arma::zeros<Matrix>(data->d,1);
	Matrix Sigma_tmp = arma::zeros<Matrix>(data->d, data->d);
		
	par->Sigma.ones();
	
	par->w.zeros();

    clock_t begin, end;
    int index_tick = 0;
	time_t t;
	srand((unsigned)time(&t));
    begin = clock();	
    for(int t=0; t<data->n; t++)
    {
        int index = par->permutation[t];
		x_t = data->Get_x(index);
		//x_t =data->X.col(index);
        y_t = data->y[index]; 
        if(norm(x_t,2) == 0.0)
            continue;     
		//$(A+uv^T)^{-1} = A^{-1} - \frac{ A^{-1} u v^T A^{-1} }{ 1 + v^T A^{-1} u }$
		double xAx = arma::accu(x_t%par->Sigma%x_t);
		if (xAx == 0.0)
			continue;
		Sigma_tmp.zeros();	
		Sigma_tmp = 1 / (par->Sigma + x_t%x_t);
		p_t = arma::accu(x_t.t()*(Sigma_tmp%par->w));       

		hat_y_t = p_t>=0.0?1.0:-1.0;

		if (y_t == hat_y_t)	{
			if (y_t == 1.0)	t_p += 1;
			else t_n += 1;			
		}
		else{
			if (hat_y_t == 1.0)	f_p += 1;			
			else f_n += 1;
		}
		if (!this->alg_name.compare("SOPD")){
			que_num += 1;
			if (y_t != hat_y_t)
			{
				par->w = par->w + y_t*x_t;
				par->Sigma = Sigma_tmp;
			}
		}else if (!this->alg_name.compare("ASOPD")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			double b_2 = par->b / (par->b + fabs(p_t));
			if (rand_num <= b_2)
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					par->w = par->w + y_t*x_t;
					par->Sigma = Sigma_tmp;
				}
			}
		}else if (!this->alg_name.compare("RSOPD")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b)
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					par->w = par->w + y_t*x_t;
					par->Sigma = Sigma_tmp;
				}
			}
		}else if (!this->alg_name.compare("AASOPD")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + 2 * fabs(p_t) + p_t*p_t*(xAx + 1)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					par->w = par->w + y_t*x_t;
					par->Sigma = Sigma_tmp;
				}
			}
		}
		else if (!this->alg_name.compare("AASOPD2")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + fabs(p_t) / sqrt(xAx)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					par->w = par->w + y_t*x_t;
					par->Sigma = Sigma_tmp;
				}
			}
		}
		else{
			cout << "Unkonw algorithm name:\t" << this->alg_name << endl;
		}
        
        end = clock();          
        if ( ((t+1) % (par->time_ticks) == 0) && index_tick < (par->num_ticks-1))
        {                       
            result->f_n[index_tick] = f_n;
            result->f_p[index_tick] = f_p;
            result->t_n[index_tick] = t_n;
            result->t_p[index_tick] = t_p;           
            result->que_num[index_tick] = que_num;
            result->time_[index_tick]=(double)(end-begin)/CLOCKS_PER_SEC;
            index_tick++;
        }
    }
    end = clock();
    result->f_n[index_tick] = f_n;
    result->f_p[index_tick] = f_p;
    result->t_n[index_tick] = t_n;
    result->t_p[index_tick] = t_p;           
    result->que_num[index_tick] = que_num;
    result->time_[index_tick]=(double)(end-begin)/CLOCKS_PER_SEC;

	result->que = (double)que_num / data->n;
	result->acc = 1 - (double)(f_n + f_p) / data->n;
	result->time = (double)(end - begin) / CLOCKS_PER_SEC;
	result->F1 = (2.0*t_p) / (2.0*t_p + f_n + f_p);
}
