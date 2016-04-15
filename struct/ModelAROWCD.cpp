#include "ModelAROWCD.h"


CModelAROWCD::CModelAROWCD()
{
}


CModelAROWCD::~CModelAROWCD()
{
}



void CModelAROWCD::Learning(CResult *result, CData *data, CParameter *par)
{
	int t_p = 0, t_n = 0, f_p = 0, f_n = 0, que_num = 0;
	double hat_y_t = 0.0, p_t = 0.0, y_t = 0.0;
	par->Reset(data->d,data->n);
	int d = data->d;
	struct CFeature_node *x_t, *x_t_tmp; 
	double *Sigma = new double[d], *SigmaX = new double[d];
	for (int i = 0; i < d; i++)
	{
		Sigma[i] = 1.0;
	}

	clock_t begin, end;
	int index_tick = 0;
	
	srand((unsigned)time(NULL));

	begin = clock();
	for (int t = 0; t<data->n; t++)
	{
		int index = par->permutation[t];
        x_t = data->x[index];		
        y_t = data->y[index];
		double x_t_norm = CVector::Norm2(x_t);
        if(x_t_norm== 0.0)
			continue;
		p_t = CVector::WDotXt(par->w,x_t);

		if (p_t >= 0) { hat_y_t = 1.0; }
		else { hat_y_t = -1.0; }

		if (y_t == hat_y_t){ 
			if (y_t == 1.0){ t_p += 1;	}
			else{ t_n += 1; }			
		}
		else{
			if (hat_y_t == 1.0)	{ f_p += 1;	}
			else { f_n += 1; }
		}
		
		if (!this->alg_name.compare("AROWCD")){
			que_num += 1;
			double ell_t = 1 - y_t*p_t;
			if (ell_t > 0.0)
			{	
				x_t_tmp = x_t;			
				while (x_t_tmp->index != -1)
				{
					SigmaX[x_t_tmp->index-1] = Sigma[x_t_tmp->index-1] * x_t_tmp->value;
					x_t_tmp++;
				}
				
				double xSx = 0.0;
				x_t_tmp = x_t;
				while (x_t_tmp->index != -1)
				{
					xSx += SigmaX[x_t_tmp->index-1] * x_t_tmp->value;
					x_t_tmp++;
				}				
				double beta_t = 1 / (xSx + par->AROW_r);
				double alpha_t = max(0.0, ell_t)*beta_t;
				
				x_t_tmp = x_t;					
				while (x_t_tmp->index != -1)
				{
					par->w[x_t_tmp->index-1] += alpha_t*y_t*Sigma[x_t_tmp->index-1] * x_t_tmp->value;
					x_t_tmp++;
				}		

				x_t_tmp = x_t;
				while (x_t_tmp->index!=-1)
				{
					int i = x_t_tmp->index-1;
					Sigma[i] = Sigma[i] - beta_t* SigmaX[i] * SigmaX[i];
					x_t_tmp++;
				}			
			}
		}	
		else if(!this-alg_name.compare("AAROWCD")){
		
		}
		else{
			cout << "Unkonw algorithm name:\t" << this->alg_name << endl;
		}

		end = clock();
		if (((t + 1) % (par->time_ticks) == 0) && index_tick < (par->num_ticks - 1))
		{
			result->f_n[index_tick] = f_n;
			result->f_p[index_tick] = f_p;
			result->t_n[index_tick] = t_n;
			result->t_p[index_tick] = t_p;
			result->que_num[index_tick] = que_num;
			result->time_[index_tick] = (double)(end - begin) / CLOCKS_PER_SEC;
			index_tick++;
		}
	}
	end = clock();
	result->f_n[index_tick] = f_n;
	result->f_p[index_tick] = f_p;
	result->t_n[index_tick] = t_n;
	result->t_p[index_tick] = t_p;
	result->que_num[index_tick] = que_num;
	result->time_[index_tick] = (double)(end - begin) / CLOCKS_PER_SEC;

	result->que = (double)que_num / data->n;
	result->acc = 1 - (double)(f_n + f_p) / data->n;
	result->time = (double)(end - begin) / CLOCKS_PER_SEC;
	result->F1 = (2.0*t_p) / (2.0*t_p + f_n + f_p);

	delete[] Sigma;
	delete[] SigmaX;
}
