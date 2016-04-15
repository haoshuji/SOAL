#include "ModelAROWC.h"


CModelAROWC::CModelAROWC()
{
}


CModelAROWC::~CModelAROWC()
{
}



void CModelAROWC::Learning(CResult *result, CData *data, CParameter *par)
{
	int t_p = 0, t_n = 0, f_p = 0, f_n = 0, que_num = 0;
	double hat_y_t = 0.0, p_t = 0.0, y_t = 0.0;
	par->Reset(data->d,data->n);
	int d = data->d;
	double **Sigma=new double*[d];
	for (int i = 0; i < d; i++)
	{
		Sigma[i] = new double[d];	
		for (int j = 0; j < d; j++)
		{
			Sigma[i][j] = 0.0;			
		}
		Sigma[i][i] = 1.0;
	}
	double *SigmaX = new double[d];
	struct CFeature_node *x_t, *x_t_tmp; 

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
		
		if (!this->alg_name.compare("AROWC")){
			que_num += 1;
			double ell_t = 1 - y_t*p_t;
			if (ell_t > 0.0)
			{				
				for (int i = 0; i < d; i++)
				{
					SigmaX[i] = 0.0;
					x_t_tmp = x_t;			
					while (x_t_tmp->index != -1)
					{
						SigmaX[i] += Sigma[i][x_t_tmp->index-1] * x_t_tmp->value;
						x_t_tmp++;
					}
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
				for (int i = 0; i < d; i++)
				{
					x_t_tmp = x_t;
					double tmp = 0.0;
					while (x_t_tmp->index != -1)
					{
						tmp += Sigma[i][x_t_tmp->index-1] * x_t_tmp->value;
						x_t_tmp++;
					}
					par->w[i] = par->w[i] + alpha_t*y_t*tmp;
				}			
				//par->w = par->w + alpha_t*y_t*par->Sigma*x_t;
				for (int i = 0; i < d; i++)	{
				for (int j = 0; j < d; j++)	{
					Sigma[i][j] = Sigma[i][j] - beta_t* SigmaX[i] * SigmaX[j];
					}
				}								
			}
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

	for (int i = 0; i < d; i++)
	{
		delete[] Sigma[i];
	}
	delete[] Sigma;
	delete[] SigmaX;
}
