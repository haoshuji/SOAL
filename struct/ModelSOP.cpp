#include "ModelSOP.h"

CModelSOP::CModelSOP(){    
}

CModelSOP::~CModelSOP(){}

void CModelSOP::Learning(CResult *result, CData *data, CParameter *par)
{
	int t_p=0, t_n=0, f_p=0, f_n = 0, que_num = 0;
    double hat_y_t = 0.0, p_t = 0.0, y_t = 0.0;
	int d = data->d;
	par->Reset(data->d,data->n);
	struct CFeature_node *x_t, *x_t_tmp; 

	double **Sigma=NULL, **Sigma_tmp = NULL;
	double *SigmaX = new double[d];
	Sigma = new double*[d]; Sigma_tmp = new double*[d];
	for (int i = 0; i < d; i++)
	{
		Sigma[i] = new double[d];		
		Sigma_tmp[i] = new double[d];
		for (int j = 0; j < d; j++)
		{
			Sigma[i][j] = 0.0;		
			Sigma_tmp[i][j] = 0.0;
		}
		Sigma[i][i] = 1.0;
	}

    clock_t begin, end;
    int index_tick = 0;
	time_t t;
	srand((unsigned)time(&t));
    begin = clock();	
    for(int t=0; t<data->n; t++)
    {
        int index = par->permutation[t];
        x_t = data->x[index];		
        y_t = data->y[index];
		double x_t_norm = CVector::Norm2(x_t);
        if(x_t_norm== 0.0)
			continue;
		
		//compute A^{-1} * x_t
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
		
		//compute x_t^T * A^{-1} * x_t
		double xSx = 0.0;
		x_t_tmp = x_t;
		while (x_t_tmp->index != -1)
		{
			xSx += SigmaX[x_t_tmp->index-1] * x_t_tmp->value;
			x_t_tmp++;
		}

		//compute 1/(1+xSx)
		double factor = 1.0/(1.0+xSx);
		for (int i = 0; i < d; i++)	{
			for (int j = 0; j < d; j++)	{
				Sigma_tmp[i][j] = Sigma[i][j] - factor * SigmaX[i] * SigmaX[j];
			}
		}
		x_t_tmp = x_t;		p_t = 0.0;
		while (x_t_tmp->index != -1)
		{
			int j = x_t_tmp->index - 1;
			double tmp=0.0; //w^T * Sigma[:][j]
			for (int i = 0; i < d; i++)
			{
				tmp += par->w[i]*Sigma_tmp[i][j];
			}
			p_t += tmp * x_t_tmp->value;
			x_t_tmp++;
		}
        if(p_t >= 0) {hat_y_t = 1.0;}
        else {hat_y_t = -1.0;}

		if (y_t == hat_y_t){
			if (y_t == 1.0){ t_p += 1; }
			else{ t_n += 1; }
		}
		else{
			if (hat_y_t == 1.0)	{ f_p += 1;	}
			else{ f_n += 1; }
		}

		if (!this->alg_name.compare("SOP")){
			que_num += 1;
			if (y_t != hat_y_t)
			{
				struct CFeature_node *x_t_tmp = x_t;
				while(x_t_tmp->index != -1)
				{
					par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
					x_t_tmp++;
				}
				for (int i = 0; i < d; i++)
				{
					for (int j = 0; j < d; j++)
					{
						Sigma[i][j] = Sigma_tmp[i][j];
					}
				}				
			}
		}else if (!this->alg_name.compare("ASOP")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			double b_2 = par->b / (par->b + fabs(p_t));
			if (rand_num <= b_2)
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						x_t_tmp++;
					}
					for (int i = 0; i < d; i++)
					{
						for (int j = 0; j < d; j++)
						{
							Sigma[i][j] = Sigma_tmp[i][j];
						}
					}			
				}
			}
		}else if (!this->alg_name.compare("RSOP")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b)
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						x_t_tmp++;
					}
					for (int i = 0; i < d; i++)
					{
						for (int j = 0; j < d; j++)
						{
							Sigma[i][j] = Sigma_tmp[i][j];
						}
					}			
				}
			}
		}else if (!this->alg_name.compare("AASOP")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + 2 * fabs(p_t) + p_t*p_t*(xSx + 1)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						x_t_tmp++;
					}
					for (int i = 0; i < d; i++)
					{
						for (int j = 0; j < d; j++)
						{
							Sigma[i][j] = Sigma_tmp[i][j];
						}
					}			
				}
			}
		}
		else if (!this->alg_name.compare("AASOP2")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + fabs(p_t) / sqrt(xSx)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						x_t_tmp++;
					}
					for (int i = 0; i < d; i++)
					{
						for (int j = 0; j < d; j++)
						{
							Sigma[i][j] = Sigma_tmp[i][j];
						}
					}			
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

	for (int i = 0; i < d; i++)
	{
		delete[] Sigma[i];
		delete[] Sigma_tmp[i];
	}
	delete[] Sigma;
	delete[] Sigma_tmp;
	delete[] SigmaX;
}
