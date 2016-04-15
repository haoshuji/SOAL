#include "ModelSOPD.h"

CModelSOPD::CModelSOPD(){    
}

CModelSOPD::~CModelSOPD(){}

void CModelSOPD::Learning(CResult *result, CData *data, CParameter *par)
{
	int t_p=0, t_n=0, f_p=0, f_n = 0, que_num = 0;
    double hat_y_t = 0.0, p_t = 0.0, y_t = 0.0;
    int d = data->d;    
	struct CFeature_node *x_t, *x_t_tmp; 
	par->Reset(data->d,data->n);
	double *Sigma=new double[d], *Sigma_tmp = new double[d];
	for (int i = 0; i < d; i++)
	{
		Sigma[i] = 1.0;	
		Sigma_tmp[i] = 1.0;
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

		x_t_tmp = x_t;
		while (x_t_tmp->index != -1)
		{
			int i = x_t_tmp->index - 1;
			double value = x_t_tmp->value;
			Sigma_tmp[i] = Sigma[i] + value * value;
			x_t_tmp++;
		}

		p_t = 0.0;
		double v_t = 0.0;
		x_t_tmp = x_t;
		while (x_t_tmp->index != -1)
		{
			int i = x_t_tmp->index - 1;
			p_t += par->w[i]*(1.0/Sigma_tmp[i])*x_t_tmp->value;
			v_t += x_t_tmp->value*(1.0/Sigma[i])*x_t_tmp->value;
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

		if (!this->alg_name.compare("SOPD")){
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
					Sigma[i] = Sigma_tmp[i];
				}
			}
		}else if (!this->alg_name.compare("ASOPD")){
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
						Sigma[x_t_tmp->index-1] += x_t_tmp->value*x_t_tmp->value;
						x_t_tmp++;
					}
				}
			}
		}else if (!this->alg_name.compare("RSOPD")){
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
						Sigma[x_t_tmp->index-1] += x_t_tmp->value*x_t_tmp->value;
						x_t_tmp++;
					}
				}
			}
		}else if (!this->alg_name.compare("AASOPD")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + 2 * fabs(p_t) + p_t*p_t*(v_t + 1)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						Sigma[x_t_tmp->index-1] += x_t_tmp->value*x_t_tmp->value;
						x_t_tmp++;
					}
				}
			}
		}
		else if (!this->alg_name.compare("AASOPD2")){
			double rand_num = (double)rand() / (double)RAND_MAX;
			if (rand_num <= par->b / (par->b + fabs(p_t) / sqrt(v_t)))
			{
				que_num += 1;
				if (y_t != hat_y_t)
				{
					struct CFeature_node *x_t_tmp = x_t;
					while(x_t_tmp->index != -1)
					{
						par->w[x_t_tmp->index-1] += y_t*x_t_tmp->value;
						Sigma[x_t_tmp->index-1] += x_t_tmp->value*x_t_tmp->value;
						x_t_tmp++;
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

	delete[] Sigma;
	delete[] Sigma_tmp;
}
