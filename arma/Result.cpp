#include "Result.h"
#include <stdio.h>

CResult::CResult()
{
    num_ticks = 0;
    f_n = NULL;
    f_p = NULL;
    t_n = NULL;
    t_p = NULL;
    que_num = NULL;
    time_ = NULL;

	que = 0.0;
	acc = 0.0;
	time = 0.0;
	F1 = 0.0;
	this->std_acc = 0.0;
	this->std_F1 = 0.0;
	this->std_que = 0.0;
	this->std_time = 0.0;
}

CResult::CResult(int m_num_ticks)
{
    num_ticks = m_num_ticks;
    f_n = new double[num_ticks];
    f_p = new double[num_ticks];
    t_n = new double[num_ticks];
    t_p = new double[num_ticks];
    que_num = new double[num_ticks];
    time_ = new double[num_ticks];
}

void CResult::Initialize(int m_num_ticks)
{
    // if(!f_n) delete []f_n;
    // if(!f_p) delete []f_p;
    // if(!t_n) delete []t_n;
    // if(!t_p) delete []t_p;
    // if(!que_num) delete []que_num;
    // if(!time_) delete []time_;
    if(!f_n)
    {
        num_ticks = m_num_ticks;
        f_n = new double[num_ticks];
        f_p = new double[num_ticks];
        t_n = new double[num_ticks];
        t_p = new double[num_ticks];
        que_num = new double[num_ticks];
        time_ = new double[num_ticks];
    }
     for (int i = 0; i < num_ticks; ++i)
    {
        f_n[i] = 0;
        f_p[i] = 0;
        t_n[i] = 0;
        t_p[i] = 0;
        que_num[i] = 0;
        time_[i] = 0;
    }

	 que = 0.0;
	 acc = 0.0;
	 time = 0.0;
	 F1 = 0.0;
	 this->std_acc = 0.0;
	 this->std_F1 = 0.0;
	 this->std_que = 0.0;
	 this->std_time = 0.0;
}

CResult& CResult::operator=(const CResult &result)
{
    if (num_ticks == 0)
    {
        Initialize(result.num_ticks);
    }

    for (int i = 0; i < num_ticks; ++i)
    {
        f_n[i] = result.f_n[i];
        f_p[i] = result.f_p[i];
        t_n[i] = result.t_n[i];
        t_p[i] = result.t_p[i];
        que_num[i] = result.que_num[i];
        time_[i] = result.time_[i];
    }
	que = result.que;
	acc = result.acc;
	time = result.time;
	F1 = result.F1;

	this->std_acc = result.std_acc;
	this->std_F1 = result.std_F1;
	this->std_que = result.std_que;
	this->std_time = result.std_time;

	return *this;
}

CResult& CResult::operator+(const CResult &result)
{
    if (num_ticks == 0)
    {
        Initialize(result.num_ticks);
    }

    for(int i = 0; i<num_ticks; ++i)
    {
        f_n[i] += result.f_n[i];
        f_p[i] += result.f_p[i];
        t_n[i] += result.t_n[i];
        t_p[i] += result.t_p[i];
        que_num[i] += result.que_num[i];
        time_[i] += result.time_[i];
    }
	que += result.que;
	acc += result.acc;
	time += result.time;
	F1 += result.F1;

	this->std_acc += result.std_acc;
	this->std_F1 += result.std_F1;
	this->std_que += result.std_que;
	this->std_time += result.std_time;
    return *this;    
}

CResult& CResult::operator/(const double num_fold)
{
    
    for(int i = 0; i<num_ticks; ++i)
    {
        f_n[i] /= num_fold;
        f_p[i] /= num_fold;
        t_n[i] /= num_fold;
        t_p[i] /= num_fold;
        que_num[i] /= num_fold;
        time_[i] /= num_fold;	
    }
	que /= num_fold;
	acc /= num_fold;
	time /= num_fold;
	F1 /= num_fold;

	this->std_que /= num_fold;
	this->std_acc /= num_fold;
	this->std_time /= num_fold;
	this->std_F1 /= num_fold;
    return *this;    
}

void CResult::Reset()
{
    for(int t =1; t<num_ticks; t++)
    {
        f_n[t] = 0;
        f_p[t] = 0;
        t_n[t] = 0;
        t_p[t] = 0;
        que_num[t] = 0;
        time_[t] = 0.0;
    }
	que = 0.0;
	time = 0.0;
	acc = 0.0;
	F1 = 0.0;
	this->std_acc = 0.0;
	this->std_F1 = 0.0;
	this->std_que = 0.0;
	this->std_time = 0.0;
}
CResult::~CResult()
{
    delete []f_n;
    delete []f_p;
    delete []t_n;
    delete []t_p;
    delete []que_num;
    delete []time_;
}
