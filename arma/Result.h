#ifndef RESULT_H
#define RESULT_H

/*
 * =====================================================================================
 *        Class:  CResult
 *  Description:  store the temporal result of algorithm
 * =====================================================================================
 */
class CResult
{
    public:
        double *f_n, *f_p, *t_n, *t_p, *que_num, *time_;
		double que, acc, time, F1;
		double std_que, std_acc, std_time, std_F1;
        int num_ticks;
        /* ====================  LIFECYCLE     ======================================= */
        CResult();
        CResult(int num_ticks);
        void Initialize(int m_num_ticks); 
        void Reset();                           /* constructor */
        ~CResult();
        
        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        /* ====================  OPERATORS     ======================================= */
        CResult& operator=(const CResult &result);
        CResult& operator+(const CResult &result);
        CResult& operator/(const double num_fold);
    protected:
        /* ====================  METHODS       ======================================= */

        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  METHODS       ======================================= */

        /* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class CResult  ----- */

#endif
