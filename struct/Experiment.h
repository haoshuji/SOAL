#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "Result.h"

#include "Data.h"
#include "Result.h"
#include "Parameter.h"
#include "Model.h"
#include "ModelPE.h"
#include "ModelPA.h"
#include "ModelSOP.h"
#include "ModelAROW.h"
#include "ModelAROWD.h"
#include "ModelAROWC.h"
#include "ModelAROWCD.h"
#include "ModelSOPD.h"
#include "Global.h"

#include <algorithm> //max min function
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

class CExperiment
{
    public:
        CParameter par;
        CData data;
		string data_dir;
        string data_name;
		string data_fullpath;
		string output_dir;
		string setting_file_fullpath;
        
		
    public:
        /* ====================  LIFECYCLE     ======================================= */                             /* constructor */
        CExperiment();
        void Parse_command_line(int argc, char **argv);

        void Exit_with_help();       

        void Run_binary_experiment(int argc, char **argv);          
        
        /* ====================  ACCESSORS     ======================================= */

        /* ====================  MUTATORS      ======================================= */

        /* ====================  OPERATORS     ======================================= */

    protected:
        /* ====================  METHODS       ======================================= */

        /* ====================  DATA MEMBERS  ======================================= */

    private:
        /* ====================  METHODS       ======================================= */

        /* ====================  DATA MEMBERS  ======================================= */

public:
	void Find_best_parameter();	
}; /* -----  end of class CExperiment  ----- */

#endif
