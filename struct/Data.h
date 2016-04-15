#ifndef DATA_H
#define DATA_H

/*
 * =====================================================================================
 *
 *       Filename:  Data.h
 *
 *    Description:  head of data struction file, supporting 
 *
 *        Version:  1.0
 *        Created:  24/02/2014 19:00:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuji Hao (), hao.shuji@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
 #include <iostream>
 #include <string>
 #include <stdlib.h>
 #include <math.h>
 #include "Global.h" 


 using namespace std;

struct CFeature_node{
    int index;
    double value;
};

class CData{
public:
    string data_file_type; //
    //num_elements: number of non-zero feature in all instances
    int num_classes, d,n,num_elements;
	int num_pos, num_neg;
    int max_line_len;
    int *ntances_each_class;
    char *line;
    double *y;// the target values
    struct CFeature_node **x; //each row represents one instance,
    struct CFeature_node *x_space;
	//Matrix X; //each column is an instance
public:
    CData();
    int Import_data(string data_file_name);
    int Read_libsvm(string data_file_name);
    void Change_to_Matrix();
    int Read_arff(string data_file_name);
	void Norm2One();
    void Multiclass_to_binary (double one_class);
    char* Read_line(FILE *input);
    ~CData();
};

#endif
