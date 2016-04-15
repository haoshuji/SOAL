/*
 * =====================================================================================
 *
 *       Filename:  Data.cpp
 *
 *    Description:  file of data struction
 *
 *        Version:  1.0
 *        Created:  24/02/2014 19:26:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shuji Hao (), hao.shuji@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#define _CRT_SECURE_NO_DEPRECATE
#include "Data.h"
#include <string.h>
#include <iostream>
#include <fstream>
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

CData::CData()
{
    data_file_type = std::string();
    num_classes = 0;
    ntances_each_class = NULL;
    d = 0;
    n =0;
    max_line_len = 1024;
    line = NULL;
    x_space = NULL;
    y = NULL;
    x = NULL;    
	num_pos = 0;
	num_neg = 0;
}

int CData::Import_data(std::string data_file_name)
{    
	cout << "Importing data\t" << data_file_name << endl;
	Read_libsvm(data_file_name);
	cout << "\t# Ins.=" << n << "\t # Pos.=" << num_pos << "\t # Neg.=" << num_neg << "\t # Dim.=" << d << endl;
    return 1;
}

void CData::Change_to_Matrix()
{
   
}

int CData::Read_libsvm(std::string file_name)
{
    max_line_len=1024;
    FILE *fp = NULL;
	fp = fopen(file_name.c_str(),"r");
    line = (char*)malloc(sizeof(char)*max_line_len); 
    int elements = 0;

    if(fp == NULL){
        std::cout << "Unable to open file\t" << file_name << endl;
		exit(1);
        return 0;
    }

    while(Read_line(fp) != NULL){
        char *p = strtok(line, " \t"); //label value doesn't account
        while(1){
            p = strtok(NULL, " \t");
            if(p == NULL || *p== '\n')
                break;
            elements++;
        } 
        elements++;
        this->n++;
    }
    rewind(fp);
    
    y = Malloc(double, n);
    x = Malloc(struct CFeature_node *, n);
    x_space = Malloc(struct CFeature_node, elements);

    int inst_max_index, max_index = 0,j=0;
    char *label,*endptr,*idx, *val;
   // double label_f = 0.0;
    for(int i=0; i<n; i++){
        inst_max_index = 0; 
        Read_line(fp);
        x[i] = &x_space[j];
        label = strtok(line, " \t\n");
        if(label == NULL)
            std::cout<<i+1<<"line error"<<std::endl;
        double label_f = strtod(label, &endptr);
        y[i] = label_f;
		if (label_f > 0.0)
		{
			num_pos++;		
		}
		else
		{
			num_neg++;
		}
        if(endptr == label || *endptr != '\0')
            std::cout<<i+1<<"line error"<<std::endl;

        while(1){
           idx = strtok(NULL, ":");
           val = strtok(NULL, " \t");

           if(val == NULL)
               break;
           int f_no = 0;
           x_space[j].index = (int)strtol(idx, &endptr, 10);
           if(endptr == idx || f_no != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
                std::cout<<i+1<<"line error"<<std::endl;
           else
               inst_max_index = x_space[j].index;

           f_no = 0;
           x_space[j].value = strtod(val, &endptr);
           if(endptr == val || f_no != 0 || (*endptr != '\0' && !isspace(*endptr)))
                std::cout<<i+1<<"line error"<<std::endl;

           ++j;
        }
        
        if(inst_max_index > max_index)
            max_index = inst_max_index;
        
        x_space[j++].index = -1;

    }
    
    this->d = max_index;
    fclose(fp);

    free(line); 
    return 1;
}

void CData::Norm2One()
{	
    double length;
    struct CFeature_node *x_t;
    for(int i=0;i<n;i++)
    {
        x_t = x[i];
        length = 0.0;
        while(x_t->index != -1)
        {
            length += x_t->value * x_t->value;  
            x_t++;
        }
        length = sqrt(length);
        if(length > 0.0)
        {
            x_t = x[i];

            while(x_t->index != -1)
            {
                x_t->value = x_t->value/length;
                x_t++;
            }
        }
    }
	cout << "Normalized data" << endl;
}

char* CData::Read_line(FILE *input)
{
    int len;

    if(fgets(line, max_line_len, input) == NULL)
        return NULL;
    while(strrchr(line, '\n') == NULL)
    {
        max_line_len *= 2;
        line = (char*) realloc(line, max_line_len);
        len = (int) strlen(line);
        if(fgets(line+len, max_line_len-len, input) == NULL)
            break;
    }
    return line;
}

int CData::Read_arff(std::string file_name)
{
    return 1;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Multiclass_to_binary
 *  Description:  Change the multiclass to binary
 *  Parameters:   one_class: the index of class user defined as class with lable 1, and other classes are change to class with label -1
 * =====================================================================================
 */
void CData::Multiclass_to_binary(double one_class)
{
    int ntance_one_label=0;

    for ( int i=0; i<n; i++ ) {
        if((double)y[i] == one_class){
            y[i] = 1;
            ntance_one_label++;
        }
        else
            y[i] = -1;
    }

    if(ntance_one_label == 0)
        std::cerr<<"Without class with label \t" <<one_class << std::endl;
    else{
		std::cout << "Done: multi-classes to binary class!" << endl;
        /*ntances_each_class[0] = ntance_one_label;
        ntances_each_class[1] = n - ntance_one_label;*/
    }
}

CData::~CData()
{
	if (x!=NULL)
	{
		free(x);
	}
	if (x_space != NULL)
	{
		free(x_space);
	}
	if (y != NULL)
		free(y);
}
