/**
Created by weilin zhang on 20/02/2020.
Copyright © 2020 weiin. All rights reserved.
Description : This header file is mainly used to construct MPI_DataType,send and receive MPI_DataType. I also added the striped decomposition solver and grid decomposition solver to here. In addition, i also add output file function and other functions that are used to help me finish this assiginment.
Version : Final
*/
#pragma once
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
class DataType
{
public:
    bool periodic = false;
	int max_steps = 100;//max step 
	int test = 1;
	int rows, columns;//the number of devided grid data chunks
	DataType()
	{
		MPI_matrix_type = nullptr;
		zero = nullptr;
		zeros_horizontal = nullptr;
		zeros_vertical = nullptr;
	}
	DataType(int row_total, int column_total,int max_step): row_total(row_total),column_total(column_total),max_steps(max_step)
	{
		MPI_matrix_type = nullptr;
		zero = nullptr;
		zeros_horizontal = nullptr;
		zeros_vertical = nullptr;
		srand(time(NULL));
		data_total.resize(row_total,vector<int > (column_total));
		for (int i = 0; i < row_total; i++)
		{
			for (int j = 0; j < column_total; j++)
			{
				data_total[i][j] = rand()%2;
			}
			
		}
	}
	~DataType()
	{
		int flag;
		MPI_Finalized(&flag);
		if (!flag && MPI_matrix_type != nullptr)
			MPI_Type_free(MPI_matrix_type);
		delete MPI_matrix_type;
		delete zero;
		delete [] zeros_vertical;
		delete [] zeros_horizontal;
		
	}
	void Clear_MPI_Type()
	{
		if (MPI_matrix_type != nullptr)
			MPI_Type_free(MPI_matrix_type);
		delete MPI_matrix_type;
		MPI_matrix_type = nullptr;
	}
	void get_N_M(void);
	void devide_random_data(int id, int p);
    void Get_Devided_data(string filePath,int id, int p, int rows, int columns);
    void Send_left_edge(int destination, int tag_num, MPI_Request *request);
	void Receive_left_edge(int destination, int tag_num, MPI_Request *request);
	void Send_right_edge(int destination, int tag_num, MPI_Request *request);
	void Receive_right_edge(int destination, int tag_num, MPI_Request *request);
	void Send_top_edge(int destination, int tag_num, MPI_Request *request);
	void Receive_top_edge(int destination, int tag_num, MPI_Request *request);
	void Send_bottom_edge(int destination, int tag_num, MPI_Request *request);
	void Receive_bottom_edge(int destination, int tag_num, MPI_Request *request);
	void Send_left_top(int destination, int tag_num, MPI_Request *request);
	void Receive_left_top(int destination, int tag_num, MPI_Request *request);
	void Send_left_bottom(int destination, int tag_num, MPI_Request *request);
	void Receive_left_bottom(int destination, int tag_num, MPI_Request *request);
	void Send_right_top(int destination, int tag_num, MPI_Request *request);
	void Receive_right_top(int destination, int tag_num, MPI_Request *request);
	void Send_right_bottom(int destination, int tag_num, MPI_Request *request);
	void Receive_right_bottom(int destination, int tag_num, MPI_Request *request);


    int num_neighbours(int ii, int jj);
    int id_from_index(int id_row, int id_column);
    void print_data(void);
    void id_to_index(int id, int& id_row, int& id_column);
    void find_dimensions(int p);
    void do_iteration(void);
	void grid_to_file(int it, int process_row, int process_column);
	void param_for_video_tofile(void);
	void grid_solver(int id, int p,int id_row, int id_column);
	void striped_solver(int id,int p);
	vector<vector<int> > get_data(){return data;}

private:
	//I am using a vector to setup the matrix. It could also easily be done using dynamic allocation
	
	// int m, n;//the total rows and columns of the divied data
	int row_total, column_total;
    int M,N;
	int n,m;
    
    //the rows and columns of whole matrix
	vector<vector<int> > data;//data for 
    vector<vector<int> > new_data;
	 vector<vector<int> > data_total;
	int* zero;
	int* zeros_vertical;
	int* zeros_horizontal;

	MPI_Datatype* MPI_matrix_type;
	void setup_MPI_type_left_edge(void);
	void setup_MPI_type_left_edge_recv(void);
	void setup_MPI_type_right_edge(void);
	void setup_MPI_type_right_edge_recv(void);
	void setup_MPI_type_top_edge(void);
	void setup_MPI_type_top_edge_recv(void);
	void setup_MPI_type_bottom_edge(void);
	void setup_MPI_type_bottom_edge_recv(void);
	void setup_MPI_type_left_top(void);
	void setup_MPI_type_left_top_recv(void);
	void setup_MPI_type_left_bottom(void);
	void setup_MPI_type_left_bottom_recv(void);
	void setup_MPI_type_right_top(void);
	void setup_MPI_type_right_top_recv(void);
	void setup_MPI_type_right_bottom(void);
	void setup_MPI_type_right_bottom_recv(void);

	
	
};
