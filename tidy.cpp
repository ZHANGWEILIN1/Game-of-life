/*
Created by weilin zhang on 20/02/2020.
Copyright © 2020 weiin. All rights reserved.
Description : This is the main function , run this code you can test the solver ,generate arbitrary 
              sizes of input random data, get the final evolution result.
Version : Final
*/
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include"DataType.h"
using namespace std;

// //initialization
bool result = true;
bool verify = false;
int id, p;
int rows, columns;
int id_row, id_column;

vector< vector < int > > grid;

/** this function is used to test whether my solver is correct i 
	@param grid  the evolution data of each processor
	@param m1 the rows of grid decomposition
	@param n1 the columns of grid decomposition
	@param periodic the information of different mode(periodic or not)
	@param result the information to show us whether the result is correct
output:
    whether the result is true or false
*/
bool test_solution(vector<vector<int> > grid,int &m1,int &n1,bool periodic, bool result);
/** this function is used to judge which kind of data we want to use(arbitrary sizes random data or the data i store in the data folder that to test the solver)
	@param verify  if it is true test the data else not 
	@param row_total the total rows of whole data
	@param col_total the total columns of whole data
	@param max_step the iteration steps
output:
       DataType
*/
DataType test_or_random(bool verify,int row_total,int col_total, int max_step);

int main(int argc, char* argv[])
{
	clock_t start_time;
	DataType A = test_or_random(verify,100,100,100);

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 1000);
    //set rows and columns
    A.find_dimensions(p);

	// //find_dimensions(p, rows, columns);
	A.param_for_video_tofile();
    ///// get id_column and id_row
    A.id_to_index(id, id_row, id_column);
	if (verify)
	{
		A.Get_Devided_data("./data/2078_source.txt",id, p, rows, columns);
	}
	else
	{
		A.devide_random_data(id,p);
	}
    A.get_N_M();
	if (columns == 1)
	{
		start_time = clock();
		A.striped_solver(id,p);
		double time_different = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<time_different<<"s"<<endl;
	}
	else
	{
		start_time = clock();
		A.grid_solver(id,p,id_row,id_column);
		double time_different = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<time_different<<"s"<<endl;
	}
	
	
	if (verify)
	{
		int m1,n1;
		bool q;
		grid = A.get_data();
		q = test_solution(grid,m1,n1,A.periodic,result);
		if (q == true)
		{
			cout<<"The result is correct"<<endl;
		}
		else
		{
			cout<<"The solution is wrong"<<endl;
		}
	}
	
	
	
	
	
	MPI_Finalize();
}



bool test_solution(vector<vector<int> > grid,int &m1,int &n1,bool periodic,bool result)
{
	DataType solution;
    vector< vector < int > >  new_grid;
	if(periodic)
	{
		solution.Get_Devided_data("./data/2078_solution_p.txt", m1, n1, id, p);
	}
	else
	{
        solution.Get_Devided_data("./data/2078_solution_np.txt", m1, n1, id, p);
	}
    new_grid = solution.get_data();

	for (int i = 1; i < m1-1; i++)
	{
		for (int j = 1; j < n1-1; j++)
		{
			if(new_grid[i][j] != grid[i][j])
			{
				result = false;
				break;
			}

		}
		
	}
	return result;
}
DataType test_or_random(bool verify,int row_total,int col_total, int max_step)
{
	if (verify)
	{
		
		DataType A;
		return A;
	}
	else
	{
		DataType A(row_total,col_total,max_step);
		return A;
	}
}