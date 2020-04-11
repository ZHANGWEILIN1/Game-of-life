
/* this version is used to high performance test,because the codes contain class ,if i want to hpc i have to intergerate all the codes  into one file */
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
// #include"DataType.h"

using namespace std;
class DataType
{
public:
    bool periodic = true;
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
void DataType::setup_MPI_type_left_edge(void)
{
	if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(N);
		vector<MPI_Aint> addresses(N);
		vector<MPI_Datatype> typelist(N);

		int cnt = 0;
		if (periodic == true || test == 1)
		{
			for (int i = 0; i < N; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[i+1][1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
		}
		else
		{
			for (int i = 0; i < N; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[i+1][0], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
		}
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}


void DataType::setup_MPI_type_left_edge_recv(void)
{
	if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(N);
		vector<MPI_Aint> addresses(N);
		vector<MPI_Datatype> typelist(N);

		int cnt = 0;
		for (int i = 0; i < N; i++)
		{
			block_length[cnt] = 1;
			MPI_Aint temp;
			MPI_Get_address(&data[i+1][0], &temp);
			addresses[cnt] = temp;
			typelist[cnt] = MPI_INT;
			cnt++;
		}
		MPI_matrix_type = new MPI_Datatype;

		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_edge(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(N);
		vector<MPI_Aint> addresses(N);
		vector<MPI_Datatype> typelist(N);

		int cnt = 0;
		if (periodic == true || test == 1)
		{
			for (int i = 0; i < N; i++)
		   {
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[i+1][M], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
		   }
			
		}
		else
		{
			for (int i = 0; i < N; i++)
		   {
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[i+1][M+1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
		   }
		}
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_edge_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(N);
		vector<MPI_Aint> addresses(N);
		vector<MPI_Datatype> typelist(N);

		int cnt = 0;
		for (int i = 0; i < N; i++)
		{
			block_length[cnt] = 1;
			MPI_Aint temp;
			MPI_Get_address(&data[i+1][M+1], &temp);
			addresses[cnt] = temp;
			typelist[cnt] = MPI_INT;
			cnt++;
		}
		MPI_matrix_type = new MPI_Datatype;

		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}

void DataType::setup_MPI_type_top_edge(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(M);
		vector<MPI_Aint> addresses(M);
		vector<MPI_Datatype> typelist(M);

		int cnt = 0;
		if (periodic == true|| test ==1)
		{
			for (int i = 0; i < M; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[1][i+1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
		}
		else
		{
			for (int i = 0; i < M; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[0][i+1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
			
		}
		
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_top_edge_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(M);
		vector<MPI_Aint> addresses(M);
		vector<MPI_Datatype> typelist(M);

		int cnt = 0;
		for (int i = 0; i < M; i++)
		{
			block_length[cnt] = 1;
			MPI_Aint temp;
			MPI_Get_address(&data[0][i+1], &temp);
			addresses[cnt] = temp;
			typelist[cnt] = MPI_INT;
			cnt++;
		}
		MPI_matrix_type = new MPI_Datatype;

		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType:: setup_MPI_type_bottom_edge(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(M);
		vector<MPI_Aint> addresses(M);
		vector<MPI_Datatype> typelist(M);

		int cnt = 0;
		if (periodic == true || test == 1)
		{
			for (int i = 0; i < M; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[N][i+1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
		}
		else
		{
			for (int i = 0; i < M; i++)
			{
				block_length[cnt] = 1;
				MPI_Aint temp;
				MPI_Get_address(&data[N+1][i+1], &temp);
				addresses[cnt] = temp;
				typelist[cnt] = MPI_INT;
				cnt++;
			}
		}
		
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_bottom_edge_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(M);
		vector<MPI_Aint> addresses(M);
		vector<MPI_Datatype> typelist(M);

		int cnt = 0;

		for (int i = 0; i < M; i++)
		{
			block_length[cnt] = 1;
			MPI_Aint temp;
			MPI_Get_address(&data[N+1][i+1], &temp);
			addresses[cnt] = temp;
			typelist[cnt] = MPI_INT;
			cnt++;
		}
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(cnt, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_left_top(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
		if (periodic == true || test == 1)
		{
			MPI_Get_address(&data[1][1], &temp);
		}
		else
		{
			MPI_Get_address(&data[0][0], &temp);
		}
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}

void DataType::setup_MPI_type_left_top_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
        MPI_Get_address(&data[0][0], &temp);
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_left_bottom(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
		if (periodic == true|| test == 1)
		{
			MPI_Get_address(&data[N][1], &temp);
		}
		else
		{
			MPI_Get_address(&data[N+1][0], &temp);
		}
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_left_bottom_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
        MPI_Get_address(&data[N+1][0], &temp);
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_top(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
		if (periodic == true || test == 1)
		{
			MPI_Get_address(&data[1][M], &temp);
		}
		else
		{
			MPI_Get_address(&data[0][M+1], &temp);
		}
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_top_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
        MPI_Get_address(&data[0][M+1], &temp);
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_bottom(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
		if (periodic == true || test == 1)
		{
			MPI_Get_address(&data[N][M], &temp);
		}
		else
		{
			MPI_Get_address(&data[N+1][M+1], &temp);
		}
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}
void DataType::setup_MPI_type_right_bottom_recv(void)
{
    if (MPI_matrix_type == nullptr)
	{
		vector<int> block_length(1);
		vector<MPI_Aint> addresses(1);
		vector<MPI_Datatype> typelist(1);

        block_length[0] = 1;
        MPI_Aint temp;
        MPI_Get_address(&data[N+1][M+1], &temp);
        addresses[0] = temp;
        typelist[0] = MPI_INT;
		MPI_matrix_type = new MPI_Datatype;
		
		MPI_Type_create_struct(1, &block_length[0], &addresses[0], &typelist[0], MPI_matrix_type);
		MPI_Type_commit(MPI_matrix_type);
	}
}

void DataType::Send_left_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_edge();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
	
}

void DataType::Receive_right_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_edge_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}

void DataType::Send_right_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_edge();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_left_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_edge_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_top_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_top_edge();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_bottom_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_bottom_edge_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_bottom_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_bottom_edge();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_top_edge(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_top_edge_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_left_top(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_top();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_right_bottom(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_bottom_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_right_bottom(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_bottom();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_left_top(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_top_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_right_top(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_top();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_left_bottom(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_bottom_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Send_left_bottom(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_left_bottom();
    MPI_Isend(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request); 
}
void DataType::Receive_right_top(int destination, int tag_num, MPI_Request *request)
{
	setup_MPI_type_right_top_recv();
    MPI_Irecv(MPI_BOTTOM, 1, *MPI_matrix_type, destination, tag_num, MPI_COMM_WORLD, request);
}
void DataType::Get_Devided_data(string filePath,int id, int p, int rows, int columns)
{
    ifstream dataFile;
	dataFile.open(filePath, std::ifstream::in);
	int row_total, column_total, row_counter = 0, start_row, start_column;
	
	if (!dataFile) {
		cerr << "error message" << endl;
	}

	
	dataFile >> this->row_total >> this->column_total;
	m = this->row_total / this->rows;
	n = this->column_total / this->columns;
	start_row = (id / this->columns) * m;
	start_column = (id % this->columns) * n;
	if (id / this->columns < this->row_total % this->rows)
	{
		m++;
		start_row += (id / this->columns);
	}
	else
	{
		start_row += this->row_total % this->rows;
	}
	if (id % this->columns < this->column_total % this->columns)
	{
		n++;
		start_column += (id % this->columns);
	}
	else
	{
		start_column += this->column_total % this->columns;
	}
	m += 2;
	n += 2;
	this->data.resize(m, vector<int>(n));
	int temp;
	while (!dataFile.eof())
	{
		if (row_counter < start_row)
		{
			for (int i = 0; i < this->column_total; i++)
			{
				dataFile >> temp;
			}

		}
		else
		{
			for (int j = 0; j < this->column_total; j++)
			{
				if ((j<start_column)|| j>=(start_column + n-2))
				{
					dataFile >> temp;
				}
				else
				{
					dataFile >> this->data[row_counter - start_row + 1][j - start_column + 1];
				}
			}
		}
		row_counter++;
		if (row_counter == start_row - 2 + m)
		{
			break;
		}
	}


	
}
int DataType::num_neighbours(int ii, int jj)
{
    int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + m) % m;
				jx = (j + jj + n) % n;
				if (this->data[ix][jx]) cnt++;
			}
	return cnt;
}
int DataType::id_from_index(int id_row, int id_column)
{
    if (periodic == true) // find periodic neighbours id
	{
		return ((id_row + rows) % rows) * columns + (id_column + columns) % columns;
	}
	if (periodic == false) 
    {
		if (id_row >= rows || id_row < 0)
		{
			test = -1;
		}
		else if (id_column >= columns || id_column < 0)
		{
			test = -1;
		}
		else
		{
			test = 1;
		}
		//return id_row * columns + id_column;
		return ((id_row + rows) % rows) * columns + (id_column + columns) % columns;

	}
}
void DataType::print_data()
{
    for (int i = 0; i < m ; i++) {
		for (int j = 0; j < n; j++) {
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
}
void DataType::id_to_index(int id, int& id_row, int& id_column)
{
    id_column = id % this->columns;
	id_row = id / this->columns;
}
void DataType::find_dimensions(int p)
{
    columns = sqrt(p);
	rows = p/this->columns;
	while (p%this->columns != 0)
	{
		this->columns--;
		this->rows = p/this->columns;
	}
	
	
}
void DataType:: do_iteration(void)
{
    new_data.resize(m,vector<int>(n));
	for (int i = 1; i < m-1; i++)
		for (int j = 1; j < n-1; j++)
		{
			new_data[i][j] = data[i][j];
			int num_n = num_neighbours(i, j);
			if (data[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_data[i][j] = (int)0;
			}
			else if (num_n == 3) new_data[i][j] = (int)1;
		}
	//shrink vector to fit and swap the new_grid to grid
	data.swap(new_data);

}
void DataType::param_for_video_tofile(void)
{
	stringstream fname;
	fstream f1;

	fname <<"/Users/zhangweilin/Desktop/final/import_data/"<< "parameter_for_video"<< ".txt";
	f1.open(fname.str().c_str(), ios_base::out);
	f1<<this->max_steps<<endl;
	f1<<this->rows<<endl;
	f1<<this->columns<<endl;

	f1.close();
}
void DataType::grid_to_file(int it, int process_row, int process_column)
{
	stringstream fname;
	fstream f1;

	fname <<"/Users/zhangweilin/Desktop/final/import_data/"<< "output" << "_" << "iteration_" << it << "_processrow_" << process_row << "_processcolumn_" << process_column << ".txt";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < m - 1; i++) //dont output boundaries received from neighbour processes
	{
		for (int j = 1; j < n - 1; j++)
			f1 <<this->data[i][j] << "\t";
		f1 << endl;
	}

	f1.close();
}
void DataType::get_N_M(void)
{
	this->N = this->m-2;
	this->M = this->n-2;
}
void DataType::grid_solver(int id, int p,int id_row, int id_column)
{
	 MPI_Request* request;
    request= new MPI_Request[16];
	// do max_steps iterations
	for (int i = 0; i < max_steps; i++)
	{
        
		
		///////////////////////////1. send left/////////////////////////////////////////
		//get the present id of processor 
        id_to_index(id, id_row, id_column);
		//shift into the left neighbour of present id
        int idsend = id_from_index(id_row, id_column - 1);
        Send_left_edge(idsend,1,&request[0]);
        Clear_MPI_Type();
		///////////////////////////////////receive from right//////////////////////////
		
		int idrecv = id_from_index(id_row, id_column + 1); 
        Receive_right_edge(idrecv,1,&request[1]);
        Clear_MPI_Type();
		/////////////////////////////2. send right///////////////////////////////////
        
        idsend = id_from_index(id_row, id_column + 1);
        Send_right_edge(idsend,2,&request[2]);
        Clear_MPI_Type();
		///////////////////////////////receive from left///////////////////////////
        
        idrecv = id_from_index(id_row, id_column - 1);
        Receive_left_edge(idrecv,2,&request[3]);
        Clear_MPI_Type();
		//////////////////////////3. send top//////////////////////////////////////
		
		idsend = id_from_index(id_row - 1, id_column);
        Send_top_edge(idsend,3,&request[4]);
        Clear_MPI_Type();
	
		///////////////////receive from bottom///////////////////////
		
		idrecv = id_from_index(id_row + 1, id_column);
        Receive_bottom_edge(idrecv,3,&request[5]);
        Clear_MPI_Type();


		////////////////////////////4. send bottom////////////////////
		
		idsend = id_from_index(id_row + 1, id_column);
        Send_bottom_edge(idsend,4,&request[6]);
        Clear_MPI_Type();
		/////////////////////////receive from top/////////////////////////////////
		
		idrecv = id_from_index(id_row - 1, id_column);Receive_top_edge(idrecv,4,&request[7]);
        Clear_MPI_Type();

		///////////////////////5. send bottom right edge//////////////////////
		
		idsend = id_from_index(id_row + 1, id_column + 1);
        Send_right_bottom(idsend,5,&request[8]);
        Clear_MPI_Type();
		////////////////////receive from top left//////////////////
		
		idrecv = id_from_index(id_row - 1, id_column - 1); 
        Receive_left_top(idrecv,5,&request[9]);
        Clear_MPI_Type();
		///////////////////////6. send top left edge/////////////////
		
		idsend = id_from_index(id_row - 1, id_column - 1);
        Send_left_top(idsend,6,&request[10]);
        Clear_MPI_Type();
		/////////////////////receive from bottom right edge//////////////
		
		idrecv = id_from_index(id_row + 1, id_column + 1); 
        Receive_right_bottom(idrecv,6,&request[11]);
        Clear_MPI_Type();
	

		//////////////////////////////7. send top right edge/////////////////////////////////
		
		idsend = id_from_index(id_row - 1, id_column + 1);
        Send_right_top(idsend,7,&request[12]);
        Clear_MPI_Type();
		////////////////////////receive from bottom left edge////////////////////////////
		
		idrecv = id_from_index(id_row + 1, id_column - 1);
        Receive_left_bottom(idrecv,7,&request[13]);
        Clear_MPI_Type();
		///////////////////////////8. send bottom left edge///////////////////////////////
		
		idsend = id_from_index(id_row + 1, id_column - 1);
        Send_left_bottom(idsend,8,&request[14]);
        Clear_MPI_Type();

		//receive from top right edge
		
		
		idrecv = id_from_index(id_row - 1, id_column + 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
        Receive_right_top(idrecv,8,&request[15]);
        Clear_MPI_Type();

		//----------

		// MPI_Waitall(16, request, MPI_STATUS_IGNORE);

		

		// do iteration...
		//----------------------------------------------------------
		do_iteration();
		// id_to_index(id, id_row, id_column);
		// grid_to_file(i, id_row, id_column); // write output grid to a file
		// // //----------------------------------------------------------

		
			
	}
}
void DataType::striped_solver(int id,int p)
{
	MPI_Request* request;
    request= new MPI_Request[4];
	for (int i = 0; i < max_steps; i++)
	{
		if (periodic)
	   {
		    MPI_Isend(&data[m-2][1], n-2, MPI_INT, (id + 1) % p, 0, MPI_COMM_WORLD, &request[0]);
            MPI_Isend(&data[1][1], n-2, MPI_INT, (p + id - 1) % p, 1, MPI_COMM_WORLD, &request[1]);

            MPI_Irecv(&data[0][1], n-2, MPI_INT, (p + id - 1) % p, 0, MPI_COMM_WORLD, &request[2]);
            MPI_Irecv(&data[m-1][1], n-2, MPI_INT, (id + 1) % p, 1, MPI_COMM_WORLD, &request[3]);
			for (int i = 1; i < m-1; i++)
		   {
			   data[i][0] = data[i][n-2];
			   data[i][n-1] = data[i][1];  
		   }
		   data[0][0] = data[0][n-2];
		   data[0][n-1] = data[0][1];
		   data[m-1][0] = data[m-1][n-2];
		   data[m-1][n-1] = data[m-1][1];
	   }
	   else
	   {
		   for (int i = 0; i < m; i++)
		   {
			   data[i][0] = 0;
			   data[i][n-1] = 0;
		   }
		   if(id < p - 1) MPI_Isend(&data[(m - 2) * n], n-2, MPI_INT, id + 1, 0, MPI_COMM_WORLD, &request[0]);
            if(id > 0) MPI_Isend(&data[n], n-2, MPI_INT, id - 1, 1, MPI_COMM_WORLD, &request[1]);

            if(id > 0) MPI_Irecv(&data[0], n-2, MPI_INT, id - 1, 0, MPI_COMM_WORLD, &request[2]);
            else for(int j = 1; j < n-1; j++) data[0][j] = 0;
            if(id < p - 1) MPI_Irecv(&data[(m - 1) * n], n-2, MPI_INT, id + 1, 1, MPI_COMM_WORLD, &request[3]);
            else for(int j = 1; j < n-1; j++) data[m - 1][j] = 0;
	   }
	   
	}
	do_iteration();
}
void DataType::devide_random_data(int id, int p)
{
	int start_row,start_column;
	m = this->row_total / this->rows;
	n = this->column_total / this->columns;
	start_row = (id / this->columns) * m;
	start_column = (id % this->columns) * n;
	if (id / this->columns < this->row_total % this->rows)
	{
		m++;
		start_row += (id / this->columns);
	}
	else
	{
		start_row += this->row_total % this->rows;
	}
	if (id % this->columns < this->column_total % this->columns)
	{
		n++;
		start_column += (id % this->columns);
	}
	else
	{
		start_column += this->column_total % this->columns;
	}
	int end_row = start_row + m;
	int end_col = start_column + n;
	m += 2;
	n += 2;
	data.resize(m,vector<int> (n));
	for (int i = start_row; i < end_row; i++)
	{
		for (int j = start_column; j < end_col; j++)
		{
			data[i-start_row+1][j-start_column+1] = data_total[i][j];
		}
		
	}
}

// //initialization
bool result = true;
bool verify = false;
int id, p;
int rows, columns;
int id_row, id_column;

vector< vector < int > > grid;
//------------------- Game of Life Serial -----------------------

bool test_solution(vector<vector<int> > grid,int &m1,int &n1,bool periodic, bool result);
DataType test_or_random(bool verify,int row_total,int col_total, int max_step);

int main(int argc, char* argv[])
{
	clock_t start_time;
	DataType A = test_or_random(verify,2000,2000,100);

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
		A.Get_Devided_data("./data/30_source.txt",id, p, rows, columns);
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
		solution.Get_Devided_data("./data/30_solution_p.txt", m1, n1, id, p);
	}
	else
	{
        solution.Get_Devided_data("./data/30_solution_np.txt", m1, n1, id, p);
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