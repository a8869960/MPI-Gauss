//
// Created by varsem on 15.11.23.
//
#include "mpi.h"
#include <iostream>
#include <cstring>
#include <cstdio>
#include "math.h"

using namespace std;

struct main_block
{
    double inverse_norm = 0;
    int k = 0;

    int i_loc = 0;
    int j_loc = 0;
//    int i_glob = 0;
};

class LinearSystem
{
 public:
  double *a = nullptr;
  double *b = nullptr;
  double *x = nullptr;
  double *buf = nullptr;
  double *A = nullptr;
  double *B = nullptr;

  double *block = nullptr;
  double *block_inv = nullptr;
  double *block_h = nullptr;
  double *helper = nullptr;

  int *indi = nullptr;
  int *indj = nullptr;
  int *indi_m = nullptr;
  int *indj_m = nullptr;

  int n = 0;
  int m = 0;
  int k = 0;
  int p = 0;
  int rows = 0;

  main_block local_main_block;

  LinearSystem (int n_, int m_, int p_, int k_)
    {
      n = n_;
      m = m_;
      p = p_;
      k = k_;

      rows = (n + m - 1) / m;
      rows = (rows + p - 1) / p;

      a = new double [n * rows * m];
      b = new double [rows * m];
      x = new double [rows * m];
      buf = new double [n * m];
      A = new double [n * rows * m];
      B = new double [rows * m];

      block = new double [m * m];
      block_inv = new double [m * m];
      block_h = new double [m * m];
      helper = new double [m * m];

      indi = new int [n / m + 1];
      indj = new int [n / m + 1];
      indi_m = new int [m];
      indj_m = new int [m];

      memset (a, 0, sizeof (double) * n * m * rows);
      memset (b, 0, sizeof (double) * m * rows);

      for (int i = 0; i < n / m + 1; i++)
        indi[i] = i;
      for (int i = 0; i < n / m + 1; i++)
        indj[i] = i;
      for (int i = 0; i < m; i++)
        {
          indi_m[i] = i;
          indj_m[i] = i;
        }
      local_main_block.k = k;
    };
  ~LinearSystem ()
    {
      delete[] a;
      delete[] b;
      delete[] x;
      delete[] buf;
      delete[] A;
      delete[] B;

      delete[] block;
      delete[] block_inv;
      delete[] block_h;
      delete[] helper;

      delete[] indi;
      delete[] indj;
      delete[] indi_m;
      delete[] indj_m;
    };
};

enum class io_status
{
  success, // 0
  not_set, // 1
  no_matrix_main, // -1
  exists_matrix_main, // 2
  error_open, // -3
  error_read // -4
};

//mpi_functions
int g2l(int i_glob, int m, int p);
int l2g(int i_loc, int m, int p, int k);
int get_max_rows (int n, int m, int p);
int get_rows (int n, int m, int p, int k);
int get_k (int m, int p, int i_glob);
int get_small_rows_in_process (int n, int m, int p, int k);
int get_all_rows (int n, int m);

//matrix_input
io_status read_matrix (double *a, int n, int m, int p, int k, const char *filename, double *buf, MPI_Comm com);
void init_matrix (double *a, int n, int m, int k, int p, double (*f)(int, int, int));
double f1 (int n, int i, int j);
double f2 (int n, int i, int j);
double f3 (int n, int i, int j);
double f4 (int n, int i, int j);
void init_b (LinearSystem *S);
void get_local_block(double* A, double* block, int i_loc, int j_loc, int n, int m, int steps, int l, int k, int p);
void get_buf_block(double *A, double *block, int i, int j, int n, int m, int steps, int l);
void put_buf_block(double *A, double *block, int i, int j, int n, int m, int steps, int l);
void put_block(double* A, double* block, int i, int j, int n, int m, int k, int l, int g, int p);
void E(double* block, int m);
void get_block_b( double *B, double *block, int i, int m, int k, int l);
void put_block_b( double *B, double *block, int i, int m, int k, int l);



//matrix_output
void print_matrix (double *a, int n, int m, int p, int k, double *buf, int r, MPI_Comm com);
void print_status (io_status st);
void print_vector (double *c, int n, int m, int p, int k, double *buf, int r, MPI_Comm com);
void print_LinearSystem (LinearSystem *S, int r, MPI_Comm com);
void print_block (double *block, int m);

//for_gauss.cpp
double matrix_norm_of_system (LinearSystem *S, MPI_Comm com);
double matrix_norm(double *A, int n);
int find_local_block_main (LinearSystem *S, int step, int steps, int l, double NORM);
void rearrange_elements (LinearSystem *S, main_block most_main, int step, MPI_Comm com);
void return_to_origin_numeration (LinearSystem *S, int steps, MPI_Comm com);

//inverse_matrix.cpp
int inverse_matrix(double *a, double *A, double *B, int n, int *indi, int *indj, double norm);

//matrix_operations
int matrixSubtraction(double* A1, double *A2, double *C, int n, int m);
void matrix_product(double *A, double* B, double* C, int n, int s, int m);

//mpi_gauss.cpp
io_status mpi_gauss (LinearSystem *S, MPI_Comm com);

//results.cpp
double r_1 (LinearSystem *S, MPI_Comm com);
double r_2 (LinearSystem *S, MPI_Comm com);

