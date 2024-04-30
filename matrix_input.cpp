//
// Created by varsem on 15.11.23.
//
#include "functions.h"

int read_array (FILE *fp, double *a, int len);

io_status read_matrix (double *a, int n, int m, int p, int k, const char *filename, double *buf, MPI_Comm com)
{
  int main_k = 0; // номер главного процесса

  FILE *fp = nullptr;
  int err = 0; // номер ошибки

  if (k == main_k)
    {
      fp = fopen (filename, "r");
      if (fp == nullptr)
        err = 1;
    }

  MPI_Bcast (&err, 1, MPI_INT, main_k, com);
  if (err == 1)
    return io_status::error_open;

  memset (buf, 0, n * m * sizeof (double));

  int max_b = (n + m - 1) / m; // число блочных строк
  for (int b = 0; b < max_b; b++)
    {
      int owner = b % p; // владелец строки
      int rows = (b * m + m <= n ? m : n - b * m); // размер этой блочной строки
      int b_loc = b / p; // локальный номер блочной строки

      if (k == main_k)
        {
          err += read_array (fp, buf, n * rows);
          if (owner == main_k)
            memcpy (a + b_loc * m * n, buf, n * rows * sizeof (double));
          else
            MPI_Send (buf, n * rows, MPI_DOUBLE, owner, 0, com);
        }
      else if (k == owner) // получает строку
        {
          MPI_Status st;
          MPI_Recv (a + b_loc * m * n, n * rows, MPI_DOUBLE, main_k, 0, com, &st);
        }
    }
    if (k == main_k)
      {
        fclose (fp);
        fp = nullptr;
      }

    MPI_Bcast (&err, 1, MPI_INT, main_k, com);
    if (err != 0)
      return io_status::error_read;

  return io_status::success;
}

int read_array (FILE *fp, double *a, int len)
{
  for (int i = 0; i < len; i++)
    {
      if (fscanf (fp, "%lf", a++) != 1)
        return -2;
    }
  return 0;
}

void init_matrix (double *a, int n, int m, int k, int p, double (*f)(int, int, int))
{
  int i_loc, i_glob, j_loc, j_glob;
  int rows = get_small_rows_in_process (n, m, p, k); //количество блочных строк в процессе

  for (i_loc = 0; i_loc < rows; i_loc++)
    {
      i_glob = l2g (i_loc, m, p, k);

      for (j_loc = 0; j_loc < n; j_loc++)
        {
          j_glob = j_loc;
          a[(i_loc) * n + j_loc] = (*f) (n, (i_glob), j_glob);
        }
    }
}

double f1 (int n, int i, int j)
{
  return n - max(i + 1, j + 1) + 1;
}

double f2 ([[maybe_unused]]int n, int i, int j)
{
  return max(i + 1, j + 1);
}

double f3 ([[maybe_unused]]int n, int i, int j)
{
  return abs(i - j);
}

double f4 ([[maybe_unused]]int n, int i, int j)
{
  return 1. / (i + j + 1);
}

int max (int i, int j)
{
  return (i > j) ? i : j;
}

void init_b (LinearSystem *S)
{
  double *b = S->b, *a = S->a;
  int n = S->n, m = S->m, rows_n = S->rows * m;

  for (int i_loc = 0; i_loc < rows_n; i_loc++)
    {
      for (int j_loc = 0; j_loc < n; j_loc += 2)
        {
          b[i_loc] += a[i_loc * n + j_loc];
        }
    }
}

void get_local_block(double *A, double *block, int i, int j,
               int n, int m, int steps, int l, int k, int p)
{
  int block_m = (l2g (i * m, m, p, k) / m == steps ? l : m), block_l = (j == steps ? l : m);
  int r, s;
  int a = i * n * m + j * m; //number of first element of the block
  memset(block, 0, sizeof(double) * m * m);

//  if (block_l == 2 and block_m == 3 )
//    cout << "flag" << endl;

  for(r = 0; r < block_m; r++)
    for(s = 0; s < block_l; s++)
      {
//        if (block_l == 2 and block_m == 3)
//          cout << a + r * n + s << " " << A[a + r * n + s] << " " << r * block_l + s << endl;
        block[r * block_l + s] = A[a + r * n + s];
      }
}

void get_buf_block(double *A, double *block, int i, int j,
                     int n, int m, int steps, int l)
{
  int block_m = m, block_l = (j == steps ? l : m);
  int r, s;
  int a = i * n * m + j * m; //number of first element of the block
  memset(block, 0, sizeof(double) * m * m);

  for(r = 0; r < block_m; r++)
    for(s = 0; s < block_l; s++)
      {

        block[r * block_l + s] = A[a + r * n + s];
      }
}

void put_buf_block(double *A, double *block, int i, int j,
                   int n, int m, int steps, int l)
{
  int block_m = m, block_l = (j == steps ? l : m);
  int r, s;
  int a = i * n * m + j * m; //number of first element of the block

  for(r = 0; r < block_m; r++)
    for(s = 0; s < block_l; s++)
      {
        A[a + r * n + s] = block[r * block_l + s];
      }
}

void put_block(double* A, double* block, int i, int j,
               int n, int m, int k, int l, int g, int p)
{
  int block_m = (l2g (i * m, m, p, g) / m == k ? l : m), block_l = (j == k ? l : m);

  int r, s;
  int a = i * n * m + j * m; //number of first element of the block

  for(r = 0; r < block_m; r++)
    {
      for(s = 0; s < block_l; s++)
        {
          A[a + r * n + s] = block[r * block_l + s];
        }
    }
}

void E(double* block, int m)
{
  memset(block, 0, sizeof(double) * m * m);

  for(int i = 0; i < m; i++)
    block[i * m + i] = 1;
}

void get_block_b( double *B, double *block, int i, int m, int k, int l)
{
  int block_m = (i == k ? l : m);

  memset(block, 0, sizeof(double) * m * m);

  int r;
  int b = i * m; //number of first element of the block

  for(r = 0; r < block_m; r++)
    {
      block[r] = B[b + r];
    }
}

void put_block_b( double *B, double *block, int i, int m, int k, int l)
{
  int block_m = (i == k ? l : m);

  int r;
  int b = i * m; //number of first element of the block

  for(r = 0; r < block_m; r++)
    B[b + r] = block[r];
}