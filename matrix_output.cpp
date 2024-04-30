//
// Created by varsem on 15.11.23.
//
#include "functions.h"

int min (int r, int l);
int print_array (double *a, int n, int m, int printed_rows, int max_print);

void print_matrix (double *a, int n, int m, int p, int k, double *buf, int r, MPI_Comm com)
{
  int main_k = 0;
  int max_b = (n + m - 1) / m; // количество блочных строк
  int printed_rows = 0;

  for (int b = 0; b < max_b; b++)
    {
      int owner = b % p;
      int rows = min (m, n - b * m);
      int b_loc = b / p;

      if (k == main_k)
        {
          if (owner == main_k)
            printed_rows += print_array (a + b_loc * n * m, n, rows, printed_rows, r);
          else
            {
              MPI_Status st;
              MPI_Recv (buf, n * rows, MPI_DOUBLE, owner, 0, com, &st);
              printed_rows += print_array (buf, n, rows, printed_rows, r);
            }
        }
      else
        {
          if (owner == k)
            MPI_Send (a + b_loc * n * m, n * rows, MPI_DOUBLE, main_k, 0, com);
        }
    }
}

int print_array (double *a, int n, int m, int printed_rows, int max_print)
{
  if (printed_rows >= max_print)
    return 0;

  int p_n = (n > max_print ? max_print : n);
  int p_m = (printed_rows + m < max_print ? m : max_print - printed_rows);

  for (int i = 0; i < p_m; i++)
    {
      for (int j = 0; j < p_n; j++)
        printf (" %10.3e", a[i * n + j]);
      printf ("\n");
    }
  return p_m;
}

void print_vector (double *c, int n, int m, int p, int k, double *buf, int r, MPI_Comm com)
{
  int main_k = 0;
  int max_b = (n + m - 1) / m; // количество блочных строк
  int printed_elements = 0;

  for (int b = 0; b < max_b; b++)
    {
      int owner = b % p; // номер процесса - владельца
      int rows = min (m, n - b * m); // кол-во обычных строк в блочной строке
      int b_loc = b / p; // номер блочной строки в локальной нумерации

      if (k == main_k)
        {
          if (owner == main_k)
            {
              for (int i = b_loc * m; i < b_loc * m + rows; i++)
                {
                  if (printed_elements >= r)
                    break;
                  printf (" %10.3e", c[i]);

                  printed_elements++;
                }
            }
          else
            {
              MPI_Status st;
              MPI_Recv (buf, n, MPI_DOUBLE, owner, 0, com, &st);
              for (int i = 0; i < rows; i++)
                {
                  if (printed_elements >= r)
                    break;
                  printf (" %10.3e", buf[i]);

                  printed_elements++;
                }
            }
        }
      else
        {
          if (k == owner)
            MPI_Send (c + b_loc, rows, MPI_DOUBLE, main_k, 0, com);
        }
    }
  if(k == main_k)
    printf("\n");
}

int min (int r, int l)
{
  return (r < l) ? r : l;
}

void print_status (io_status st)
{
  switch (st)
    {
      case io_status::error_open:
        cout << "Can't open the file" << endl;
      break;

      case io_status::error_read:
        cout << "Can't read element" << endl;
      break;

      case io_status::no_matrix_main:
        cout << "No matrix main" << endl;
      break;

      case io_status::exists_matrix_main:
        cout << "Exists matrix main" << endl;
      break;

      case io_status::success:
        break;

      case io_status::not_set:
        break;
    }
}

void print_LinearSystem (LinearSystem *S, int r, MPI_Comm com)
{
  int main_k = 0;

  if (S->k == main_k)
    cout << "Matrix A:" << endl;
  print_matrix (S->a, S->n, S->m, S->p, S->k, S->buf, r, com);
  if (S->k == main_k)
    cout << "Matrix B:" << endl;
  print_vector (S->b, S->n, S->m, S->p, S->k, S->buf, r, com);
}

void print_block (double *block, int m)
{
  for (int ii = 0; ii < m; ii++)
    {
      for (int jj = 0; jj < m; jj++)
        printf(" %.3e", block[ii * m + jj]);
      cout << endl;
    }
    cout << endl;
}