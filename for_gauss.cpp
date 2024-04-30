//
// Created by varsem on 15.11.23.
//
#include "functions.h"
#define eps 1e-15

void change_rows_in_one_process (LinearSystem *S, int i1_glob, int i2_glob, bool is_all = true);
void change_rows_in_two_processes (LinearSystem *S, int i1_glob, int i2_glob, MPI_Comm com, bool is_all = true);

double matrix_norm_of_system (LinearSystem *S, MPI_Comm com)
{
  double norm = 0, local_sum = 0, global_sum;
  double *a = S->A;
  int n = S->n, m = S->m, k = S->k, p = S->p;
  int rows_n = get_small_rows_in_process (n, m, p, k);

  for (int j = 0; j < n; j++)
    {
      local_sum = 0;
      for (int i = 0; i < rows_n; i++)
        local_sum += fabs (a[i * n + j]);

      MPI_Allreduce (&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);
      if (global_sum > norm)
        norm = global_sum;
    }

  return norm;
}

double matrix_norm(double *A, int n)
{
  double norm = 0, helper = 0;

  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        helper += abs (A[i * n + j]);

      if (helper > norm)
        norm = helper;

      helper = 0;
    }
  return norm;
}

int find_local_block_main (LinearSystem *S, int step, int steps, int l, double NORM)
{
  double *a = S->a, *block = S->block, *block_inv = S->block_inv, *block_h = S->block_h;
  int *indi_m = S->indi_m, *indj_m = S->indj_m;
  int n = S->n, m = S->m, k = S->k, p = S->p;
  int i_loc, count = 0, imax = -1, jmax = -1;
  double norm, min = 1.7976931348623158e+308;

  i_loc = (step % p <= k ? step / p : step / p + 1);
  int b = get_rows (n, m, p, k);

  for (; i_loc < b; i_loc++)
    {
      for (int j_loc = step; j_loc < steps; j_loc++)
        {
          get_local_block (a, block, i_loc, S->indj[j_loc], n, m, steps, l, k, p);
          if (inverse_matrix (block, block_inv, block_h, m, indi_m, indj_m, NORM) == 0)
            {
              norm = matrix_norm (block_inv, m);
              if (norm < min)
                {
                  min = norm;
                  imax = i_loc;
                  jmax = j_loc;
                }
              else
                count++;
            }
        }
    }
  if (count == (steps - step) * (steps - step))
    return -1;

  S->local_main_block.i_loc = imax;
  S->local_main_block.j_loc = jmax;
  S->local_main_block.inverse_norm = min;

  return 0;
}

void rearrange_elements (LinearSystem *S, main_block most_main, int step, [[maybe_unused]]MPI_Comm com)
{
  int helper, m = S->m, p = S->p, k = S->k;;

  // Перестановка столбцов
  helper = S->indj[step];
  S->indj[step] = S->indj[most_main.j_loc];
  S->indj[most_main.j_loc] = helper;

  //Перестановка строк
  int i_glob = l2g (most_main.i_loc * m, m, p, most_main.k) / m;

  helper = S->indi[step];
  S->indi[step] = S->indi[i_glob];
  S->indi[i_glob] = helper;

  //Меняем строки между процессами
  //Глобальные номера строк, которе надо поменять - i_glob и step
  if (i_glob == step)
    return;

  int owner1 = i_glob % p;
  int owner2 = step % p;
  if (k != owner1 && k != owner2)
    return;
//  cout << owner1 << owner2 << endl;
///////////
  //Теперь (k == owner1) || (k == owner2)
  if (owner1 == owner2)
    change_rows_in_one_process (S, i_glob, step);
  else // owner1 != owner2
    change_rows_in_two_processes (S, i_glob, step, com);
}

void change_rows_in_one_process (LinearSystem *S, int i1_glob, int i2_glob, bool is_all)
{
  int n = S->n, m = S->m, k = S->k, p = S->p;
  double *a = S->a, *b = S->b, *block = S->block, *block_h = S->block_h;
  int *indj = S->indj;
  int step = i2_glob;
  int steps = n / m;
  int bl = (n - steps * m == 0 ? steps : steps + 1);
  int l = n - steps * m;
  if (l == 0) l = m;
  int i1_loc = g2l (i1_glob * m, m, p) / m;
  int i2_loc = g2l (i2_glob * m, m, p) / m;
  int pos1 = n * m * i1_loc, pos2 = n * m * i2_loc;
  double helper;

  if (i1_glob == i2_glob)
    return;

  if (is_all)
    {
      for (int im = step; im < bl; im++)
        {
          get_local_block (a, block, i1_loc, indj[im], n, m, steps, l, k, p);
          get_local_block (a, block_h, i2_loc, indj[im], n, m, steps, l, k, p);
          put_block (a, block, i2_loc, indj[im], n, m, steps, l, k, p);
          put_block (a, block_h, i1_loc, indj[im], n, m, steps, l, k, p);
        }
    }
  pos1 = m * i1_loc;
  pos2 = m * i2_loc;
  for (int i = 0; i < m; i++)
    {
      helper = b[pos1 + i];
      b[pos1 + i] = b[pos2 + i];
      b[pos2 + i] = helper;
    }
}

void change_rows_in_two_processes (LinearSystem *S, int i1_glob, int i2_glob, MPI_Comm com, bool is_all)
{
  int n = S->n, m = S->m, k = S->k, p = S->p;
  int *indj = S->indj;
  int step = i2_glob;
  int steps = n / m;
  int bl = (n - steps * m == 0 ? steps : steps + 1);
  int l = n - steps * m;
  if (l == 0) l = m;
  int owner2 = (i1_glob % p == k ? i2_glob % p : i1_glob % p);
  int i_glob = (i1_glob % p == k ? i1_glob : i2_glob);
  double *a = S->a, *b = S->b, *buf = S->buf, *block = S->block;
  MPI_Status st;

//  cout << i1_glob << i2_glob << endl;
  int i_loc = g2l (i_glob * m, m, p) / m;
  int pos = n * m * i_loc;
  //Change a_rows
//  is_all = false;
  if (is_all)
    {
      for (int j = step; j < bl; j++)
        {
          get_local_block (a, block, i_loc, indj[j], n, m, steps, l, k, p);
          MPI_Send (block, min (m * m, (n - m * step) * (n - m * step)), MPI_DOUBLE, owner2, 0, com);
          MPI_Recv (block, min (m * m, (n - m * step) * (n - m * step)), MPI_DOUBLE, owner2, 0, com, &st);
          put_block (a, block, i_loc, indj[j], n, m, steps, l, k, p);
        }
    }
  //Change b
  pos = i_loc * m;
  MPI_Send (b + pos, m, MPI_DOUBLE, owner2, 0, com);
  MPI_Recv (buf, m, MPI_DOUBLE, owner2, 0, com, &st);
  memcpy (b + pos, buf, m * sizeof (double));
}

void return_to_origin_numeration (LinearSystem *S, int steps, MPI_Comm com)
{
  double *b = S->b, *x = S->x;
  int n = S->n, m = S->m, k = S->k, p = S->p;
  int *indj = S->indj;
  int rows = get_small_rows_in_process (n, m, p, k);

  if (rows == 0)
    return;

  if (m == n)
    {
      for (int i = 0; i < n; i++)
        x[i] = b[i];
      return;
    }


  for (int t = 0; t < steps; t++)
    {
      int j;
      for (j = t; j < steps; j++)
        if (indj[j] == t)
          break;
      int i_glob = j;
      int owner1 = t % p, owner2 = i_glob % p;

      if (t != i_glob && (k == owner1 || k == owner2))
        {
          if (owner1 == owner2)
            change_rows_in_one_process (S, t, i_glob, false);
          else
            change_rows_in_two_processes (S, i_glob, t, com, false);
        }

      int helper = indj[t];
      indj[t] = t;
      indj[j] = helper;
    }

  for (int i = 0; i < rows; i++)
    x[i] = b[i];
}
