//
// Created by varvara.semenova on 4/3/24.
//
#include "functions.h"

double norm_1 (double *y, int n, int m, int p, int k, MPI_Comm com);
void Ax (LinearSystem *S, MPI_Comm com);

double norm_1 (double *y, int n, int m, int p, int k, MPI_Comm com)
{
  double local_sum = 0, global_sum;
  int rows = get_small_rows_in_process (n, m, p, k);

  for (int i = 0; i < rows; i++)
    local_sum += fabs (y[i]);

  MPI_Allreduce (&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);
  return global_sum;
}

void Ax (LinearSystem *S, MPI_Comm com)
{
  double *a = S->A, *x = S->x, *buf = S->buf;
  int n = S->n, m = S->m, p = S->p, k = S->k;
  int rows = get_small_rows_in_process (n, m, p, k);

  int steps = n / m;
  int l = n - steps * m;
  int bl = (l == 0 ? steps : steps + 1);

  for (int i_glob = 0; i_glob < bl; i_glob++)
    {
      int owner = i_glob % p;
      int size = (i_glob == steps ? l : m);
      int i_loc = g2l (i_glob * m, m, p);

      if (k == owner)
        memcpy (buf + i_glob * m, x + i_loc, size * sizeof (double));

      MPI_Bcast (buf + i_glob * m, size, MPI_DOUBLE, owner, com);
    }

  matrix_product (a, buf, S->b, rows, n, 1);
}

double r_1 (LinearSystem *S, MPI_Comm com)
{
  double *b = S->b, *B = S->B;
  int n = S->n, k = S->k, m = S->m, p = S->p;
  int rows = get_small_rows_in_process (n, m, p, k);
  double norm = matrix_norm_of_system (S, com);

  Ax (S, com);
  matrixSubtraction (b, B, b, rows, 1); // Ax - b

  double nominator = norm_1 (b, n, m, p, k, com);
  double denominator = norm_1 (B, n, m, p, k, com);

  if (fabs (denominator) < 1e-16 * norm)
    return -2;

  return nominator / denominator;
}

double r_2 (LinearSystem *S, MPI_Comm com)
{
  double *x = S->x;
  int n = S->n, k = S->k, m = S->m, p = S->p;
  int rows = get_small_rows_in_process (n, m, p, k);
  double local_sum = 0, global_sum;

  for (int i = 0; i < rows; i++)
    {
      int i_glob = l2g (i, m, p, k);
      local_sum += fabs (x[i] - (i_glob + 1) % 2);
    }

  MPI_Allreduce (&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);

  return global_sum;
}