//
// Created by varsem on 3/16/24.
//
#include "functions.h"

int g2l(int i_glob, int m, int p)
{
  int i_glob_m = i_glob / m;
  int i_loc_m = i_glob_m / p;
  return i_loc_m * m + i_glob % m;
}

int l2g(int i_loc, int m, int p, int k)
{
  int i_loc_m = i_loc / m;
  int i_glob_m = i_loc_m * p + k;
  return i_glob_m * m + i_loc % m;
}

// максимальное число строк на процесс
int get_max_rows (int n, int m, int p) // 3 2 2
{
  // число блочных строк
  int b = (n + m - 1) / m; // 2
  return (b + p - 1) / p; // (2 + 3 - 1) / 2 = 2
}

// количество блочных строк в процессе
int get_rows (int n, int m, int p, int k) // 3 1 3 0
{
  // Число блочных строк
  int b = (n + m - 1) / m; // 3
  return b % p <= k ? b / p : b / p + 1; // 0
}

//В каком процессе лежит строка i_glob
int get_k (int m, int p, int i_glob)
{
  int i_glob_m = i_glob / m;
  return i_glob_m % p;
}

// количество блочных строк
int get_all_rows (int n, int m)
{
  return (n + m - 1) / m;
}

// сколько обычных строк в блочной строке под номером b
//int get_rows_n_in_b (int b, int n, in)

//сколько обычных значимых строк в процессе k
int get_small_rows_in_process (int n, int m, int p, int k)
{
  int rows_n = 0;
  int all_rows = get_all_rows (n, m);

  for (int i_glob = 0; i_glob < all_rows; i_glob++)
    {
      int owner = i_glob % p;
      if (k == owner)
        rows_n += min (m, n - i_glob * m);
    }
  return rows_n;
}