//
// Created by varsem on 15.11.23.
//
#define eps 1.2e-16

#include "functions.h"

using namespace std;

int matrixMax(double *A, int k, int n, int *indi, int *indj, double norm);

int inverse_matrix(double *a, double *A, double *B, int n, int *indi, int *indj, double norm)
{
  if(n == 1)
    {
      if(abs(a[0]) < 1e-15 * norm)
        return -1;

      A[0] = 1. / a[0];
      return 0;
    }

  int i, j, ii;

  for(i = 0; i < n; i++)
    {
      indi[i] = i;
      indj[i] = i;
    }

  memcpy(A, a, sizeof(double) * n * n);

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      {
        if(i == j)
          B[i * n + j] = 1;
        else
          B[i * n + j] = 0;
      }

  //Прямой ход
  for(i = 0; i < n; i++)
    {
      if(matrixMax(A, i, n, indi, indj, norm) == -1)
        return -1;

      double Aii = A[indi[i] * n + indj[i]];

      for(j = 0; j < i + 1; j++)
        B[indi[i] * n + indj[j]] = B[indi[i] * n + indj[j]] / Aii;
      for(j = i + 1; j < n; j++)
        {
          A[indi[i] * n + indj[j]] = A[indi[i] * n + indj[j]] / Aii;
          B[indi[i] * n + indj[j]] = B[indi[i] * n + indj[j]] / Aii;
        }
      A[indi[i] * n + indj[i]] = 1;

      for(ii = i + 1; ii < n; ii++)
        {
          for(j = 0; j < i + 1; j++)
            B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - B[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
          for(j = i + 1; j < n; j++)
            {
              A[indi[ii] * n + indj[j]] = A[indi[ii] * n + indj[j]] - A[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
              B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - B[indi[i] * n + indj[j]] * A[indi[ii] * n + indj[i]];
            }
          A[indi[ii] * n + indj[i]] = 0;
        }
    }

  //Обратный ход
  for(i = n - 1; i >= 0; i--)
    {
      for(ii = i - 1; ii >= 0; ii--)
        {
          double Aa = A[indi[ii] * n + indj[i]];

          for(j = n - 1; j >= 0; j--)
            {
              B[indi[ii] * n + indj[j]] = B[indi[ii] * n + indj[j]] - Aa * B[indi[i] * n + indj[j]];
              A[indi[ii] * n + indj[j]] = A[indi[ii] * n + indj[j]] - Aa * A[indi[i] * n + indj[j]];
            }
        }
    }

  for(i = 0; i < n; i++)
    indi[i] = i;

  for(j = 0; j < n; j++)
    for(i = 0; i < n; i++)
      {
        if(abs(A[indi[i] * n + j] - 1) < eps)
          {
            int helper = indi[i];
            indi[i] = indi[j];
            indi[j] = helper;
          }
      }

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      A[i * n + j] = B[indi[i] * n + j];

  return 0;
}

int matrixMax(double *A, int k, int n, int *indi, int *indj, double norm)
{
  double max = 0;
  int imax = k, jmax = k;

  for(int i = k; i < n; i++)
    for(int j = k; j < n; j++)
      {
        if(abs(A[indi[i] * n + indj[j]]) > max)
          {
            max = abs(A[indi[i] * n + indj[j]]);
            imax = i;
            jmax = j;
          }
      }

  if(abs(A[indi[imax] * n + indj[jmax]]) < eps * norm)
    return -1;

  int helper;

  helper = indi[imax];
  indi[imax] = indi[k];
  indi[k] = helper;

  helper = indj[jmax];
  indj[jmax] = indj[k];
  indj[k] = helper;

  return 0;
}