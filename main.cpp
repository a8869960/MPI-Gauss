//
// Created by varsem on 15.11.23.
//
#include "functions.h"
#include <fenv.h>

int main (int argc, char* argv[])
{
  feenableexcept (FE_ALL_EXCEPT ^ FE_INEXACT);
  int task = 11;
  int n, m, p, r, s, k;
  int main_k = 0;

  MPI_Init (&argc, &argv);
  MPI_Comm com = MPI_COMM_WORLD;
  MPI_Comm_size (com, &p);
  MPI_Comm_rank (com, &k);

  //Check input parameters
  if (!((   argc == 6
         && sscanf (argv[1], "%d", &n) == 1
         && sscanf (argv[2], "%d", &m) == 1
         && sscanf (argv[3], "%d", &r) == 1
         && sscanf (argv[4], "%d", &s) == 1
         && s == 0)
         ||(argc == 5
         && sscanf (argv[1], "%d", &n) == 1
         && sscanf (argv[2], "%d", &m) == 1
         && sscanf (argv[3], "%d", &r) == 1
         && sscanf (argv[4], "%d", &s) == 1
         && s != 0)
  ))
    {
      cout << "Wrong input parameters: -1" << endl;
      cout << argv[0] << " n m r s filename" << endl;
      MPI_Finalize ();
      return 0;
    }

  if (n < 1 || m < 1 || n < m || s < 0 || s > 4 || r < 0)
    {
      cout << "Wrong value of input parameters" << endl;
      MPI_Finalize ();
      return 0;
    }

  //Инициализация матрицы A
  int rows = get_max_rows (n, m, p);
  LinearSystem S (n, m, p, k);
  double (*g[])(int, int, int) = {f1, f2, f3, f4};

  if(argc == 5)
    init_matrix (S.a, n, m, k, p, g[s - 1]);
  else
    {
      io_status st = read_matrix (S.a, n, m, p, k, argv[5], S.buf, com);
      if (st != io_status::success)
        {
          if (k == main_k)
            print_status (st);
          MPI_Finalize ();
          return 0;
        }
    }

  //Инициализация матрицы B
  init_b (&S);

  //Вывод матриц
//  if (k == main_k && r != 0)
//    cout << "Matrix A:" << endl;
//  print_matrix (S.a, n, m, p, k, S.buf, r, com);
//  if (k == main_k && r != 0)
//    cout << "Matrix B:" << endl;
//  print_vector (S.b, n, m, p, k, S.buf, r, com);
    if (r != 0)
      print_LinearSystem (&S, r, com);

  memcpy (S.A, S.a, sizeof (double) * rows * n * m);
  memcpy (S.B, S.b, sizeof (double) * rows * m);

  //Gauss Method
  double t1 = MPI_Wtime ();
  io_status st = mpi_gauss (&S, com);
  if (k == main_k) print_status (st);
  t1 = MPI_Wtime () - t1;

  if (k == main_k && st == io_status::success && r != 0)
    cout << "Result:" << endl;
  if (r != 0 && st == io_status::success)
    print_vector (S.x, n, m, p, k, S.buf, r, com);

  // подсчет невязки
  double r1 = -1, r2 = -1;
  double t2 = MPI_Wtime ();
  if (st == io_status::success)
    {
      r1 = r_1 (&S, com);
      r2 = r_2 (&S, com);
    }
  t2 = MPI_Wtime () - t2;

  if (k == main_k)
    printf (
        "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
        argv[0], task, r1, r2, t1, t2, s, n, m, p);

  MPI_Finalize();

  return 0;
}

