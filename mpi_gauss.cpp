//
// Created by varsem on 3/16/24.
//
#include "functions.h"

io_status mpi_gauss (LinearSystem *S, MPI_Comm com)
{
  int n = S->n, m = S->m, p = S->p, k = S->k, step;
  double *a = S->a, *b = S->b, *buf = S->buf;
  double *block = S->block, *block_inv = S->block_inv, *block_h = S->block_h;
  int *indi_m = S->indi_m, *indj_m = S->indj_m, *indj = S->indj;
  int steps, l, bl;
  int st_loc = 1, st_glob = 0, rows = get_rows (n, m, p, k);
  bool is_print = false;
//   [[maybe_unused]]double time = 0, t;

  steps = n / m; //how many blocks m*m
  l = n - steps * m; //how long last block
  bl = (l == 0 ? steps : steps + 1); //number of all blocks
//  t = bl % p;
  if (l == 0) l = m;
  double norm = matrix_norm_of_system (S, com);

  //Прямой ход
  for (step = 0; step < steps; step++)
    {
      int owner_of_main_row = step % p;
      
      //Нахождение главного элемента
      st_loc = find_local_block_main (S, step, steps, l, norm);
      if (p != 1)
        MPI_Allreduce (&st_loc, &st_glob, 1, MPI_INT, MPI_MAX, com); // если ни в одном нет главного, то вернет не 0
      else
        st_glob = st_loc;
      if (st_glob != 0)
        return io_status::no_matrix_main;

      //Нахождение "самого" главного
      main_block most_main;
      int main_block_k = 0;
      if (p != 1)
        {
          MPI_Allreduce (&(S->local_main_block), &most_main, 1, MPI_DOUBLE_INT, MPI_MINLOC, com);
          main_block_k = most_main.k;
          MPI_Bcast (&S->local_main_block.i_loc, 2, MPI_INT, main_block_k, com);
        }
      most_main.i_loc = S->local_main_block.i_loc;
      most_main.j_loc = S->local_main_block.j_loc;

      if (is_print)
        {
          if (k == most_main.k)
            cout << "MAIN BLOCK: " << most_main.i_loc << " " << most_main.j_loc << " "
                 << a[most_main.i_loc * m * n + indj[most_main.j_loc]] << " " << main_block_k << endl;
        }

      if (most_main.i_loc < 0 || most_main.j_loc < 0)
        return io_status::no_matrix_main;

      //Переставляем элементы
      rearrange_elements (S, most_main, step, com);
     

      //Делим "самую" главную строку
      if (k == owner_of_main_row)
        {
          //Обращение главного элемента
          int i_loc = g2l (step * m, m, p) / m;
          get_local_block (a, block, i_loc, indj[step], n, m, steps, l, k, p);
          inverse_matrix (block, block_inv, block_h, m, indi_m, indj_m, norm);

          //Деление строчки матрицы А
          E (block, m);
          put_block (a, block, i_loc, indj[step], n, m, steps, l, k, p);
          for (int j = step + 1; j < bl; j++)
            {
              int size_l = (j != steps ? m : l);

              get_local_block(a, block, i_loc, indj[j], n, m, steps, l, k, p);
              matrix_product(block_inv, block, block_h, m, m, size_l);
              put_block(a, block_h, i_loc, indj[j], n, m, steps, l, k, p);
            }

          //Деление строчки матрицы В
          get_block_b(b, block, i_loc, m, steps, l);
          matrix_product(block_inv, block, block_h, m, m, 1);
          put_block_b(b, block_h, i_loc, m, steps, l);
        }

      if (is_print)
        {
          MPI_Allreduce (&(S->local_main_block), &most_main, 1, MPI_DOUBLE_INT, MPI_MINLOC, com);
          if (k == 0)
            cout << endl << "DIVISION STEP: " << step << endl;
          print_LinearSystem (S, 10, com);
        }

      //ЗАНУЛЕНИЕ СТОЛБЦА
//       t = MPI_Wtime ();
      //Обмениваемся "самой" главной строчкой
      int index = g2l (step * m, m, p) / m;
     if (p == 1)
       buf = a + index * n * m;
     else
       {
          for (int j = step + 1; j < bl; j++)
            {
              if (k == owner_of_main_row)
                get_local_block (a, block, index, indj[j], n, m, steps, l, k, p);
              MPI_Bcast (block, m * m, MPI_DOUBLE, owner_of_main_row, com);
              put_buf_block (buf, block, 0, indj[j], n, m, steps, l);
            }
       }

      int ind = (step % p < k ? step / p : step / p + 1);

      //Обмениваемся присоединенным элементом
      double *bb = S->helper;
      if (k == step % p)
        {
          get_block_b (b, bb, g2l (step * m, m, p) / m, m, steps, l);
        }
      if (p != 1)
        MPI_Bcast (bb, m * m, MPI_DOUBLE, step % p, com);
      
      for(int i = ind; i < rows; i++)
        { 
            get_local_block(a, block, i, indj[step], n, m, steps, l, k, p); // size = size_m * m

            for(int j = step + 1; j < bl; j++)
            {
                int size_m = (l2g (i * m, m, p, k) / m != steps ? m : l);
                int size_l = (j == steps ? l : m);

                get_buf_block(buf, block_h, 0, indj[j], n, m, steps, l); //size = m * size_l
                matrix_product(block, block_h, block_inv, size_m, m, size_l);
                get_local_block(a, block_h, i, indj[j], n, m, steps, l, k, p);
                matrixSubtraction(block_h, block_inv, block_h, size_m, size_l);
                put_block(a, block_h, i, indj[j], n, m, steps, l, k, p);
            }
            matrix_product(block, bb, block_inv, m, m, 1);

            get_block_b(b, block_h, i, m, steps, l);
            matrixSubtraction(block_h, block_inv, block_h, 1, m);
            put_block_b(b, block_h, i, m, steps, l);
        }

      memset(block, 0, sizeof(double) * m * m);
      for(int i = ind; i < get_rows (n, m, p, k); i++)
        put_block(a, block, i, indj[step], n, m, steps, l, k, p);
//  t = MPI_Wtime () - t;
//       time += t;
    if (is_print)
      {
        MPI_Allreduce (&(S->local_main_block), &most_main, 1, MPI_DOUBLE_INT, MPI_MINLOC, com);
        if (k == 0)
          cout << endl << "DIFFERENCE STEP: " << step << endl;
        print_LinearSystem (S, 10, com);
      }
    }
//   MPI_Bcast (&st_loc, 1, MPI_INT, 0, com);

  //Делим последний блок
  if(bl == steps + 1 && k == steps % p)
    {
      get_local_block(a, block, rows - 1, indj[steps], n, m, steps, l, k, p);
      if(inverse_matrix(block, block_inv, block_h, l, indi_m, indj_m, norm) == -1)
        st_loc = -1;
      else
        {
          E(block, l);
          put_block(a, block, rows - 1, indj[steps], n, m, steps, l, k, p);

          get_block_b(b, block, rows - 1, m, k, l);
          matrix_product(block_inv, block, block_h, l, l, 1);
          put_block_b(b, block_h, rows - 1, m, k, l);
        }
    }
  if (bl == steps + 1)
    {
      MPI_Bcast (&st_loc, 1, MPI_INT, steps % p, com);
      if (st_loc == -1)
        return io_status::no_matrix_main;
    }

  if (is_print)
    {
      if (k == 0)
        cout << endl << "DIV STEP: " << step << endl;
      print_LinearSystem (S, 10, com);
      if (k == 0)
        cout << endl << endl << endl << endl << endl << endl;
    }
    
    
  //ОБРАТНЫЙ ХОД
  int ind = rows;
  if(bl == steps + 1)
    {
      int owner = steps % p;
      if (k == owner)
        {
          ind--;
          get_block_b (b, block, rows - 1, m, steps, l); //size = 1 * m or m * 1
        }
      if (p != 1)
        MPI_Bcast (block, m * l, MPI_DOUBLE, owner, com);

      for(int i = ind - 1; i >= 0; i--)
        {
          get_local_block(a, block_h, i, indj[steps], n, m, steps, l, k, p); //size = m * l
          matrix_product(block_h, block, block_inv, m, l, 1); //size = m * 1 or 1 * m
          get_block_b(b, block_h, i, m, steps, l); //size = 1 * m or m * 1
          matrixSubtraction(block_h, block_inv, block_h, 1, m);
          put_block_b(b, block_h, i, m, steps, l);
        }
    }

  for(step = steps - 1; step >= 0; step--)
    {
      int owner = step % p;
      if (k == owner)
        {
          ind--;
          get_block_b (b, block, g2l (step * m, m, p) / m, m, steps, l); //size = 1 * m or m * 1
        }
      if (p != 1)
        MPI_Bcast (block, m * m, MPI_DOUBLE, owner, com);

      for(int i = ind - 1; i >= 0; i--)
        {
          get_local_block(a, block_h, i, indj[step], n, m, steps, l, k, p); //size = m * m
          matrix_product(block_h, block, block_inv, m, m, 1); //size = m * 1 or 1 * m
          get_block_b(b, block_h, i, m, steps, l); //size = 1 * or m * 1
          matrixSubtraction(block_h, block_inv, block_h, 1, m); //
          put_block_b(b, block_h, i, m, steps, l);
        }

      if (is_print)
        {
//          int q = 0, t;
//          MPI_Allreduce (&t, &q, 1, MPI_DOUBLE, MPI_SUM, com);
          if (k == owner)
            cout << "INV " << owner << endl;
          print_LinearSystem (S, 10, com);
        }
    }

  if (is_print)
    {
      if (k == 0)
        cout << endl << "END: " << step << endl;
      print_LinearSystem (S, 10, com);
    }

  return_to_origin_numeration (S, steps, com);
//   cout << time << endl;

  return io_status::success;
}
