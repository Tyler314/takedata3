
inline void update_D_F(PyArrayObject *py_y, PyArrayObject *py_p, 
                       long i, long N, long Nz, double *D, double *F) {
  long j;
  y_assert(i - 1);
  y_assert(i + N - 1);

  *D += -y(i - 1) + y(i + N - 1);

  *F = 0;
  for(j = Nz; j < N; ++j) {
    y_assert(i + j);
    p_assert(j);
    *F += y(i + j) * p(j);
  }
  
}

inline void calc_A_b(double P, double Q, double D, double F, long N, 
                     double *A, double *b) {
  *b = (D*Q - F*P) / (N*Q - P*P);
  *A = (F - P*(*b)) / Q;
}

inline void calc_err_flag(PyArrayObject *py_y, PyArrayObject *py_p,
                          long i, long N, long Nz, double th, 
                          double A, double b, double *E, long *flag) {
  long j;
  double res;

  *flag = 0;
  *E = 0;

  for(j = 0; j < N; ++j) {
    y_assert(i + j); 
    p_assert(j);
    res = A * p(j) + b - y(i + j);
    *E += res * res;
    if(j > Nz && *flag == 0 && (res > th || res < -th)) {
      *flag = 1;
    }
  }
}


inline void calc_next(PyArrayObject *py_y, PyArrayObject *py_p, 
                      long i, long N, long Nz, double P, double Q,
                      double *D, double *F, double *A, double *b, double *E, 
                      long *flag, double th) {
  update_D_F(py_y, py_p, i, N, Nz, D, F);
  calc_A_b(P, Q, *D, *F, N, A, b);
  calc_err_flag(py_y, py_p, i, N, Nz, th, *A, *b, E, flag);
}


void find_next_pulse(PyArrayObject *py_y, PyArrayObject *py_p,
                     PyArrayObject *py_ret_i_flag, PyArrayObject *py_ret_A_b,
                     long i, long i_max, double P, double Q, double th, 
                     long N, long Nz) {
  long j; // Loop index. 
  
  // Fitting variables. 
  double D = 0, F = 0;
  
  double A = 0, Ap = 0, An = 0;
  double b = 0, bp = 0, bn = 0;
  double E = 0, Ep = 0, En = 0;
  
  long flag = 0, flagn = 0, flagp = 0;
  
  // Initialize D and F.
  for(j = 0; j < N; ++j) {
    y_assert(i + j);
    p_assert(j);
    
    D += y(i + j);
    F += y(i + j) * p(j);
  }
  
  calc_A_b(P, Q, D, F, N, &Ap, &bp);
  
  // Loop through indices until Ap > th. 
  while(Ap < th && i < i_max - 3) {
    ++i;
    update_D_F(py_y, py_p, i, N, Nz, &D, &F);
    calc_A_b(P, Q, D, F, N, &Ap, &bp);
  }
  
  calc_err_flag(py_y, py_p, i, N, Nz, th, Ap, bp, &Ep, &flagp);
  
  // At this point we have the previous point. Now we need the current and next.
  ++i; 
  calc_next(py_y, py_p, i, N, Nz, P, Q, &D, &F, &A, &b, &E, &flag, th);
  ++i; 
  calc_next(py_y, py_p, i, N, Nz, P, Q, &D, &F, &An, &bn, &En, &flagn, th);
  
  while(i < i_max - 1) {
    
    // Check break condition. 
    if(E < En && E < Ep) {
      ret_i_flag(0) = i - 1;
      ret_i_flag(1) = flag;
      ret_A_b(0) = A;
      ret_A_b(1) = b;
      break;
    }
    
    // Not a local minimum, so move on. 
    Ap = A;  bp = b;  Ep = E;  flagp = flag;
    A  = An;  b = bn;  E = En; flag  = flagn;
    
    ++i;
    calc_next(py_y, py_p, i, N, Nz, P, Q, &D, &F, &An, &bn, &En, &flagn, th);
  }
}                     
                     
