long   N  = p_shape(0);
long   Nz = 0;

while(p(Nz) == 0) {
  ++Nz;
  p_assert(Nz);
 }

find_next_pulse(py_y, py_p, py_ret_i_flag, py_ret_A_b, i, i_max, 
                P, Q, th, N, Nz);

