DATA_SECTION

  init_int n
  init_vector x(1,40)

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number r(1e-05,100)
  init_bounded_number k(1e-05,100)
  init_bounded_number c(1e-05,100)
  init_bounded_number s(1e-05,100)
  vector mu(1,n) // per capita mort prob
PROCEDURE_SECTION
  mu = x + r * elem_prod((1 - x / k), (x - c) / k);
  f = 0.5 * n * log(2 * M_PI) + n * log(s) + 0.5 * norm2(x - mu) / square(s);

