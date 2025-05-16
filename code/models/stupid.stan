//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  vector[2] d13Ca;
  vector[2] s1;
  vector[2] s2;
  vector[2] s3;
  vector[2] s4;
  vector[2] s5;
  vector[2] Ci0;
  vector[2] A0;
  vector[2] gb;
  real b;
  vector[2] d13Cp;
  vector[2] Dab;
  vector[2] GCLab;
  vector[2] GCWab;
  real gam;
}

transformed data{
  real d_v = 0.000940096;
  real s2_beta = s2[1] / s2[2] ^ 2;
  real amax_v = (s3[1] * (1 - s3[1])) / s3[2] ^ 2 - 1;
  real gc_v = (s4[1] * (1 - s4[1])) / s4[2] ^ 2 - 1;
  real meso_v = (s5[1] * (1 - s5[1])) / s5[2] ^ 2 - 1;
  real Ci0_beta = Ci0[1] / Ci0[2] ^ 2;
  real A0_beta = A0[1] / A0[2] ^ 2;
  real gb_beta = gb[1] / gb[2] ^ 2;
  real s1_beta = s1[1] / s1[2] ^ 2;
}

// The parameters accepted by the model.
parameters {
  real<lower=100, upper=8000> ca;
  real<lower=-12, upper=-4> d13Ca_m;
  real<lower=0> s2_m;
  real<lower=0, upper=1> amax_scale;
  real<lower=0, upper=1> gc_scale;
  real<lower=0, upper=1> meso_scale;
  real<lower=1> Ci0_m;
  real<lower=0> A0_m;
  real<lower=0> gb_m;
  real<lower=1e-7> l;
  real<lower=0> gcl_m;
  real<lower=0, upper=1> SA_gc;
  real<lower=0> s1_m;
}

transformed parameters{
  real Pl = gcl_m * 0.5e-6;
  real SA = SA_gc * 0.5 ^ 2;
  real D = SA / (pi() * (Pl / 2) ^ 2);

  real amax = SA * amax_scale / D ;
  real gcmax = (d_v * D * amax) / (l + ((pi() / 2) * sqrt(amax / pi()))) / 1.6;
  real gcop = (1 / (gcmax * gc_scale) + 1 / gb_m) ^ -1;

  real q_a = 1 / gcop * (Ci0_m - gam);
  real q_b = gam * (-2 * A0_m / gcop + ca - 1 / meso_scale - 2 * Ci0_m + 2 * gam) 
    + Ci0_m * (-A0_m / gcop - ca + 1 / meso_scale);
  real q_c = A0_m * (gam * (2 * ca - 2 / meso_scale - Ci0_m - 2 * gam) + Ci0_m 
    * (ca - 1 / meso_scale));

  real A = (-q_b - sqrt(q_b ^ 2 - 4 * q_a * q_c)) / (2 * q_a);

  real ci = ca - A / gcop - 1 / meso_scale;
  real D13C = 4.4 + (b - 4.4) * ci / ca;
  real d13C_m = d13Ca_m - D13C;
}

// The model to be estimated.
model {
  d13Ca_m ~ normal(d13Ca[1], d13Ca[2]);
  
  s2_m ~ gamma(s2[1] * s2_beta, s2_beta);
  
  amax_scale ~ beta(s3[1] * amax_v, (1 - s3[1]) * amax_v); // aka s3
  
  gc_scale ~ beta(s4[1] * gc_v, (1 - s4[1]) * gc_v); // aka s4
  
  meso_scale ~ beta(s5[1] * meso_v, (1 - s5[1]) * meso_v); // aka s5
  
  Ci0_m ~ gamma(Ci0[1] * Ci0_beta, Ci0_beta);
  
  A0_m ~ gamma(A0[1] * A0_beta, A0_beta);
  
  gb_m ~ gamma(gb[1] * gb_beta, gb_beta);
  
  l ~ uniform(1e-6, 1e-4);
  
  gcl_m ~ gamma(2, 0.08);

  SA_gc ~ beta(1.5, 20);
  
  s1_m ~ gamma(s1[1] * s1_beta, s1_beta);

  d13Cp[1] ~ normal(d13C_m, d13Cp[2]);
  Dab[1] ~ normal(D, Dab[2]);
  GCLab[1] ~ normal(Pl / s1_m, GCLab[2]);
  GCWab[1] ~ normal(l / s2_m, GCWab[2]);
}

