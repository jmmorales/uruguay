functions { 
    
  /* compute the kronecker product
  * Args: 
  *   A,B: matrices 
  * Returns: 
  *   kronecker product of A and B
  */ 
  matrix kronecker(matrix A, matrix B) { 
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
    for (i in 1:cols(A)) { 
      for (j in 1:rows(A)) { 
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
      } 
    } 
    return kron; 
  } 
   
   
  // function for covariance based on Gaussian prior (from McElreath)
  // removed self variance as we assume no change between visits
   
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  } 

} 
    
data { 
  int<lower=1> N;                 // total number of observations 
  int<lower=1> T;                 // Number of (potential) replicates at each site
  int<lower=1> E;                 // Number of establecimientos
  array[N,T] int<lower=0,upper=1> Y;    // response variable 
  array[N] int<lower=1, upper=T> visits;
  int<lower=1> K;                 // number of sample-level predictors (potrero + establecimiento) parapoder construir una Z ?nica
  int<lower=1> N_J;               // num of groups (bird spp)
  int<lower=1> L_J;               // num group level predictors (traits)
  int<lower=1> N_1;               // num sites (transectas)
  array[N] int<lower=1,upper=N_J> J;    // group id 
  array[N]int<lower=1,upper=N_1> J_1;  // site id
  int<lower=1> Kx;
  matrix[N, Kx] X;                // obs-level design matrix 
  matrix[N_J, L_J] TT;       // group-level traits presence
  matrix[N_J, N_J] C;             // phylogenetic correlation matrix
  vector[N_J] ones;               // vector on 1s
  int<lower=1> Kp;                // va a ser 1 que es la probabilidad de detecci?n
  int<lower=1> L_Jp;              // num group level predictors (traits) for detection probability
  matrix[N_J, L_Jp] TTp;          // group-level traits for detection probability
  array[N] int<lower=1, upper=E> est;   // identifica el establecimiento de la observacion
  int<lower=1> K_e;               // number of sample-level predictors para establecimiento
  matrix[E, K_e] XE;              // establecimiento-level design matrix 
  matrix[N_1,N_1] Dmat;           // sites distance matrix
  int<lower=1> NY;                 // Number of years
  array[N] int<lower=1, upper=NY> years;   // identifica el a?o de la observacion
}

parameters {
  corr_matrix[K] Omega;     // correlation matrix for var-covar of betas (potrero + establecimiento)
  vector<lower=0>[K] tau;   // scales for the variance covariance of betas (potrero + establecimiento)
  vector[L_J * K] z;        // coeffs for traits
  vector[N_J * K] betas;
  real<lower=0,upper=1> rho;  // correlation between phylogeny and betas
  real<lower=0> taup;         // scales for the variance covariance of betas
  vector[N_J*Kp] ps;          //vector con las p.d de las especies
  vector[L_Jp * Kp] zp;       // coeffs for traits for detection probability
  real<lower=0,upper=1> rhop; // correlation between phylogeny and betas
  real<lower=0> sigmae;       // varianza de la normal de la jerarqu?a establecimiento
  matrix[N_J, E] beta0;       // Intercepto 1 (por establecimiento y especie)
  vector[NY] beta00;     // Intercepto 2 (por a?o)
  vector[N_1] r_1_1;          // site level random effects   
  real<lower=0> etasq;        // max covariance between sites
  real<lower=0> rhosq;        // covariance decay rate
//  real<lower=0> delta;        // variance for sites
  real<lower=0> sigmaee;       // varianza de la normal de la jerarqu?a a?o
}
    
transformed parameters { 
  matrix[K, K] Sigma = quad_form_diag(Omega, tau);
  matrix[N_J*K, N_J*K] S = kronecker(Sigma, rho * C + (1-rho) * diag_matrix(ones));
  matrix[L_J, K] Z = to_matrix(z, L_J, K); // potrero+estb   
  vector[N_J * K] m = to_vector(TT * Z);          // mean of coeffs potrero+estb  
  matrix[N_J, K] b_m = to_matrix(betas, N_J, K);  // coeffs potrero+estb  
  matrix[N_J, N_J] Sp = taup *(rhop * C + (1-rhop) * diag_matrix(ones));
  matrix[L_Jp, Kp] Zp = to_matrix(zp, L_Jp, Kp);  // mean probab. of detection   
  vector[N_J * Kp] mp = to_vector(TTp * Zp);      // mean of coeffs
  matrix[N_J, Kp] p_m = to_matrix(ps, N_J, Kp);   // coeffs
} 
    
model {
  matrix[N_1,N_1] SIGMA;
  // priors
  rhosq ~ exponential( 0.5 );
  etasq ~ exponential( 2 );
//  delta ~ student_t(3,0,10);
  SIGMA = cov_GPL2(Dmat, etasq, rhosq);
  r_1_1 ~ multi_normal( rep_vector(0,N_1) , SIGMA );

  Omega ~ lkj_corr(2);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
  betas ~ multi_normal(m, S);
  //rho ~ beta(2,2);
  z ~ normal(0,1);
  taup ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()
  ps ~ multi_normal(mp, Sp);
  zp ~ normal(0,1);
  //rhop ~ beta(2,2);
  sigmae ~ student_t(3,0,10); // cauchy(0, 2.5); // lognormal()//varianza de la normal para establecimiento
  sigmaee ~ student_t(3,0,10); 
  
  target += log_sum_exp(log(0.5) +  beta_lpdf(rho|1, 100), log(0.5) +  beta_lpdf(rho|2, 2));
  target += log_sum_exp(log(0.5) +  beta_lpdf(rhop|1, 100), log(0.5) +  beta_lpdf(rhop|2, 2));
  
  for (ee in 1:E){
    for(s in 1:N_J){
      beta0[s,ee] ~ normal(XE[ee,1]*b_m[s,5]+XE[ee,2]*b_m[s,6], sigmae);// sampleo para cada especie una probabilidad basal por establecimiento
    }
  }
  
  for(ii in 1:NY)
  {
      beta00[ii] ~ normal(0, sigmaee);
  }
    
  {
    vector[N] mu;
    for (n in 1:N){
      mu[n] = beta0[J[n], est[n]]+ beta00[years[n]] +  X[n,1] * b_m[J[n],1] + X[n,2] * b_m[J[n],2]+ X[n,3] * b_m[J[n],3] + X[n,4] * b_m[J[n],4]+ r_1_1[J_1[n]];
      if (sum(Y[n, 1:visits[n]]) > 0)
      target += log(inv_logit(mu[n])) + bernoulli_lpmf(Y[n,1:visits[n]]  | inv_logit(p_m[J[n],1]));
      else
      target += log_sum_exp(log(inv_logit(mu[n])) + bernoulli_lpmf(Y[n,1:visits[n]]  | inv_logit(p_m[J[n],1])), log(1 - inv_logit(mu[n])));
    }
  }
}

generated quantities{
      vector[N] mus;
      array[N] int pa;
      
    for (n in 1:N){
      mus[n] = beta0[J[n], est[n]] + beta00[years[n]] 
      +  X[n,1] * b_m[J[n],1] + X[n,2] * b_m[J[n],2]+ X[n,3] * b_m[J[n],3] + X[n,4] * b_m[J[n],4]
      + r_1_1[J_1[n]];
      pa[n] = bernoulli_logit_rng(mus[n]);

    }
}
