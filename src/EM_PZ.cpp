// EM algorithm for Latent fixed-effect Model (LFM)
// Input: p-values, fixed effects(Z)

#include "EM_PZ.hpp"

void EM_PZ(vec Pvalue, mat Z, uword M, uword L, double& alpha, vec& b, vec& pi1, uword& iter_times,
          uword maxiter, double tol){

  double Lq = 0;
  double Lq_old = Lq;
  mat ZT = Z.t();
  rowvec one = ones<rowvec>(L);

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec Zb = Z*b;
    vec sigma_Zb = 1/(1+exp(-Zb));
    // sigma_Zb.elem(find(sigma_Zb>0.999)).fill(0.999);
    // sigma_Zb.elem(find(sigma_Zb<0.001)).fill(0.001);
    vec comp_pos = exp(Zb) * alpha % pow(Pvalue,(alpha-1));
    pi1 = comp_pos/(comp_pos+1);

    // compute incomplete log likelihood
    Lq = (sum(pi1)*log(alpha) + (alpha-1)*dot(pi1, log(Pvalue)) + dot(pi1, Zb-log(pi1+(pi1==0)))
            - dot(1-pi1, log(1-pi1+(pi1==1))) + sum(log(1-sigma_Zb)));
    // printf("Stage 2 iter=%d\n", iter+1);

    // M step
    alpha = -sum(pi1)/dot(pi1,log(Pvalue));
    vec g = ZT*(sigma_Zb-pi1);
    mat comp_H = sqrt(sigma_Zb%(1-sigma_Zb))*one%Z;
    mat H = comp_H.t()*comp_H;
    b -= solve(H, g);

    // check convergence
    if (iter != 0){
      if (Lq < Lq_old){
        printf("L is not increasing. \n");
        break;
      }
      if ((Lq-Lq_old)/abs(Lq) < tol){
        iter_times = iter+1;
        break;
      }
      else{
        iter_times = maxiter;
      }
    }
    Lq_old = Lq;
  }
}
