// EM algorithm for Two Groups Model
// Input: p-values

#include "EM_P.hpp"

void EM_P(vec Pvalue, uword M, double& alpha, double& pi1_, vec& pi1, uword& iter_times,
          uword maxiter, double tol){

  double L = 0;
  double L_old = L;

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec comp_pos = (pi1_*alpha) * pow(Pvalue,(alpha-1));
    pi1 = comp_pos/(comp_pos+1-pi1_);

    // compute incomplete log likelihood
    L = (sum(pi1)*(log(alpha)+log(pi1_)-log(1-pi1_)) + (alpha-1)*dot(pi1, log(Pvalue))
      - dot(pi1, log(pi1+(pi1==0))) + M*log(1-pi1_) - dot(1-pi1, log(1-pi1+(pi1==1))));
    // printf("Stage 1 iter=%d L=%f\n", iter+1, L);

    // M step
    alpha = -sum(pi1)/dot(pi1,log(Pvalue));
    pi1_ = sum(pi1)/M;

    // check convergence
    if (iter != 0){
      if (L < L_old){
        printf("L is not increasing. \n");
        break;
      }
      if ((L-L_old)/abs(L) < tol){
        iter_times = iter+1;
        break;
      }
      else{
        iter_times = maxiter;
      }
    }
    L_old = L;
  }
}
