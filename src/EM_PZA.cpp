// EM algorithm for Latent Sparse Mixed Model (LSMM)
// Input: p-values, fixed effects(Z), random effects(A)

#include "EM_PZA.hpp"

void EM_PZA(vec Pvalue, mat Z, mat A, uword M, uword L, uword K, double& alpha, vec& b, vec& pi1,
            double& sigma2, double& omega, vec& omegak, vec& mu, vec& xi, vec&beta, vec& LQ,
            uword& iter_times, uword maxiter, double tol){

  double Lq = 0;
  double Lq_old = Lq;

  mat ZT = Z.t();
  mat A2 = square(A);
  mat A2T = A2.t();
  rowvec one = ones<rowvec>(L);
  vec u = zeros<vec>(K);
  vec y_tilde = A*(omegak%mu);
  vec sigma_xi = 1/(1+exp(-xi));
  vec lambda_xi = (sigma_xi - 0.5)/xi/2;

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec s2 = sigma2/(1 + 2*sigma2*A2T*lambda_xi);
    vec Zb = Z*b;

    for (uword i = 0; i < K; i++){
      vec y_tilde_i = y_tilde - A.col(i)*omegak[i]*mu[i];
      mu[i] = s2[i] * dot(pi1 - 0.5 - 2*lambda_xi%(Zb + y_tilde_i), A.col(i));
      u[i] = log(omega/(1-omega)) + log(s2[i]/sigma2)/2 + mu[i]*mu[i]/s2[i]/2;
      omegak[i] = 1/(1+exp(-u[i]));
      y_tilde = y_tilde_i + A.col(i)*omegak[i]*mu[i];
    }
    vec v = log(alpha) + (alpha-1)*log(Pvalue) + Zb + y_tilde;
    pi1 = 1/(1+exp(-v));
    xi = sqrt(square(Zb + y_tilde) + A2*(omegak%(s2+mu%mu)-omegak%omegak%mu%mu));

    // compute Lq
    sigma_xi = 1/(1+exp(-xi));
    lambda_xi = (sigma_xi - 0.5)/xi/2;
    Lq = ( (alpha-1)*dot(pi1, log(Pvalue)) + sum(pi1)*log(alpha)
      - dot(pi1, log(pi1+(pi1==0))) - dot(1-pi1, log(1-pi1+(pi1==1)))
      + dot(pi1, Zb+y_tilde) + sum(log(sigma_xi)-(Zb+y_tilde+xi)/2)
      - dot(omegak, s2+mu%mu)*0.5/sigma2 + sum(omegak)*(log(omega)+0.5-log(1-omega)-0.5*log(sigma2))
      + K*log(1-omega) + dot(omegak, 0.5*log(s2)-log(omegak+(omegak==0)))
      - dot(1-omegak, log(1-omegak+(omegak==1)))
    );
    LQ[iter] = Lq;
    // printf("Stage 4 iter=%d\n", iter+1);

    // M step
    alpha = -sum(pi1)/dot(pi1,log(Pvalue));
    vec g = ZT*(-pi1+2*lambda_xi%(Zb+y_tilde)+0.5);
    mat comp_H = sqrt(lambda_xi)*one%Z;
    mat H = 2*comp_H.t()*comp_H;
    b -= solve(H, g);
    sigma2 = dot(omegak, s2+mu%mu)/sum(omegak);
    omega = sum(omegak)/K;

    // check convergence
    if (iter != 0){
      if (Lq < Lq_old){
        printf("L is not increasing. \n");
        break;
      }
      if ((Lq-Lq_old)/abs(Lq) < tol){
        iter_times = iter+1;
        beta = mu%omegak;
        break;
      }
      else{
        iter_times = maxiter;
        beta = mu%omegak;
      }
    }
    Lq_old = Lq;
  }
}
