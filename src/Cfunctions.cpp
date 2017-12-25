# include <RcppArmadillo.h>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This file contains all the C++ functions used

///////////////////////////////////////////////////////////////////////////
//              Huber estimator for entry-wise Cov Matrix               //
//////////////////////////////////////////////////////////////////////////
////See section 2.4.2 for more details

///////////////  Huber loss function ///////////////////////////
//Calculate Huber Loss l_{tau}(X_{it}'X_{ik} -\sigma)
// [[Rcpp::export]]
arma::mat Huber_loss (arma::mat Vi, arma::mat Vj, float Z, float CT)
{
  using namespace arma;
  int t;       float v1;
  mat Loss;    Loss.zeros(1,1);
  mat Huber;   Huber.zeros(1,1);
  mat M1;      M1.zeros(1,1);
  int T=Vi.n_cols;


  for(t=0; t<T; t++){
    M1=(Vi(t)*Vj(t)-Z);
    v1=fabs(M1(0));

    if(v1> fabs(CT)) {Huber(0)=CT*(2*M1(0)-CT);}
    else {Huber(0)= M1(0)*M1(0);}
    //printf("\n v1=%f  ct=%f Huber=%f", v1, fabs(CT), Huber(0));
    Loss+=Huber(0)/T;
  }

  return Loss;
}


///////////////  Huber gradient function ///////////////////////////
//Calculate the Gradient of Huber Loss
// [[Rcpp::export]]
arma::mat Huber_gradient (arma::mat Vi, arma::mat Vj, float Z, float CT)
{
  using namespace arma;
  int t,T;       float v1;
  mat Grad;        Grad.zeros(1,1);
  mat Huber_dot;   Huber_dot.zeros(1,1);
  mat M1;          M1.zeros(1,1);
  T=Vi.n_cols;


  for(t=0; t<T; t++){
    M1=Vi(t)*Vj(t)-Z;
    v1=M1(0);
    if (v1> fabs(CT))     {Huber_dot(0)=2*CT;}
    else if (v1 < -1*fabs(CT)) {Huber_dot(0)= -2*CT;}
    else {Huber_dot(0)=2*M1(0);}
    //printf("\n x=%f   Huber_dot=%f", M1(0), Huber_dot(0));
    Grad-= 1*Huber_dot(0)/T;
  }

  return Grad;
}



///////////////  Gradient descent of Huber loss  ///////////////////////////
//Minimize Huber loss with gradient descent
// [[Rcpp::export]]
float Huber_descent (arma::mat Vi, arma::mat Vj, float Z, float CT)
{
  using namespace arma;

  int k;    float v1=0, v2=0;
  float Z_1=Z, Z_2=0;
  mat test; test.zeros(1,1);
  for(k=1; k<500; k++){
    v1= as_scalar(Huber_loss (Vi, Vj, Z_1, CT));
    test=Huber_gradient (Vi, Vj, Z_1, CT);
    Z_2=Z_1;      Z_1-=0.5* test(0)/sqrt(static_cast<double>(k));
    v2=as_scalar(Huber_loss (Vi, Vj, Z_1, CT));
    //printf("\n %dth v1=%f    v2=%f   A1=%f    A2=%f    \n", k, v1, v2,fabs(v1-v2));
    if(fabs(v1-v2)<1.0e-8 )k=500;
  }

  return Z_2;

}


///////////////  5-fold CV for tau selection  ///////////////////////////
//Select tuning parameter of Huber loss via 5-fold cross validation
// [[Rcpp::export]]
float Robust_CV (arma::mat Vi, arma::mat Vj)
{
  using namespace arma;
  int i,k,T, T_vali=0;
  float Z=0, Z_hat=0, MSE_vali, MSE_small, ct_o=5, ct, range;;


  //mat vx;        vx=;
  T=Vi.n_cols;
  //Z=as_scalar(mean(vx,1));
  Z=0;
  T_vali=T/5;

  mat Vi_1;      mat Vi_2;
  mat Vi_train;  mat Vi_vali;
  mat Vj_1;      mat Vj_2;
  mat Vj_train;  mat Vj_vali;
  mat validation;

  range=as_scalar(sqrt(T*cov(Vi)/2));


  for(i=5,MSE_small=1.0e8, ct_o=1.0e8; i<=25; i++){

    //ct=(abs(vx).max())*i*0.2;
    ct=range*i*0.1;
    if(i==25)ct=1.0e7;
    //printf("\n---------------  %dth CT=%f ------------------\n",i, ct);

    MSE_vali=0;
    for(k=0; k<5; k++){
      //printf("\n---------------  %dth ------------------\n",k);
      Vi_1.resize(0,0);    Vi_2.resize(0,0);
      Vj_1.resize(0,0);    Vj_2.resize(0,0);
      Vi_vali=Vi.cols(span(k*T_vali, (k+1)*T_vali-1));
      Vj_vali=Vj.cols(span(k*T_vali, (k+1)*T_vali-1));

      if(k > 0){
        Vi_1=Vi.cols(span(0,k*T_vali-1));
        Vj_1=Vj.cols(span(0,k*T_vali-1));
      }

      if(k < 4){
        Vi_2=Vi.cols(span((k+1)*T_vali,T-1));
        Vj_2=Vj.cols(span((k+1)*T_vali,T-1));
      }

      Vi_train=join_rows(Vi_1,Vi_2);
      Vj_train=join_rows(Vj_1,Vj_2);


      Z_hat = Huber_descent (Vi_train, Vj_train, Z, ct);


      MSE_vali+= as_scalar(Huber_loss ( Vi_vali, Vj_vali, Z_hat, ct));


    }
    if(MSE_vali<MSE_small){MSE_small=MSE_vali; ct_o=ct;}


  }

  return ct_o;



}


///////////////////////////////////////////////////////////////////////////
//              Huber estimator for robust factor estimation            //
//////////////////////////////////////////////////////////////////////////
////See section 2.4.2 for more details
//See equaiton (2.11)


///////////////  Huber loss function ///////////////////////////
//Calculate Huber \sum_{j=1}^p Loss l_{tau}(Z_j -\hat{b}_j'f)
// [[Rcpp::export]]
arma::mat Huber_loss_F (arma::mat X, arma::mat phi, arma::mat B, float CT, int T)
{
  using namespace arma;
  int t;       float v1;
  mat Loss;    Loss.zeros(1,1);
  mat Huber;   Huber.zeros(1,1);
  mat M1;      M1.zeros(1,1);

  for(t=0; t<T; t++){
    M1=X(t)-phi.row(t)*B;
    v1=fabs(M1(0));
    if(v1> fabs(CT)) {Huber(0)=CT*(2*M1(0)-CT);}
    else {Huber(0)= M1(0)*M1(0);}
    Loss+=Huber/T;
  }

  return Loss;
}


///////////////  Huber gradient function ///////////////////////////
//Calculate the Gradient of Huber Loss
// [[Rcpp::export]]
arma::mat Huber_gradient_F (arma::mat X, arma::mat phi, arma::mat B, float CT, int T)
{
  using namespace arma;
  int t, J=B.n_rows;       float v1;
  mat Grad;        Grad.zeros(J,1);
  mat Huber_dot;   Huber_dot.zeros(1,1);
  mat M1;          M1.zeros(1,1);
  mat phi_t;       phi_t=phi.t();

  for(t=0; t<T; t++){
    M1=X(t)-phi.row(t)*B;
    v1=M1(0);
    if (v1> fabs(CT))     {Huber_dot(0)=2*CT;}
    else if (v1 < -1*fabs(CT)) {Huber_dot(0)= -2*CT;}
    else {Huber_dot(0)=2*M1(0);}

    Grad-= 1*Huber_dot(0)*phi_t.col(t)/T;
  }

  return Grad;
}

///////////////  Gradient descent of Huber loss ///////////////////////////
//Minimize Huber loss with gradient descent
// [[Rcpp::export]]
arma::mat Huber_descent_F (arma::mat X, arma::mat phi, arma::mat B, float CT)
{
  using namespace arma;

  int k, T=phi.n_rows, J=phi.n_cols;
  float v1=0, v2=0;
  mat b_1;  b_1.zeros(J,1);
  mat b_2;  b_2.zeros(J,1);
  mat test; test.zeros(J,1);
  b_1=B;

  for(k=1; k<500; k++){
    v1= as_scalar(Huber_loss_F (X, phi, b_1, CT, T));
    test=Huber_gradient_F (X, phi, b_1, CT, T);
    b_2=b_1;      b_1-=0.5* test/sqrt(static_cast<double>(k));
    v2=as_scalar(Huber_loss_F (X, phi, b_1, CT, T));

    if(fabs(v1-v2)<1.0e-10 || v1<v2+1.0e-8)k=500;
  }

  return b_2;

}


///////////////  5-fold CV for tau selection  ///////////////////////////
//Select tuning parameter of Huber loss via 5-fold cross validation
// [[Rcpp::export]]
float Robust_CV_F (arma::mat vx, arma::mat phi)
{
  using namespace arma;
  int i,k,T=phi.n_rows, J=phi.n_cols, T_vali=0;
  float MSE_vali, MSE_small, ct_o=5, ct;
  T_vali=T/5;

  mat vx_1;      mat vx_2;
  mat vx_train;  mat vx_vali;    mat vx_hat;
  mat phi_1;     mat phi_2;
  mat phi_vali;  mat phi_train;
  vec b_2;       vec b_sol;
  vec b_1;       b_1.zeros(J,1);

  b_1=solve(phi, vx);

  for(i=1,MSE_small=1.0e8, ct_o=1.0e8; i<=25; i++){

    ct=(abs(vx).max())*i/20;
    if(i==25)ct=1.0e7;


    for(k=0, MSE_vali=0; k<5; k++){
      //printf("\n---------------  %dth ------------------\n",k);
      vx_1.resize(0,0);    vx_2.resize(0,0);
      vx_vali=vx.rows(span(k*T_vali, (k+1)*T_vali-1));
      if(k > 0)vx_1=vx.rows(span(0,k*T_vali-1));
      if(k < 4)vx_2=vx.rows(span((k+1)*T_vali,T-1));
      vx_train=join_cols(vx_1,vx_2);
      phi_1.resize(0,0);    phi_2.resize(0,0);
      phi_vali=phi.rows(span(k*T_vali, (k+1)*T_vali-1));
      if(k > 0)phi_1=phi.rows(span(0,k*T_vali-1));
      if(k < 4)phi_2=phi.rows(span((k+1)*T_vali,T-1));
      phi_train=join_cols(phi_1,phi_2);

      b_2 = Huber_descent_F (vx_train, phi_train, b_1, ct);

      vx_hat=phi_vali*b_2;

      MSE_vali+=as_scalar((vx_hat-vx_vali).t()*(vx_hat-vx_vali)/T_vali);
      //MSE_vali+= as_scalar(Huber_loss_F ( vx_vali, phi_vali, b_2, ct,T_vali));


    }
    if(MSE_vali<MSE_small){MSE_small=MSE_vali; ct_o=ct;}


  }

  return ct_o;

}




///////////////////////////////////////////////////////////////////////////
//            Robust estimate of mu                                     //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X,  constant term of the tuning parameter C_tau
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust(float C_tau, arma::mat X)
{
  using namespace arma;
  int i, P, N;
  //Initial value of Huber descent
  P=X.n_rows;    N=X.n_cols;

  float Z=0.5;
  //The order of Tau see Theorem 2.7
  float Tau;

  mat Xi;
  mat mu_hat; mu_hat.zeros(P);
  mat mu_one; mu_one.ones(N);


  for(i=0;i<P;i++){
      Rcpp::checkUserInterrupt();
    Xi=X.row(i);
    //Un-comment the following line if you want to choose Tau via 5 fold CV
    Tau= Robust_CV ((Xi),trans(mu_one));
    mu_hat(i)=Huber_descent (Xi, mu_one, Z, Tau);

  }


  return mu_hat;

}

///////////////////////////////////////////////////////////////////////////
//            Robust estimate of mu and factor coefficients              //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X,  constant term of the tuning parameter C_tau
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust_F(float C_tau, arma::mat X, arma::mat phi)
{
  using namespace arma;
  int i, P, N, K;
  //Initial value of Huber descent
  P=X.n_rows;    N=X.n_cols;
  K=phi.n_cols;
  //The order of Tau see Theorem 2.7
  float Tau;
  mat F_H_0; F_H_0.ones(K);

  mat Xi;
  mat mu_hat; mu_hat.zeros(K,P);

  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    Xi=X.row(i);
    F_H_0=solve(phi,trans(Xi));
    Tau= Robust_CV_F (trans(Xi),(phi));
   mu_hat.col(i)=Huber_descent_F(Xi, phi, F_H_0, Tau);
  }
 return mu_hat;

}
///////////////////////////////////////////////////////////////////////////
//       Entry-wise Huber robust eatimaiton of covariance matrix        //
//////////////////////////////////////////////////////////////////////////
//Input: Data matrix X, estimated mean mu_hat, constant term of the tuning parameter C_tau
//Output: Estimated cov matrix Sigma_hat
// [[Rcpp::export]]
arma::mat Cov_Huber(float C_tau, arma::mat X, arma::mat mu_hat)
{
  using namespace arma;
  int i, j, P=X.n_rows;

  //Initial value of Huber descent
  float Z=0.5;

  //Tuning parameter
  float Tau;

  //Define the matrices
  mat Xi, Xj;
  mat Sigma_hat; Sigma_hat.zeros(P,P);


  //Entry-wise Huber method
  for(i=0;i<P;i++){
      Rcpp::checkUserInterrupt();
    for(j=0;j<=i;j++){
      Xi=X.row(i); Xj=X.row(j);
      //Un-comment the following line if you want to choose Tau via 5 fold CV
      Tau= Robust_CV (Xi,Xj);
      //printf("\n (%d, %d) th CT= %f", i,j,CT);
      Sigma_hat(i,j)=Huber_descent (Xi, Xj, Z, Tau);
      Sigma_hat(i,j)-=mu_hat(i)*mu_hat(j);
      Sigma_hat(j,i)=Sigma_hat(i,j);
    }
  }


  return Sigma_hat;

}





///////////////////////////////////////////////////////////////////////////
//             Return all the eigenvalues and eigenvectors      //
//////////////////////////////////////////////////////////////////////////

//Input:  Covariance matrix M
//Output: all
// [[Rcpp::export]]
arma::mat Eigen_Decomp( arma::mat M)
{
  using namespace arma;
  int  P=M.n_rows;

  //Define matrices for eigenvalues and eigen-vectors
  vec eigval_cov;  eigval_cov.zeros(P);
  mat eigvec_cov;  eigvec_cov.zeros(P,P);
  mat eigall_cov;  eigall_cov.zeros(P,P+1);

  eig_sym(eigval_cov, eigvec_cov, M);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);
  eigall_cov = join_rows(eigvec_cov, eigval_cov);


  return eigall_cov;

}




///////////////////////////////////////////////////////////////////////////
//             Estimate the loadings via Robust cov and PCA              //
//////////////////////////////////////////////////////////////////////////

//Input:  Covariance matrix M and the number of factors K
//Output: Eigen-vectors corresponding to top K eigenvalues
// [[Rcpp::export]]
arma::mat Loading_Robust(int K, arma::mat M)
{
  using namespace arma;
  int i, P=M.n_rows;

  //Define matrices for eigenvalues and eigen-vectors
  vec eigval_cov;  eigval_cov.zeros(P);
  mat eigvec_cov;  eigvec_cov.zeros(P,P);
  mat Lambda_hat;  Lambda_hat.zeros(P,K);

  eig_sym(eigval_cov, eigvec_cov, M);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);

  for(i=0;i<K;i++)Lambda_hat.col(i)=eigvec_cov.col(i)*sqrt(eigval_cov(i));


  return Lambda_hat;

}




