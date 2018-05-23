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
  float Z=0, Z_hat=0, MSE_vali, MSE_small, ct_o=5, ct, range;


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
//      Infulence function used for U-type robust cov estimation         //
//////////////////////////////////////////////////////////////////////////
////See section 2.4.1 for more details

//Infulence Matrix
// [[Rcpp::export]]
arma::mat Influence_Huber (arma::mat X, float tau)
{
  using namespace arma;
  int i;       float v1;
  int d; d=X.n_rows;

  vec eigval_cov;  eigval_cov.zeros(d);
  vec eigval_sign; eigval_sign.zeros(d);
  vec eigval_phi;  eigval_phi.zeros(d);
  mat eigvec_cov;  eigvec_cov.zeros(d,d);
  mat X_return;    X_return.zeros(d,d);

  eig_sym(eigval_cov, eigvec_cov, X);

  eigval_phi=eigval_cov;

  eigval_sign=sign(eigval_cov);

  for(i=0; i<d; i++){
    v1=as_scalar(eigval_cov[i]);
    if(std::abs(v1)>tau)eigval_phi[i]=eigval_sign[i]*tau;

  }

  mat eig_lambda; eig_lambda.zeros(d,d);
  eig_lambda.diag()=eigval_phi;


  X_return=eigvec_cov * eig_lambda * eigvec_cov.t();


  return X_return;

}










///////////////////////////////////////////////////////////////////////////
//            Robust estimate of mu                                     //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust(arma::mat X)
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
//            Robust estimate of mu: no cross validation                               //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust_noCV(arma::mat X, arma::mat tau)
{
  using namespace arma;
  int i, P, N;
  //Initial value of Huber descent
  P=X.n_rows;    N=X.n_cols;
  float tau_select;

  float Z=0.5;


  mat Xi;
  mat mu_hat; mu_hat.zeros(P);
  mat mu_one; mu_one.ones(N);


  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    Xi=X.row(i);
    tau_select  = tau(i);
    mu_hat(i)=Huber_descent (Xi, mu_one, Z,tau_select);

  }


  return mu_hat;

}
///////////////////////////////////////////////////////////////////////////
//            Robust estimate of mu and factor coefficients              //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust_F(arma::mat X, arma::mat phi)
{
  using namespace arma;
  int i, P, K;
  //Initial value of Huber descent
  P=X.n_rows;
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
//            Robust estimate of mu and factor coefficients: no CV        //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X
//Output: Estimated mean mu_hat
// [[Rcpp::export]]
arma::mat mu_robust_F_noCV(arma::mat X, arma::mat phi,arma::mat tau)
{
  using namespace arma;
  int i, P, K;
  //Initial value of Huber descent
  P=X.n_rows;
  K=phi.n_cols;
  //The order of Tau see Theorem 2.7
  float tau_select;
  mat F_H_0; F_H_0.ones(K);

  mat Xi;
  mat mu_hat; mu_hat.zeros(K,P);

  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    Xi=X.row(i);
    F_H_0=solve(phi,trans(Xi));
    tau_select  = tau(i);
    mu_hat.col(i)=Huber_descent_F(Xi, phi, F_H_0, tau_select);
  }
  return mu_hat;

}
///////////////////////////////////////////////////////////////////////////
//       Entry-wise Huber robust eatimaiton of covariance matrix        //
//////////////////////////////////////////////////////////////////////////
//Input: Data matrix X, estimated mean mu_hat
//Output: Estimated cov matrix Sigma_hat
// [[Rcpp::export]]
arma::mat Cov_Huber(arma::mat X, arma::mat mu_hat)
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


//Output: Estimated tuning parameter for covriance matrix
// [[Rcpp::export]]
arma::mat Cov_Huber_tune( arma::mat X, float tau)
{
  using namespace arma;
  int i, j, P=X.n_rows, N=X.n_cols;

  //Tuning parameter
  float tautemp;

  //Define the matrices
  mat Xi, Xj;
  mat CT; CT.zeros(P,P);
  mat Z; Z.ones(N);
  mat sd;

  // Entry-wise Huber method
  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    for(j=0;j<=i;j++){
      Xi=X.row(i); Xj=X.row(j);
      Z = Xi%Xj;
      sd = arma::stddev(Z,1,1);
      tautemp  = N/log(static_cast<double>((P^2)*N));
      CT(i,j) = tau*sd(0,0)*sqrt(static_cast<double>(tautemp));
      //cout << CT(i,j);
      CT(j,i)=CT(i,j);
    }
  }

  return CT;

}

///////////////////////////////////////////////////////////////////////////
//       Entry-wise Huber robust eatimaiton of covariance matrix   : noCV     //
//////////////////////////////////////////////////////////////////////////
//Input: Data matrix X, estimated mean mu_hat
//Output: Estimated cov matrix Sigma_hat
// [[Rcpp::export]]
arma::mat Cov_Huber_noCV(arma::mat X, arma::mat mu_hat, arma::mat tau)
{
  using namespace arma;
  int i, j, P=X.n_rows;

  //Initial value of Huber descent
  float Z=0.5;

  //Tuning parameter
  float tau_select;

  //Define the matrices
  mat Xi, Xj;
  mat Sigma_hat; Sigma_hat.zeros(P,P);


  //Entry-wise Huber method
  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    for(j=0;j<=i;j++){
      Xi=X.row(i); Xj=X.row(j);
      tau_select = tau(i,j);
      //printf("\n (%d, %d) th CT= %f", i,j,CT);
      Sigma_hat(i,j)=Huber_descent(Xi, Xj, Z, tau_select);
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
//             Estimate the loadings via Sample COV and PCA              //
//////////////////////////////////////////////////////////////////////////

//Input:  Covariance matrix M and the number of factors K
//Output: Eigen-vectors corresponding to top K eigenvalues
// [[Rcpp::export]]
arma::mat Loading_Sample(int K, arma::mat M)
{
  using namespace arma;
  int P=M.n_rows;

  //Define matrices for eigenvalues and eigen-vectors
  vec eigval_cov;  eigval_cov.zeros(P);
  mat eigvec_cov;  eigvec_cov.zeros(P,P);
  mat Lambda_hat;  Lambda_hat.zeros(P,K);

  eig_sym(eigval_cov, eigvec_cov, M);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);

  Lambda_hat=eigvec_cov.cols(0,K-1)*sqrt(static_cast<double>(P));



  return Lambda_hat;

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




///////////////////////////////////////////////////////////////////////////
//            U-type robust eatimaiton of covariance matrix             //
//////////////////////////////////////////////////////////////////////////

//Input: Data matrix X and tuning parameter Tau
//Output: Estimated cov matrix Sigma_U
// [[Rcpp::export]]
arma::mat Cov_U(float C_tau, arma::mat X)
{
  using namespace arma;
  int i, j, P=X.n_rows, N=X.n_cols;
  float v1=0, v2=0, v3=0;
  float Tau=C_tau*P*sqrt(N/log(N));
  //Define the matrices
  mat A, B;
  mat Sigma_U; Sigma_U.zeros(P,P);


  //U-type COV estimate method
  for(i=1;i<N;i++){
    for(j=0;j<i;j++){
      A=X.col(j)-X.col(i);
      B=A*A.t();
      v1=as_scalar(A.t()*A);
      v2=v1/2;
      if(v2>Tau)v2=Tau;
      v3=v2/v1/2;
      Sigma_U+=v3*B;
    }
  }

  Sigma_U=Sigma_U/N/(N-1);


  return Sigma_U;

}

