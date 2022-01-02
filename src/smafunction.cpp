#include <RcppArmadillo.h>
#include <iostream>
#include <armadillo>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;


///////////////////////general function 
double sum(arma::vec x){
  double sumx = 0;
  int p = x.size();
  for(int i = 0;i < p;i++){
    sumx+=x(i);
  }
  return sumx;
}

arma::mat Deleteone(arma::mat Data, int DeleteNum){
  int N = Data.n_rows;
  arma::mat newData;
  arma::mat newData1;
  arma::mat newData2;
  if(DeleteNum == 0){
    newData = Data.rows(1, N-1);
  }
  if(0 < DeleteNum && DeleteNum < N-1){
    newData1 = Data.rows(0, DeleteNum-1);
    newData2 = Data.rows(DeleteNum+1, N-1);
    newData = join_cols(newData1, newData2);
  }
  if(DeleteNum == N-1){
    newData = Data.rows(0, N-2);
  }
  return newData;
}

arma::vec Deletevec(arma::vec x,int u){
  int p = x.size();
  arma::vec x_new(p-1);
  for(int i =0;i < u ;i++){
    x_new(i) = x(i);
  }
  for(int j = u;j < p-1;j++){
    x_new(j) = x(j+1);
  }
  return x_new;
}

double prasadraoest(int m, int p, arma::vec y, arma::mat X, arma::vec D){
  
  arma::vec coef = inv(X.t() * X) * X.t() * y;
  arma::vec resid = y - X * coef;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sumy = 0.0;
  arma::mat tXX_inv = inv(X.t() * X);
  for(int i=0; i<m; i++){
    sum1 = sum1 + pow(resid(i),2);
    arma::vec temp2 = D(i)*(1 - X.row(i)* tXX_inv * X.row(i).t());
    sum2 = sum2 + temp2(0);
    sumy = sumy + pow(y(i) - mean(y),2);
  }
  double Ahat = 0.0001;
  if(p>1){
    double tempA = (sum1 - sum2)/(m - p);
    if(tempA > 0.0001){
      Ahat = tempA;
    }
  }else{
    double tempA = sumy/(m-1) - sum(D)/m;
    if(tempA > 0.0001){
      Ahat = tempA;
    }
  }
  return Ahat;
}

arma::mat invervariance(double scalar, arma::vec D){
  int m = D.size();
  arma::mat tempD = arma::zeros(m,m);
  for(int i=0; i<m; i++){
    tempD(i,i) = 1.0/(scalar + D(i));
  }
  return tempD;
}

arma::mat Pmat(double scalar, arma::mat X, arma::vec D){
  arma::mat tempp = invervariance(scalar,D) - (invervariance(scalar,D) * X * inv(X.t() * 
    invervariance(scalar,D) * X) * X.t() * invervariance(scalar,D));
  return tempp;
}

double resimaxilikelihood(arma::vec y, arma::mat X, arma::vec D, int maxiter){
  int m = y.size();
  int p = X.n_cols;
  double psihat = prasadraoest(m, p, y, X, D);
  for(int it=0; it <maxiter; it++){
    if(psihat < 0.0001){
      psihat = 0.0001;
    }else{
      double trace1 = 0.0;
      double trace2 = 0.0;
      arma::mat pmatpre = Pmat(psihat,X,D);
      arma::mat pmatpre2 = Pmat(psihat,X,D) * Pmat(psihat,X,D);
      for(int t=0; t<m; t++){
        trace1 = trace1 + pmatpre(t,t);
        trace2 = trace2 + pmatpre2(t,t);
      }
      arma::vec tempbias = (y.t() * pmatpre2 * y - trace1)/trace2;
      double bias = tempbias(0);
      double psi_temp = psihat + bias;
      psihat = psi_temp;
      if(bias < sqrt(m)){
        if(psihat < 0){psihat = 0.0001;}
        break; }
      if(psihat < 0){psihat = 0.0001;}
    }
  }
  return psihat;
}

List empiricalbayes(arma::mat x, arma::vec y, arma::vec D, int m, int p){
  arma::mat tempx = arma::zeros(p,p);
  arma::mat tempxy = arma::zeros(p, 1);
  for(int i=0; i<m; i++){
    tempx = tempx + x.row(i).t() * x.row(i);  
    tempxy = tempxy + y(i) * x.row(i).t();
  }
  
  arma::mat tempx_inv = inv(tempx);
  arma::vec beta_ols = tempx_inv * tempxy;
  double sum1 = 0.0; 
  double sum2 = 0.0; 
  for(int ii=0; ii<m; ii++){
    arma::vec tempsum1 = pow(y(ii) - x.row(ii) * beta_ols, 2);
    sum1 = sum1 + tempsum1(0);
    arma::vec tempsum2 = (1 - (x.row(ii) * tempx_inv * x.row(ii).t())) * D(ii);
    sum2 = sum2 + tempsum2(0);
  }
  double Ahat = (sum1 - sum2)/(m-p);
  if(Ahat < 0){
    Ahat = 0;
  }
  arma::mat tempb1 = arma::zeros(p, p);
  arma::mat tempb2 = arma::zeros(p, 1);;
  for(int j=0; j<m; j++){
    tempb1 = tempb1 + (x.row(j).t() * x.row(j))/(Ahat + D(j));
    tempb2 = tempb2 + x.row(j).t() * y(j)/((Ahat + D(j)));
  }
  arma::vec bhat = inv(tempb1) * tempb2;
  return List::create(Named("bhat") = bhat, Named("Ahat") = Ahat);
}

//////////////////// Mcjack method 
double thetafun(arma::vec beta, double A, arma::mat XX, double xi){
  double theta = dot(beta,XX) + pow(A,0.5) *xi;
  return theta;
}

double thetahatfun(arma::vec beta, double A, double D, arma::mat XX, double YY){
  double temp = dot(beta,XX);
  double thetahat = temp + A * (YY - temp)/(A+D);
  return thetahat;
}

List Jackfun(arma::mat X, arma::vec Y, arma::vec D, int m){
  arma::vec Ajack(m);
  int p = X.n_cols;
  arma::mat bjack = arma::zeros(m,p);
  for(int j = 0;j < m;j++){
    arma::mat X_j = Deleteone(X,j);
    arma::vec Y_j = Deletevec(Y,j);
    arma::vec D_j = Deletevec(D,j);
    Ajack(j) = resimaxilikelihood(Y_j, X_j, D_j,200);
    arma::mat Vj_inv = arma::zeros(m-1,m-1);
    for(int i = 0;i < m - 1;i++){
      Vj_inv(i,i) = 1/(Ajack(j) + D_j(i));
    }
    arma::mat temp_vec = inv(X_j.t() * Vj_inv * X_j) * X_j.t() * Vj_inv * Y_j;
    bjack.row(j) = temp_vec.t();
  }
  return List::create(Named("Aj") = Ajack, Named("bj") = bjack);
}

//[[Rcpp::export]]
arma::vec mspe_FH_McJack(int m, arma::mat X,arma::vec Y, arma::vec D, arma::vec bhat, double Ahat, const int K, int maxiter){
  
  List temp = Jackfun(X,Y,D,m);
  arma::vec Aj = temp["Aj"];
  arma::mat bhatj = temp["bj"];
  arma::vec bpsi(m);
  for(int t = 0;t < m;t++){
    arma::vec dif(K);
    arma::mat difmat(m,K);
    for(int k = 0;k < K;k++){
      arma::mat xi_mat = arma::randn(1,1);
      double xi = xi_mat(0,0);
      double theta_ik = thetafun(bhat,Ahat,X.row(t),xi);
      arma::mat temp_mat = arma::randn(1,1);
      double Y_ik = theta_ik + pow(D(t),0.5) * temp_mat(0,0);
      arma::vec Y_k = Y;
      Y_k(t) = Y_ik;
      double Ahatk = resimaxilikelihood(Y_k,X,D,maxiter);
      arma::mat V_invk = arma::zeros(m,m);
      for(int i = 0;i < m;i++){
        V_invk(i,i) = 1/(Ahatk + D(i));
      }
      arma::vec bhatk = inv(X.t() * V_invk * X) *  X.t() * V_invk * Y_k;
      
      double thetahat_ik = thetahatfun(bhatk ,Ahatk ,D(t) ,X.row(t) ,Y_ik);
      dif(k) = pow(thetahat_ik - theta_ik,2);
      List tempk = Jackfun(X, Y_k, D,m);
      arma::vec Ajk = tempk["Aj"];
      arma::mat bhatjk = tempk["bj"];
      for(int ii = 0;ii < m;ii++){
        double theta_kj = thetafun(bhatj.row(ii).t() ,Aj(ii) ,X.row(t) ,xi);
        double thetahat_kj = thetahatfun(bhatjk.row(ii).t() ,Ajk(ii) ,D(t) ,X.row(t), Y_ik);
        difmat(ii,k) = pow(thetahat_kj - theta_kj, 2);
      }
    }
    arma::vec meandifmat = mean(difmat,1);
    arma::vec logvec(m);
    for(int i=0; i<m; i++){
      logvec(i) = log(meandifmat(i));
    }
    bpsi(t) = m * log(sum(dif)/K) - (m - 1) * sum(logvec)/ m ;

  }
  arma::vec bpsi_inv = exp(bpsi);
  return bpsi_inv;
}

///////////////////// Jackknife method
List mseu(arma::vec D,int index, int m, arma::mat X, arma::mat Y, int maxiter){
  arma::vec bi_u(m);
  arma::vec thetaHat_u(m);
  for(int u = 0;u < m;u++){
    arma::vec D_new = Deletevec(D, u);
    arma::mat dataX = Deleteone(X, u);
    arma::mat dataY = Deleteone(Y, u);
    arma::vec D_u = Deletevec(D,u);
    
    double A_REML_u = resimaxilikelihood(dataY,dataX,D_new,maxiter);
    arma::mat v_inv = arma::zeros(m-1,m-1);
    for(int i=0; i<(m-1); i++){
      v_inv(i,i) = 1.0/(A_REML_u+D_u(i));
    }
    arma::vec bhat_u = inv(dataX.t() * v_inv * dataX) * dataX.t() * v_inv * dataY;
    double B_u = D(index) / (A_REML_u+D(index));
    bi_u(u) = A_REML_u * B_u;
    arma::mat temp = X.row(index) * bhat_u + (1 - B_u) * (Y(index) - X.row(index) * bhat_u);
    thetaHat_u(u) = temp(0,0);
  }
  return List::create(Named("bhatu") = bi_u, Named("thetahatu") = thetaHat_u);
}

// [[Rcpp::export]]
arma::vec mspe_FH_Jackknife(int p, int m, arma::mat X,arma::mat Y,arma::vec D, arma::vec bhat, double A_REML, int maxiter){
  arma::vec mseFH(m);
  arma::vec B = D / (D+A_REML);
  arma::vec b = A_REML * B;
  arma::vec thetahat(m);
  arma::vec temp = X * bhat;
  
  for(int i = 0;i<m;i++){
    thetahat(i) = temp(i) +(1 - B(i))*(Y(i) - temp(i));
  }
  for(int i = 0;i < m;i++){
    List uhat = mseu(D,i,m,X,Y,maxiter);
    arma::vec bhatu =uhat["bhatu"];
    arma::vec thetahatu = uhat["thetahatu"];
    mseFH(i) = b(i) - (m - 1) * sum(bhatu - b(i))/m + (m - 1) * sum(pow(thetahatu - thetahat(i),2))/m;
  }
  return mseFH;
}

//////////////////// Parameter Bootstrap method
List FH_bootstrap(int m, int p, arma::mat X, arma::vec Y_star, arma::vec D, arma::mat Y){
  List Est_FH = empiricalbayes(X, Y_star, D, m, p);
  double psi_FH = Est_FH["Ahat"];
  arma::vec bhat = Est_FH["bhat"];
  arma::vec g1 = arma::vec(m);
  arma::vec g2 = arma::vec(m);
  arma::vec g3 = arma::vec(m);
  arma::vec theta = arma::vec(m);
  arma::mat inv_V = arma::zeros(m,m);
  for(int v=0; v<m; v++){
    inv_V(v,v) = 1/(psi_FH + D(v));
  }
  for(int i=0; i<m; i++){
    g1(i) = psi_FH * D(i) / (psi_FH + D(i));
    arma::vec tempg2 = ((pow(D(i), 2) / pow(psi_FH + D(i), 2)) * X.row(i) * 
      inv(X.t() * inv_V * X) * X.row(i).t());
    g2(i) = tempg2(0);
    arma::vec temptheta = X.row(i) * bhat + (Y(i) - X.row(i) * bhat) * psi_FH / (psi_FH + D(i));
    theta(i) = temptheta(0);
  }
  return List::create(Named("g1") = g1, Named("g2") = g2, Named("theta") = theta,Named("bhat") = bhat, Named("psi_FH") = psi_FH);
}

//[[Rcpp::export]]
arma::vec mspe_FH_PB(arma::mat X, arma::vec Y, arma::vec D, int K){
  int m = Y.size();
  int p = X.n_cols;
  
  List temp2 = FH_bootstrap(m, p, X, Y, D, Y);
  arma::vec g1 = temp2["g1"];
  arma::vec g2 = temp2["g2"];
  arma::vec theta = temp2["theta"];
  arma::vec bhat = temp2["bhat"];
  double psi_FH = temp2["psi_FH"];
  arma::mat g1_star = arma::zeros(m, K);
  arma::mat g2_star = arma::zeros(m, K);
  arma::mat theta_star = arma::zeros(m, K);
  
  for(int k=0; k<K; k++){
    arma::vec eistar(m);
    for(int e=0; e<m; e++){
      arma::vec tempe = pow(D(e), 0.5)* arma::randn(1,1);
      eistar(e) = tempe(0);
    }
    arma::vec Y_star = X*bhat + pow(psi_FH, 0.5) * arma::randn(m,1) + eistar; 
    List temp1 = FH_bootstrap(m, p, X, Y_star, D, Y);
    arma::vec g1_stark = temp1["g1"];
    g1_star.col(k) = g1_stark;
    arma::vec g2_stark = temp1["g2"];
    g2_star.col(k) = g2_stark;
    arma::vec thetak_star = temp1["theta"];
    theta_star.col(k) = thetak_star;
  }
  
  arma::mat summat = arma::zeros(m,K);
  for(int i=0; i<K; i++){
    for(int j=0; j<m; j++){
      summat(j,i) = pow(theta_star(j, i) - theta(j), 2);
    }
  }
  arma::mat tempmat = g1_star + g2_star;
  arma::vec tempmean1(m);
  arma::vec tempmean2(m);
  for(int i=0; i<m; i++){
    tempmean1(i) = mean(tempmat.row(i));
    tempmean2(i) = mean(summat.row(i));
  }
  
  arma::vec mspe_theta = 2.0 * (g1 + g2) - tempmean1 + tempmean2;
  return mspe_theta;
}

/////////////////// Double parameter Bootstrap method
List FH_DB_REML(int m, int p, arma::mat X, arma::vec Y, arma::vec D, int maxiter){
  double Ahat = resimaxilikelihood(Y, X, D, maxiter);
  arma::mat V_inv = arma::zeros(m,m);
  for(int i=0; i<m; i++){
    V_inv(i,i) = 1/(Ahat + D(i));
  }
  arma::vec bhat = arma::inv(X.t() * V_inv * X) * X.t() * V_inv * Y;
  arma::vec theta(m);
  for(int j=0; j<m; j++){
    arma::vec temp_theta = X.row(j) * bhat + (Y(j) - X.row(j) * bhat) * (Ahat/(Ahat + D(j)));
    theta(j) = temp_theta(0);
  }
  return List::create(Named("theta") = theta, Named("bhat") = bhat, Named("Ahat") = Ahat);
}

//[[Rcpp::export]]
arma::vec mspe_FH_DB(arma::mat X, arma::vec y, arma::vec D,int K, int C, arma::vec b, double Ahat, int maxiter){
  int m = y.size();
  int p = X.n_cols;
  arma::mat umat = arma::zeros(m, K);
  arma::mat vmat = arma::zeros(m, K);
  arma::mat temp_mat = arma::zeros(m, K);
  for(int k=0; k<K; k++){
    arma::vec thetak = X * b +  pow(Ahat, 0.5) * arma::randn(m,1);
    arma::vec eistar1(m);
    for(int d=0; d<m; d++){
      arma::vec tempe1 = pow(D(d), 0.5)* arma::randn(1,1);
      eistar1(d) = tempe1(0);
    }
    arma::vec Y_star1 = thetak + eistar1;
    List temp1 = FH_DB_REML(m, p, X, Y_star1, D, maxiter);
    arma::vec bhatk = temp1["bhat"];
    double Ahatk = temp1["Ahat"];
    arma::vec thetahatk = temp1["theta"];
    for(int u=0; u<m; u++){
      umat(u , k) = pow(thetahatk(u) - thetak(u), 2);
    }
    for(int c=0; c<C; c++){
      arma::vec thetac = X * bhatk +  pow(Ahatk, 0.5) * arma::randn(m,1);
      arma::vec eistar2(m);
      for(int d=0; d<m; d++){
        arma::vec tempe2 = pow(D(d), 0.5)* arma::randn(1,1);
        eistar2(d) = tempe2(0);
      }
      arma::vec Y_star2 = thetac + eistar2;
      
      List temp2 = FH_DB_REML(m, p, X, Y_star2, D, maxiter);
      arma::vec thetahatkc = temp2["theta"];
      for(int v=0; v<m; v++){
        temp_mat(v,c) = pow(thetahatkc(v) - thetac(v), 2);
      }
    }
    for(int tt=0; tt<m; tt++){
      vmat(tt, k)  = mean(temp_mat.row(tt));
    }
  }
  arma::vec u(m);
  arma::vec v(m);
  for(int uu=0; uu<m; uu++){
    u(uu) = mean(umat.row(uu));
    v(uu) = mean(vmat.row(uu));
  }
  arma::vec mspe_DB(m);
  for(int ii=0; ii<m; ii++){
    if(u(ii) >= v(ii)){
      mspe_DB(ii) = u(ii) + atan(m * (u(ii) - v(ii)))/m;
    }else{
      mspe_DB(ii) = pow(u(ii),2)/(u(ii) + atan(m * (v(ii) - u(ii)))/m);
    }
  }  
  return mspe_DB;
}

////////////////// Sucma method 
arma::vec AhatK(int m ,arma::mat X, arma::vec bhat, double A_REML, arma::vec D, int K,int p){
  arma::vec A_REML_K(K); 
  for(int i = 0;i < K;i++){
    arma::vec Y_K = X * bhat + sqrt(A_REML) * arma::randn(m,1) + sqrt(D) % arma::randn(m,1);
    A_REML_K(i) = resimaxilikelihood(Y_K,X,D,1000);
  }
  return A_REML_K;
}

// [[Rcpp::export]]
arma::vec mspe_FH_Sucma(int p, int m, arma::mat X,arma::vec Y,arma::vec D, arma::vec bhat, double A_REML, int K){
  arma::vec a_yAhat = A_REML * D/(A_REML + D);
  arma::vec a_ykAhat = a_yAhat;
  arma::vec A_REML_K = AhatK(m, X ,bhat ,A_REML, D, K, p );
  arma::mat a_ykAhatk = arma::mat(m, K);
  for(int j = 0;j < K;j++){
    a_ykAhatk.col(j) = A_REML_K(j) * D/(A_REML_K(j)+D);
  }
  arma::vec mseFH = 2 * a_yAhat - mean(a_ykAhatk,1);
  return mseFH;
}
