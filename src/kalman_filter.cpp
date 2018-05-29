#include "kalman_filter.h"
#include <math.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;
#include <iostream>

using namespace std;
// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  #ifdef PRINT_H_
  std::cout << "--------------------------------------" << std::endl << std::endl;
  std::cout << "Start Predict" << std::endl << std::endl;
  std::cout << "x_:" <<x_<< std::endl << std::endl;
  std::cout << "F_:" <<F_<< std::endl << std::endl;
  std::cout << "P_:" <<P_<< std::endl << std::endl;
  std::cout << "Q_:" <<Q_<< std::endl << std::endl;
  #endif //PRINT_H_
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_= F_*P_*Ft + Q_;
  #ifdef PRINT_H_
  std::cout << "x_pred:" <<x_<< std::endl << std::endl;
   #endif //PRINT_H_  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  VectorXd y = z - H_ * x_; 
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  #ifdef PRINT_HEADERS_H_
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Update LASER" << std::endl;
  #endif //PRINT_H_
  #ifdef PRINT_H_
    std::cout << "H_:" <<H_<< std::endl<< std::endl;
    std::cout << "x_:" <<x_<< std::endl<< std::endl;
    std::cout << "z:" <<z<< std::endl<< std::endl;
    std::cout << "y:" <<y<< std::endl<< std::endl;
    std::cout << "Ht:" <<Ht<< std::endl<< std::endl;
    std::cout << "S:" <<S<< std::endl<< std::endl;
    std::cout << "K:" <<K<< std::endl<< std::endl;
  #endif //PRINT_H_

  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;

  #ifdef PRINT_H_
    std::cout << "--------------------------------------"<< std::endl;
    std::cout << "x_update:" <<x_<< std::endl << std::endl;
    std::cout << "P_update:" <<x_<< std::endl << std::endl;
  #endif //PRINT_H_
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd h = VectorXd(3);

  


  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);


  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float rho_dot = (px*vx + py*vy) / rho;

  if ((fabs(px)<0.0001) || (fabs(py)<0.0001))
  {
    theta = 0;
    rho_dot = 0;
  }
  h << rho, theta, rho_dot;

  
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  VectorXd y = z - h;
  
  while ( y(1) > M_PI || y(1) < -M_PI ) 
  {
    if ( y(1) > M_PI )
    {
      y(1) -= (2*M_PI);
    } 
    else 
    {
      y(1) += (2*M_PI);
    }
  }
  #ifdef PRINT_HEADERS_H_
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "UpdateEKF RADAR" << std::endl;
    std::cout << "theta:" <<theta<< std::endl<< std::endl;
  #endif //PRINT_H_
  #ifdef PRINT_H_
    
    std::cout << "z:" <<z<< std::endl<< std::endl;
    std::cout << "h(x'):" <<h<< std::endl<< std::endl;
    std::cout << "x_:" <<x_<< std::endl<< std::endl;
    std::cout << "y:" <<y<< std::endl<< std::endl<< std::endl;
    std::cout << " Ht:" << Ht<< std::endl<< std::endl<< std::endl;
    std::cout << " R_:" << R_<< std::endl<< std::endl<< std::endl;
    std::cout << " P_:" << P_<< std::endl<< std::endl<< std::endl;
    std::cout << "K:" <<K<< std::endl<< std::endl;
    std::cout << "S:" <<S<< std::endl<< std::endl;
  #endif //PRINT_H_

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;

  #ifdef PRINT_H
  std::cout << " P_updated:" << P_<< std::endl<< std::endl;
  std::cout << " x_updated:" << x_<< std::endl<< std::endl;
  #endif //PRINT_H_
}
