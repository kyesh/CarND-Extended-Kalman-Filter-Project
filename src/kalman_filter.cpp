#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
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
    * predict the state
  */
	cout << "Line1" << endl;
	x_ = F_ * x_;
	cout << "Line2" << endl;
	MatrixXd Ft = F_.transpose();
	cout << "Line3" << endl;
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
        VectorXd z_pred = ToPolar(x_);
        VectorXd y = z - z_pred;
        MatrixXd Ht = H_.transpose();
        MatrixXd S = H_ * P_ * Ht + R_;
        MatrixXd Si = S.inverse();
        MatrixXd PHt = P_ * Ht;
        MatrixXd K = PHt * Si;

        //new estimate
        x_ = x_ + (K * y);
        long x_size = x_.size();
        MatrixXd I = MatrixXd::Identity(x_size, x_size);
        P_ = (I - K * H_) * P_;

}

VectorXd KalmanFilter::ToPolar(const VectorXd& x_state){

    VectorXd h_x(3,1);
    h_x << 1,0,0;
    
    
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    float c1 = sqrt(px*px+py*py);
    
    if ((fabs(px) < 0.0001)||(fabs(c1) < 0.0001)){
        cout << "ToPolar Error - Division by zero" << endl;
        return h_x;
    }
    else {
        h_x << c1, atan2(py, px), (vx*px+vy*py)/c1;
    }
    
    return h_x;
}