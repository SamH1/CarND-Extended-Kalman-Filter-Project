#include "kalman_filter.h"
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

// prep output file for debug
string out_file_name_1 = "LogFile1.txt";
std::ofstream out_file_1(out_file_name_1.c_str(), std::ofstream::out);

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

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
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

// prep output file for debug
// string out_file_name_ = "LogFile1.txt";
// std::ofstream out_file_(out_file_name_.c_str(), std::ofstream::out);

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	float x = x_(0);
	float y = x_(1);
	float vx = x_(2);
	float vy = x_(3);


	float rho = sqrt(x*x + y*y);
	float theta = atan2(y,x);
	float ro_dot = (x*vx + y*vy) / rho;
	VectorXd z_pred = VectorXd(3);
	z_pred << rho, theta, ro_dot;

	VectorXd y_ = z - z_pred;

	// normalize angle to be betwwen -Pi and Pi

	while (y_(1) > M_PI) {
		y_(1) -= 2 * M_PI;
	}

	while (y_(1) < -M_PI) {
		y_(1) += 2 * M_PI;
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	out_file_1 << "x_ = " << x_ << endl;
	//new state
	x_ = x_ + (K * y_);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

	// print the output
	out_file_1 << "\nx'_ = " << x_ << endl;
	out_file_1 << "\ny_ = " << y_ << endl;
	out_file_1 << "\nH_ = " << H_ << endl;
	out_file_1 << "\nS_ = " << S << endl;
	out_file_1 << "\nI_ = " << I << endl;
	out_file_1 << "\nK_ = " << K << endl;
	out_file_1 << "-----------------------------" << endl;

}