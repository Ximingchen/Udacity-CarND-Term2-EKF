#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

// the elements of kalman filter are listed above;
void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_; // update x (k+1|k) 
	P_ = F_ * P_ * F_.transpose() + Q_; // update to get P(k+1|k)
}


void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	// here P is represented by P(k|k)
	
	VectorXd y = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si; // calculate Kalman gain

	int dim = x_.size();
	MatrixXd I = MatrixXd::Identity(dim, dim); // define the identity matrix
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

/*
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  
  // the same procedure as above but replacing the H function with Hj
	// first compute the predicted x in z coordinate
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	// z = rho phi rho dot
	float rho = sqrt(px*px + py * py);
	float phi = atan2(py, px);
	float rhodot = (rho > 0.0000001) ? ((px*vy + py * vx) / rho) : rho;

	VectorXd z_predict(3);
	z_predict << rho, phi, rhodot;

	VectorXd y = z - z_predict; // the measurement difference in the new coordinate (dimension 3)
	// consider normalizing the vector y since phi should be within 0 - 2pi
	
	while (y(1) > M_PI || y(1) < -M_PI) {
		if (y(1) > M_PI)
			y(1) -= 2.* M_PI;
		else
			y(1) += 2.* M_PI;
	}
	

	// notice that here we should use Hj instead, this will be handled in FusionEKF.cpp
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_* P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si; // calculate Kalman gain

	int dim = x_.size();
	MatrixXd I = MatrixXd::Identity(dim, dim);
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}
*/
void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// update the state by using Extended Kalman Filter equations

	// get components of predicted state
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	// get components of radar measurement space
	float rho = sqrt(px * px + py * py);
	float phi = atan2(py, px);
	float rho_dot = (rho != 0 ? (px * vx + py * vy) / rho : 0);

	// define predicted position and speed
	VectorXd z_pred(3);
	z_pred << rho, phi, rho_dot;

	// measurement update
	VectorXd y = z - z_pred;

	// normalize angle
	double width = 2 * M_PI;   //
	double offsetValue = y(1) + M_PI;   // value relative to 0
	y(1) = (offsetValue - (floor(offsetValue / width) * width)) - M_PI;

	MatrixXd PHt = P_ * H_.transpose();
	MatrixXd S = H_ * PHt + R_;
	MatrixXd K = PHt * S.inverse();

	// new state
	x_ = x_ + (K * y);
	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}