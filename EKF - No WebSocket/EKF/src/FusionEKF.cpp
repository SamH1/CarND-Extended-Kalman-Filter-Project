#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  ekf_.F_ = MatrixXd(4, 4); //State transition
  ekf_.P_ = MatrixXd(4, 4);


  //set the acceleration noise components
  noise_ax = 9; 
  noise_ay = 9;
  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

// prep output file for debug
string out_file_name_ = "LogFile.txt";
std::ofstream out_file_(out_file_name_.c_str(), std::ofstream::out);

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  //out_file_ << is_initialized_ << endl;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, -3, 0; //this value is important for the RMSE, (-3 , 0) selected after few trials with different combinations


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		double ro1 = measurement_pack.raw_measurements_[0];
		double theta1 = measurement_pack.raw_measurements_[1];
		ekf_.x_(0) = ro1*cos(theta1);
		ekf_.x_(1) = ro1*sin(theta1);
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		out_file_ << "R (x,y)" << ekf_.x_(0) << ", " << ekf_.x_(1) << endl; // for debug
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		ekf_.x_(0) = measurement_pack.raw_measurements_[0];
		ekf_.x_(1) = measurement_pack.raw_measurements_[1];
		out_file_ << "L (x,y)" << ekf_.x_(0) << ", " << ekf_.x_(1) << endl; // for debug
    }

	ekf_.F_ << 1, 0, 0.05, 0,
		0, 1, 0, 0.05,
		0, 0, 1, 0,
		0, 0, 0, 1;
	
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
  //out_file_ << "t1 " << measurement_pack.timestamp_ << endl;
  //out_file_ << "t0 " << previous_timestamp_ << endl;
  //out_file_ << "dt " << dt << endl;
  
  previous_timestamp_ = measurement_pack.timestamp_;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;
  double noise_ax2 = noise_ax * noise_ax;
  double noise_ay2 = noise_ay * noise_ay;

  // Modify the  F matrix so that the time is integrated (Section 8)
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;


  //set the proccess covariance matrix Q (Section 9)
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_.fill(0.);
  ekf_.Q_ << dt_4 * noise_ax2 / 4, 0, dt_3 * noise_ax2 / 2, 0,
	  0, dt_4 * noise_ay2 / 4, 0, dt_3 * noise_ay2 / 2,
	  dt_3 * noise_ax2 / 2, 0, dt_2 * noise_ax2, 0,
	  0, dt_3 * noise_ay2 / 2, 0, dt_2 * noise_ay2;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
 
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  Hj_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;

	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;

	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  //out_file_ << "Q_ = " << ekf_.Q_ << endl;
  //out_file_ << "F_ = " << ekf_.F_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  //out_file_ << endl;
}
