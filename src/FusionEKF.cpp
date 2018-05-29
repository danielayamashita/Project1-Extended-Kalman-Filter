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
  //set the acceleration noise components
  
  noise_ax = 9;
  noise_ay = 9;
  //the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);

  //the initial transition matrix P_
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			       0, 1, 0, 0,
			       0, 0, 1000, 0,
			       0, 0, 0, 1000;
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


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
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,ekf_.H_,R_radar_,ekf_.Q_);
      ekf_.x_(0) = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_(1) = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);

      ekf_.x_(2) =  measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_(3) =  measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);

      if(ekf_.x_(0) <0.00001)
      {
        ekf_.x_(0) = 0.00001;
      }
      if(ekf_.x_(1) <0.00001)
      {
        ekf_.x_(1) = 0.00001;
      }

      #ifdef PRINT_HEADERS_H_
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "INIT RADAR Measurement:" << std::endl;
      #endif
      
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
     //ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,ekf_.H_,R_laser_,ekf_.Q_);
     ekf_.x_<< measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
     #ifdef PRINT_HEADERS_H_
     std::cout << "--------------------------------------" << std::endl;
     std::cout << "INIT LASER Measurement:"<< std::endl;
     #endif
    }

    previous_timestamp_ = measurement_pack.timestamp_ ;
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
  //compute the time elapsed between the current and previous measurements
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  #ifdef PRINT_H_
  std::cout << "dt:" <<  dt << std::endl;
  #endif
	ekf_.F_ << 1, 0, dt, 0,
			       0, 1, 0, dt,
			       0, 0, 1, 0,
			       0, 0, 0, 1;
  //ekf_.F_(0,2) = dt;
  //ekf_.F_(1,3) = dt;

  float dt_4;
	float dt_3;
	float dt_2;
	VectorXd z;
	z=VectorXd(2);
  //std::cout << "--------------------------------------" << std::endl;
  //std::cout << "Raw Meas:"<< std::endl;
  //std::cout << measurement_pack.raw_measurements_[0]<<measurement_pack.raw_measurements_[1] <<  z << std::endl;
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];

  //std::cout << "Z:" <<  z << std::endl;

  ekf_.Q_ = MatrixXd(4, 4);
	dt_4 = dt*dt*dt*dt/4; 
	dt_3 = dt*dt*dt/2;
	dt_2 = dt*dt;
  ekf_.Q_   << dt_4*noise_ax    , 0                 ,dt_3*noise_ax      , 0,
		            0               , dt_4*noise_ay     ,0                  , dt_3*noise_ay,
		            dt_3*noise_ax   , 0                 ,dt_2*noise_ax      , 0,
		            0               , dt_3*noise_ay     , 0                 , dt_2*noise_ay;

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
    // Radar update
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    //ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,ekf_.H_,R_radar_,ekf_.Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    //ekf_.Init(ekf_.x_,ekf_.P_,ekf_.F_,ekf_.H_,R_laser_,ekf_.Q_);
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  #ifdef PRINT_OUTPUTS_H_
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  #endif
}
