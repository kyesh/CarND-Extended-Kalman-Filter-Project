#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  //cout << "Start Fusion EKF contructor" <<  endl;
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
	//create a 4D state vector, we don't know yet the values of the x state

	//kf_.x_ = VectorXd(4);

	//state covariance matrix P

	ekf_.P_ = MatrixXd(4, 4);

	ekf_.P_ << 1, 0, 0, 0,
		  0, 1, 0, 0,
		  0, 0, 1000, 0,
		  0, 0, 0, 1000;





	//measurement covariance

//	kf_.R_ = MatrixXd(2, 2);

//	kf_.R_ << 0.0225, 0,

//			  0, 0.0225;



	//measurement matrix

	//kf_.H_laser_ = MatrixXd(2, 4);

	H_laser_ << 1, 0, 0, 0,

			  0, 1, 0, 0;

	//TODO: Hj_



	//the initial transition matrix F_

	//kf_.F_ = MatrixXd(4, 4);

	//kf_.F_ << 1, 0, 1, 0,

	//		  0, 1, 0, 1,

	//		  0, 0, 1, 0,

	//		  0, 0, 0, 1;



	//set the acceleration noise components

	//noise_ax = 5;

	//noise_ay = 5;

	printf("Constructor Complete\n");
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

        float noise_ax = 29;

        float noise_ay = 29;


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: Not Initilized " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
		  0, 1, 0, 1,
		  0, 0, 1, 0,
		  0, 0, 0, 1;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	//cout << "Radar!" << endl;
	float x,y;
	//y = sin(phi)*rho;
	y = sin(measurement_pack.raw_measurements_[1])*measurement_pack.raw_measurements_[0];
	//x = cos(phi)*rho;
	x = cos(measurement_pack.raw_measurements_[1])*measurement_pack.raw_measurements_[0];	
	//cout << "x: " << x << " y: " << y << endl;
	ekf_.x_ << x,y,0,0;
	//cout << "ekf_.x_:" << ekf_.x_ << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	//cout << "Laser!" << endl;
		//set the state with the initial location and zero velocity

		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;






    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
	previous_timestamp_ = measurement_pack.timestamp_;
	//cout << "EKF Initilized" << endl;
    return;
  }
 //cout << "Start Prediction" << endl;
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	//compute the time elapsed between the current and previous measurements

	float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0; // 1000000.0;	//dt - expressed in seconds

	previous_timestamp_ = measurement_pack.timestamp_;



	float dt_2 = dt * dt;

	float dt_3 = dt_2 * dt;

	float dt_4 = dt_3 * dt;



	//Modify the F matrix so that the time is integrated
 //cout << "before setting F" << endl;
	ekf_.F_(0, 2) = dt;

	ekf_.F_(1, 3) = dt;



	//set the process covariance matrix Q
//cout << "before initilizging Q" << endl;
	ekf_.Q_ = MatrixXd(4, 4);
//cout << "before setting Q" << endl;
	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,

			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,

			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,

			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
//cout << "before calle predict" << endl;
  ekf_.Predict();
//cout << "End Prediction" << endl;

//cout << "ekf_.x_:" << ekf_.x_ << endl;

//cout << "Start Update" << endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	//cout << "Radar!" << endl;
	Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;

	//ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
	//cout << "Laser!" << endl;
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;

	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
 //cout << "End Update" << endl;
 //cout << "ekf_.x_:" << ekf_.x_ << endl;
}
