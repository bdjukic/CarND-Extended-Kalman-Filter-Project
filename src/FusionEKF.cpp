#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// State vector
	ekf_.x_ = VectorXd(4);

	// Initializing matrices
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.P_ = MatrixXd(4, 4);

	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	// The initial transition matrix F_
	ekf_.F_ << 1, 0, 1, 0, 
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;

	// State covariance matrix P
	ekf_.P_ << 1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1000, 0,
				0, 0, 0, 1000;

	// Measurement covariance matrix - Laser
	R_laser_ << 0.0225, 0,
				0, 0.0225;

	// Measurement covariance matrix - Radar
	R_radar_ << 0.09, 0, 0,
				0, 0.0009, 0,
				0, 0, 0.09;
				
	H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;

	Hj_ << 1, 1, 0, 0,
			1, 1, 0, 0,
			1, 1, 1, 1; 

	// Set the acceleration noise components
	acceleration_noise_ax = 9;
	acceleration_noise_ay = 9;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	if (!is_initialized_) {
		cout << "EKF: " << endl;
		
		float position_x;
		float position_y;
		float velocity_x;
		float velocity_y;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			// Convert from polar to cartesian coordinates
			position_x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
			position_y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
			
			velocity_x = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
			velocity_y = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			position_x = measurement_pack.raw_measurements_[0];
			position_y = measurement_pack.raw_measurements_[1];
			
			// We don't get velocty from the LIDAR measurement
			velocity_x = 0;
			velocity_y = 0;
		}
		
		ekf_.x_ << position_x, position_y, velocity_x, velocity_y;

		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;

		return;
	}

	// Compute the time elapsed between the current and previous measurements
	float dt_in_seconds = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	float dt_2 = dt_in_seconds * dt_in_seconds;
	float dt_3 = dt_2 * dt_in_seconds;
	float dt_4 = dt_3 * dt_in_seconds;

	// Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt_in_seconds;
	ekf_.F_(1, 3) = dt_in_seconds;

 	// Set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);

	ekf_.Q_ << dt_4 / 4 * acceleration_noise_ax, 0, dt_3 / 2 * acceleration_noise_ax, 0,
				0, dt_4 / 4 * acceleration_noise_ay, 0, dt_3 / 2 * acceleration_noise_ay,
				dt_3 / 2 * acceleration_noise_ax, 0, dt_2 * acceleration_noise_ax, 0,
				0, dt_3/2 * acceleration_noise_ay, 0, dt_2 * acceleration_noise_ay;

	ekf_.Predict();

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		
		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;
		
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} else {
	    ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
