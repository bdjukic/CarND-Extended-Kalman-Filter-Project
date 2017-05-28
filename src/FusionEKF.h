#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
	FusionEKF();
	virtual ~FusionEKF();

	void ProcessMeasurement(const MeasurementPackage &measurement_pack);

	KalmanFilter ekf_;

private:
	bool is_initialized_;
	long long previous_timestamp_;

	Tools tools;
	
	Eigen::MatrixXd R_laser_;
	Eigen::MatrixXd R_radar_;
	Eigen::MatrixXd H_laser_;
	Eigen::MatrixXd Hj_;

	float acceleration_noise_ax;
	float acceleration_noise_ay;
};

#endif /* FusionEKF_H_ */
