//Unscented Kalman Filter
#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = .75;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/

	n_x_ = 5;

	n_aug_ = n_x_ + 2;

	lambda_ = 3 - n_x_;

	
}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	if (!is_initialized_) {

		//initializeUKF(meas_package);
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			// x and y positions, known from laser
			x_(0) = meas_package.raw_measurements_(0);
			x_(1) = meas_package.raw_measurements_(1);

			// No data on velocity, yaw and yaw rate.
			x_(2) = 0;
			x_(3) = 0;
			x_(4) = 0;

			
		}
		else
		{
			double radial_distance = meas_package.raw_measurements_(0);
			double angle = meas_package.raw_measurements_(1);
			double radial_speed = meas_package.raw_measurements_(2);

			x_(0) = radial_distance * std::cos(angle);
			x_(1) = radial_distance * std::sin(angle);

			x_(2) = radial_speed;
			x_(3) = angle;
			x_(4) = 0;

			
		}

		// Initialize covariance matrix P
		P_.setIdentity(5, 5);


		VectorXd weights = VectorXd(2 * n_aug_ + 1);

		weights(0) = double(lambda_) / double(lambda_ + n_aug_);

		for (int index = 1; index < 2 * n_aug_ + 1; ++index)
		{
			weights(index) = 0.5 / (lambda_ + n_aug_);
		}
		
		weights_ = weights;
		Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;

		// After initialization return, no need to do prediction and update cycle on initialization

		return;

	}

	if (meas_package.sensor_type_ == MeasurementPackage::LASER) {


		if (use_laser_)
		{
			// Predict new state 
			double time_delta = (meas_package.timestamp_ - time_us_) / 1000000.0;
			time_us_ = meas_package.timestamp_;

			Prediction(time_delta);
			UpdateLidar(meas_package);
		}
	}
	else
	{
		if (use_radar_)
		{

			// Predict new state 
			double time_delta = (meas_package.timestamp_ - time_us_) / 1000000.0;
			time_us_ = meas_package.timestamp_;

			Prediction(time_delta);
			UpdateRadar(meas_package);
		}
	}

}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/



	MatrixXd Xsig_ = MatrixXd(5, 11);

	Xsig_.col(0) = x_;

	//calculate square root of P
	MatrixXd P_root = P_.llt().matrixL();

	double lambda_r = std::sqrt(lambda_ + n_x_);

	for (int index = 0; index < n_x_; ++index)
	{
		Xsig_.col(index + 1) = x_ + (lambda_r * P_root.col(index));
		Xsig_.col(index + n_x_ + 1) = x_ - (lambda_r * P_root.col(index));

	}

	/////Create Xsig_aug matrix
	// Create aug vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0);
	x_aug.head(5) = x_;

	// Create aug covariance matrix
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	//create aug square root matrix
	MatrixXd P_aug_root = P_aug.llt().matrixL();

	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.col(0) = x_aug;

	double lambda_r = std::sqrt(lambda_ + n_aug_);

	for (int index = 0; index < n_aug_; ++index)
	{
		Xsig_aug.col(index + 1) =
			x_aug + (lambda_r * P_aug_root.col(index));

		Xsig_aug.col(index + n_aug_ + 1) =
			x_aug - (lambda_r * P_aug_root.col(index));
	}
	/////aug_sigma_point matrix created





	////create sigma_point_prediction matrix
	MatrixXd predictions = MatrixXd(n_x_, 2 * n_aug_ + 1);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		VectorXd current_state = Xsig_aug.col(index);

		float longitudinal_speed = current_state(2);
		float yaw = current_state(3);
		float yaw_speed = current_state(4);

		float random_linear_acceleration = current_state(5);
		float random_yaw_acceleration = current_state(6);

		float speed_ratios = longitudinal_speed / yaw_speed;
		float interpolated_yaw = yaw + (delta_t * yaw_speed);

		VectorXd change_vector = VectorXd(n_x_);

		if (std::abs(yaw_speed) < 0.001)
		{
			change_vector(0) = delta_t * longitudinal_speed * std::cos(yaw);
			change_vector(1) = delta_t * longitudinal_speed * std::sin(yaw);

		}
		else
		{
			change_vector(0) = speed_ratios * (std::sin(interpolated_yaw) - std::sin(yaw));
			change_vector(1) = speed_ratios * (std::cos(yaw) - std::cos(interpolated_yaw));
		}

		change_vector(2) = 0;
		change_vector(3) = yaw_speed * delta_t;
		change_vector(4) = 0;

		VectorXd noise_vector = VectorXd(n_x_);
		double timeproduct = 0.5 * delta_t * delta_t;

		noise_vector(0) = timeproduct * std::cos(yaw) * random_linear_acceleration;
		noise_vector(1) = timeproduct * std::sin(yaw) * random_linear_acceleration;
		noise_vector(2) = delta_t * random_linear_acceleration;
		noise_vector(3) = timeproduct * random_yaw_acceleration;
		noise_vector(4) = delta_t * random_yaw_acceleration;

		predictions.col(index) = current_state.head(5) + change_vector + noise_vector;

	}
	Xsig_pred_ = predictions;
	////Xsig_pred_ is predicted_sigma_point matrix

	//x_ = mean of predictions

	x_ << 0, 0, 0, 0, 0;
	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		x_ += weights_(index) * Xsig_pred_.col(index);
	}

	// Normalize yaw angle
	while (x_(3) < -M_PI) { x_(3) += 2.0 * M_PI; }
	while (x_(3) >  M_PI) { x_(3) -= 2.0 * M_PI; }


	/////P_ = Predicted Coveriance Matrix

	MatrixXd covariance_matrix = MatrixXd(n_x_, n_x_);
	covariance_matrix.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		VectorXd diff = Xsig_pred_.col(index) - x_;

		// Normalize yaw angle 
		while (diff(3) < -M_PI) { diff(3) += 2.0 * M_PI; }
		while (diff(3) >  M_PI) { diff(3) -= 2.0 * M_PI; }

		covariance_matrix += weights_(index) * diff * diff.transpose();
	}
	P_ = covariance_matrix;
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/

	// Laser measurements have two dimensions, px and py
	int n_z = 2;

	//MatrixXd measurements_predictions = getLaserMeasurementsPredictions(Xsig_pred_);
	MatrixXd measurements_predictions = MatrixXd(2, 2 * n_aug_ + 1);
	//cout << "hallo" << endl;
	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		// predictions to laser measurements 
		measurements_predictions.col(index) = Xsig_pred_.col(index).head(2);
	}


	//prediction of measurement mean
	VectorXd mean_measurement_prediction = VectorXd(2);
	mean_measurement_prediction.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		mean_measurement_prediction += weights_(index) * measurements_predictions.col(index);
	}

	//	S = Laser_Measurement_Prediction_covariance matrix
	MatrixXd S = MatrixXd(2, 2);
	S.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		VectorXd diff = measurements_predictions.col(index) - mean_measurement_prediction;
		S += weights_(index) * diff * diff.transpose();
	}

	MatrixXd measurement_noise_covariance_matrix(2, 2);
	measurement_noise_covariance_matrix.fill(0);
	measurement_noise_covariance_matrix(0, 0) = std_laspx_ * std_laspx_;
	measurement_noise_covariance_matrix(1, 1) = std_laspy_ * std_laspy_;

	S += measurement_noise_covariance_matrix;

	//T = Laser_Cross_Correlation matrix
	MatrixXd T = MatrixXd(n_x_, 2);
	T.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{

		VectorXd state_diff = Xsig_pred_.col(index) - x_;
		VectorXd measurement_diff = measurements_predictions.col(index) - mean_measurement_prediction;

		T += weights_(index) * state_diff * measurement_diff.transpose();

	}


	MatrixXd K = T * S.inverse();

	VectorXd measurement(2);
	measurement(0) = meas_package.raw_measurements_(0);
	measurement(1) = meas_package.raw_measurements_(1);

	VectorXd measurement_diff = measurement - mean_measurement_prediction;

	// Update
	x_ = x_ + (K * measurement_diff);
	P_ = P_ - (K * S * K.transpose());

	// Normalize yaw angle
	while (x_(3) < -M_PI) { x_(3) += 2.0 * M_PI; }
	while (x_(3) >  M_PI) { x_(3) -= 2.0 * M_PI; }
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/

	// Radar measurements has three dimensions; radial distance, angle, velocity
	int n_z = 3;

	//create measurements_predictions matrix
	MatrixXd measurements_predictions = MatrixXd(3, 2 * n_aug_ + 1);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		double px = Xsig_pred_(0, index);
		double py = Xsig_pred_(1, index);
		double v = Xsig_pred_(2, index);
		double yaw = Xsig_pred_(3, index);

		double vx = v * std::cos(yaw);
		double vy = v * std::sin(yaw);

		double radial_distance = std::sqrt((px * px) + (py * py));

		measurements_predictions(0, index) = radial_distance;
		measurements_predictions(1, index) = std::atan2(py, px);
		measurements_predictions(2, index) = ((px * vx) + (py * vy)) / radial_distance;
	}



	//create vector with mean of predictions
	VectorXd mean_measurement_prediction = VectorXd(3);
	mean_measurement_prediction.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		mean_measurement_prediction += weights_(index) * measurements_predictions.col(index);
	}

	//	S = Radar_Measurement_Prediction_Covariance matrix
	MatrixXd S = MatrixXd(3, 3);
	S.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{
		VectorXd diff = measurements_predictions.col(index) - mean_measurement_prediction;
		//normalize angle
		while (diff(1) < -M_PI) { diff(1) += 2.0 * M_PI; }
		while (diff(1) >  M_PI) { diff(1) -= 2.0 * M_PI; }


		S += weights_(index) * diff * diff.transpose();
	}

	MatrixXd measurement_noise_covariance_matrix(3, 3);
	measurement_noise_covariance_matrix.fill(0);
	measurement_noise_covariance_matrix(0, 0) = std_radr_ * std_radr_;
	measurement_noise_covariance_matrix(1, 1) = std_radphi_ * std_radphi_;
	measurement_noise_covariance_matrix(2, 2) = std_radrd_ * std_radrd_;

	S += measurement_noise_covariance_matrix;



	//T = Radar_Cross_Correlation matrix

	MatrixXd T = MatrixXd(n_x_, 3);
	T.fill(0);

	for (int index = 0; index < 2 * n_aug_ + 1; ++index)
	{

		VectorXd state_diff = Xsig_pred_.col(index) - x_;
		//normalize angle
		while (state_diff(3) < -M_PI) { state_diff(3) += 2.0 * M_PI; }
		while (state_diff(3) >  M_PI) { state_diff(3) -= 2.0 * M_PI; }

		VectorXd measurement_diff = measurements_predictions.col(index) - mean_measurement_prediction;
		//normalize angle
		while (measurement_diff(1) < -M_PI) { measurement_diff(1) += 2.0 * M_PI; }
		while (measurement_diff(1) >  M_PI) { measurement_diff(1) -= 2.0 * M_PI; }


		T += weights_(index) * state_diff * measurement_diff.transpose();

	}


	MatrixXd K = T * S.inverse();

	VectorXd measurement(3);
	measurement(0) = meas_package.raw_measurements_(0);
	measurement(1) = meas_package.raw_measurements_(1);
	measurement(2) = meas_package.raw_measurements_(2);

	VectorXd measurement_diff = measurement - mean_measurement_prediction;
	//normalize angle
	while (measurement_diff(1) < -M_PI) { measurement_diff(1) += 2.0 * M_PI; }
	while (measurement_diff(1) >  M_PI) { measurement_diff(1) -= 2.0 * M_PI; }


	// Update
	x_ = x_ + (K * measurement_diff);
	P_ = P_ - (K * S * K.transpose());

	// Normalize yaw angle
	while (x_(3) < -M_PI) { x_(3) += 2.0 * M_PI; }
	while (x_(3) >  M_PI) { x_(3) -= 2.0 * M_PI; }
}










