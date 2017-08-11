#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using namespace std::placeholders;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
#ifdef UKF_DEBUG
  std::cout << "Init..." << std::endl;
#endif

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  n_x_ = 5;

  n_aug_ = 7;

  // spreading parameter
  lambda_ = 3 - n_aug_;

  // sigma points predicted
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.setConstant(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // augmented covariance matrix
  Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

  // matrices of measurements std deviations
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  // laser measurement matrix
  H_ = MatrixXd(2, n_x_);
	H_ << 1, 0, 0, 0, 0,
			  0, 1, 0, 0, 0;

  Ht_ = H_.transpose();

  // identity matrix
  I_ = MatrixXd::Identity(n_x_, n_x_);

  // tools
  tools = Tools();

#ifdef UKF_DEBUG
  std::cout << "...done!" << std::endl;
#endif
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::RADAR: {
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double px, py;
        tools.CalculateCartesian(rho, phi, px, py);
        x_ << px, py, 0, 0, 0;
        break;
      }
      case MeasurementPackage::LASER: {
        double px = meas_package.raw_measurements_[0];
        double py = meas_package.raw_measurements_[1];
        x_ << px, py, 0, 0, 0;
        break;
      }
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
 	time_us_ = meas_package.timestamp_;

  Prediction(dt);
  switch (meas_package.sensor_type_) {
    case MeasurementPackage::RADAR: {
      if (!use_radar_) return;
      UpdateRadar(meas_package);
      break;
    }
    case MeasurementPackage::LASER: {
      if (!use_laser_) return;
      UpdateLidar(meas_package);
      break;
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
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
#ifdef UKF_DEBUG
  std::cout << "Prediction step..." << std::endl;
#endif

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.setZero();
  x_aug.head(n_x_) = x_;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;

  // calculate square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  MatrixXd addSigmas = sqrt(lambda_ + n_aug_) * A;
  for (int i = 0; i < n_aug_; ++i) {
      Xsig_aug.col(i + 1) = x_aug + addSigmas.col(i);
      Xsig_aug.col(n_aug_ + i + 1) = x_aug - addSigmas.col(i);
  }

  // predicting sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    //extract values
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double delta_t_sq = 0.5 * delta_t * delta_t;
    double siny = sin(yaw);
    double cosy = cos(yaw);

    double px_p, py_p;
    double yawd_p = yawd + delta_t * nu_yawdd;
    double yaw_p = yaw + delta_t * yawd + delta_t_sq * nu_yawdd;
    double v_p = v + delta_t * nu_a;
    if (fabs(yawd) < 1e-2) {
      py_p = py + v * siny * delta_t + delta_t_sq * siny * nu_a;
      px_p = px + v * cosy * delta_t + delta_t_sq * cosy * nu_a;
    } else {
      double v_y = v / yawd;
      py_p = py + v_y * (-cos(yaw + yawd * delta_t) + cosy) + delta_t_sq * siny * nu_a;
      px_p = px + v_y * (sin(yaw + yawd * delta_t) - siny) + delta_t_sq * cosy * nu_a;
    }

    VectorXd x_out = VectorXd(n_x_);
    x_out << px_p, py_p, v_p, yaw_p, yawd_p;

    Xsig_pred_.col(i) = x_out;
  }

  // predicted state mean
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    x_ +=  weights_(i) * Xsig_pred_.col(i);
  }
  x_(3) = tools.NormalizeAngle(x_(3));

  // predicted state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.NormalizeAngle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

#ifdef UKF_DEBUG
  std::cout << "...done!" << std::endl;
#endif

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
#ifdef UKF_DEBUG
  std::cout << "Update laser step..." << std::endl;
#endif

  int n_z = 2;
  VectorXd measurement = VectorXd(n_z);
  measurement << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  VectorXd y = measurement - H_ * x_;

  // support matrices
	MatrixXd S = H_ * P_ * Ht_ + R_laser_;
  MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht_ * Si;

	// new estimate
	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;

  NIS_laser_ = y.transpose() * Si * y;

#ifdef UKF_SHOW_NIS
  std::cout << "NIS laser = " << NIS_laser_ << std::endl;
#endif

#ifdef UKF_DEBUG
  std::cout << "...done!" << std::endl;
#endif
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

#ifdef UKF_DEBUG
  std::cout << "Update radar step..." << std::endl;
#endif

  int n_z = 3;
  VectorXd measurement = VectorXd(n_z);
  measurement << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],
                  meas_package.raw_measurements_[2];

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double ksi = Xsig_pred_(3, i);
    // double ksi_dot = Xsig_pred_(4, i);

    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rho_dot = (fabs(rho) < 1e-2) ? 0 : v * (px * cos(ksi) + py * sin(ksi)) / rho;
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rho_dot;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.setZero();

  // matrix of differences (sigma points - z_pred)
  MatrixXd z_diff = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd vec_diff = Zsig.col(i) - z_pred;
    vec_diff(1) = tools.NormalizeAngle(vec_diff(1));
    S += weights_(i) * vec_diff * vec_diff.transpose();

    // store diff value - will be used later too
    z_diff.col(i) = vec_diff;
  }

  // update S using R matrix
  S += R_radar_;

  // calculate the matrix for cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.NormalizeAngle(x_diff(3));
    Tc += weights_(i) * x_diff * z_diff.col(i).transpose();
  }

  // storing inversed S - will be used twice
  MatrixXd Si = S.inverse();

  // calculate Kalman gain K
  MatrixXd K = Tc * Si;

  //update state mean and covariance matrix

  VectorXd zvec_diff = measurement - z_pred;
  zvec_diff(1) = tools.NormalizeAngle(zvec_diff(1));
  x_ += K * zvec_diff;
  x_(3) = tools.NormalizeAngle(x_(3));
  P_ -= K * S * K.transpose();

  // calculate NIS
  NIS_radar_ = zvec_diff.transpose() * Si * zvec_diff;

#ifdef UKF_SHOW_NIS
  std::cout << "NIS radar = " << NIS_radar_ << std::endl;
#endif

#ifdef UKF_DEBUG
  std::cout << "...done!" << std::endl;
#endif

}
