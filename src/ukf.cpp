#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // Augmented state dimension
  n_x_ = 5;

  // Sigma point spreading parameter
  n_aug_ = n_x_ + 2;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0; //0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0; //0.15;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  //NIS for Laser
  NIS_laser_ = 0;

  //NIS for Radar
  NIS_radar_ = 0;
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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "UKF: " << endl;

    x_ << 0.0, //px [m]
          0.0, //py [m]
          0.0, //v [m/s]
          0.0, //yaw [rad]
          0.0; //yawd [rad/s]

    P_ << 1.0, 0, 0, 0, 0,
          0, 1.0, 0, 0, 0,
          0, 0, 1.0, 0, 0,
          0, 0, 0, 1.0, 0,
          0, 0, 0, 0, 1.0;

    double px, py;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      px = rho * cos(phi);
      py = rho * sin(phi);
    }

    x_(0) = px;
    x_(1) = py;

    // time when the state is true, in us
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   *  Updating
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  /*******************************************************************************
   *  Generating Augmented Sigma Points
   *******************************************************************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = pow(std_a_, 2);
  P_aug(n_x_+1, n_x_+1) = pow(std_yawdd_, 2);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  //calculate square root of P_aug
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //set remaining sigma points
  for (int i=0; i<n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_)*L.col(i);
  }

  /*******************************************************************************
   *  Sigma Point Prediction
   *******************************************************************************/
  for (int i=0; i<2*n_aug_+1; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predict sigma points
    double px_pred, py_pred, v_pred, yaw_pred, yawd_pred;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_pred = px + v/yawd * (sin(yaw+yawd*delta_t) - sin(yaw));
      py_pred = py + v/yawd * (-cos(yaw+yawd*delta_t) + cos(yaw));
    } else {
      px_pred = px + v * cos(yaw) * delta_t;
      py_pred = py + v * sin(yaw) * delta_t;
    }

    v_pred = v + 0;
    yaw_pred = yaw + yawd * delta_t;
    yawd_pred = yawd + 0;

    px_pred += 0.5 * pow(delta_t,2) * cos(yaw) * nu_a;
    py_pred += 0.5 * pow(delta_t,2) * sin(yaw) * nu_a;
    v_pred += delta_t * nu_a;
    yaw_pred += 0.5 * pow(delta_t,2) * nu_yawdd;
    yawd_pred += delta_t * nu_yawdd;

    VectorXd Xsig_pred_col = VectorXd(n_x_);
    Xsig_pred_col << px_pred,
                     py_pred,
                     v_pred,
                     yaw_pred,
                     yawd_pred;

    //write predicted sigma points into right column
    Xsig_pred_.col(i) = Xsig_pred_col;
  }

  /*******************************************************************************
   *  Predicted Mean and Covariance
   *******************************************************************************/
  for (int i=0; i<2*n_aug_+1; i++) {
    //set weights
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_+n_aug_);
    } else {
      weights_(i) = 1 / (2 * (lambda_+n_aug_));
    }
  }

  //predict state mean
  x_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2.0 * M_PI;
    }
    while (M_PI <= x_diff(3)) {
      x_diff(3) -= 2.0 * M_PI;
    }

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
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

  /*******************************************************************************
   *  Predict Laser Measurement
   *******************************************************************************/
  int n_z = 2;

  MatrixXd H = MatrixXd(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  /******************************

  z_pred =      H_     *  x_

                           /  px  \
  / px \   /1, 0, 0, 0, 0\ |  py  |
  |    | = |             | |  v   |
  \ py /   \0, 1, 0, 0, 0/ | yaw  |
                           \ yawd /

  *******************************/

  VectorXd z_pred = H * x_;

  /*******************************************************************************
   *  Linear KF Update for Laser
   *******************************************************************************/
  MatrixXd Ht = H.transpose();
  MatrixXd R = MatrixXd(n_z, n_z);
  R << pow(std_laspx_,2),                 0,
                       0, pow(std_laspy_,2);
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();

  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;

  //calculate NIS for Laser
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
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

  /*******************************************************************************
   *  Predict Radar Measurement
   *******************************************************************************/
  int n_z = 3;

  //create sigma point matrix in measurment space
  MatrixXd Zsig_aug = MatrixXd(n_z, 2*n_aug_+1);

  for (int i=0; i<2*n_aug_+1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double rho = sqrt(pow(px,2) + pow(py,2));
    double phi = atan2(py, px);
    double rhod = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;

    VectorXd Zsig_aug_col = VectorXd(n_z);
    Zsig_aug_col << rho,
                    phi,
                    rhod;

    Zsig_aug.col(i) = Zsig_aug_col;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //calculate mean predicted measurement
  for (int i=0; i<2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig_aug.col(i);
  }

  /*******************************************************************************
   *  Unscented KF Update for Radar
   *******************************************************************************/
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //calculate innovation covariance matrix S
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig_aug.col(i) - z_pred;
    // angle normalization for phi [rad]
    while (z_diff(1) < -M_PI) {
      z_diff(1) += 2.0 * M_PI;
    }
    while (M_PI <= z_diff(1)) {
      z_diff(1) -= 2.0 * M_PI;
    }

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << pow(std_radr_,2),                  0,                 0,
                      0, pow(std_radphi_,2),                 0,
                      0,                  0, pow(std_radrd_,2);
  S += R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization for yaw [rad]
    while (x_diff(3) < -M_PI) {
      x_diff(3) += 2.0 * M_PI;
    }
    while (M_PI <= x_diff(3)) {
      x_diff(3) -= 2.0 * M_PI;
    }

    VectorXd z_diff = Zsig_aug.col(i) - z_pred;
    // angle normalization for phi [rad]
    while (z_diff(1) < -M_PI) {
      z_diff(1) += 2.0 * M_PI;
    }
    while (M_PI <= z_diff(1)) {
      z_diff(1) -= 2.0 * M_PI;
    }

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  // angle normalization for phi [rad]
  while (z_diff(1) < -M_PI) {
    z_diff(1) += 2.0 * M_PI;
  }
  while (M_PI <= z_diff(1)) {
    z_diff(1) -= 2.0 * M_PI;
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //calculate NIS for Radar
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
