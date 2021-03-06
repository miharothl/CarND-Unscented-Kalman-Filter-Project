#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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

  verbose_= false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  //Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  //radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  //radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  //radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //set state dimension
  n_x_ = 5;

  //define spreading parameter
  lambda_ = 3 - n_x_;

  //set augmented dimension
  n_aug_ = 7;


  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.);
  x_pred_ = VectorXd(5);
  x_pred_.fill(0.);
  P_pred_ = MatrixXd(5, 5);
  P_pred_.fill(0.);


  //vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
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

  if (!this->is_initialized_) {

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];

      float px = ro * cos(phi);
      float py = ro * sin(phi);

      x_ << px, py, ro_dot, 0., 0. ;

      is_initialized_ = true;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];

      x_ << px, py, 0., 0., 0. ;

      is_initialized_ = true;
    }

    previous_timestamp_ = meas_package.timestamp_;

    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  this->Prediction(dt);

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    this->UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    this->UpdateLidar(meas_package);
  }

  return;
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

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  this->AugmentedSigmaPoints(&Xsig_aug);
  this->SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred_);
  this->PredictMeanAndCovariance(Xsig_pred_, &x_pred_, &P_pred_);
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
  int n_z = 2;

  VectorXd z_pred = VectorXd(2);
  MatrixXd S_pred = MatrixXd(2, 2);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S_pred.fill(0.0);
  this->PredictLidarMeasurement(Xsig_pred_, &z_pred, &S_pred, &Zsig);

  VectorXd x_out = VectorXd(5);
  MatrixXd P_out = MatrixXd(5, 5);
  this->UpdateStateFromLidarMeasurement(meas_package, z_pred, S_pred, Zsig, &x_out, &P_out);
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

  int n_z = 3;

  VectorXd z_pred = VectorXd(3);
  MatrixXd S_pred = MatrixXd(3, 3);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  this->PredictRadarMeasurement(Xsig_pred_, &z_pred, &S_pred, &Zsig);

  VectorXd x_out = VectorXd(5);
  MatrixXd P_out = MatrixXd(5, 5);
  this->UpdateStateFromRadarMeasurement(meas_package, z_pred, S_pred, Zsig, &x_out, &P_out);
}


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  if (verbose_)
    cout << "AugmentedSigmaPoints..." << endl;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.fill(0.);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  if (verbose_)
    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;
}


void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double dt, MatrixXd* Xsig_out) {

  if (verbose_)
    cout << "SigmaPointPrediction..." << endl;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double delta_t = dt;

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  if (verbose_)
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;
}


void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {

  if (verbose_)
    cout << "PredictMeanAndCovariance..." << endl;

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    x_diff(3)  = this->NormalizeAngle(x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  if (verbose_) {
    std::cout << "Predicted state" << std::endl;
    std::cout << x << std::endl;
    std::cout << "Predicted covariance matrix" << std::endl;
    std::cout << P << std::endl;
  }

  //write result
  *x_out = x;
  *P_out = P;
}


void UKF::PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {

  if (verbose_)
    cout << "PredictRadarMeasurement..." << endl;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = this->NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  if (verbose_) {
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;
  }

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}


void UKF::PredictLidarMeasurement(MatrixXd Xsig_pred, VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {

  if (verbose_)
    cout << "PredictLidarMesurement..." << endl;

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<   std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  S = S + R;

  if (verbose_) {
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;
  }

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}


void UKF::UpdateStateFromRadarMeasurement(MeasurementPackage meas_package, VectorXd z_pred, MatrixXd S_pred,
                                          MatrixXd Zsig, VectorXd *x_out, MatrixXd *P_out) {

  if (verbose_)
  cout << "UpdateStateFromRadarMeasurement..." << endl;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1),
          meas_package.raw_measurements_(2);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = this->NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;
    //angle normalization
    x_diff(3) = this->NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd S_pred_inverse = S_pred.inverse();
  MatrixXd K = Tc * S_pred_inverse;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = this->NormalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_pred_ + K * z_diff;
  P_ = P_pred_ - K*S_pred*K.transpose();

  if (verbose_) {
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  }

  NIS_radar_ = z_diff.transpose() * S_pred_inverse * z_diff;

  //write result
  *x_out = x_;
  *P_out = P_;
}


void UKF::UpdateStateFromLidarMeasurement(MeasurementPackage meas_package, VectorXd z_pred, MatrixXd S_pred,
                                          MatrixXd Zsig, VectorXd *x_out, MatrixXd *P_out) {

  if (verbose_)
    cout << "UpdateStateFromLidarMeasurement..." << endl;

  //set measurement dimension, laser can measure px, py
  int n_z = 2;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;
    //angle normalization
    x_diff(3) = this->NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd S_pred_inverse = S_pred.inverse();
  MatrixXd K = Tc * S_pred_inverse;

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_pred_ + K * z_diff;
  P_ = P_pred_ - K*S_pred*K.transpose();

  NIS_laser_ = z_diff.transpose() * S_pred_inverse * z_diff;

  if (verbose_) {
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  }

  //write result
  *x_out = x_;
  *P_out = P_;
}


double UKF::NormalizeAngle(double angle) {
  return atan2(sin(angle), cos(angle));
}