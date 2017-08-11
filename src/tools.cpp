#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  * Calculate the RMSE here.
  */
  int size = estimations.size();
  assert(size != 0 && ground_truth.size() == size);

  VectorXd rmse(4);
	rmse.setZero();

	for (int i = 0; i < size; ++i) {
    VectorXd res = estimations[i] - ground_truth[i];
    res = res.array() * res.array();
    rmse += res;
	}

	//calculating the mean
	rmse = rmse.array() / size;

	//calculating the squared root
	rmse = rmse.array().sqrt();

	return rmse;
}

double Tools::NormalizeAngle(double angle) {
  double result = angle;
  double twopi = 2 * M_PI;
  while (result < -M_PI) result += twopi;
  while (result > M_PI) result -= twopi;
  return result;
}

void Tools::CalculateCartesian(const double &rho, const double &phi, double &px, double &py) {
  px = rho * sin(phi);
  py = -rho * cos(phi);
}
