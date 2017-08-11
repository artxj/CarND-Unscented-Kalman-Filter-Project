#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * Returns normalized angle in range [-pi; pi]
  */
  double NormalizeAngle(double angle);

  /**
  * A helper method to calculate Cartesian coordinates from Polar ones
  */
  void CalculateCartesian(const double &rho, const double &phi, double &px, double &py);

};

#endif /* TOOLS_H_ */
