#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  //lesson 5, section 23
	VectorXd error_sum(4);
	error_sum << 0, 0, 0, 0;

	for (int i = 0; i < estimations.size(); ++i)
	{
		VectorXd difference = estimations[i] - ground_truth[i];
		VectorXd squared_difference = difference.array() * difference.array();

		error_sum += squared_difference;
	}

	VectorXd MSE = error_sum / estimations.size();
	return MSE.array().sqrt();*/
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	for (int i = 0; i < estimations.size(); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient wise multiplication
		residual = residual.array() * residual.array();

		rmse += residual;
	}
	// calculate the mean
	rmse = rmse / estimations.size();
	//calculate the squared root
	rmse = rmse.array().sqrt();
	cout << rmse << "=rmse" << endl;
	return rmse;
}