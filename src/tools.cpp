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
  */
  VectorXd res;
  VectorXd rmse(4);
  VectorXd sum(4);
  rmse << 0,0,0,0;
  sum << 0,0,0,0;

   // check the validity of the following inputs:
  if ((estimations.size()==0) || (estimations.size() != ground_truth.size()))
  {
	    cout<< "Error: estimation or ground_truth does not have a good size."<< endl;
	    return rmse;
	}
  //accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code 
		res = (estimations[i] - ground_truth[i]);
		res = res.array()*res.array();
		rmse = rmse+res; 
	}
  rmse/=estimations.size();
	rmse = rmse.array().sqrt();
  return rmse;
 
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

 	MatrixXd Hj(3,4);
	//recover state parameters


	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
 /*
  if ((pow(px,2) + pow(py,2)) >= 0.0001)
	{
	     Hj << px/sqrt(pow(px,2) + pow(py,2)),py/sqrt(pow(px,2) + pow(py,2)),0,0,
            -py/(pow(px,2) + pow(py,2)),px/(pow(px,2) + pow(py,2)),0,0,
            (px*(vx*py-vy*px))/pow(pow(px,2) + pow(py,2),3/2),(py*(vy*px-vx*py))/pow(pow(px,2) + pow(py,2),3/2),px/sqrt(pow(px,2) + pow(py,2)),py/sqrt(pow(px,2) + pow(py,2));
	}
	else
	{
	    cout<<"Error: Division by zero"<<endl;
	}
  return Hj;
  */


 //pre-compute a set of terms to avoid repeated calculation

	double c1 = px*px+py*py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);



	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "ERROR - CalculateJacobian () - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;


}
