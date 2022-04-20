/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#include "Rwg_SO3_local_parameterization.h"
RwgSO3LocalParameterization::RwgSO3LocalParameterization(const std::vector<int>& constant_parameters):constancy_mask_(4, 0){
    if(constant_parameters.empty())
        return;
    std::vector<int> constant = constant_parameters;
    std::sort(constant.begin(), constant.end());
//    CHECK_GE(constant.front(), 0)
//        << "Indices indicating constant parameter must be greater than zero.";
//    CHECK_LT(constant.back(), 4)
//        << "Indices indicating constant parameter must be less than the size "
//        << "of the parameter block.";
//    CHECK(std::adjacent_find(constant.begin(), constant.end()) == constant.end())
//    << "The set of constant parameters cannot contain duplicates";
    for (size_t i = 0; i < constant_parameters.size(); ++i) {
        constancy_mask_[constant_parameters[i]] = 1;
    }
}
bool RwgSO3LocalParameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
{
    Eigen::Map<const Sophus::SO3d> R(x);

    double delta_[LocalSize()];
    mempcpy(delta_, delta, sizeof(double) * LocalSize());
    for (int i = 0; i < 3; ++i) {
        delta_[i] = constancy_mask_[i] ? 0 : delta_[i];
    }

//    std::cout << "delta_  " << delta_[0] << " " << delta_[1] << " " << delta_[2] << std::endl;

    Sophus::SO3d R_delta = Sophus::SO3d::exp(Eigen::Map<const Eigen::Vector3d>(delta_));

    Eigen::Map<Sophus::SO3d> R_plus_delta(x_plus_delta);

    R_plus_delta = R * R_delta;

    std::cout << R.unit_quaternion().coeffs().transpose() << std::endl;

    return true;
}
bool RwgSO3LocalParameterization::ComputeJacobian(const double *x, double *jacobian) const
{
    Eigen::Map<const Sophus::SO3d> R(x);
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>> j(jacobian);//关于该jacobian的说明参加 https://fzheng.me/2018/01/23/ba-demo-ceres/
//    j.topRows<3>().setIdentity();
//    j.bottomRows<1>().setZero();
    j = R.Dx_this_mul_exp_x_at_0();

//    std::cout << "j " << j << std::endl;

    return true;
}
