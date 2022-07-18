/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#pragma once
#include <ros/assert.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "integration_base.h"

#include "../utility/utility.h"
#include "../estimator/parameters.h"
#include "../utility/sophus_utils.hpp"
#include <ceres/ceres.h>

class ZUPTFactor : public ceres::SizedCostFunction<6, 7, 9, 7, 9>
{
  public:
    ZUPTFactor(){};
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {

        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);

        Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);


        Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
        residual.block<3, 1>(0, 0) = dR * R0 * RIC[0].transpose() * Qi.inverse() * Vi;
        residual.block<3, 1>(3, 0) = dR * R0 * RIC[0].transpose() * Qj.inverse() * Vj;

        Matrix<double, 6, 6> sqrt_info = 100000 * Eigen::Matrix<double, 6, 6>::Identity();
        residual = sqrt_info * residual;
        ROS_WARN_STREAM("residual " << residual.transpose());

        if (jacobians)
        {

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                jacobian_pose_i.setZero();

                jacobian_pose_i.block<3, 3>(0, 3) = dR * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qi.inverse() * Vi);

                jacobian_pose_i = sqrt_info * jacobian_pose_i;

//                ROS_ASSERT(fabs(jacobian_pose_j.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_pose_j.minCoeff()) < 1e8);
            }

            if(jacobians[1]){
                Eigen::Map<Eigen::Matrix<double, 6, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                jacobian_speedbias_i.setZero();

                jacobian_speedbias_i.block<3, 3>(0, 0) = dR * R0 * RIC[0].transpose() * Qi.inverse();

                jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

//                ROS_ASSERT(fabs(jacobian_ix_sx.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_ix_sx.minCoeff()) < 1e8);
            }

            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                jacobian_pose_j.setZero();

                jacobian_pose_j.block<3, 3>(3, 3) = -dR * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qj.inverse() * Vj);

                jacobian_pose_j = sqrt_info * jacobian_pose_j;

//                ROS_ASSERT(fabs(jacobian_pose_j.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_pose_j.minCoeff()) < 1e8);
            }

            if(jacobians[3]){
                Eigen::Map<Eigen::Matrix<double, 6, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                jacobian_speedbias_j.setZero();

                jacobian_speedbias_j.block<3, 3>(3, 0) = dR * R0 * RIC[0].transpose() * Qj.inverse();

                jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

//                ROS_ASSERT(fabs(jacobian_ix_sx.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_ix_sx.minCoeff()) < 1e8);
            }

        }

        return true;
    }

};

