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

class NonholonomicFactor : public ceres::SizedCostFunction<3, 7, 9, 7, 9, 7>
{
  public:
    NonholonomicFactor(IntegrationBase* _pre_integration): pre_integration(_pre_integration){};
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {

        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);

        Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);

        Eigen::Vector3d TIO(parameters[4][0], parameters[4][1], parameters[4][2]);
        Eigen::Quaterniond dR(parameters[4][6], parameters[4][3], parameters[4][4], parameters[4][5]);

#if 0
        if ((Bai - pre_integration->linearized_ba).norm() > 0.10 ||
            (Bgi - pre_integration->linearized_bg).norm() > 0.01)
        {
            pre_integration->repropagate(Bai, Bgi);
        }
#endif

        Eigen::Map<Eigen::Vector3d> residual(residuals);
        Vector3d dV = pre_integration->vel_1 - pre_integration->linearized_vel;
        Matrix3d R_w_0;
        Vector3d w_x_0 = pre_integration->linearized_gyr - pre_integration->linearized_bg;
        R_w_0 << 0, -w_x_0(2), w_x_0(1),
                w_x_0(2), 0, -w_x_0(0),
                -w_x_0(1), w_x_0(0), 0;

        Matrix3d R_w_1;
        Vector3d w_x_1 = pre_integration->gyr_1 - pre_integration->linearized_bg;
        R_w_1 << 0, -w_x_1(2), w_x_1(1),
                w_x_1(2), 0, -w_x_1(0),
                -w_x_1(1), w_x_1(0), 0;
        residual = dV - dR * R0 * RIC[0].transpose() * (Qj.inverse() * Vj - Qi.inverse() * Vi - R_w_1 * TIO + R_w_0 * TIO);

        Matrix3d A = dR * R0 * RIC[0].transpose() * Utility::skewSymmetric(TIO);
        double dt = pre_integration->sum_dt;
        Matrix3d cov = 2 * (VEL_N_wheel * VEL_N_wheel) * Eigen::Matrix3d::Identity() + 2 * A * ((GYR_N * GYR_N) * Eigen::Matrix3d::Identity() + dt * (GYR_W * GYR_W) * Eigen::Matrix3d::Identity()) * A.transpose();
//        Matrix3d cov = 2 * (VEL_N_wheel * VEL_N_wheel) * Eigen::Matrix3d::Identity() + 2 * (GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
        Eigen::Matrix<double, 3, 3> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 3, 3>>(cov.inverse()).matrixL().transpose();
        residual = sqrt_info * residual;

        if (jacobians)
        {

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]); //jacobians[i] is a row-major array of size num_residuals x parameter_block_sizes_[i]
                jacobian_pose_i.setZero();

                jacobian_pose_i.block<3, 3>(0, 0) =  dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qi.inverse() * Vi);

                jacobian_pose_i = sqrt_info * jacobian_pose_i;

                if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
                {
                    ROS_WARN("numerical unstable in preintegration");
                    //std::cout << sqrt_info << std::endl;
                    //ROS_BREAK();
                }
            }
            if (jacobians[1])
            {
                Eigen::Map<Eigen::Matrix<double, 3, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                jacobian_speedbias_i.setZero();

                jacobian_speedbias_i.block<3, 3>(0, 0) = dR.toRotationMatrix() * R0 * RIC[0].transpose() * Qi.inverse();

                jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

//                ROS_ASSERT(fabs(jacobian_pose_j.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_pose_j.minCoeff()) < 1e8);
            }
            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                jacobian_pose_j.setZero();

                jacobian_pose_j.block<3, 3>(0, 0) = -dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qj.inverse() * Vj);

                jacobian_pose_j = sqrt_info * jacobian_pose_j;

//                ROS_ASSERT(fabs(jacobian_ex_pose.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_ex_pose.minCoeff()) < 1e8);
            }

            if(jacobians[3]){
                Eigen::Map<Eigen::Matrix<double, 3, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                jacobian_speedbias_j.setZero();

                jacobian_speedbias_j.block<3, 3>(0, 0) = -dR.toRotationMatrix() * R0 * RIC[0].transpose() * Qj.inverse();


                jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

//                ROS_ASSERT(fabs(jacobian_ix_sx.maxCoeff()) < 1e8);
//                ROS_ASSERT(fabs(jacobian_ix_sx.minCoeff()) < 1e8);
            }
            if (jacobians[4])
            {
                Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_tio(jacobians[4]);
                jacobian_tio.setZero();

                jacobian_tio.block<3, 3>(0, 0) = dR * R0 * RIC[0].transpose() * (R_w_1 - R_w_0);

                jacobian_tio.block<3, 3>(0, 3) = dR * Utility::skewSymmetric(R0 * RIC[0].transpose() * (Qj.inverse() * Vj - Qi.inverse() * Vi - R_w_1 * TIO + R_w_0 * TIO));

                jacobian_tio = sqrt_info * jacobian_tio;

                //ROS_ASSERT(fabs(jacobian_speedbias_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_j.minCoeff()) < 1e8);
            }

        }

        return true;
    }

    IntegrationBase* pre_integration;

};

