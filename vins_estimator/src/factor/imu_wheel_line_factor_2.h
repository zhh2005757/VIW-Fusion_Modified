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

#include "../utility/utility.h"
#include "../estimator/parameters.h"
#include "integration_base.h"

#include <ceres/ceres.h>
#include <ceres/gradient_checker.h>
#include "pose_local_parameterization.h"

class IMUWheelLineFactor : public ceres::SizedCostFunction<18, 7, 9, 7, 9, 7>
{
  public:
    IMUWheelLineFactor() = delete;
    IMUWheelLineFactor(IntegrationBase* _pre_integration): pre_integration(_pre_integration)
    {
    }
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
    {

        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
        Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

        Eigen::Vector3d Pj(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond Qj(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
        Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
        Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);

        Eigen::Vector3d TIO(parameters[4][0], parameters[4][1], parameters[4][2]);
        Eigen::Quaterniond dR(parameters[4][6], parameters[4][3], parameters[4][4], parameters[4][5]);

//Eigen::Matrix<double, 15, 15> Fd;
//Eigen::Matrix<double, 15, 12> Gd;

//Eigen::Vector3d pPj = Pi + Vi * sum_t - 0.5 * g * sum_t * sum_t + corrected_delta_p;
//Eigen::Quaterniond pQj = Qi * delta_q;
//Eigen::Vector3d pVj = Vi - g * sum_t + corrected_delta_v;
//Eigen::Vector3d pBaj = Bai;
//Eigen::Vector3d pBgj = Bgi;

//Vi + Qi * delta_v - g * sum_dt = Vj;
//Qi * delta_q = Qj;

//delta_p = Qi.inverse() * (0.5 * g * sum_dt * sum_dt + Pj - Pi);
//delta_v = Qi.inverse() * (g * sum_dt + Vj - Vi);
//delta_q = Qi.inverse() * Qj;

#if 0
        if ((Bai - pre_integration->linearized_ba).norm() > 0.10 ||
            (Bgi - pre_integration->linearized_bg).norm() > 0.01)
        {
            pre_integration->repropagate(Bai, Bgi);
        }
#endif

        Eigen::Map<Eigen::Matrix<double, 18, 1>> residual(residuals);
        residual = pre_integration->evaluate(Pi, Qi, Vi, Bai, Bgi,
                                            Pj, Qj, Vj, Baj, Bgj, TIO, dR);

        Eigen::Matrix<double, 18, 18> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 18, 18>>(pre_integration->covariance.inverse()).matrixL().transpose();
//        Eigen::Matrix<double, 18, 18> sqrt_info_1 = Eigen::Matrix<double, 18, 18>::Zero();
//        sqrt_info_1.block<15,15>(0,0) = Eigen::LLT<Eigen::Matrix<double, 15, 15>>(pre_integration->covariance.block<15,15>(0,0).inverse()).matrixL().transpose();
        //sqrt_info.setIdentity();
//        std::cout<<"sqrt_info_1 :\n"<<sqrt_info_1<<std::endl;
//        std::cout<<"sqrt_info :\n"<<sqrt_info<<std::endl;

        residual = sqrt_info * residual;

        Matrix3d R_w_0;
//        Vector3d w_x_0 = -R0 * RIC[0].transpose() * (pre_integration->gyr_0 - pre_integration->linearized_bg);
        Vector3d w_x_0 = pre_integration->linearized_gyr - Bgi;
//        ROS_WARN_STREAM("w_x_0 " << w_x_0.transpose());
        R_w_0 << 0, -w_x_0(2), w_x_0(1),
                w_x_0(2), 0, -w_x_0(0),
                -w_x_0(1), w_x_0(0), 0;


        Matrix3d R_w_1;
//        Vector3d w_x_1 = -R0 * RIC[0].transpose() * (pre_integration->gyr_1 - pre_integration->linearized_bg);
        Vector3d w_x_1 = pre_integration->gyr_1 - Bgj;
//        ROS_WARN_STREAM("w_x_1 " << w_x_1.transpose());
        R_w_1 << 0, -w_x_1(2), w_x_1(1),
                w_x_1(2), 0, -w_x_1(0),
                -w_x_1(1), w_x_1(0), 0;

        Vector3d vel_0 = pre_integration->linearized_vel;
//        ROS_WARN_STREAM("vel_0 " << vel_0.transpose());
        Vector3d vel_1 = pre_integration->vel_1;


        if (jacobians)
        {
            double sum_dt = pre_integration->sum_dt;
            Eigen::Matrix3d dp_dba = pre_integration->jacobian.template block<3, 3>(O_P, O_BA);
            Eigen::Matrix3d dp_dbg = pre_integration->jacobian.template block<3, 3>(O_P, O_BG);

            Eigen::Matrix3d dq_dbg = pre_integration->jacobian.template block<3, 3>(O_R, O_BG);

            Eigen::Matrix3d dv_dba = pre_integration->jacobian.template block<3, 3>(O_V, O_BA);
            Eigen::Matrix3d dv_dbg = pre_integration->jacobian.template block<3, 3>(O_V, O_BG);

            if (pre_integration->jacobian.maxCoeff() > 1e8 || pre_integration->jacobian.minCoeff() < -1e8)
            {
                ROS_WARN("numerical unstable in preintegration");
                //std::cout << pre_integration->jacobian << std::endl;
///                ROS_BREAK();
            }

            if (jacobians[0])
            {
                Eigen::Map<Eigen::Matrix<double, 18, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
                jacobian_pose_i.setZero();

                jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
                jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

#if 0
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Qj.inverse() * Qi).toRotationMatrix();
#else
                Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q)).bottomRightCorner<3, 3>();
#endif

                jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (G * sum_dt + Vj - Vi));

                jacobian_pose_i.block<3, 3>(O_VO, O_R) = dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qi.inverse() * Vi);

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
                Eigen::Map<Eigen::Matrix<double, 18, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
                jacobian_speedbias_i.setZero();
                jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
                jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
                jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;

#if 0
            jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -dq_dbg;
#else
                //Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                //jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * corrected_delta_q).bottomRightCorner<3, 3>() * dq_dbg;
                jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration->delta_q).bottomRightCorner<3, 3>() * dq_dbg;
#endif

                jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
                jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
                jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

                jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

                jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

                jacobian_speedbias_i.block<3, 3>(O_VO, O_V - O_V) = dR.toRotationMatrix() * R0 * RIC[0].transpose() * Qi.inverse();
                jacobian_speedbias_i.block<3, 3>(O_VO, O_BG - O_V) = -dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(TIO);

                jacobian_speedbias_i = sqrt_info * jacobian_speedbias_i;

                //ROS_ASSERT(fabs(jacobian_speedbias_i.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_i.minCoeff()) < 1e8);
            }
            if (jacobians[2])
            {
                Eigen::Map<Eigen::Matrix<double, 18, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[2]);
                jacobian_pose_j.setZero();

                jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();

#if 0
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Eigen::Matrix3d::Identity();
#else
                Eigen::Quaterniond corrected_delta_q = pre_integration->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration->linearized_bg));
                jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

                jacobian_pose_j.block<3, 3>(O_VO, O_R) = -dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(Qj.inverse() * Vj);
#endif

                jacobian_pose_j = sqrt_info * jacobian_pose_j;

                //ROS_ASSERT(fabs(jacobian_pose_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_pose_j.minCoeff()) < 1e8);
            }
            if (jacobians[3])
            {
                Eigen::Map<Eigen::Matrix<double, 18, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
                jacobian_speedbias_j.setZero();

                jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

                jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

                jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

                jacobian_speedbias_j.block<3, 3>(O_VO, O_V - O_V) = -dR.toRotationMatrix() * R0 * RIC[0].transpose() * Qj.inverse();
                jacobian_speedbias_j.block<3, 3>(O_VO, O_BG - O_V) = dR.toRotationMatrix() * R0 * RIC[0].transpose() * Utility::skewSymmetric(TIO);

                jacobian_speedbias_j = sqrt_info * jacobian_speedbias_j;

                //ROS_ASSERT(fabs(jacobian_speedbias_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_j.minCoeff()) < 1e8);
            }
            if (jacobians[4])
            {
                Eigen::Map<Eigen::Matrix<double, 18, 7, Eigen::RowMajor>> jacobian_tio(jacobians[4]);
                jacobian_tio.setZero();

                jacobian_tio.block<3, 3>(O_VO, O_P) = dR * R0 * RIC[0].transpose() * (R_w_1 - R_w_0);

                jacobian_tio.block<3, 3>(O_VO, O_R) = dR * Utility::skewSymmetric(R0 * RIC[0].transpose() * (Qj.inverse() * Vj - Qi.inverse() * Vi - R_w_1 * TIO + R_w_0 * TIO));

//                Matrix3d tmp_jaco = jacobian_tio.block<3, 3>(O_P, O_R);

//                ROS_WARN_STREAM("jacobian_tio.block<3, 3>(O_P, O_R) " << tmp_jaco);

                jacobian_tio = sqrt_info * jacobian_tio;

//                ROS_WARN_STREAM("jacobian  " << jacobian_tio);

                //ROS_ASSERT(fabs(jacobian_speedbias_j.maxCoeff()) < 1e8);
                //ROS_ASSERT(fabs(jacobian_speedbias_j.minCoeff()) < 1e8);
            }
        }

        return true;
    }

    //bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

    //void checkCorrection();
    //void checkTransition();
    //void checkJacobian(double **parameters);
    IntegrationBase* pre_integration;

};

