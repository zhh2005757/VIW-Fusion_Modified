/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Qin Tong (qintonguav@gmail.com)
 *******************************************************/

#pragma once
#include <eigen3/Eigen/Dense>
#include <iostream>
//#include "../factor/imu_factor.h"
#include "../factor/imu_wheel_factor.h"
#include "../utility/utility.h"
#include <ros/ros.h>
#include <map>
#include "../estimator/feature_manager.h"
#include "../factor/wheel_integration_base.h"

using namespace Eigen;
using namespace std;

class ImageFrame
{
    public:
        ImageFrame(){};
        ImageFrame(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& _points, double _t):t{_t},is_key_frame{false}
        {
            points = _points;
        };
        map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>> > > points;
        double t;
        Matrix3d R;
        Vector3d T;
        IntegrationBase *pre_integration;
        WheelIntegrationBase *pre_integration_wheel;
        bool is_key_frame;
};
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs);
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &wheel_s);
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x);
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x, Vector3d &wheel_s, Vector3d* Bas);
bool LinearAlignment_NoVelocity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x);
bool LinearAlignment_Stereo(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x);
void WheelExtrisincInitialize(map<double, ImageFrame> &all_image_frame, MatrixXd &r_A, Matrix3d &R0, Matrix3d &rio, Vector3d &tio, VectorXd &x);
void WheelExtrisincInitialize_Road(map<double, ImageFrame> &all_image_frame, MatrixXd &r_A, Matrix3d &R0, Matrix3d &rio, Vector3d &tio, VectorXd &x);
void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x);
