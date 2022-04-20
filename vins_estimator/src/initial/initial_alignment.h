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
#include "../factor/imu_factor.h"
//#include "../factor/imu_wheel_factor.h"
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
        // copy structer
//        ImageFrame(const ImageFrame &imageFrame){
//            cout << "ImageFrame(ImageFrame &imageFrame) " << endl;
//            if (nullptr != pre_integration) {
//                delete pre_integration;
//            }
//                pre_integration = new IntegrationBase(*imageFrame.pre_integration);
//            if (nullptr != pre_integration_wheel) {
//                delete pre_integration_wheel;
//                pre_integration_wheel = new WheelIntegrationBase(*imageFrame.pre_integration_wheel);
//            }
//            R = imageFrame.R;
//            T = imageFrame.T;
//            points = imageFrame.points;
//            t = imageFrame.t;
//            is_key_frame = imageFrame.is_key_frame;
//        }
        //copy operator overload
        ImageFrame& operator=(const ImageFrame& imageFrame){
//            cout << "ImageFrame& operator=(const ImageFrame& imageFrame) " << endl;
            if (this != &imageFrame) {
                if (nullptr != imageFrame.pre_integration) {
//                    cout << "imageFrame.pre_integration->acc_0" << imageFrame.pre_integration->acc_0 << endl;
                    pre_integration = new IntegrationBase(imageFrame.pre_integration->acc_0,imageFrame.pre_integration->gyr_0,imageFrame.pre_integration->linearized_ba,imageFrame.pre_integration->linearized_bg);
                    pre_integration->dt_buf = imageFrame.pre_integration->dt_buf;
                    pre_integration->acc_buf = imageFrame.pre_integration->acc_buf;
                    pre_integration->gyr_buf = imageFrame.pre_integration->gyr_buf;
//                    copy = true;
//                    cout << "Constructed___copy    " << copy << endl;
                }
//                if (nullptr != imageFrame.pre_integration_wheel) {
//                    pre_integration_wheel = new WheelIntegrationBase(imageFrame.pre_integration_wheel->vel_0,imageFrame.pre_integration_wheel->gyr_0,imageFrame.pre_integration_wheel->linearized_sx,imageFrame.pre_integration_wheel->linearized_sy,imageFrame.pre_integration_wheel->linearized_sw,imageFrame.pre_integration_wheel->linearized_td);
//                    pre_integration_wheel->dt_buf = imageFrame.pre_integration_wheel->dt_buf;
//                    pre_integration_wheel->vel_buf = imageFrame.pre_integration_wheel->vel_buf;
//                    pre_integration_wheel->gyr_buf = imageFrame.pre_integration_wheel->gyr_buf;
//                }
                R = imageFrame.R;
                T = imageFrame.T;
                points = imageFrame.points;
                t = imageFrame.t;
                is_key_frame = imageFrame.is_key_frame;
            }
            return *this;
        }
//        ~ImageFrame(){
//            cout << "~ImageFrame" << endl;
//            delete pre_integration;
//            pre_integration = nullptr;
//            delete pre_integration_wheel;
//            pre_integration_wheel = nullptr;
//        }

        ImageFrame(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>>& _points, double _t):t{_t},is_key_frame{false}
        {
            points = _points;
        };

//        ~ImageFrame(){
//            cout << "Deconstructed___copy    " << copy << endl;
//            if (copy){
//                delete pre_integration;
//                pre_integration = nullptr;
//                copy = false;
//            }
//        }

        void ImageFrame_Deconstruced(){
            if (nullptr != pre_integration) {
                delete pre_integration;
                pre_integration = nullptr;
            }
        }

        map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>> > > points;
        double t;
        Matrix3d R;
        Vector3d T;
        IntegrationBase *pre_integration = nullptr;
        WheelIntegrationBase *pre_integration_wheel = nullptr;
        bool is_key_frame;
//        bool copy = false;

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
