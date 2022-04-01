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

#include "initial_alignment.h"


void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs)
{
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);
        tmp_A.setZero();
        VectorXd tmp_b(3);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        A += tmp_A.transpose() * tmp_A;
        b += tmp_A.transpose() * tmp_b;
    }
    delta_bg = A.ldlt().solve(b);
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());

    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta_bg;

    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    }

}

/*
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &wheel_s )
{
//    Matrix3d A;
//    Vector3d b;
    MatrixXd A(4, 4);
    VectorXd b(4);
//    Vector3d delta_bg;
    VectorXd delta(4);
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
//        MatrixXd tmp_A(3, 3);
        MatrixXd tmp_A(6, 4);
        tmp_A.setZero();
//        VectorXd tmp_b(3);
        VectorXd tmp_b(6);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        Eigen::Quaterniond qio(RIO);
//        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_A.block<3,3>(0,0) = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_A.block<3,1>(3,3) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(3, 2);
//        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        tmp_b.segment<3>(0) = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        tmp_b.segment<3>(3) = 2 * (frame_j->second.pre_integration_wheel->delta_q.inverse() * qio.inverse() * q_ij * qio).vec();
        A += tmp_A.transpose() * tmp_A;
        b += tmp_A.transpose() * tmp_b;
    }
    delta = A.ldlt().solve(b);
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta.segment<3>(0).transpose());


    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta.segment<3>(0);

    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    }
    double delta_sw = delta(3);
    wheel_s(2) += delta_sw;
    double sw=wheel_s(2);
    ROS_WARN_STREAM("wheel intrinsic sw initial calibration " << sw);
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration_wheel->repropagate(wheel_s(0),wheel_s(1),wheel_s(2));
    }
}
*/

void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &wheel_s )
{
//    Matrix3d A;
//    Vector3d b;
    MatrixXd A(1, 1);
    VectorXd b(1);
//    Vector3d delta_bg;
    VectorXd delta_sw(1);
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
//        MatrixXd tmp_A(3, 3);
        MatrixXd tmp_A(3, 1);
        tmp_A.setZero();
//        VectorXd tmp_b(3);
        VectorXd tmp_b(3);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        Eigen::Quaterniond qio(RIO);
        tmp_A = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(3, 2);
//        tmp_A.block<3,3>(0,0) = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
//        tmp_A.block<3,1>(3,3) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(3, 2);
//        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
//        tmp_b.segment<3>(0) = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        tmp_b = 2 * (frame_j->second.pre_integration_wheel->delta_q.inverse() * qio.inverse() * q_ij * qio).vec();
        A += tmp_A.transpose() * tmp_A;
        b += tmp_A.transpose() * tmp_b;
    }
    delta_sw = A.ldlt().solve(b);
//    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_sw(0));


//    for (int i = 0; i <= WINDOW_SIZE; i++)
//        Bgs[i] += delta.segment<3>(0);
//
//    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
//    {
//        frame_j = next(frame_i);
//        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
//    }
//    double delta_sw = delta(3);
    wheel_s(2) += delta_sw(0);
    double sw=wheel_s(2);
    ROS_WARN_STREAM("wheel intrinsic sw initial calibration " << sw);
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration_wheel->repropagate(wheel_s(0),wheel_s(1),wheel_s(2));
    }
}

MatrixXd TangentBasis(Vector3d &g0)
{
    Vector3d b, c;
    Vector3d a = g0.normalized();
    Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}

void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d* Bas)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 12);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
//            tmp_A.block<3, 3>(0, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
//            tmp_A.block<3, 3>(3, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BA);
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
            b.tail<6>() += r_b.tail<6>();

            A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
            A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
        }
            A = A * 1000.0;
            b = b * 1000.0;
            x = A.ldlt().solve(b);
            VectorXd dg = x.segment<2>(n_state - 6);
            g0 = (g0 + lxly * dg).normalized() * G.norm();
            //double s = x(n_state - 1);
    }   
    g = g0;

//    Bas[0] = x.tail<3>();
//    ROS_WARN_STREAM("Bas " << Bas[0].transpose());
//    for (int i = 0; i <= WINDOW_SIZE; i++)
//        Bas[i] += delta_ba;
//
//    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
//    {
//        frame_j = next(frame_i);
//        frame_j->second.pre_integration->repropagate(delta_ba, Vector3d::Zero());
//    }

}

bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d* Bas)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);

//    double cond = sqrt(1.0/A.ldlt().rcond());
//    ROS_WARN_STREAM("condition num: %f " << cond);

    double s = x(n_state - 1) / 100.0;
    ROS_WARN_STREAM("estimated scale: %f " << s);
    g = x.segment<3>(n_state - 4);
//    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
    {
        return false;
    }

    RefineGravity(all_image_frame, g, x, Bas);
    s = (x.tail<4>())(0) / 100.0;
    (x.tail<4>())(0) = s;
    ROS_WARN_STREAM("refine scale: %f " << s);
//    ROS_WARN_STREAM(" refine  g   " << g.norm() << " " << g.transpose());
    if(s < 0.0 )
        return false;
    else
        return true;
//    return false;
}

void RefineGravity_Stereo(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 12);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
//            tmp_A.block<3, 3>(0, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0 - frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T);

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
//            tmp_A.block<3, 3>(3, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BA);
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
            b.tail<6>() += r_b.tail<6>();

            A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
            A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 6);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;

//    Bas[0] = x.tail<3>();
//    ROS_WARN_STREAM("Bas " << Bas[0].transpose());
//    for (int i = 0; i <= WINDOW_SIZE; i++)
//        Bas[i] += delta_ba;
//
//    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
//    {
//        frame_j = next(frame_i);
//        frame_j->second.pre_integration->repropagate(delta_ba, Vector3d::Zero());
//    }

}

bool LinearAlignment_Stereo(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) ;
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);

//    double cond = sqrt(1.0/A.ldlt().rcond());
//    ROS_WARN_STREAM("condition num: %f " << cond);

//    double s = x(n_state - 1) / 100.0;
//    ROS_WARN_STREAM("estimated scale: %f " << s);
    g = x.segment<3>(n_state - 4);
//    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
//    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
//    {
//        return false;
//    }

    RefineGravity_Stereo(all_image_frame, g, x);
//    s = (x.tail<4>())(0) / 100.0;
//    (x.tail<4>())(0) = s;
//    ROS_WARN_STREAM("refine scale: %f " << s);
//    ROS_WARN_STREAM(" refine  g   " << g.norm() << " " << g.transpose());
//    if(s < 0.0 )
//        return false;
//    else
        return true;
//    return false;
}

void RefineGravity_Joint(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d* Bas)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 12);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            tmp_A.block<3, 3>(0, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 3>(3, 9) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BA);
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            Matrix<double, 6, 6> cov = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
//            cov_inv.setIdentity();

            cov.block<3, 3>(0, 0) = frame_j->second.pre_integration->covariance.block<3, 3>(0, 0);
            cov.block<3, 3>(0, 3) = frame_j->second.pre_integration->covariance.block<3, 3>(0, 6);
            cov.block<3, 3>(3, 0) = frame_j->second.pre_integration->covariance.block<3, 3>(6, 0);
            cov.block<3, 3>(3, 3) = frame_j->second.pre_integration->covariance.block<3, 3>(6, 6);
            cov_inv = cov.inverse();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
            b.tail<6>() += r_b.tail<6>();

            A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
            A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 6);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;

}

bool LinearAlignment_Joint(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d* Bas, Vector3d* Bgs)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 7;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;

    double s=0;

        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++) {
            frame_j = next(frame_i);

            MatrixXd tmp_A(9, 16);
            tmp_A.setZero();
            VectorXd tmp_b(9);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;
            Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);

            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
            tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
//            tmp_A.block<3, 3>(0, 10) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA);
            tmp_A.block<3, 3>(0, 13) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BG);
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];

            tmp_A.block<3, 3>(6, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(6, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 3>(6, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
//            tmp_A.block<3, 3>(6, 10) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BA);
            tmp_A.block<3, 3>(6, 13) = -frame_j->second.pre_integration->jacobian.block<3, 3>(O_V, O_BG);
            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration->delta_v;

            tmp_A.block<3, 3>(3, 13) = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
            tmp_b.block<3, 1>(3, 0) = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();

            Matrix<double, 9, 9> cov_inv = Matrix<double, 9, 9>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
//            cov_inv.setIdentity();
            cov_inv = frame_j->second.pre_integration->covariance.block<9, 9>(0, 0).inverse();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<10, 10>() += r_A.bottomRightCorner<10, 10>();
            b.tail<10>() += r_b.tail<10>();

            A.block<6, 10>(i * 3, n_state - 10) += r_A.topRightCorner<6, 10>();
            A.block<10, 6>(n_state - 10, i * 3) += r_A.bottomLeftCorner<10, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);

//    double cond = sqrt(1.0/A.ldlt().rcond());
//    ROS_WARN_STREAM("condition num: %f " << cond);

        s = x(n_state - 7) / 100.0;
        ROS_WARN_STREAM("estimated scale: %f " << s);
        g = x.segment<3>(n_state - 10);
        ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());

        Vector3d delta_bg = x.tail<3>();
        ROS_WARN_STREAM("delta_bg " << delta_bg.transpose());

        if (fabs(g.norm() - G.norm()) > 0.5 || s < 0) {
            return false;
        }

        RefineGravity_Joint(all_image_frame, g, x, Bas);
        s = (x.tail<4>())(0) / 100.0;
        (x.tail<4>())(0) = s;
        ROS_WARN_STREAM("refine scale: %f " << s);
        ROS_WARN_STREAM(" refine     " << g.norm() << " " << g.transpose());

        Vector3d delta_ba = x.tail<3>();
        ROS_WARN_STREAM("delta_ba " << delta_ba.transpose());

        if (delta_ba.norm() > 0.2)
            return false;

        for (int i = 0; i <= WINDOW_SIZE; i++) {
            Bgs[i] += delta_bg;
        }

        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++) {
            frame_j = next(frame_i);
            frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
        }

//
        for (int i = 0; i <= WINDOW_SIZE; i++) {
            Bas[i] += delta_ba;
        }

        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++) {
            frame_j = next(frame_i);
            frame_j->second.pre_integration->repropagate(Bas[0], Vector3d::Zero());
        }

    if(s < 0.1 )
        return false;
    else {
        return true;
    }
//    return false;
}

void RefineGravity_Sv(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = 6 * (all_frame_count - 1);

    MatrixXd A{6, 4};
    A.setZero();
    VectorXd b{6};
    b.setZero();
    MatrixXd r_A{n_state, 4};
    r_A.setZero();
    VectorXd r_b{n_state};
    r_b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 4);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


//            tmp_A.block<3, 1>(0, 0) = -RIO * frame_j->second.pre_integration_wheel->vel_0 * dt;
//            tmp_A.block<3, 2>(0, 1) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
//            tmp_A.block<3, 1>(0, 3) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) ;
//            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 1>(3, 0) = -RIO * frame_j->second.pre_integration_wheel->vel_0 + frame_i->second.R.transpose() * frame_j->second.R * RIO * frame_j->second.pre_integration_wheel->vel_1;
            tmp_A.block<3, 2>(3, 1) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0 ;


//            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
//            cov_inv.setIdentity();

            r_A.block<6, 4>(i * 6, 0) = tmp_A;
            r_b.segment<6>(i * 6) = tmp_b;

        }
        A = r_A.transpose() * r_A;
        A = A * 1000.0;
        b = r_A.transpose() * r_b;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(1);
        g0 = (g0 + lxly * dg).normalized() * G.norm();

//        double s = x(3) / 100.0 ;
//        cout << "v is :" << frame_j->second.pre_integration->delta_v.transpose() << s << endl;
    }
    g = g0;
}

bool LinearAlignment_Sv(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d &wheel_s)
{
    int all_frame_count = all_image_frame.size();
    int n_state = 6 * (all_frame_count - 1);

    MatrixXd A{5, 5};
    A.setZero();
    VectorXd b{5};
    b.setZero();
    MatrixXd r_A{n_state, 5};
    r_A.setZero();
    VectorXd r_b{n_state};
    r_b.setZero();
    MatrixXd cov_inv{n_state, n_state};
    cov_inv.setIdentity();
    VectorXd w{all_frame_count-1};
    w.setOnes();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k=0;k<30;k++) {
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++) {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 5);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;

            tmp_A.block<3, 1>(0, 0) = -RIO * frame_j->second.pre_integration_wheel->vel_0 * dt;
            tmp_A.block<3, 3>(0, 1) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
            tmp_A.block<3, 1>(0, 4) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T);
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p +
                                      frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
            //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
            tmp_A.block<3, 1>(3, 0) = -RIO * frame_j->second.pre_integration_wheel->vel_0 +
                                      frame_i->second.R.transpose() * frame_j->second.R * RIO *
                                      frame_j->second.pre_integration_wheel->vel_1;
            tmp_A.block<3, 3>(3, 1) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
            //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

//        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
//        cov_inv.setIdentity();

            r_A.block<6, 5>(i * 6, 0) = w(i) * tmp_A;
            r_b.segment<6>(i * 6) = w(i) * tmp_b;

//        cout << "v is :" << frame_j->second.pre_integration->delta_v.transpose() <<  endl;
//        cout << "p is :" << frame_j->second.pre_integration->delta_p.transpose() <<  endl;
//        cout << "p1 is :" << frame_j->second.pre_integration_wheel->delta_p.transpose() <<  endl;
//        cout << "RIO is :" << RIO << endl;
//        cout << "frame.T :" << frame_i->second.T.transpose() << endl;

        }
        A = r_A.transpose() * r_A;
        A = A * 1000.0;
        b = r_A.transpose() * r_b;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        cout << "s is :" << x(4) <<  " sx is :" << x(0) << endl;
        for (i = 0; i < all_frame_count - 1; i++) {
            w(i) = exp(-fabs((r_A * x - r_b)(i)) * 100);
        }
    }
//    for (int i=0;i<30;i++) {
//        A = r_A.transpose() * cov_inv * r_A;
//        A = A * 1000.0;
//        b = r_A.transpose() * cov_inv * r_b;
//        b = b * 1000.0;
//        x = A.ldlt().solve(b);
//        cout << "s is :" << x(4) <<  " sx is :" << x(0) << endl;
//        cov_inv = ((r_A * x - r_b) * (r_A * x - r_b).transpose()).inverse();
//    }

//    cout << "residual : " << endl << r_A * x - r_b << endl;

    double cond = sqrt(1.0/A.ldlt().rcond());
    ROS_WARN_STREAM("condition num: %f " << cond);

    double s = x(4) ;

    if (s < 0.0) {}
    else {
        ROS_WARN_STREAM("estimated scale: %f " << s);
        g = x.segment<3>(1);
        ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
        if (fabs(g.norm() - G.norm()) > 0.5 || s < 0) {
            return false;
        }
        double sx = wheel_s(0);
        double sy = wheel_s(1);
        double sw = wheel_s(2);
        sx = x(0);
        if (sx < 0)
            return false;
        sy = sx;
        ROS_WARN_STREAM("sx: %f " << sx << " " << "sy: %f " << sy);
    }

    RefineGravity_Sv(all_image_frame, g, x);
    s = (x.tail<1>())(0) ;
    (x.tail<1>())(0) = s;
    if (s < 0.0) {}
    else {
        ROS_WARN_STREAM("refine scale: %f " << s);
        ROS_WARN_STREAM(" refine     " << g.norm() << " " << g.transpose());

        double sx = wheel_s(0);
        double sy = wheel_s(1);
        double sw = wheel_s(2);
        sx = x(0);
        if (sx < 0)
            return false;
        sy = sx;
        ROS_WARN_STREAM("sx: %f " << sx << " " << "sy: %f " << sy);

        wheel_s(0) = sx;
        wheel_s(1) = sy;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++) {
            frame_j = next(frame_i);
            frame_j->second.pre_integration_wheel->repropagate(sx, sy, sw);
        }
    }
//    if(s < 0.0 )
//        return false;
//    else
//        return true;
    return false;

}

// Stereo
void RefineGravity_NoVelocity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = 3 * (all_frame_count-2);

    MatrixXd A{3, 3};
    A.setZero();
    VectorXd b{3};
    b.setZero();
    MatrixXd r_A{n_state, 3};
    r_A.setZero();
    VectorXd r_b{n_state};
    r_b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    map<double, ImageFrame>::iterator frame_k;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(next(frame_i)) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);
            frame_k = next(frame_j);

            MatrixXd tmp_A(3, 3);
            tmp_A.setZero();
            VectorXd tmp_b(3);
            tmp_b.setZero();

            double dt12 = frame_j->second.pre_integration->sum_dt;
            double dt23 = frame_k->second.pre_integration->sum_dt;

            tmp_A.block<3, 2>(0, 0) = 0.5 * frame_i->second.R.transpose() * (dt23 * dt23 * dt12 + dt12 * dt12 * dt23) * lxly;
//            tmp_A.block<3, 1>(0, 2) = frame_i->second.R.transpose() * ((frame_k->second.T - frame_j->second.T) * dt12 - (frame_j->second.T - frame_i->second.T) * dt23) ;
            tmp_b.block<3, 1>(0, 0) = frame_i->second.R.transpose() * frame_j->second.R * frame_k->second.pre_integration->delta_p * dt12 - frame_j->second.pre_integration->delta_p * dt23 + frame_j->second.pre_integration->delta_v * dt12 * dt23 - 0.5 * frame_i->second.R.transpose() * (dt23 * dt23 * dt12 + dt12 * dt12 * dt23) * g0 - frame_i->second.R.transpose() * (frame_j->second.R - frame_k->second.R) * TIC[0] * dt12 + frame_i->second.R.transpose() * (frame_i->second.R - frame_j->second.R) * TIC[0] * dt23 - frame_i->second.R.transpose() * ((frame_k->second.T - frame_j->second.T) * dt12 - (frame_j->second.T - frame_i->second.T) * dt23);

            r_A.block<3, 3>(i * 3, 0) = tmp_A;
            r_b.segment<3>(i * 3) = tmp_b;

        }
        A = r_A.transpose() * r_A;
        A = A * 1000.0;
        b = r_A.transpose() * r_b;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(0);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignment_NoVelocity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    int n_state = 6 * (all_frame_count - 2);

    MatrixXd A{5, 5};
    A.setZero();
    VectorXd b{5};
    b.setZero();
    MatrixXd r_A{n_state, 5};
    r_A.setZero();
    VectorXd r_b{n_state};
    r_b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    map<double, ImageFrame>::iterator frame_k;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(next(frame_i)) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);
        frame_k = next(frame_j);

        MatrixXd tmp_A(6, 5);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt12 = frame_j->second.pre_integration->sum_dt;
        double dt23 = frame_k->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = 0.5 * frame_i->second.R.transpose() * (dt23 * dt23 * dt12 + dt12 * dt12 * dt23);
//        tmp_A.block<3, 1>(0, 3) = frame_i->second.R.transpose() * ((frame_k->second.T - frame_j->second.T) * dt12 - (frame_j->second.T - frame_i->second.T) * dt23) / 100 ;
        tmp_b.block<3, 1>(0, 0) = frame_i->second.R.transpose() * frame_j->second.R * frame_k->second.pre_integration->delta_p * dt12 - frame_j->second.pre_integration->delta_p * dt23 + frame_j->second.pre_integration->delta_v * dt12 * dt23 - frame_i->second.R.transpose() * (frame_j->second.R - frame_k->second.R) * TIC[0] * dt12 + frame_i->second.R.transpose() * (frame_i->second.R - frame_j->second.R) * TIC[0] * dt23 - frame_i->second.R.transpose() * ((frame_k->second.T - frame_j->second.T) * dt12 - (frame_j->second.T - frame_i->second.T) * dt23);

        //wheel preintegration w.r.t s,sx
//        tmp_A.block<3, 1>(3, 3) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
//        tmp_A.block<3, 1>(3, 4) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
//        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);


        r_A.block<6, 5>(i * 6, 0) = tmp_A;
        r_b.segment<6>(i * 6) = tmp_b;

    }
//    ROS_WARN_STREAM("i: " << i);
    A = r_A.transpose() * r_A;
    A = A * 1000.0;
//    ROS_WARN_STREAM("A size:  " << A);
    b = r_A.transpose() * r_b;
    b = b * 1000.0;
//    ROS_WARN_STREAM("b size:  " << b);
    x = A.ldlt().solve(b);

    double cond = sqrt(1.0/A.ldlt().rcond());
    ROS_WARN_STREAM("condition num: %f " << cond);

//    double s = x(3) / 100 ;
//    ROS_WARN_STREAM("estimated scale: %f " << s);
//    double delta_sx=x(4);
//    ROS_WARN_STREAM("delta_sx: %f " << delta_sx);
    g = x.segment<3>(0);
    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
//    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
//    {
//        return false;
//    }

    RefineGravity_NoVelocity(all_image_frame, g, x);
//    s = (x.tail<1>())(0) ;
//    (x.tail<1>())(0) = s;
//    ROS_WARN_STREAM("refine scale: %f " << s);
    ROS_WARN_STREAM(" refine g    " << g.norm() << " " << g.transpose());
//    if(s < 0.0 )
//        return false;
//    else
//        return true;
//    return true;

}

//original wheel 

/*
void RefineGravityWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(9, 9);
            tmp_A.setZero();
            VectorXd tmp_b(9);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

            tmp_A.block<3,1>(6, 8) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);

            Matrix<double, 9, 9> cov_inv = Matrix<double, 9, 9>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 3);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignmentWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(9, 10);
        tmp_A.setZero();
        VectorXd tmp_b(9);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;

        tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 9, 9> cov_inv = Matrix<double, 9, 9>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
    double s = x(n_state - 1) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    ROS_WARN_STREAM("estimated scale: %f " << s);
    g = x.segment<3>(n_state - 4);
    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
    {
        return false;
    }

    RefineGravityWithWheel(all_image_frame, g, x);
    s = (x.tail<1>())(0) / 100.0;
    (x.tail<1>())(0) = s;
    ROS_WARN_STREAM("refine scale: %f " << s);
    ROS_WARN_STREAM(" refine     " << g.norm() << " " << g.transpose());
    return s >= 0.0;
}
*/

// modified with respect to joint estimation （ Sx Sy Sw)

/*
void RefineGravityWithIntrWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 2 + 1;
    int n_state = all_frame_count * 3 + 2 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

//            MatrixXd tmp_A(9, 9);
            MatrixXd tmp_A(12, 12);
            tmp_A.setZero();
            VectorXd tmp_b(12);
            tmp_b.setZero();

            Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
            Eigen::Quaterniond qio(RIO);

            double dt = frame_j->second.pre_integration->sum_dt;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

            tmp_A.block<3, 1>(6, 8) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
            tmp_A.block<3, 1>(6, 9) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
            // tmp_A.block<3, 1>(6, 10) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);
            // tmp_A.block<3, 1>(6, 11) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 2);
            // tmp_A.block<3, 1>(9,11) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(3, 2);
//            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
            // tmp_b.block<3, 1>(9,0) = 2 * (frame_j->second.pre_integration_wheel->delta_q.inverse() * qio.inverse() * q_ij * qio).vec();
            Matrix<double, 12, 12> cov_inv = Matrix<double, 12, 12>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

//            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
//            b.tail<3>() += r_b.tail<3>();

            A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
            b.tail<6>() += r_b.tail<6>();

//            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
//            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();

            A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
            A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 6);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignmentWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d &wheel_s)
{
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 3 + 1;
    int n_state = all_frame_count * 3 + 3 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

//        MatrixXd tmp_A(9, 10);
        MatrixXd tmp_A(12, 13);
        tmp_A.setZero();
        VectorXd tmp_b(12);
        tmp_b.setZero();

        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        Eigen::Quaterniond qio(RIO);

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;

//        tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
//        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
        tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
        // tmp_A.block<3, 1>(6, 10) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
        // tmp_A.block<3, 1>(6, 11) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);
        tmp_A.block<3, 1>(6, 12) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 2);
        tmp_A.block<3, 1>(9,12) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(3, 2);
        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p  - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
        tmp_b.block<3, 1>(9,0) = 2 * (frame_j->second.pre_integration_wheel->delta_q.inverse() * qio.inverse() * q_ij * qio).vec();
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 12, 12> cov_inv = Matrix<double, 12, 12>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

//        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
//        b.tail<4>() += r_b.tail<4>();

        A.bottomRightCorner<7, 7>() += r_A.bottomRightCorner<7, 7>();
        b.tail<7>() += r_b.tail<7>();

//        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
//        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();

        A.block<6, 7>(i * 3, n_state - 7) += r_A.topRightCorner<6, 7>();
        A.block<7, 6>(n_state - 7, i * 3) += r_A.bottomLeftCorner<7, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
//    double s = x(n_state - 1) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    double s = x(n_state - 4) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    ROS_WARN_STREAM("estimated scale: %f " << s);
//    g = x.segment<3>(n_state - 4);
    g = x.segment<3>(n_state - 7);
    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
    {
        return false;
    }

    RefineGravityWithIntrWheel(all_image_frame, g, x);
    s = (x.tail<4>())(0) / 100.0;
    (x.tail<4>())(0) = s;
    ROS_WARN_STREAM("refine scale: %f " << s);
    ROS_WARN_STREAM(" refine     " << g.norm() << " " << g.transpose());

    double sx=wheel_s(0);
    double sy=wheel_s(1);
    double sw=wheel_s(2);
    sx += (x.tail<3>())(0);
    if (sx < 0)
        return false;
   sy=sx;
   sw += (x.tail<3>())(2);
   wheel_s(0)=sx;
    wheel_s(1)=sy;
    wheel_s(2)=sw;
    ROS_WARN_STREAM("sx: %f " <<  sx << " "  << "sy: %f " << sy << " " << "sw: %f " << sw);
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration_wheel->repropagate(sx,sy,sw);
    }

    return s >= 0.0;

}
*/

// modified with respect to disjoint estimation (mono)

/*
void RefineGravityWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 2 + 1;
    int n_state = all_frame_count * 3 + 2 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

//            MatrixXd tmp_A(9, 9);
            MatrixXd tmp_A(12, 12);
            tmp_A.setZero();
            VectorXd tmp_b(12);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;

//            double s= (x.tail<3>())(0) / 100.0;

            // visual and IMU
            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;
//            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0-s*frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) ;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

             // wheel and visual
//             tmp_A.block<3, 1>(6, 8) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
////            tmp_A.block<3, 1>(6, 9) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
////             tmp_A.block<3, 1>(6, 10) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);
//             tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
////            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO) -s* (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) ;

              // wheel and imu
            tmp_A.block<3, 3>(9, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(9, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;

            // pbo
            tmp_A.block<3, 3>(9, 9) = Matrix3d::Identity() - frame_i->second.R.transpose() * frame_j->second.R;

            // sx sy
//            tmp_A.block<3, 1>(9, 9) = RIO * frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
//            tmp_A.block<3, 1>(9, 10) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);

            // pbo
            tmp_b.block<3, 1>(9, 0) = frame_j->second.pre_integration->delta_p - RIO * frame_j->second.pre_integration_wheel->delta_p;

//            tmp_b.block<3, 1>(9, 0) = frame_j->second.pre_integration->delta_p - RIO * frame_j->second.pre_integration_wheel->delta_p - TIO + frame_i->second.R.transpose() * frame_j->second.R * TIO- frame_i->second.R.transpose() * dt * dt / 2 * g0;

            Matrix<double, 12, 12> cov_inv = Matrix<double, 12, 12>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();


            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

//            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
//            b.tail<3>() += r_b.tail<3>();

            A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
            b.tail<6>() += r_b.tail<6>();

//            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
//            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();

            A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
            A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 6);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignmentWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d &wheel_s)
{
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 3 + 1;
    int n_state = all_frame_count * 3 + 3 + 4;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

//        MatrixXd tmp_A(9, 10);
        MatrixXd tmp_A(12, 13);
        tmp_A.setZero();
        VectorXd tmp_b(12);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        //IMU P Q preintegration w.r.t s,v,g
        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];

        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;




        //wheel preintegration w.r.t s,sx
//        tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
////        tmp_A.block<3, 1>(6, 10) = n*(-frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0));
//        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);




//        //wheel IMU P preintegration w.r.t sx,v,g
        tmp_A.block<3, 3>(9, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(9, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();

        // estimate Pbo
        tmp_A.block<3, 3>(9, 10) = Matrix3d::Identity() - frame_i->second.R.transpose() * frame_j->second.R;

        // sx sy
//        tmp_A.block<3, 1>(9, 10) = RIO * frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
//        tmp_A.block<3, 1>(9, 11) = frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);

        // estimate Pbo
        tmp_b.block<3, 1>(9, 0) = frame_j->second.pre_integration->delta_p - RIO * frame_j->second.pre_integration_wheel->delta_p;
//        tmp_b.block<3, 1>(9, 0) = frame_j->second.pre_integration->delta_p - RIO * frame_j->second.pre_integration_wheel->delta_p - TIO + frame_i->second.R.transpose() * frame_j->second.R * TIO;

        Matrix<double, 12, 12> cov_inv = Matrix<double, 12, 12>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        //TODO think about the covariance matrix
        cov_inv.setIdentity();


//        cout << "wheel Jacabian: " << tmp_A.block<3, 1>(6, 10).transpose() << endl;
//        cout <<"wheel delta_q: " << frame_j->second.pre_integration_wheel->delta_q.toRotationMatrix() << endl;
//        cout << "ba Jacabian: " << frame_j->second.pre_integration->jacobian.block<3, 3>(O_P, O_BA).transpose();
//        cout << "delta_q & frame i j: " << frame_j->second.pre_integration->delta_q << " " << frame_i->second.R.transpose() * frame_j->second.R << endl;


        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

//        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
//        b.tail<4>() += r_b.tail<4>();

        A.bottomRightCorner<7, 7>() += r_A.bottomRightCorner<7, 7>();
        b.tail<7>() += r_b.tail<7>();

//        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
//        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();

        A.block<6, 7>(i * 3, n_state - 7) += r_A.topRightCorner<6, 7>();
        A.block<7, 6>(n_state - 7, i * 3) += r_A.bottomLeftCorner<7, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);


    double cond = sqrt(1.0/A.ldlt().rcond());
    ROS_WARN_STREAM("condition num: %f " << cond);
//    double s = x(n_state - 1) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    double s = x(n_state - 4) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    ROS_WARN_STREAM("estimated scale: %f " << s);
//    g = x.segment<3>(n_state - 4);
    g = x.segment<3>(n_state - 7);
    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
//    if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
//    {
//        return false;
//    }

    Vector3d pbo = x.tail<3>();
    ROS_WARN_STREAM(" result pbo     " << pbo.transpose());

//    double sx=wheel_s(0);
//    double sy=wheel_s(1);
//    double sw=wheel_s(2);
//    sx += (x.tail<3>())(1);
//    if (sx < 0)
//        return false;
//    sy=sx;
////     sy += (x.tail<3>())(2);
////     if (sy < 0)
////         return false;
//    ROS_WARN_STREAM("sx: %f " << sx << " " << "sy: %f " << sy);
//
//
//    Vector3d tmp_V=x.segment<3>(0);
//    ROS_WARN_STREAM("tmp velocity :" << tmp_V.transpose());


    RefineGravityWithWheel(all_image_frame, g, x);
     s = (x.tail<4>())(0) / 100.0;
     (x.tail<4>())(0) = s;
     ROS_WARN_STREAM("refine scale: %f " << s);
    ROS_WARN_STREAM(" refine g    " << g.norm() << " " << g.transpose());
    pbo = x.tail<3>();
    ROS_WARN_STREAM(" refine pbo     " << pbo.transpose());
//    sx=wheel_s(0);
//    sy=wheel_s(1);
//    sw=wheel_s(2);
//    sx += (x.tail<3>())(1);
//    sy=sx;
////    sy += (x.tail<3>())(2);
//    ROS_WARN_STREAM("refine sx: %f " << sx << " " << "sy: %f " << sy);
//
//    tmp_V=x.segment<3>(0);
//    ROS_WARN_STREAM("refine tmp velocity :" << tmp_V.transpose());

//    wheel_s(0)=sx;
//    wheel_s(1)=sy;
//    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
//    {
//        frame_j = next(frame_i);
//        frame_j->second.pre_integration_wheel->repropagate(sx,sy,sw);
//    }

    return s >= 0.1;
//    return false;
}
*/

// modified with respect to disjoint estimation (stereo) 
/*
 void m_RefineGravityWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 2 + 1;
    int n_state = all_frame_count * 3 + 2 + 3;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

//            MatrixXd tmp_A(9, 9);
            MatrixXd tmp_A(9, 11);
            tmp_A.setZero();
            VectorXd tmp_b(9);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;

            // double s= (x.tail<3>())(0) / 100.0;


            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            // tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            // tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0-frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) ;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

            // tmp_A.block<3, 1>(6, 8) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
            tmp_A.block<3, 1>(6, 9) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
            // tmp_A.block<3, 1>(6, 10) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);
//            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
            tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO) - (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) ;
            // + frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 2) * delta_sw
            Matrix<double, 9, 9> cov_inv = Matrix<double, 9, 9>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            // ROS_WARN_STREAM("tmp_A is : \n" << tmp_A);
            // ROS_WARN_STREAM("tmp_b is : \n" << tmp_b.transpose());


            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

//            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
//            b.tail<3>() += r_b.tail<3>();

            A.bottomRightCorner<5, 5>() += r_A.bottomRightCorner<5, 5>();
            b.tail<5>() += r_b.tail<5>();

//            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
//            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();

            A.block<6, 5>(i * 3, n_state - 5) += r_A.topRightCorner<6, 5>();
            A.block<5, 6>(n_state - 5, i * 3) += r_A.bottomLeftCorner<5, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 5);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignmentWithWheel(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d &wheel_s)
{
    int all_frame_count = all_image_frame.size();
//    int n_state = all_frame_count * 3 + 3 + 1;
    int n_state = all_frame_count * 3 + 3 + 3;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

//        MatrixXd tmp_A(9, 10);
        MatrixXd tmp_A(9, 12);
        tmp_A.setZero();
        VectorXd tmp_b(9);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        // tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
        // tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
         tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0]-frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T);
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;

//        tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
//        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO);
        // tmp_A.block<3, 1>(6, 9) = (frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T) / 100;
        tmp_A.block<3, 1>(6, 10) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 0);
        // tmp_A.block<3, 1>(6, 11) = -frame_j->second.pre_integration_wheel->jacobian.template block<3, 1>(0, 1);
        tmp_b.block<3, 1>(6, 0) = frame_j->second.pre_integration_wheel->delta_p - RIO.transpose() * frame_i->second.R.transpose() * frame_j->second.R * TIO + (frame_i->second.R * RIO).transpose() * frame_j->second.R * TIC[0] - RIO.transpose() * (TIC[0] - TIO)-(frame_i->second.R * RIO).transpose() * (frame_j->second.T - frame_i->second.T);
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 9, 9> cov_inv = Matrix<double, 9, 9>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

//        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
//        b.tail<4>() += r_b.tail<4>();

        A.bottomRightCorner<6, 6>() += r_A.bottomRightCorner<6, 6>();
        b.tail<6>() += r_b.tail<6>();

//        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
//        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();

        A.block<6, 6>(i * 3, n_state - 6) += r_A.topRightCorner<6, 6>();
        A.block<6, 6>(n_state - 6, i * 3) += r_A.bottomLeftCorner<6, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
//    double s = x(n_state - 1) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    // double s = x(n_state - 3) / 100.0;//TODO:尺度因子处理，除以了100,（PS：100在tmp_A.block<3, 1>(0, 9)有所体现，但具体为何如此处理还未知
    // ROS_WARN_STREAM("estimated scale: %f " << s);
//    g = x.segment<3>(n_state - 4);
    g = x.segment<3>(n_state - 6);
    ROS_WARN_STREAM(" result g     " << g.norm() << " " << g.transpose());
    // if(fabs(g.norm() - G.norm()) > 0.5 || s < 0)
    // {
    //     return false;
    // }

    m_RefineGravityWithWheel(all_image_frame, g, x);
    // s = (x.tail<3>())(0) / 100.0;
    // (x.tail<3>())(0) = s;
    // ROS_WARN_STREAM("refine scale: %f " << s);
    ROS_WARN_STREAM(" refine     " << g.norm() << " " << g.transpose());

    double sx=wheel_s(0);
    double sy=wheel_s(1);
    double sw=wheel_s(2);
    sx += (x.tail<3>())(1);
    if (sx < 0 || sx > 2) 
        return false;
   sy=sx;
    // sy += (x.tail<3>())(2);
    // if (sy < 0)
        // return false;
    ROS_WARN_STREAM("sx: %f " << sx << "sy: %f " << sy);
    wheel_s(0)=sx;
    wheel_s(1)=sy;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration_wheel->repropagate(sx,sy,sw);
    }

    // return s >= 0.0;
    return true;

}
*/

/*
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x)
{
    solveGyroscopeBias(all_image_frame, Bgs);
    if(USE_WHEEL)
        return LinearAlignmentWithWheel(all_image_frame, g, x);
    else
        return LinearAlignment(all_image_frame, g, x);
}
*/

void WheelExtrisincInitialize(map<double, ImageFrame> &all_image_frame, MatrixXd &r_A, Matrix3d &R0, Matrix3d &rio, Vector3d &tio, VectorXd &x){

    MatrixXd A{(all_image_frame.size() ) * 3, 6};
    A.setZero();
    VectorXd b{(all_image_frame.size() )* 3};
    b.setZero();
    VectorXd err{all_image_frame.size() * 3 };
    err.setZero();
    VectorXd x_ep;
    for (int n=0;n<1;n++) {
        Vector3d        Ps[(all_image_frame.size())];
        Vector3d        Vs[(all_image_frame.size())];
        Matrix3d        Rs[(all_image_frame.size())];
        map<double, ImageFrame>::iterator frame_i;
        map<double, ImageFrame>::iterator frame_j;

        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++ ) {
            frame_j = next(frame_i);
            Vs[i] = R0 * frame_i->second.R * x.segment<3>(i * 3);
            Rs[i] = R0 * frame_i->second.R;
            Matrix3d R_w_x;
            Vector3d w_x = R0 * RIC[0].transpose() * (frame_j->second.pre_integration->linearized_gyr - frame_j->second.pre_integration->linearized_bg);
            R_w_x << 0, -w_x(2), w_x(1),
                    w_x(2), 0, -w_x(0),
                    -w_x(1), w_x(0), 0;
            A.block<3, 3>(i * 3, 0) = exp(-err.segment<3>(i * 3).norm() * 100) * RIC[0] * R0.transpose();
            A.block<3, 3>(i * 3, 3) = exp(-err.segment<3>(i * 3).norm() * 100) * RIC[0] * R0.transpose() * R_w_x;
            b.segment<3>(i * 3) = exp(-err.segment<3>(i * 3).norm() * 100) * Rs[i].transpose() * Vs[i];

        }
        {
            Vs[all_image_frame.size()-1] = R0 * frame_j->second.R * x.segment<3>((all_image_frame.size()-1) * 3);
            Rs[all_image_frame.size()-1] = R0 * frame_j->second.R;
            Matrix3d R_w_x;
            Vector3d w_x = R0 * RIC[0].transpose() * (frame_j->second.pre_integration->gyr_1 - frame_j->second.pre_integration->linearized_bg);
            R_w_x << 0, -w_x(2), w_x(1),
                    w_x(2), 0, -w_x(0),
                    -w_x(1), w_x(0), 0;
            A.block<3, 3>((all_image_frame.size()-1) * 3, 0) = exp(-err.segment<3>((all_image_frame.size()-1) * 3).norm() * 100) * RIC[0] * R0.transpose();
            A.block<3, 3>((all_image_frame.size()-1) * 3, 3) = exp(-err.segment<3>((all_image_frame.size()-1) * 3).norm() * 100) * RIC[0] * R0.transpose() * R_w_x;
            b.segment<3>((all_image_frame.size()-1) * 3) = exp(-err.segment<3>((all_image_frame.size()-1) * 3).norm() * 100) * Rs[all_image_frame.size()-1].transpose() * Vs[all_image_frame.size()-1];

        }
        r_A = A.transpose() * A;
        VectorXd r_b = A.transpose() * b;
        x_ep = r_A.ldlt().solve(r_b);
        err = (A * x_ep - b).cwiseAbs();
//            ROS_WARN_STREAM("err     " << err);
        ROS_WARN_STREAM("max err     " << err.maxCoeff());
    }
    r_A = r_A.inverse();
//        for (int i = 0; i <= frame_count; i++)
//        {
//            Vb[i] = x_ep.head<3>();
//        }

    Matrix3d tmp_R = Eigen::Quaterniond::FromTwoVectors(x_ep.head<3>(), Vector3d{1,0,0}).toRotationMatrix();
    rio = RIC[0] * R0.transpose() * Utility::ypr2R(Eigen::Vector3d{Utility::R2ypr(tmp_R).x() , 0, 0}).transpose();
    tio = x_ep.tail<3>();
    ROS_WARN_STREAM("rio     " << rio);
    ROS_WARN_STREAM("rio ypr     " << Utility::R2ypr_m(rio).transpose());
    ROS_WARN_STREAM("tio     " << tio.transpose());

}

bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x, Vector3d &wheel_s, Vector3d* Bas)
{
    solveGyroscopeBias(all_image_frame, Bgs);
//    solveGyroscopeBias(all_image_frame, Bgs, wheel_s);
//    if(USE_WHEEL) {
////        solveGyroscopeBias(all_image_frame, Bgs, wheel_s);
////        for(int i=0;i<20;i++){
////            if(!LinearAlignmentWithWheel(all_image_frame, g, x,wheel_s))
////                return false;
////        }
////        return true;
//        return LinearAlignment_NoVelocity(all_image_frame, g, x);
////        return LinearAlignmentWithWheel(all_image_frame, g, x);
//        return LinearAlignmentWithWheel(all_image_frame, g, x,wheel_s);
////        return LinearAlignment_NoVelocity(all_image_frame, g, x);
////        return LinearAlignment_Sv(all_image_frame, g, x, wheel_s);
//    }
//    else{
//        return LinearAlignment_Joint(all_image_frame, g, x, Bas, Bgs);
        return LinearAlignment(all_image_frame, g, x, Bas);
//    }

}
