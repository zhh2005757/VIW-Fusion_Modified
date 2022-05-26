/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#include "visualization.h"

using namespace ros;
using namespace Eigen;
ros::Publisher pub_odometry, pub_latest_odometry, pub_latest_wheel_odometry, pub_latest_pure_wheel_odometry;
ros::Publisher pub_path;
ros::Publisher pub_groundtruth;
ros::Publisher pub_wheel_preintegration;
ros::Publisher pub_point_cloud, pub_margin_cloud;
ros::Publisher pub_key_poses;
ros::Publisher pub_camera_pose;
ros::Publisher pub_camera_pose_visual;
nav_msgs::Path path;
nav_msgs::Path groundtruth_path;

ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

ros::Publisher pub_image_track;

CameraPoseVisualization cameraposevisual(1, 0, 0, 1);
static double sum_of_path = 0;
static Vector3d last_path(0.0, 0.0, 0.0);

size_t pub_counter = 0;

void registerPub(ros::NodeHandle &n)
{
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 1000);
    pub_latest_wheel_odometry = n.advertise<nav_msgs::Odometry>("wheel_propagate", 1000);
    pub_latest_pure_wheel_odometry = n.advertise<nav_msgs::Odometry>("pure_wheel_propagate", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);
    pub_groundtruth = n.advertise<nav_msgs::Path>("groundtruth", 1000);
    pub_wheel_preintegration = n.advertise<nav_msgs::Path>("wheel_preintegration", 1000);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    pub_point_cloud = n.advertise<sensor_msgs::PointCloud>("point_cloud", 1000);
    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("margin_cloud", 1000);
    pub_key_poses = n.advertise<visualization_msgs::Marker>("key_poses", 1000);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 1000);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 1000);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 1000);
    pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 1000);
    pub_image_track = n.advertise<sensor_msgs::Image>("image_track", 1000);

    cameraposevisual.setScale(0.1);
    cameraposevisual.setLineWidth(0.01);
}

void pubLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_odometry.publish(odometry);
}
void pubWheelLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_wheel_odometry.publish(odometry);
}
void pubPureWheelLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_pure_wheel_odometry.publish(odometry);
}
void pubTrackImage(const cv::Mat &imgTrack, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgTrackMsg = cv_bridge::CvImage(header, "bgr8", imgTrack).toImageMsg();
    pub_image_track.publish(imgTrackMsg);
}


void printStatistics(const Estimator &estimator, double t)
{
//    if(ESTIMATE_INTRINSIC_WHEEL) {
//        std::ofstream ofs(INTRINSIC_ITERATE_PATH, ios::app);
//        if (!ofs.is_open()) {
//            ROS_WARN("cannot open %s", INTRINSIC_ITERATE_PATH.c_str());
//        }
//        ofs << estimator.sx << " " << estimator.sy << " " << estimator.sw << std::endl;
//    }
//    if(ESTIMATE_EXTRINSIC_WHEEL) {
//        std::ofstream ofs(EXTRINSIC_WHEEL_ITERATE_PATH, ios::app);
//        if (!ofs.is_open()) {
//            ROS_WARN("cannot open %s", EXTRINSIC_WHEEL_ITERATE_PATH.c_str());
//        }
//        ofs << estimator.tio.transpose() <<" "<< Utility::R2ypr(estimator.rio).transpose()<<std::endl;
//    }
//    if(ESTIMATE_EXTRINSIC) {
//        std::ofstream ofs(EXTRINSIC_CAM_ITERATE_PATH, ios::app);
//        if (!ofs.is_open()) {
//            ROS_WARN("cannot open %s", EXTRINSIC_CAM_ITERATE_PATH.c_str());
//        }
//        ofs << estimator.tic[0].transpose() <<" "<< Utility::R2ypr(estimator.ric[0]).transpose()<<std::endl;
//    }

    if (estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    //printf("position: %f, %f, %f\r", estimator.Ps[WINDOW_SIZE].x(), estimator.Ps[WINDOW_SIZE].y(), estimator.Ps[WINDOW_SIZE].z());
    ROS_DEBUG_STREAM("position: " << estimator.Ps[WINDOW_SIZE].transpose());
    ROS_DEBUG_STREAM("orientation: " << estimator.Vs[WINDOW_SIZE].transpose());
    if (ESTIMATE_EXTRINSIC || ESTIMATE_EXTRINSIC_WHEEL || USE_PLANE)
    {
        cv::FileStorage fs(EX_CALIB_RESULT_PATH, cv::FileStorage::WRITE);
        if(ESTIMATE_EXTRINSIC){
            for (int i = 0; i < NUM_OF_CAM; i++)
            {
                //ROS_DEBUG("calibration result for camera %d", i);
                ROS_DEBUG_STREAM("extirnsic tic: " << estimator.tic[i].transpose());
                ROS_DEBUG_STREAM("extrinsic ric: " << Utility::R2ypr(estimator.ric[i]).transpose());

                Eigen::Matrix4d eigen_T = Eigen::Matrix4d::Identity();
                eigen_T.block<3, 3>(0, 0) = estimator.ric[i];
                eigen_T.block<3, 1>(0, 3) = estimator.tic[i];
                cv::Mat cv_T;
                cv::eigen2cv(eigen_T, cv_T);
                if(i == 0)
                    fs << "body_T_cam0" << cv_T ;
                else
                    fs << "body_T_cam1" << cv_T ;
            }
        }

        if(ESTIMATE_EXTRINSIC_WHEEL){
            //ROS_DEBUG("calibration result for camera %d", i);
            ROS_DEBUG_STREAM("extirnsic tio: " << estimator.tio.transpose());
            ROS_WARN_STREAM("extrinsic rio: " << Utility::R2ypr(estimator.rio).transpose());

            Eigen::Matrix4d eigen_T = Eigen::Matrix4d::Identity();
            eigen_T.block<3, 3>(0, 0) = estimator.rio;
            eigen_T.block<3, 1>(0, 3) = estimator.tio;
            cv::Mat cv_T;
            cv::eigen2cv(eigen_T, cv_T);
            fs << "body_T_wheel" << cv_T ;
        }

        if(USE_PLANE){
            ROS_WARN_STREAM("plane zpw: " << estimator.zpw);
            ROS_WARN_STREAM("plane rpw: " << Utility::R2ypr(estimator.rpw).transpose());

            Eigen::Matrix3d eigen_T = estimator.rpw;
            cv::Mat cv_T;
            cv::eigen2cv(eigen_T, cv_T);
            fs << "plane_R_world" << cv_T ;
            fs << "plane_Z_world" <<estimator.zpw;
        }

        fs.release();
    }
    if(ESTIMATE_INTRINSIC_WHEEL){
        cv::FileStorage fs(IN_CALIB_RESULT_PATH, cv::FileStorage::WRITE);

        if(ESTIMATE_INTRINSIC_WHEEL){
            //ROS_DEBUG("calibration result for camera %d", i);
            ROS_WARN_STREAM("intirnsic sx: %f " << " " <<  estimator.sx << " " << "sy: %f " << estimator.sy << " " << "sw: %f" << estimator.sw);
            fs << "sx" << estimator.sx;
            fs << "sy" << estimator.sy;
            fs << "sw" << estimator.sw;
        }


    }

    static double sum_of_time = 0;
    static int sum_of_calculation = 0;
    sum_of_time += t;
    sum_of_calculation++;
    ROS_DEBUG("vo solver costs: %f ms", t);
    ROS_DEBUG("average of time %f ms", sum_of_time / sum_of_calculation);

    sum_of_path += (estimator.Ps[WINDOW_SIZE] - last_path).norm();
    last_path = estimator.Ps[WINDOW_SIZE];
    ROS_DEBUG("sum of path %f", sum_of_path);
    if (ESTIMATE_TD){
        ROS_INFO("td %f", estimator.td);
        std::ofstream ofs(TD_PATH, ios::app);
        if (!ofs.is_open()) {
            ROS_WARN("cannot open %s", TD_PATH.c_str());
        }
        ofs << estimator.td<<std::endl;
    }

    if (ESTIMATE_TD_WHEEL){
        ROS_INFO("td_wheel %f", estimator.td_wheel);
        std::ofstream ofs(TD_WHEEL_PATH, ios::app);
        if (!ofs.is_open()) {
            ROS_WARN("cannot open %s", TD_WHEEL_PATH.c_str());
        }
        ofs << estimator.td_wheel<<std::endl;
    }

}

void pubOdometry(const Estimator &estimator, const std_msgs::Header &header)
{
//    if (estimator.solver_flag == Estimator::SolverFlag::INITIAL)
//    {
//        nav_msgs::Odometry odometry;
//        odometry.header = header;
//        odometry.header.frame_id = "world";
//        odometry.child_frame_id = "world";
//        Quaterniond tmp_Q;
//        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
//        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
//        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
//        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
//        odometry.pose.pose.orientation.x = tmp_Q.x();
//        odometry.pose.pose.orientation.y = tmp_Q.y();
//        odometry.pose.pose.orientation.z = tmp_Q.z();
//        odometry.pose.pose.orientation.w = tmp_Q.w();
//        // odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
//        odometry.twist.twist.linear.x = estimator.Bas[0].x();
//        // odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
//        odometry.twist.twist.linear.y = estimator.Bas[0].y();
//        // odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
//        odometry.twist.twist.linear.z = estimator.Bas[0].z();
//        pub_odometry.publish(odometry);
//    }
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR) //NON_LINEAR
    {
        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();
//        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
//        odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
//        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
//        odometry.twist.twist.linear.x = estimator.Bas[WINDOW_SIZE].x();
//        odometry.twist.twist.linear.y = estimator.Bas[WINDOW_SIZE].y();
//        odometry.twist.twist.linear.z = estimator.Bas[WINDOW_SIZE].z();
        odometry.twist.twist.linear.x = estimator.x_test;
        odometry.twist.twist.linear.y = estimator.y_test;
        odometry.twist.twist.linear.z = estimator.z_test;
//        odometry.twist.twist.linear.x = estimator.yaw_test;
//        odometry.twist.twist.linear.y = estimator.pitch_test;
//        odometry.twist.twist.linear.z = estimator.roll_test;
        pub_odometry.publish(odometry);

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose = odometry.pose.pose;
        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);

//        // write result to file
//        ofstream foutC(VINS_RESULT_PATH, ios::app);
//        foutC.setf(ios::fixed, ios::floatfield);
//        foutC.precision(0);
//        foutC << header.stamp.toSec() * 1e9 << ",";
//        foutC.precision(5);
//        foutC << estimator.Ps[WINDOW_SIZE].x() << ","
//              << estimator.Ps[WINDOW_SIZE].y() << ","
//              << estimator.Ps[WINDOW_SIZE].z() << ","
//              << tmp_Q.w() << ","
//              << tmp_Q.x() << ","
//              << tmp_Q.y() << ","
//              << tmp_Q.z() << ","
//              << estimator.Vs[WINDOW_SIZE].x() << ","
//              << estimator.Vs[WINDOW_SIZE].y() << ","
//              << estimator.Vs[WINDOW_SIZE].z() << "," << endl;
//        foutC.close();
//        Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
//        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
//                                                          tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());

        // write result to file
        ofstream foutC(VINS_RESULT_PATH, ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC << std::setprecision(0)
              << header.stamp.toSec()*1e9<< " "
              << std::setprecision(9)
              << pose_stamped.pose.position.x << " "
              << pose_stamped.pose.position.y << " "
              << pose_stamped.pose.position.z << " "
              << pose_stamped.pose.orientation.x << " "
              << pose_stamped.pose.orientation.y << " "
              << pose_stamped.pose.orientation.z << " "
              << pose_stamped.pose.orientation.w << std::endl;
        foutC.close();
        auto tmp_T = pose_stamped.pose.position;
        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x, tmp_T.y, tmp_T.z,
               tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());

        // write result to file --RIO
        ofstream foutRIO(OUTPUT_FOLDER + "/calib_rio_nonlinear" + to_string(estimator.now) + ".csv", ios::app);
        foutRIO.setf(ios::fixed, ios::floatfield);
        foutRIO << std::setprecision(0)
              << header.stamp.toSec() * 1e9 << ","
              << std::setprecision(9)
              << estimator.yaw_test << ","
              << estimator.pitch_test << ","
              << estimator.roll_test << ","
              << estimator.tmpR_test << ","
              << std::endl;
        foutRIO.close();

        // write result to file --TIO
        ofstream foutTIO(OUTPUT_FOLDER + "/calib_tio_nonlinear" + to_string(estimator.now) + ".csv", ios::app);
        foutTIO.setf(ios::fixed, ios::floatfield);
        foutTIO << std::setprecision(0)
              << header.stamp.toSec() * 1e9 << ","
              << std::setprecision(9)
              << estimator.x_test << ","
              << estimator.y_test << ","
              << estimator.z_test << ","
              << std::endl;
        foutTIO.close();


//        if (estimator.line_start) {
//            if (!estimator.rio_finish) {
//                // write result to file --RIO
//                ofstream foutC(OUTPUT_FOLDER + "/calib_rio_" + to_string(estimator.now) + ".csv", ios::app);
//                foutC.setf(ios::fixed, ios::floatfield);
//                foutC << std::setprecision(0)
//                      << header.stamp.toSec() * 1e9 << ","
//                      << std::setprecision(9)
//                      << estimator.yaw_test << ","
//                      << estimator.pitch_test << ","
//                      << estimator.roll_test << ","
//                      << std::endl;
//                foutC.close();
//            }
//        }else {
//            if (estimator.rio_finish) {
//                // write result to file --TIO
//                ofstream foutC(OUTPUT_FOLDER + "/calib_tio_" + to_string(estimator.now) + ".csv", ios::app);
//                foutC.setf(ios::fixed, ios::floatfield);
//                foutC << std::setprecision(0)
//                      << header.stamp.toSec() * 1e9 << ","
//                      << std::setprecision(9)
//                      << estimator.x_test << ","
//                      << estimator.y_test << ","
//                      << estimator.z_test << ","
//                      << std::endl;
//                foutC.close();
//            }
//        }

    }
}
void pubGroundTruth(Estimator &estimator, const std_msgs::Header &header, Eigen::Matrix<double, 7, 1>& pose, const double td)
{
    estimator.mGTBuf.lock();
    if(estimator.groundtruthBuf.empty()){
        estimator.mGTBuf.unlock();
        return;
    }
    double groundtruth_time = estimator.groundtruthBuf.front().first;
    double header_time = header.stamp.toSec() - OFFSET_SIM;
    pose = estimator.groundtruthBuf.front().second;
    while(groundtruth_time < (header_time - 1e-5)){
        estimator.groundtruthBuf.pop();
        if(estimator.groundtruthBuf.empty())
            break;
        groundtruth_time = estimator.groundtruthBuf.front().first;
        pose = estimator.groundtruthBuf.front().second;
    }
    if(!estimator.groundtruthBuf.empty() && groundtruth_time > (header_time + 1e-5))
    {

        ROS_INFO("wait for new frame");
        ROS_INFO("groundtruth_time: %f, header_time: %f", groundtruth_time, header_time);
        estimator.mGTBuf.unlock();
        return;
    }
    if(estimator.groundtruthBuf.empty() || abs(groundtruth_time - header_time)>1e-5)
    {
        ROS_ERROR("can not find corresponding groundtruth");
    }
    estimator.mGTBuf.unlock();


    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {
        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.orientation.w = pose[0];
        pose_stamped.pose.orientation.x = pose[1];
        pose_stamped.pose.orientation.y = pose[2];
        pose_stamped.pose.orientation.z = pose[3];
        pose_stamped.pose.position.x = pose[4];
        pose_stamped.pose.position.y = pose[5];
        pose_stamped.pose.position.z = pose[6];
        groundtruth_path.header = header;
        groundtruth_path.header.frame_id = "world";
        groundtruth_path.poses.push_back(pose_stamped);
        pub_groundtruth.publish(groundtruth_path);

        // write result to file
        ofstream foutC(GROUNDTRUTH_PATH, ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC << std::setprecision(0)
        << header.stamp.toSec()*1e9<< " "
        << std::setprecision(9)
        << pose_stamped.pose.position.x << " "
        << pose_stamped.pose.position.y << " "
        << pose_stamped.pose.position.z << " "
        << pose_stamped.pose.orientation.x << " "
        << pose_stamped.pose.orientation.y << " "
        << pose_stamped.pose.orientation.z << " "
        << pose_stamped.pose.orientation.w << std::endl;

//        foutC.setf(ios::fixed, ios::floatfield);
//        foutC.precision(0);
//        foutC << header.stamp.toSec() * 1e9 << ",";
//        foutC.precision(5);
//        foutC << estimator.Ps[WINDOW_SIZE].x() << ","
//              << estimator.Ps[WINDOW_SIZE].y() << ","
//              << estimator.Ps[WINDOW_SIZE].z() << ","
//              << tmp_Q.w() << ","
//              << tmp_Q.x() << ","
//              << tmp_Q.y() << ","
//              << tmp_Q.z() << ","
//              << estimator.Vs[WINDOW_SIZE].x() << ","
//              << estimator.Vs[WINDOW_SIZE].y() << ","
//              << estimator.Vs[WINDOW_SIZE].z() << "," << endl;
        foutC.close();
//        Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
//        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
//               tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());

        auto tmp_T = pose_stamped.pose.position;
        auto tmp_Q = pose_stamped.pose.orientation;
        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x, tmp_T.y, tmp_T.z,
               tmp_Q.w, tmp_Q.x, tmp_Q.y, tmp_Q.z);

    }
}
void pubWheelPreintegration(const Eigen::Vector3d& P, const Eigen::Quaterniond& Q,const std_msgs::Header &header)
{
    static nav_msgs::Path preintegration_path;
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header = header;
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose.orientation.w =Q.w();
    pose_stamped.pose.orientation.x = Q.x();
    pose_stamped.pose.orientation.y = Q.y();
    pose_stamped.pose.orientation.z = Q.z();
    pose_stamped.pose.position.x = P[0];
    pose_stamped.pose.position.y = P[1];
    pose_stamped.pose.position.z = P[2];
    preintegration_path.header = header;
    preintegration_path.header.frame_id = "world";
    preintegration_path.poses.push_back(pose_stamped);
    pub_wheel_preintegration.publish(preintegration_path);


}
void pubKeyPoses(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.key_poses.size() == 0)
        return;
    visualization_msgs::Marker key_poses;
    key_poses.header = header;
    key_poses.header.frame_id = "world";
    key_poses.ns = "key_poses";
    key_poses.type = visualization_msgs::Marker::SPHERE_LIST;
    key_poses.action = visualization_msgs::Marker::ADD;
    key_poses.pose.orientation.w = 1.0;
    key_poses.lifetime = ros::Duration();

    //static int key_poses_id = 0;
    key_poses.id = 0; //key_poses_id++;
    key_poses.scale.x = 0.05;
    key_poses.scale.y = 0.05;
    key_poses.scale.z = 0.05;
    key_poses.color.r = 1.0;
    key_poses.color.a = 1.0;

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        geometry_msgs::Point pose_marker;
        Vector3d correct_pose;
        correct_pose = estimator.key_poses[i];
        pose_marker.x = correct_pose.x();
        pose_marker.y = correct_pose.y();
        pose_marker.z = correct_pose.z();
        key_poses.points.push_back(pose_marker);
    }
    pub_key_poses.publish(key_poses);
}

void pubCameraPose(const Estimator &estimator, const std_msgs::Header &header)
{
    int idx2 = WINDOW_SIZE - 1;

    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {
        int i = idx2;
        Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[0]);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();

        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
        if(STEREO)
        {
            Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[1];
            Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[1]);
            cameraposevisual.add_pose(P, R);
        }
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}


void pubPointCloud(const Estimator &estimator, const std_msgs::Header &header)
{
    sensor_msgs::PointCloud point_cloud, loop_point_cloud;
    point_cloud.header = header;
    loop_point_cloud.header = header;


    for (auto &it_per_id : estimator.f_manager.feature)
    {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        if (it_per_id.start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id.solve_flag != 1)
            continue;
        int imu_i = it_per_id.start_frame;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
        Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);
        point_cloud.points.push_back(p);
    }
    ROS_INFO("good point size: %d", point_cloud.points.size());
    pub_point_cloud.publish(point_cloud);


    // pub margined potin
    sensor_msgs::PointCloud margin_cloud;
    margin_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    { 
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        //if (it_per_id->start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id->solve_flag != 1)
        //        continue;

        if (it_per_id.start_frame == 0 && it_per_id.feature_per_frame.size() <= 2 
            && it_per_id.solve_flag == 1 )
        {
            int imu_i = it_per_id.start_frame;
            Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
            Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

            geometry_msgs::Point32 p;
            p.x = w_pts_i(0);
            p.y = w_pts_i(1);
            p.z = w_pts_i(2);
            margin_cloud.points.push_back(p);
        }
    }
    pub_margin_cloud.publish(margin_cloud);
}


void pubTF(const Estimator &estimator, const std_msgs::Header &header)
{
    if( estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion q;
    // body frame
    Vector3d correct_t;
    Quaterniond correct_q;
    correct_t = estimator.Ps[WINDOW_SIZE];
    correct_q = estimator.Rs[WINDOW_SIZE];

    transform.setOrigin(tf::Vector3(correct_t(0),
                                    correct_t(1),
                                    correct_t(2)));
    q.setW(correct_q.w());
    q.setX(correct_q.x());
    q.setY(correct_q.y());
    q.setZ(correct_q.z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "body"));

    // camera frame
    transform.setOrigin(tf::Vector3(estimator.tic[0].x(),
                                    estimator.tic[0].y(),
                                    estimator.tic[0].z()));
    q.setW(Quaterniond(estimator.ric[0]).w());
    q.setX(Quaterniond(estimator.ric[0]).x());
    q.setY(Quaterniond(estimator.ric[0]).y());
    q.setZ(Quaterniond(estimator.ric[0]).z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "camera"));

    // wheel frame
    if(USE_WHEEL && !ONLY_INITIAL_WITH_WHEEL){
        transform.setOrigin(tf::Vector3(estimator.tio.x(),
                                        estimator.tio.y(),
                                        estimator.tio.z()));
        q.setW(Quaterniond(estimator.rio).w());
        q.setX(Quaterniond(estimator.rio).x());
        q.setY(Quaterniond(estimator.rio).y());
        q.setZ(Quaterniond(estimator.rio).z());
        transform.setRotation(q);
        br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "wheel"));
    }


    // plane frame
    if(USE_PLANE){
        transform.setOrigin(tf::Vector3(0,
                                        0,
                                        estimator.zpw));
        q.setW(Quaterniond(estimator.rpw).w());
        q.setX(Quaterniond(estimator.rpw).x());
        q.setY(Quaterniond(estimator.rpw).y());
        q.setZ(Quaterniond(estimator.rpw).z());
        transform.setRotation(q);
        br.sendTransform(tf::StampedTransform(transform, header.stamp, "plane", "world"));
//        std::cout<<"plane rpw: "<< Eigen::Quaterniond(estimator.rpw).coeffs().transpose()<<" zpw: "<<estimator.zpw<<std::endl;
        ROS_DEBUG_STREAM("plane rpw: "<< Eigen::Quaterniond(estimator.rpw).coeffs().transpose()<<" zpw: "<<estimator.zpw);
    }


    nav_msgs::Odometry odometry;
    odometry.header = header;
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = estimator.tic[0].x();
    odometry.pose.pose.position.y = estimator.tic[0].y();
    odometry.pose.pose.position.z = estimator.tic[0].z();
    Quaterniond tmp_q{estimator.ric[0]};
    odometry.pose.pose.orientation.x = tmp_q.x();
    odometry.pose.pose.orientation.y = tmp_q.y();
    odometry.pose.pose.orientation.z = tmp_q.z();
    odometry.pose.pose.orientation.w = tmp_q.w();
    pub_extrinsic.publish(odometry);



}

void pubKeyframe(const Estimator &estimator)
{
    // pub camera pose, 2D-3D points of keyframe
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR && estimator.marginalization_flag == 0)
    {
        int i = WINDOW_SIZE - 2;
        //Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Vector3d P = estimator.Ps[i];
        Quaterniond R = Quaterniond(estimator.Rs[i]);

        nav_msgs::Odometry odometry;
        odometry.header.stamp = ros::Time(estimator.Headers[WINDOW_SIZE - 2]);
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();
        //printf("time: %f t: %f %f %f r: %f %f %f %f\n", odometry.header.stamp.toSec(), P.x(), P.y(), P.z(), R.w(), R.x(), R.y(), R.z());

        pub_keyframe_pose.publish(odometry);


        sensor_msgs::PointCloud point_cloud;
        point_cloud.header.stamp = ros::Time(estimator.Headers[WINDOW_SIZE - 2]);
        point_cloud.header.frame_id = "world";
        for (auto &it_per_id : estimator.f_manager.feature)
        {
            int frame_size = it_per_id.feature_per_frame.size();
            if(it_per_id.start_frame < WINDOW_SIZE - 2 && it_per_id.start_frame + frame_size - 1 >= WINDOW_SIZE - 2 && it_per_id.solve_flag == 1)
            {

                int imu_i = it_per_id.start_frame;
                Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
                Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0])
                                      + estimator.Ps[imu_i];
                geometry_msgs::Point32 p;
                p.x = w_pts_i(0);
                p.y = w_pts_i(1);
                p.z = w_pts_i(2);
                point_cloud.points.push_back(p);

                int imu_j = WINDOW_SIZE - 2 - it_per_id.start_frame;
                sensor_msgs::ChannelFloat32 p_2d;
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.y());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.y());
                p_2d.values.push_back(it_per_id.feature_id);
                point_cloud.channels.push_back(p_2d);
            }

        }
        pub_keyframe_point.publish(point_cloud);
    }
}