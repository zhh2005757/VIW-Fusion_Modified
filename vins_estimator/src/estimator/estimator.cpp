/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *******************************************************/

#include "estimator.h"
#include "../utility/visualization.h"
#include "../factor/pose_subset_parameterization.h"
#include "../factor/orientation_subset_parameterization.h"
#include "../factor/Rwg_SO3_local_parameterization.h"
#include "../factor/lidarFactor.hpp"
#define DISTORTION 1

constexpr double SCAN_PERIOD = 0.1;
constexpr double DISTANCE_SQ_THRESHOLD = 25;
constexpr double NEARBY_SCAN = 2.5;

double para_q[4] = {0, 0, 0, 1};
double para_t[3] = {0, 0, 0};
Eigen::Map<Eigen::Quaterniond> q_last_curr(para_q); //TODO Use IMU to deskew (interpolation)
Eigen::Map<Eigen::Vector3d> t_last_curr(para_t);


Estimator::Estimator(): f_manager{Rs}
{
    ROS_INFO("init begins");

    kdtreeCornerLast.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    kdtreeSurfLast.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());
    cornerPointsSharp.reset(new pcl::PointCloud<PointType>());
    cornerPointsLessSharp.reset(new pcl::PointCloud<PointType>());
    surfPointsFlat.reset(new pcl::PointCloud<PointType>());
    surfPointsLessFlat.reset(new pcl::PointCloud<PointType>());
    laserCloudCornerLast.reset(new pcl::PointCloud<PointType>());
    laserCloudSurfLast.reset(new pcl::PointCloud<PointType>());
    laserCloudFullRes.reset(new pcl::PointCloud<PointType>());

    //for rtabmap
    for(int i=0; i<WINDOW_SIZE+1; ++i)
    {
        pre_integrations[i] = nullptr;
        pre_integrations_wheel[i] = nullptr;
    }
    tmp_pre_integration = nullptr;
    tmp_wheel_pre_integration = nullptr;
    last_marginalization_info = nullptr;

    initThreadFlag = false;
    clearState();
}

Estimator::~Estimator()
{
    if (MULTIPLE_THREAD)
    {
        processThread.join();
        cout << "join thread \n" << endl;
    }
}

void Estimator::clearState()
{
    mProcess.lock();
    while(!accBuf.empty())
        accBuf.pop();
    while(!gyrBuf.empty())
        gyrBuf.pop();
    while(!featureBuf.empty())
        featureBuf.pop();

    while(!wheelVelBuf.empty())
        wheelVelBuf.pop();
    while(!wheelGyrBuf.empty())
        wheelGyrBuf.pop();

    while(!angular_buf.empty())
        angular_buf.pop();

    while(!wheel_velocity_buf.empty())
        wheel_velocity_buf.pop();

    while(!yaw_sum_vec.empty())
        yaw_sum_vec.pop_front();

    while(!xy_sum_vec.empty())
        xy_sum_vec.pop_front();

    while(!LS_list.empty())
        LS_list.pop_front();

    prevTime = -1;
    prevTime_wheel = -1;
    prevTime_lidar = -1;
    curTime = 0;
    curTime_wheel = 0;
    curTime_lidar = -1;
    openExEstimation = 0;
    openPlaneEstimation = 0;
    openExWheelEstimation = 0;
    openIxEstimation = 0;
    initP = Eigen::Vector3d(0, 0, 0);
    initR = Eigen::Matrix3d::Identity();
    inputImageCnt = 0;
    initFirstPoseFlag = false;
    tmp_R.setZero();

    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        Vb[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();
        vel_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
        {
            delete pre_integrations[i];
        }
        pre_integrations[i] = nullptr;

        dt_buf_wheel[i].clear();
        linear_velocity_buf_wheel[i].clear();
        angular_velocity_buf_wheel[i].clear();
        if (pre_integrations_wheel[i] != nullptr)
        {
            delete pre_integrations_wheel[i];
        }
        pre_integrations_wheel[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }
    rwg = Matrix3d::Identity();

    tio = Vector3d::Zero();
    rio = Matrix3d::Identity();
    til = Vector3d::Zero();;
    ril = Matrix3d::Identity();
    sx = 0;
    sy = 0;
    sw = 0;

    rpw = Matrix3d::Identity();
    zpw = 0;

    first_imu = false,
    first_wheel = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();
    all_image_frame_init.clear();
    all_image_frame_init_now.clear();
    yaw_sum_vec.clear();
    xy_sum_vec.clear();
    static_flag = true;
    static_init_flag = false;


    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (tmp_wheel_pre_integration != nullptr)
        delete tmp_wheel_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    tmp_wheel_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();
    imu_data_list.clear();

    f_manager.clearState();

    failure_occur = 0;

    mProcess.unlock();
}

void Estimator::setParameter()
{
    mProcess.lock();
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
        cout << " exitrinsic cam " << i << endl  << ric[i] << endl << tic[i].transpose() << endl;
    }
    tio = TIO;
    tio_0 = TIO;
    rio = RIO;
    rio_0 = RIO;
    cout << " exitrinsic wheel " << endl  << rio << endl << tio.transpose() << endl;
    sx = SX;
    sy = SY;
    sw = SW;
    cout << " initrinsic wheel " << endl  << sx << " " << sy << " " << sw << endl;
    ril = RIL;
    til = TIL;
    cout << " exitrinsic lidar " << endl  << ril << endl << til.transpose() << endl;

    f_manager.setRic(ric);
    ProjectionTwoFrameOneCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionTwoFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionOneFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    td = TD;
    td_wheel = TD_WHEEL;
    g = G;
    cout << "set g " << g.transpose() << endl;
    featureTracker.readIntrinsicParameter(CAM_NAMES);

    std::cout << "MULTIPLE_THREAD is " << MULTIPLE_THREAD << '\n';
    if (MULTIPLE_THREAD && !initThreadFlag)
    {
        initThreadFlag = true;
        processThread = std::thread(&Estimator::processMeasurements, this);
    }
    mProcess.unlock();
}

void Estimator::changeSensorType(int use_imu, int use_stereo)
{
    bool restart = false;
    mProcess.lock();
    if(!use_imu && !use_stereo)
        printf("at least use two sensors! \n");
    else
    {
        if(USE_IMU != use_imu)
        {
            USE_IMU = use_imu;
            if(USE_IMU)
            {
                // reuse imu; restart system
                restart = true;
            }
            else
            {
                if (last_marginalization_info != nullptr)
                    delete last_marginalization_info;

                tmp_pre_integration = nullptr;
                last_marginalization_info = nullptr;
                last_marginalization_parameter_blocks.clear();
            }
        }
        
        STEREO = use_stereo;
        printf("use imu %d use stereo %d\n", USE_IMU, STEREO);
    }
    mProcess.unlock();
    if(restart)
    {
        clearState();
        setParameter();
    }
}

void Estimator::inputImage(double t, const cv::Mat &_img, const cv::Mat &_img1)
{
    inputImageCnt++;
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    TicToc featureTrackerTime;

    if(_img1.empty())
        featureFrame = featureTracker.trackImage(t, _img);
    else
        featureFrame = featureTracker.trackImage(t, _img, _img1);
    //printf("featureTracker time: %f\n", featureTrackerTime.toc());

    if (SHOW_TRACK)
    {
        cv::Mat imgTrack = featureTracker.getTrackImage();
        pubTrackImage(imgTrack, t);
    }
    
    if(MULTIPLE_THREAD)
    {
#if 1
        if(inputImageCnt % 2 == 0) // 30 Hz to 15 Hz
        {
            mBuf.lock();
            featureBuf.push(make_pair(t, featureFrame));
            mBuf.unlock();
        }
#else
        mBuf.lock();
        featureBuf.push(make_pair(t, featureFrame));
        mBuf.unlock();
#endif
    }
    else
    {
        mBuf.lock();
        featureBuf.push(make_pair(t, featureFrame));
        mBuf.unlock();
        TicToc processTime;
        processMeasurements();
        printf("process time: %f\n", processTime.toc());
    }
    
}
void Estimator::inputFeature(double t, const vector<cv::Point2f>& _features0, const vector<cv::Point2f>& _features1)
{
    inputImageCnt++;
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    TicToc featureTrackerTime;

    if(_features1.empty())
        featureFrame = featureTracker.trackFeature(t, _features0);
    else
        featureFrame = featureTracker.trackFeature(t, _features0, _features1);
    //printf("featureTracker time: %f\n", featureTrackerTime.toc());

    if (SHOW_TRACK)
    {
        cv::Mat imgTrack = featureTracker.getTrackImage();
        pubTrackImage(imgTrack, t);
    }

    if(MULTIPLE_THREAD)
    {
        if(inputImageCnt % 2 == 0)
        {
            mBuf.lock();
            featureBuf.push(make_pair(t, featureFrame));
            mBuf.unlock();
        }
    }
    else
    {
        mBuf.lock();
        featureBuf.push(make_pair(t, featureFrame));
        mBuf.unlock();
        TicToc processTime;
        processMeasurements();

//        if(solver_flag == SolverFlag::NON_LINEAR){
//            std::ofstream ofs(PROCESS_TIME_PATH, ios::app);
//            if (!ofs.is_open()) {
//                ROS_WARN("cannot open %s", PROCESS_TIME_PATH.c_str());
//            }
//            ofs << processTime.toc()<<std::endl;
//        }

        printf("process time: %f\n", processTime.toc());
    }

}
void Estimator::inputPointCloud(double t, const pcl::PointCloud<pcl::PointXYZ> &laserCloudIn)
{
    FeatureExtractor featureExtractor;
    featureExtractor.processPointCLoud(laserCloudIn);
    if(MULTIPLE_THREAD)
    {
        mBuf.lock();
        cornerSharpBuf.push(make_pair(t,featureExtractor.cornerPointsSharp));
        cornerLessSharpBuf.push(make_pair(t,featureExtractor.cornerPointsLessSharp));
        surfFlatBuf.push(make_pair(t,featureExtractor.surfPointsFlat));
        surfLessFlatBuf.push(make_pair(t,featureExtractor.surfPointsLessFlat));
        fullPointsBuf.push(make_pair(t,featureExtractor.laserCloud));
        mBuf.unlock();
    }
    else
    {
        mBuf.lock();
        cornerSharpBuf.push(make_pair(t,featureExtractor.cornerPointsSharp));
        cornerLessSharpBuf.push(make_pair(t,featureExtractor.cornerPointsLessSharp));
        surfFlatBuf.push(make_pair(t,featureExtractor.surfPointsFlat));
        surfLessFlatBuf.push(make_pair(t,featureExtractor.surfPointsLessFlat));
        fullPointsBuf.push(make_pair(t,featureExtractor.laserCloud));
        mBuf.unlock();
        TicToc processTime;
        processMeasurements();
        printf("process time: %f\n", processTime.toc());
    }

}
void Estimator::inputGroundtruth(double t, Eigen::Matrix<double,7,1>& _data)
{

        mGTBuf.lock();
        groundtruthBuf.push(make_pair(t, _data));
        mGTBuf.unlock();

}
void Estimator::inputIMU(double t, const Vector3d &linearAcceleration, const Vector3d &angularVelocity)
{
    mBuf.lock();
    accBuf.push(make_pair(t, linearAcceleration));
    gyrBuf.push(make_pair(t, angularVelocity));
    //printf("input imu with time %f \n", t);
    mBuf.unlock();

    if (solver_flag == NON_LINEAR)
    {
        mPropagate.lock();
        fastPredictIMU(t, linearAcceleration, angularVelocity);
        pubLatestOdometry(latest_P, latest_Q, latest_V, t);
        mPropagate.unlock();
    }
}
void Estimator::inputWheel(double t, const Vector3d &linearVelocity, const Vector3d &angularVelocity)
{
    mWheelBuf.lock();
    wheelVelBuf.push(make_pair(t, linearVelocity));
    wheelGyrBuf.push(make_pair(t, angularVelocity));
    //printf("input imu with time %f \n", t);
    mWheelBuf.unlock();

    if (solver_flag == NON_LINEAR)
    {
        mWheelPropagate.lock();
        fastPredictWheel(t, linearVelocity, angularVelocity);
        pubWheelLatestOdometry(latest_P_wheel, latest_Q_wheel, latest_V_wheel, t);
        Eigen::Quaterniond q;
        Eigen::Vector3d p;
        Eigen::Vector3d v;
        fastPredictPureWheel(t, linearVelocity, angularVelocity, p, q, v);
        pubPureWheelLatestOdometry(p, q, v, t);
        mWheelPropagate.unlock();
    }
}
void Estimator::inputFeature(double t, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &featureFrame)
{
    mBuf.lock();
    featureBuf.push(make_pair(t, featureFrame));
    mBuf.unlock();

    if(!MULTIPLE_THREAD)
        processMeasurements();
}


bool Estimator::getIMUInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &accVector,
                                vector<pair<double, Eigen::Vector3d>> &gyrVector)
{
    if(accBuf.empty())
    {
        printf("not receive imu\n");
        return false;
    }
    //printf("get imu from %f %f\n", t0, t1);
    //printf("imu front time %f   imu end time %f\n", accBuf.front().first, accBuf.back().first);
    if(t1 <= accBuf.back().first)
    {
        while (accBuf.front().first <= t0)
        {
            accBuf.pop();
            gyrBuf.pop();
        }
        while (accBuf.front().first < t1)
        {
            accVector.push_back(accBuf.front());
            accBuf.pop();
            gyrVector.push_back(gyrBuf.front());
            gyrBuf.pop();
        }
        accVector.push_back(accBuf.front());
        gyrVector.push_back(gyrBuf.front());
    }
    else
    {
        printf("wait for imu\n");
        return false;
    }
    return true;
}

bool Estimator::getWheelInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &velVector,
                               vector<pair<double, Eigen::Vector3d>> &gyrVector)
{
    if(wheelVelBuf.empty())
    {
        printf("not receive wheel\n");
        return false;
    }
    //printf("get imu from %f %f\n", t0, t1);
    //printf("imu front time %f   imu end time %f\n", wheelVelBuf.front().first, wheelVelBuf.back().first);
    if(t1 <= wheelVelBuf.back().first)
    {
        while (wheelVelBuf.front().first <= t0)
        {
            wheelVelBuf.pop();
            wheelGyrBuf.pop();
        }
        while (wheelVelBuf.front().first < t1)
        {
            velVector.push_back(wheelVelBuf.front());
            wheelVelBuf.pop();
            gyrVector.push_back(wheelGyrBuf.front());
            wheelGyrBuf.pop();
        }
        velVector.push_back(wheelVelBuf.front());
        gyrVector.push_back(wheelGyrBuf.front());
    }
    else
    {
        printf("wait for wheel\n");
        return false;
    }
//    ROS_INFO("velVector.size: %d", static_cast<int>(velVector.size()));
    return true;
}

bool Estimator::IMUAvailable(double t)
{
    if(!accBuf.empty() && t <= accBuf.back().first)
        return true;
    else
        return false;
}
bool Estimator::WheelAvailable(double t)
{
    if(!wheelVelBuf.empty() && t <= wheelVelBuf.back().first)
        return true;
    else
        return false;
}
bool Estimator::LidarAvailable(double t)
{
    if(!cornerSharpBuf.empty() && t >= cornerSharpBuf.front().first)
        return true;
    else
        return false;
}
void Estimator::processMeasurements()
{
    while (1)
    {
        //printf("process measurments\n");
        // 时间 特征点ID 图像id xyz_uv_vel
        pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > feature;
        vector<pair<double, Eigen::Vector3d>> accVector, gyrVector;
        vector<pair<double, Eigen::Vector3d>> velWheelVector, gyrWheelVector;
        bool isLidar = false;
        if(!featureBuf.empty())
        {
            feature = featureBuf.front();
            curTime = feature.first + td;
//            ROS_WARN_STREAM("curTime " << curTime);
//            if(!cornerLessSharpBuf.empty()){
//                double lidarTime = cornerLessSharpBuf.front().first;
//                ROS_WARN_STREAM("lidarTime " << lidarTime);
//                cornerLessSharpBuf.pop();
//            }
            curTime_wheel = curTime - td_wheel;
            while(1)
            {
                if ((!USE_IMU  || IMUAvailable(feature.first + td)))
                    break;
                else
                {
                    printf("wait for imu ... \n");
                    if (! MULTIPLE_THREAD)
                        return;
                    std::chrono::milliseconds dura(5);
                    std::this_thread::sleep_for(dura);
                }
            }
            while(1)
            {
                if ((!USE_WHEEL  || WheelAvailable(feature.first + td - td_wheel)))
                    break;
                else
                {
                    printf("wait for wheel ... \n");
                    if (! MULTIPLE_THREAD)
                        return;
                    std::chrono::milliseconds dura(5);
                    std::this_thread::sleep_for(dura);
                }
            }
            if (USE_LIDAR && LidarAvailable(feature.first)){
                isLidar = true;
                cornerPointsSharpFeature = cornerSharpBuf.front();
                cornerSharpBuf.pop();
                cornerPointsLessSharpFeature = cornerLessSharpBuf.front();
                cornerLessSharpBuf.pop();
                surfPointsFlatFeature = surfFlatBuf.front();
                surfFlatBuf.pop();
                surfPointsLessFlatFeature = surfLessFlatBuf.front();
                surfLessFlatBuf.pop();
                curTime_lidar = cornerPointsSharpFeature.first;
                *cornerPointsSharp = cornerPointsSharpFeature.second;
                *cornerPointsLessSharp = cornerPointsLessSharpFeature.second;
                *surfPointsFlat = surfPointsFlatFeature.second;
                *surfPointsLessFlat = surfPointsLessFlatFeature.second;
                //TODO add here
                int idx = 0;
                while(camera_in_lidarFrame_num > 1) {
//                    s_1 = (feature.first - cornerPointsSharpFeature.first) / SCAN_PERIOD;
                    double t0 = camera_in_lidarFrame[idx].first;
                    double t1 = camera_in_lidarFrame[idx + 1].first;
                    LaserOdometryProcess(t0, t1);
                    idx ++;
                    camera_in_lidarFrame_num--;
                }
                camera_in_lidarFrame.clear();
                camera_in_lidarFrame_num = 0;
            }
            camera_in_lidarFrame.emplace_back(curTime,camera_in_lidarFrame_num);
            camera_in_lidarFrame_num++;
            mBuf.lock();
            if(USE_IMU)
                //获取两图像帧时间戳之间的加速度和陀螺仪数据
                getIMUInterval(prevTime, curTime, accVector, gyrVector);// 对imu的时间进行判断，将队列里的imu数据放入到accVector和gyrVector中，完成之后返回true
            featureBuf.pop();//每次运行完之后都删除featureBuf中的元素，直到为空，已经把要删除的这个值给了feature
            mBuf.unlock();
            //TODO确认下是否统一成mBuf
            mWheelBuf.lock();
            if(USE_WHEEL)
                getWheelInterval(prevTime_wheel, curTime_wheel, velWheelVector, gyrWheelVector);
            mWheelBuf.unlock();
#if !WHEEL
            if(USE_IMU)
            {
                if(!initFirstPoseFlag)
                    initFirstIMUPose(accVector);
                for(size_t i = 0; i < accVector.size(); i++)
                {
                    //获取两帧IMU数据之间的dt
                    double dt;
                    if(i == 0)
                        dt = accVector[i].first - prevTime;
                    else if (i == accVector.size() - 1)
                        dt = curTime - accVector[i - 1].first;
                    else
                        dt = accVector[i].first - accVector[i - 1].first;
                    //预积分,并通过惯性解算得到Rs, Ps, Vs初值，用于后续processImage
                    processIMU(accVector[i].first, dt, accVector[i].second, gyrVector[i].second);
                }
            }
            if(USE_WHEEL)
            {
                for(size_t i = 0; i < velWheelVector.size(); i++)
                {
                    //获取两帧wheel数据之间的dt
                    double dt;
                    if(i == 0)
                        dt = velWheelVector[i].first - prevTime_wheel;
                    else if (i == velWheelVector.size() - 1)
                        dt = curTime_wheel - velWheelVector[i - 1].first;
                    else
                        dt = velWheelVector[i].first - velWheelVector[i - 1].first;
                    //预积分
                    processWheel(velWheelVector[i].first, dt, velWheelVector[i].second, gyrWheelVector[i].second);
                }
            }
#else
            if(USE_IMU && USE_WHEEL)
            {
                /* IMU */
                size_t k=0;
                /* IMU */

                for(size_t i = 0; i < velWheelVector.size(); i++)
                {
                    //获取两帧wheel数据之间的dt
                    double dt;
                    if(i == 0) {
                        dt = velWheelVector[i].first - prevTime_wheel;

                        /* IMU */
                        while(accVector[k].first < velWheelVector[i].first && k < accVector.size()){
//                            ROS_WARN_STREAM("k-a is " << k);
                            Vector3d tmpGyrVec = gyrVector[k].second;
                            //获取两帧IMU数据之间的dt
                            double dt_0;
                            if(k == 0)
                                dt_0 = accVector[k].first - prevTime;
                            else if (k == accVector.size() - 1)
                                dt_0 = curTime - accVector[k - 1].first;
                            else
                                dt_0 = accVector[k].first - accVector[k - 1].first;
                            //预积分,并通过惯性解算得到Rs, Ps, Vs初值，用于后续processImage
                            processIMU(accVector[k].first, dt_0, accVector[k].second, tmpGyrVec, velWheelVector[i].second);
                            k++;
                        }
                        /* IMU */
                    }
                    else if (i == velWheelVector.size() - 1) {
                        dt = curTime_wheel - velWheelVector[i - 1].first;

                        /* IMU */
                        Vector3d inter_vel = (velWheelVector[i-1].second + velWheelVector[i].second) / 2.0;
                        while(accVector[k].first > velWheelVector[i-1].first && accVector[k].first < velWheelVector[i].first && k < accVector.size()){
                            Vector3d tmpGyrVec = gyrVector[k].second;
                            //获取两帧IMU数据之间的dt
                            double dt_0;
                            if(k == 0)
                                dt_0 = accVector[k].first - prevTime;
                            else if (k == accVector.size() - 1)
                                dt_0 = curTime - accVector[k - 1].first;
                            else
                                dt_0 = accVector[k].first - accVector[k - 1].first;
                            //预积分,并通过惯性解算得到Rs, Ps, Vs初值，用于后续processImage
                            processIMU(accVector[k].first, dt_0, accVector[k].second, tmpGyrVec, inter_vel);
                            k++;
                        }
                        /* IMU */

                        /* IMU */
                        while(accVector[k].first > velWheelVector[i].first && k < accVector.size()){
                            Vector3d tmpGyrVec = gyrVector[k].second;
                            //获取两帧IMU数据之间的dt
                            double dt_0;
                            if(k == 0)
                                dt_0 = accVector[k].first - prevTime;
                            else if (k == accVector.size() - 1)
                                dt_0 = curTime - accVector[k - 1].first;
                            else
                                dt_0 = accVector[k].first - accVector[k - 1].first;
                            //预积分,并通过惯性解算得到Rs, Ps, Vs初值，用于后续processImage
                            processIMU(accVector[k].first, dt_0, accVector[k].second, tmpGyrVec, velWheelVector[i].second);
                            k++;
                        }
                        /* IMU */
                    }
                    else {
                        dt = velWheelVector[i].first - velWheelVector[i - 1].first;

                        /* IMU */
                        Vector3d inter_vel = (velWheelVector[i-1].second + velWheelVector[i].second) / 2.0;
                        while(accVector[k].first > velWheelVector[i-1].first && accVector[k].first < velWheelVector[i].first && k < accVector.size()){
                            Vector3d tmpGyrVec = gyrVector[k].second;
                            //获取两帧IMU数据之间的dt
                            double dt_0;
                            if(k == 0)
                                dt_0 = accVector[k].first - prevTime;
                            else if (k == accVector.size() - 1)
                                dt_0 = curTime - accVector[k - 1].first;
                            else
                                dt_0 = accVector[k].first - accVector[k - 1].first;
                            //预积分,并通过惯性解算得到Rs, Ps, Vs初值，用于后续processImage
                            processIMU(accVector[k].first, dt_0, accVector[k].second, tmpGyrVec, inter_vel);
                            k++;
                        }
                        /* IMU */
                    }
                    //预积分
                    processWheel(velWheelVector[i].first, dt, velWheelVector[i].second, gyrWheelVector[i].second);
                }
            }
#endif

//            static_init_flag = static_initialize();
            if (static_init_flag) {
                mProcess.lock();
//                static_initialize();
//                static_flag = false;
                processImage(feature.second, feature.first, isLidar);
                prevTime = curTime;
                prevTime_wheel = curTime_wheel;
                prevTime_lidar = curTime_lidar;

                printStatistics(*this, 0);

                std_msgs::Header header;
                header.frame_id = "world";
                header.stamp = ros::Time(feature.first);

                pubOdometry(*this, header);
                Eigen::Matrix<double, 7, 1> pose;
                pubGroundTruth(*this, header, pose, td);
                pubKeyPoses(*this, header);
                pubCameraPose(*this, header);
//                pubLidarPointCloud(*this,header);
                pubPointCloud(*this, header);
                pubKeyframe(*this);
                pubTF(*this, header);

                //可视化预积分积分的轨迹
//                if (USE_WHEEL && solver_flag == NON_LINEAR) {
//                    Eigen::Vector3d P;
//                    Eigen::Quaterniond Q;
//                    integrateWheelPreintegration(feature.first, P, Q, pose);
//                    pubWheelPreintegration(P, Q, header);
//                }

                mProcess.unlock();
            }else {
                static_init_flag = static_initialize();
            }
        }

        if (! MULTIPLE_THREAD)
            break;

        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }
}

//利用重力信息，初始化最开始状态中的旋转矩阵
void Estimator::initFirstIMUPose(vector<pair<double, Eigen::Vector3d>> &accVector)
{
    printf("init first imu pose\n");
    initFirstPoseFlag = true;
    //return;
    Eigen::Vector3d averAcc(0, 0, 0);
    int n = (int)accVector.size();
    for(size_t i = 0; i < accVector.size(); i++)
    {
        averAcc = averAcc + accVector[i].second;
    }
    averAcc = averAcc / n;
    printf("averge acc %f %f %f\n", averAcc.x(), averAcc.y(), averAcc.z());
    // 主要利用基于重力方向，得到的roll pitch信息，由于yaw信息通过重力方向，并不能恢复出来，因此减去yaw
    Matrix3d R0 = Utility::g2R(averAcc);
    double yaw = Utility::R2ypr(R0).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    Rs[0] = R0;
    cout << "init R0 " << endl << Rs[0] << endl;
    //Vs[0] = Vector3d(5, 0, 0);
}

void Estimator::initFirstPose(Eigen::Vector3d p, Eigen::Matrix3d r)
{
    Ps[0] = p;
    Rs[0] = r;
    initP = p;
    initR = r;
}

#if !WHEEL
void Estimator::processIMU(double t, double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)//第一帧IMU
    {
        first_imu = true;
        acc_0 = linear_acceleration;//始终为中值积分里的第一个数据
        gyr_0 = angular_velocity;
        oldest_time = t;
    }
    if (!static_init_flag && static_flag) {
        newest_time = t;
        imu_count++;
        ImuData imudata;
        imudata.timestamp = t;
        imudata.am = linear_acceleration;
        imudata.wm = angular_velocity;
        imu_data_list.push_back(imudata);
    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    if (frame_count != 0)//frame_count+1指代的是frame_cout帧的预积分， 可以从optimization中添加预积分约束中发现
    {
        //预积分
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
        //if(solver_flag != NON_LINEAR)
            tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);

        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;
        //这个是IMU惯性解算, Ps[j], Rs[j], Vs[j]的初值在processImage里通过上一帧的状态量进行初始化
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;

    // checkLine and checkZeroV
    if(angular_buf.size() < 5) {
        angular_buf.push(gyr_0);
        angular_v_sum += gyr_0.norm();
    }else{
        angular_v_sum -= angular_buf.front().norm();
        angular_buf.pop();
    }
}
#else
void Estimator::processIMU(double t, double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity, const Vector3d &wheel_velocity)
{
    if (!first_imu)//第一帧IMU
    {
        first_imu = true;
        acc_0 = linear_acceleration;//始终为中值积分里的第一个数据
        gyr_0 = angular_velocity;
        vel_0 = wheel_velocity;
        oldest_time = t;
    }
//    if (!static_init_flag && static_flag) {
        newest_time = t;
        imu_count++;
        ImuData imudata;
        imudata.timestamp = t;
        imudata.am = linear_acceleration;
        imudata.wm = angular_velocity;
        imu_data_list.push_back(imudata);
//    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, vel_0, Bas[frame_count], Bgs[frame_count]};//重建一个新的，为下一次processIMU做准备
    }
    if (frame_count != 0)//frame_count+1指代的是frame_cout帧的预积分， 可以从optimization中添加预积分约束中发现
    {
        //预积分
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity, wheel_velocity);
        //if(solver_flag != NON_LINEAR)
        tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity, wheel_velocity);

        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);
        vel_velocity_buf[frame_count].push_back(wheel_velocity);

        int j = frame_count;
        //这个是IMU惯性解算, Ps[j], Rs[j], Vs[j]的初值在processImage里通过上一帧的状态量进行初始化
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
    vel_0 = wheel_velocity;

    // checkLine and checkZeroV
    if(angular_buf.size() < 5) {
        angular_buf.push(gyr_0);
        angular_v_sum += gyr_0.norm();
    }else{
        angular_v_sum -= angular_buf.front().norm();
        angular_buf.pop();
    }
}
#endif

void Estimator::processWheel(double t, double dt, const Vector3d &linear_velocity, const Vector3d &angular_velocity)
{
    if (!first_wheel)//第一帧IMU
    {
        first_wheel = true;
        vel_0_wheel = linear_velocity;//始终为中值积分里的第一个数据
        gyr_0_wheel = angular_velocity;
    }

    if (!pre_integrations_wheel[frame_count])
    {
        pre_integrations_wheel[frame_count] = new WheelIntegrationBase{vel_0_wheel, gyr_0_wheel, sx, sy, sw, td_wheel};
    }
    if (frame_count != 0)
    {
        //预积分
        pre_integrations_wheel[frame_count]->push_back(dt, linear_velocity, angular_velocity);
        //if(solver_flag != NON_LINEAR)
        tmp_wheel_pre_integration->push_back(dt, linear_velocity, angular_velocity);

        dt_buf_wheel[frame_count].push_back(dt);
        linear_velocity_buf_wheel[frame_count].push_back(linear_velocity);
        angular_velocity_buf_wheel[frame_count].push_back(angular_velocity);

//        int j = frame_count;
//        //这个是IMU惯性解算, Ps[j], Rs[j], Vs[j]的初值在processImage里通过上一帧的状态量进行初始化
//        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
//        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
//        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
//        Vector3d un_acc_1 = Rs[j] * (linear_velocity - Bas[j]) - g;
//        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
//        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
//        Vs[j] += dt * un_acc;
    }
    vel_0_wheel = linear_velocity;
    gyr_0_wheel = angular_velocity;

    // checkLine and checkZeroV
    if(wheel_velocity_buf.size() < 5) {
        wheel_velocity_buf.push(vel_0_wheel);
        wheel_v_sum += vel_0_wheel.norm();
    }else{
        wheel_v_sum -= wheel_velocity_buf.front().norm();
        wheel_velocity_buf.pop();
    }
}
void Estimator::integrateWheelPreintegration( double t, Eigen::Vector3d& P, Eigen::Quaterniond& Q, const Eigen::Matrix<double, 7, 1>& pose){
    static bool firstPreint = false;
    static Eigen::Vector3d Pwo;
    static Eigen::Quaterniond Qwo;
    if(!firstPreint)
    {
        Pwo = Eigen::Vector3d(pose[4], pose[5], pose[6]);
        Qwo.w() = pose[0];
        Qwo.x() = pose[1];
        Qwo.y() = pose[2];
        Qwo.z() = pose[3];
        std::cout<<"integrateWheelPreintegration initial pose: \n"<<Pwo.transpose()<<std::endl<<Qwo.coeffs().transpose()<<std::endl;
        firstPreint = true;
    }else{
        Pwo = Qwo * all_image_frame[t].pre_integration_wheel->delta_p + Pwo.eval();
        Qwo = Qwo * all_image_frame[t].pre_integration_wheel->delta_q;
    }
    P = Pwo;
    Q = Qwo;
}
void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const double header, const bool isLidar)
{
    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());
    if (f_manager.addFeatureCheckParallax(frame_count, image, td))
    {
        marginalization_flag = MARGIN_OLD;
//        ROS_WARN_STREAM("KEY-FRAME");
        //printf("keyframe\n");
    }
    else
    {
        marginalization_flag = MARGIN_SECOND_NEW;
        //printf("non-keyframe\n");
    }

    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_INFO("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());
    Headers[frame_count] = header;

    ImageFrame imageframe(image, header);
    imageframe.pre_integration = tmp_pre_integration;//tmp_pre_intergration在前面processIMU计算得到
    imageframe.pre_integration_wheel = tmp_wheel_pre_integration;//tmp_wheel_pre_intergration在前面processWheel计算得到
    all_image_frame.insert(make_pair(header, imageframe));
#if WHEEL
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, vel_0, Bas[frame_count], Bgs[frame_count]};//重建一个新的，为下一次processIMU做准备
#else
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};//重建一个新的，为下一次processIMU做准备
#endif
    tmp_wheel_pre_integration = new WheelIntegrationBase{vel_0_wheel, gyr_0_wheel, sx, sy, sw, td_wheel};//重建一个新的，为下一次processWheel做准备



    if(ESTIMATE_EXTRINSIC == 2)
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");
        if (frame_count != 0)
        {
            //获取frame_count - 1 与 frame_count两个图像帧匹配的特征点，特征点在归一化相机坐标系
            vector<pair<Vector3d, Vector3d>> corres = f_manager.getCorresponding(frame_count - 1, frame_count);
            Matrix3d calib_ric;
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric))
            {
                ROS_WARN("initial extrinsic rotation calib success");
                ROS_WARN_STREAM("initial cam extrinsic rotation: " << endl << calib_ric);
                ric[0] = calib_ric;
                RIC[0] = calib_ric;
                ESTIMATE_EXTRINSIC = 1;
            }
        }
    }

    if (solver_flag == INITIAL)
    {
        // monocular + IMU initilization
        if (!STEREO && USE_IMU)
        {
//            ROS_WARN_STREAM("frame count: " << frame_count);
            if (frame_count == WINDOW_SIZE)
            {
                bool result = false;
                if(ESTIMATE_EXTRINSIC != 2 && ESTIMATE_EXTRINSIC_WHEEL != 2 && (header - initial_timestamp) > 0.1)
                {
                    result = initialStructure();
                    initial_timestamp = header;
                }
                if(result)
                {
                    static_flag = false;
                    optimization();
                    updateLatestStates();
                    solver_flag = NON_LINEAR;
                    slideWindow();
                    ROS_INFO("Initialization finish!");
                }
                else
                    slideWindow();
            }
        }

        // stereo + IMU initilization
        if(STEREO && USE_IMU)
        {
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
            f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
            if (frame_count == WINDOW_SIZE)
            {
                map<double, ImageFrame>::iterator frame_it;
                int i = 0;
                for (frame_it = all_image_frame.begin(); frame_it != all_image_frame.end(); frame_it++)
                {
                    frame_it->second.R = Rs[i];
                    frame_it->second.T = Ps[i];
                    i++;
                }
                if (!static_flag){
                    solveGyroscopeBias(all_image_frame, Bgs);
                    for (int i = 0; i <= WINDOW_SIZE; i++)
                    {
                        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
                    }
                }else{
                    for (int i = 0; i <= WINDOW_SIZE; i++)
                    {
                        pre_integrations[i]->repropagate(Bas[0], Bgs[0]); // use the ba and bg from static initialization to propagete the preintegration
                    }
                }
                optimization();
                updateLatestStates();
                solver_flag = NON_LINEAR;
                slideWindow();
                ROS_INFO("Initialization finish!");
            }
        }

        // stereo only initilization
        if(STEREO && !USE_IMU)
        {
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
            f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
            optimization();

            if(frame_count == WINDOW_SIZE)
            {
                optimization();
                updateLatestStates();
                solver_flag = NON_LINEAR;
                slideWindow();
                ROS_INFO("Initialization finish!");
            }
        }

        if(frame_count < WINDOW_SIZE)//滑窗未满，下一图像帧时刻的状态量用上一图像帧时刻的状态量进行初始化
        {                               //（PS：processIMU函数中，会对Ps、Vs、Rs用中值积分进行更新）
            frame_count++;
            int prev_frame = frame_count - 1;
            Ps[frame_count] = Ps[prev_frame];
            Vs[frame_count] = Vs[prev_frame];
            Rs[frame_count] = Rs[prev_frame];
            Bas[frame_count] = Bas[prev_frame];
            Bgs[frame_count] = Bgs[prev_frame];
            if (isLidar) {
                EdgePointsFrameVec[frame_count] = EdgePointsSum;
                PlanePointsFrameVec[frame_count] = PlanePointsSum;
            }
        }

        //设定平面参数的初始值
        if(solver_flag == NON_LINEAR && USE_PLANE)
            initPlane();


    }
    else//初始化成功，优化环节
    {
        TicToc t_solve;
        if(!USE_IMU) //当不存在imu时，使用pnp方法进行位姿预测;存在imu时使用imu积分进行预测
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
        f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
        optimization();
        set<int> removeIndex;
        outliersRejection(removeIndex);//基于重投影误差，检测外点
        if(FEATURE0_TOPIC.empty())
            f_manager.removeOutlier(removeIndex);
        if (! MULTIPLE_THREAD)
        {
//            ROS_INFO("removeIndex size: %d", removeIndex.size());
            if(FEATURE0_TOPIC.empty())
                featureTracker.removeOutliers(removeIndex); //若路标点为外点，则对前端图像跟踪部分的信息进行剔除更新;主要包括prev_pts, ids， track_cnt
            predictPtsInNextFrame(); //预测路标点在下一时刻左图中的坐标，基于恒速模型
        }
            
        ROS_DEBUG("solver costs: %fms", t_solve.toc());

        if (failureDetection())//默认直接返回false，具体判定失败的条件可根据具体场景修正
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            clearState(); //清除状态、重新设置参数，相当于重新开启vio
            setParameter();
            ROS_WARN("system reboot!");
            return;
        }

        slideWindow();
        f_manager.removeFailures();
        // prepare output of VINS
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
        updateLatestStates();



//        Eigen::Matrix3d rpw_temp = (Rs[WINDOW_SIZE] * rio).transpose();
//        Eigen::AngleAxisd qpw_temp(rpw_temp);
//        double zpw_temp = -(rpw_temp * (Ps[WINDOW_SIZE] + Rs[WINDOW_SIZE] * tio))[2];
//        std::cout<<"qpw: "<<qpw_temp.axis().transpose()<<" zpw: "<<zpw_temp<<std::endl;
    }
}
void Estimator::initPlane(){
    double zpws = 0.0;
    std::vector<Eigen::Quaterniond> qpws;
    for (int i = 0; i < frame_count; ++i) {
        Eigen::Matrix3d rpw_temp = (Rs[i] * rio).transpose();
        qpws.emplace_back(Eigen::Quaterniond(rpw_temp));
        zpws += -(rpw_temp * (Ps[i] + Rs[i] * tio))[2];
    }
    Eigen::Quaterniond qpw = Utility::quaternionAverage(qpws);
    rpw = qpw.toRotationMatrix();
    //对于平面，yaw是多余的，设为0
    double yaw = Utility::R2ypr(rpw).x();
    rpw = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * rpw;
    zpw = zpws / frame_count;
    std::cout<<"Init plane:  rpw: "<<Eigen::AngleAxisd(rpw).axis().transpose()<<" zpw: "<<zpw<<std::endl;
}

bool Estimator::checkLine()
{
    double aver_w = angular_v_sum / angular_buf.size();
    ROS_WARN_STREAM("average angular velocity  %f!" << aver_w << " buf size" << angular_buf.size());
    if (aver_w > 0.2) {  // 0.05 kaist  0.2 realsensed455
        ROS_INFO("Corner checked!");
        return false;
    }
    else {
        ROS_INFO("Line checked!");
        return true;
    }
}

bool Estimator::checkZeroV()
{
    double aver_w = wheel_v_sum / wheel_velocity_buf.size();
    ROS_WARN_STREAM("average wheel velocity  %f!" << aver_w << " buf size" << wheel_velocity_buf.size());
    if (aver_w > 0.01) {
        return false;
    }
    else {
        ROS_INFO("Zero Velocity checked!");
        return true;
    }
}

bool Estimator::static_initialize()
{

    // First lets collect a window of IMU readings from the newest measurement to the oldest
    vector<ImuData> window_1to0, window_2to1;
    for (const ImuData &data : imu_data_list) {
        if (data.timestamp > newest_time - 0.5 * init_window_time && data.timestamp <= newest_time - 0.0 * init_window_time) {
            window_1to0.push_back(data);
        }
        if (data.timestamp > newest_time - 1.0 * init_window_time && data.timestamp <= newest_time - 0.5 * init_window_time) {
            window_2to1.push_back(data);
        }
    }

    // Return if both of these failed
    if (window_1to0.size() < 2 || window_2to1.size() < 2) {
        ROS_WARN("[init-s]: unable to select window of IMU readings, not enough readings\n");
        return false;
    }

    // Calculate the sample variance for the newest window from 1 to 0
    Eigen::Vector3d a_avg_1to0 = Eigen::Vector3d::Zero();
    for (const ImuData &data : window_1to0) {
        a_avg_1to0 += data.am;
    }
    a_avg_1to0 /= (int)window_1to0.size();
    double a_var_1to0 = 0;
    for (const ImuData &data : window_1to0) {
        a_var_1to0 += (data.am - a_avg_1to0).dot(data.am - a_avg_1to0);
    }
    a_var_1to0 = std::sqrt(a_var_1to0 / ((int)window_1to0.size() - 1));

    // Calculate the sample variance for the second newest window from 2 to 1
    Eigen::Vector3d a_avg_2to1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d w_avg_2to1 = Eigen::Vector3d::Zero();
    for (const ImuData &data : window_2to1) {
        a_avg_2to1 += data.am;
        w_avg_2to1 += data.wm;
    }
    a_avg_2to1 = a_avg_2to1 / window_2to1.size();
    w_avg_2to1 = w_avg_2to1 / window_2to1.size();
    double a_var_2to1 = 0;
    for (const ImuData &data : window_2to1) {
        a_var_2to1 += (data.am - a_avg_2to1).dot(data.am - a_avg_2to1);
    }
    a_var_2to1 = std::sqrt(a_var_2to1 / ((int)window_2to1.size() - 1));
    // ROS_WARN(BOLDGREEN "[init-s]: IMU excitation, %.4f,%.4f\n" RESET, a_var_1to0, a_var_2to1);

    // If it is below the threshold and we want to wait till we detect a jerk
    if (a_var_1to0 < init_imu_thresh && a_var_2to1 < init_imu_thresh) {
        ROS_WARN("[init-s]: no IMU excitation, below threshold %.4f < %.4f\n", a_var_1to0, init_imu_thresh);
        return false;
    }

    // We should also check that the old state was below the threshold!
    // This is the case when we have started up moving, and thus we need to wait for a period of stationary motion
    if (a_var_2to1 > init_imu_thresh && a_var_1to0 > init_imu_thresh) {
        ROS_WARN("[init-s]: to much IMU excitation, above threshold %.4f > %.4f\n", a_var_2to1, init_imu_thresh);
        static_flag = false;
        return true;
    }

    if (a_var_2to1 > zupt_thresh && a_var_1to0 < zupt_thresh) {
        ROS_WARN("[init-s]: ZUPT Checked!");
        static_flag = true;
        return true;
    }

    if (!static_init_flag) {
        gb = a_avg_2to1;
        G.z() = gb.norm();
        Vector3d estimate_g = gb / gb.norm();
        ROS_WARN_STREAM("estimate g " << gb.transpose());
        ROS_WARN_STREAM("estimate g mag " << gb.norm());
        R0 = Utility::g2R(estimate_g);
        ROS_WARN_STREAM("estimate R0 ypr " << Utility::R2ypr(R0).transpose());
        ROS_WARN_STREAM("R0 * RIC     " << Utility::R2ypr(R0 * RIC[0]).transpose());

        Bgs[0] = w_avg_2to1;
        Bas[0] = a_avg_2to1 - R0.transpose() * G;
        ROS_WARN_STREAM("estimate ba " << Bas[0].transpose());
        ROS_WARN_STREAM("estimate bg " << Bgs[0].transpose());
        static_flag = true;
    }

    if (solver_flag == NON_LINEAR)
        static_flag = false;

    return true;

}



bool Estimator::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility
//    if (checkZeroV())
//        return false;
    // global sfm
    Quaterniond Q[frame_count + 1]; //R_w_c  from camera frame to world frame. tzhang
    Vector3d T[frame_count + 1]; // t_w_c
    map<int, Vector3d> sfm_tracked_points; //观测到的路标点的在世界坐标系的位置，索引为路标点的编号
    vector<SFMFeature> sfm_f;
    for (auto &it_per_id : f_manager.feature) //对所有路标点进行遍历
    {
        int imu_j = it_per_id.start_frame - 1;
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id; //路标点的编号设置tmp_feature的ID
        for (auto &it_per_frame : it_per_id.feature_per_frame) //对观测到路标点j的所有图像帧进行遍历
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.point;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()})); //构建路标点在左相机图像坐标系下的观测（去畸变）
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    if (!relativePose(relative_R, relative_T, l)) //通过本质矩阵求取滑窗最后一帧（WINDOW_SIZE）到图像帧l的旋转和平移变换
    {//共视点大于20个、视差足够大，才进行求取
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    GlobalSFM sfm; // 通过global sfm求取滑窗中的图像帧位姿，以及观测到的路标点的位置
    if(!sfm.construct(frame_count + 1, Q, T, l, //只有frame_count == WINDOW_SIZE才会调用initialStructure，此时frame_count即为WINDOW_SIZE
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
        ROS_INFO("global SFM failed!");
        marginalization_flag = MARGIN_OLD; //global sfm求解失败，对老的图像帧进行边缘化
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        if((frame_it->first) == Headers[i])
        {
            frame_it->second.is_key_frame = true;
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose(); //得到R_w_i 也即R_c0_i
            frame_it->second.T = T[i]; //TODO(tzhang):没有对imu与cam平移进行处理，是因为两者之间的平移量较小，影响不大？ 可优化
            i++;
            continue; //若图像帧在滑窗中，直接利用上述global sfm的结果乘以R_c_i，得到imu到world的变换矩阵即R_w_i
        }
        if((frame_it->first) > Headers[i]) // 时间戳比较，实现逐帧判断是不是滑窗内的帧
        {
            i++;
        }
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if(it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_INFO("Not enough points for solve pnp !");
            return false;
        }
        /**
         *bool cv::solvePnP(    求解pnp问题
         *   InputArray  objectPoints,   特征点的3D坐标数组
         *   InputArray  imagePoints,    特征点对应的图像坐标
         *   InputArray  cameraMatrix,   相机内参矩阵
         *   InputArray  distCoeffs,     畸变系数的输入向量
         *   OutputArray     rvec,       旋转向量, r_c_w
         *   OutputArray     tvec,       平移向量, t_c_w
         *   bool    useExtrinsicGuess = false, 为真则使用提供的初始估计值
         *   int     flags = SOLVEPNP_ITERATIVE 采用LM优化
         *)
         */
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_INFO("solve pnp fail!");
            return false;
        }
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp); // 通过PnP求解得到t_c0_c
        T_pnp = R_pnp * (-T_pnp);
        frame_it->second.R = R_pnp * RIC[0].transpose();//得到R_c0_i              R_c0_c to R_c0_i
        frame_it->second.T = T_pnp; //t_c0_c （PS：未考虑camera与imu之间的平移,由于尺度因子未知，此时不用考虑；在求解出尺度因子后，会考虑）
    }
    if (visualInitialAlign())
        return true;
    else
    {
        ROS_INFO("misalign visual structure with IMU");
        return false;
    }

}

bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
//    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);  //IMU与camera对准，计算更新了陀螺仪偏置bgs、重力向量gc0、尺度因子s、速度vk
    Vector3d wheel_s(sx,sy,sw);
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x, wheel_s, Bas);  //IMU与camera对准，计算更新了陀螺仪偏置bgs、重力向量gc0、尺度因子s、速度vk
//    sx = (x.tail<4>())(0);
//    sx=wheel_s(0);
//    sy=wheel_s(1);
//    sw=wheel_s(2);
//    ROS_WARN_STREAM("initial wheel_s: " << wheel_s.transpose());
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // save the initialization data for the following bas refinement
//    if (ONLY_INITIAL_WITH_WHEEL) {
//        for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != all_image_frame.end(); it++) {
//            ImageFrame imageFrame;
//            imageFrame = it->second;
//            all_image_frame_init.insert(make_pair(it->first, imageFrame));
//        }
//        g_init = g;
//    }

    // change state
    //Headers[i]存储图像帧时间戳，初始化过程中，仅边缘化老的图像帧，因此留在滑窗中的均为关键帧
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i]].R;
        Vector3d Pi = all_image_frame[Headers[i]].T;
        Ps[i] = Pi;//t_c0_c
        Rs[i] = Ri;//R_c0_b
        all_image_frame[Headers[i]].is_key_frame = true;
    }

//    double s = (x.tail<1>())(0);
//    double s = (x.tail<3>())(0);
//    double s=0.65;
    double s = (x.tail<4>())(0);

//    TIO = x.tail<3>();
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        //得到陀螺仪bias新值，预积分重新传播；此处的预积分为Estimator类中的预积分；
        //图像帧中预积分的重新传播，在计算bg后已经完成
        //TODO(tzhang)：Estimator与图像帧中均维护预积分，重复，是否有必要？待优化
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    //利用尺度因子，对t_c0_b的更新；注意此时用到了camera与imu之间的平移量
    for (int i = frame_count; i >= 0; i--)//根据该式可以得出Ps[0] = Pc0b0 = 0;
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);//sPc0ck - Rc0bk*tbc  sP_c0_bk
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            // update initial Vs
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3); // V_c0_bk
        }
    }
    ROS_WARN_STREAM("g0     " << g.transpose());  // g_c0
    if (static_flag){
        g = R0 * gb;  // g_w
        R0 = R0 * RIC[0];
        ROS_WARN_STREAM("g = R0 * gb " << g.transpose());
    }else {
        R0 = Utility::g2R(g);//根据基于当前世界坐标系计算得到的重力方向与实际重力方向差异，计算当前世界坐标系的修正量；
//    R0 = RIC[0];
        //注意：由于yaw不可观，修正量中剔除了yaw影响，也即仅将世界坐标系的z向与重力方向对齐
        yaw = Utility::R2ypr(R0 * Rs[0]).x();
        R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
        ROS_WARN_STREAM("RIC * R0.transpose()     " << Utility::R2ypr(RIC[0] * R0.transpose()).transpose());
        g = R0 * g;  // g_w
        ROS_WARN_STREAM("g = R0 * g " << g.transpose());
        //Matrix3d rot_diff = R0 * Rs[0].transpose();
    }
    Matrix3d rot_diff = R0;//将世界坐标系与重力方向对齐，之前的世界坐标系Rs[0]根据图像帧c0定义得到，并未对准到重力方向  R_w_c0
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];//t_w_bi
        Rs[i] = rot_diff * Rs[i];//R_w_bi
        Vs[i] = rot_diff * Vs[i];//V_w_bi
//        Vb[i] = R0 * RIC[0].transpose() * Rs[i].transpose() * Vs[i];
    }
    ROS_WARN_STREAM("g0     " << g.transpose());
    ROS_WARN_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());
    ROS_WARN_STREAM("R0  " << Utility::R2ypr(R0).transpose());
//    if(!ESTIMATE_RIO) {
//        dR = rio.transpose() * RIC[0] * R0.transpose();
        dR = rio_0.transpose();
//    }
//    ROS_WARN_STREAM("dR " << Utility::R2ypr(dR).transpose());
//    R0 = RIC[0];

    // end




    f_manager.clearDepth();//清除路标点状态，假定所有路标点逆深度均为估计；注意后续优化中路标点采用逆深度，而初始化过程中路标点采用三维坐标
    f_manager.triangulate(frame_count, Ps, Rs, tic, ric); //基于SVD的路标点三角化，双目情形：利用左右图像信息； 非双目情形：利用前后帧图像信息

    return true;
}

bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            double threshold = FEATURE0_TOPIC.empty()?30:8;
//            ROS_INFO("parallax threshold: %f",threshold);
            if(average_parallax * FOCAL_LENGTH > threshold && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * FOCAL_LENGTH, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        if(USE_IMU)
        {
            para_SpeedBias[i][0] = Vs[i].x();
            para_SpeedBias[i][1] = Vs[i].y();
            para_SpeedBias[i][2] = Vs[i].z();

            para_SpeedBias[i][3] = Bas[i].x();
            para_SpeedBias[i][4] = Bas[i].y();
            para_SpeedBias[i][5] = Bas[i].z();

            para_SpeedBias[i][6] = Bgs[i].x();
            para_SpeedBias[i][7] = Bgs[i].y();
            para_SpeedBias[i][8] = Bgs[i].z();
        }
    }

    para_Ex_NHC[0][0] = tio_0.x();
    para_Ex_NHC[0][1] = tio_0.y();
    para_Ex_NHC[0][2] = tio_0.z();
    Quaterniond q1{dR};
    para_Ex_NHC[0][3] = q1.x();
    para_Ex_NHC[0][4] = q1.y();
    para_Ex_NHC[0][5] = q1.z();
    para_Ex_NHC[0][6] = q1.w();



    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }

    para_Ex_Pose_wheel[0][0] = tio.x();
    para_Ex_Pose_wheel[0][1] = tio.y();
    para_Ex_Pose_wheel[0][2] = tio.z();
    Quaterniond q{rio};
    para_Ex_Pose_wheel[0][3] = q.x();
    para_Ex_Pose_wheel[0][4] = q.y();
    para_Ex_Pose_wheel[0][5] = q.z();
    para_Ex_Pose_wheel[0][6] = q.w();

    para_Ix_sx_wheel[0][0] = sx;
    para_Ix_sy_wheel[0][0] = sy;
    para_Ix_sw_wheel[0][0] = sw;

    Quaterniond q2{rpw};
    para_plane_R[0][0] = q.x();
    para_plane_R[0][1] = q.y();
    para_plane_R[0][2] = q.z();
    para_plane_R[0][3] = q.w();
    para_plane_Z[0][0] = zpw;

    para_Ex_Pose_lidar[0][0] = til.x();
    para_Ex_Pose_lidar[0][1] = til.y();
    para_Ex_Pose_lidar[0][2] = til.z();
    Quaterniond q3{ril};
    para_Ex_Pose_lidar[0][3] = q3.x();
    para_Ex_Pose_lidar[0][4] = q3.y();
    para_Ex_Pose_lidar[0][5] = q3.z();
    para_Ex_Pose_lidar[0][6] = q3.w();

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);

    para_Td[0][0] = td;
    para_Td_wheel[0][0] = td_wheel;
}

void Estimator::double2vector()
{
    Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Vector3d origin_P0 = Ps[0];

    if (failure_occur)
    {
        origin_R0 = Utility::R2ypr(last_R0);
        origin_P0 = last_P0;
        failure_occur = 0;
    }

    if(USE_IMU)
    {
        Vector3d origin_R00 = Utility::R2ypr(Quaterniond(para_Pose[0][6],
                                                          para_Pose[0][3],
                                                          para_Pose[0][4],
                                                          para_Pose[0][5]).toRotationMatrix());
        double y_diff = origin_R0.x() - origin_R00.x(); //regard c0 as the world coordinate so the following estimated position need to eliminate the effect of the yaw of Rs[0] (caculate the y_diff w.r.t Rs[0])   2022.03.28 by zhh
        // TODO
        Matrix3d rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));
        if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
        {
            ROS_DEBUG("euler singular point!");
            rot_diff = Rs[0] * Quaterniond(para_Pose[0][6],
                                           para_Pose[0][3],
                                           para_Pose[0][4],
                                           para_Pose[0][5]).toRotationMatrix().transpose();
        }

        tio_0 = Vector3d(para_Ex_NHC[0][0],
                         para_Ex_NHC[0][1],
                         para_Ex_NHC[0][2]);   // In fact here tio_0 is the toi according to the residual formula

        dR = Quaterniond(para_Ex_NHC[0][6],
                         para_Ex_NHC[0][3],
                         para_Ex_NHC[0][4],
                         para_Ex_NHC[0][5]).normalized().toRotationMatrix().transpose();
//        dR = Utility::ypr2R(Vector3d{Utility::R2ypr(dR).x(), 0 , 0});

        for (int i = 0; i <= WINDOW_SIZE; i++)
        {

            Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
            
            Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                    para_Pose[i][1] - para_Pose[0][1],
                                    para_Pose[i][2] - para_Pose[0][2]) + origin_P0;


            Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                        para_SpeedBias[i][1],
                                        para_SpeedBias[i][2]);



            Bas[i] = Vector3d(para_SpeedBias[i][3],
                              para_SpeedBias[i][4],
                              para_SpeedBias[i][5]);


            Bgs[i] = Vector3d(para_SpeedBias[i][6],
                              para_SpeedBias[i][7],
                              para_SpeedBias[i][8]);

        }

    }
    else
    {
        for (int i = 0; i <= WINDOW_SIZE; i++)
        {
            Rs[i] = Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
            
            Ps[i] = Vector3d(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);
        }
    }

    if(USE_IMU)
    {
        for (int i = 0; i < NUM_OF_CAM; i++)
        {
            tic[i] = Vector3d(para_Ex_Pose[i][0],
                              para_Ex_Pose[i][1],
                              para_Ex_Pose[i][2]);
            ric[i] = Quaterniond(para_Ex_Pose[i][6],
                                 para_Ex_Pose[i][3],
                                 para_Ex_Pose[i][4],
                                 para_Ex_Pose[i][5]).normalized().toRotationMatrix();
        }
    }

    if(USE_WHEEL)
    {
        tio = Vector3d(para_Ex_Pose_wheel[0][0],
                       para_Ex_Pose_wheel[0][1],
                       para_Ex_Pose_wheel[0][2]);
        rio = Quaterniond(para_Ex_Pose_wheel[0][6],
                          para_Ex_Pose_wheel[0][3],
                          para_Ex_Pose_wheel[0][4],
                          para_Ex_Pose_wheel[0][5]).normalized().toRotationMatrix();
#if !WHEEL
        ROS_WARN_STREAM("calib RIO " << endl << rio);
        ROS_WARN_STREAM("calib RIO ypr " << Utility::R2ypr_m(rio).transpose());  // atan
        ROS_WARN_STREAM("calib TIO     " << tio.transpose());
        ROS_WARN_STREAM("tio T    " << (-rio.transpose() * tio).transpose());
#endif
        sx = para_Ix_sx_wheel[0][0];
        sy = para_Ix_sy_wheel[0][0];
        sw = para_Ix_sw_wheel[0][0];

        td_wheel = para_Td_wheel[0][0];

    }

    if(USE_PLANE){
        zpw = para_plane_Z[0][0];
        rpw = Quaterniond(para_plane_R[0][3], para_plane_R[0][0], para_plane_R[0][1], para_plane_R[0][2]).normalized().toRotationMatrix();
    }

    if(USE_LIDAR){
        til = Vector3d(para_Ex_Pose_lidar[0][0],
                       para_Ex_Pose_lidar[0][1],
                       para_Ex_Pose_lidar[0][2]);
        ril = Quaterniond(para_Ex_Pose_lidar[0][6],
                          para_Ex_Pose_lidar[0][3],
                          para_Ex_Pose_lidar[0][4],
                          para_Ex_Pose_lidar[0][5]).normalized().toRotationMatrix();
        ROS_WARN_STREAM("calib RIL     " << endl << ril);  // atan
        ROS_WARN_STREAM("calib TIL     " << til.transpose());
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);

    if(USE_IMU)
        td = para_Td[0][0];

}

bool Estimator::failureDetection()
{
    return false;
    if (f_manager.last_track_num < 2)
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;
    }
    if (Bas[WINDOW_SIZE].norm() > 2.5)
    {
        ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    Vector3d tmp_P = Ps[WINDOW_SIZE];
    if ((tmp_P - last_P).norm() > 5)
    {
        //ROS_INFO(" big translation");
        //return true;
    }
    if (abs(tmp_P.z() - last_P.z()) > 1)
    {
        //ROS_INFO(" big z translation");
        //return true; 
    }
    Matrix3d tmp_R = Rs[WINDOW_SIZE];
    Matrix3d delta_R = tmp_R.transpose() * last_R;
    Quaterniond delta_Q(delta_R);
    double delta_angle;
    delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;
    if (delta_angle > 50)
    {
        ROS_INFO(" big delta_angle ");
        //return true;
    }
    return false;
}

void Estimator::optimization()
{
    TicToc t_whole, t_prepare;
    vector2double();

    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = NULL;
    loss_function = new ceres::HuberLoss(1.0);
    //loss_function = new ceres::CauchyLoss(1.0 / FOCAL_LENGTH);
    //ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
    for (int i = 0; i < frame_count + 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
//        if (solver_flag == INITIAL)
//            problem.SetParameterBlockConstant(para_Pose[i]);    // This is a stop button !!!
        if(USE_IMU) {
            problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);
//            problem.AddParameterBlock(para_SpeedBias[i], SIZE_BIAS);
        }
    }
//    ceres::LocalParameterization *local_parameterization = new RwgSO3LocalParameterization({2});
//    problem.AddParameterBlock(para_Gravity[0], SIZE_G, local_parameterization);
////    if (solver_flag == NON_LINEAR)
//        problem.SetParameterBlockConstant(para_Gravity[0]);    // This is a stop button !!!
    //TODO 为什么有IMU就不需要固定第一帧位姿
    if(!USE_IMU)
        problem.SetParameterBlockConstant(para_Pose[0]);

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
//        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        ceres::LocalParameterization *local_parameterization;
        if(ESTIMATE_EXTRINSIC) {
            // 固定某些量
            switch (CAM_EXT_ADJ_TYPE) {
                case CameraExtrinsicAdjustType::ADJUST_CAM_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 6});
                    break;
                case CameraExtrinsicAdjustType::ADJUST_CAM_ROTATION:
                    local_parameterization = new PoseSubsetParameterization({0, 1, 2, 6});
                    break;
                case CameraExtrinsicAdjustType::ADJUST_CAM_TRANSLATION:
                    local_parameterization = new PoseSubsetParameterization({3, 4, 5, 6});
                    break;
                case CameraExtrinsicAdjustType::ADJUST_CAM_NO_ROTATION_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 3, 4, 5, 6});
                    break;
                default:
                    local_parameterization = new PoseSubsetParameterization({});
            }
        } else
            local_parameterization = new PoseLocalParameterization();

        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);
        if ((ESTIMATE_EXTRINSIC && frame_count == WINDOW_SIZE && Vs[0].norm() > 0.2) || openExEstimation)
        {
            //ROS_INFO("estimate extinsic param");
            openExEstimation = 1;
        }
        else
        {
            //ROS_INFO("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);
        }
    }

    if(USE_WHEEL && !ONLY_INITIAL_WITH_WHEEL){
//        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();

        ceres::LocalParameterization *local_parameterization;
        if(ESTIMATE_EXTRINSIC_WHEEL) {
            // 固定某些量
            switch (WHEEL_EXT_ADJ_TYPE) {
                case WheelExtrinsicAdjustType::ADJUST_WHEEL_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 6});
                    break;
                case WheelExtrinsicAdjustType::ADJUST_WHEEL_ROTATION:
                    local_parameterization = new PoseSubsetParameterization({0, 1, 2, 6});
                    break;
                case WheelExtrinsicAdjustType::ADJUST_WHEEL_TRANSLATION:
                    local_parameterization = new PoseSubsetParameterization({3, 4, 5, 6});
                    break;
                case WheelExtrinsicAdjustType::ADJUST_WHEEL_NO_ROTATION_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 3, 4, 5, 6});
                    break;
                default:
                    local_parameterization = new PoseSubsetParameterization({});
            }
        } else
            local_parameterization = new PoseLocalParameterization();

        problem.AddParameterBlock(para_Ex_Pose_wheel[0], SIZE_POSE, local_parameterization);
        if ((ESTIMATE_EXTRINSIC_WHEEL && frame_count == WINDOW_SIZE && Vs[0].norm() > 0.2) || openExWheelEstimation)
        {
            //ROS_INFO("estimate extrinsic param");
            openExWheelEstimation = 1;
        }
        else
        {
            //ROS_INFO("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose_wheel[0]);
        }
        problem.AddParameterBlock(para_Ix_sx_wheel[0],1);
        problem.AddParameterBlock(para_Ix_sy_wheel[0],1);
        problem.AddParameterBlock(para_Ix_sw_wheel[0],1);
        if ((ESTIMATE_INTRINSIC_WHEEL && frame_count == WINDOW_SIZE && Vs[0].norm() > 0.2) || openIxEstimation)
        {
            //ROS_INFO("estimate intrinsic param");
            openIxEstimation = 1;
        }
        else
        {
            //ROS_INFO("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ix_sx_wheel[0]);
            problem.SetParameterBlockConstant(para_Ix_sy_wheel[0]);
            problem.SetParameterBlockConstant(para_Ix_sw_wheel[0]);
        }
    }

//    if(USE_PLANE){
//        ceres::LocalParameterization *local_parameterization = new OrientationSubsetParameterization(std::vector<int>{2});
//        problem.AddParameterBlock(para_plane_R[0], SIZE_ROTATION, local_parameterization);
//        problem.AddParameterBlock(para_plane_Z[0], 1);
//        if (frame_count == WINDOW_SIZE || openPlaneEstimation)
//        {
//            //ROS_INFO("estimate extrinsic param");
//            openPlaneEstimation = 1;
//        }
//        else
//        {
//            //ROS_INFO("fix extinsic param");
//            problem.SetParameterBlockConstant(para_plane_R[0]);
//            problem.SetParameterBlockConstant(para_plane_Z[0]);
//        }
//    }

    if(USE_LIDAR){
//        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();

        ceres::LocalParameterization *local_parameterization;
        if(ESTIMATE_EXTRINSIC_LIDAR) {
            // 固定某些量
            switch (LIDAR_EXT_ADJ_TYPE) {
                case LidarExtrinsicAdjustType::ADJUST_LIDAR_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 6});
                    break;
                case LidarExtrinsicAdjustType::ADJUST_LIDAR_ROTATION:
                    local_parameterization = new PoseSubsetParameterization({0, 1, 2, 6});
                    break;
                case LidarExtrinsicAdjustType::ADJUST_LIDAR_TRANSLATION:
                    local_parameterization = new PoseSubsetParameterization({3, 4, 5, 6});
                    break;
                case LidarExtrinsicAdjustType::ADJUST_LIDAR_NO_ROTATION_NO_Z:
                    local_parameterization = new PoseSubsetParameterization({2, 3, 4, 5, 6});
                    break;
                default:
                    local_parameterization = new PoseSubsetParameterization({});
            }
        } else
            local_parameterization = new PoseLocalParameterization();

        problem.AddParameterBlock(para_Ex_Pose_lidar[0], SIZE_POSE, local_parameterization);
        if ((ESTIMATE_EXTRINSIC_LIDAR && frame_count == WINDOW_SIZE && Vs[0].norm() > 0.2) || openExLidarEstimation)
        {
            //ROS_INFO("estimate extrinsic param");
            openExLidarEstimation = 1;
        }
        else
        {
            //ROS_INFO("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose_lidar[0]);
        }
    }


    problem.AddParameterBlock(para_Td[0], 1);
    problem.AddParameterBlock(para_Td_wheel[0], 1);
    //TODO 为什么判断 Vs[0].norm() < 0.2
    if (!ESTIMATE_TD || Vs[0].norm() < 0.2)
        problem.SetParameterBlockConstant(para_Td[0]);
    if (!ESTIMATE_TD_WHEEL || Vs[0].norm() < 0.2)
        problem.SetParameterBlockConstant(para_Td_wheel[0]);
    if (last_marginalization_info && last_marginalization_info->valid)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }
    if(USE_IMU)
    {
        for (int i = 0; i < frame_count; i++)
        {
            int j = i + 1;
            if (pre_integrations[j]->sum_dt > 10.0) //两图像帧之间时间过长，不使用中间的预积分
                continue;
#if !WHEEL
            IMUFactor* imu_factor = new IMUFactor(pre_integrations[j]);
            problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
#else
            IMUWheelLineFactor* imu_factor = new IMUWheelLineFactor(pre_integrations[j]);
            problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j], para_Ex_NHC[0]);
#endif
//            std::vector<const double *> parameters(4);
//            parameters[0] = para_Pose[i];
//            parameters[1] = para_SpeedBias[i];
//            parameters[2] = para_Pose[j];
//            parameters[3] = para_SpeedBias[j];
//            imu_factor->check(const_cast<double **>(parameters.data()));
        }
        ceres::LocalParameterization *local_parameterization;
        if (ESTIMATE_TIO)
        {
            if (ESTIMATE_RIO){
                local_parameterization = new PoseSubsetParameterization({6});
            }else{
                local_parameterization = new PoseSubsetParameterization({3,4,5,6});
            }
        }else{
            if (ESTIMATE_RIO){
                local_parameterization = new PoseSubsetParameterization({0,1,2,6});
            }else{
                local_parameterization = new PoseLocalParameterization();
            }
        }
        problem.AddParameterBlock(para_Ex_NHC[0], SIZE_NHC, local_parameterization);
        if(!ESTIMATE_TIO && !ESTIMATE_RIO || solver_flag == INITIAL){
            problem.SetParameterBlockConstant(para_Ex_NHC[0]);
        }
    }
    if(USE_WHEEL && !ONLY_INITIAL_WITH_WHEEL)
    {
        for (int i = 0; i < frame_count; i++)
        {
            int j = i + 1;
            if (pre_integrations_wheel[j]->sum_dt > 10.0) //两图像帧之间时间过长，不使用中间的预积分
                continue;
            // TODO 添加卡方检验
            WheelFactor* wheel_factor = new WheelFactor(pre_integrations_wheel[j]);
            problem.AddResidualBlock(wheel_factor, NULL, para_Pose[i], para_Pose[j], para_Ex_Pose_wheel[0], para_Ix_sx_wheel[0], para_Ix_sy_wheel[0], para_Ix_sw_wheel[0], para_Td_wheel[0]);

//            std::vector<const double *> parameters(7);
//            parameters[0] = para_Pose[i];
//            parameters[1] = para_Pose[j];
//            parameters[2] = para_Ex_Pose_wheel[0];
//            parameters[3] = para_Ix_sx_wheel[0];
//            parameters[4] = para_Ix_sy_wheel[0];
//            parameters[5] = para_Ix_sw_wheel[0];
//            parameters[6] = para_Td_wheel[0];
//            wheel_factor->check(const_cast<double **>(parameters.data()));
        }
    }

//    if(USE_WHEEL && ONLY_INITIAL_WITH_WHEEL)
//    {
//        for (int i = 0; i < frame_count; i++)
//        {
//            int j = i + 1;
//            if (pre_integrations[j]->sum_dt > 10.0) //两图像帧之间时间过长，不使用中间的预积分
//                continue;
//            // TODO 添加卡方检验
//            NonholonomicFactor* nonholonomic_factor = new NonholonomicFactor(pre_integrations[j]);
//            problem.AddResidualBlock(nonholonomic_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j], para_Ex_NHC[0]);
//
//        }
//    }

    if(USE_PLANE)
    {
        for (int i = 0; i < frame_count; i++)
        {
            PlaneFactor* plane_factor = new PlaneFactor();
            problem.AddResidualBlock(plane_factor, NULL, para_Pose[i], para_Ex_Pose_wheel[0], para_plane_R[0], para_plane_Z[0]);

//            std::vector<const double *> parameters(4);
//            parameters[0] = para_Pose[i];
//            parameters[1] = para_Ex_Pose_wheel[0];
//            parameters[2] = para_plane_R[0];
//            parameters[3] = para_plane_Z[0];
//            plane_factor->check(const_cast<double **>(parameters.data()));
        }
    }

    if(static_flag){
        for (int i = 0; i < frame_count; i++)
        {
            int j = i + 1;
            ZUPTFactor* zupt_factor = new ZUPTFactor();
            problem.AddResidualBlock(zupt_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
        }
    }



    if(USE_LIDAR)
    {
        for (int i = 0; i < frame_count; i++)
        {
            if (!EdgePointsFrameVec[i].empty()) {
                vector<EdgePoints> tmpEdgePoints = EdgePointsFrameVec[i];
//                ROS_WARN_STREAM("lidar i " << i);
                for (int j = 0; j < tmpEdgePoints.size(); ++j) {
                    ceres::CostFunction *lidar_edge_factor = LidarEdgeFactor::Create(tmpEdgePoints[i].curr_point,
                                                                                 tmpEdgePoints[i].last_point_a,
                                                                                 tmpEdgePoints[i].last_point_b,
                                                                                 tmpEdgePoints[i].s, s_1);
                    problem.AddResidualBlock(lidar_edge_factor, loss_function, para_Pose[i], para_Pose[i+1], para_Ex_Pose_lidar[0]);
                }
                vector<PlanePoints> tmpPlanePoints = PlanePointsFrameVec[i];
                for (int j = 0; j < tmpPlanePoints.size(); ++j) {
                    ceres::CostFunction *lidar_plane_factor = LidarPlaneFactor::Create(tmpPlanePoints[i].curr_point,
                                                                                  tmpPlanePoints[i].last_point_a,
                                                                                  tmpPlanePoints[i].last_point_b,
                                                                                  tmpPlanePoints[i].last_point_c,
                                                                                  tmpPlanePoints[i].s, s_1);
                    problem.AddResidualBlock(lidar_plane_factor, loss_function, para_Pose[i], para_Pose[i+1], para_Ex_Pose_lidar[0]);
                }
            }
//
//            std::vector<const double *> parameters(4);
//            parameters[0] = para_Pose[i];
//            parameters[1] = para_Ex_Pose_wheel[0];
//            parameters[2] = para_plane_R[0];
//            parameters[3] = para_plane_Z[0];
//            plane_factor->check(const_cast<double **>(parameters.data()));
        }
    }

    int f_m_cnt = 0;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < 4)
            continue;
 
        ++feature_index;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            //每个具有共视关系的帧都会形成因子
            if (imu_i != imu_j)
            {
                Vector3d pts_j = it_per_frame.point;
                ProjectionTwoFrameOneCamFactor *f_td = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);
                //TODO add lidar factor here (only add once)
//                ROS_WARN_STREAM("imu_i " << imu_i);
//                ROS_WARN_STREAM("imu_j " << imu_j);
//                ROS_WARN_STREAM("feature_index " << feature_index);

//                std::vector<const double *> parameters(6);
//                parameters[0] = para_Pose[imu_i];
//                parameters[1] = para_Pose[imu_j];
//                parameters[2] = para_Ex_Pose[0];
//                parameters[3] = para_Feature[feature_index];
//                parameters[4] = para_Td[0];
//                f_td->check(const_cast<double **>(parameters.data()));
            }
            //双目情况，还会与具有共视关系的帧的右目形成因子
            if(STEREO && it_per_frame.is_stereo)
            {                
                Vector3d pts_j_right = it_per_frame.pointRight;
                if(imu_i != imu_j)
                {
                    ProjectionTwoFrameTwoCamFactor *f = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]);
                }
                else
                {
                    ProjectionOneFrameTwoCamFactor *f = new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f, loss_function, para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]);
                }
               
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    //printf("prepare for ceres: %f \n", t_prepare.toc());

    ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
//    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    if (marginalization_flag == MARGIN_OLD)
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    else
        options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    //cout << summary.BriefReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    //printf("solver costs: %f \n", t_solver.toc());

    double2vector();
    //printf("frame_count: %d \n", frame_count);

    if(frame_count < WINDOW_SIZE) //滑窗未满
        return;
    /*滑窗满了，进行边缘化处理*/
    TicToc t_whole_marginalization;
    if (marginalization_flag == MARGIN_OLD)
    {
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        vector2double();
        // 先验部分，基于先验残差，边缘化滑窗中第0帧时刻的状态向量
        if (last_marginalization_info && last_marginalization_info->valid)
        {
            vector<int> drop_set;
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);//先验部分如果存在与要被边缘化变量相同的变量，也要对先验部分进行边缘化
            }
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                           last_marginalization_parameter_blocks,
                                                                           drop_set);
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        //imu 预积分部分，基于第0帧与第1帧之间的预积分残差，边缘化第0帧状态向量
        if(USE_IMU)
        {
            if (pre_integrations[1]->sum_dt < 10.0)
            {
#if !WHEEL
                IMUFactor* imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});//边缘化 para_Pose[0], para_SpeedBias[0]
#else
                IMUWheelLineFactor* imu_factor = new IMUWheelLineFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1], para_Ex_NHC[0]},
                                                                           vector<int>{0, 1});//边缘化 para_Pose[0], para_SpeedBias[0]
#endif
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }
        //wheel 预积分部分，基于第0帧与第1帧之间的预积分残差，边缘化第0帧状态向量
        if(USE_WHEEL && !ONLY_INITIAL_WITH_WHEEL)
        {
            if (pre_integrations_wheel[1]->sum_dt < 10.0)
            {
                WheelFactor* wheel_factor = new WheelFactor(pre_integrations_wheel[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(wheel_factor, NULL,
                                                                               vector<double *>{para_Pose[0], para_Pose[1], para_Ex_Pose_wheel[0], para_Ix_sx_wheel[0], para_Ix_sy_wheel[0], para_Ix_sw_wheel[0], para_Td_wheel[0]},
                                                                               vector<int>{0});//边缘化 para_Pose[0]
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }

//        if(USE_WHEEL && ONLY_INITIAL_WITH_WHEEL)
//        {
//            if (pre_integrations[1]->sum_dt < 10.0)
//            {
//                NonholonomicFactor* nonholonomic_factor = new NonholonomicFactor(pre_integrations[1]);
//                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(nonholonomic_factor, NULL,
//                                                                               vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1], para_Ex_NHC[0]},
//                                                                               vector<int>{0, 1});//边缘化 para_Pose[0], para_SpeedBias[0]
//                marginalization_info->addResidualBlockInfo(residual_block_info);
//            }
//        }

        //平面部分，边缘化第0帧状态向量
        if(USE_PLANE)
        {
            PlaneFactor* plane_factor = new PlaneFactor();
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(plane_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_Ex_Pose_wheel[0], para_plane_R[0], para_plane_Z[0]},
                                                                           vector<int>{0});//边缘化 para_Pose[0]
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        if(static_flag){
            ZUPTFactor* zupt_factor = new ZUPTFactor();
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(zupt_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});//边缘化 para_SpeedBias[0]
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        if(USE_LIDAR)
        {
            if (!EdgePointsFrameVec[0].empty()) {
                vector<EdgePoints> tmpEdgePoints = EdgePointsFrameVec[0];
                for (int j = 0; j < tmpEdgePoints.size(); ++j) {
                    ceres::CostFunction *lidar_edge_factor = LidarEdgeFactor::Create(tmpEdgePoints[0].curr_point,
                                                                                 tmpEdgePoints[0].last_point_a,
                                                                                 tmpEdgePoints[0].last_point_b,
                                                                                 tmpEdgePoints[0].s, s_1);
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(lidar_edge_factor, NULL,
                                                                                   vector<double *>{para_Pose[0], para_Pose[1], para_Ex_Pose_lidar[0]},
                                                                                   vector<int>{0});//边缘化 para_Pose[0]
                    marginalization_info->addResidualBlockInfo(residual_block_info);
                }
                vector<PlanePoints> tmpPlanePoints = PlanePointsFrameVec[0];
                for (int j = 0; j < tmpPlanePoints.size(); ++j) {
                    ceres::CostFunction *lidar_plane_factor = LidarPlaneFactor::Create(tmpPlanePoints[0].curr_point,
                                                                                  tmpPlanePoints[0].last_point_a,
                                                                                  tmpPlanePoints[0].last_point_b,
                                                                                  tmpPlanePoints[0].last_point_c,
                                                                                  tmpPlanePoints[0].s, s_1);
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(lidar_plane_factor, NULL,
                                                                                   vector<double *>{para_Pose[0], para_Pose[1], para_Ex_Pose_lidar[0]},
                                                                                   vector<int>{0});//边缘化 para_Pose[0]
                    marginalization_info->addResidualBlockInfo(residual_block_info);
                }
            }
        }

        //图像部分，基于与第0帧相关的图像残差，边缘化第一次观测的图像帧为第0帧的路标点和第0帧
        {
            int feature_index = -1;
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (it_per_id.used_num < 4)
                    continue;

                ++feature_index;

                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                if (imu_i != 0)
                    continue;

                Vector3d pts_i = it_per_id.feature_per_frame[0].point;//主帧

                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if(imu_i != imu_j)
                    {
                        Vector3d pts_j = it_per_frame.point;
                        ProjectionTwoFrameOneCamFactor *f_td = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, loss_function,
                                                                                        vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                                                                                        vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    if(STEREO && it_per_frame.is_stereo)
                    {
                        Vector3d pts_j_right = it_per_frame.pointRight;
                        if(imu_i != imu_j)
                        {
                            ProjectionTwoFrameTwoCamFactor *f = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                           vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]},
                                                                                           vector<int>{0, 4});
                            marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                        else
                        {
                            ProjectionOneFrameTwoCamFactor *f = new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                           vector<double *>{para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]},
                                                                                           vector<int>{2});
                            marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                    }
                }
            }
        }

        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());
        
        TicToc t_margin;
        marginalization_info->marginalize();
        ROS_DEBUG("marginalization %f ms", t_margin.toc());
        //仅仅改变滑窗double部分地址映射，具体值的改变通过slideWindow和vector2double函数完成；记住边缘化仅仅改变A和b，不改变状态向量
        //由于第0帧观测到的路标点全被边缘化，即边缘化后保存的状态向量中没有路标点;因此addr_shift无需添加路标点
        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)//最老图像帧数据丢弃，从i=1开始遍历
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];// i数据保存到1-1指向的地址，滑窗向前移动一格
            if(USE_IMU)
                addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
        }
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
        addr_shift[reinterpret_cast<long>(para_Ex_Pose_wheel[0])] = para_Ex_Pose_wheel[0];
        addr_shift[reinterpret_cast<long>(para_Ix_sx_wheel[0])] = para_Ix_sx_wheel[0];
        addr_shift[reinterpret_cast<long>(para_Ix_sy_wheel[0])] = para_Ix_sy_wheel[0];
        addr_shift[reinterpret_cast<long>(para_Ix_sw_wheel[0])] = para_Ix_sw_wheel[0];

        addr_shift[reinterpret_cast<long>(para_plane_R[0])] = para_plane_R[0];
        addr_shift[reinterpret_cast<long>(para_plane_Z[0])] = para_plane_Z[0];

        addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
        addr_shift[reinterpret_cast<long>(para_Td_wheel[0])] = para_Td_wheel[0];

        addr_shift[reinterpret_cast<long>(para_Ex_NHC[0])] = para_Ex_NHC[0];

        addr_shift[reinterpret_cast<long>(para_Ex_Pose_lidar[0])] = para_Ex_Pose_lidar[0];

        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = parameter_blocks;
        
    }
    else //将次新的图像帧数据边缘化
    {
        if (last_marginalization_info && //存在先验边缘化信息时才进行次新帧边缘化;否则仅仅通过slidewindow，丢弃次新帧
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();
            if (last_marginalization_info && last_marginalization_info->valid)
            {
                vector<int> drop_set; //记录需要丢弃的变量在last_marginalization_parameter_blocks中的索引
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // construct new marginlization_factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }

            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());
            
            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    if(USE_IMU)
                        addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    if(USE_IMU)
                        addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
            addr_shift[reinterpret_cast<long>(para_Ex_Pose_wheel[0])] = para_Ex_Pose_wheel[0];
            addr_shift[reinterpret_cast<long>(para_Ix_sx_wheel[0])] = para_Ix_sx_wheel[0];
            addr_shift[reinterpret_cast<long>(para_Ix_sy_wheel[0])] = para_Ix_sy_wheel[0];
            addr_shift[reinterpret_cast<long>(para_Ix_sw_wheel[0])] = para_Ix_sw_wheel[0];

            addr_shift[reinterpret_cast<long>(para_plane_R[0])] = para_plane_R[0];
            addr_shift[reinterpret_cast<long>(para_plane_Z[0])] = para_plane_Z[0];
            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
            addr_shift[reinterpret_cast<long>(para_Td_wheel[0])] = para_Td_wheel[0];

            addr_shift[reinterpret_cast<long>(para_Ex_NHC[0])] = para_Ex_NHC[0];

            addr_shift[reinterpret_cast<long>(para_Ex_Pose_lidar[0])] = para_Ex_Pose_lidar[0];


            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;
            
        }
    }
    //printf("whole marginalization costs: %f \n", t_whole_marginalization.toc());
    //printf("whole time for ceres: %f \n", t_whole.toc());
}

void Estimator::slideWindow()
{
    TicToc t_margin;
    if (marginalization_flag == MARGIN_OLD) // 边缘化最老的图像帧，即次新的图像帧为关键帧
    {
        double t_0 = Headers[0];
        back_R0 = Rs[0];
        back_P0 = Ps[0];
        if (frame_count == WINDOW_SIZE) //仅在滑窗满时，进行滑窗边缘化处理
        {
            //1、滑窗中的数据往前移动一帧；运行结果就是WINDOW_SIZE位置的状态为之前0位置对应的状态
            // 0,1,2...WINDOW_SIZE——>1,2...WINDOW_SIZE,0
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Headers[i] = Headers[i + 1];
                Rs[i].swap(Rs[i + 1]);
                Ps[i].swap(Ps[i + 1]);
#if !WHEEL
                if(USE_IMU)
                {
                    std::swap(pre_integrations[i], pre_integrations[i + 1]);

                    dt_buf[i].swap(dt_buf[i + 1]);
                    linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                    angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                    Vs[i].swap(Vs[i + 1]);
                    Bas[i].swap(Bas[i + 1]);
                    Bgs[i].swap(Bgs[i + 1]);
                }
#else
                if(USE_IMU && USE_WHEEL)
                {
                    std::swap(pre_integrations[i], pre_integrations[i + 1]);

                    dt_buf[i].swap(dt_buf[i + 1]);
                    linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                    angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);
                    vel_velocity_buf[i].swap(vel_velocity_buf[i + 1]);

                    Vs[i].swap(Vs[i + 1]);
                    Bas[i].swap(Bas[i + 1]);
                    Bgs[i].swap(Bgs[i + 1]);
                }
#endif
                if(USE_WHEEL)
                {
                    std::swap(pre_integrations_wheel[i], pre_integrations_wheel[i + 1]);

                    dt_buf_wheel[i].swap(dt_buf_wheel[i + 1]);
                    linear_velocity_buf_wheel[i].swap(linear_velocity_buf_wheel[i + 1]);
                    angular_velocity_buf_wheel[i].swap(angular_velocity_buf_wheel[i + 1]);
                }
            }
            //2、处理前，WINDOW_SIZE位置的状态为之前0位置对应的状态；处理后，WINDOW_SIZE位置的状态为之前WINDOW_SIZE位置对应的状态;之前0位置对应的状态被剔除
            // 0,1,2...WINDOW_SIZE——>1,2...WINDOW_SIZE,WINDOW_SIZE
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];

#if WHEEL
            if(USE_IMU && USE_WHEEL)
            {
                Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
                Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
                Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, vel_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
                vel_velocity_buf[WINDOW_SIZE].clear();
            }
#else
            if(USE_IMU)
            {
                Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
                Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
                Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
            }
#endif

            if(USE_WHEEL)
            {

                delete pre_integrations_wheel[WINDOW_SIZE];
                pre_integrations_wheel[WINDOW_SIZE] = new WheelIntegrationBase{vel_0_wheel, gyr_0_wheel, sx, sy, sw, td_wheel};

                dt_buf_wheel[WINDOW_SIZE].clear();
                linear_velocity_buf_wheel[WINDOW_SIZE].clear();
                angular_velocity_buf_wheel[WINDOW_SIZE].clear();
            }

            //3、对时刻t_0(对应滑窗第0帧)之前的所有数据进行剔除；即all_image_frame中仅保留滑窗中图像帧0与图像帧WINDOW_SIZE之间的数据
            if (true || solver_flag == INITIAL)
            {
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                delete it_0->second.pre_integration_wheel;
                all_image_frame.erase(all_image_frame.begin(), it_0);
            }
            slideWindowOld();
        }
    }
    else //边缘化次新的图像帧，主要完成的工作是数据衔接
    {
        if (frame_count == WINDOW_SIZE) //仅在滑窗满时，进行滑窗边缘化处理
        {//0,1,2...WINDOW_SIZE-2, WINDOW_SIZE-1, WINDOW_SIZE——>0,,1,2...WINDOW_SIZE-2,WINDOW_SIZE, WINDOW_SIZE
            Headers[frame_count - 1] = Headers[frame_count];
            Ps[frame_count - 1] = Ps[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];
#if !WHEEL
            if(USE_IMU) //IMU数据衔接，预积分的传播
            {
                for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
                {
                    double tmp_dt = dt_buf[frame_count][i];
                    Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                    Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];

                    pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity); //预积分的传播

                    dt_buf[frame_count - 1].push_back(tmp_dt); //TODO(tzhang): 数据保存有冗余，integration_base中也保存了同样的数据
                    linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                    angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
                }

                Vs[frame_count - 1] = Vs[frame_count];
                Bas[frame_count - 1] = Bas[frame_count];
                Bgs[frame_count - 1] = Bgs[frame_count];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
            }
#else
            if(USE_IMU && USE_WHEEL) //IMU数据衔接，预积分的传播
            {
                for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
                {
                    double tmp_dt = dt_buf[frame_count][i];
                    Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                    Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];
                    Vector3d tmp_wheel_velocity = vel_velocity_buf[frame_count][i];

//                    ROS_WARN_STREAM("tmp_dt " << tmp_dt);
                    pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity, tmp_wheel_velocity); //预积分的传播

                    dt_buf[frame_count - 1].push_back(tmp_dt); //TODO(tzhang): 数据保存有冗余，integration_base中也保存了同样的数据
                    linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                    angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
                    vel_velocity_buf[frame_count - 1].push_back(tmp_wheel_velocity);
                }

                Vs[frame_count - 1] = Vs[frame_count];
                Bas[frame_count - 1] = Bas[frame_count];
                Bgs[frame_count - 1] = Bgs[frame_count];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, vel_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
                vel_velocity_buf[WINDOW_SIZE].clear();
            }
#endif
            if(USE_WHEEL) //WHEEL数据衔接，预积分的传播
            {
                for (unsigned int i = 0; i < dt_buf_wheel[frame_count].size(); i++)
                {
                    double tmp_dt = dt_buf_wheel[frame_count][i];
                    Vector3d tmp_linear_velocity = linear_velocity_buf_wheel[frame_count][i];
                    Vector3d tmp_angular_velocity = angular_velocity_buf_wheel[frame_count][i];

                    pre_integrations_wheel[frame_count - 1]->push_back(tmp_dt, tmp_linear_velocity, tmp_angular_velocity); //预积分的传播

                    dt_buf_wheel[frame_count - 1].push_back(tmp_dt); //TODO(tzhang): 数据保存有冗余，integration_base中也保存了同样的数据
                    linear_velocity_buf_wheel[frame_count - 1].push_back(tmp_linear_velocity);
                    angular_velocity_buf_wheel[frame_count - 1].push_back(tmp_angular_velocity);
                }

                delete pre_integrations_wheel[WINDOW_SIZE];
                pre_integrations_wheel[WINDOW_SIZE] = new WheelIntegrationBase{vel_0_wheel, gyr_0_wheel, sx, sy, sw, td_wheel};

                dt_buf_wheel[WINDOW_SIZE].clear();
                linear_velocity_buf_wheel[WINDOW_SIZE].clear();
                angular_velocity_buf_wheel[WINDOW_SIZE].clear();
            }
            slideWindowNew(); //更新第一次观测到路标点的图像帧的索引
        }
    }
}

void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}

void Estimator::slideWindowOld()
{
    sum_of_back++;

    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth) //非线性优化阶段，除了更新第一次观测到路标点的图像帧的索引，还需更新路标点的逆深度
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        R0 = back_R0 * ric[0]; //R_w_cl  老的第0帧左相机坐标系与世界坐标系之间的相对旋转
        R1 = Rs[0] * ric[0];   //R_w_cl  新的第0帧左相机坐标系与世界坐标系之间的相对旋转
        P0 = back_P0 + back_R0 * tic[0];
        P1 = Ps[0] + Rs[0] * tic[0];
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();  //初始化未完成，只是更新第一次观测到路标点的图像帧的索引
}


void Estimator::getPoseInWorldFrame(Eigen::Matrix4d &T)
{
    T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = Rs[frame_count];
    T.block<3, 1>(0, 3) = Ps[frame_count];
}

void Estimator::getPoseInWorldFrame(int index, Eigen::Matrix4d &T)
{
    T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = Rs[index];
    T.block<3, 1>(0, 3) = Ps[index];
}

void Estimator::predictPtsInNextFrame()
{
    //printf("predict pts in next frame\n");
    if(frame_count < 2)
        return;
    // predict next pose. Assume constant velocity motion
    Eigen::Matrix4d curT, prevT, nextT;
    getPoseInWorldFrame(curT);
    getPoseInWorldFrame(frame_count - 1, prevT);
    nextT = curT * (prevT.inverse() * curT);
    map<int, Eigen::Vector3d> predictPts;

    for (auto &it_per_id : f_manager.feature)
    {
        if(it_per_id.estimated_depth > 0)
        {
            int firstIndex = it_per_id.start_frame;
            int lastIndex = it_per_id.start_frame + it_per_id.feature_per_frame.size() - 1;
            //printf("cur frame index  %d last frame index %d\n", frame_count, lastIndex);
            if((int)it_per_id.feature_per_frame.size() >= 2 && lastIndex == frame_count)
            {
                double depth = it_per_id.estimated_depth;
                Vector3d pts_j = ric[0] * (depth * it_per_id.feature_per_frame[0].point) + tic[0];
                Vector3d pts_w = Rs[firstIndex] * pts_j + Ps[firstIndex];
                Vector3d pts_local = nextT.block<3, 3>(0, 0).transpose() * (pts_w - nextT.block<3, 1>(0, 3));
                Vector3d pts_cam = ric[0].transpose() * (pts_local - tic[0]);
                int ptsIndex = it_per_id.feature_id;
                predictPts[ptsIndex] = pts_cam;
            }
        }
    }
    featureTracker.setPrediction(predictPts);
    //printf("estimator output %d predict pts\n",(int)predictPts.size());
}

double Estimator::reprojectionError(Matrix3d &Ri, Vector3d &Pi, Matrix3d &rici, Vector3d &tici,
                                 Matrix3d &Rj, Vector3d &Pj, Matrix3d &ricj, Vector3d &ticj, 
                                 double depth, Vector3d &uvi, Vector3d &uvj)
{
    Vector3d pts_w = Ri * (rici * (depth * uvi) + tici) + Pi;
    Vector3d pts_cj = ricj.transpose() * (Rj.transpose() * (pts_w - Pj) - ticj);
    Vector2d residual = (pts_cj / pts_cj.z()).head<2>() - uvj.head<2>();
    double rx = residual.x();
    double ry = residual.y();
    return sqrt(rx * rx + ry * ry);
}

void Estimator::outliersRejection(set<int> &removeIndex)
{
    //return;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature) // 遍历所有路标点
    {
        double err = 0;
        int errCnt = 0;
        it_per_id.used_num = it_per_id.feature_per_frame.size();//即被观察到该路标点的图像帧数目
        if (it_per_id.used_num < 4) // 对观测少于4次的路标点，不进行外点判断
            continue;
        feature_index ++;
        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;
        double depth = it_per_id.estimated_depth;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i != imu_j) //不同时刻，左相机在不同帧之间的重投影误差计算
            {
                Vector3d pts_j = it_per_frame.point;             
                double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                    Rs[imu_j], Ps[imu_j], ric[0], tic[0],
                                                    depth, pts_i, pts_j);
                err += tmp_error;
                errCnt++;
                //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
            }
            // need to rewrite projecton factor.........
            if(STEREO && it_per_frame.is_stereo) // 双目情形
            {
                
                Vector3d pts_j_right = it_per_frame.pointRight;
                if(imu_i != imu_j) //不同时刻，左右图像帧之间的重投影误差
                {            
                    double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                        Rs[imu_j], Ps[imu_j], ric[1], tic[1],
                                                        depth, pts_i, pts_j_right);
                    err += tmp_error;
                    errCnt++;
                    //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
                }
                else //相同时刻，左右图像帧之间的重投影误差 TODO(tzhang)：此处不同时刻判断没啥用，代码冗余
                {
                    double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                        Rs[imu_j], Ps[imu_j], ric[1], tic[1],
                                                        depth, pts_i, pts_j_right);
                    err += tmp_error;
                    errCnt++;
                    //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
                }       
            }
        }
        double ave_err = err / errCnt;
        if(ave_err * FOCAL_LENGTH > 3) // 若平均的重投影均方根过大，则判定该路标点为外点; 添加该路标点编号至removeIndex中
            removeIndex.insert(it_per_id.feature_id);

    }
}
//中值积分IMU惯性解算
void Estimator::fastPredictIMU(double t, Eigen::Vector3d linear_acceleration, Eigen::Vector3d angular_velocity)
{
    double dt = t - latest_time;
    latest_time = t;
    //latest_P, latest_Q, latest_V的初值将会在视觉初始化中从updateLatestStates()中获得
    Eigen::Vector3d un_acc_0 = latest_Q * (latest_acc_0 - latest_Ba) - g;
    Eigen::Vector3d un_gyr = 0.5 * (latest_gyr_0 + angular_velocity) - latest_Bg;
    latest_Q = latest_Q * Utility::deltaQ(un_gyr * dt);
    Eigen::Vector3d un_acc_1 = latest_Q * (linear_acceleration - latest_Ba) - g;
    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    latest_P = latest_P + dt * latest_V + 0.5 * dt * dt * un_acc;
    latest_V = latest_V + dt * un_acc;
    latest_acc_0 = linear_acceleration;
    latest_gyr_0 = angular_velocity;
}

//中值积分Wheel航迹解算
void Estimator::fastPredictWheel(double t, Eigen::Vector3d linear_velocity, Eigen::Vector3d angular_velocity)
{
    double dt = t - latest_time_wheel;
    latest_time_wheel = t;
    Eigen::Vector3d un_gyr = 0.5 * latest_sw * (latest_gyr_0 + angular_velocity);
    Eigen::Vector3d un_vel_0 = latest_Q_wheel * latest_vel_wheel_0;
    Eigen::Matrix3d latest_sv = Eigen::Vector3d(latest_sx, latest_sy, 1).asDiagonal();
    //在初始化完成后会进行updateLatestState，对这些latest量进行更新
    latest_Q_wheel = latest_Q_wheel * Utility::deltaQ(un_gyr * dt);
    latest_V_wheel = 0.5 * latest_sv * (latest_Q_wheel * linear_velocity + un_vel_0);
    latest_P_wheel = latest_P_wheel + dt * latest_V_wheel ;
    latest_vel_wheel_0 = linear_velocity;
    latest_gyr_wheel_0 = angular_velocity;

}
//中值积分Wheel航迹解算
void Estimator::fastPredictPureWheel(double t, Eigen::Vector3d linear_velocity, Eigen::Vector3d angular_velocity, Eigen::Vector3d &P, Eigen::Quaterniond &Q, Eigen::Vector3d &V)
{
    static bool first_time = false;
    static Eigen::Quaterniond Q_latest;
    static Eigen::Vector3d V_latest = Eigen::Vector3d::Zero();
    static Eigen::Vector3d P_latest = Eigen::Vector3d::Zero();
    static Eigen::Vector3d vel_latest_0 = Eigen::Vector3d::Zero();
    static Eigen::Vector3d gyr_latest_0 = Eigen::Vector3d::Zero();
    static double t_latest;
    if(!first_time){
        first_time = true;
        Q_latest = latest_Q_wheel;
        V_latest = latest_V_wheel;
        P_latest = latest_P_wheel;
        vel_latest_0 = latest_vel_wheel_0;
        gyr_latest_0 = latest_gyr_wheel_0;
        t_latest = latest_time_wheel;
        std::cout<<"fastPredictPureWheel initial pose: \n"<<P_latest.transpose()<<std::endl<<Q_latest.coeffs().transpose()<<std::endl;
    }

    double dt = t - t_latest;
    t_latest = t;
    Eigen::Vector3d un_gyr = 0.5 * latest_sw * (gyr_latest_0 + angular_velocity);
    Eigen::Vector3d un_vel_0 = Q_latest * vel_latest_0;
    Eigen::Matrix3d latest_sv = Eigen::Vector3d(latest_sx, latest_sy, 1).asDiagonal();
    //在初始化完成后会进行updateLatestState，对这些latest量进行更新
    Q_latest = Q_latest * Utility::deltaQ(un_gyr * dt);
    V_latest = 0.5 * latest_sv * (Q_latest * linear_velocity + un_vel_0);
    P_latest = P_latest + dt * V_latest ;
    vel_latest_0 = linear_velocity;
    gyr_latest_0 = angular_velocity;

    P = P_latest;
    Q = Q_latest;
    V = V_latest;
}
//获取滑窗中最新帧时刻的状态，并在世界坐标系下进行中值积分；初始化完成后，最新状态在inputIMU函数中发布
void Estimator::updateLatestStates()
{
    //IMU
    mPropagate.lock();
    latest_time = Headers[frame_count] + td;
    latest_P = Ps[frame_count];
    latest_Q = Rs[frame_count];
    latest_V = Vs[frame_count];
    latest_Ba = Bas[frame_count];
    latest_Bg = Bgs[frame_count];
    latest_acc_0 = acc_0;
    latest_gyr_0 = gyr_0;
    mBuf.lock();
    queue<pair<double, Eigen::Vector3d>> tmp_accBuf = accBuf;
    queue<pair<double, Eigen::Vector3d>> tmp_gyrBuf = gyrBuf;
    mBuf.unlock();
    while(!tmp_accBuf.empty())
    {
        double t = tmp_accBuf.front().first;
        Eigen::Vector3d acc = tmp_accBuf.front().second;
        Eigen::Vector3d gyr = tmp_gyrBuf.front().second;
        fastPredictIMU(t, acc, gyr);
        tmp_accBuf.pop();
        tmp_gyrBuf.pop();
    }

    if(ESTIMATE_RIO) {
//        Matrix3d tmp_R = Utility::ypr2R(Vector3d{Utility::R2ypr(dR).x(), 0 , 0});
        Matrix3d tmp_R = dR;
//        rio_0 = (tmp_R * R0 * RIC[0].transpose()).transpose();
        rio_0 = dR.transpose();
        yaw_test = Utility::R2ypr(rio_0).x();
        pitch_test = Utility::R2ypr(rio_0).y();
        roll_test = Utility::R2ypr(rio_0).z();
        tmpR_test = Utility::R2ypr(tmp_R).x();
        ROS_WARN_STREAM("rio  ypr   " << Utility::R2ypr_m(rio_0).transpose());
        ROS_WARN_STREAM("rio     " << endl << rio_0);
        ROS_WARN_STREAM("tmpR ypr     " << endl << Utility::R2ypr_m(tmp_R).transpose());
//        yaw_sum_vec.push_back(Utility::R2ypr_m(tmp_R).norm());
//        if (yaw_sum_vec.size() == 100) {
//            double sum = std::accumulate(std::begin(yaw_sum_vec), std::end(yaw_sum_vec), 0.0);
//            double mean = sum / yaw_sum_vec.size();
//
//            double accum = 0.0;
//            std::for_each(std::begin(yaw_sum_vec), std::end(yaw_sum_vec), [&](const double d) {
//                accum += (d - mean) * (d - mean);
//            });
//
//            double stdev = sqrt(accum / (yaw_sum_vec.size() - 1));
////            x_test = stdev;
//            ROS_WARN_STREAM("std rotation: " << stdev);
//            if (stdev < 0.05) {
//                ESTIMATE_RIO = 0;
//            } else
//                yaw_sum_vec.pop_front();
//        }
    }else{
        ROS_WARN_STREAM("RIO Calibration finished! " << Utility::R2ypr_m(rio_0).transpose());
        ROS_WARN_STREAM("RIO Calibration  " << endl << rio_0);
    }

    if(ESTIMATE_TIO) {
        x_test = tio_0.x();
        y_test = tio_0.y();
        z_test = tio_0.z();
        ROS_WARN_STREAM("tio     " << tio_0.transpose());
//        xy_sum_vec.push_back(tio_0.norm());
//        if (xy_sum_vec.size() == 100) {
//            double sum = std::accumulate(std::begin(xy_sum_vec), std::end(xy_sum_vec), 0.0);
//            double mean = sum / xy_sum_vec.size();
//
//            double accum = 0.0;
//            std::for_each(std::begin(xy_sum_vec), std::end(xy_sum_vec), [&](const double d) {
//                accum += (d - mean) * (d - mean);
//            });
//
//            double stdev = sqrt(accum / (xy_sum_vec.size() - 1));
////            y_test = stdev;
//            ROS_WARN_STREAM("std translation: " << stdev);
//            if (stdev < 0.05) {
//                ESTIMATE_TIO = 0;
//            } else
//                xy_sum_vec.pop_front();
//        }
    }else{
        ROS_WARN_STREAM("TIO Calibration finished! " << tio_0.transpose());
    }


    mPropagate.unlock();

    //my code
    //Wheel
    mWheelPropagate.lock();
    latest_time_wheel = Headers[frame_count] + td - td_wheel;
    latest_Q_wheel = Rs[frame_count] * RIO;
    latest_P_wheel = Rs[frame_count] * TIO + Ps[frame_count];
    latest_sx = sx;
    latest_sy = sy;
    latest_sw = sw;
    //latest_V_wheel可以在fastPredictWheel算出来，不需要初值
//    latest_V_wheel = Vs[frame_count];
    latest_vel_wheel_0 = vel_0_wheel;

    latest_gyr_wheel_0 = gyr_0_wheel;
    mWheelBuf.lock();
    queue<pair<double, Eigen::Vector3d>> tmp_wheel_velBuf = wheelVelBuf;
    queue<pair<double, Eigen::Vector3d>> tmp_wheel_gyrBuf = wheelGyrBuf;
    mWheelBuf.unlock();
    //这里更新了latest_p v q，因此要buf里剩余的测量值都进行航迹推算，这样才能保证下次测量数据输入时inputWheel，是使用的最新的状态进行航迹推算
    while(!tmp_wheel_velBuf.empty())
    {
        double t = tmp_wheel_velBuf.front().first;
        Eigen::Vector3d vel = tmp_wheel_velBuf.front().second;
        Eigen::Vector3d gyr = tmp_wheel_gyrBuf.front().second;
        fastPredictWheel(t, vel, gyr);
        tmp_wheel_velBuf.pop();
        tmp_wheel_gyrBuf.pop();
    }
    mWheelPropagate.unlock();
}

// undistort lidar point
void Estimator::TransformToStart(PointType const *const pi, PointType *const po)
{
    //interpolation ratio
    double s;
    if (DISTORTION)
        s = (pi->intensity - int(pi->intensity)) / SCAN_PERIOD;
    else
        s = 1.0;
    //s = 1;
    Eigen::Quaterniond q_point_last = Eigen::Quaterniond::Identity().slerp(s, q_last_curr);
    Eigen::Vector3d t_point_last = s * t_last_curr;
    Eigen::Vector3d point(pi->x, pi->y, pi->z);
    Eigen::Vector3d un_point = q_point_last * point + t_point_last;

    po->x = un_point.x();
    po->y = un_point.y();
    po->z = un_point.z();
    po->intensity = pi->intensity;
}

// transform all lidar points to the start of the next frame

void Estimator::TransformToAnyTime(PointType const *const pi, PointType *const po, double s)
{

    // undistort point first
    pcl::PointXYZI un_point_tmp;
    TransformToStart(pi, &un_point_tmp);

    Eigen::Vector3d un_point{un_point_tmp.x, un_point_tmp.y, un_point_tmp.z};
    Eigen::Vector3d point_end = q_last_curr.inverse() * (un_point - t_last_curr);

    //transform points from the end of the current frame to the next closest image frame  --by zhh
    Eigen::Quaterniond q_next_image_point = Eigen::Quaterniond::Identity().slerp(s, q_last_curr); // constant velocity assumption
    Eigen::Vector3d t_next_image_point = s * t_last_curr;

    point_end = q_next_image_point.inverse() * (point_end - t_next_image_point);



#if USE_IMU_FOR_DESKEW
    Eigen::Vector3d imuShift{imuShiftFromStartX, imuShiftFromStartY, imuShiftFromStartZ};

    Eigen::Quaterniond imuQstart = Eigen::AngleAxisd(imuYawStart, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(imuPitchStart, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(imuRollStart, Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond imuQlast = Eigen::AngleAxisd(-imuRollLast, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(-imuPitchLast, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-imuYawLast, Eigen::Vector3d::UnitY());

//    Eigen::Quaterniond imuQstart = Eigen::AngleAxisd(imuRollStart, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(imuPitchStart, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(imuYawStart, Eigen::Vector3d::UnitY());
//    Eigen::Quaterniond imuQlast = Eigen::AngleAxisd(-imuYawLast, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(-imuPitchLast, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-imuRollLast, Eigen::Vector3d::UnitZ());

    point_end = imuQlast * imuQstart * (point_end - imuShift);
#endif

    po->x = point_end.x();
    po->y = point_end.y();
    po->z = point_end.z();

    //Remove distortion time info
    po->intensity = int(pi->intensity);
}

void Estimator::TransformToEnd(PointType const *const pi, PointType *const po)
{
    // undistort point first
    pcl::PointXYZI un_point_tmp;
    TransformToStart(pi, &un_point_tmp);

    Eigen::Vector3d un_point(un_point_tmp.x, un_point_tmp.y, un_point_tmp.z);
    Eigen::Vector3d point_end = q_last_curr.inverse() * (un_point - t_last_curr);

#if USE_IMU_FOR_DESKEW
    Eigen::Vector3d imuShift{imuShiftFromStartX, imuShiftFromStartY, imuShiftFromStartZ};

    Eigen::Quaterniond imuQstart = Eigen::AngleAxisd(imuYawStart, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(imuPitchStart, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(imuRollStart, Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond imuQlast = Eigen::AngleAxisd(-imuRollLast, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(-imuPitchLast, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-imuYawLast, Eigen::Vector3d::UnitY());

//    Eigen::Quaterniond imuQstart = Eigen::AngleAxisd(imuRollStart, Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(imuPitchStart, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(imuYawStart, Eigen::Vector3d::UnitY());
//    Eigen::Quaterniond imuQlast = Eigen::AngleAxisd(-imuYawLast, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(-imuPitchLast, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-imuRollLast, Eigen::Vector3d::UnitZ());

    point_end = imuQlast * imuQstart * (point_end - imuShift);
#endif

    po->x = point_end.x();
    po->y = point_end.y();
    po->z = point_end.z();

    //Remove distortion time info
    po->intensity = int(pi->intensity);
}

void Estimator::LaserOdometryProcess(double t0, double t1){
    int cornerPointsSharpNum = cornerPointsSharp->points.size();
    int surfPointsFlatNum = surfPointsFlat->points.size();

    pcl::PointXYZI pointSel;
    std::vector<int> pointSearchInd;
    std::vector<float> pointSearchSqDis;

    // initializing
    if (!systemInited)
    {
        systemInited = true;
        std::cout << "Initialization finished \n";
    }
    else {
        EdgePointsSum.clear();
        PlanePointsSum.clear();
        // transform last points to current frame (match at the time of current frame)
        int cornerPointsLastNum = laserCloudCornerLast->points.size();
        for (int i = 0; i < cornerPointsLastNum; i++)
        {
            TransformToEnd(&laserCloudCornerLast->points[i], &laserCloudCornerLast->points[i]);
        }

        int surfPointsLastNum = laserCloudSurfLast->points.size();
        for (int i = 0; i < surfPointsLastNum; i++)
        {
            TransformToEnd(&laserCloudSurfLast->points[i], &laserCloudSurfLast->points[i]);
        }
        // find correspondence for corner features
        for (int i = 0; i < cornerPointsSharpNum; ++i) {
            double pointTime = prevTime_lidar + cornerPointsSharp->points[i].intensity;
            if (pointTime < t0 || pointTime > t1)
                continue;
            TransformToStart(&(cornerPointsSharp->points[i]), &pointSel);
            // transform to the end for matching on the end
            TransformToEnd(&pointSel, &pointSel);
            kdtreeCornerLast->nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);

            int closestPointInd = -1, minPointInd2 = -1;
            if (pointSearchSqDis[0] < DISTANCE_SQ_THRESHOLD) {
                closestPointInd = pointSearchInd[0];
                int closestPointScanID = int(laserCloudCornerLast->points[closestPointInd].intensity);

                double minPointSqDis2 = DISTANCE_SQ_THRESHOLD;
                // search in the direction of increasing scan line
                for (int j = closestPointInd + 1; j < (int) laserCloudCornerLast->points.size(); ++j) {
                    // if in the same scan line, continue
                    if (int(laserCloudCornerLast->points[j].intensity) <= closestPointScanID)
                        continue;

                    // if not in nearby scans, end the loop
                    if (int(laserCloudCornerLast->points[j].intensity) > (closestPointScanID + NEARBY_SCAN))
                        break;

                    double pointSqDis = (laserCloudCornerLast->points[j].x - pointSel.x) *
                                        (laserCloudCornerLast->points[j].x - pointSel.x) +
                                        (laserCloudCornerLast->points[j].y - pointSel.y) *
                                        (laserCloudCornerLast->points[j].y - pointSel.y) +
                                        (laserCloudCornerLast->points[j].z - pointSel.z) *
                                        (laserCloudCornerLast->points[j].z - pointSel.z);

                    if (pointSqDis < minPointSqDis2) {
                        // find nearer point
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                }

                // search in the direction of decreasing scan line
                for (int j = closestPointInd - 1; j >= 0; --j) {
                    // if in the same scan line, continue
                    if (int(laserCloudCornerLast->points[j].intensity) >= closestPointScanID)
                        continue;

                    // if not in nearby scans, end the loop
                    if (int(laserCloudCornerLast->points[j].intensity) < (closestPointScanID - NEARBY_SCAN))
                        break;

                    double pointSqDis = (laserCloudCornerLast->points[j].x - pointSel.x) *
                                        (laserCloudCornerLast->points[j].x - pointSel.x) +
                                        (laserCloudCornerLast->points[j].y - pointSel.y) *
                                        (laserCloudCornerLast->points[j].y - pointSel.y) +
                                        (laserCloudCornerLast->points[j].z - pointSel.z) *
                                        (laserCloudCornerLast->points[j].z - pointSel.z);

                    if (pointSqDis < minPointSqDis2) {
                        // find nearer point
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                }
            }
            if (minPointInd2 >= 0) // both closestPointInd and minPointInd2 is valid
            {
                Eigen::Vector3d curr_point(cornerPointsSharp->points[i].x,
                                           cornerPointsSharp->points[i].y,
                                           cornerPointsSharp->points[i].z);
                Eigen::Vector3d last_point_a(laserCloudCornerLast->points[closestPointInd].x,
                                             laserCloudCornerLast->points[closestPointInd].y,
                                             laserCloudCornerLast->points[closestPointInd].z);
                Eigen::Vector3d last_point_b(laserCloudCornerLast->points[minPointInd2].x,
                                             laserCloudCornerLast->points[minPointInd2].y,
                                             laserCloudCornerLast->points[minPointInd2].z);
                double s;
                if (DISTORTION)
                    s = (cornerPointsSharp->points[i].intensity - int(cornerPointsSharp->points[i].intensity)) /
                        SCAN_PERIOD;
                else
                    s = 1.0;
                EdgePoints edgePoints = {curr_point, last_point_a, last_point_b, s};
                EdgePointsSum.push_back(edgePoints);
//                ceres::CostFunction *cost_function = LidarEdgeFactor::Create(curr_point, last_point_a, last_point_b, s);
//                problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
            }
        }

        // find correspondence for plane features
        for (int i = 0; i < surfPointsFlatNum; ++i) {
            TransformToStart(&(surfPointsFlat->points[i]), &pointSel);
            kdtreeSurfLast->nearestKSearch(pointSel, 1, pointSearchInd, pointSearchSqDis);

            int closestPointInd = -1, minPointInd2 = -1, minPointInd3 = -1;
            if (pointSearchSqDis[0] < DISTANCE_SQ_THRESHOLD) {
                closestPointInd = pointSearchInd[0];

                // get closest point's scan ID
                int closestPointScanID = int(laserCloudSurfLast->points[closestPointInd].intensity);
                double minPointSqDis2 = DISTANCE_SQ_THRESHOLD, minPointSqDis3 = DISTANCE_SQ_THRESHOLD;

                // search in the direction of increasing scan line
                for (int j = closestPointInd + 1; j < (int) laserCloudSurfLast->points.size(); ++j) {
                    // if not in nearby scans, end the loop
                    if (int(laserCloudSurfLast->points[j].intensity) > (closestPointScanID + NEARBY_SCAN))
                        break;

                    double pointSqDis = (laserCloudSurfLast->points[j].x - pointSel.x) *
                                        (laserCloudSurfLast->points[j].x - pointSel.x) +
                                        (laserCloudSurfLast->points[j].y - pointSel.y) *
                                        (laserCloudSurfLast->points[j].y - pointSel.y) +
                                        (laserCloudSurfLast->points[j].z - pointSel.z) *
                                        (laserCloudSurfLast->points[j].z - pointSel.z);

                    // if in the same or lower scan line
                    if (int(laserCloudSurfLast->points[j].intensity) <= closestPointScanID &&
                        pointSqDis < minPointSqDis2) {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    }
                        // if in the higher scan line
                    else if (int(laserCloudSurfLast->points[j].intensity) > closestPointScanID &&
                             pointSqDis < minPointSqDis3) {
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }

                // search in the direction of decreasing scan line
                for (int j = closestPointInd - 1; j >= 0; --j) {
                    // if not in nearby scans, end the loop
                    if (int(laserCloudSurfLast->points[j].intensity) < (closestPointScanID - NEARBY_SCAN))
                        break;

                    double pointSqDis = (laserCloudSurfLast->points[j].x - pointSel.x) *
                                        (laserCloudSurfLast->points[j].x - pointSel.x) +
                                        (laserCloudSurfLast->points[j].y - pointSel.y) *
                                        (laserCloudSurfLast->points[j].y - pointSel.y) +
                                        (laserCloudSurfLast->points[j].z - pointSel.z) *
                                        (laserCloudSurfLast->points[j].z - pointSel.z);

                    // if in the same or higher scan line
                    if (int(laserCloudSurfLast->points[j].intensity) >= closestPointScanID &&
                        pointSqDis < minPointSqDis2) {
                        minPointSqDis2 = pointSqDis;
                        minPointInd2 = j;
                    } else if (int(laserCloudSurfLast->points[j].intensity) < closestPointScanID &&
                               pointSqDis < minPointSqDis3) {
                        // find nearer point
                        minPointSqDis3 = pointSqDis;
                        minPointInd3 = j;
                    }
                }

                if (minPointInd2 >= 0 && minPointInd3 >= 0) {

                    Eigen::Vector3d curr_point(surfPointsFlat->points[i].x,
                                               surfPointsFlat->points[i].y,
                                               surfPointsFlat->points[i].z);
                    Eigen::Vector3d last_point_a(laserCloudSurfLast->points[closestPointInd].x,
                                                 laserCloudSurfLast->points[closestPointInd].y,
                                                 laserCloudSurfLast->points[closestPointInd].z);
                    Eigen::Vector3d last_point_b(laserCloudSurfLast->points[minPointInd2].x,
                                                 laserCloudSurfLast->points[minPointInd2].y,
                                                 laserCloudSurfLast->points[minPointInd2].z);
                    Eigen::Vector3d last_point_c(laserCloudSurfLast->points[minPointInd3].x,
                                                 laserCloudSurfLast->points[minPointInd3].y,
                                                 laserCloudSurfLast->points[minPointInd3].z);
                    double s;
                    if (DISTORTION)
                        s = (surfPointsFlat->points[i].intensity - int(surfPointsFlat->points[i].intensity)) /
                            SCAN_PERIOD;
                    else
                        s = 1.0;
                    PlanePoints planePoints = {curr_point, last_point_a, last_point_b, last_point_c, s};
                    PlanePointsSum.push_back(planePoints);
//                    ceres::CostFunction *cost_function = LidarPlaneFactor::Create(curr_point, last_point_a, last_point_b, last_point_c, s);
//                    problem.AddResidualBlock(cost_function, loss_function, para_q, para_t);
                }
            }
        }
    }

    int cornerPointsLessSharpNum = cornerPointsLessSharp->points.size();
    for (int i = 0; i < cornerPointsLessSharpNum; i++)
    {
        TransformToEnd(&cornerPointsLessSharp->points[i], &cornerPointsLessSharp->points[i]);
    }

    int surfPointsLessFlatNum = surfPointsLessFlat->points.size();
    for (int i = 0; i < surfPointsLessFlatNum; i++)
    {
        TransformToEnd(&surfPointsLessFlat->points[i], &surfPointsLessFlat->points[i]);
    }

    int laserCloudFullResNum = laserCloudFullRes->points.size();
    for (int i = 0; i < laserCloudFullResNum; i++)
    {
        TransformToEnd(&laserCloudFullRes->points[i], &laserCloudFullRes->points[i]);
    }

    pcl::PointCloud<PointType>::Ptr laserCloudTemp = cornerPointsLessSharp;
    cornerPointsLessSharp = laserCloudCornerLast;
    laserCloudCornerLast = laserCloudTemp;

    laserCloudTemp = surfPointsLessFlat;
    surfPointsLessFlat = laserCloudSurfLast;
    laserCloudSurfLast = laserCloudTemp;

    kdtreeCornerLast->setInputCloud(laserCloudCornerLast);
    kdtreeSurfLast->setInputCloud(laserCloudSurfLast);

}
