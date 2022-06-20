//
// Created by zhh2005757 on 22-6-13.
//

#ifndef VINS_FEATURE_EXTRACTOR_H
#define VINS_FEATURE_EXTRACTOR_H

#pragma once
#include "../estimator/parameters.h"
#include "../utility/tic_toc.h"

class FeatureExtractor {
public:
    void processPointCLoud(const pcl::PointCloud<pcl::PointXYZ> &laserCloudIn);

    int N_SCANS = 16;
    const double scanPeriod = 0.1;
//    static float cloudCurvature[400000];
    int cloudSortInd[400000];
    int cloudNeighborPicked[400000];
    int cloudLabel[400000];
//    bool comp (int i,int j) { return (cloudCurvature[i]<cloudCurvature[j]); }
//    struct {
//        bool operator()(int a, int b) const
//        {
//            return (cloudCurvature[a]<cloudCurvature[b]);
//        }
//    } comp;

    pcl::PointCloud<PointType> cornerPointsSharp;
    pcl::PointCloud<PointType> cornerPointsLessSharp;
    pcl::PointCloud<PointType> surfPointsFlat;
    pcl::PointCloud<PointType> surfPointsLessFlat;
    pcl::PointCloud<PointType>::Ptr laserCloud;
};



#endif //VINS_FEATURE_EXTRACTOR_H
