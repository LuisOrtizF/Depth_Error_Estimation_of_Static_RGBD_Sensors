//standard includes
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>

//Opencv includes
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>

//PCL includes
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/registration/icp.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <boost/thread/thread.hpp>
#include <pcl/console/time.h>

//libcbdetect
#include "rt_nonfinite.h"
#include "fce.h"

using namespace std;
using namespace cv;

static void CornersViewSave();
static void CreateIdealCLoud(int);
static void FindCorners(cv::Mat);
static void RealCloud_Reg_RMS();
boost::shared_ptr<pcl::visualization::PCLVisualizer> MultiViewPorts(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr,
                                                                    pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr,
                                                                    pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr);
std::string getpath();

//Color and gray left image
cv::Mat src,src_gray;

//Deph Map
cv::Mat depthMap;

//Camera intrinsic parameters
cv::Mat intrinsics = cv::Mat::zeros(3, 3, CV_32F);

//Corners Coordinates in Image System (u,v) pixels
vector<Point2f> corners_uv,corners1,corners2,corners3;

//Corners Coordinates in Camera System (X,Y,Z) mm
vector<Point3f> corners_world;

bool flag=true,ZED=false,KINECT=false;

int aux1=0,aux4[6],aux5,resolution;
const char *aux7;

//Color and depth filenames
char direccion1[200],direccion2[200],direccion3[200],direccion4[200],direccion5[200];
const char *buildPath;

//Chessboards Sizes
Size chess1,chess2,chess3;

//Results File
ofstream DepthError;

//Ideal, real, final point clouds
pcl::PointCloud<pcl::PointXYZRGB>::Ptr ideal_cloud(new pcl::PointCloud<pcl::PointXYZRGB>),
                                       real_cloud(new pcl::PointCloud<pcl::PointXYZRGB>),
                                       Final(new pcl::PointCloud<pcl::PointXYZRGB>);

static void CornersViewSave(){

    //Visualization of Corners and write image
    sprintf(direccion1,"%s/CornersDetected/%s/c_%dm.png",buildPath,aux7,aux5);
    imwrite(direccion1,src);
    aux5=aux5+1;

    namedWindow("Corners Detection", WINDOW_NORMAL);
    moveWindow("Corners Detection", 0, 0);
    resizeWindow("Corners Detection", 700, 700);
    imshow("Corners Detection", src);
    cout << "---Press 'q' key to continue---\n"<<endl;
    waitKey(0);
    destroyWindow("Corners Detection");
}

static void CreateIdealCLoud(int aux){

    //CREATE IDEAL AND REAL CLOUD WITH COODINATES OF THE CORNERS IN THE WORLD-OBJECT SYSTEM (mm)
    //(0,0) in the center of board
    //FIND IDEAL COORDINATES OF THE CORNERS IN THE CAMERA SYSTEM (mm)
    //Chess board Square size(mm)
    const float lado_h = 197, lado_w = 197;

    int x,y,z;

    if(aux==1){
        //number of internal corners of the plane (x-y) 3D chessboard
        //(3 en x(cols-width), 4 en y(rows-height))
        if(aux4[0]==3 && aux4[1]==4){
            //PLANE XY
            for (x = aux4[0]-1; x >= 0; x--){
               for (y = aux4[1]-1; y >= 0; y--)
                   corners_world.push_back(Point3f(290+x*lado_w, 148+y*lado_h, 0));
            }
        }
        else{
            //PLANE YZ
            for (y = aux4[0]-1; y >= 0; y--){
               for (z = aux4[1]-1; z >= 0; z--)
                   corners_world.push_back(Point3f(0,145+y*lado_h, 288+z*lado_w));
            }
        }
    }

    if(aux==2){
        if(aux4[0]==3 && aux4[2]==4){
            //PLANE YZ (1)
            //(3 en z(cols-width), 4 en y(rows-height))
            for (z=0; z < aux4[0]; z++){
               for (y = aux4[2]-1; y >= 0; y--)
                   corners_world.push_back(Point3f(0, 145+y*lado_h,288+z*lado_w));
            }
        }
        else{
            //PLANE YZ (1)
            //(3 en z(cols-width), 4 en y(rows-height))
            for (y = aux4[0]-1; y >= 0; y--){
               for (z = aux4[2]-1; z >= 0; z--)
                   corners_world.push_back(Point3f(0,145+y*lado_h, 288+z*lado_w));
            }
        }
        //PLANE XY (2)
        //(3 en x(cols-width), 4 en y(rows-height))
        for (x = aux4[1]-1; x >= 0; x--){
           for (y = aux4[3]-1; y >= 0; y--)
               corners_world.push_back(Point3f(195+x*lado_w, 148+y*lado_h, 0));
        }
    }

    if(aux==3){
        //Big Checkerboard
        ////PLANE YZ (1)
        //for (z = 0; z <aux4[3]; z++)
        //{
        //   for (y = aux4[0]-1; y >= 0; y--)
        //       corners_world.push_back(Point3f(0,45+y*lado_h, 185+z*lado_w));
        //}
        ////PLANE XZ (2)

        //for (z = 0; z <aux4[4]; z++)
        //{
        //   for (x = 0; x<aux4[1]; x++)
        //       corners_world.push_back(Point3f(185+x*lado_w, 0, 320+z*lado_h));
        //}
        ////PLANE XY (3)
        //for (y = aux4[5]-1; y >= 0; y--)
        //{
        //   for (x = 0; x < aux4[2]; x++)
        //       corners_world.push_back(Point3f(210+x*lado_w, 45+y*lado_h, 0));
        //}

        //Small Checkerboard
        //Chess board Square size(mm)
        const float lado1 = 50,lado2 = 50;

        //PLANE YZ (1)
        for (z = 2; z < aux4[0]+2; z++){
           for (y = aux4[3]+1; y > 1; y--)
               corners_world.push_back(Point3f(0,y*lado1, z*lado2));
        }
        //PLANE XZ (2)
        for (z = 2; z < aux4[1]+2; z++){
            for (x = 2; x < aux4[4]+2; x++)
               corners_world.push_back(Point3f(x*lado1, 0, z*lado2));
        }
        //PLANE XY (3)
        for (y = aux4[2]+1; y > 1; y--){
           for (x = 2; x < aux4[5]+2; x++)
               corners_world.push_back(Point3f(x*lado1, y*lado2, 0));
        }
    }
}

static void FindCorners(cv::Mat src_aux)
{
  emxArray_real_T *corners;
  emxArray_real_T *tam_chess;
  emxArray_uint8_T *image;

  //  static unsigned char* uv1 = NULL;
  //  uv1=new unsigned char[2742336];
  //  uv1=new unsigned char[2073600];

  emxInitArray_real_T(&corners, 2);
  emxInitArray_real_T(&tam_chess, 2);

  int w=src_aux.cols, h=src_aux.rows;

  static int iv2[2] = { h, w };

  int idx0,idx1;

  image = emxCreateND_uint8_T(2, *(int (*)[2])&iv2[0]);

  for (idx0 = 0; idx0 < h; idx0++) {
    for (idx1 = 0; idx1 < w; idx1++) {
      image->data[idx0 + h * idx1] = src_aux.at<uchar>(idx0,idx1);
    }
  }

  fce(image, corners, tam_chess);

  int aux1,aux2=0,aux3[tam_chess->size[0]];

  flag=true;

  for (idx0 = 0; idx0 < corners->size[0]; idx0++) {
    corners_uv.push_back(Point2f(corners->data[idx0], corners->data[idx0 + corners->size[0]]));
  }

  for (idx0 = 0; idx0 < tam_chess->size[0]; idx0++) {
    aux1=1;
    for (idx1 = 0; idx1 < tam_chess->size[1]; idx1++) {
      aux1=tam_chess->data[idx0 + tam_chess->size[0] * idx1]*aux1;
      aux4[idx0 + tam_chess->size[0] * idx1]=tam_chess->data[idx0 + tam_chess->size[0] * idx1];
    }
    aux2=aux1+aux2;
    aux3[idx0]=aux1;
  }

  if(aux2==corners->size[0]){

      if(tam_chess->size[0]==1){

        flag=false;

        for (idx0 = 0; idx0 < aux3[0]; idx0++){
          corners1.push_back(corners_uv[idx0]);
        }

        chess1=Size(tam_chess->data[1], tam_chess->data[0]);

        cornerSubPix(src_aux, corners1, chess1, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));

        corners_uv.clear();

        for (idx0 = 0; idx0 < corners1.size(); idx0++) {
            corners_uv.push_back(corners1[idx0]);
        }

        drawChessboardCorners(src, chess1, cv::Mat(corners1), true);

        cout<<"3. Find Corners: "<<tam_chess->size[0]<<" Chessboards, "<<corners->size[0]<<" Corners."<<endl;

        CreateIdealCLoud(1);
      }

      if(tam_chess->size[0]==2){

        flag=false;

        for (idx0 = 0; idx0 < aux3[0]; idx0++) {
          corners1.push_back(corners_uv[idx0]);
        }

        for (idx0 = aux3[0]; idx0 < (aux3[0]+aux3[1]); idx0++) {
          corners2.push_back(corners_uv[idx0]);
        }

        chess1=Size(tam_chess->data[2], tam_chess->data[0]);
        chess2=Size(tam_chess->data[3], tam_chess->data[1]);

        cornerSubPix(src_aux, corners1, chess1, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));
        cornerSubPix(src_aux, corners2, chess2, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));

        corners_uv.clear();

        for (idx0 = 0; idx0 < corners1.size(); idx0++) {
            corners_uv.push_back(corners1[idx0]);
            //cout<<corners_uv[idx0]<<endl;
        }

        for (idx0 = 0; idx0 < corners2.size(); idx0++) {
            corners_uv.push_back(corners2[idx0]);
            //cout<<corners_uv[idx0+corners1.size()]<<endl;
        }

        drawChessboardCorners(src, chess1, cv::Mat(corners1), true);
        drawChessboardCorners(src, chess2, cv::Mat(corners2), true);

        cout<<"3. Find Corners: "<<tam_chess->size[0]<<" Chessboards, "<<corners->size[0]<<" Corners."<<endl;

        CreateIdealCLoud(2);
      }

      if(tam_chess->size[0]==3){

        flag=false;

        for (idx0 = 0; idx0 < aux3[0]; idx0++) {
          corners1.push_back(corners_uv[idx0]);
        }

        for (idx0 = aux3[0]; idx0 < (aux3[0]+aux3[1]); idx0++) {
          corners2.push_back(corners_uv[idx0]);
        }

        for (idx0 = (aux3[0]+aux3[1]); idx0 < aux2; idx0++) {
          corners3.push_back(corners_uv[idx0]);
        }

        chess1=Size(tam_chess->data[3], tam_chess->data[0]);
        chess2=Size(tam_chess->data[4], tam_chess->data[1]);
        chess3=Size(tam_chess->data[5], tam_chess->data[2]);

        cornerSubPix(src_aux, corners1, chess1, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));
        cornerSubPix(src_aux, corners2, chess2, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));
        cornerSubPix(src_aux, corners3, chess3, Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 0.001));

        corners_uv.clear();

        for (idx0 = 0; idx0 < corners1.size(); idx0++) {
            corners_uv.push_back(corners1[idx0]);
        }

        for (idx0 = 0; idx0 < corners2.size(); idx0++) {
            corners_uv.push_back(corners2[idx0]);
        }

        for (idx0 = 0; idx0 < corners3.size(); idx0++) {
            corners_uv.push_back(corners3[idx0]);
        }

        drawChessboardCorners(src, chess1, cv::Mat(corners1), true);
        drawChessboardCorners(src, chess2, cv::Mat(corners2), true);
        drawChessboardCorners(src, chess3, cv::Mat(corners3), true);

        cout<<"3. Find Corners: "<<tam_chess->size[0]<<" Chessboards, "<<corners->size[0]<<" Corners."<<endl;

        CreateIdealCLoud(3);
     }
  }

  //  delete [] uv1;  // When done, free memory pointed to by uv1.
  //  uv1 = NULL;     // Clear uv1 to prevent using invalid memory reference.

  emxDestroyArray_real_T(tam_chess);
  emxDestroyArray_real_T(corners);
  emxDestroyArray_uint8_T(image);
}

static void RealCloud_Reg_RMS(){

    // FIND REAL COORDINATES OF THE CORNERS IN THE CAMERA SYSTEM (mm)
    // CONVERT PIXEL COORDINATES TO CAMERA COORDINATES

    switch (resolution) {
    case 2742336:
            intrinsics.at<float>(0,0)=1400;    //fx
            intrinsics.at<float>(1,1)=1400;    //fy
            intrinsics.at<float>(0,2)=1208;    //Cx
            intrinsics.at<float>(1,2)=565.867; //Cy
            intrinsics.at<float>(2,2)=1;
        break;
    case 2073600:
        if (ZED){
            intrinsics.at<float>(0,0)=1400;    //fx
            intrinsics.at<float>(1,1)=1400;    //fy
            intrinsics.at<float>(0,2)=1064;    //Cx
            intrinsics.at<float>(1,2)=484.867; //Cy
            intrinsics.at<float>(2,2)=1;
        }
        else{
            intrinsics.at<float>(0,0)=1056.395; //fx
            intrinsics.at<float>(1,1)=1057.192; //fy
            intrinsics.at<float>(0,2)=963.771;  //Cx
            intrinsics.at<float>(1,2)=535.469;  //Cy
            intrinsics.at<float>(2,2)=1;
        }
        break;
    case 921600:
            intrinsics.at<float>(0,0)=699.999; //fx
            intrinsics.at<float>(1,1)=699.999; //fy
            intrinsics.at<float>(0,2)=690.502; //Cx
            intrinsics.at<float>(1,2)=330.933; //Cy
            intrinsics.at<float>(2,2)=1;
        break;
    case 252672:
            intrinsics.at<float>(0,0)=350;      //fx
            intrinsics.at<float>(1,1)=350;      //fy
            intrinsics.at<float>(0,2)=360.751;  //Cx
            intrinsics.at<float>(1,2)=172.967;  //Cy
            intrinsics.at<float>(2,2)=1;
        break;
    case 307200:
            intrinsics.at<float>(0,0)=522.9552; //fx
            intrinsics.at<float>(1,1)=523.8785; //fy
            intrinsics.at<float>(0,2)=329.0866; //Cx
            intrinsics.at<float>(1,2)=254.7133; //Cy
            intrinsics.at<float>(2,2)=1;
        break;
    case 217088:
            intrinsics.at<float>(0,0)=367.736373303424; //fx
            intrinsics.at<float>(1,1)=367.498516141750; //fy
            intrinsics.at<float>(0,2)=264.409624145552; //Cx
            intrinsics.at<float>(1,2)=201.559776273164; //Cy
            intrinsics.at<float>(2,2)=1;
        break;
    default:
        break;
    }

    vector<Point3f> corners_camera; //Point coordinates on camera coodinates system
    Point3f p3d;

    float coord_Z;
    int i;

    for (i = 0; i < corners_uv.size(); i++)
    {
      coord_Z = depthMap.at<unsigned short>(cvRound(corners_uv[i].y), cvRound(corners_uv[i].x));
      //cout<<coord_Z<<endl;
      p3d=Point3f((corners_uv[i].x-intrinsics.at<float>(0,2))*coord_Z/intrinsics.at<float>(0,0), (corners_uv[i].y-intrinsics.at<float>(1,2))*coord_Z/intrinsics.at<float>(1,1), coord_Z);
      corners_camera.push_back(p3d);
    }

    // Create ideal point cloud

    uint8_t r(0), g(0), b(0);

    for (i = 0; i < corners_world.size(); i++)
    {
        pcl::PointXYZRGB point;
        point.x = corners_world[i].x;
        point.y = corners_world[i].y;
        point.z = corners_world[i].z;
        uint32_t rgb = (static_cast<uint32_t>(r) << 16 | static_cast<uint32_t>(g) << 8 | static_cast<uint32_t>(b));
        point.rgb = *reinterpret_cast<float*>(&rgb);
        ideal_cloud->points.push_back (point);
    }

    ideal_cloud->width = (int) ideal_cloud->points.size ();
    ideal_cloud->height = 1;
    //pcl::io::savePCDFileASCII("ideal_cloud.pcd", *ideal_cloud);

    //Create real point cloud

    uint8_t rr(255), gg(0), bb(0);

    for (int i = 0; i < corners_camera.size(); i++)
    {
        pcl::PointXYZRGB point;
        point.x = corners_camera[i].x;
        point.y = corners_camera[i].y;
        point.z = corners_camera[i].z;
        uint32_t rgb = (static_cast<uint32_t>(rr) << 16 | static_cast<uint32_t>(gg) << 8 | static_cast<uint32_t>(bb));
        point.rgb = *reinterpret_cast<float*>(&rgb);
        real_cloud->points.push_back (point);
    }

    real_cloud->width = (int) real_cloud->points.size ();
    real_cloud->height = 1;
    //pcl::io::savePCDFileASCII("real_cloud.pcd", *real_cloud);

    // REGISTRATION
    //    pcl::IterativeClosestPoint<pcl::PointXYZRGB, pcl::PointXYZRGB> icp;
    //    icp.setInputCloud (real_cloud);
    //    icp.setInputTarget (ideal_cloud);
    //    icp.setMaximumIterations (1);
    //    icp.align (*Final);

    pcl::registration::TransformationEstimationSVD<pcl::PointXYZRGB,pcl::PointXYZRGB> ASVD;
    pcl::registration::TransformationEstimationSVD<pcl::PointXYZRGB,pcl::PointXYZRGB>::Matrix4 transform_matrix;
    ASVD.estimateRigidTransformation (*real_cloud,*ideal_cloud,transform_matrix);
    pcl::transformPointCloud(*real_cloud, *Final, transform_matrix);

    //pcl::io::savePCDFileASCII("Final.pcd", *Final);

    cout<<"4. Create Ideal and Actual Point Clouds"<<endl;

    // DEPTH VS ERROR

    cv::Mat corners_cameraZ(1, corners_camera.size(), CV_32FC1);

    for (i = 0; i < corners_camera.size(); i++) {
        corners_cameraZ.at<float>(0, i) = corners_camera[i].z;
    }

    double minVal;
    double maxVal;

    minMaxLoc( corners_cameraZ, &minVal, &maxVal);

    double dist, dist_x, dist_y, dist_z, valx, valy, valz, error;

    for (size_t i = 0; i < Final->points.size (); i++)
    {
      valx=Final->points[i].x-ideal_cloud->points[i].x;
      valy=Final->points[i].y-ideal_cloud->points[i].y;
      valz=Final->points[i].z-ideal_cloud->points[i].z;

      dist_x = pow(valx,2);
      dist_y = pow(valy,2);
      dist_z = pow(valz,2);

      dist = dist_x+dist_y+dist_z;
      error = pow(dist,0.5);

      if(corners_camera[i].z > 500 && corners_camera[i].z < 20001)
          DepthError <<corners_camera[i].z<< " " <<error<<endl;
    }

    cout <<"5. Obtem Depth Error Succesfull"<<"\n"<<endl;
}

// Multiple Point Clouds Visualizer
boost::shared_ptr<pcl::visualization::PCLVisualizer> MultiViewPorts(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud_ic,
                                                                    pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud_cc,
                                                                    pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud_res)
{
    // Open 3D viewer
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("A_Versatile_Method_for_Depth_Data_Error_Estimation_in_RGB-D_Sensors"));

    // Create three separated viewports
    int v1 (0), v2 (1), v3 (2);

    // View a Ideal Cube Point Cloud
    viewer->createViewPort (0.0, 0.5, 0.5, 1, v1);
    viewer->setBackgroundColor(1.0, 1.0, 1.0, v1);
    viewer->addText("Ideal Corners Location", 10, 10, 20, 0.0, 0.0, 0.0, "v1 text", v1);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb1(cloud_ic);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_ic, rgb1,"cloud1", v1);

    // View a Captured Cube Point Cloud
    viewer->createViewPort (0.5, 0.5, 1.0, 1.0, v2);
    viewer->setBackgroundColor(1.0, 1.0, 1.0, v2);
    viewer->addText("Actual Corners Location", 10, 10, 20, 0.0, 0.0, 0.0, "v2 text", v2);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb2(cloud_cc);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_cc, rgb2,"cloud2", v2);

    //View a result of the SVD Algorithm
    viewer->createViewPort (0.0, 0.0, 1.0, 0.5, v3);
    viewer->setBackgroundColor(1.0, 1.0, 1.0, v3);
    viewer->addText("SVD Registration-Ideal and Actual Corners Location", 10, 10, 20, 0.0, 0.0, 0.0, "v3 text", v3);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb3(cloud_res);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_res, rgb3,"cloud3",   v3);
    viewer->addPointCloud<pcl::PointXYZRGB> (cloud_ic, rgb1,"cloud1v3",  v3);

    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5,"cloud1",v1);
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5,"cloud2",v2);
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5,"cloud3",v3);
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5,"cloud1v3",v3);

    viewer->addCoordinateSystem(900);
    viewer->setSize(700,700);
    viewer->setPosition(800,0);

    //viewer->initCameraParameters();

    viewer->setCameraPosition(900, 900, 900, 0, 0, 0, v1);
    viewer->setCameraPosition(900, 900, 900, 0, 0, 0, v2);
    viewer->setCameraPosition(900, 900, 900, 0, 0, 0, v3);

    pcl::PointXYZRGBA point;

    point.getArray3fMap () << 1000, 0, 0;
    viewer->addText3D ("Xw", point, 60, 1, 0, 0, "x_");

    point.getArray3fMap () << 0, 1000, 0;
    viewer->addText3D ("Yw", point, 60, 0, 1, 0, "y_");

    point.getArray3fMap () << 0, 0, 1100;
    viewer->addText3D ("Zw", point, 60, 0, 0, 1, "z_");

    return(viewer);
}

std::string getpath() {
    char buf[PATH_MAX + 1];
    if(readlink("/proc/self/exe", buf, sizeof(buf) - 1) == -1)
        throw std::string("readlink() failed");
    std::string str(buf);
    return str.substr(0, str.rfind('/'));
}

int main(int argc, char* argv[])
{
    cout<<"\nA_Versatile_Method_for_Depth_Data_Error_Estimation_in_RGB-D_Sensors"<<endl;

    if(argc != 4)
    {
        cout << "Use: './depth_error <save_windows> <device> <devide_resolution>'\n"
             << "      <save_windows>  : '0' off\n"
             << "                      : '1' on\n"
             << "      <device>        : '1' for ZED\n"
             << "                      : '2' for KINECT\n"
            //  << "      <devide_resolution>: '1' 2208x1242\n"
            //  << "                         : '2' 1920x1080\n"
            //  << "                         : '3' 1280x720\n"
             << "      <devide_resolution>: '4' 672x376 (for ZED)\n"
             << "                         : '5' 640x480 (for KINECT 1)\n"
             << "                         : '6' 512x424 (for KINECT 2)\n"
             << endl;
        return 1;
    }

    int param1, param2, param3;

    sscanf (argv[1],"%d",&param1);
    sscanf (argv[2],"%d",&param2);
    sscanf (argv[3],"%d",&param3);

    switch(param2){
    case 1:
        ZED=true;
        break;
    case 2:
        KINECT=true;
        break;
    default:
        return 1;
        break;
    }

    switch(param3){
    // case 1:
    //     aux7= "2208x1242";
    //     break;
    // case 2:
    //     if(ZED)
    //         aux7= "1920x1080";
    //     if(KINECT)
    //         aux7= "1920x1080_kinect";
    //     break;
    // case 3:
    //     aux7= "1280x720";
    //     break;
    case 4:
        aux7= "672x376";
        break;
    case 5:
        aux7= "640x480";
        break;
    case 6:
        aux7= "512x424";
        break;
    default:
        return 1;
        break;
    }

    string aux8 = getpath();

    buildPath = aux8.c_str();

    sprintf(direccion5,"%s/Results/%s/Table_%s.txt",buildPath,aux7,aux7);

    DepthError.open(direccion5);

    fce_initialize();

    //    double tiempo=0;
    //    double tiempo2=0;

    for(int i = 1; i < 6; i++)
    {
       cout<<"\n1. Init"<<endl;

       //Load Undistorted images because the equations for find depth
       //are not include distortions coeficients
       //In stereo devices image is left image

       cout<<"2. Load Image and Depth Map"<<endl;

       sprintf(direccion1,"../Data/2208x1242/%d_l.png",i);

       switch(param3){
        //    case 1:
        //        sprintf(direccion1,"../Data/2208x1242/%d_l.png",i);
        //        sprintf(direccion2,"../Data/2208x1242/%d_d.png",i);
        //        sprintf(direccion3,"%s/Results/2208x1242/%d_r.png",buildPath,i);
        //    break;
        //    case 2:
        //        if(ZED){
        //         sprintf(direccion1,"../Data/1920x1080/%d_l.png",i);
        //         sprintf(direccion2,"../Data/1920x1080/%d_d.png",i);
        //         sprintf(direccion3,"%s/Results/1920x1080/%d_r.png",buildPath,i);
        //        }
        //        if(KINECT){
        //         sprintf(direccion1,"../Data/1920x1080_kinect/%d_l.png",i);
        //         sprintf(direccion2,"../Data/1920x1080_kinect/%d_d.png",i);
        //         sprintf(direccion3,"%s/Results/1920x1080_kinect/%d_r.png",buildPath,i);
        //        }
        //        break;
        //    case 3:
        //        sprintf(direccion1,"../Data/1280x720/%d_l.png",i);
        //        sprintf(direccion2,"../Data/1280x720/%d_d.png",i);
        //        sprintf(direccion3,"%s/Results/1280x720/%d_r.png",buildPath,i);
        //        break;
           case 4:
               sprintf(direccion1,"../Data/672x376/%d_l.png",i);
               sprintf(direccion2,"../Data/672x376/%d_d.png",i);
               sprintf(direccion3,"%s/Results/672x376/%d_r.png",buildPath,i);
               break;
           case 5:
               sprintf(direccion1,"../Data/640x480/%d_l.png",i);
               sprintf(direccion2,"../Data/640x480/%d_d.png",i);
               sprintf(direccion3,"%s/Results/640x480/%d_r.png",buildPath,i);
               break;
           case 6:
               sprintf(direccion1,"../Data/512x424/%d_l.png",i);
               sprintf(direccion2,"../Data/512x424/%d_d.png",i);
               sprintf(direccion3,"%s/Results/512x424/%d_r.png",buildPath,i);
               break;
           default:
               return 1;
               break;
       }

       sprintf(direccion4,"%d_l.png",i);

       //_COLOR If set, always convert image to the 3 channel BGR color image.
       //_ANYDEPTH If set, return 16-bit/32-bit image when the input has the corresponding depth, otherwise convert it to 8-bit.

       //pcl::console::TicToc t_captura;
       //t_captura.tic();

       src = imread(direccion1, CV_LOAD_IMAGE_COLOR);

       //tiempo+=t_captura.toc ();
       //t_captura.tic();

       depthMap = imread(direccion2, CV_LOAD_IMAGE_ANYDEPTH);

       //tiempo2+=t_captura.toc ();

       resolution = src.rows*src.cols;

       cvtColor(src, src_gray, CV_BGR2GRAY);

       FindCorners(src_gray);

       if(flag)
       {
           cout<<"---No Corners Detected on Image '"<<direccion4<<"'"<<endl;
       }
       else
       {
            RealCloud_Reg_RMS();

            if(param1==1)
            {
                aux5=i;
                CornersViewSave();

                boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;

                viewer = MultiViewPorts(ideal_cloud, real_cloud, Final);

                while(!viewer->wasStopped ())
                {
                  viewer->spinOnce ();
                }

                viewer->saveScreenshot(direccion3);
                viewer->close();
              }
       }

       corners_uv.clear();
       corners1.clear();
       corners2.clear();
       corners3.clear();
       corners_world.clear();
       ideal_cloud->points.clear();
       Final->points.clear();
       real_cloud->points.clear();
    }
    //    cout<<"Promedio1: "<< tiempo/20<<"ms"<<endl;
    //    cout<<"Promedio2: "<< tiempo2/20<<"ms"<<endl;
    DepthError.close();
    fce_terminate();
    return 0;
}