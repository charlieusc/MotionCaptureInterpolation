#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "motion.h"
#include "interpolator.h"
#include "types.h"
using namespace std;

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);
    
      if(startKeyframe>=200 && startKeyframe<=500){
          cout<<startPosture->bone_rotation[0].p[2]<<endl;
      }
    
    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
        for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++){
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;
         
        if (bone == 0 && startKeyframe+frame>=200 && startKeyframe+frame<=500) {
            cout<<interpolatedPosture.bone_rotation[bone].p[2]<<endl;
            //cout<<pInputMotion->GetPosture(startKeyframe+frame)->bone_rotation[bone].p[2]<<endl;
        }
        
        
        
    }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}
double degree2radian(double degree){
    return degree / 180 * M_PI;
}
void matrixMultip(double R1[9], double R2[9], double RR[9]){
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
                RR[i*3+j] = 0;
        }
    }
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                RR[i*3+j] += R1[i*3+k]*R2[k*3+j];
            }
            
        }
    }
}
void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
  // students should implement this
    double theta1 = degree2radian(angles[0]);
    double theta2 = degree2radian(angles[1]);
    double theta3 = degree2radian(angles[2]);
    
    double Rz[] = {cos(theta3), -sin(theta3), 0,
                   sin(theta3), cos(theta3),  0,
                   0,           0,            1};
    
    double Ry[] = {cos(theta2), 0,  sin(theta2),
                   0,           1,            0,
                   -sin(theta2), 0, cos(theta2)};
    
    double Rx[] = {1,           0,            0,
                   0,  cos(theta1), -sin(theta1),
                   0,  sin(theta1), cos(theta1)};
    double Rr[9];
    
    matrixMultip(Rz, Ry, Rr);
    matrixMultip(Rr, Rx, R);
}

void calculateBEA1B2(vector &p0, vector &p1, vector &p2, vector &p3, vector &a1, vector &b2, bool fp0, bool fp3){
    double scaler = 1.0/3;
    if (fp0 == false) {
        //first key frame
        a1 = p1 + (p2*2 - p3 - p1) * scaler;
        b2 = p2 - ((p2*2 - p1 + p3) * 0.5 - p2) * scaler;
    }else if(fp3 == false){
        //last key frame
        a1 = p1 + ((p1*2 - p0 + p2) * 0.5 - p1) * scaler;
        b2 = p2 + (p1*2 - p0 - p2) * scaler;
    }else{
        //inter frame
        a1 = p1 + ((p1*2 - p0 + p2) * 0.5 - p1) * scaler;
        b2 = p2 - ((p2*2 - p1 + p3) * 0.5 - p2) * scaler;
    }
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  // students should implement this
    int inputLength = pInputMotion->GetNumFrames();// frames are indexed 0, ..., inputLength-1
    int startKeyframe = 0;
    vector p0,p1,p2,p3,a1,b2;
    
    vector bp0,bp1,bp2,bp3,ba1,bb2;
    bool fp0, fp3, fbp0, fbp3;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;
        int preKeyframe = startKeyframe - N - 1;
        int nextKeyframe = endKeyframe + N + 1;
        
        Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
        Posture * prePosture = NULL;
        Posture * nextPosture = NULL;
        
        
        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);
        
        if(preKeyframe >= 0){
            prePosture = pInputMotion->GetPosture(preKeyframe);
        }
        if(nextKeyframe<inputLength){
            nextPosture = pInputMotion->GetPosture(nextKeyframe);
        }
        
        //calculate p1, a1, b2, p2 for root interpolation
        p1 = startPosture->root_pos;
        p2 = endPosture->root_pos;

        if (startKeyframe == 0) {
            //first key frame
            fp0 = false;
            p3 = nextPosture->root_pos;
        }else if(nextKeyframe>inputLength){
            //last key frame
            p0 = prePosture->root_pos;
            fp3 = false;
        }else{
            //inter frame
            p0 = prePosture->root_pos;
            p3 = nextPosture->root_pos;
            fp0 = true;
            fp3 = true;
        }
        calculateBEA1B2(p0, p1, p2, p3, a1, b2, fp0, fp3);
        
        
        if(startKeyframe>=200 && startKeyframe<=500){
            cout<<startPosture->bone_rotation[0].p[2]<<endl;
        }
        
        
        // interpolate in between
        for(int frame=1; frame<=N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N+1);

            // interpolate root position
            interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, a1, b2, p2);
            
            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++){
                //calculate bp1, ba1, bb2, bp2 for bone interpolation
                bp1 = startPosture->bone_rotation[bone];
                bp2 = endPosture->bone_rotation[bone];
                
                if (startKeyframe == 0) {
                    //first key frame
                    fbp0 = false;
                    bp3 = nextPosture->bone_rotation[bone];
                }else if(nextKeyframe>inputLength){
                    //last key frame
                    bp0 = prePosture->bone_rotation[bone];
                    fp3 = false;
                }else{
                    //inter frame
                    bp0 = prePosture->bone_rotation[bone];
                    bp3 = nextPosture->bone_rotation[bone];
                    fp0 = true;
                    fbp3 = true;
                }

                calculateBEA1B2(bp0, bp1, bp2, bp3, ba1, bb2, fbp0, fbp3);
                interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, bp1, ba1, bb2, bp2);
                
                if (bone == 0 && startKeyframe+frame>=200 && startKeyframe+frame<=500) {
                    cout<<interpolatedPosture.bone_rotation[bone].p[2]<<endl;
                    //cout<<pInputMotion->GetPosture(startKeyframe+frame)->bone_rotation[bone].p[0]<<endl;
                }
            }
            
            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }
        
        startKeyframe = endKeyframe;
    }
    
    for(int frame=startKeyframe+1; frame<inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  // students should implement this
    int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1
    
    int startKeyframe = 0;
    //initiate quaternion for bone interpulation
    Quaternion<double>qStart, qEnd, qInter;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;
        
        Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
        
        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);
        
        if(startKeyframe>=200 && startKeyframe<=500){
            cout<<startPosture->bone_rotation[0].p[2]<<endl;
        }
        
        // interpolate in between
        for(int frame=1; frame<=N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N+1);
            
            // interpolate root position
            interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;
            
            
            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++){
                //give value to start and end quaternion
                Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
                Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);
                qInter = Slerp(t, qStart, qEnd);
                Quaternion2Euler(qInter, interpolatedPosture.bone_rotation[bone].p);
                if (bone == 0 && startKeyframe+frame>=200 && startKeyframe+frame<=500) {
                    cout<<interpolatedPosture.bone_rotation[bone].p[2]<<endl;
                    //cout<<pInputMotion->GetPosture(startKeyframe+frame)->bone_rotation[bone].p[0]<<endl;
                }

            }
            
            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }
        
        startKeyframe = endKeyframe;
    }
    
    for(int frame=startKeyframe+1; frame<inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::calculateBQAnBn(Quaternion<double> &q0, Quaternion<double> &q1, Quaternion<double> &q2, Quaternion<double> &q3, Quaternion<double> &an, Quaternion<double> &bn, bool fq0, bool fq3){
    double scaler = 1.0/3;
    Quaternion<double> a1;
    Quaternion<double> anp;
    Quaternion<double> bN;
    Quaternion<double> temp;

    
    if (fq0 == false) {
        //first key frame
        temp = Double(q3, q2);
        a1 = Slerp(scaler, q1, temp);
        an = a1;
        
        temp = Double(q1, q2);
        anp = Slerp(0.5, temp, q3);
        bn = Slerp(-scaler, q2, anp);
    }else if(fq3 == false){
        //last key frame
        temp = Double(q0, q1);
        anp = Slerp(0.5, temp, q2);
        an = Slerp(scaler, q1, anp);
        
        temp = Double(q0, q1);
        bN = Slerp(scaler, q2, temp);
        bn = bN;
        
    }else{
        //inter frame
        temp = Double(q0, q1);
        anp = Slerp(0.5, temp, q2);
        an = Slerp(scaler, q1, anp);
        
        temp = Double(q1, q2);
        anp = Slerp(0.5, temp, q3);
        bn = Slerp(-scaler, q2, anp);
    }
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
    int inputLength = pInputMotion->GetNumFrames();// frames are indexed 0, ..., inputLength-1
    int startKeyframe = 0;
    vector p0,p1,p2,p3,a1,b2;
    
    vector bp0,bp1,bp2,bp3,ba1,bb2;
    Quaternion<double> q0,q1,q2,q3,an,bn,qInter;
    bool fp0, fp3, fq0, fq3;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;
        int preKeyframe = startKeyframe - N - 1;
        int nextKeyframe = endKeyframe + N + 1;
        
        Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
        Posture * prePosture = NULL;
        Posture * nextPosture = NULL;
        
        
        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);
        
        if(preKeyframe >= 0){
            prePosture = pInputMotion->GetPosture(preKeyframe);
        }
        if(nextKeyframe<inputLength){
            nextPosture = pInputMotion->GetPosture(nextKeyframe);
        }
        
        //calculate p1, a1, b2, p2 for root interpolation
        p1 = startPosture->root_pos;
        p2 = endPosture->root_pos;
        
        if (startKeyframe == 0) {
            //first key frame
            fp0 = false;
            p3 = nextPosture->root_pos;
        }else if(nextKeyframe>inputLength){
            //last key frame
            p0 = prePosture->root_pos;
            fp3 = false;
        }else{
            //inter frame
            p0 = prePosture->root_pos;
            p3 = nextPosture->root_pos;
            fp0 = true;
            fp3 = true;
        }
        calculateBEA1B2(p0, p1, p2, p3, a1, b2, fp0, fp3);
        
        if(startKeyframe>=600 && startKeyframe<=800){
            cout<<startPosture->bone_rotation[2].p[0]<<endl;
        }
        
        // interpolate in between
        for(int frame=1; frame<=N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N+1);
            
            // interpolate root position
            interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, a1, b2, p2);
            
            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++){
                //calculate q1, an, bn, q2 for bone interpolation
                Euler2Quaternion(startPosture->bone_rotation[bone].p, q1);
                Euler2Quaternion(endPosture->bone_rotation[bone].p, q2);
                
                if (startKeyframe == 0) {
                    //first key frame
                    fq0 = false;
                    Euler2Quaternion(nextPosture->bone_rotation[bone].p, q3);
                }else if(nextKeyframe>inputLength){
                    //last key frame
                    Euler2Quaternion(prePosture->bone_rotation[bone].p, q0);
                    fq3 = false;
                }else{
                    //inter frame
                    Euler2Quaternion(prePosture->bone_rotation[bone].p, q0);
                    Euler2Quaternion(nextPosture->bone_rotation[bone].p, q3);
                    fq0 = true;
                    fq3 = true;
                }
                
                calculateBQAnBn(q0, q1, q2, q3, an, bn, fq0, fq3);
                qInter = DeCasteljauQuaternion(t, q1, an, bn, q2);
                Quaternion2Euler(qInter, interpolatedPosture.bone_rotation[bone].p);
                
                if (bone == 2 && startKeyframe+frame>=600 && startKeyframe+frame<=800) {
                    cout<<interpolatedPosture.bone_rotation[bone].p[0]<<endl;
                    //cout<<pInputMotion->GetPosture(startKeyframe+frame)->bone_rotation[bone].p[0]<<endl;
                }
                
            }
            
            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }
        
        startKeyframe = endKeyframe;
    }
    
    for(int frame=startKeyframe+1; frame<inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
    
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
    double R[9];
    Euler2Rotation(angles, R);
    q = Quaternion<double>::Matrix2Quaternion(R);
    q.Normalize();
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
    double R[9];
    q.Quaternion2Matrix(R);
    Rotation2Euler(R, angles);
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_)
{
  Quaternion<double> result;
    double cosTheta = qStart.Gets()*qEnd_.Gets() + qStart.Getx()*qEnd_.Getx() + qStart.Gety()*qEnd_.Gety() + qStart.Getz()*qEnd_.Getz();
    double theta = acos(cosTheta);

    //always chose the shortest angle
    if (cosTheta<0) {
        theta = M_PI - acos(-cosTheta);
        qEnd_ = qEnd_ * (-1.0);
    }
    
    //if theta equals to 0 reutrn qStart
    if (sin(theta) == 0.0) {
        result = qStart;
        return result;
    }
    
    result = sin((1-t)*theta)/sin(theta) * qStart + sin(t*theta)/sin(theta) * qEnd_;
    result.Normalize();
  return result;
    
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  Quaternion<double> result;
    /*
    double pdotq = p.Gets()*q.Gets() + p.Getx()*q.Getx()+ p.Gety()*q.Gety() + p.Getz()*q.Getz();
    result = 2*pdotq*q - p;
     */
    result = Slerp(2, p, q);
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
    vector result;
    vector Q0, Q1, Q2, R0, R1;
    Q0 = p0 + (p1 - p0) * t;
    Q1 = p1 + (p2 - p1) * t;
    Q2 = p2 + (p3 - p2) * t;
    R0 = Q0 + (Q1 - Q0) * t;
    R1 = Q1 + (Q2 - Q1) * t;
    result = R0 + (R1 - R0) * t;
    return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
  Quaternion<double> result;
    Quaternion<double> Q0, Q1, Q2, R0, R1;
    Q0 = Slerp(t, p0, p1);
    Q1 = Slerp(t, p1, p2);
    Q2 = Slerp(t, p2, p3);
    R0 = Slerp(t, Q0, Q1);
    R1 = Slerp(t, Q1, Q2);
    result = Slerp(t,R0, R1);
  return result;
}

