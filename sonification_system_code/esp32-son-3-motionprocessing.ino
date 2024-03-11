void compute_motion (void) {

  // check for out-out-range gyro readings ( > 490 dps ): 
  #ifdef SENSOR3
    float high_gyro_reading = 
    max(max(max(max(Sensor1_Gyro.x,Sensor1_Gyro.y),Sensor1_Gyro.z),
            max(max(Sensor2_Gyro.x,Sensor2_Gyro.y),Sensor2_Gyro.z)),
            max(max(Sensor3_Gyro.x,Sensor3_Gyro.y),Sensor3_Gyro.z));
  #else 
    float high_gyro_reading = 
    max(max(max(Sensor1_Gyro.x,Sensor1_Gyro.y),Sensor1_Gyro.z),
        max(max(Sensor2_Gyro.x,Sensor2_Gyro.y),Sensor2_Gyro.z));
  #endif

  if ( high_gyro_reading > high_gyro_readings ) {
    high_gyro_readings = high_gyro_reading;
  }

  // find the ROT vectors / total rotation about each axis from start (in radians)
  Sensor1_Rot = vadd( Sensor1_Rot, vscaler_mult(Sensor1_Gyro,dt) );
  Sensor2_Rot = vadd( Sensor2_Rot, vscaler_mult(Sensor2_Gyro,dt) );
  #ifdef SENSOR3
    Sensor3_Rot = vadd( Sensor3_Rot, vscaler_mult(Sensor3_Gyro,dt) );
  #endif

  #ifdef SUBTRACTGYROBIAS
    Sensor1_Gyro = vsubtract(Sensor1_Gyro,Sensor1_GyroBias); // recall that the GyroBias vectors have already had each component multipled by dt
    Sensor2_Gyro = vsubtract(Sensor2_Gyro,Sensor2_GyroBias);
    #ifdef SENSOR3
      Sensor3_Gyro = vsubtract(Sensor3_Gyro,Sensor3_GyroBias);
    #endif
  #endif 
  
  // Find quats describing sensor rotation during sample period
  Vectors q1 = rot_quat ( vscaler_mult(Sensor1_Gyro,dt), Sensor1_Xaxis, Sensor1_Yaxis, Sensor1_Zaxis );
  Vectors q2 = rot_quat ( vscaler_mult(Sensor2_Gyro,dt), Sensor2_Xaxis, Sensor2_Yaxis, Sensor2_Zaxis );
  #ifdef SENSOR3
    Vectors q3 = rot_quat ( vscaler_mult(Sensor3_Gyro,dt), Sensor3_Xaxis, Sensor3_Yaxis, Sensor3_Zaxis );
  #endif

  // update axis position relative relative to their initial position / in coordinates defined by their initial position
  Sensor1_Xaxis = qvq(q1,Sensor1_Xaxis);
  Sensor1_Yaxis = qvq(q1,Sensor1_Yaxis);
  Sensor1_Zaxis = qvq(q1,Sensor1_Zaxis);
  Sensor2_Xaxis = qvq(q2,Sensor2_Xaxis);
  Sensor2_Yaxis = qvq(q2,Sensor2_Yaxis);
  Sensor2_Zaxis = qvq(q2,Sensor2_Zaxis);
  #ifdef SENSOR3
    Sensor3_Xaxis = qvq(q3,Sensor3_Xaxis);
    Sensor3_Yaxis = qvq(q3,Sensor3_Yaxis);
    Sensor3_Zaxis = qvq(q3,Sensor3_Zaxis);
  #endif

  // Compute quats describing total sensor rotation from initial position 
  if (FIRSTSAMPLE) {
    FIRSTSAMPLE = false;
    Sensor1_Quat = q1;
    Sensor2_Quat = q2;
    #ifdef SENSOR3
      Sensor3_Quat = q3;
    #endif
  } else {
    Sensor1_Quat = quat_mult( q1, Sensor1_Quat );
    Sensor2_Quat = quat_mult( q2, Sensor2_Quat );
    #ifdef SENSOR3
      Sensor3_Quat = quat_mult( q3, Sensor3_Quat );
    #endif
  }
  
}

Vectors rot_quat ( Vectors r, Vectors ax, Vectors ay, Vectors az ) {
  ax = normalize(ax); // we must ensure this is normalized, although it should be as well (if increased speed needed, drop this step).
  ay = normalize(ay); // we must ensure this is normalized, although it should be as well (if increased speed needed, drop this step).
  az = normalize(az); // we must ensure this is normalized, although it should be as well (if increased speed needed, drop this step).
  Vectors qx = formquat(r.x,ax);
  Vectors qy = formquat(r.y,ay);
  Vectors qz = formquat(r.z,az);
  Vectors q = quat_mult( qx, quat_mult( qy, qz) );
  q = normalize4(q); // we must ensure this is normalized, although it should be as well (if increased speed needed, drop this step).
  return(q);
}

Vectors formquat ( float theta, Vectors v ) {
  Vectors q;
  float cos2theta = cosf(theta/2.0);
  float sin2theta = sinf(theta/2.0);
  q.r = cos2theta;
  q.x = sin2theta * v.x;
  q.y = sin2theta * v.y;
  q.z = sin2theta * v.z;
  return (q); 
}

Vectors quat_mult ( Vectors q1, Vectors q2 ) {
  Vectors q = vadd( vscaler_mult(q2,q1.r) , vadd( vscaler_mult(q1,q2.r) , crossp(q1,q2) ) );
  q.r = q1.r * q2.r - dotp(q1,q2); 
  return(q);
}

Vectors normalize ( Vectors v ) {
  Vectors w;
  w = vdivide(v,vmag(v));
  return(w); 
}

Vectors normalize4 ( Vectors v ) {
  Vectors w;
  w = vdivide4(v,vmag4(v));
  return(w); 
}

Vectors vscaler_mult ( Vectors thisVector1, float thisfloat ) {
  Vectors output;
  output.x = thisVector1.x * thisfloat;
  output.y = thisVector1.y * thisfloat;
  output.z = thisVector1.z * thisfloat;
  return(output);
}

float dotp ( Vectors thisVector1, Vectors thisVector2 ) {
  float dotproduct;
  dotproduct = thisVector1.x * thisVector2.x + thisVector1.y * thisVector2.y + thisVector1.z * thisVector2.z ;
  return dotproduct;
}

Vectors crossp ( Vectors thisVector1, Vectors thisVector2 ) {
  Vectors crossproduct;
  crossproduct.x = thisVector1.y * thisVector2.z - thisVector1.z * thisVector2.y ;
  crossproduct.y = thisVector1.z * thisVector2.x - thisVector1.x * thisVector2.z ;
  crossproduct.z = thisVector1.x * thisVector2.y - thisVector1.y * thisVector2.x ;
  return(crossproduct);
}

float vdst ( Vectors thisVector1, Vectors thisVector2 ) {
  float x = thisVector1.x - thisVector2.x;
  float y = thisVector1.y - thisVector2.y;
  float z = thisVector1.z - thisVector2.z;
  float distance = sqrtf ( x*x + y*y + z*z );
  return distance;
}

float vdst4 ( Vectors thisVector1, Vectors thisVector2 ) {
  float x = thisVector1.x - thisVector2.x;
  float y = thisVector1.y - thisVector2.y;
  float z = thisVector1.z - thisVector2.z;
  float r = thisVector1.r - thisVector2.r;
  float distance = sqrtf ( x*x + y*y + z*z + r*r );
  return distance;
}

float vdstSq4 ( Vectors thisVector1, Vectors thisVector2 ) {
  float x = thisVector1.x - thisVector2.x;
  float y = thisVector1.y - thisVector2.y;
  float z = thisVector1.z - thisVector2.z;
  float r = thisVector1.r - thisVector2.r;
  float distance = x*x + y*y + z*z + r*r ;
  return distance;
}

float vmag ( Vectors thisVector ) {
  float MagV;
  MagV = sqrtf ( thisVector.x * thisVector.x + thisVector.y * thisVector.y + thisVector.z * thisVector.z );
  return MagV;
}

float vmag4 ( Vectors thisVector ) {
  float MagV;
  MagV = sqrtf ( thisVector.r * thisVector.r + thisVector.x * thisVector.x + thisVector.y * thisVector.y + thisVector.z * thisVector.z );
  return MagV;
}

float vmagSqrd ( Vectors thisVector ) {
  float MagV;
  MagV = thisVector.x * thisVector.x + thisVector.y * thisVector.y + thisVector.z * thisVector.z ;
  return MagV;
}

Vectors vdivide ( Vectors v, float r ) {
  Vectors w;
  w.x = v.x /r;
  w.y = v.y /r;
  w.z = v.z /r;
  return(w);
}

Vectors vdivide4 ( Vectors v, float r ) {
  Vectors w;
  w.r = v.r /r;
  w.x = v.x /r;
  w.y = v.y /r;
  w.z = v.z /r;
  return(w);
}

Vectors vsubtract ( Vectors v1, Vectors v2 ) {
  Vectors output;
  output.x = v1.x - v2.x;
  output.y = v1.y - v2.y;
  output.z = v1.z - v2.z;
  return output;
}

Vectors vadd ( Vectors v1, Vectors v2 ) {
  Vectors output;
  output.x = v1.x + v2.x;
  output.y = v1.y + v2.y;
  output.z = v1.z + v2.z;
  return output;
}

Vectors vrunning_avg ( Vectors thisVector1, Vectors thisVector2, int thisNum ) {
  Vectors ravg;
  ravg.x = (thisVector1.x * thisNum + thisVector2.x) / (1 + thisNum);
  ravg.y = (thisVector1.y * thisNum + thisVector2.y) / (1 + thisNum);
  ravg.z = (thisVector1.z * thisNum + thisVector2.z) / (1 + thisNum);
  return(ravg);
}

Vectors qvq ( Vectors thisQuat, Vectors thisPosition ) {
  // Fewer computations than the other definition
  // Also, note that we're assuming our quaternions are unit quaternions (they should be, or should be very close). 

  thisPosition.r = 0;
  float qv_r = thisQuat.r * thisPosition.r - thisQuat.x * thisPosition.x - thisQuat.y * thisPosition.y - thisQuat.z * thisPosition.z ;
  float qv_x = thisQuat.r * thisPosition.x + thisQuat.x * thisPosition.r + thisQuat.y * thisPosition.z - thisQuat.z * thisPosition.y ; 
  float qv_y = thisQuat.r * thisPosition.y - thisQuat.x * thisPosition.z + thisQuat.y * thisPosition.r + thisQuat.z * thisPosition.x ;
  float qv_z = thisQuat.r * thisPosition.z + thisQuat.x * thisPosition.y - thisQuat.y * thisPosition.x + thisQuat.z * thisPosition.r ; 

  // This is the same operation as above, expect the sign is flipped everywhere we have a non-real quat component (since the second q in qvq is the inverse of q). 
  //float qvq_r = qv_r * thisQuat.r + qv_x * thisQuat.x + qv_y * thisQuat.y + qv_z * thisQuat.z ; // no point in computing r value of the output
  float qvq_x = qv_r * thisQuat.x * -1.0 + qv_x * thisQuat.r - qv_y * thisQuat.z + qv_z * thisQuat.y ;
  float qvq_y = qv_r * thisQuat.y * -1.0 + qv_x * thisQuat.z + qv_y * thisQuat.r - qv_z * thisQuat.x ;
  float qvq_z = qv_r * thisQuat.z * -1.0 - qv_x * thisQuat.y + qv_y * thisQuat.x + qv_z * thisQuat.r ;

  Vectors thisFinalPosition;
  thisFinalPosition.r = 0; // qvq_r ; // physically meaningful vectors have r = 0
  thisFinalPosition.x = qvq_x ;
  thisFinalPosition.y = qvq_y ;
  thisFinalPosition.z = qvq_z ;

  return(thisFinalPosition);

}