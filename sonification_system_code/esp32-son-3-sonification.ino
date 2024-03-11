int Son_pitch( float d ) {

  int TF_temp = (int) ( (float)PitchValueMax * expf( -1.0 * ( d / a_pitch ) ) ); 

  return(TF_temp);

}

int Son_amp( float er ) {

  int amp_temp = (int) ( (float)amp_max * expf( -1.0 * ( er / a_amp ) ) );

  return(amp_temp); 

}

void Playback_model( void ) {

  for ( int i = 0; i <= MOTsamplecount_max; i++ ) {
    #ifdef SENSOR3
      float to_end = vdst4(model_q1[i],model_q1[MOTsamplecount_max]) + 
                     vdst4(model_q2[i],model_q2[MOTsamplecount_max]) + 
                     vdst4(model_q3[i],model_q3[MOTsamplecount_max]);
    #else 
      float to_end = vdst4(model_q1[i],model_q1[MOTsamplecount_max]) + 
                     vdst4(model_q2[i],model_q2[MOTsamplecount_max]);
    #endif
    amplitude_update(amp_max);
    UpdateTF(Son_pitch(to_end));
    delay(1); 
  }

  amplitude_update(0);
  //delaywithFIFOreset(1); // Not needed since we're not using FIFO
  
}

void Playback_reach( void ) {

  for ( int i = 0; i < MaxReadingsToSave; i++ ) {
        if ( i <= (Realsamplecount - 1) ) {

          float error1 = vdst4(quat1_saved[i],model_q1[MOTsamplecount_saved[i]]); // max error = 2
          float error2 = vdst4(quat2_saved[i],model_q2[MOTsamplecount_saved[i]]); // max error = 2
          #ifdef SENSOR3
            float error3 = vdst4(quat3_saved[i],model_q3[MOTsamplecount_saved[i]]); // max error = 2
            float error_total = error1 + error2 + erro3; // max of 6
            float to_end = vdst4(quat1_saved[i],model_q1[MOTsamplecount_max]) + 
                          vdst4(quat2_saved[i],model_q2[MOTsamplecount_max]) + 
                          vdst4(quat3_saved[i],model_q3[MOTsamplecount_max]);
          #else 
            float error_total = error1 + error2; // max of 4
            float to_end = vdst4(quat1_saved[i],model_q1[MOTsamplecount_max]) + 
                          vdst4(quat2_saved[i],model_q2[MOTsamplecount_max]);
          #endif

          amplitude_update(Son_amp(error_total));
          UpdateTF(Son_pitch(to_end));
          delay(1); 
          
        }
      }
  
  amplitude_update(0);

}

void UpdateTF( int new_pitch ) {
  new_pitch = constrain(new_pitch,PitchValueMin,PitchValueMax);
  timerSW_top = (int) timerSizeprescaler / (new_pitch * 2);
  timerAlarmWrite(timerSW, timerSW_top, true); // true = reload automatically
}

void amplitude_update( int thisamplitude ){
  thisamplitude = constrain(thisamplitude,0,amp_max);
  amp = thisamplitude;
}