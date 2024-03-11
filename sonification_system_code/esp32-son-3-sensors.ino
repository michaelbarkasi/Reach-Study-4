void FetchReadingsSensor (ICM_20948_SPI &myICM, Vectors &storeAccReadingInThisVector, Vectors &storeGyroReadingInThisVector) { 
  
  myICM.getAGMT(); // get the A, G, M, and T readings
  // Scaled data (1000dps gyro, 16g accel range)

  storeAccReadingInThisVector.x = myICM.agmt.acc.axes.x * Acc_multiple; // convert first to g by dividing by 2048, then to m/s^2 by multiplying by 9.80665, i.e. (myICM.agmt.acc.axes.x / 2048.0) * 9.80665
  storeAccReadingInThisVector.y = myICM.agmt.acc.axes.y * Acc_multiple;
  storeAccReadingInThisVector.z = myICM.agmt.acc.axes.z * Acc_multiple;

  storeGyroReadingInThisVector.x = myICM.agmt.gyr.axes.x * Gyro_multiple; // convert to degrees/s by dividing by 32.8, then to radians by multiplying by * (PI/180.0), i.e. (myICM.agmt.gyr.axes.x / 32.8) * (PI/180.0)
  storeGyroReadingInThisVector.y = myICM.agmt.gyr.axes.y * Gyro_multiple;
  storeGyroReadingInThisVector.z = myICM.agmt.gyr.axes.z * Gyro_multiple; 

  /* storeMagReadingInThisVector.x = myICM.agmt.mag.axes.x;
  storeMagReadingInThisVector.y = myICM.agmt.mag.axes.y;
  storeMagReadingInThisVector.z = myICM.agmt.mag.axes.z; */
 
}

void initializeSensor_SPI (ICM_20948_SPI &myICM, int thisPin, int thisNum ) {

  bool initialized = false;
  while (!initialized) {

    myICM.begin(thisPin, SPI_PORT, SPI_FREQ); // Here we are using the user-defined SPI_FREQ as the clock speed of the SPI bus

    SERIAL_PORT.print(F("Initialization of sensor "));
    SERIAL_PORT.print(thisNum,1);
    SERIAL_PORT.print(F(" returned: "));
    SERIAL_PORT.println(myICM.statusString());
    if (myICM.status != ICM_20948_Stat_Ok)
    {
      SERIAL_PORT.println("Trying again...");
      delay(500);
    }
    else
    {
      initialized = true;
    }
  }
  SERIAL_PORT.print(F("Device "));
  SERIAL_PORT.print(thisNum,1);
  SERIAL_PORT.println(F(" connected!"));

}

void configureInterrupts (ICM_20948_SPI &myICM) {

  // Now we're going to set up interrupts. There are a lot of options, but for this test we're just configuring the interrupt pin and enabling interrupts to tell us when new data is ready
  /*
    ICM_20948_Status_e  cfgIntActiveLow         ( bool active_low );
    ICM_20948_Status_e  cfgIntOpenDrain         ( bool open_drain );
    ICM_20948_Status_e  cfgIntLatch             ( bool latching );                          // If not latching then the interrupt is a 50 us pulse
    ICM_20948_Status_e  cfgIntAnyReadToClear    ( bool enabled );                           // If enabled, *ANY* read will clear the INT_STATUS register. So if you have multiple interrupt sources enabled be sure to read INT_STATUS first
    ICM_20948_Status_e  cfgFsyncActiveLow       ( bool active_low );
    ICM_20948_Status_e  cfgFsyncIntMode         ( bool interrupt_mode );                    // Can ue FSYNC as an interrupt input that sets the I2C Master Status register's PASS_THROUGH bit
    ICM_20948_Status_e  intEnableI2C            ( bool enable );
    ICM_20948_Status_e  intEnableDMP            ( bool enable );
    ICM_20948_Status_e  intEnablePLL            ( bool enable );
    ICM_20948_Status_e  intEnableWOM            ( bool enable );
    ICM_20948_Status_e  intEnableWOF            ( bool enable );
    ICM_20948_Status_e  intEnableRawDataReady   ( bool enable );
    ICM_20948_Status_e  intEnableOverflowFIFO   ( uint8_t bm_enable );
    ICM_20948_Status_e  intEnableWatermarkFIFO  ( uint8_t bm_enable );
 */
  myICM.cfgIntActiveLow(true);  // Active low to be compatible with the breakout board's pullup resistor
  myICM.cfgIntOpenDrain(false); // Push-pull, though open-drain would also work thanks to the pull-up resistors on the breakout
  myICM.cfgIntLatch(true);      // Latch the interrupt until cleared; Not sure why, but with high speed DMP does not work if we latch; must set to false. 
  SERIAL_PORT.print(F("cfgIntLatch returned: "));
  SERIAL_PORT.println(myICM.statusString());

  myICM.intEnableRawDataReady(true); // enable interrupts on the DMP
  SERIAL_PORT.print(F("intEnableRawDataReady returned: "));
  SERIAL_PORT.println(myICM.statusString());

  //  // Note: weirdness with the Wake on Motion interrupt being always enabled.....
  //  uint8_t zero_0 = 0xFF;
  //  ICM_20948_execute_r( &myICM._device, AGB0_REG_INT_ENABLE, (uint8_t*)&zero_0, sizeof(uint8_t) );
  //  SERIAL_PORT.print("INT_EN was: 0x"); SERIAL_PORT.println(zero_0, HEX);
  //  zero_0 = 0x00;
  //  ICM_20948_execute_w( &myICM._device, AGB0_REG_INT_ENABLE, (uint8_t*)&zero_0, sizeof(uint8_t) );

}

void configureSensor (ICM_20948_SPI &myICM) {

// Here we are doing a SW reset to make sure the device starts in a known state
  myICM.swReset();
  if (myICM.status != ICM_20948_Stat_Ok)
  {
    SERIAL_PORT.print(F("Software Reset returned: "));
    SERIAL_PORT.println(myICM.statusString());
  }
  delay(250);

  // Now wake the sensor up
  myICM.sleep(false);
  myICM.lowPower(false);

  // The next few configuration functions accept a bit-mask of sensors for which the settings should be applied.

  // Set Gyro and Accelerometer to a particular sample mode
  // options: ICM_20948_Sample_Mode_Continuous
  //          ICM_20948_Sample_Mode_Cycled
  myICM.setSampleMode((ICM_20948_Internal_Acc | ICM_20948_Internal_Gyr), ICM_20948_Sample_Mode_Continuous);
  SERIAL_PORT.print(F("setSampleMode returned: "));
  SERIAL_PORT.println(myICM.statusString());

  //  ICM_20948_smplrt_t mySmplrt;
  //  mySmplrt.g = 54;
  //  myICM.setSampleRate(ICM_20948_Internal_Gyr, mySmplrt);
  //  SERIAL_PORT.print(F("setSampleRate returned: "));
  //  SERIAL_PORT.println(myICM.statusString());

  // Set full scale ranges for acc and gyr; limit gyr to dps500 (users are not rotating their body faster than that)
  ICM_20948_fss_t myFSS; // This uses a "Full Scale Settings" structure that can contain values for all configurable sensors

  myFSS.a = gpm16; // (ICM_20948_ACCEL_CONFIG_FS_SEL_e)
                  // gpm2
                  // gpm4
                  // gpm8
                  // gpm16

  myFSS.g = dps1000; // (ICM_20948_GYRO_CONFIG_1_FS_SEL_e)
                    // dps250
                    // dps500
                    // dps1000
                    // dps2000

  myICM.setFullScale((ICM_20948_Internal_Acc | ICM_20948_Internal_Gyr), myFSS);
  if (myICM.status != ICM_20948_Stat_Ok)
  {
    SERIAL_PORT.print(F("setFullScale returned: "));
    SERIAL_PORT.println(myICM.statusString());
  }

  // Set gyro sample rate divider with GYRO_SMPLRT_DIV
  // Set accel sample rate divider with ACCEL_SMPLRT_DIV_2
  ICM_20948_smplrt_t mySmplrt;
  //mySmplrt.g = 19; // ODR is computed as follows: 1.1 kHz/(1+GYRO_SMPLRT_DIV[7:0]). 19 = 55Hz. InvenSense Nucleo example uses 19 (0x13).
  //mySmplrt.a = 19; // ODR is computed as follows: 1.125 kHz/(1+ACCEL_SMPLRT_DIV[11:0]). 19 = 56.25Hz. InvenSense Nucleo example uses 19 (0x13).
  #ifdef RUNHALFSPEED
    mySmplrt.g = 4; // 225Hz
    mySmplrt.a = 4; // 225Hz
  #endif
  //mySmplrt.g = 1; // 550Hz // I would have assumed this should be 2, but with 2, only get 377 or 378 updates/sec; with 1, get 565, 566 --- Is it because this number is an exponent on divider, e.g. 1100 / 2^x?
  //mySmplrt.a = 1; // 562.5Hz // I would have assumed this should be 2, but with 2, only get 377 or 378 updates/sec; with 1, get 565, 566
  //mySmplrt.g = 8; // 112Hz
  //mySmplrt.a = 8; // 112Hz
  #ifdef RUNFULLSPEED 
    mySmplrt.g = 0; // 1100Hz 
    mySmplrt.a = 0; // 1125Hz
  #endif
  myICM.setSampleRate((ICM_20948_Internal_Acc | ICM_20948_Internal_Gyr), mySmplrt); 
  SERIAL_PORT.print(F("setSampleRate (Acc and Gyro) returned: "));
  SERIAL_PORT.println(myICM.statusString());

  // ... don't want or need the low-pass filter (we want to detect fast changes), at least, crank the filter down as low as it goes (499bw/376.5)
  // Assume these filters are running averages of some sort, e.g. nyquist bandwidth of 265Hz means a running average (of faster samples) is returned 265 times/sec
  // Set up Digital Low-Pass Filter configuration
  ICM_20948_dlpcfg_t myDLPcfg;    // Similar to FSS, this uses a configuration structure for the desired sensors
  myDLPcfg.a = acc_d473bw_n499bw; // (ICM_20948_ACCEL_CONFIG_DLPCFG_e)
                                  // acc_d246bw_n265bw        - means 3db bandwidth is 246 hz and nyquist bandwidth is 265 hz
                                  // acc_d111bw4_n136bw    *
                                  // acc_d50bw4_n68bw8
                                  // acc_d23bw9_n34bw4
                                  // acc_d11bw5_n17bw
                                  // acc_d5bw7_n8bw3        - means 3 db bandwidth is 5.7 hz and nyquist bandwidth is 8.3 hz
                                  // acc_d473bw_n499bw

  myDLPcfg.g = gyr_d361bw4_n376bw5; // (ICM_20948_GYRO_CONFIG_1_DLPCFG_e)
                                    // gyr_d196bw6_n229bw8   
                                    // gyr_d151bw8_n187bw6   *
                                    // gyr_d119bw5_n154bw3
                                    // gyr_d51bw2_n73bw3
                                    // gyr_d23bw9_n35bw9
                                    // gyr_d11bw6_n17bw8
                                    // gyr_d5bw7_n8bw9
                                    // gyr_d361bw4_n376bw5

  myICM.setDLPFcfg((ICM_20948_Internal_Acc | ICM_20948_Internal_Gyr), myDLPcfg);
  if (myICM.status != ICM_20948_Stat_Ok)
  {
    SERIAL_PORT.print(F("setDLPcfg returned: "));
    SERIAL_PORT.println(myICM.statusString());
  }

  // Choose whether or not to use DLPF
  // Here we're also showing another way to access the status values, and that it is OK to supply individual sensor masks to these functions
  ICM_20948_Status_e accDLPEnableStat = myICM.enableDLPF(ICM_20948_Internal_Acc, false);
  ICM_20948_Status_e gyrDLPEnableStat = myICM.enableDLPF(ICM_20948_Internal_Gyr, false);
  SERIAL_PORT.print(F("Enable/disable DLPF for Accelerometer returned: "));
  SERIAL_PORT.println(myICM.statusString(accDLPEnableStat));
  SERIAL_PORT.print(F("Enable/disable DLPF for Gyroscope returned: "));
  SERIAL_PORT.println(myICM.statusString(gyrDLPEnableStat));

  // Choose whether or not to start the magnetometer
  //  myICM.startupMagnetometer();
  //  if (myICM.status != ICM_20948_Stat_Ok)
  //  {
  //    SERIAL_PORT.print(F("startupMagnetometer returned: "));
  //    SERIAL_PORT.println(myICM.statusString());
  //  }

}