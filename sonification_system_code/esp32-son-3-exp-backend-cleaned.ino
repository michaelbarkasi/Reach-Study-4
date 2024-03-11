/*
 * Code for Wave APS Prototype 4
 *  
 * 3-Sensor Motion Sonification via ttPWM and 3D position processing for multi-pivot body movements
 *
 * by Michael Barkasi
 * copyright 2023
 * 
 * For use with three TDK InvenSense ICM-20948 inertial sensors (SparkFun breakouts ID 15335) (SPI), and 
 *  ESP32 feather (Adafruit HUZZAH32, ID 3591), mounted on a custom PCB feather.
 *
 * Requires the SparkFun libary ICM_20948.h and board manager for Adafruit ESP32 feather.
 * 
 **************************************************************
 * Some code related to SPI/Interrupts for Sen20948 taken from: 
 * 
 * Example3_Interrupts.ino
 * ICM 20948 Arduino Library Demo
 * Builds on Example2_Advanced.ino to set up interrupts when data is ready
 * Owen Lyke @ SparkFun Electronics
 * Original Creation Date: June 5 2019
 *
 * Some code related to the setup for Sen20948 taken from (and now modified): 
 * 
 * Example2_Advanced.ino
 * ICM 20948 Arduino Library Demo
 * Shows how to use granular configuration of the ICM 20948
 * Owen Lyke @ SparkFun Electronics
 * Original Creation Date: April 17 2019
 * 
 * Example7_DMP_Quat6_EulerAngles.ino
 * ICM 20948 Arduino Library Demo
 * Initialize the DMP based on the TDK InvenSense ICM20948_eMD_nucleo_1.0 example-icm20948
 * Paul Clark, April 25th, 2021
 * Based on original code by:
 * Owen Lyke @ SparkFun Electronics
 * Original Creation Date: April 17 2019
 * 
 * Also, the redefine of "CM_20948_Status_e ICM_20948::initializeDMP(void)" in SensorSetup is from: 
 * 
 * Example10_DMP_FastMultipleSensors.ino
 * ICM 20948 Arduino Library Demo
 * Initialize the DMP based on the TDK InvenSense ICM20948_eMD_nucleo_1.0 example-icm20948
 * Paul Clark, April 25th, 2021
 * Based on original code by:
 * Owen Lyke @ SparkFun Electronics
 * Original Creation Date: April 17 2019
 * 
 *  From SparkFun on examples 7 and 10: 
 *  ** This example is based on InvenSense's _confidential_ Application Note "Programming Sequence for DMP Hardware Functions".
 *  ** We are grateful to InvenSense for sharing this with us.
 *  
 * License for this code contained in License.md. This license is the MIT License, which is found here: https://github.com/sparkfun/SparkFun_ICM-20948_ArduinoLibrary/blob/main/License.md
 *      and here: https://opensource.org/licenses/MIT
 *      
 *      The relevant part: "Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
 *                          files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
 *                          modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
 *                          is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included 
 *                          in all copies or substantial portions of the Software."
 * 
 * ***************************************************************
 * Sound Generation: 
 *
 * Setup base wave form @ 152,380 Hz using the ESP32's hardware PWM
 *  152,380 Hz is the fastest oscilation its hardware PWM can achieve with at least 9 bit of resolution
 *  Get this with: clk_src = LEDC_APB_CLK (80 MHz).
 *  ... then duty resolution is integer (log 2 (LEDC_APB_CLK / frequency))
 *     see: 
 *       https://docs.espressif.com/projects/esp-idf/en/latest/esp32/api-reference/peripherals/ledc.html
 *       https://github.com/espressif/esp-idf/blob/5893797bf068ba6f72105fff289ead370b4591a3/examples/peripherals/ledc/ledc_basic/main/ledc_basic_example_main.c
 * Initial code to get PWM working from: https://randomnerdtutorials.com/esp32-pwm-arduino-ide/
 * Note: We really want this PWM to run in LEDC_HIGH_SPEED_MODE; the code from the esp-idf gives this control, but not the code for arduino;
 *  Assuming arduino puts PWM into high-speed mode by default --- I have confirmed with this with oscilloscope readings.
 *
 *  How the sound is generated: Two timers are used. The first (PWM) controls a base wave form of digital pulses with ultrasonic (e.g., 152 kHz) oscilation. 
 *   E.g., 152 kHz is too fast for a speaker, which (acting as a low-pass filter) only "sees" the average voltage, which is digital HIGH * duty cycle. 
 *   (Duty cycle is pulse width, e.g. pulse on for 70% of the period and off for 30%.) The second timer oscilates pulse width (duty cycle) of this wave at a hearable frequency, 
 *   flipping between 0 and some second nonzero value. This second nonzero value essentially ends up controlling volume, since it covaries with the voltage "seen" on each pulse by the speaker.
 *   
 * ***************************************************************
 * BlueTooth: 
 *  
 * Based on the Arduino SerialToSerialBT example, which came with the following preface:
 * 
 * //This example code is in the Public Domain (or CC0 licensed, at your option.)
 * //By Evandro Copercini - 2018
 * //
 * //This example creates a bridge between Serial and Classical Bluetooth (SPP)
 * //and also demonstrate that SerialBT have the same functionalities of a normal Serial
 * 
*/

// **************************************************************
// Basic Parameters: 

const int update_rate = 1000; // Controls both sound updates and when to process sensor data/ update motion models
const float dt = 1.0 / update_rate;

const bool TERMINAL = true; // play feedback online (terminal = false) or after reach is completed (terminal = true)?
const int main_block_count = 75;

//#define DEMO // uncomment to switch to just 5 random motions before sonification (first is model)
//#define SENSOR3 // uncomment to use sensor 3 as well (currently, unit cannot process all 3 sensors at update_rate = 1000)
//#define RUNHALFSPEED // uncomment to run sensors at 225 Hz (does not change sonification update rate)
#define RUNFULLSPEED // uncomment to run sensors at 1100 Hz (gyro) and 1125 Hz (accel) 
#define SUBTRACTGYROBIAS // uncomment to subtract estimated gyro bias from raw gyro readings
//#define USEBT // uncomment to communicate over BlueTooth instead of wired serial port
//#define DACwave // uncomment to use DAC instead of ttPWM for sound generation

// **************************************************************
// Sensor Communications (SPI and Interrupts): 

#include "ICM_20948.h" // Click here to get the library: http://librarymanager/All#SparkFun_ICM_20948_IMU
// Definitions for Interrupts
#define INT_PIN_myICM1 33//13 
#define INT_PIN_myICM2 32//12
// Definitions for SPI
#define SPI_PORT SPI     
#define SPI_FREQ 4500000 // 7000000 is max for Sen20948 (see datasheet)? other TDK sheets suggest 2.5M. 4500000 seems to work. 
#define CS_PIN_myICM1 27//33     
#define CS_PIN_myICM2 21//15
ICM_20948_SPI myICM1; // Arbitrary name for the sensor
ICM_20948_SPI myICM2; // Arbitrary name for the sensor
// Variables for interrupts (pulling sensor data)
volatile bool isrFired1 = false;
volatile bool isrFired2 = false;

#ifdef SENSOR3 
  #define INT_PIN_myICM3 14//27
  #define CS_PIN_myICM3 15//32  
  ICM_20948_SPI myICM3; // Arbitrary name for the sensor
  volatile bool isrFired3 = false;
#endif

// **************************************************************
// Variables controlling base pulse wave production / volume control: 

#ifdef DACwave
  const int soundPin = A1;  // DAC is set to use DAC1, which is A1; DAC2 is A0
  const int amp_max_true = 255; // DAC is 8-bit, so only 256 levels
#else 
  const int soundPin = 25; // pin G (GPIO) 25 is actually the same pin on breakout board as A1
  const int base_wave_freq = 152380; // highest availabe with 9 bit resolution is 152380
  const int pwm_Channel = 0;
  const int pwm_resolution = (int) log2( 80000000 / base_wave_freq ) ;
  const int amp_max_true = pow(2,pwm_resolution) - 1; // for 9 bit resolution: 2^9 - 1 = 511; for 8 bit 2^8 - 1 = 255
#endif

#ifdef DEMO
  const int amp_max = (int) (amp_max_true*0.5); // Make it quiter for the demo
#else
  const int amp_max = amp_max_true;
#endif

// **************************************************************
// More Sound Control: 

// Basic audio parameters
volatile int amp = (int) amp_max;
const int amp_default = (int) amp_max;
const int pitchinitial = 440;
//volatile int TF = pitchinitial;
const int PitchValueMax = 6000;
const int PitchValueMin = 220;
float a_amp = 1.0;
float a_pitch = 1.0;

// setting the timer for generating hearable oscilation 
hw_timer_t * timerSW = NULL;
portMUX_TYPE timerSWMux = portMUX_INITIALIZER_UNLOCKED;
const int timerSWprescaler = 2; // Timers are 80Mhz; counters seem to be at least 32bit, so, don't need to prescale down low
const long timerSize = 80000000;
const int timerSizeprescaler = timerSize/timerSWprescaler;
volatile int timerSW_top = (int) timerSizeprescaler / (pitchinitial * 2); // controls frequency of hearable oscilation
volatile bool up = true;

// Variables for the timer controlling sound updates / motion processing
hw_timer_t * timerSU = NULL;
portMUX_TYPE timerSUMux = portMUX_INITIALIZER_UNLOCKED;
const int timerSU_top = (int) timerSizeprescaler / update_rate; 

// **************************************************************
// Setup for Bluetooth:

#ifdef USEBT 
  #include "BluetoothSerial.h"
  #if !defined(CONFIG_BT_ENABLED) || !defined(CONFIG_BLUEDROID_ENABLED)
  #error Bluetooth is not enabled! Please run `make menuconfig` to and enable it
  #endif
  BluetoothSerial SerialBT;
  #define SERIAL_PORT SerialBT 
#else
  #define SERIAL_PORT Serial
#endif

// **************************************************************
// Motion Processing

// Core Vector Variables
#ifdef SENSOR3 
  struct Vectors {
    float r = 0.0; // Uncomment if we need to use quaternions, in which case this is q0, hence, it should initialize as zero, in case it doesn't get set elsewhere and we're dealing with a vector in R^3
    float x = 0.0; // q1; if ever encoding quaternions, think of x/y/z as q1/q2/q3.  
    float y = 0.0; // q2; it doesn't matter what we initialize these values as, they all get overwritten. 0/1/0/0 is a holdover from previous programs that used quaternions.
    float z = 0.0; // q3
  } Sensor1_Acc, Sensor1_Gyro, // Sensor1_Mag, 
    Sensor2_Acc, Sensor2_Gyro, // Sensor2_Mag,
    Sensor3_Acc, Sensor3_Gyro, // Sensor3_Mag,
    Sensor1_Rot, Sensor2_Rot, Sensor3_Rot,
    Sensor1_Quat, Sensor2_Quat, Sensor3_Quat,
    Sensor1_GravityIP, Sensor2_GravityIP, Sensor3_GravityIP,
    Sensor1_Xaxis, Sensor1_Yaxis, Sensor1_Zaxis,
    Sensor2_Xaxis, Sensor2_Yaxis, Sensor2_Zaxis,
    Sensor3_Xaxis, Sensor3_Yaxis, Sensor3_Zaxis,
    Sensor1_XaxisIP, Sensor1_YaxisIP, Sensor1_ZaxisIP,
    Sensor2_XaxisIP, Sensor2_YaxisIP, Sensor2_ZaxisIP,
    Sensor3_XaxisIP, Sensor3_YaxisIP, Sensor3_ZaxisIP,
    Sensor1_GyroBias, Sensor2_GyroBias, Sensor3_GyroBias,
    zerov;
#else
  struct Vectors {
    float r = 0.0; // Uncomment if we need to use quaternions, in which case this is q0, hence, it should initialize as zero, in case it doesn't get set elsewhere and we're dealing with a vector in R^3
    float x = 0.0; // q1; if ever encoding quaternions, think of x/y/z as q1/q2/q3.  
    float y = 0.0; // q2; it doesn't matter what we initialize these values as, they all get overwritten. 0/1/0/0 is a holdover from previous programs that used quaternions.
    float z = 0.0; // q3
  } Sensor1_Acc, Sensor1_Gyro, // Sensor1_Mag, 
    Sensor2_Acc, Sensor2_Gyro, // Sensor2_Mag
    Sensor1_Rot, Sensor2_Rot, 
    Sensor1_Quat, Sensor2_Quat,
    Sensor1_GravityIP, Sensor2_GravityIP,
    Sensor1_Xaxis, Sensor1_Yaxis, Sensor1_Zaxis,
    Sensor2_Xaxis, Sensor2_Yaxis, Sensor2_Zaxis,
    Sensor1_XaxisIP, Sensor1_YaxisIP, Sensor1_ZaxisIP,
    Sensor2_XaxisIP, Sensor2_YaxisIP, Sensor2_ZaxisIP,
    Sensor1_GyroBias, Sensor2_GyroBias,
    zerov;
#endif

// Motion Processing Variables
const float Gyro_multiple = PI / (32.8 * 180.0); // puts gyro readings in radians / s; depends on sensor full-range settings!!
const float Acc_multiple = 9.80665 / 2048.0; // puts accel readings in m/s^2; depends on sensor full-range settings!!
float high_gyro_readings = 0;

// Reach Number
int reachnum = 1; 
int randomsamplenum = 0; // which of the first random motions is to be treated as the model during sonification? (Set randomly)
#ifdef DEMO
 int randomsamplesize = 5; // How many random samples should be collected before sonifying? 
#else
 int randomsamplesize = 25; // How many random samples should be collected before sonifying? 
#endif

// Initial Position Detection 
bool IPfound = false;
const float IPtolerancefraction = 0.1; // This fraction used to set the above value, e.g. ".1 = within 10 percent of grav vector magnitude". 
float IPtolerance = 1.0; // How close to the gravity vectors recorded at the start do we need to be for the unit to think user is in the IP? (Set below.)

// Tracking Unit Sample Rate
unsigned long time1 = 0; // time1, time2, and time3 used to time how long we've been recording, during live sonification
unsigned long time2 = 0;
unsigned long time3 = 0;
float samplerate = 1.0;

// Tracking Sample Number and Model (MOT) Number
int Realsamplecount = 0; // because of the compensation for different velocity, this comes apart from the MOTsamplecount
int MOTsamplecount = 0;
int MOTsamplecount_max = 0;
int MOTadvance = 1; // used in real-time time-warping algorithm

// Error Measurements 
float to_end_start = 1.0;
float error1 = 0;
float error2 = 0;
#ifdef SENSOR3
  float error3 = 0;
#endif
float error_total = 0;

// Variables to Store Saved Data
const int MaxReadingsToSave = 2500;
struct Vectors model_q1[MaxReadingsToSave];
struct Vectors model_q2[MaxReadingsToSave];
Vectors* quat1_saved;
Vectors* quat2_saved;
#ifdef SENSOR3
  struct Vectors model_q3[MaxReadingsToSave];
  Vectors* quat3_saved;
#endif
int MOTsamplecount_saved[MaxReadingsToSave] = {0};

// Unit State Tracking Variables 
bool FIRSTSAMPLE = true;
bool IN_MOTION = false;
bool WAITING_FOR_REST = true;
bool WAITING_FOR_INITATION = true;
bool MOTION_JUST_ENDED = false;
bool MOTION_JUST_STARTED = false;

// Motion Detection Variables
int jostle_buff_count = 0;
int jostle_buff_count_saved = 1;
int jostle_noise_filter = 100;
int rest_buff_count = 0;
int rest_noise_filter = 100; 
#define MOTIONCHECKSENSOR_Gyro Sensor2_Gyro // What sensor are we using for the motion check?
const float motionthreshold_gyro = (12.0 * (PI/180.0)) * (12.0 * (PI/180.0)); // How fast (in radians per second) must a sensor be moving (i.e., rotating) for it to be considered "in motion"?
  // Number here ^ is magnetude of the *total* rotation <x,y,z> vector, not rotation about any one axis. (3 degrees per axis is standard noise)
const float motionthreshold_gyro_high = (motionthreshold_gyro * 3.0) * (motionthreshold_gyro * 3.0);

// **************************************************************
// Interrupt functions: 

bool UPDATE_SOUND_NOW = true; // The timer (TC2) controlling sound update rate will flip this on when it's time to read.

void IRAM_ATTR onTimerSU() { // interrupt handler for timer generating sound updates / motion processing: 
  portENTER_CRITICAL_ISR(&timerSWMux);
  UPDATE_SOUND_NOW = true;
  portEXIT_CRITICAL_ISR(&timerSWMux);
}
void IRAM_ATTR onTimerSW() { // for timer generating hearable oscilation
  portENTER_CRITICAL_ISR(&timerSWMux);
  #ifdef DACwave
    if (up) {
      dacWrite(DAC1, amp);
    } else {
      dacWrite(DAC1, 0);
    }
  #else 
    if (up) {
      ledcWrite(pwm_Channel, amp);
    } else {
      ledcWrite(pwm_Channel, 0);
    }
  #endif
  up = !up;
  portEXIT_CRITICAL_ISR(&timerSWMux);
}
void icmISR1(void) { // These three for grabbing sensor data
  isrFired1 = true; // Can't use I2C within ISR on 328p, so just set a flag to know that data is available
}
void icmISR2(void) {
  isrFired2 = true; // Can't use I2C within ISR on 328p, so just set a flag to know that data is available
}
#ifdef SENSOR3
  void icmISR3(void) {
    isrFired3 = true; // Can't use I2C within ISR on 328p, so just set a flag to know that data is available
  }
#endif

// **************************************************************

void setup() {

  // Setup serial port / BlueTooth
  Serial.begin(115200);
  #ifdef USEBT
    SerialBT.begin("Wave APS (ESP32) Sonification Unit"); // BlueTooth device name
    Serial.println("The device started, now you can pair it with bluetooth!");
  #else 
    Serial.println("The device started!");
  #endif

  // Allocate memory to store online readings: 
  quat1_saved = (Vectors*)malloc(MaxReadingsToSave * sizeof(Vectors));
  if (quat1_saved == nullptr) {
    Serial.println("Heap memory allocation failure!");
    while (true) {
      // An error occurred; loop indefinitely
    }
  }
  quat2_saved = (Vectors*)malloc(MaxReadingsToSave * sizeof(Vectors));
  if (quat2_saved == nullptr) {
    Serial.println("Heap memory allocation failure!");
    while (true) {
      // An error occurred; loop indefinitely
    }
  }
  #ifdef SENSOR3
    quat3_saved = (Vectors*)malloc(MaxReadingsToSave * sizeof(Vectors));
    if (quat3_saved == nullptr) {
      Serial.println("Heap memory allocation failure!");
      while (true) {
        // An error occurred; loop indefinitely
      }
    }
  #endif

  // Pause all setup until connected to serial / BlueTooth
  while (!SERIAL_PORT.available()) {
    SERIAL_PORT.println(F("Press any key to begin."));
    delay(6000);
  }
  #ifdef USEBT 
    SERIAL_PORT.println(F("Connected to serial port through Bluetooth!"));
  #else 
    SERIAL_PORT.println(F("Connected to serial port!"));
  #endif

  if ( TERMINAL ) {
    SERIAL_PORT.println();
    SERIAL_PORT.println(F("Program: Terminal Sonification (Online PostBlock)"));
  } else {
    SERIAL_PORT.println();
    SERIAL_PORT.println(F("Program: Online Sonification (Terminal PostBlock)"));
  }

  SERIAL_PORT.println();
  SERIAL_PORT.print(F("Update rate selected: "));
  SERIAL_PORT.println(update_rate,1);

  SERIAL_PORT.println();
  SERIAL_PORT.print(F("Setting up base pulse wave: "));
  
  // Setup base pulse wave
  #ifdef DACwave
    SERIAL_PORT.println(F("DAC selected"));
    SERIAL_PORT.print(F("Volume resolution: "));
    SERIAL_PORT.println(amp_max,1);
  #else
    SERIAL_PORT.println(F("ttPWM selected"));
    SERIAL_PORT.print(F("Volume resolution: "));
    SERIAL_PORT.println(amp_max);
    SERIAL_PORT.print(F("Base wave frequency: "));
    SERIAL_PORT.println(base_wave_freq,1);
    ledcSetup(pwm_Channel, base_wave_freq, pwm_resolution); // configure PWM functionalitites
    ledcAttachPin(soundPin, pwm_Channel); // attach the channel to the pin generating the wave
  #endif
  SERIAL_PORT.println();
  
  // Setup timer used for hearable oscilation
  timerSW = timerBegin(0, timerSWprescaler, true); // true indicates counter goes up
  timerAttachInterrupt(timerSW, &onTimerSW, true); // true indicates edge interrupt (false for level)
  timerAlarmWrite(timerSW, timerSW_top, true); // true for reload automatically
  timerAlarmEnable(timerSW);
  
  // Setup timer for sound updates/ motion processing
  timerSU = timerBegin(1, timerSWprescaler, true); // true indicates counter goes up
  timerAttachInterrupt(timerSU, &onTimerSU, true); // true indicates edge interrupt (false for level)
  timerAlarmWrite(timerSU, timerSU_top, true); // true for reload automatically
  timerAlarmEnable(timerSU);

  // Enable SPI communication 
  SPI_PORT.begin();

  // Configure Sensors: 

  myICM1.enableDebugging(); // Uncomment this line to enable helpful debug messages on Serial about ICM1
  myICM2.enableDebugging(); // Uncomment this line to enable helpful debug messages on Serial about ICM2

  // Setup Interrupts
  // myICM1
  pinMode(INT_PIN_myICM1, INPUT_PULLUP); // Using a pullup b/c ICM-20948 Breakout board has an onboard pullup as well and we don't want them to compete
  attachInterrupt(digitalPinToInterrupt(INT_PIN_myICM1), icmISR1, FALLING); // Set up a falling interrupt; 
  // myICM2
  pinMode(INT_PIN_myICM2, INPUT_PULLUP);                                   
  attachInterrupt(digitalPinToInterrupt(INT_PIN_myICM2), icmISR2, FALLING); 

  // Initialize the sensors over SPI
  initializeSensor_SPI (myICM1, CS_PIN_myICM1, 1);
  initializeSensor_SPI (myICM2, CS_PIN_myICM2, 2);

  // Configure the DMP for each sensor
  // configureDMP (myICM1);
  // configureDMP (myICM2);

  // Configure sensors for raw data pulls
  configureSensor (myICM1);
  configureSensor (myICM2);

  // Configure interrupts for each sensor
  configureInterrupts (myICM1);
  configureInterrupts (myICM2);

  #ifdef SENSOR3
  
    myICM3.enableDebugging(); // Uncomment this line to enable helpful debug messages on Serial about ICM3
  
    // Setup Interrupts
    // myICM3
    pinMode(INT_PIN_myICM3, INPUT_PULLUP);                                   
    attachInterrupt(digitalPinToInterrupt(INT_PIN_myICM3), icmISR3, FALLING); 
  
    // Initialize the sensors over SPI
    initializeSensor_SPI (myICM3, CS_PIN_myICM3, 3);
  
    // Configure the DMP for each sensor
    // configureDMP (myICM3);
  
    // Configure sensors for raw data pulls
    configureSensor (myICM3);
  
    // Configure interrupts for each sensor
    configureInterrupts (myICM3);

  #endif

  // Set axes: 
  Sensor1_Xaxis.x = 1.0;
  Sensor1_Yaxis.y = 1.0;
  Sensor1_Zaxis.z = 1.0;
  Sensor2_Xaxis.x = 1.0;
  Sensor2_Yaxis.y = 1.0;
  Sensor2_Zaxis.z = 1.0;
  #ifdef SENSOR3
    Sensor3_Xaxis.x = 1.0;
    Sensor3_Yaxis.y = 1.0;
    Sensor3_Zaxis.z = 1.0;
  #endif

  Sensor1_XaxisIP = Sensor1_Xaxis;
  Sensor1_YaxisIP = Sensor1_Yaxis;
  Sensor1_ZaxisIP = Sensor1_Zaxis;
  Sensor2_XaxisIP = Sensor2_Xaxis;
  Sensor2_YaxisIP = Sensor2_Yaxis;
  Sensor2_ZaxisIP = Sensor2_Zaxis;
  #ifdef SENSOR3
    Sensor3_XaxisIP = Sensor3_Xaxis;
    Sensor3_YaxisIP = Sensor3_Yaxis;
    Sensor3_ZaxisIP = Sensor3_Zaxis;
  #endif

  SERIAL_PORT.println();
  SERIAL_PORT.println(F("Configuration complete!"));
  SERIAL_PORT.println();

  int IPcounter = 0;
  int IPcounterMax = 2000;
  bool IPSAVED = false;

  // Pause all setup to calibrate sound level (initiates at max)
  SERIAL_PORT.println(F("Press any key once sound level is calibrated."));
  while (SERIAL_PORT.available()) {
    SERIAL_PORT.read();
    delay(1);
  }
  while (!SERIAL_PORT.available()) {
    delay(1);
  }

  amplitude_update(0);

  // Cue for syncing unit recordings with external OptiTrack recording
  SERIAL_PORT.println();
  SERIAL_PORT.println(F("Begin OPTITRACK reading in 5 seconds ... "));

  for ( int i = 5 ; i >= 0 ; i-- ) {

    if ( i > 0 ) {
      SERIAL_PORT.println();
      SERIAL_PORT.println(i,1);
      delay(1000);
    } else {
      SERIAL_PORT.println();
      double time_now = (double) millis() / 1000;
      SERIAL_PORT.print(time_now,3);
      SERIAL_PORT.println(F(", NOW!"));
    }
    
  }

  delay(500);

  // Pause until user is still
  SERIAL_PORT.println();
  SERIAL_PORT.println(F("Must Capture IP and gyro drift; press any key when still."));
  while (SERIAL_PORT.available()) {
    SERIAL_PORT.read();
    delay(1);
  }
  while (!SERIAL_PORT.available()) {
    delay(1);
  }

  delay(50);

  while (!IPSAVED) { // Grab IP and estimate gyro drift

    if (isrFired1) { // If our isr flag is set then clear the interrupts on the ICM
      isrFired1 = false;
      FetchReadingsSensor (myICM1, Sensor1_Acc, Sensor1_Gyro); //, Sensor1_Mag); // Not pulling Mag data
    }
    if (isrFired2) { // If our isr flag is set then clear the interrupts on the ICM
      isrFired2 = false;
      FetchReadingsSensor (myICM2, Sensor2_Acc, Sensor2_Gyro); //, Sensor2_Mag); 
    }
    #ifdef SENSOR3
      if (isrFired3) { // If our isr flag is set then clear the interrupts on the ICM
        isrFired3 = false;
        FetchReadingsSensor (myICM3, Sensor3_Acc, Sensor3_Gyro); //, Sensor3_Mag);
      }
    #endif

    if (IPcounter == 0) {
      Sensor1_GravityIP = Sensor1_Acc; // Grab the acceleration vector as the gravity vector at IP
      Sensor2_GravityIP = Sensor2_Acc;
      #ifdef SENSOR3
        Sensor3_GravityIP = Sensor3_Acc;
      #endif
      Sensor1_GyroBias = Sensor1_Gyro; // the gyro should be at rest (0,0,0), so whatever readings we get here need to be subtracted as we integrate
      Sensor2_GyroBias = Sensor2_Gyro;
      #ifdef SENSOR3
        Sensor3_GyroBias = Sensor3_Gyro;
      #endif
      SERIAL_PORT.println();
      SERIAL_PORT.println(F("grabbing IP."));
    }
    if (IPcounter > 0) {
      Sensor1_GravityIP = vrunning_avg(Sensor1_GravityIP,Sensor1_Acc,IPcounter); // Now take the vectors a bunch more and average them
      Sensor2_GravityIP = vrunning_avg(Sensor2_GravityIP,Sensor2_Acc,IPcounter);
      #ifdef SENSOR3
        Sensor3_GravityIP = vrunning_avg(Sensor3_GravityIP,Sensor3_Acc,IPcounter);
      #endif
      Sensor1_GyroBias = vrunning_avg(Sensor1_GyroBias,Sensor1_Gyro,IPcounter);
      Sensor2_GyroBias = vrunning_avg(Sensor2_GyroBias,Sensor2_Gyro,IPcounter);
      #ifdef SENSOR3
        Sensor3_GyroBias = vrunning_avg(Sensor3_GyroBias,Sensor3_Gyro,IPcounter);
      #endif
    }
    if (IPcounter == IPcounterMax) {
      IPSAVED = true;
      IPtolerance = IPtolerancefraction * vmag(Sensor1_GravityIP);           
      SERIAL_PORT.println();
      SERIAL_PORT.println(F("IP Saved! Gravity Vectors:"));
      print_vector3ln(Sensor1_GravityIP);
      print_vector3ln(Sensor2_GravityIP);
      #ifdef SENSOR3
        print_vector3ln(Sensor3_GravityIP);
      #endif
      float GMag1 = vmag (Sensor1_GravityIP);
      float GMag2 = vmag (Sensor2_GravityIP);
      #ifdef SENSOR3
        float GMag3 = vmag (Sensor3_GravityIP);
      #endif
      SERIAL_PORT.println();
      SERIAL_PORT.println(F("Gravity Vector Magnitudes:"));
      SERIAL_PORT.println(GMag1,3);
      SERIAL_PORT.println(GMag2,3);
      #ifdef SENSOR3
        SERIAL_PORT.println(GMag3,3);
      #endif
      Sensor1_GyroBias = vscaler_mult(Sensor1_GyroBias,dt);
      Sensor2_GyroBias = vscaler_mult(Sensor2_GyroBias,dt);
      #ifdef SENSOR3
        Sensor3_GyroBias = vscaler_mult(Sensor3_GyroBias,dt);
      #endif
      SERIAL_PORT.println();
      SERIAL_PORT.println(F("Gyro Bias Vectors: (How many degrees will the ROT drift per sample at rest?)"));
      print_vector3ln(Sensor1_GyroBias);
      print_vector3ln(Sensor2_GyroBias);
      #ifdef SENSOR3
        print_vector3ln(Sensor3_GyroBias);
      #endif 
      //Sensor1_GyroBias = vscaler_mult(Sensor1_GyroBias,dt);
      //Sensor2_GyroBias = vscaler_mult(Sensor2_GyroBias,dt);
      //Sensor3_GyroBias = vscaler_mult(Sensor3_GyroBias,dt);
      //SERIAL_PORT.println();
      //SERIAL_PORT.println(F("Gyro Bias / ROT drift Vectors: (How many degrees will the ROT drift per sample at rest?)"));
      //print_vector3(Sensor1_GyroBias);
      //print_vector3(Sensor2_GyroBias);
      //#ifdef SENSOR3
      //  print_vector3ln(Sensor3_GyroBias);
      //#endif 
    }
    
    IPcounter++;

    myICM1.clearInterrupts();
    myICM2.clearInterrupts();
    #ifdef SENSOR3
      myICM3.clearInterrupts();
    #endif

  }

  randomSeed(millis());
  randomsamplenum = random(1,(randomsamplesize -4)); 

  SERIAL_PORT.print(F("Random Sample Number (Model): "));
  SERIAL_PORT.println(randomsamplenum,1);

  SERIAL_PORT.println();
  SERIAL_PORT.println(F("Time, Reach, Sample, amp, pitch, total-error, MOTnum, qax, qay, qaz, qar, qbx, qby, qbz, qbr, high-gyros, total-samples, total-time, rate"));

  // Sanity Check! FIFO should not be turned on, or used, but make sure it's turned off for sure!
  myICM1.enableFIFO(false); // Need to keep FIFO current?
  myICM2.enableFIFO(false);
  #ifdef SENSOR3
    myICM3.enableFIFO(false);
  #endif 
 
}

void loop() {

  // Always fetch readings as soon as interrupts are noticed
  if (isrFired1) { // If our isr flag is set then clear the interrupts on the ICM
    isrFired1 = false;
    FetchReadingsSensor (myICM1, Sensor1_Acc, Sensor1_Gyro); //, Sensor1_Mag); // Not pulling Mag data
  }
  if (isrFired2) { // If our isr flag is set then clear the interrupts on the ICM
    isrFired2 = false;
    FetchReadingsSensor (myICM2, Sensor2_Acc, Sensor2_Gyro); //, Sensor2_Mag); 
  }
  #ifdef SENSOR3
    if (isrFired3) { // If our isr flag is set then clear the interrupts on the ICM
      isrFired3 = false;
      FetchReadingsSensor (myICM3, Sensor3_Acc, Sensor3_Gyro); //, Sensor3_Mag);
    }
  #endif

  //  Run Motion Check (always keep tabs on if we're in motion or still) 
  float vmag_gyro = vmagSqrd(MOTIONCHECKSENSOR_Gyro);
  if (vmag_gyro < motionthreshold_gyro) { //make sure to select a sensor whose movement correlates well with the whole system
    if (jostle_buff_count_saved == jostle_buff_count) { // don't wait around if jostling isn't being detected
      rest_buff_count += 6;
    } else {
      rest_buff_count++;
    }
    if (rest_buff_count > rest_noise_filter) { // if not moving 
      IN_MOTION = false; // Signal we're at rest     
      if (WAITING_FOR_REST) MOTION_JUST_ENDED = true; // signal that motion has just ended
    }
    jostle_buff_count_saved = jostle_buff_count;
  } else { // else the gyro readings indicate movement
    if (vmag_gyro > motionthreshold_gyro_high) { // don't wait around if we're definitely moving
      jostle_buff_count += 11; 
    } else {
      jostle_buff_count++;
    }
    if (jostle_buff_count > jostle_noise_filter) { // If we're moving
      IN_MOTION = true; // Signal we're in motion
      if (WAITING_FOR_INITATION) MOTION_JUST_STARTED = true; // signal that we have initiated movement
    }
  }

  // There are certain routines that should run only once, as soon as the unit decides it's moving
  if (MOTION_JUST_STARTED) {
    Realsamplecount = 0; // Reset the two variables keeping track of the number of samples since motion has initiated
    MOTsamplecount = 0;
    MOTadvance = 1;
    rest_buff_count = 0; // Reset the rest buffer
    jostle_buff_count = 0; // Reset the jostle buffer
    WAITING_FOR_INITATION = false; // By flipping to false, unit is remembering (for next time through loop) that the motion hasn't just started
    MOTION_JUST_STARTED = false; // This status should turn itself off once it's run
    WAITING_FOR_REST = true;
    //PRINTED_RUN_SUMMARY_DATA = false; 
  }

  // There are certain routines that should run only once, as soon as the unit decides its still
  if (MOTION_JUST_ENDED) {
    time2 = micros(); 
    amplitude_update(0); // End to motion should always turn the sound off
    jostle_buff_count = 0; // Reset the jostle buffer
    rest_buff_count = 0; // Reset the rest buffer
    WAITING_FOR_REST = false; // By flipping to false, unit is remembering (for next time through loop) that the motion hasn't just ended
    MOTION_JUST_ENDED = false; // This status should turn itself off once it's run 
    WAITING_FOR_INITATION = true;

    if ( IPfound && Realsamplecount > 250 ) { // don't count short accidental jerks under 250ms

      if ( reachnum > randomsamplesize ) {
        if ( TERMINAL && reachnum <= main_block_count ) {
          Playback_reach(); 
        } else if ( !TERMINAL && reachnum > main_block_count ) {
          Playback_reach(); 
        }
        delay(500);
        Playback_model();
      }
       
      time3 = time2 - time1;
      float samplingtime = (float)( (float)time3 / (float)1000000 );
  
      samplerate = Realsamplecount / samplingtime; 

      for ( int i = 0; i < MaxReadingsToSave; i++ ) {
        if ( i <= (Realsamplecount - 1) ) {

          error1 = vdst4(quat1_saved[i],model_q1[MOTsamplecount_saved[i]]); // max error = 2
          error2 = vdst4(quat2_saved[i],model_q2[MOTsamplecount_saved[i]]); // max error = 2
          #ifdef SENSOR3
            error3 = vdst4(quat3_saved[i],model_q3[MOTsamplecount_saved[i]]); // max error = 2
            error_total = error1 + error2 + erro3; // max of 6
            float to_end = vdst4(quat1_saved[i],model_q1[MOTsamplecount_max]) + 
                           vdst4(quat2_saved[i],model_q2[MOTsamplecount_max]) + 
                           vdst4(quat3_saved[i],model_q3[MOTsamplecount_max]);
          #else 
            error_total = error1 + error2; // max of 4
            float to_end = vdst4(quat1_saved[i],model_q1[MOTsamplecount_max]) + 
                           vdst4(quat2_saved[i],model_q2[MOTsamplecount_max]);
          #endif

          int amp_reconstructed = Son_amp(error_total); 
          int TF_reconstructed = Son_pitch(to_end);

          amp_reconstructed = constrain(amp_reconstructed,0,amp_max);
          TF_reconstructed = constrain(TF_reconstructed,PitchValueMin,PitchValueMax);

          double time_now = (double) millis() / 1000;
          
          if ( i < (Realsamplecount - 1) ) {
            SERIAL_PORT.print(time_now,3);
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.print(reachnum,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(i,1); // print real sample count
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(amp_reconstructed,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(TF_reconstructed,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(error_total,5);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(MOTsamplecount_saved[i],1);
            SERIAL_PORT.print(F(", "));
            print_vector4(quat1_saved[i]);
            SERIAL_PORT.print(F(", "));
            print_vector4ln(quat2_saved[i]);
          }
          if ( i == (Realsamplecount -1) ) {
            SERIAL_PORT.print(time_now,3);
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.print(reachnum,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(i,1); // print real sample count
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(amp_reconstructed,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(TF_reconstructed,1);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(error_total,5);
            SERIAL_PORT.print(F(", ")); 
            SERIAL_PORT.print(MOTsamplecount_saved[i],1);
            SERIAL_PORT.print(F(", "));
            print_vector4 (quat1_saved[i]);
            SERIAL_PORT.print(F(", "));
            print_vector4 (quat2_saved[i]);
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.print(high_gyro_readings,3);
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.print(Realsamplecount); 
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.print(samplingtime,3); 
            SERIAL_PORT.print(F(", "));
            SERIAL_PORT.println(samplerate,1);
          }
        }
      }

      if ( reachnum == randomsamplenum ) {

        MOTsamplecount_max = Realsamplecount - 1; // Gets taken one to far, but okay, since starts at zero. 
        
        for ( int i = 0; i < MaxReadingsToSave; i++ ) {
          if ( i < Realsamplecount ) { // Realsamplecount is taken one too far
            model_q1[i] = quat1_saved[i];
            model_q2[i] = quat2_saved[i]; 
            #ifdef SENSOR3
              model_q3[i] = quat3_saved[i];
            #endif
          }
        }

        #ifdef SENSOR3
        to_end_start = vdst4(quat1_saved[0],model_q1[MOTsamplecount_max]) + 
                       vdst4(quat2_saved[0],model_q2[MOTsamplecount_max]) + 
                       vdst4(quat3_saved[0],model_q3[MOTsamplecount_max]);
        #else
        to_end_start = vdst4(quat1_saved[0],model_q1[MOTsamplecount_max]) + 
                       vdst4(quat2_saved[0],model_q2[MOTsamplecount_max]);
        #endif

        a_amp = to_end_start / ( -1.0 * logf( 1.0 / (float)amp_max ) );
        a_pitch = to_end_start / ( -1.0 * logf( (float)pitchinitial / (float)PitchValueMax ) );

        //float loga = ( -1.0 * logf( 1.0 / (float)amp_max ) );
        //float logp = ( -1.0 * logf( (float)pitchinitial / (float)PitchValueMax ) );

        // SERIAL_PORT.print(F("a_amp: "));
        //     SERIAL_PORT.println(a_amp,5);
        //     SERIAL_PORT.print(F("a_pitch: "));
        //     SERIAL_PORT.println(a_pitch,5);
        //     SERIAL_PORT.print(F("to_end_start: "));
        //     SERIAL_PORT.println(to_end_start,5);
        //     SERIAL_PORT.print(F("loga: "));
        //     SERIAL_PORT.println(loga,5);
        //     SERIAL_PORT.print(F("logp: "));
        //     SERIAL_PORT.println(logp,5);
        
      }

      reachnum++;
      
    }

    IPfound = false; // Signal we're now out of the IP

    high_gyro_readings = 0;
    Sensor1_Rot = zerov;
    Sensor2_Rot = zerov;
    Sensor1_Quat = zerov;
    Sensor2_Quat = zerov;
    Sensor1_Xaxis = Sensor1_XaxisIP;
    Sensor1_Yaxis = Sensor1_YaxisIP;
    Sensor1_Zaxis = Sensor1_ZaxisIP;
    Sensor2_Xaxis = Sensor2_XaxisIP;
    Sensor2_Yaxis = Sensor2_YaxisIP;
    Sensor2_Zaxis = Sensor2_ZaxisIP;
    #ifdef SENSOR3
      Sensor3_Rot = zerov;
      Sensor3_Quat = zerov;
      Sensor3_Xaxis = Sensor3_XaxisIP;
      Sensor3_Yaxis = Sensor3_YaxisIP;
      Sensor3_Zaxis = Sensor3_ZaxisIP;
    #endif

    //delaywithFIFOreset(1); // Not needed since we're not using FIFO
    
  }

  if (!IN_MOTION && !IPfound) {
        
    if ( vdst(Sensor1_GravityIP,Sensor1_Acc) < IPtolerance ){ // We are near the IP
      if ( vdst(Sensor2_GravityIP,Sensor2_Acc) < IPtolerance ) {
        #ifdef SENSOR3
          if ( vdst(Sensor3_GravityIP,Sensor3_Acc) < IPtolerance ) {
        #endif

          IPfound = true;
          UpdateTF(pitchinitial);
          amplitude_update((int)amp_default*0.15);
          rest_buff_count = 0; // Reset the rest buffer so the unit response doesn't bounce around too much
          jostle_buff_count = 0; // Reset the jostle buffer so the unit response doesn't bounce around too much
          
        #ifdef SENSOR3
          }
        #endif
      }
    }

  }

  if ( UPDATE_SOUND_NOW && IN_MOTION && IPfound ) {

    if ( Realsamplecount == 0 ) {
      time1 = micros();
      FIRSTSAMPLE = true;
    }

    compute_motion();

    if ( Realsamplecount < MaxReadingsToSave ) {
      quat1_saved[Realsamplecount] = Sensor1_Quat;
      quat2_saved[Realsamplecount] = Sensor2_Quat;
      #ifdef SENSOR3
        quat3_saved[Realsamplecount] = Sensor3_Quat;
      #endif
    }

    if ( reachnum > randomsamplesize ) {

      // compute error
      error1 = vdst4(Sensor1_Quat,model_q1[MOTsamplecount]); // max error = 2
      error2 = vdst4(Sensor2_Quat,model_q2[MOTsamplecount]); // max error = 2
      #ifdef SENSOR3
        error3 = vdst4(Sensor3_Quat,model_q3[MOTsamplecount]); // max error = 2
        error_total = error1 + error2 + erro3; // max of 6
        float to_end = vdst4(Sensor1_Quat,model_q1[MOTsamplecount_max]) + 
                        vdst4(Sensor2_Quat,model_q2[MOTsamplecount_max]) + 
                        vdst4(Sensor3_Quat,model_q3[MOTsamplecount_max]);
      #else 
        error_total = error1 + error2; // max of 4
        float to_end = vdst4(Sensor1_Quat,model_q1[MOTsamplecount_max]) + 
                        vdst4(Sensor2_Quat,model_q2[MOTsamplecount_max]);
      #endif

      if ( TERMINAL && reachnum > main_block_count ) {
          // compute and update sound 
          amplitude_update(Son_amp(error_total));
          UpdateTF(Son_pitch(to_end));
        } else if ( !TERMINAL && reachnum <= main_block_count ) {
          // compute and update sound 
          amplitude_update(Son_amp(error_total));
          UpdateTF(Son_pitch(to_end));
        }

      // save data for printing at end of motion 
      if ( Realsamplecount < MaxReadingsToSave ) {
        MOTsamplecount_saved[Realsamplecount] = MOTsamplecount;
      }
      
      // determine what MOTsample we should be at, based on current QUAT vectors
      if ( MOTsamplecount < 100 ) { // first 100 MOTsamplecount are a special case; just advance once (initiation too noisey)
        MOTsamplecount++;
      }
      else {
        #ifdef SENSOR3
          float distn1 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount-1]) +
                         vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount-1]) + 
                         vdstSq4(Sensor3_Quat,model_q3[MOTsamplecount-1]);
          float dist0 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount]) +
                        vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount]) + 
                        vdstSq4(Sensor3_Quat,model_q3[MOTsamplecount]);
          float dist1 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount+1]) +
                        vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount+1]) + 
                        vdstSq4(Sensor3_Quat,model_q3[MOTsamplecount+1]);
        #else
          float distn1 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount-1]) +
                         vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount-1]);
          float dist0 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount]) +
                        vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount]);
          float dist1 = vdstSq4(Sensor1_Quat,model_q1[MOTsamplecount+1]) +
                        vdstSq4(Sensor2_Quat,model_q2[MOTsamplecount+1]);
        #endif 
        if ( MOTadvance == 0 ) {
          if ( dist1 < dist0 ) {
            MOTadvance = 2;
          } else {
            MOTadvance = 1;
          }
        } else {
          if ( distn1 < dist0 ) {
            MOTadvance = 0;
          } else if ( dist1 < dist0 ) {
            MOTadvance = 2;
          } else {
            MOTadvance = 1;
          }
        }
        MOTsamplecount = MOTsamplecount + MOTadvance; 
      }
      if (MOTsamplecount > MOTsamplecount_max) MOTsamplecount = MOTsamplecount_max; // don't let us go beyond MOTsamplecount_saved (returns gibberish for error next update)

    }
    
    Realsamplecount++;
        
    UPDATE_SOUND_NOW = false;
    
  }

  myICM1.clearInterrupts();
  myICM2.clearInterrupts();
  #ifdef SENSOR3
    myICM3.clearInterrupts();
  #endif

}

void print_vector3 ( Vectors v ) {

  SERIAL_PORT.print(v.x,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.y,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.z,5); 
  SERIAL_PORT.print(F(", "));
  
}

void print_vector3ln ( Vectors v ) {

  SERIAL_PORT.print(v.x,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.y,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.println(v.z,5); 
  
}

void print_vector4 ( Vectors v ) {

  SERIAL_PORT.print(v.x,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.y,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.z,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.r,5); 
  
}

void print_vector4ln ( Vectors v ) {

  SERIAL_PORT.print(v.x,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.y,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.print(v.z,5); 
  SERIAL_PORT.print(F(", "));
  SERIAL_PORT.println(v.r,5); 
  
}

// void delaywithFIFOreset(int milli) {
//   delay(milli);
//   myICM1.resetFIFO(); // Need to keep FIFO current?
//   myICM2.resetFIFO();
//   #ifdef SENSOR3
//     myICM3.resetFIFO();
//   #endif
// }
