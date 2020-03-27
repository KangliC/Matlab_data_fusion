# Matlab_data_fusion
IMU data analysis and data fusion

File motion1 and static1 content data from real gyro sensor(mpu9150). User should load one of them before applying any analysis function.

Notice that the sensor data in motion1 and static1 are not calibrated (The data just being read out from the sensor.)!

It is easy to do bias calculation of gyro sensor because we have accel sensor data which is 'stable data' 
while the gyro data is 'unstable'. Users can see it from DCM model. But it is not that easy to do accel bias 
calculation. Usually the bias can be find out with manual method. And here we pressent a method to do that online.

One possible method to do online calibration of Accel sensor bias pressented below. 
We do prediction on object velocity and accel sensor bias and config the measurement 
as velocity. Since the system for testing is static we can image the velocity is 0 m/s. 
A linear kalman model can be build and do calculation. We can see the result of velocity of the 
system keeps pretty stable and stop increasing with the new accel bias kalman model. For the 
limitation of lacking extra sensor data like velocity and position, we can only do the model in static. But it works pretty well, we can see from the .jpg result. A better model can 
be built include the position data which is also 0 for our static system.

Apply cof_betaTest(sensorData, 0.001, 0, 0.5, 'alg_name') function to start the analysis, alg_name is the name of the algorithm you want.
