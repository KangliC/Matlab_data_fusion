function [obj, Qret, roll, pitch, heading] = cof_betaTest(data, deltat, zeta, beta, algorithm_pick)

tArray = [];
EulerAngArray_roll = [];
EulerAngArray_pitch = [];
EulerAngArray_heading = [];
accelArray = [];
distArray = [];
velArray = [];

EulerAng.roll = 0;
EulerAng.pitch = 0;
EulerAng.heading = 0;
[m,n] = size(data);
Q.q1 = 1;
Q.q2 = 0;
Q.q3 = 0;
Q.q4 = 0;

%cov(data(1))
dcmobj = DCM_IMU;%'MeasurementVariance', var(data(3,1:n))
dcm9obj = DCM_IMU9;%('DCMVariance', var(data(3,1:n)));
pi_imu6obj = PI_IMU6;
eskf_imu9obj = ESKF_IMU9;
fkf_imu9obj = FKF_IMU9;

tic;
for i = 1:n
    accel = [data(1,i) data(2,i) data(3,i)];
    gyro = [data(4,i) data(5,i) data(6,i)];
    mag = [data(7,i) data(8,i) data(9,i)];
    
    if strcmp(algorithm_pick, 'FKF')
        fkf_imu9obj = fkf_imu9obj.Update_IMU9(gyro,accel,mag,deltat);
        accelArray = [accelArray fkf_imu9obj.acc];
        distArray = [distArray fkf_imu9obj.dist];
        velArray = [velArray fkf_imu9obj.vel];
        EulerAngArray_roll = [EulerAngArray_roll fkf_imu9obj.roll];
        EulerAngArray_pitch = [EulerAngArray_pitch fkf_imu9obj.pitch];
        EulerAngArray_heading = [EulerAngArray_heading fkf_imu9obj.yaw];
    elseif strcmp(algorithm_pick, 'DCM')
        dcmobj = dcmobj.UpdateIMU(gyro,accel,mag,deltat);
        accelArray = [accelArray dcmobj.a];
        distArray = [distArray dcmobj.dist];
        velArray = [velArray dcmobj.vel];
        EulerAngArray_roll = [EulerAngArray_roll dcmobj.roll/pi*180];
        EulerAngArray_pitch = [EulerAngArray_pitch dcmobj.pitch/pi*180];
        EulerAngArray_heading = [EulerAngArray_heading dcmobj.yaw/pi*180];
    elseif strcmp(algorithm_pick, 'ESKF_eCOMPASS')
        eskf_imu9obj = eskf_imu9obj.Update_eCompass(accel, mag);
        EulerAngArray_roll = [EulerAngArray_roll eskf_imu9obj.roll/pi*180];
        EulerAngArray_pitch = [EulerAngArray_pitch eskf_imu9obj.pitch/pi*180];
        EulerAngArray_heading = [EulerAngArray_heading eskf_imu9obj.yaw/pi*180];
    elseif strcmp(algorithm_pick, 'PI_IMU6')
        pi_imu6obj = pi_imu6obj.UpdateIMU6(gyro,accel,deltat);
        accelArray = [accelArray pi_imu6obj.a];
        EulerAngArray_roll = [EulerAngArray_roll pi_imu6obj.roll];
        EulerAngArray_pitch = [EulerAngArray_pitch pi_imu6obj.pitch];
        EulerAngArray_heading = [EulerAngArray_heading pi_imu6obj.yaw];
    elseif strcmp(algorithm_pick, 'ESKF_IMU9')
        eskf_imu9obj = eskf_imu9obj.Update_IMU9(gyro,accel,mag,deltat);
        accelArray = [accelArray eskf_imu9obj.AccGl];
        EulerAngArray_roll = [EulerAngArray_roll eskf_imu9obj.roll];
        EulerAngArray_pitch = [EulerAngArray_pitch eskf_imu9obj.pitch];
        EulerAngArray_heading = [EulerAngArray_heading eskf_imu9obj.yaw];
    elseif strcmp(algorithm_pick, 'DCM9')
        dcm9obj = dcm9obj.UpdateIMU9(gyro,accel,mag,deltat);
        accelArray = [accelArray dcm9obj.a];
        distArray = [distArray dcm9obj.dist];
        velArray = [velArray dcm9obj.vel];
        EulerAngArray_roll = [EulerAngArray_roll dcm9obj.roll/pi*180];
        EulerAngArray_pitch = [EulerAngArray_pitch dcm9obj.pitch/pi*180];
        EulerAngArray_heading = [EulerAngArray_heading dcm9obj.yaw/pi*180];
    elseif strcmp(algorithm_pick, 'ahrs')
        Q = ahrs_filter(Q,data(4,i),data(5,i),data(6,i),data(1,i),data(2,i),data(3,i),data(7,i),data(8,i),data(9,i),deltat,zeta,beta);
        Qx.q0 = Q.q1;
        Qx.q1 = Q.q2;
        Qx.q2 = Q.q3;
        Qx.q3 = Q.q4;
        EulerAng = computeAngles(Qx);
        EulerAngArray_roll = [EulerAngArray_roll EulerAng.roll];
        EulerAngArray_pitch = [EulerAngArray_pitch EulerAng.pitch];
        EulerAngArray_heading = [EulerAngArray_heading EulerAng.heading];
    else
        error('invalid argument');
    end
    tArray = [tArray i];
end
ct = toc;
roll = EulerAng.roll;
pitch = EulerAng.pitch;
heading = EulerAng.heading;
Qret = Q;
obj = eskf_imu9obj;

figure;
subplot(4,1,1);
set(gcf, 'Color', 'White');
hold on; box on;

ylim([-180 180]);
plot(tArray, EulerAngArray_roll(1,:), 'r-', 'LineWidth', 1);
plot(tArray, EulerAngArray_pitch(1,:), 'g-', 'LineWidth', 1);
plot(tArray, EulerAngArray_heading(1,:), 'b-', 'LineWidth', 1);
set(gca,'FontSize',12); 
xlabel('Time (ms)');
ylabel('Euler Ang (deg)');
legend('Roll', 'Pitch', 'Heading');
title(['Consume ' num2str(ct) ' Seconds'],'FontSize',10);
if strcmp(algorithm_pick, 'ahrs')
    title(['Beta = ' num2str(beta), '  Zeta = ' num2str(zeta), '  t = ' num2str(ct)], 'FontSize', 10);
    return
elseif strcmp(algorithm_pick, 'PI_IMU6')
    title(['t = ' num2str(ct)], 'FontSize', 10);
    return
end

subplot(4,1,2);
hold on;box on;
set(gcf, 'Color', 'White');
plot(tArray, accelArray(1,:), 'r-', 'LineWidth', 1);
plot(tArray, accelArray(2,:), 'g-', 'LineWidth', 1);
plot(tArray, accelArray(3,:), 'b-', 'LineWidth', 1);
set(gca,'FontSize',12); 
xlabel('Time (ms)');
ylabel('Acceleration (m/s^2)');
legend('x', 'y', 'z');

subplot(4,1,3);
hold on;box on;
set(gcf, 'Color', 'White');
plot(tArray, velArray(1,:), 'r-', 'LineWidth', 1);
plot(tArray, velArray(2,:), 'g-', 'LineWidth', 1);
plot(tArray, velArray(3,:), 'b-', 'LineWidth', 1);
set(gca,'FontSize',12); 
xlabel('Time (ms)');
ylabel('Velocity (m/s)');
legend('x', 'y', 'z');

subplot(4,1,4);
hold on;box on;
set(gcf, 'Color', 'White');
plot(tArray, distArray(1,:), 'r-', 'LineWidth', 1);
plot(tArray, distArray(2,:), 'g-', 'LineWidth', 1);
plot(tArray, distArray(3,:), 'b-', 'LineWidth', 1);
set(gca,'FontSize',12); 
xlabel('Time (ms)');
ylabel('Distance (m)');
legend('x', 'y', 'z');

return
