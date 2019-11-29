classdef FKF_IMU9 < handle
        %% Public properties
    properties (Access = public)
        g0 = 9.8;                   % gravitation (change according to your area)
        gt = [0 0 0]';            % estimated non-gravitational accelerations on earth frame
        acc = zeros(3,1);           % gt in m/s^2
        q_gyro_bias2 = 0.001^2;     % very small number to make bias change slowly
        r_acc2 = 0.1^2;             % variance of calibrated accelerometer (g-component)
        r_a2 = 0.2;                % large variance for some unknown acceleration (acc = a + g)
        r_mag2 = 0.5^2;             % magnetmeter variance
        q_gyro_bias2_init = 0.3^2;  % initial variance of bias states (for bias estimator)
        
        vel = zeros(3,1);           % estimated velocity on earth frame
        dist = zeros(3,1);          % estimated distance from origin on earth frame;
        yaw = 0;                    % Yaw angle around z axis (in ZYX convention)
        pitch = 0;                  % Pitch angle around y axis
        roll = 0;                   % Roll angle around x axis
        P = [];                     % estimate covariance (these are initialized in constructor below)
        %H = [];                    % observation model (static), here is 1
        %Qv = [];                    % measure noise covariance 
        %Qw = [];                    % process noise covariance
        cov_gyro = [];
        cov_acc_mag = [];
        qt = [1 0 0 0]';            % quaternion
        
        mag_sensor = [0 0 0]';
        inited = false;
        cnt = 0;
        gg = [0 0 0]';
        gbias = [0 0 0.01468]';
        mnsum = [0 0];
    end
    
    methods (Access = public)
        
        function obj = FKF_IMU9(varargin)
            obj.cov_gyro = obj.q_gyro_bias2_init*eye(3);
            obj.cov_acc_mag = [obj.r_acc2*eye(3) zeros(3,3);
                                zeros(3,3) obj.r_mag2*eye(3)];
            obj.P = eye(4);
        end
        
        function obj = Update_IMU9(obj, Gyro, Accelerometer, Magnetometer, dt)
            %Data adapt%            
            % control input (angular velocities from gyroscopes)
            if (size(Gyro,1) == 3), u = Gyro;
            else u = Gyro';
            end
            % measurements/observations (acceleromeres)
            if (size(Accelerometer,1) == 3), a = Accelerometer;
            else a = Accelerometer';
            end            
            % measurements/observations (magnetometer)
            if (size(Magnetometer,1) == 3), mag_tmp = Magnetometer;
            else mag_tmp = Magnetometer';
            end
            % mag data fill in case the mag data is not available
            if mag_tmp'*mag_tmp == 0
                if obj.mag_sensor'*obj.mag_sensor == 0
                    return;
                else
                    mag_tmp = obj.mag_sensor;
                end
            else
                obj.mag_sensor = mag_tmp;
            end
            accel = a;
            mag = mag_tmp;
            mag = mag/norm(mag);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            q0 = obj.qt(1);q1 = obj.qt(2);q2 = obj.qt(3);q3 = obj.qt(4);
            qp = [-q1 -q2 -q3;
                  q0  -q3 q2;
                  q3  q0  -q1;
                  -q2 q1  q0];
%             qp = [q1 q2 q3;
%                   -q0  -q3 q2;
%                   q3  -q0  -q1;
%                   -q2 q1  -q0];
              
            Omega = [0   -u(1) -u(2)  -u(3);
                     u(1) 0     u(3)  -u(2);
                     u(2) -u(3) 0     u(1);
                     u(3) u(2)  -u(1) 0];
%             Omega = [0   u(1) u(2)  u(3);
%                     -u(1) 0   u(3)  -u(2);
%                     -u(2) -u(3) 0     u(1);
%                     -u(3) u(2)  -u(1) 0];
                 
            qt_ = (eye(4) + dt/2*Omega)*obj.qt;  % quaternion prediction
            qt_ = qt_./norm(qt_);
            
            if obj.inited
                q_0 = qt_(1);q_1 = qt_(2);q_2 = qt_(3);q_3 = qt_(4);
                %q_0 = obj.qt(1);q_1 = obj.qt(2);q_2 = obj.qt(3);q_3 = obj.qt(4);
                P1 = [q_0 q_1 -q_2 -q_3; -q_3 q_2 q_1 -q_0; q_2 q_3 q_0 q_1];
                P2 = [q_3 q_2 q_1 q_0; q_0 -q_1 q_2 -q_3; -q_1 -q_0 q_3 q_2];
                P3 = [-q_2 q_3 -q_0 q_1; q_1 q_0 q_3 q_2; q_0 -q_1 -q_2 q_3];
            
                DCM = [P1*qt_,P2*qt_,P3*qt_];
                dg_b = DCM*obj.gt; % pure accel in predicted body frame
                %a = a - dg_b;% real accel - pure accel = accel caused by gravity
                val = P3*qt_;
                a_predict = a - val;
                a_len = sqrt(a_predict'*a_predict);
                obj.cov_acc_mag = [(a_len*obj.r_a2 + obj.r_acc2)*eye(3) zeros(3,3);
                                    zeros(3,3) obj.r_mag2*eye(3)];              
            end

            a = a./norm(a);
            ax = a(1);ay = a(2);az = a(3);
            mx = mag(1);my = mag(2);mz = mag(3);
            % mD magnet down, mN magnet north
            mD = -0.8047;%dot(a,mag);%
            mN = 0.5937;%sqrt(1-mD^2);%
            obj.mnsum = obj.mnsum + [mD mN];
            
            Wm = [mN*mag(1)+mD*mag(3) mD*mag(2) mN*mag(3)-mD*mag(1) -mN*mag(2);
                mD*mag(2) mN*mag(1)-mD*mag(3) mN*mag(2) mN*mag(3)+mD*mag(1);
                mN*mag(3)-mD*mag(1) mN*mag(2) -mN*mag(1)-mD*mag(3) mD*mag(2);
                -mN*mag(2) mN*mag(3)+mD*mag(1) mD*mag(2) -mN*mag(1)+mD*mag(3)];
            Wa = [az ay -ax 0;
                  ay -az 0 ax;
                 -ax 0 -az ay;
                   0 ax ay az];
            
            qm = (Wa+eye(4))*(Wm+eye(4))*obj.qt;%/4;  % quaternion from measurement
            qm = qm./norm(qm);
            
            if obj.inited
                Qw = (dt^2)/4*qp*obj.cov_gyro*qp'; % cov for process

                Jacob = zeros(4,6);
                Jacob(1,1)= -q0*(mN*mz-mD*mx)-q1*mN*my-q2*(1-mN*mx-mD*mz)-q3*mD*my;%-q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
                Jacob(1,2)= q0*mD*my+q1*(mN*mx-mD*mz+1)+q2*mN*my+q3*(mN*mz+mD*mx);%q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
                Jacob(1,3)= q0*(mN*mx+mD*mz+1)+q1*mD*my+q2*(mN*mz-mD*mx)-q3*mN*my;%q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
                Jacob(1,4)= (ax*mD + mN + az*mN)*q0 + ay*mN*q1 + (-((1 + az)*mD) + ax*mN)*q2 + ay*mD*q3;
                Jacob(1,5)= ay*mD*q0 + (mD + az*mD - ax*mN)*q1 + ay*mN*q2 - (ax*mD + mN + az*mN)*q3;
                Jacob(1,6)= mD*(q0 + az*q0 - ay*q1 + ax*q2) + mN*(-(ax*q0) + q2 + az*q2 + ay*q3);

                Jacob(2,1)= -q0*mN*my+q1*(mN*mz+mD*mx)+q2*mD*my+q3*(1-mN*mx+mD*mz);%q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
                Jacob(2,2)= q0*(mN*mx+mD*mz+1)+q1*mD*my+q2*(mN*mz-mD*mx)-q3*mN*my;%q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
                Jacob(2,3)= -q0*mD*my-q1*(mN*mx-mD*mz+1)-q2*mN*my-q3*(mN*mz+mD*mx);%-((1 + mN*mx)*q1) - mD*(my*q0 - mz*q1 + mx*q3) - mN*(my*q2 + mz*q3);
                Jacob(2,4)= ay*mN*q0+(1-az)*mN*q1+ax*mD*q1-ay*mD*q2+(1-az)*mD*q3-ax*mN*q3;%ay*(mN*q0 - mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(mD*q1 - mN*q3);
                Jacob(2,5)= mD*(q0 - az*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + (-1 + az)*q2 + ay*q3);
                Jacob(2,6)= ay*(mD*q0 + mN*q2) + mD*((-1 + az)*q1 + ax*q3) + mN*(ax*q1 + q3 - az*q3);

                Jacob(3,1)= -q0*(mN*mx+mD*mz+1)-q1*mD*my-q2*(mN*mz-mD*mx)+q3*mN*my;%-((1 + mN*mx + mD*mz)*q0) - mD*my*q1 + mD*mx*q2 - mN*mz*q2 + mN*my*q3;
                Jacob(3,2)= -q0*mN*my+q1*(mN*mz+mD*mx)+q2*mD*my+q3*(1-mN*mx+mD*mz);%q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
                Jacob(3,3)= -q0*(mN*mz-mD*mx)-q1*mN*my-q2*(1-mN*mx-mD*mz)-q3*mD*my;%-q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
                Jacob(3,4)= mD*((-1 + az)*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + q2 - az*q2 + ay*q3);
                Jacob(3,5)= ay*(-(mN*q0) + mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);
                Jacob(3,6)= mN*(q0 - az*q0 + ay*q1) - ax*(mD*q0 + mN*q2) + mD*((-1 + az)*q2 + ay*q3);

                Jacob(4,1)= q0*mD*my+q1*(mN*mx-mD*mz+1)+q2*mN*my+q3*(mN*mz+mD*mx);%q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
                Jacob(4,2)= q0*(mN*mz-mD*mx)+q1*mN*my+q2*(1-mN*mx-mD*mz)+q3*mD*my;%q2 + mN*(mz*q0 + my*q1 - mx*q2) - mD*(mx*q0 + mz*q2 - my*q3);
                Jacob(4,3)= -q0*mN*my+q1*(mN*mz+mD*mx)+q2*mD*my+q3*(1-mN*mx+mD*mz);%q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
                Jacob(4,4)= -(ay*(mD*q0 + mN*q2)) + ax*(mN*q1 + mD*q3) + (1 + az)*(mD*q1 - mN*q3);
                Jacob(4,5)= (1 + az)*(-(mN*q0) + mD*q2) + ax*(mD*q0 + mN*q2) + ay*(mN*q1 + mD*q3);
                Jacob(4,6)= ay*(mN*q0 - mD*q2) + (1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);
                Jacob = 0.25*Jacob;

                Qv = Jacob*obj.cov_acc_mag*Jacob'; % cov for measure
                P0 = (eye(4)+dt/2*Omega)*obj.P*(eye(4)+dt/2*Omega)' + Qw; % covariance prediction
                K = P0/(P0+Qv);
                qt = qt_+K*(qm - qt_);
                obj.qt = qt./norm(qt);
                obj.P = (eye(4) - K)*P0;%*(eye(4) - K)' + K*Qv*K';
            else
                obj.cnt = obj.cnt + 1;
                if obj.cnt > 800
                    obj.inited = true;
                end
                obj.qt = qm;
            end
            
            q_0 = obj.qt(1);q_1 = obj.qt(2);q_2 = obj.qt(3);q_3 = obj.qt(4);
            P1 = [q_0 q_1 -q_2 -q_3; -q_3 q_2 q_1 -q_0; q_2 q_3 q_0 q_1];
            P2 = [q_3 q_2 q_1 q_0; q_0 -q_1 q_2 -q_3; -q_1 -q_0 q_3 q_2];
            P3 = [-q_2 q_3 -q_0 q_1; q_1 q_0 q_3 q_2; q_0 -q_1 -q_2 q_3];
            DCM = [P1*obj.qt,P2*obj.qt,P3*obj.qt];
            obj.gt = DCM'*accel - [0 0 1]';% cal global pure accel eliminate gravity
            %obj.gg = obj.gg + obj.gt;
            obj.gt = obj.gt - obj.gbias;
            
            % acc, vel, dist
            obj.acc = obj.gt*obj.g0;
            obj.vel = obj.vel + obj.acc*dt;
            obj.dist = obj.dist + obj.vel*dt + 0.5*obj.acc*dt^2;
            
            Q.q0 = obj.qt(1);
            Q.q1 = obj.qt(2);
            Q.q2 = obj.qt(3);
            Q.q3 = obj.qt(4);
            euler = computeAngles(Q);
            obj.yaw = euler.heading;
            obj.roll = euler.roll;
            obj.pitch = euler.pitch;
        end
    end
end