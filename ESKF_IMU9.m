classdef ESKF_IMU9 < handle
    
    properties (Access = public)
        g0 = 9.8;	% gravitation around Helsinki, Finland (change according to your area)
        B0 = 71;    % calculated geomagnetic sphere
        Qvy = 0.0006;  % gyro sensor noise variance units (rad/s)^2
        Qwb = 1e-4; % gyro offset random walk units (rad/s)^2
        Qvb = 8;	% magnetometer sensor noise variance units uT^2 defining minimum deiation from geomagnetic sphere.
        Qvg = 0.002;     % accelerometer sensor noise variance units g^2 defining minimum deviation from 1g sphere.
        Gyro_Offset_MIN = -0.2; % minimum permissible power on gyro offsets (rad/s)
        Gyro_Offset_MAX = 0.2; % maximum permissible power on gyro offsets (rad/s)
        
        Rot = eye(3,3);
        Q = struct('q0',1,'q1',0,'q2',0,'q3',0);    %[1,0,0,0];
        yaw = 0;	% Yaw angle around z axis (in ZYX convention)
        pitch = 0;	% Pitch angle around y axis
        roll = 0;	% Roll angle around x axis distortion

        B_theta = 0;            % in radian
        berr = [0; 0; 0];
        qgerr = [0;0;0];
        qmerr = [0;0;0];
        gyro_offset = [0 0 0]';	% gyro bias
        AccGl = [0 0 0]';       %linear acceleration in global frame
    end
    properties (Access = private)
        gyro_sum = 0;
        gyro_cnt = 0;
        mag_prev = [0 0 0]';
        mag_first_lock = 0;
    end
    
    methods (Access = private)
        % Bdelta is in radian
        function [R, Bdelta, B, G] = eCompass(obj, accel, magnat)
            if (size(accel,1) == 3), a = accel;
            else a = accel';
            end
            if (size(magnat,1) == 3), m = magnat;
            else m = magnat';
            end
            
            Bdelta = 0;
            Rz = a;
            Rx = m;
            Ry = cross(Rz,Rx);
            Rx = cross(Ry,Rz);
            
            mod = [norm(Rx) norm(Ry) norm(Rz)];
            if ~(mod(1)==0 || mod(2)==0 || mod(3)==0)
                R = [Rx/mod(1) Ry/mod(2) Rz/mod(3)];
            else
                R = eye(3,3);
                return;
            end
            
            B = norm(m);
            G = mod(3);
            gdotb = dot(a,m);
            if ~(B == 0 || G == 0)
               Bdelta = asin(gdotb/(B*G)); 
            end
        end
        
        function [Q] = qfromGyro(obj, gyro, dt)
            deltaRad = dt*norm(gyro);
            sinhalf = sin(0.5*deltaRad);
            if deltaRad ~= 0
                tmp = sinhalf * dt / deltaRad;
                Q.q1 = tmp*gyro(1);
                Q.q2 = tmp*gyro(2);
                Q.q3 = tmp*gyro(3);
            else
                Q.q1 = 0;
                Q.q2 = 0;
                Q.q3 = 0;
            end
            tmpq0 = Q.q1^2+Q.q2^2+Q.q3^2;
            if tmpq0 <= 1
                Q.q0 = sqrt(1-tmpq0);
            else
                Q.q0 = 0;
            end
        end
        function [R] = RfromQ(obj, q)
            R(1,1) = 2*q.q0^2 + 2*q.q1^2 - 1;
            R(1,2) = 2*q.q1*q.q2 + 2*q.q0*q.q3;
            R(1,3) = 2*(q.q1*q.q3 - q.q0*q.q2);
            R(2,1) = 2*(q.q1*q.q2 - q.q0*q.q3);
            R(2,2) = 2*(q.q0^2 + q.q2^2) - 1;
            R(2,3) = 2*(q.q2*q.q3 + q.q0*q.q1);
            R(3,1) = 2*(q.q1*q.q3 + q.q0*q.q2);
            R(3,2) = 2*(q.q2*q.q3 - q.q0*q.q1);
            R(3,3) = 2*(q.q0^2 + q.q3^2) - 1;
        end
        function [Q] = QfromR(obj, r)
            q0sq = (1 + r(1,1) + r(2,2) + r(3,3))*0.25;
            %q0 is always non-negative
            Q.q0 = sqrt(abs(q0sq));
            if Q.q0 > 1e-4
                tmp = 0.25/Q.q0;
                Q.q1 = (r(2,3) - r(3,2))*tmp;
                Q.q2 = (r(3,1) - r(1,3))*tmp;
                Q.q3 = (r(1,2) - r(2,1))*tmp;
            else
                Q.q1 = sqrt(abs(0.5 + 0.5*r(1,1) - q0sq));
                Q.q2 = sqrt(abs(0.5 + 0.5*r(2,2) - q0sq));
                Q.q3 = sqrt(abs(0.5 + 0.5*r(3,3) - q0sq));
                if r(2,3) < r(3,2) 
                    Q.q1 = -Q.q1; 
                end
                if r(3,1) < r(1,3) 
                    Q.q2 = -Q.q2; 
                end
                if r(1,2) < r(2,1) 
                    Q.q3 = -Q.q3; 
                end
            end
            qnorm = norm([Q.q0 Q.q1 Q.q2 Q.q3]);
            Q.q0 = Q.q0/qnorm;
            Q.q1 = Q.q1/qnorm;
            Q.q2 = Q.q2/qnorm;
            Q.q3 = Q.q3/qnorm;
        end
        
        % v = Q* x u x Q
        % q = 1/sqrt(2) * {sqrt(1 + u.v) - u x v / sqrt(1 + u.v)}
        function [Q] = QrotU2V(obj,u,v)
            uxv = cross(u,v);
            squdotv = sqrt(abs(1+dot(u,v)));
            Q.q0 = squdotv/sqrt(2);
            if squdotv ~= 0
                tmp = 1/sqrt(2)/squdotv;
                Q.q1 = -uxv(1)*tmp;
                Q.q2 = -uxv(2)*tmp;
                Q.q3 = -uxv(3)*tmp;
            else
                y_z = u(2) - u(3);
                z_x = u(3) - u(1);
                x_y = u(1) - u(2);
                tmp = sqrt(y_z^2 + z_x^2 + x_y^2);
                if tmp ~= 0
                    Q.q1 = y_z/tmp;
                    Q.q2 = z_x/tmp;
                    Q.q3 = x_y/tmp;
                else
                    Q.q1 = 1/sqrt(2);
                    Q.q2 = -1/sqrt(2);
                    Q.q3 = 0;
                end
            end
        end
    end
    
    methods (Access = public)
        function obj = ESKF_IMU9(varargin)
            %updateP = true;
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'Gravity'), obj.g0 = varargin{i+1};
                %elseif  strcmp(varargin{i}, 'Covariance'), obj.P = varargin{i+1}; updateP = false;
                else
                    error('Invalid argument');
                end
            end
            
        end
        
        function obj = Update_eCompass(obj, accel, magnat)
            if norm(magnat) ~= 0
                [R, Bdelta, B, G] = obj.eCompass(accel, magnat);
                obj.yaw = atan2(R(1,2),R(1,1));%atan2(Ry(1),Rx(1));
                obj.roll = atan2(R(2,3),R(3,3));%atan2(Rz(2),Rz(3));
                obj.pitch = asin(-R(1,3));%asin(-Rz(1));
                obj.B_theta = Bdelta;
            end
        end
        
        function obj = Update_IMU9(obj, gyro, accel, magnat, period)
            % adapt data
            if (size(gyro,1) ~= 3)
                gyro = gyro';
            end
            if (size(accel,1) ~= 3)
                accel = accel';
            end
            if (size(magnat,1) ~= 3)
                magnat = magnat';
            end
            
            if magnat'*magnat == 0
                if obj.mag_prev'*obj.mag_prev == 0
                    return;
                else
                    magnat = obj.mag_prev;
                end
            else
                obj.mag_prev = magnat;
            end
            qPri = obj.Q;
            % eleminate bias from gyro data
            gyro = gyro - obj.gyro_offset;
            [Q] = obj.qfromGyro(gyro, period);
            qtmp.q0 = qPri.q0*Q.q0 - qPri.q1*Q.q1 - qPri.q2*Q.q2 - qPri.q3*Q.q3;
            qtmp.q1 = qPri.q0*Q.q1 + qPri.q1*Q.q0 + qPri.q2*Q.q3 - qPri.q3*Q.q2;
            qtmp.q2 = qPri.q0*Q.q2 - qPri.q1*Q.q3 + qPri.q2*Q.q0 + qPri.q3*Q.q1;
            qtmp.q3 = qPri.q0*Q.q3 + qPri.q1*Q.q2 - qPri.q2*Q.q1 + qPri.q3*Q.q0;
            qPri = qtmp;        %priori q with new gyro data
            [RMi] = obj.RfromQ(qPri);%rotation matrix from priori q
            [R6DOF, Bdelta, B, G] = obj.eCompass(accel, magnat);%calculate R6dof from accel and magnatmeter
            % calibrated mag data available, and first time lock the
            % initial attitude to 6dof.
            if obj.mag_first_lock == 0
                obj.B_theta = Bdelta;
                RMi = R6DOF;
                [obj.Q] = obj.QfromR(R6DOF);%get q6dof from R6dof
                obj.mag_first_lock = 1;
            end
            tmp = G - 1;
            QvgQa = 3*tmp*tmp;
            if QvgQa < obj.Qvg
                QvgQa = obj.Qvg;
            end
            tmp = B - obj.B0;
            QvbQd = 3*tmp*tmp;
            if QvbQd < obj.Qvb
                QvbQd = obj.Qvb;
            end
            atmp = [R6DOF(1,3); R6DOF(2,3); R6DOF(3,3)];%geo accel from 6dof
            % priori estimate gravity(sensor frame)
            gPri = [RMi(1,3); RMi(2,3); RMi(3,3)]; %geo accel from gyro
            % vector err between 6dof and gyro in geo accel direction.
            [qtmp] = obj.QrotU2V(atmp,gPri);
            Z_err = [0 0 0 0 0 0]';
            Z_err(1) = qtmp.q1;
            Z_err(2) = qtmp.q2;
            Z_err(3) = qtmp.q3;
            %geo magnetor from 6dof
            atmp = [R6DOF(1,1)*cos(Bdelta)+R6DOF(1,3)*sin(Bdelta); R6DOF(2,1)*cos(Bdelta)+R6DOF(2,3)*sin(Bdelta); R6DOF(3,1)*cos(Bdelta)+R6DOF(3,3)*sin(Bdelta)];
            %atmp = [R6DOF(1,1)*cos(obj.B_theta)+R6DOF(1,3)*sin(obj.B_theta); R6DOF(2,1)*cos(obj.B_theta)+R6DOF(2,3)*sin(obj.B_theta); R6DOF(3,1)*cos(obj.B_theta)+R6DOF(3,3)*sin(obj.B_theta)];
            % priori estimate geo magnetor from gyro(sensor frame)
            bPri = [RMi(1,1)*cos(obj.B_theta)+RMi(1,3)*sin(obj.B_theta); RMi(2,1)*cos(obj.B_theta)+RMi(2,3)*sin(obj.B_theta); RMi(3,1)*cos(obj.B_theta)+RMi(3,3)*sin(obj.B_theta)];
            %bPri = [RMi(1,1)*cos(Bdelta)+RMi(1,3)*sin(Bdelta); RMi(2,1)*cos(Bdelta)+RMi(2,3)*sin(Bdelta); RMi(3,1)*cos(Bdelta)+RMi(3,3)*sin(Bdelta)];
            % vector err between 6dof and gyro in geo magnetor direction.
            [qtmp] = obj.QrotU2V(atmp,bPri);
            Z_err(4) = qtmp.q1;
            Z_err(5) = qtmp.q2;
            Z_err(6) = qtmp.q3;
            % Qw,k = E[wk wk'] 
            Qw = zeros(9,9);
            Qw(7,7) = obj.berr(1)^2 + obj.Qwb/3;
            Qw(8,8) = obj.berr(2)^2 + obj.Qwb/3;
            Qw(9,9) = obj.berr(3)^2 + obj.Qwb/3;
            Qvy_Qwb_OVER3 = obj.Qvy/3 + obj.Qwb/3;
            Qw(1,1) = obj.qgerr(1)^2 + 0.25 * period^2 * (Qw(7,7) + Qvy_Qwb_OVER3);
            Qw(2,2) = obj.qgerr(2)^2 + 0.25 * period^2 * (Qw(8,8) + Qvy_Qwb_OVER3);
            Qw(3,3) = obj.qgerr(3)^2 + 0.25 * period^2 * (Qw(9,9) + Qvy_Qwb_OVER3);
            Qw(4,4) = obj.qmerr(1)^2 + 0.25 * period^2 * (Qw(7,7) + Qvy_Qwb_OVER3);
            Qw(5,5) = obj.qmerr(2)^2 + 0.25 * period^2 * (Qw(8,8) + Qvy_Qwb_OVER3);
            Qw(6,6) = obj.qmerr(3)^2 + 0.25 * period^2 * (Qw(9,9) + Qvy_Qwb_OVER3);
            Qw(1,7) = (obj.qgerr(1)*obj.berr(1) - period/2 * (Qw(7,7)));
            Qw(2,8) = (obj.qgerr(2)*obj.berr(2) - period/2 * (Qw(8,8)));
            Qw(3,9) = (obj.qgerr(3)*obj.berr(3) - period/2 * (Qw(9,9)));
            Qw(4,7) = (obj.qmerr(1)*obj.berr(1) - period/2 * (Qw(7,7)));
            Qw(5,8) = (obj.qmerr(2)*obj.berr(2) - period/2 * (Qw(8,8)));
            Qw(6,9) = (obj.qmerr(3)*obj.berr(3) - period/2 * (Qw(9,9)));
            Qw(7,1) = Qw(1,7);
            Qw(8,2) = Qw(2,8);
            Qw(9,3) = Qw(3,9);
            Qw(7,4) = Qw(4,7);
            Qw(8,5) = Qw(5,8);
            Qw(9,6) = Qw(6,9);
            %Qv
            Qv = zeros(6,6);
            Qv(1,1) = (QvgQa + (obj.Qwb + obj.Qvy) * period^2) / 12;
            Qv(2,2) = Qv(1,1);
            Qv(3,3) = Qv(1,1);
            Qv(4,4) = (QvbQd / (obj.B0^2) + (obj.Qwb + obj.Qvy) * period^2) / 12;
            Qv(5,5) = Qv(4,4);
            Qv(6,6) = Qv(4,4);
            %Ck
            Ck = [eye(3) zeros(3,3) -0.5*period*eye(3);
                  zeros(3,3) eye(3) -0.5*period*eye(3)];
            %P- = Qwk
            % Kalman gain
            QwCkT = Qw*Ck';
            Kk = QwCkT / (Ck*QwCkT + Qv);
            X_err = Kk*Z_err;
            % set qgerr, qmerr, berr
            obj.qgerr = X_err(1:3);
            obj.qmerr = X_err(4:6);
            obj.berr = X_err(7:9);
            qtmp.q1 = -obj.qgerr(1);
            qtmp.q2 = -obj.qgerr(2);
            qtmp.q3 = -obj.qgerr(3);
            qtmp.q0 = sqrt(abs(1-qtmp.q1^2-qtmp.q2^2-qtmp.q3^2));
            [rtmp] = obj.RfromQ(qtmp);
            % posteriori gravity(sensor frame)
            gPos = rtmp*gPri;
            qtmp.q1 = -obj.qmerr(1);
            qtmp.q2 = -obj.qmerr(2);
            qtmp.q3 = -obj.qmerr(3);
            qtmp.q0 = sqrt(abs(1-qtmp.q1^2-qtmp.q2^2-qtmp.q3^2));
            [rtmp] = obj.RfromQ(qtmp);
            % posteriori geomagnatic(sensor frame)
            bPos = rtmp*bPri;
            %calculate R6dof from posteriori gravity and geomagnatmeter
            [obj.Rot, obj.B_theta, ~] = obj.eCompass(gPos, bPos);
            [obj.Q] = obj.QfromR(obj.Rot);
            offset_change_max = sqrt(obj.Qwb)*period; 
            for i = 1:3
                if obj.berr(i) > offset_change_max
                    obj.gyro_offset(i) = obj.gyro_offset(i) - offset_change_max;
                elseif obj.berr(i) < -offset_change_max
                    obj.gyro_offset(i) = obj.gyro_offset(i) + offset_change_max;
                else
                    obj.gyro_offset(i) = obj.gyro_offset(i) - obj.berr(i);
                end
                %obj.gyro_offset(i) = obj.gyro_offset(i) - obj.berr(i);
                
                if obj.gyro_offset(i) > obj.Gyro_Offset_MAX
                    obj.gyro_offset(i) = obj.Gyro_Offset_MAX;
                elseif obj.gyro_offset(i) < obj.Gyro_Offset_MIN
                    obj.gyro_offset(i) = obj.Gyro_Offset_MIN;
                end
            end
            accgl = obj.Rot'*accel;
            obj.AccGl = obj.g0*([0; 0; 1] - accgl);
            
            Q.q1 = obj.Q.q0;
            Q.q2 = obj.Q.q1;
            Q.q3 = obj.Q.q2;
            Q.q4 = obj.Q.q3;
            euler = computeAngles(Q);
            obj.yaw = euler.heading;
            obj.roll = euler.roll;
            obj.pitch = euler.pitch;
        end
    end
end
