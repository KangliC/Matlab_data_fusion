classdef UKF_DCM_IMU < handle
% UKF_DCM_IMU Implementation 
%
%
%   Date          Author           
%   03/09/2020    Kangli Chu     

    %% Public properties
    properties (Access = public)
        g0 = 9.8;                   % gravitation (change according to your area)
        
        mag_bias = [0 0 0]';
        state = [0 0 1 0 0 0]';     % States are lowest row of rotation matrix and gyroscope x y and z biases
                                    % (C_31, C_32, C_33, w_b1, w_b2, w_b3)
        q_dcm2 = 0.0003;             % estimated variance of dcm states (gyro variance per second)
        q_gyro_bias2 = 0.001^2;    % very small number to make bias change slowly
        r_acc2 = 0.5^2;             % variance of calibrated accelerometer (g-component)
        r_a2 = 10^2;                % large variance for some unknown acceleration (acc = a + g)
        q_dcm2_init = 1^2;          % initial variance of dcm states (for attitude estimation)
        q_gyro_bias2_init = 0.1^2;  % initial variance of bias states (for bias estimator)
        a = zeros(3,1);             % estimated non-gravitational accelerations on earth frame
        vel = zeros(3,1);           % estimated velocity on earth frame
        dist = zeros(3,1);          % estimated distance from origin on earth frame;
        yaw = 0;                    % Yaw angle around z axis (in ZYX convention)
        pitch = 0;                  % Pitch angle around y axis
        roll = 0;                   % Roll angle around x axis
        P = [];                     % estimate covariance (these are initialized in constructor below)
        H = [];                     % observation model (static)
        Q = [];                     % proces noise covariance (static part)

        Dcm = [1 0 0; 0 1 0; 0 0 1];
        Qt = [1 0 0 0]';
        first_row = [1 0 0]';       % first row of of the rotation matrix (for yaw angle estimate)
        
        mag_initial = false;
        mag_vector = [1 0 0]';
        mag_err_int = [0;0;0]; 
        
        cnt = 0;
        gbias = [0 0 0.01468]';
        err_b = 0;
        Ki = 0.001;
        Kp = 0.01;
        
        % UKF parameter
        L = 6;%numel(x);          %numer of states  
        m = 3;%numel(z);          %numer of measurements  
        alpha = 1e0;              %default, tunable  
        ki = 7;                   %default, tunable  
        beta = 2;                 %best for gaussian, tunable 
    end

    %% Public methods
    methods (Access = public)
        function obj = UKF_DCM_IMU(varargin)
            updateP = true;
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'Gravity'), obj.g0 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'State'), obj.state = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Covariance'), obj.P = varargin{i+1}; updateP = false;
                elseif  strcmp(varargin{i}, 'DCMVariance'), obj.q_dcm2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'BiasVariance'), obj.q_gyro_bias2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'InitialDCMVariance'), obj.q_dcm2_init = varargin{i+1};
                elseif  strcmp(varargin{i}, 'InitialBiasVariance'), obj.q_gyro_bias2_init = varargin{i+1};                    
                elseif  strcmp(varargin{i}, 'MeasurementVariance'), obj.r_acc2 = varargin{i+1};
                elseif  strcmp(varargin{i}, 'MeasurementVarianceVariableGain'), obj.r_a2 = varargin{i+1};
                else
                    error('Invalid argument');
                end
            end
            
            if (updateP), obj.P = [obj.q_dcm2_init*eye(3), zeros(3,3); zeros(3,3), obj.q_gyro_bias2_init*eye(3)]; end
            obj.H = [eye(3)*obj.g0, zeros(3,3)];
            obj.Q = [obj.q_dcm2*eye(3), zeros(3,3); zeros(3,3) obj.q_gyro_bias2*eye(3)];
        end
        function [Q] = QfromDCM(obj, dcm)
            r = dcm';
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
        function [dcm] = DCMfromQ(obj, q)
            R(1,1) = 2*q.q0^2 + 2*q.q1^2 - 1;
            R(1,2) = 2*q.q1*q.q2 + 2*q.q0*q.q3;
            R(1,3) = 2*(q.q1*q.q3 - q.q0*q.q2);
            R(2,1) = 2*(q.q1*q.q2 - q.q0*q.q3);
            R(2,2) = 2*(q.q0^2 + q.q2^2) - 1;
            R(2,3) = 2*(q.q2*q.q3 + q.q0*q.q1);
            R(3,1) = 2*(q.q1*q.q3 + q.q0*q.q2);
            R(3,2) = 2*(q.q2*q.q3 - q.q0*q.q1);
            R(3,3) = 2*(q.q0^2 + q.q3^2) - 1;
            dcm = R';
        end
        function [Q] = Qmultiply(obj,p,q)
            Q.q0 = p.q0*q.q0 - p.q1*q.q1 - p.q2*q.q2 - p.q3*q.q3;
            Q.q1 = p.q0*q.q1 + p.q1*q.q0 + p.q2*q.q3 - p.q3*q.q2;
            Q.q2 = p.q0*q.q2 - p.q1*q.q3 + p.q2*q.q0 + p.q3*q.q1;
            Q.q3 = p.q0*q.q3 + p.q1*q.q2 - p.q2*q.q1 + p.q3*q.q0;
        end
        function [y,Y,P,Y1] = ut(obj,X,Wm,Wc,n,Q,u,A,B) 
            L = size(X,2);  
            y = zeros(n,1);  
            Y = zeros(n,L);  
            for k=1:L  
                Y(:,k) = X(:,k) + A*X(:,k) + B*u;%f(X(:,k));  
                y = y+Wm(k)*Y(:,k);  
            end  
            Y1 = Y-y(:,ones(1,L));  
            P = Y1*diag(Wc)*Y1' + Q;
        end
        function [y,Y,P,Y1] = ut2(obj,X,Wm,Wc,n,R,H) 
            L = size(X,2);  
            y = zeros(n,1);  
            Y = zeros(n,L);  
            for k=1:L  
                Y(:,k) = H*X(:,k);%f(X(:,k));  
                y = y+Wm(k)*Y(:,k);  
            end  
            Y1 = Y-y(:,ones(1,L));  
            P = Y1*diag(Wc)*Y1' + R;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Sigma points around reference point  
        % Inputs:  
        % x: reference point  
        % P: covariance  
        % c: coefficient  
        % Output:  
        % X: Sigma points  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function X = sigmas(obj,x,P,c) 
              A = c*chol(P)';  %tranditional
%             [U,S,~] = svd(P);   %svd method
%             A = c*U*(S^0.5);    %svd method
%             L = numel(x);
%             A = zeros(L,L);
%             for i = 1:L
%                 A(i,i) = c*sqrt(P(i,i));
%             end          
            Y = x(:,ones(1,numel(x)));  
            X = [x Y+A Y-A];  
        end
        
        function obj = UpdateIMU(obj, Gyroscope, Accelerometer, Magnetometer, SamplePeriod)
            % control input (angular velocities from gyroscopes)
            if (size(Gyroscope,1) == 3), u = Gyroscope;
            else u = Gyroscope';
            end
            % measurements/observations (acceleromeres)
            if (size(Accelerometer,1) == 3), z = Accelerometer*obj.g0;
            else z = Accelerometer'*obj.g0;
            end
            % measurements/observations (magnetormeter)
            if (size(Magnetometer,1) == 3), mag = Magnetometer;
            else mag = Magnetometer';
            end   
            
            realz = z;
            x = obj.state;
            x_last = x;
            Q_ = SamplePeriod^2 * obj.Q; %Process noise covariance with time dependent noise
            
            % "rotation operators, wx wy wz"
            C3X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

            % Model generation
            A = [zeros(3,3) -SamplePeriod*C3X; zeros(3,6)];
            B = [SamplePeriod*C3X; zeros(3,3)];
            %%%%%%%%%%%%%%%%%%
            %   UKF
            %%%%%%%%%%%%%%%%%%
            lambda = obj.alpha^2*(obj.L+obj.ki)-obj.L;        %scaling factor  
            c = obj.L+lambda;                                 %scaling factor  
            Wm = [lambda/c 0.5/c+zeros(1,2*obj.L)];           %weights for sigmas  
            Wc = Wm;  
            Wc(1) = Wc(1)+(1-obj.alpha^2+obj.beta);           %weights for covariance  
            c = sqrt(c);  
            % Kalman a priori prediction
            X = obj.sigmas(x,obj.P,c);  %X sigmas here is a list of vector
            [x_predict,~,P1,~] = obj.ut(X,Wm,Wc,obj.L,Q_,u,A,B);%unscented transformation of prediction
            X1=obj.sigmas(x_predict,P1,c);
            X2=X1-x_predict(:,ones(1,size(X1,2)));        %deviation of X1
            
            % recompute R using the error between acceleration and the model of g 
            % (estimate of the magnitude of a0 in a = a0 + g)
            a_predict = z - x_predict(1:3)*obj.g0;
            a_len = sqrt(a_predict'*a_predict);
            R = (a_len*obj.r_a2 + obj.r_acc2)*eye(3);
            [z1,~,P2,Z2] = obj.ut2(X1,Wm,Wc,obj.m,R,obj.H);    %unscented transformation of measurments 
            P12 = X2*diag(Wc)*Z2';                        %transformed cross-covariance  
            K = P12/P2;  
            % update a posteriori
            x = x_predict + K * (z-z1);
            % update a posteriori covariance
            obj.P = P1 - K * P2 * K'; % for using any K

            % normalization of x & P (divide by DCM vector length)
            dcm_vector_length = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
%             J_33 = [x(2)^2 + x(3)^2,    -x(1)*x(2),         -x(1)*x(3); ...
%                     -x(1)*x(2),         x(1)^2 + x(3)^2,    -x(2)*x(3); ...
%                     -x(1)*x(3),         -x(2)*x(3),         x(1)^2 + x(2)^2];        
%             J = [ J_33 / (dcm_vector_length^3), zeros(3,3); 
%                     zeros(3,3), eye(3)];

            % Laplace approximation of normalization function for x to P, J = Jacobian(f,x)
            % P_new = E[J*(x-x0)*(x-x0)'*J'] = J*E[(x-x0)*(x-x0)']*J' = J*P*J'
            %obj.P = J*obj.P*J';
            x(1:3) = x(1:3) ./ dcm_vector_length;
            obj.state = x;

            % update UX / C3X
            UX = [0 -(u(3)-x(6)) u(2)-x(5); 
                u(3)-x(6) 0 -(u(1)-x(4)); 
                -(u(2)-x(5)) u(1)-x(4) 0];
            C3X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];            
            % compute Euler angles (not exactly a part of the extended Kalman filter)
            % yaw integration through full rotation matrix
            %u_nb = u - x(4:6);
            if (true)         
                % alternative method estimating the whole rotation matrix
                % integrate full rotation matrix (using first row estimate in memory)
                x1 = obj.first_row - SamplePeriod*UX*obj.first_row; %rotate x1 by x1 x u_nb
                x1 = x1/norm(x1);
                x2 = C3X * x1; %second row x2 = (state x x1)
%                 x1 = C3X' * x2; %recalculate first row x1 = (x2 * state) (ensure perpendicularity)
%                 x1 = x1 / norm(x1); % normalize length
                dcm = [x1';x2';x(1:3)'];
                [Qdcm] = obj.QfromDCM(dcm);
                if (mag'*mag) ~= 0 
                    qtmp.q0 = 0;
                    qtmp.q1 = mag(1);
                    qtmp.q2 = mag(2);
                    qtmp.q3 = mag(3);
                    [qtmp] = obj.Qmultiply(Qdcm,qtmp);
                    Qdcm_cong.q0 = Qdcm.q0;
                    Qdcm_cong.q1 = -Qdcm.q1;
                    Qdcm_cong.q2 = -Qdcm.q2;
                    Qdcm_cong.q3 = -Qdcm.q3;
                    [qmag] = obj.Qmultiply(qtmp,Qdcm_cong);
                    mag_E = [qmag.q1;qmag.q2;qmag.q3];
                    mag_E = mag_E/norm(mag_E);%normalize, in earth frame
                    mag_E = [mag_E(1);mag_E(2);0];
                    mag_err = cross(mag_E,obj.mag_vector);%cal err
                    qtmp.q0 = 0;
                    qtmp.q1 = mag_err(1);
                    qtmp.q2 = mag_err(2);
                    qtmp.q3 = mag_err(3);
                    [qtmp] = obj.Qmultiply(Qdcm_cong,qtmp);
                    [mag_err_q] = obj.Qmultiply(qtmp,Qdcm);
                    mag_err_B = [mag_err_q.q1;mag_err_q.q2;mag_err_q.q3];%err in body frame
                    %complimentary filter
                    obj.mag_err_int = obj.mag_err_int + obj.Ki*mag_err_B*SamplePeriod;
                    mag_change = obj.mag_err_int + obj.Kp*mag_err_B;
                    
                    q_mag_change.q0 = sqrt(1 - norm(mag_change)^2);
                    q_mag_change.q1 = -mag_change(1);
                    q_mag_change.q2 = -mag_change(2);
                    q_mag_change.q3 = -mag_change(3);
                    
                    [qtmp] = obj.Qmultiply(q_mag_change,Qdcm);
                    q_mag_change.q1 = mag_change(1);
                    q_mag_change.q2 = mag_change(2);
                    q_mag_change.q3 = mag_change(3);
                    [Qdcm] = obj.Qmultiply(qtmp,q_mag_change);
                end
                obj.Qt = [Qdcm.q0 Qdcm.q1 Qdcm.q2 Qdcm.q3]';
                [euler] = computeAngles(Qdcm);
                obj.yaw = euler.heading/180*pi;
                obj.roll = euler.roll/180*pi;
                obj.pitch = euler.pitch/180*pi;
                [obj.Dcm] = obj.DCMfromQ(Qdcm);
                obj.first_row = obj.Dcm(1,:)';
                obj.a = obj.Dcm*realz - [0;0;obj.g0];
                obj.a = obj.a - obj.gbias*obj.g0;% test eliminate gbias.
                
            else
                % alternative method estimating the whole rotation matrix
                % integrate full rotation matrix (using first row estimate in memory)
                x1 = obj.first_row - SamplePeriod*UX*obj.first_row; %rotate x1 by x1 x u_nb
                x1 = x1/norm(x1);
                x2 = C3X * x1; %second row x2 = (state x x1)
                x1 = C3X' * x2; %recalculate first row x1 = (x2 * state) (ensure perpendicularity)
                x1 = x1/norm(x1); % normalize length
                yaw_g = atan2(x2(1),x1(1));
                obj.pitch = asin(-x(1));
                obj.roll = atan2(x(2),x(3));  
                cos_roll = cos(obj.roll);
                sin_roll = cos_roll*x(2)/x(3);
                %cos_pitch = x(2)/sin_roll;
                if (mag'*mag) ~= 0
                    %init to no rotate on yaw to find the mx,my,mz in the flat plane.
                    x2 = [0 cos_roll -sin_roll]';
                    x1 = -C3X * x2; %second row x2 = (state x x1)
                    magx = x1'*mag;
                    magy = x2'*mag;
                    yaw_m = atan2(-magy, magx);%'-' here for the reason that rotate frame not the vector.
                    %complementary filter
                    delta_yaw = yaw_m-obj.yaw;
                    if delta_yaw > pi
                        delta_yaw = delta_yaw - 2*pi;
                    elseif delta_yaw < -pi
                        delta_yaw = delta_yaw + 2*pi;
                    end
                    obj.err_b = obj.err_b + obj.Ki*(delta_yaw)*SamplePeriod;
                    obj.yaw = yaw_g + obj.Kp*(delta_yaw) + obj.err_b;
                else
                    obj.yaw = yaw_g;
                end
                %compute new pitch and roll angles from a posteriori states 
                cos_yaw = cos(obj.yaw);
                sin_yaw = sin(obj.yaw);
                cos_pitch = cos(obj.pitch);
                sin_pitch = -x(1);
                %[cos(obj.pitch)*cos(obj.yaw) -cos(obj.roll)*sin(obj.yaw)+sin(obj.roll)*(-x(1))*cos(obj.yaw) sin(obj.roll)*sin(obj.yaw)+cos(obj.roll)*(-x(1))*cos(obj.yaw)]';
                obj.first_row = [cos_pitch*cos_yaw -cos_roll*sin_yaw+sin_roll*(sin_pitch)*cos_yaw sin_roll*sin_yaw+cos_roll*(sin_pitch)*cos_yaw]';
                x2 = C3X*obj.first_row;
                % estimated non-gravitational acceleration in global frame
                obj.a = [obj.first_row';x2';x(1:3)']*z - [0;0;obj.g0]; % acceleration estimate (g reduced)  
            end
            obj.vel = obj.vel + obj.a*SamplePeriod;
            obj.dist = obj.dist + obj.vel*SamplePeriod + 0.5*obj.a*SamplePeriod^2;
        end
    end
end
