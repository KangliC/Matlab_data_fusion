%============================================================================
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%============================================================================

classdef DCM_IMU9 < handle
%   DCM_IMU9 Implementation of IMU algorithm
%
%   Date          Author         Source           Notes    
%   09/12/2019    Kangli Chu     Heikki Hyyti     Initial release

    %% Public properties
    properties (Access = public)
        g0 = 9.8;                   % gravitation around Helsinki, Finland (change according to your area)
        B0 = 70;                    % earth magnetic value in uT
        state = [1 0 0 0 0 1 0 0 0]'; % States are lowest row of rotation matrix and gyroscope x y and z biases
                                    % (C_11, C_12, C_13, C_31, C_32, C_33, w_b1, w_b2, w_b3)
        q_dcm2 = 0.00003;            % estimated variance of dcm states (gyro variance per second)
        q_gyro_bias2 = 0.001^2;     % very small number to make bias change slowly
        r_acc2 = 0.5^2;             % variance of calibrated accelerometer (g-component)
        r_a2 = 10^2;                % large variance for some unknown acceleration (acc = a + g)
        r_mag2 = 0.806;
        r_m2 = 250;
        q_dcm2_init = 1^2;          % initial variance of dcm states (for attitude estimation)
        q_gyro_bias2_init = 0.1^2;  % initial variance of bias states (for bias estimator)
        a = zeros(3,1);             % estimated non-gravitational accelerations
        vel = zeros(3,1);           % estimated velocity on earth frame
        dist = zeros(3,1);          % estimated distance from origin on earth frame;
        m = zeros(3,1);             % estimated non-earth magnetical disturbance
        yaw = 0;                    % Yaw angle around z axis (in ZYX convention)
        pitch = 0;                  % Pitch angle around y axis
        roll = 0;                   % Roll angle around x axis distortion
        P = [];                     % estimate covariance (these are initialized in constructor below)
        H = [];                     % observation model (static)
        Q = [];                     % proces noise covariance (static part)
        Rot = [1 0 0; 0 1 0; 0 0 1];% initial rotation matrix
        first_row = [1 0 0]';       % first row of of the rotation matrix (for yaw angle estimate)
        mag_sensor = [0 0 0]';
    end

    %% Public methods
    methods (Access = public)
        function obj = DCM_IMU9(varargin)
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
            
            if (updateP), obj.P = [obj.q_dcm2_init*eye(3), zeros(3,3), zeros(3,3); zeros(3,3), obj.q_dcm2_init*eye(3), zeros(3,3); zeros(3,3), zeros(3,3), obj.q_gyro_bias2_init*eye(3)]; 
            end
            obj.H = [eye(3)*40.05, eye(3)*(-57.22), zeros(3,3); zeros(3,3), eye(3)*obj.g0, zeros(3,3)];
            obj.Q = [obj.q_dcm2*eye(3), zeros(3,3), zeros(3,3); zeros(3,3), obj.q_dcm2*eye(3), zeros(3,3); zeros(3,3), zeros(3,3), obj.q_gyro_bias2*eye(3)];
        end
        
        function obj = UpdateIMU9(obj, Gyroscope, Accelerometer, Magnetometer, SamplePeriod)
            %Data adapt%
            % measurements/observations (acceleromeres)
            if (size(Accelerometer,1) == 3), z = Accelerometer*obj.g0;
            else z = Accelerometer'*obj.g0;
            end            
            z = z*0.990257;%calibration simulation
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
                obj.mag_sensor = mag_tmp;%obj.r_m2 = 25;
            end
            mag = mag_tmp;
            
            % Kalman Filter Calculation %
            x = obj.state;
            %x_last = x;
            Q_ = SamplePeriod^2 * obj.Q; %Process noise covariance with time dependent noise
            
            % control input (angular velocities from gyroscopes)
            if (size(Gyroscope,1) == 3), u = Gyroscope;
            else u = Gyroscope';
            end
            
            % "rotation operators"
            C1X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
            C3X = [0 -x(6) x(5); x(6) 0 -x(4); -x(5) x(4) 0];
            
            UX = [0         -(u(3)-x(9))    u(2)-x(8); 
                  u(3)-x(9)     0         -(u(1)-x(7)); 
                -(u(2)-x(8))  u(1)-x(7)     0];

            % Model generation
            A = [zeros(3,3) zeros(3,3) -SamplePeriod*C1X; 
                 zeros(3,3) zeros(3,3) -SamplePeriod*C3X; 
                 zeros(3,9)];
            B = [SamplePeriod*C1X; 
                 SamplePeriod*C3X; 
                 zeros(3,3)];
            F = eye(9) + [-SamplePeriod*UX, zeros(3,3), -SamplePeriod*C1X; zeros(3,3), -SamplePeriod*UX, -SamplePeriod*C3X; zeros(3,9)];

            % Kalman a priori prediction
            x_predict = x + A*x + B*u;
            P_predict = F * obj.P * F' + Q_;

            % recompute R using the error between acceleration and the model of g 
            % (estimate of the magnitude of a0 in a = a0 + g)
            a_predict = z - x_predict(4:6)*obj.g0;
            a_len = sqrt(a_predict'*a_predict);
            m_len = abs(sqrt(mag'*mag)- obj.B0);   % deviation from earth magnetical value in uT
            
            R = [(m_len*obj.r_m2 + obj.r_mag2)*eye(3), zeros(3,3); zeros(3,3), (a_len*obj.r_a2 + obj.r_acc2)*eye(3)];

            % update magnetic related value in H
%             mag_earth = obj.Rot*mag;
%             mag_earth_x = sqrt(mag_earth(1)^2 + mag_earth(2)^2);%31.85;%
%             mag_earth_z = mag_earth(3);%-38.22;%
            %obj.H = [eye(3)*mag_earth_x, eye(3)*mag_earth_z, zeros(3,3); zeros(3,3), eye(3)*obj.g0, zeros(3,3)];
            % Kalman innovation
            y = [mag;z] - obj.H*x_predict;
            S = obj.H * P_predict * obj.H' + R;

            % Kalman gain
            K = P_predict * obj.H' / S;

            % update a posteriori
            x = x_predict + K * y;

            % update a posteriori covariance
            IKH = eye(9) - K*obj.H;
            obj.P = IKH * P_predict * IKH' + K * R * K'; % for using any K

            % normalization of x & P (divide by DCM vector length)
            dcm_vector_length = sqrt(x(1)^2 + x(2)^2 + x(3)^2);     % C1
            dcm_vector_length2 = sqrt(x(4)^2 + x(5)^2 + x(6)^2);    % C3
            J_33 = [ x(2)^2+x(3)^2,   -x(1)*x(2),       -x(1)*x(3); ...
                    -x(1)*x(2),        x(1)^2+x(3)^2,   -x(2)*x(3); ...
                    -x(1)*x(3),       -x(2)*x(3),        x(1)^2+x(2)^2];    
                
            J_33_2 = [x(5)^2+x(6)^2,    -x(4)*x(5),        -x(4)*x(6); ...
                     -x(4)*x(5),         x(4)^2+x(6)^2,    -x(5)*x(6); ...
                     -x(4)*x(6),        -x(5)*x(6),         x(4)^2+x(5)^2];   
            J = [J_33/(dcm_vector_length^3), zeros(3,3), zeros(3,3); ...
                zeros(3,3), J_33_2/(dcm_vector_length2^3), zeros(3,3); ...
                zeros(3,3), zeros(3,3), eye(3)];
            % Laplace approximation of normalization function for x to P, J = Jacobian(f,x)
            % P_new = E[J*(x-x0)*(x-x0)'*J'] = J*E[(x-x0)*(x-x0)']*J' = J*P*J'
            obj.P = J*obj.P*J';
            x(1:3) = x(1:3) ./ dcm_vector_length;
            x(4:6) = x(4:6) ./ dcm_vector_length2;           
            obj.state = x;
            
            % compute Euler angles (not exactly a part of the extended Kalman filter)
            % yaw integration through full rotation matrix
            u_nb = u - x(7:9);
            C3X = [0 -x(6) x(5); x(6) 0 -x(4); -x(5) x(4) 0];

            % alternative method estimating the whole rotation matrix
            % integrate full rotation matrix (using first row estimate in memory)
            x1 = obj.state(1:3);%obj.first_row + SamplePeriod*UX'*obj.first_row; %rotate x1 by x1 x u_nb
            x2 = C3X * x1; %second row x2 = (state x x1)
            x2 = x2 ./ sqrt(x2(1)^2 + x2(2)^2 + x2(3)^2); % normalize length of the second row
            %x1 = C3X' * x2; %recalculate first row x1 = (x2 * state) (ensure perpendicularity)
            obj.first_row = x1;% ./ sqrt(x1(1)^2 + x1(2)^2 + x1(3)^2); % normalize length
            x3 = obj.state(4:6);
            obj.Rot = [x1';x2';x3'];

            % update magnetic related value in H
%             mag_earth = obj.Rot*mag;
%             mag_earth_x = sqrt(mag_earth(1)^2 + mag_earth(2)^2);
%             mag_earth_z = mag_earth(3);%-38.22;%
%             obj.H = [eye(3)*mag_earth_x, eye(3)*mag_earth_z, zeros(3,3); zeros(3,3), eye(3)*obj.g0, zeros(3,3)];
            
            %compute new pitch and roll angles from a posteriori states
            obj.pitch = asin(-x(4));
            obj.roll = atan2(x(5),x(6));          
            obj.yaw = atan2(x2(1),x1(1));
            %x2 = C3X*x(1:3);
            % estimated non-gravitational acceleration in global frame
            obj.a = obj.Rot*z - [0;0;obj.g0]; % acceleration estimate (g reduced) 
            obj.vel = obj.vel + obj.a*SamplePeriod;
            obj.dist = obj.dist + obj.vel*SamplePeriod + 0.5*obj.a*SamplePeriod^2;
            % save the estimated non-gravitational acceleration
            %obj.a = z - x(4:6)*obj.g0; % acceleration estimate (g reduced) 
            %obj.m = mag - x(1:3);
        end
    end
end
