classdef ACCEL_KF < handle
    properties (Access = public)
        state = [0 0 0 0 0 0]';
        vcov_init = 0.01^2;
        abcov_init = [0.001^2 0 0;
                      0 0.001^2 0;
                      0 0 0.002^2];%0.005^2;
        vcov = 0.01^2;
        abcov = [0.0001 0 0;
                 0 0.0001 0;
                 0 0 0.0001];%0.01^2;
        vmcov = 0.120^2;
        P = [];
        Q = [];
        R = [];
    end
    
    methods (Access = public)
        function obj = ACCEL_KF(varargin)
            obj.P = [eye(3)*obj.vcov_init zeros(3,3);
                     zeros(3,3) eye(3)*obj.abcov_init];
            obj.Q = [eye(3)*obj.vcov zeros(3,3);
                     zeros(3,3) eye(3)*obj.abcov];
            obj.R = eye(3)*obj.vmcov;
        end
        % x = A * x + B * u; x = [Vx Vy Vz abx aby abz]'
        % y = H * x
        function obj = UpdateKF(obj, accel, dt)
             A = [1 0 0 -dt 0 0;
                 0 1 0 0 -dt 0;
                 0 0 1 0  0 -dt;
                 0 0 0 1  0  0;
                 0 0 0 0  1  0;
                 0 0 0 0  0  1];
             B = [dt 0  0;
                  0  dt 0;
                  0  0  dt;
                  0  0  0;
                  0  0  0;
                  0  0  0];
             H = [1 0 0 0 0 0;
                  0 1 0 0 0 0;
                  0 0 1 0 0 0];
             Q_ = dt^2 * obj.Q;
             x = A * obj.state + B * accel;
             y_ = H * x;
             y = zeros(3,1);%0*randn(3,1);%true measurement, it is 0 matrix here.
             Pk = A*obj.P*A' + Q_;
             K = Pk*H'/(H*Pk*H'+obj.R);
             obj.state = x + K*(y-y_);%measurement y is 0 always in this case
             % obj.P = IKH * P_predict * IKH' + K * R * K'; % for using any K
             IKH = eye(6) - K * H;
             obj.P = IKH * Pk * IKH' + K * obj.R * K';
        end
    end
end
