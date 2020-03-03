classdef PI_IMU6 < handle
    
    properties (Access = public)
        rot3 = [0 0 1];
        a = zeros(3,1);
        q0 = 1;
        q1 = 0;
        q2 = 0;
        q3 = 0;        
        Ki = 0.001;
        Kp = 0.3;
        errInt = [0 0 0];
        yaw = 0;
        pitch = 0;
        roll = 0;
    end
    
    methods (Access = public)
        function obj = PI_IMU6(varargin)
        
        end
        function obj = UpdateIMU6(obj, gyro, accel, dt)
            acc = accel / norm(accel);
            err = cross(acc, obj.rot3);
            obj.errInt = obj.errInt + obj.Ki*err*dt;
            gy = gyro + obj.Kp*err + obj.errInt;
            
            q0Last = obj.q0;
            q1Last = obj.q1;
            q2Last = obj.q2;
            q3Last = obj.q3;
            
            obj.q0 = obj.q0 + 0.5*(-q1Last * gy(1) - q2Last * gy(2) - q3Last * gy(3))*dt;
            obj.q1 = obj.q1 + 0.5*( q0Last * gy(1) + q2Last * gy(3) - q3Last * gy(2))*dt;
            obj.q2 = obj.q2 + 0.5*( q0Last * gy(2) - q1Last * gy(3) + q3Last * gy(1))*dt;
            obj.q3 = obj.q3 + 0.5*( q0Last * gy(3) + q1Last * gy(2) - q2Last * gy(1))*dt;
            
            normval = norm([obj.q0 obj.q1 obj.q2 obj.q3]);
            obj.q0 = obj.q0/normval;
            obj.q1 = obj.q1/normval;
            obj.q2 = obj.q2/normval;
            obj.q3 = obj.q3/normval;
            
            obj.rot3(1) = 2*(obj.q1*obj.q3 - obj.q0*obj.q2);
            obj.rot3(2) = 2*(obj.q2*obj.q3 + obj.q0*obj.q1);
            obj.rot3(3) = 1 - 2*(obj.q1)^2 - 2*(obj.q2)^2;
            Q.q0 = obj.q0;
            Q.q1 = obj.q1;
            Q.q2 = obj.q2;
            Q.q3 = obj.q3;
            euler = computeAngles(Q);
            obj.yaw = euler.heading;
            obj.roll = euler.roll;
            obj.pitch = euler.pitch;
        end
    end
end
