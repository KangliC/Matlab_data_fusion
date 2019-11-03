function [EulerAng]=computeAngles(Qi)
    Q.q4 = Qi.q3;
    Q.q3 = Qi.q2;
    Q.q2 = Qi.q1;
    Q.q1 = Qi.q0;
	roll = atan2(Q.q1*Q.q2 + Q.q3*Q.q4, 0.5 - Q.q2*Q.q2 - Q.q3*Q.q3);
	pitch = asin(-2.0 * (Q.q2*Q.q4 - Q.q1*Q.q3));
	heading = atan2(Q.q2*Q.q3 + Q.q1*Q.q4, 0.5 - Q.q3*Q.q3 - Q.q4*Q.q4);
    
% 	heading = atan2(2*Q.q2*Q.q3 - 2*Q.q1*Q.q4, 2*Q.q2*Q.q2 + 2*Q.q1*Q.q1-1);
% 	pitch = -asin(2*Q.q2*Q.q4 + 2*Q.q1*Q.q3);
% 	roll = atan2(2*Q.q3*Q.q4 - 2*Q.q1*Q.q2, 2*Q.q1*Q.q1 + 2*Q.q4*Q.q4-1);
    EulerAng.roll = roll/pi*180;
    EulerAng.pitch = pitch/pi*180;
    EulerAng.heading = heading/pi*180;
return
