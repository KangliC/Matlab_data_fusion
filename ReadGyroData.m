delete(instrfindall);
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf, 'Color', 'White');
hold on; box on;
at_roll = animatedline('MaximumNumPoints',1000,'Color','r', 'LineWidth', 1);
at_pitch = animatedline('MaximumNumPoints',1000,'Color','g', 'LineWidth', 1);
at_heading = animatedline('MaximumNumPoints',1000,'Color','b', 'LineWidth', 1);
ylim([-180 180]);
xlim([0 100]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;
T = 1000000000;
A = [];
cnt = 1;
ax = [];
ay = [];
az = [];
gx = [];
gy = [];
gz = [];
mx = [];
my = [];
mz = [];
Q.q1 = 1;
Q.q2 = 0;
Q.q3 = 0;
Q.q4 = 0;
pitch = 0;
roll = 0;
heading = 0;
EulerAng.roll = 0;
EulerAng.pitch = 0;
EulerAng.heading = 0;
gyroMeasError = pi*2/180;
gyroMeasDrift = pi*0.5/180;
beta = sqrt(3/4)*gyroMeasError;
zeta = sqrt(3/4)*gyroMeasDrift;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = serial('COM3');
s.BaudRate = 3000000;
s.DataBits = 8;
s.ByteOrder = 'littleEndian';
s.Parity = 'even';
s.StopBits = 1;
s.Terminator = 'CR/LF';
s.FlowControl = 'none';
s.InputBufferSize = 10485760;
s.Timeout = 1; 
s.BytesAvailableFcnMode = 'terminator';
fopen(s);
%fprintf(s,'hello');
s.TransferStatus;
get(s,{'InputBufferSize','BytesAvailable'});
while(t < T)
	data = fscanf(s,'%f %f %f %f %f %f %f %f %f\r\n');
    if length(data) ~= 9
        t = t+1;
        disp(data);
        continue;
    end
	Q = ahrs_filter(Q,data(4,1),data(5,1),data(6,1),data(1,1),data(2,1),data(3,1),data(7,1),data(8,1),data(9,1),0.001,zeta,beta);
    EulerAng = computeAngles(Q);
    if t>500
        xlim([t-100 t]); 
    end
    addpoints(at_roll,t,EulerAng.roll);
    addpoints(at_pitch,t,EulerAng.pitch);
    addpoints(at_heading,t,EulerAng.heading);
    drawnow limitrate;
    roll = EulerAng.roll;
    pitch = EulerAng.pitch;
    heading = EulerAng.heading;
    t = t+1;
     
%      if length(ax) == 50000
%          break;
%      end
end

fclose(s);
delete(s);
clear s;

