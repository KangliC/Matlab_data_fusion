delete(instrfindall);
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf, 'Color', 'White');

subplot(2,2,1); hold on; box on;
accel_1 = animatedline('MaximumNumPoints',1000,'Color','r', 'LineWidth', 1);
accel_2 = animatedline('MaximumNumPoints',1000,'Color','g', 'LineWidth', 1);
accel_3 = animatedline('MaximumNumPoints',1000,'Color','b', 'LineWidth', 1);
ylim([-4 4]);
xlim([0 500]);

subplot(2,2,2); hold on; box on;
gyro_1 = animatedline('MaximumNumPoints',1000,'Color','r', 'LineWidth', 1);
gyro_2 = animatedline('MaximumNumPoints',1000,'Color','g', 'LineWidth', 1);
gyro_3 = animatedline('MaximumNumPoints',1000,'Color','b', 'LineWidth', 1);
ylim([-10 10]);
xlim([0 500]);

subplot(2,2,3); hold on; box on;
mag_1 = animatedline('MaximumNumPoints',1000,'Color','r', 'LineWidth', 1);
mag_2 = animatedline('MaximumNumPoints',1000,'Color','g', 'LineWidth', 1);
mag_3 = animatedline('MaximumNumPoints',1000,'Color','b', 'LineWidth', 1);
ylim([-100 100]);
xlim([0 500]);

subplot(2,2,4); hold on; box on;
temp = animatedline('MaximumNumPoints',1000,'Color','r', 'LineWidth', 1);
ylim([20 40]);
xlim([0 500]);
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
a = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = serial('COM6');
s.BaudRate = 3000000;
s.DataBits = 8;
s.ByteOrder = 'littleEndian';
s.Parity = 'even';
s.StopBits = 1;
s.Terminator = 'CR/LF';
s.FlowControl = 'none';
s.InputBufferSize = 10485760;
s.Timeout = 30; 
s.BytesAvailableFcnMode = 'terminator';
fopen(s);
%fprintf(s,'hello');
s.TransferStatus;
get(s,{'InputBufferSize','BytesAvailable'});
while(t < T)
     data = fscanf(s,'%f %f %f %f %f %f %f %f %f %f\r\n');

    if length(data) ~= 10
        t = t+1;
        disp(data);
        continue;
    end
     b = toc(a);
     
     subplot(2,2,1);
      if t>500
        xlim([t-500 t]); 
     end    
     addpoints(accel_1,t,data(1,1));
     addpoints(accel_2,t,data(2,1));
     addpoints(accel_3,t,data(3,1));
     if b > 1/10
        drawnow;
        a = tic;
     end
     subplot(2,2,2);
     if t>500
        xlim([t-500 t]); 
     end     
     addpoints(gyro_1,t,data(4,1));
     addpoints(gyro_2,t,data(5,1));
     addpoints(gyro_3,t,data(6,1));
     if b>1/10
        drawnow;
     end
     subplot(2,2,3);
     if t>500
        xlim([t-500 t]); 
     end     
     addpoints(mag_1,t,data(7,1));
     addpoints(mag_2,t,data(8,1));
     addpoints(mag_3,t,data(9,1));
     if b > 1/10
        drawnow;
     end
%      subplot(2,2,4);
%      if t>500
%         xlim([t-500 t]); 
%      end
%      addpoints(temp,t,data(10,1));
%      if b>1/10
%         drawnow;
%      end
     t = t+1;
%      ax = [ax data(1,1)];
%      ay = [ay data(2,1)];
%      az = [az data(3,1)];
%      gx = [gx data(4,1)];
%      gy = [gy data(5,1)];
%      gz = [gz data(6,1)];
%      
%      if length(ax) == 500000
%          break;
%      end
end

fclose(s);
delete(s);
clear s;
