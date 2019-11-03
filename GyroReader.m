%%close all;
delete(instrfindall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt_limt = 30000;
sensorData = [];
i = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while( i < cnt_limt)
    data = fscanf(s,'%f %f %f %f %f %f %f %f %f\r\n');
    
    if length(data) ~= 9
        disp(data);
        continue;
    end
    sensorData = [sensorData data];
    i = i+1;
end

fclose(s);
delete(s);
clear s;
delete(instrfindall);

