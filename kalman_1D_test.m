dArray = [];
dhArray = [];
dkArray = [];
tArray = [];
X = 2;
dt = 0.1;
T = 10;
P = 0;
K = 0;
M = 0;
Var_v = 3;
Var_w = 5;
x = 0;
y = 0;
close all;
figure;
set(gcf, 'Color', 'White');
box on;
h = animatedline('MaximumNumPoints',100, 'Color','g');
h2 = animatedline('MaximumNumPoints',100, 'Color','r');
%h = plot(x,y);
for t = 0:dt:T
    M = 10 + 2*rand() - 1;
    P = P + Var_w;
    K = P/(P + Var_v);
    X = X + K*(M - X);
    P = (1-K)*P;
    dhArray = [dhArray M];
    dArray = [dArray X];
    dkArray = [dkArray K];
    tArray = [tArray t];
    y = X;
    x = t;
    y2 = X+4;
    addpoints(h,x,y);
    addpoints(h2,x,y2);
    %refreshdata(h);
    drawnow;
    %X = int32(X);
end

% close all;
% figure;
% set(gcf, 'Color', 'White');
% hold on; box on;
% 
% subplot(3,1,1);hold on; box on;
% plot(tArray, dArray(1,:), 'b-', 'LineWidth', 2);
% subplot(3,1,2);hold on; box on;
% plot(tArray, dhArray(1,:), 'g-', 'LineWidth', 2);
% subplot(3,1,1);hold on; box on;
% plot(tArray, dkArray(1,:), 'r-', 'LineWidth', 2);

