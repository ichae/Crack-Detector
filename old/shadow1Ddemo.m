%% shadow removal demo
clc
clear all
close all

x = 0:0.5:100;
mu = 80;
sigma = 60;
bias = 0.6;
bg = bias + .5*exp(-((x-mu)/sigma).^2);   % actual shadow function
signal = 1- 0.15*(heaviside(x-10)-heaviside(x-15) + heaviside(x-50)-heaviside(x-60) + heaviside(x-80)-heaviside(x-85)); % features
f = bg.*(signal);

%% Case 1: no noise
% Estimate shadow
f = f/max(f);
s = 20;
sz = 6*s;
h = fspecial('gaussian',sz,s);
shadow = imfilter(f,h,'same','replicate');

u = f./shadow;

%Plottings
close all;
figure(5); 
plot(x,f,'-r','LineWidth',2);  hold on;
plot(x,shadow,'-b','LineWidth',2); hold on;
plot(x,u,'-k','LineWidth',2); hold on;
legend('f(x)','s(x)','u(x)');
hold off;
%% Case 2: Noisy
f = bg.*(signal) + 0.06*rand(1,length(f));
f = f/max(f);

% s = 10;
% sz = 6*s;
% h = fspecial('gaussian',sz,s);
% shadow = imfilter(f,h,'same','replicate');
shadow = 1;
iter = 500;
u = f;
lambda = 0.6;
figure(6)
for ii = 1 : iter
    u(1) = u(2); u(end) = u(end-1);
    uxx = gradient(gradient(u));
    
    term1 = (f-u.*shadow).*shadow;
    term2 = uxx;
    
    term1 = term1./(max(abs(term1(:))));
    term2 = term2./(max(abs(term2(:))));
    
    chng = (1-lambda)*term1+ lambda*term2;
    dt = 0.01/max(abs(chng(:)));
%     dt = 0.1;
    u = u + dt*chng;
%     u = u/max(u(:));
%     plot(x,u,'-k','LineWidth',3);drawnow; 
end

plot(x,f,'-r','LineWidth',2);drawnow; hold on; 
plot(x,shadow,'-b','LineWidth',2);drawnow; hold on;
plot(x,u,'-k','LineWidth',2);
% plot(x,f./shadow,'-m','LineWidth',2);
legend('f(x)','s(x)','u(x)','f(x)/s(x)');
hold off;
