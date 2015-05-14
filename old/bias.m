
X = 0:0.1:100;

signal = 1+3*exp(-X.^2/35) + 1*rand(1,length(X))+ 2*(1-exp(-(X-480)/500)) ;
% plot(X,signal);

% b = 0.5+zeros(1,length(X));
% b(1,1:560) = 3;

b = 0.1*X+2;

biased_signal = b+signal;

subplot(1,3,1); plot(X,signal);
subplot(1,3,2); plot(X,b);
subplot(1,3,3); plot(X,biased_signal);

%%
h = 