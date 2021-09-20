% Euler RGBM
clear all
rng('default');
%parameters
mu = 0.02; sigma= 0.15; X0 = 25; mu_bar = mu-0.5*sigma^2;
T = 20; N = 10000; M=2500; dt = T/M; tau= 0.2; %reallocation 
q        = 0.01:0.01:1;
n        = 100     %quantiles %deciles in my case
X        = ones(M,N)*X0;
%Z       = ones(M,N)*X0;

dW = sqrt(dt)*randn(M,N);
for i  = 2:M
X(i,:) = X(i-1,:)+X(i-1,:).*(mu*dt+sigma*dW(i,:))-tau*(X(i-1,:)-mean(X(i-1,:)))*dt; %RGBM
%Z(i,:) = Z(i-1,:)+Z(i-1,:).*(mu*dt+sigma*dW(i,:)); %GBM
end
tnew   = ones(M,1)*dt;
tc     = cumsum(tnew);
plot(tc, X); 
l=10;

Xini       = X(l*(M/T), :)  % wealth at time t=1  (as at t=0 everyone has X0 money)
Xend       = X(T*(M/T), :)    % wealth at tim1 t=2 (final time)
sXini = sort(Xini);
sXend = sort(Xend);
dataSec_ini = reshape(sXini,n,[]);
dataSec_end = reshape(sXend,n,[]);
Mini = max(dataSec_ini)
Mend = max(dataSec_end)

yini       = quantile (Xini, q) 
yend       = quantile (Xend, q) 

% Zini       = Z(l*(M/T), :)  % wealth at time t_0  (as at t=0 everyone has X0 money)
% Zend       = Z(T*(M/T), :)    % wealth at tim1 t_end (final time)
% uini       = quantile (Zini, q) 
% uend       = quantile (Zend, q) 
% sZini = sort(Zini);
% sZend = sort(Zend);
% dataSec_zini = reshape(sZini,n,[]);
% dataSec_zend = reshape(sZend,n,[]);
% Uini = mean(dataSec_zini)
% Uend = mean(dataSec_zend)

[Xini_sort id1] = sort (Xini)     % for NAGIC 
Xend_sort = Xend(id1)
dataSections = reshape(Xend_sort,n,[]); 
M = mean(dataSections)

%{
mean_ini   =  mean(log(Xini)) % for theor Nagic
mean_end   =  mean(log(Xend))
std_ini = std(log(Xini))
std_end = std(log(Xend))
rho = corrcoef(log(Xini),log(Xend))
%}

for p      = 1:n
%gu(p)       = Uend(p)/Uini(p) -1; % for GBM GIC overlay
g(p)       = Mend(p)/Mini(p) -1;
ng(p)      = M(p)/Mini(p)-1;          %Xend_sort(p*(N/(n+1)))
end

%plot(q*100,gu, '-','LineWidth', 2);
plot(q*100,g,'-','LineWidth', 2);
%set(gca,'YLim',[-1.5 1.5])
hold on
plot(q*100,ng,'-','LineWidth', 2);
%set(gca,'YLim',[-1.5 1.5])
%plot(q*100,gt,'-','LineWidth', 2);
legend({'GIC in RGBM','NAGIC in RGBM'},'Location','best', 'FontSize', 13)
hold off
xlabel('Quantiles', 'FontSize', 13);
ylabel('Relative Chnage in Wealth ', 'FontSize', 13);
title(['t = ' , num2str(l), ' and t^{\prime} = ', num2str(T)],'FontSize', 13);
