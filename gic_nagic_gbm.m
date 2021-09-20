% Euler GBM
clear all
rng('default');
%parameters
mu = 0.02; sigma= 0.15; X0 = 25;
T = 20; N = 10000; M=2500; dt = T/M; mu_bar = mu-0.5*sigma^2;
q        = 0.01:0.01:1 
n        = 100     %quantiles %deciles in my case
X = ones(M,N)*X0;
dW = sqrt(dt)*randn(M,N);
for i = 2:M
X(i,:) = X(i-1, :).*(1 + dt*mu + sigma*dW(i,:));
end
tnew=ones(M,1)*dt;
tc=cumsum(tnew);
EX=X0*exp(tc*mu);
plot(tc, X); 
%plot(tc,X(:,1),'r');
%   hold on
%   plot(tc,X(:,2),'b');
%   plot(tc,X(:,3),'r');
%   plot(tc,X(:,4),'c');
%   plot(tc,X(:,5),'g');
%   hold off
 xlabel('Time', 'FontSize', 14);
 ylabel('Wealth', 'FontSize', 14);
 title('Simulations of GBM - N=10^4, 5 randomly chosen trajectories ', 'FontSize', 15);
%plot(tc, EX, 'r-*')
%hold off
l=10;
Xini       = X(l*(M/T), :)  % wealth at time t_0  (as at t=0 everyone has X0 money)
Xend       = X(T*(M/T), :)    % wealth at time 1  (final time)
sXini = sort(Xini);
sXend = sort(Xend);
dataSec_ini = reshape(sXini,n,[]);
dataSec_end = reshape(sXend,n,[]);
Mini = mean(dataSec_ini)
Mend = mean(dataSec_end)

[Xini_sort id1] = sort (Xini)
Xend_sort = Xend(id1)
dataSections = reshape(Xend_sort,n,[]); 
M = mean(dataSections)

yini       = quantile (Xini, q) 
yend       = quantile (Xend, q) 

for p    = 1:n
g(p)     = Mend(p)/Mini(p) -1; %gic
%lg(p)     = log(yend(p))-log(yini(p));
ng(p)      = M(p)/Mini(p)-1;        %Xend_sort(p*(N/(n))) %nagic
gt(p)     = exp(mu_bar*(T)+norminv(q(p))*sigma*sqrt(T))/exp(mu_bar*l+norminv(q(p))*sigma*sqrt(l)) -1; %analytical gic
end
plot(q*100,gt*100,'-','LineWidth', 3);
hold on
plot(q*100,g*100,'-','LineWidth', 3);
hold on
%plot(q*100,lg,'-','LineWidth', 2);
%hold on 
plot(q*100,ng*100,'-','LineWidth', 2);
legend( {'Analytical GIC','GIC in GBM', 'NaGIC in GBM'}, 'Location','northwest', 'FontSize', 13)
hold off
xlabel('Quantile', 'FontSize', 13);
ylabel('Relative change in wealth (%)', 'FontSize', 13);
title(['t = ' , num2str(l), ' year and t^{\prime} = ', num2str(T), ' year'],  'FontSize', 15);
%}
