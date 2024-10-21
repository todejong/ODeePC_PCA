close all; clear all; clc; 

Range = [-2,2]
SineData = [25, 40, 1]
Band = [0, 1]
NumPeriod = 1
Nu = 1
TrainTest = 2/3
N = 10;
Ttrain = 50+N

Ttest = 1000
Period = Ttrain+Ttest


u_data = idinput([Period 1 NumPeriod],'sine',Band,Range,SineData)';
u_train = u_data(1:Ttrain)
u_test = u_data(Ttrain+1:Ttrain+Ttest)



mu = 1;
Ts = 0.1; %Sampling time


xk = [0;0]
xvec = []

for k = 1:length(u_data)
xk = [xk(1) + Ts*xk(2);  xk(2)+Ts*(mu*(1-xk(1)^2)*xk(2)-xk(1)+(u_data(k)))];
xvec = [xvec xk];

end 

figure()
subplot(3,1,1)
plot(xvec(1,:))
subplot(3,1,2)
plot(xvec(2,:))
subplot(3,1,3)
plot(u_data)


save('data/udata','u_test','u_train',"N",'Ttrain','Ttest')