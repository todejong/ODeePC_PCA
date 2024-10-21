close all; clear all; clc; 

load('kernel.mat')

w = warning ('off','all');


nx = 2
nu = 1;


%% System parameters

mu = 1;
Ts = 0.1; %Sampling time


%--------------------------------------------------------------------------
% 0.01<=x1<=2, -20<=x2<=0, 0<=u<=1
y_min = -3;
y_max = 3;
u_min = -5;
u_max = 5;

N = 10;

wQ = 1000;
wX = 10*eye(nx);
wR = 10;
lambda = 10;



%% Implicit Yalmip

u = sdpvar(N,1);
y = sdpvar(N,1);
x0 = sdpvar(nx,1);
ref = sdpvar(N,1);
gvar=sdpvar(10,1);


objective = (y-ref)'*wQ*(y-ref)+(u)'*wR*(u) + lambda*gvar'*gvar;

constraints = [[PHI;ones(1,length(PHI(1,:)))]*pinv(Y)*gvar == zeros(q+1,1)];
constraints = [constraints, Y*pinv([PHI;ones(1,length(PHI(1,:)))])*[Phi(u,x0,U,X0,sigma_u,sigma_x,K_mean_col,K_mean_all,Vk(:,1:q));1] + gvar == y  ];


    
% constraints = [ PhiY*gvar == [Phi(u,x0,U,X0,V,sigma_u,sigma_x);y]];

for k = 1:N
    % constraints = [constraints,  -3<=u(k)<=3];
end
Parameters = {x0,ref};
Outputs = {u};

options = sdpsettings('solver', 'fmincon', 'verbose', 0, 'debug', 0);
controller = optimizer(constraints, objective, options, Parameters, Outputs);





%% Simulation
x0 = [0.5;0];

fs = 10;                    % Sampling frequency (samples per second)
Ts = 1/fs;                   % seconds per sample
StopTime = 25;                % seconds
t = (0:Ts:StopTime-Ts)';        % seconds
F = 0.1;
R = 1*sin(2*pi*F*t)';% works 


% R = [zeros(1,50),2*ones(1,150)];
% D = 0.1*[zeros(1,length(R)-50),ones(1,50)]

k_sim = length(R)-N;

xk = x0*ones(1,k_sim+1);
uk = zeros(nu,k_sim);
yk = xk(1,:);
tk = zeros(1,k_sim);
r_plot = 0;

t = 0;

u0 = x0(1);

for k = 1:k_sim
    disp(string(k) + '/' + string(k_sim))
    t(k+1) =  Ts*k;
    xini = xk(:,k);

    r = R(k:k+N-1)';
    r_plot(k+1) = r(1);
    tic;
    %%% mpc sol
    
    Uk = controller({xk(:,k),r});
    u_mpc =  Uk(1:nu);
    tk(:,k) = toc;
    uk(:,k) = u_mpc;

    xk(1,k+1) = xk(1,k) + Ts*xk(2,k);
    xk(2,k+1) = xk(2,k) + Ts*(mu*(1-xk(1,k)^2)*xk(2,k)-xk(1,k)+(uk(:,k))); % + D(k);
    yk(k+1) = xk(1,k+1);

end

%%

curr_fig = figure;
curr_axes1=axes('Parent',curr_fig,'FontSize',11,'FontName','Times New Roman');
box(curr_axes1,'on');
hold(curr_axes1,'all');
%your plots
subplot(2,1,1)
hold on
plot(r_plot(1:end),LineWidth=1.5,LineStyle="--",Color=[0 0 0])
plot(yk,LineWidth=1.5)
axis tight
grid on
ylabel('$y$',Interpreter='latex')
subplot(2,1,2)
plot(uk,LineWidth=1.5)
ylabel('$u$',Interpreter='latex')
xlabel('$t [s]$',Interpreter='latex')
grid on
axis tight
%your plots
set(gca,'TickLabelInterpreter','Latex')
set(curr_fig,'Units','centimeters','PaperSize',[20.98 29.68],'PaperUnits','centimeters','PaperPosition',[0 0 12 8])
% you can change 9 and 6.3 to change the ratios of the plot...
savefig('product_kernel_control.fig') %change it with the name you want to give to the .fig plot
print -depsc product_kernel_control %change it with the name you want to give to the .eps figure


mean(tk)

save('plot_data','r_plot', 'yk','uk') 