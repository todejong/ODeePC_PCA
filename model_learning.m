close all; clear all; clc; 

load('data/udata.mat')
Tu = Ttrain-N;
gamma = 0;                  % regularization parameter


%% System parameters
mu = 1;
Ts = 0.1; %Sampling time
q = 200; % number of principal components

%% Generate points 
u = u_train

x10 = [-3:0.5:3];
x20 = [-3:0.5:3];

Tx = length(x10)*length(x20);


k = 0
for i = 1:length(x10)
    for j = 1:length(x20)
        k = k+1;
        xx0{k} = [x10(i);x20(j)];
    end 
end
U = [];
Y = [];
X0 = []

for i = 1:Tu
U = [U,u(i:i+N-1)'];

end


for kk = 1:length(xx0)
    disp('kk = '+string(kk)+'/'+string(length(xx0)))
    x0 = xx0{kk};
    X0 = [X0,x0];

    for j = 1:Tu
        u = U(:,j);
        for k = 1:N
           if k == 1
                y(1,k) = x0(1) + Ts*x0(2);
                y(2,k) = x0(2) + Ts*(mu*(1-x0(1)^2)*x0(2)-x0(1)+(u(k)));
            else 
                y(1,k) = y(1,k-1) + Ts*y(2,k-1);
                y(2,k) = y(2,k-1) + Ts*(mu*(1-y(1,k-1)^2)*y(2,k-1)-y(1,k-1)+(u(k)));
            end 
        end 
    Y = [Y,y(1,:)'];
    end 

end 

%% Kernel width
 
sigma_x = sqrt(0.8)
sigma_u = sqrt(200)

%% define the Gram matrix
for i = 1:length(xx0)
        disp('i = '+string(i)+'/'+string(length(xx0)))

    % Gaussian kernel
    Kx(i,:) = exp(-sum((X0(:,i)-X0).^2,1)/sigma_x^2);
   
    % Hardy multiquadrics
    % Kx(i,:) = sqrt(1+sum((X0(:,i)-X0).^2,1)/sigma_x^2);

    % reverse multiquadrics
    % Kx(i,:) = 1./sqrt(1+sum((X0(:,i)-X0).^2,1)/sigma_x^2);

    % Matern32
    % Kx(i,:) =  Matern_32(X0(:,i),X0,sigma_x);

    % Matern52
    % Kx(i,:) =  Matern_52(X0(:,i),X0,sigma_x);
end



for i = 1:Tu
            disp('i = '+string(i)+'/'+string(Tu))
    % Gaussian kernel
    Ku(i,:) = exp(-sum((U(:,i)-U).^2,1)/sigma_u^2);
   
    % Hardy multiquadrics
    % Ku(i,:) = sqrt(1+sum((U(:,i)-U).^2,1)/sigma_u^2);

    % reverse multiquadrics
    % Ku(i,:) = 1./sqrt(1+sum((U(:,i)-U).^2,1)/sigma_u^2);

    % Matern32
    % Ku(i,:) =  Matern_32(U(:,i),U,sigma_u);

    % Matern52
    % Ku(i,:) =  Matern_52(U(:,i),U,sigma_u);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('test 1')

%% simulate system

CKU = cond(Ku)
CKX = cond(Kx)


yker = []
Ku_inv = inv(chol(Ku))*inv(chol(Ku)');
Kx_inv = inv(chol(Kx))*inv(chol(Kx)');
Kinv = kron(Kx_inv,Ku_inv);
% Kinv = kron(inv(Ku),inv(Kx));



K_prod = kron(Kx,Ku);
disp('test 2')

% K_centered = centerKernelMatrix(K_prod);

[K_centered,K_mean_row,K_mean_col,K_mean_all] = centerKernelMatrix(K_prod);
disp('test 3')

tic;
[Uk,Sk,Vk] = svds(K_centered,q) ;
t1 = toc;
disp('test 4')




%% check PCA
knew_x = k_new(U(:,1),X0(:,1),U,X0,sigma_u,sigma_x,K_mean_col,K_mean_all);

Phi = Phi(U(:,1),X0(:,1),U,X0,sigma_u,sigma_x,K_mean_col,K_mean_all,Vk(:,1:q));
PHI = Vk(:,1:q)'*K_centered;

PHIinv = pinv(PHI);
%% 
xini0 = [2;2]
xini = xini0

clear y

for k = 1:length(u_test)
    i
    if k == 1
        y(1,k) = xini(1) + Ts*xini(2);
        y(2,k) = xini(2) + Ts*(mu*(1-xini(1)^2)*xini(2)-xini(1)+(u_test(k)));
    else 
        y(1,k) = y(1,k-1) + Ts*y(2,k-1);
        y(2,k) = y(2,k-1) + Ts*(mu*(1-y(1,k-1)^2)*y(2,k-1)-y(1,k-1)+(u_test(k)));
    end 
end 



disp('test 3')


Yker = []


for j = 1:length(u_test)/N
disp('j = '+string(j)+'/'+string(length(u_test)/N))

xker = xini;
uker = u_test((j-1)*N+1:(j)*N)';

clear k
for i = 1:length(X0(1,:))
    % Gaussian:
    kx(i,1) = exp(-sum((X0(:,i)-xker).^2,1)/sigma_x^2);   

    % Hardy multiquadrics
    % kx(i,1) = sqrt(1+sum((X0(:,i)-xker).^2,1)/sigma_x^2);   

    % reverse multiquadrics
    % kx(i,1) = 1/sqrt(1+sum((X0(:,i)-xker).^2,1)/sigma_x^2);

    % Matern32
    % kx(i,:) =  Matern_32(X0(:,i),xker,sigma_x);

    % Matern52
    % kx(i,:) =  Matern_52(X0(:,i),xker,sigma_x);
end 
for i = 1:length(U(1,:))
    % Gaussian:
    ku(i,1) = exp(-sum((U(:,i)-uker).^2,1)/sigma_u^2);   

    % Hardy multiquadrics
    % ku(i,1) = sqrt(1+sum((U(:,i)-uker).^2,1)/sigma_u^2);   

    % reverse multiquadrics
    % ku(i,1) = 1/sqrt(1+sum((U(:,i)-uker).^2,1)/sigma_u^2);

    % Matern32
    % ku(i,:) =  Matern_32(U(:,i),uker,sigma_u);

    % Matern52
    % ku(i,:) =  Matern_52(U(:,i),uker,sigma_u);


end 

k = kron(kx,ku);
yker = Y*Kinv*k;

Yker = [Yker,[NaN(N*(j-1),1);xini(1);yker;NaN(length(u_test)-(j*N),1)]];

xini = [y(1,j*N);y(2,j*N)];

end 
disp('test 4')

%% make figures

curr_fig = figure;
curr_axes1=axes('Parent',curr_fig,'FontSize',11,'FontName','Times New Roman');
box(curr_axes1,'on');
hold(curr_axes1,'all');
%your plots
subplot(2,1,1)
hold on
% yline(10,'r--')
% yline(-10,'r--')
plot([0:length(y(1,:))],[xini0(1),y(1,:)],'k',LineWidth=2)
for i = 1:length(Yker(1,:))
    plot([0:length(Yker(:,i))-1],Yker(:,i),LineWidth=2)
    xline((i-1)*N,'k--')
end 
ylabel('$x_1$',Interpreter='latex')
axis tight
legend('$x_1$','$x_{1,pred}$',Interpreter='latex')
subplot(2,1,2)
plot([0:length(u_test)-1],u_test,'k',LineWidth=2)
ylabel('$u$',Interpreter='latex')
xlabel('$k$',Interpreter='latex')
grid on
 %your plots
set(gca,'TickLabelInterpreter','Latex')
set(curr_fig,'Units','centimeters','PaperSize',[20.98 29.68],'PaperUnits','centimeters','PaperPosition',[0 0 12 8])
% you can change 9 and 6.3 to change the ratios of the plot...
savefig('figures/kernel_smi_product.fig') %change it with the name you want to give to the .fig plot
print -depsc figures/kernel_smi_product %change it with the name you want to give to the .eps figure

YKinv = Y*Kinv;

save('control_DeePC\kernel','YKinv','X0','U','Tx','Tu','sigma_x','sigma_u','K_mean_col','K_mean_all')

YPHIinv = Y*PHIinv;

Yker = [];
xini0 = [2;2]
xini = xini0


for j = 1:length(u_test)/N
disp('j = '+string(j)+'/'+string(length(u_test)/N))

xker = xini;
uker = u_test((j-1)*N+1:(j)*N)';

clear Phi
Phi = Phi(uker,xker,U,X0,sigma_u,sigma_x,K_mean_col,K_mean_all,Vk(:,1:q));


yker = YPHIinv*Phi;

Yker = [Yker,[NaN(N*(j-1),1);xini(1);yker;NaN(length(u_test)-(j*N),1)]];

xini = [y(1,j*N);y(2,j*N)];

end 
disp('test 4')

%% make figures

curr_fig = figure;
curr_axes1=axes('Parent',curr_fig,'FontSize',11,'FontName','Times New Roman');
box(curr_axes1,'on');
hold(curr_axes1,'all');
%your plots
subplot(2,1,1)
hold on
% yline(10,'r--')
% yline(-10,'r--')
plot([0:length(y(1,:))],[xini0(1),y(1,:)],'k',LineWidth=2)
for i = 1:length(Yker(1,:))
    plot([0:length(Yker(:,i))-1],Yker(:,i),LineWidth=2)
    xline((i-1)*N,'k--')
end 
ylabel('$x_1$',Interpreter='latex')
axis tight
legend('$x_1$','$x_{1,pred}$',Interpreter='latex')
subplot(2,1,2)
plot([0:length(u_test)-1],u_test,'k',LineWidth=2)
ylabel('$u$',Interpreter='latex')
xlabel('$k$',Interpreter='latex')
grid on
 %your plots
set(gca,'TickLabelInterpreter','Latex')
set(curr_fig,'Units','centimeters','PaperSize',[20.98 29.68],'PaperUnits','centimeters','PaperPosition',[0 0 12 8])
% you can change 9 and 6.3 to change the ratios of the plot...
savefig('figures/kernel_smi_product_PCA.fig') %change it with the name you want to give to the .fig plot
print -depsc figures/kernel_smi_product_PCA %change it with the name you want to give to the .eps figure

YKinv = Y*Kinv;

Theta = YPHIinv;

save('control_DeePC\kernel','YKinv','X0','U','Tx','Tu','sigma_x','sigma_u','K_mean_col','K_mean_all','q','Vk','PHI','Y','Theta')
