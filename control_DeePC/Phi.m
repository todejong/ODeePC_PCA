function Phi = Phi(uker,xker,U,X0,sigma_u,sigma_x,K_mean_col,K_mean_all,V)

Tu = length(U(1,:));
Tx = length(X0(1,:));
k = [];
for i = 1:Tu
    % Gaussian kernel
    ku(i,1) = exp(-sum((U(:,i)-uker).^2,1)/sigma_u^2);

    % reverse multiquadrics
    % ku(i,1) = 1/sqrt(1+sum((U(:,i)-uker).^2,1)/sigma_u^2);
end 

for i = 1:Tx
    % Gaussian kernel
    kx(i,1) = exp(-sum((X0(:,i)-xker).^2,1)/sigma_x^2);

    % reverse multiquadrics
    % kx(i,1) = 1/sqrt(1+sum((X0(:,i)-xker).^2,1)/sigma_x^2);
end 
for i = 1:length(kx)
    for j = 1:length(ku)
        k = [k;kx(i)*ku(j)];
    end 
end 

knew_x = k-mean(k)-K_mean_col+K_mean_all;
Phi = V'*knew_x;
end

