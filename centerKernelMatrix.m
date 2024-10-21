function [K_centered,K_mean_row,K_mean_col,K_mean_all] = centerKernelMatrix(K)
    % Input:
    % K - n x n Gram matrix (kernel matrix)
    
    % Get the number of data points (size of the Gram matrix)
    n = size(K, 1);
    
    % Create a column vector of ones
    one_n = ones(n, 1);
    
    % Compute the mean of rows and columns
    K_mean_row = mean(K, 1);  % Row-wise mean (1 x n)
    K_mean_col = mean(K, 2);  % Column-wise mean (n x 1)
    K_mean_all = mean(K_mean_row);  % Grand mean (scalar)
    
    % Center the Gram matrix
    K_centered = K - K_mean_row - K_mean_col + K_mean_all;
end
