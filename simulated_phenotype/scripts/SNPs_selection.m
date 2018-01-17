function [chosen_SNP_ids, chosen_LDscores] = SNPs_selection(X_data, n)

% INPUT:
% - genotypic data [samples x SNPs]
% - number of SNPs to be chosen

A           = corr(X_data);
variance    = ones(size(A));

current_ids = 1:size(A,1);


for s = 1:n
    
    % compute LD-score of each SNP
    l = sum( (A.^2).* variance  ,2);
    
    % find a SNP with the highest score
    [chosen_LDscores(s), temp_id] = max(l);
    
    chosen_SNP_ids(s) = current_ids(temp_id);
    current_ids(temp_id) = [];
    
    
    % PARTIAL CORRELATION
    
    % chosen SNP(s)
    Z = X_data(:, chosen_SNP_ids);
    
    X = X_data;
    X(:,chosen_SNP_ids) = []; % remove chosen variables (SNPs)
 
    [A, ~, resid, variance] = partialcorr_my(X, X, Z); 
    variance = repmat(variance(1:size(A,1)),[size(A,1), 1]);

    
% // the code below produces the same output (A, variance) but it is slower //  
%     A = nan(size(X,2), size(X,2));
%     variance = nan(size(A));
%     for i = 1:size(X,2)
%         for j = i:size(X,2)
%             i ,j 
%             [~,~,stats_x]  = glmfit(Z, X(:,i) );  
%             [~,~,stats_y]  = glmfit(Z, X(:,j) );
%             
%             % partial correlation -> correlation between residuals
%             A(i,j) = corr(stats_x.resid, stats_y.resid);
%             A(j,i) = A(i,j);
%             
%             % variance 
%             variance(i,j) = var(stats_y.resid);
%             variance(j,i) = var(stats_x.resid);
% 
%         end
%     end
% // 

end