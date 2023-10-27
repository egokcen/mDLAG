function plotCrossValPerf_gfa(res)
%
% plotCrossValPerf_gfa(res)
%
% Description: Plot cross-validated performance for the cross-validated
%              results in res.
% 
% Arguments:
%
%     res -- structure whose i-th entry (corresponding to the i-th latent
%            dimensionality) has fields
%              R       -- rank of factorization
%              R2      -- average cross-validated leave-group-out R^2
%              R2_sem  -- standard error of R2 across CV folds
%              MSE     -- average cross-validated leave-group-out 
%                         mean-squared error
%              MSE_sem -- standard error of MSE across CV folds
%              model   -- GFA model estimated using all data
%
% Outputs:
%     None.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     02 Jul 2022 -- Initial full revision.
    
RList = [res.R];

% Plot R^2 versus model rank
R2 = [res.R2];
R2_sem = [res.R2_sem];
[best_performance, maxIdx] = max(R2);
best_sem = res(maxIdx).R2_sem;
% NOTE: max used in this manner selects the first instance that
%       satisfies the input condition.
[~, bestModel] = max(R2 >= best_performance - best_sem);

figure;
subplot(1,2,1);
hold on;
errorbar(RList, R2, R2_sem, 'o-', 'MarkerFaceColor', 'b');
plot(RList(bestModel), R2(bestModel), 'g*', 'markersize', 5);
plot(RList, ones(1,length(RList)).*best_performance,'k--');
plot(RList, ones(1,length(RList)).*best_performance + best_sem,'r--');
plot(RList, ones(1,length(RList)).*best_performance - best_sem,'r--');
xlabel('Model rank');
ylabel('Cross-validated R^2');
hold off;

% Plot MSE versus model rank
MSE = [res.MSE];
MSE_sem = [res.MSE_sem];
[best_performance, minIdx] = min(MSE);
best_sem = res(minIdx).MSE_sem;
% NOTE: min used in this manner selects the first instance that
%       satisfies the input condition.
[~, bestModel] = max(MSE <= best_performance + best_sem); 

subplot(1,2,2);
hold on;
errorbar(RList, MSE, MSE_sem, 'o-', 'MarkerFaceColor', 'b');
plot(RList(bestModel), MSE(bestModel), 'g*', 'markersize', 5);
plot(RList, ones(1,length(RList)).*best_performance,'k--');
plot(RList, ones(1,length(RList)).*best_performance + best_sem,'r--');
plot(RList, ones(1,length(RList)).*best_performance - best_sem,'r--');
xlabel('Model rank');
ylabel('Cross-validated MSE');
hold off;
