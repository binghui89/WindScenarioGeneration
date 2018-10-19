function  [copModel] = func_choose_best_copula_model(U)

% X1 = up_ramps(:,3);
% X2 = up_ramps(:,8);
% X3 = data_objTransformed;

copulaModelName = {'gaussian', 't', 'clayton', 'gumbel', 'frank'};

% % % copModels = [];

% initialize statistical metrics
ll = zeros(1, length(copulaModelName));
aic = zeros(1, length(copulaModelName));
bic = zeros(1, length(copulaModelName));
aks = zeros(1, length(copulaModelName));
snc = zeros(1, length(copulaModelName));
copModels = cell(1,5);

for j = 1 : numel(copulaModelName)
    copulaparams = copula.fit(copulaModelName{j}, U);
    dbg('fitcopulas', 3, 'Computing statistics.\n');    
    [ll(1,j), aic(1,j), bic(1,j), aks(1,j), snc(1,j)] = copula.fitStatistics(copulaparams, U);
    % copModels(j).data_obj = copulaparams;
    copModels{j} = copulaparams;
end

% Sort fitted copula models using ll, aic, and bic
% detailed information about ll, aic, and bic, see the function 'copula.fitStatistics'
% largest ll, smallest aic and bic
[~, ll_max_idx] = max(ll);
[~, aic_min_idx] = min(aic);
[~, bic_min_idx] = min(bic);

temp = tabulate([ll_max_idx, aic_min_idx, bic_min_idx]);
[~, idxxxxx] = max(temp(:,3));

if isnumeric(idxxxxx)
    copModel = copModels{idxxxxx};
    copModel.aic = aic(idxxxxx);
    copModel.bic = bic(idxxxxx);
    copModel.ll = ll(idxxxxx);
else 
    error('the best copula model is not chosen!!!')   
end


figure
plot(ll); hold on
plot(aic); hold on
plot(bic); hold on
close















    