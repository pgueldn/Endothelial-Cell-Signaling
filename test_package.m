function test_package
% Check archive consistency and execute one baseline ODE integration.
root = fileparts(mfilename('fullpath'));
original = pwd;
cleanup = onCleanup(@() cd(original));
cd(root);
M = load('refcase_EC_Network_AngII_v5.mat');
E = CreateValidationMatrix(M.speciesNames,4);
assert(nnz(abs(E)~=0.5)==79);
folders = {'constant','reaction_specific'};
expected = [43.5846202531646 18.428696331097;68.3385569620254 10.8222335857063];
for k=1:2
    T = readtable(fullfile('data',folders{k},'validation_scores_by_run.csv'));
    values = T.ValidationAgreementPct(~logical(T.Failed));
    assert(numel(values)==100000 && all(isfinite(values)));
    assert(abs(mean(values)-expected(k,1))<1e-9);
    assert(abs(std(values)-expected(k,2))<1e-9);
end
W = readtable(fullfile('data','figures','Figure_5_wss_values.csv'));
assert(height(W)==121 && W.Wss(1)==0 && W.Wss(end)==1);
rpar = readmatrix(fullfile('data','constant','best_rpar_matrix_rows_w_n_EC50.csv'));
assert(isequal(size(rpar),[3 M.numEdges]));
y0 = zeros(1,numel(M.speciesNames));
y0(1:5)=0.2;
y0(1)=0.5;
[~,y] = ODE_master(rpar,y0,M.speciesNames,M.tau,M.ymax,M.numEdges,...
    fullfile(root,'ODElist_EC_Network_AngII_v5.mat'));
assert(all(isfinite(y(:))));
fprintf('Package data checks and baseline ODE integration passed.\n');
end
