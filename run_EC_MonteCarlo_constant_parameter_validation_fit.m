%% run_EC_MonteCarlo_constant_parameter_validation_fit.m

close all; clear; clc;

%% SECTION 1. USER SETTINGS
% Full search: 100000 samples. Use a small numMC for a trial run.

version = 'AngII_v5';

numMC = 100000;

rngSeed = 1;

topNToSave = 50;

saveOutputTimeCourses = false;
maxTimeCourseRunsToStore = 2000;

validationDirectionThreshold = 0.02;
outputFolder = fullfile('results','constant');

globalNominal.w = 0.81;
globalNominal.n = 1.68;
globalNominal.EC50 = 0.577;

globalSamplingRange.wAbs    = [0.25 1.00];
globalSamplingRange.nAbs    = [1.00 3.00];
globalSamplingRange.EC50Abs = [0.30 0.80];

numInputs = 6;
ActiveInputs = 1:5;
inputValue = 0.2;
stressInputValue = 0.5;
perturbValue = 0.5;

inputOverrides = {
    'Stress', stressInputValue
};

validationPerturbSpecies = {'Stress','Wss','Hypoxia','TNFa'};

validationSpeciesSubset = {'VEGF','VEGFR','IL6','ppET1','ET1','eNOS','FAK','TGFb','pSmad23',...
'pSmad15','p38','ERK','JNK','PKC','DAG','PI3K','Akt','mTOR','p70S6K','RhoA',...
'ROCK','MLCK','Actin','NFkB','HIF1a','MCP1','TSP1','Proliferation','Migration','Apoptosis','EndMT'};

outputSpeciesNames = {
    'NO','ET1','ppET1','eNOS','Klf2','NFkB','HIF1a', ...
    'TGFb','ALK1','ALK5','pSmad15','pSmad23', ...
    'VEGF','VEGFR','IL6','MCP1', ...
    'Proliferation','Migration','Apoptosis','EndMT'};

%% SECTION 2. LOCATE MODEL FILES AND LOAD MODEL

scriptPath = mfilename('fullpath');
if ~isempty(scriptPath)
    modelDir = fileparts(scriptPath);
else
    modelDir = pwd;
end
cd(modelDir);

refFile = ['refcase_EC_Network_', version, '.mat'];
ODEfilename = ['ODElist_EC_Network_', version];
odeFileMat = [ODEfilename, '.mat'];

assert(exist(refFile,'file') == 2, 'Could not find %s in %s.', refFile, modelDir);
assert(exist(odeFileMat,'file') == 2, 'Could not find %s in %s.', odeFileMat, modelDir);
assert(exist('ODE_master.m','file') == 2, 'Could not find ODE_master.m in %s.', modelDir);
assert(exist('CreateValidationMatrix.m','file') == 2, 'Could not find CreateValidationMatrix.m in %s.', modelDir);

load(refFile, 'reactionRules', 'tau', 'ymax', 'speciesNames', 'numEdges');

speciesNames = speciesNames(:);
numSpecies = numel(speciesNames);

if ~exist('reactionRules','var')
    reactionRules = cell(0,1);
else
    reactionRules = reactionRules(:);
end

rng(rngSeed, 'twister');

resultsDir = fullfile(modelDir, outputFolder);
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fprintf('Loaded EC model version %s with %d species and %d rpar edges.\n', version, numSpecies, numEdges);
fprintf('Running CONSTANT-PARAMETER Monte Carlo search.\n');
fprintf('Results will be saved to:\n  %s\n', resultsDir);

%% SECTION 3. VALIDATION TARGETS AND OUTPUT SPECIES

validationPerturbIdx = findSpeciesIndices(speciesNames, validationPerturbSpecies, true);
validationSpeciesIdx = findSpeciesIndices(speciesNames, validationSpeciesSubset, true);
outputIdx = findSpeciesIndices(speciesNames, outputSpeciesNames, true);
outputSpeciesNamesFound = speciesNames(outputIdx);

Exp_validation = CreateValidationMatrix(speciesNames, numel(validationPerturbSpecies));
knownValidationMask = abs(Exp_validation) ~= 0.5;
knownValidationN = nnz(knownValidationMask);

fprintf('Validation matrix contains %d known directional checks.\n', knownValidationN);

SamplingRangeTable = table( ...
    {'w'; 'n'; 'EC50'}, ...
    [globalSamplingRange.wAbs(1); globalSamplingRange.nAbs(1); globalSamplingRange.EC50Abs(1)], ...
    [globalSamplingRange.wAbs(2); globalSamplingRange.nAbs(2); globalSamplingRange.EC50Abs(2)], ...
    [globalNominal.w; globalNominal.n; globalNominal.EC50], ...
    'VariableNames', {'Parameter','Low','High','OriginalNominal'});
writetable(SamplingRangeTable, fullfile(resultsDir, 'constant_parameter_sampling_ranges.csv'));

%% SECTION 4. BUILD INITIAL CONDITION

y0_base = zeros(1, numSpecies);
y0_base(ActiveInputs) = inputValue;
y0_base = applyInputOverrides(y0_base, speciesNames, inputOverrides);

%% SECTION 5. RUN MONTE CARLO SEARCH AND SCORE VALIDATION AGREEMENT

globalParamSamples = nan(numMC, 3);   % columns: W, HillN, EC50

finalSpecies = nan(numMC, numSpecies);
finalOutputs = nan(numMC, numel(outputIdx));

validationScores = nan(numMC,1);
validationMatches = nan(numMC,1);
validationKnownCounts = nan(numMC,1);
validationDiffStore = nan(numMC, numel(validationPerturbSpecies), numSpecies);
validationQualStore = nan(numMC, numel(validationPerturbSpecies), numSpecies);
validationMatchStore = false(numMC, numel(validationPerturbSpecies), numSpecies);

failedRuns = false(numMC,1);
failureMessages = cell(numMC,1);

if saveOutputTimeCourses
    outputTimeCourses = [];
    storedTimeCourseRuns = nan(maxTimeCourseRunsToStore,1);
    timeCourseStoreCount = 0;
    tSim = [];
else
    outputTimeCourses = [];
    storedTimeCourseRuns = [];
    timeCourseStoreCount = 0;
    tSim = [];
end

fprintf('Running %d constant-parameter Monte Carlo validation simulations...\n', numMC);

for runIdx = 1:numMC

    [w_run, n_run, EC50_run] = sampleGlobalConstantParameters(globalSamplingRange);
    globalParamSamples(runIdx,:) = [w_run, n_run, EC50_run];

    rpar = buildConstantRpar(w_run, n_run, EC50_run, numEdges);

    try
        [t, y] = ODE_master(rpar, y0_base, speciesNames, tau, ymax, numEdges, ODEfilename);

        if any(~isfinite(y(:)))
            error('Non-finite values were produced by baseline ODE_master.');
        end

        yBaselineFinal = y(end,:);
        finalSpecies(runIdx,:) = yBaselineFinal;
        finalOutputs(runIdx,:) = yBaselineFinal(outputIdx);

        if saveOutputTimeCourses && timeCourseStoreCount < maxTimeCourseRunsToStore
            if isempty(tSim)
                tSim = t(:);
                outputTimeCourses = nan(numel(tSim), numel(outputIdx), maxTimeCourseRunsToStore);
            end

            if numel(t) == numel(tSim)
                timeCourseStoreCount = timeCourseStoreCount + 1;
                storedTimeCourseRuns(timeCourseStoreCount) = runIdx;
                outputTimeCourses(:,:,timeCourseStoreCount) = y(:,outputIdx);
            else
                warning('Run %d had a different time vector length and was not stored for time-course plotting.', runIdx);
            end
        end

        [scorePct, nMatch, nKnown, Diff_store, Qual_store, Match_store] = scoreExperimentalValidation( ...
            rpar, yBaselineFinal, validationPerturbIdx, perturbValue, Exp_validation, ...
            validationDirectionThreshold, speciesNames, tau, ymax, numEdges, ODEfilename);

        validationScores(runIdx) = scorePct;
        validationMatches(runIdx) = nMatch;
        validationKnownCounts(runIdx) = nKnown;
        validationDiffStore(runIdx,:,:) = Diff_store;
        validationQualStore(runIdx,:,:) = Qual_store;
        validationMatchStore(runIdx,:,:) = Match_store;

    catch ME
        failedRuns(runIdx) = true;
        failureMessages{runIdx} = ME.message;
        warning('Run %d failed: %s', runIdx, ME.message);
    end

    if mod(runIdx, max(1,round(numMC/10))) == 0 || runIdx == numMC
        currentValidScores = validationScores(~failedRuns & isfinite(validationScores));
        if isempty(currentValidScores)
            bestSoFar = NaN;
        else
            bestSoFar = max(currentValidScores);
        end
        fprintf('  Completed %d of %d. Failed: %d. Best score so far: %.2f%%\n', ...
            runIdx, numMC, nnz(failedRuns), bestSoFar);
    end
end

validRuns = find(~failedRuns & isfinite(validationScores));
fprintf('Monte Carlo validation complete. Valid runs: %d / %d.\n', numel(validRuns), numMC);

if isempty(validRuns)
    error('No valid runs were completed. Check failureMessages in the workspace.');
end

[bestScore, localBestIdx] = max(validationScores(validRuns));
bestRun = validRuns(localBestIdx);

bestW = globalParamSamples(bestRun,1);
bestN = globalParamSamples(bestRun,2);
bestEC50 = globalParamSamples(bestRun,3);
bestRpar = buildConstantRpar(bestW, bestN, bestEC50, numEdges);

fprintf('Best run: %d with %.2f%% agreement (%d of %d known checks).\n', ...
    bestRun, bestScore, validationMatches(bestRun), validationKnownCounts(bestRun));
fprintf('Best constant parameters: w = %.6g, n = %.6g, EC50 = %.6g\n', bestW, bestN, bestEC50);

if saveOutputTimeCourses && ~isempty(outputTimeCourses)
    outputTimeCourses = outputTimeCourses(:,:,1:timeCourseStoreCount);
    storedTimeCourseRuns = storedTimeCourseRuns(1:timeCourseStoreCount);
end

%% SECTION 6. SAVE NUMERIC RESULTS AND BEST PARAMETER SETS

runInfo = struct();
runInfo.version = version;
runInfo.numMC = numMC;
runInfo.rngSeed = rngSeed;
runInfo.modelDir = modelDir;
runInfo.resultsDir = resultsDir;
runInfo.globalNominal = globalNominal;
runInfo.globalSamplingRange = globalSamplingRange;
runInfo.ActiveInputs = ActiveInputs;
runInfo.inputValue = inputValue;
runInfo.inputOverrides = inputOverrides;
runInfo.perturbValue = perturbValue;
runInfo.validationPerturbSpecies = validationPerturbSpecies;
runInfo.validationDirectionThreshold = validationDirectionThreshold;
runInfo.validationSpeciesSubset = validationSpeciesSubset;
runInfo.outputSpeciesNames = outputSpeciesNames;
runInfo.failedRuns = failedRuns;
runInfo.failureMessages = failureMessages;
runInfo.bestRun = bestRun;
runInfo.bestScore = bestScore;
runInfo.bestW = bestW;
runInfo.bestN = bestN;
runInfo.bestEC50 = bestEC50;
runInfo.saveOutputTimeCourses = saveOutputTimeCourses;
runInfo.maxTimeCourseRunsToStore = maxTimeCourseRunsToStore;
runInfo.storedTimeCourseRuns = storedTimeCourseRuns;

bestValidationDiff = squeeze(validationDiffStore(bestRun,:,:));
bestValidationQual = squeeze(validationQualStore(bestRun,:,:));
bestValidationMatch = squeeze(validationMatchStore(bestRun,:,:));

save(fullfile(resultsDir, 'MonteCarlo_constant_parameter_validation_fit_all_results.mat'), ...
    'runInfo', 'speciesNames', 'reactionRules', 'Exp_validation', 'knownValidationMask', ...
    'globalParamSamples', 'finalSpecies', 'finalOutputs', 'outputIdx', ...
    'outputSpeciesNamesFound', 'failedRuns', 'failureMessages', 'validationScores', ...
    'validationMatches', 'validationKnownCounts', 'validationDiffStore', ...
    'validationQualStore', 'validationMatchStore', 'bestRun', 'bestScore', ...
    'bestW', 'bestN', 'bestEC50', 'bestRpar', 'bestValidationDiff', ...
    'bestValidationQual', 'bestValidationMatch', 'tSim', 'outputTimeCourses', ...
    'storedTimeCourseRuns', '-v7.3');

GlobalParameterTable = table((1:numMC)', failedRuns, globalParamSamples(:,1), globalParamSamples(:,2), globalParamSamples(:,3), ...
    validationScores, validationMatches, validationKnownCounts, ...
    'VariableNames', {'Run','Failed','W','HillN','EC50','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks'});
writetable(GlobalParameterTable, fullfile(resultsDir, 'sampled_global_constant_parameters.csv'));

FinalSpeciesTable = [table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    globalParamSamples(:,1), globalParamSamples(:,2), globalParamSamples(:,3), ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks','W','HillN','EC50'}), ...
    array2table(finalSpecies, 'VariableNames', matlab.lang.makeValidName(speciesNames))];
writetable(FinalSpeciesTable, fullfile(resultsDir, 'final_species_values_with_validation_scores.csv'));

FinalOutputTable = [table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    globalParamSamples(:,1), globalParamSamples(:,2), globalParamSamples(:,3), ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks','W','HillN','EC50'}), ...
    array2table(finalOutputs, 'VariableNames', matlab.lang.makeValidName(outputSpeciesNamesFound))];
writetable(FinalOutputTable, fullfile(resultsDir, 'final_output_species_values_with_validation_scores.csv'));

SummaryStats = summarizeFinalOutputs(finalOutputs(validRuns,:), outputSpeciesNamesFound);
writetable(SummaryStats, fullfile(resultsDir, 'summary_statistics_outputs.csv'));

ValidationScoreTable = table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    globalParamSamples(:,1), globalParamSamples(:,2), globalParamSamples(:,3), ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks','W','HillN','EC50'});
writetable(ValidationScoreTable, fullfile(resultsDir, 'validation_scores_by_run.csv'));

BestParameterTable = table(bestRun, bestScore, validationMatches(bestRun), validationKnownCounts(bestRun), bestW, bestN, bestEC50, ...
    'VariableNames', {'Run','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks','W','HillN','EC50'});
writetable(BestParameterTable, fullfile(resultsDir, 'best_constant_parameter_set.csv'));
writematrix(bestRpar, fullfile(resultsDir, 'best_rpar_matrix_rows_w_n_EC50.csv'));

scoreForSort = validationScores;
scoreForSort(failedRuns | ~isfinite(scoreForSort)) = -Inf;
[~, sortedIdx] = sort(scoreForSort, 'descend');
sortedIdx = sortedIdx(scoreForSort(sortedIdx) > -Inf);
topRuns = sortedIdx(1:min(topNToSave, numel(sortedIdx)));
TopScoreTable = ValidationScoreTable(topRuns,:);
writetable(TopScoreTable, fullfile(resultsDir, 'top_validation_scores.csv'));

BestValidationTable = buildValidationResultTable(bestValidationDiff, bestValidationQual, bestValidationMatch, ...
    Exp_validation, knownValidationMask, validationPerturbSpecies, speciesNames);
writetable(BestValidationTable, fullfile(resultsDir, 'best_validation_direction_matrix_long.csv'));

fprintf('Saved constant-parameter validation fit results and best parameter set.\n');

%% SECTION 7. PLOTS

makeSummaryPlots(resultsDir, finalOutputs, outputSpeciesNamesFound, validRuns, tSim, outputTimeCourses, storedTimeCourseRuns);
makeValidationPlots(resultsDir, validationScores, failedRuns, bestRun, bestScore, ...
    bestValidationDiff, bestValidationQual, Exp_validation, knownValidationMask, ...
    validationPerturbSpecies, speciesNames, validationSpeciesIdx, validationSpeciesSubset);

fprintf('Saved plots in %s.\n', resultsDir);
fprintf('Done. Best constant-parameter rpar is saved in best_rpar_matrix_rows_w_n_EC50.csv.\n');

%% LOCAL FUNCTIONS

function [w, n, EC50] = sampleGlobalConstantParameters(globalSamplingRange)

w = uniformSample(globalSamplingRange.wAbs(1), globalSamplingRange.wAbs(2));
n = uniformSample(globalSamplingRange.nAbs(1), globalSamplingRange.nAbs(2));

validEC50Upper = 0.98 * (0.5 .^ (1 ./ n));
ec50Low = globalSamplingRange.EC50Abs(1);
ec50High = min(globalSamplingRange.EC50Abs(2), validEC50Upper);

if ec50High <= ec50Low
    ec50High = validEC50Upper;
    ec50Low = min(ec50Low, 0.90 * validEC50Upper);
    ec50Low = max(0.01, ec50Low);
end

EC50 = uniformSample(ec50Low, ec50High);
end

function rpar = buildConstantRpar(w, n, EC50, numEdges)
rpar = zeros(3, numEdges);
rpar(1,:) = w;
rpar(2,:) = n;
rpar(3,:) = EC50;
rpar(1,6) = 0;

end

function x = uniformSample(lo, hi)
if hi < lo
    tmp = lo; lo = hi; hi = tmp;
end
if abs(hi - lo) < eps
    x = lo;
else
    x = lo + (hi - lo) .* rand();
end
end

function y0 = applyInputOverrides(y0, speciesNames, inputOverrides)
if isempty(inputOverrides)
    return;
end
for i = 1:size(inputOverrides,1)
    species = inputOverrides{i,1};
    value = inputOverrides{i,2};
    idx = find(strcmp(speciesNames, species), 1);
    if isempty(idx)
        warning('Input override species not found: %s', species);
    else
        y0(idx) = value;
    end
end
end

function idx = findSpeciesIndices(speciesNames, requestedSpecies, warnIfMissing)
idx = [];
for i = 1:numel(requestedSpecies)
    thisName = requestedSpecies{i};
    thisIdx = find(strcmp(speciesNames, thisName), 1);
    if isempty(thisIdx)
        if warnIfMissing
            warning('Requested species not found: %s', thisName);
        end
    else
        idx(end+1) = thisIdx; %#ok<AGROW>
    end
end
end

function [scorePct, nMatch, nKnown, Diff_store, Qual_store, Match_store] = scoreExperimentalValidation( ...
    rpar, yBaselineFinal, validationPerturbIdx, perturbValue, Exp_validation, validationDirectionThreshold, ...
    speciesNames, tau, ymax, numEdges, ODEfilename)

numPert = numel(validationPerturbIdx);
numSpecies = numel(speciesNames);
Diff_store = nan(numPert, numSpecies);

for p = 1:numPert
    y0_pert = yBaselineFinal;
    y0_pert(validationPerturbIdx(p)) = y0_pert(validationPerturbIdx(p)) + perturbValue;

    [~, y_perturbed] = ODE_master(rpar, y0_pert, speciesNames, tau, ymax, numEdges, ODEfilename);

    if any(~isfinite(y_perturbed(:)))
        error('Non-finite values were produced by perturbed ODE_master.');
    end

    Diff_store(p,:) = y_perturbed(end,:) - yBaselineFinal;
end

Diff_for_sign = Diff_store;
Diff_for_sign(abs(Diff_for_sign) < validationDirectionThreshold) = 0;
Qual_store = sign(Diff_for_sign);

knownMask = abs(Exp_validation) ~= 0.5;
Match_store = (Qual_store == Exp_validation) & knownMask;

nMatch = nnz(Match_store);
nKnown = nnz(knownMask);
scorePct = 100 * nMatch / nKnown;
end

function SummaryStats = summarizeFinalOutputs(finalOutputsValid, outputSpeciesNamesFound)
numOutputs = numel(outputSpeciesNamesFound);
Mean = nan(numOutputs,1);
Median = nan(numOutputs,1);
Std = nan(numOutputs,1);
P05 = nan(numOutputs,1);
P25 = nan(numOutputs,1);
P75 = nan(numOutputs,1);
P95 = nan(numOutputs,1);

for j = 1:numOutputs
    x = finalOutputsValid(:,j);
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    Mean(j) = mean(x);
    Median(j) = median(x);
    Std(j) = std(x);
    pct = prctile(x, [5 25 75 95]);
    P05(j) = pct(1);
    P25(j) = pct(2);
    P75(j) = pct(3);
    P95(j) = pct(4);
end

Species = outputSpeciesNamesFound(:);
SummaryStats = table(Species, Mean, Median, Std, P05, P25, P75, P95);
end

function T = buildValidationResultTable(Diff_store, Qual_store, Match_store, Exp_validation, knownMask, validationPerturbSpecies, speciesNames)
numPert = size(Diff_store,1);
numSpecies = size(Diff_store,2);
N = numPert * numSpecies;
Perturbation = cell(N,1);
Species = cell(N,1);
ExpectedDirection = nan(N,1);
ObservedDifference = nan(N,1);
ObservedDirection = nan(N,1);
KnownLiteratureCheck = false(N,1);
Matched = false(N,1);

row = 0;
for p = 1:numPert
    for j = 1:numSpecies
        row = row + 1;
        Perturbation{row} = validationPerturbSpecies{p};
        Species{row} = speciesNames{j};
        ExpectedDirection(row) = Exp_validation(p,j);
        ObservedDifference(row) = Diff_store(p,j);
        ObservedDirection(row) = Qual_store(p,j);
        KnownLiteratureCheck(row) = knownMask(p,j);
        Matched(row) = Match_store(p,j);
    end
end

T = table(Perturbation, Species, ExpectedDirection, ObservedDifference, ObservedDirection, KnownLiteratureCheck, Matched);
end

function makeSummaryPlots(resultsDir, finalOutputs, outputSpeciesNamesFound, validRuns, tSim, outputTimeCourses, storedTimeCourseRuns)
if isempty(validRuns)
    warning('No valid runs available for plotting.');
    return;
end

plotData = finalOutputs(validRuns,:);
numOutputs = size(plotData,2);

fig1 = figure('Name','Monte Carlo final output distributions','Color','w');
nRows = ceil(sqrt(numOutputs));
nCols = ceil(numOutputs / nRows);
for j = 1:numOutputs
    subplot(nRows, nCols, j);
    histogram(plotData(:,j));
    title(outputSpeciesNamesFound{j}, 'Interpreter','none');
    xlabel('Final value');
    ylabel('Count');
end
set(fig1, 'Position', [100 100 1200 800]);
saveas(fig1, fullfile(resultsDir, 'final_output_histograms.png'));
savefig(fig1, fullfile(resultsDir, 'final_output_histograms.fig'));

fig2 = figure('Name','Monte Carlo final output heatmap','Color','w');
imagesc(plotData);
colorbar;
xlabel('Output species');
ylabel('Valid Monte Carlo run');
xticks(1:numOutputs);
xticklabels(outputSpeciesNamesFound);
xtickangle(45);
title('Baseline final output values across constant-parameter Monte Carlo runs');
set(fig2, 'Position', [100 100 1200 500]);
saveas(fig2, fullfile(resultsDir, 'final_output_heatmap.png'));
savefig(fig2, fullfile(resultsDir, 'final_output_heatmap.fig'));

if ~isempty(tSim) && ~isempty(outputTimeCourses)
    validStored = storedTimeCourseRuns(isfinite(storedTimeCourseRuns));
    if isempty(validStored)
        return;
    end

    fig3 = figure('Name','Monte Carlo output time course envelopes','Color','w');
    nRows = ceil(sqrt(numOutputs));
    nCols = ceil(numOutputs / nRows);

    for j = 1:numOutputs
        subplot(nRows, nCols, j);
        y = squeeze(outputTimeCourses(:, j, :));
        med = median(y, 2, 'omitnan');
        lo = prctile(y, 5, 2);
        hi = prctile(y, 95, 2);

        fill([tSim; flipud(tSim)], [lo; flipud(hi)], [0.8 0.8 0.8], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5);
        hold on;
        plot(tSim, med, 'LineWidth', 1.5);
        title(outputSpeciesNamesFound{j}, 'Interpreter','none');
        xlabel('Time');
        ylabel('Model value');
    end
    set(fig3, 'Position', [100 100 1200 800]);
    saveas(fig3, fullfile(resultsDir, 'output_timecourse_envelopes.png'));
    savefig(fig3, fullfile(resultsDir, 'output_timecourse_envelopes.fig'));
end
end

function makeValidationPlots(resultsDir, validationScores, failedRuns, bestRun, bestScore, ...
    bestValidationDiff, bestValidationQual, Exp_validation, knownValidationMask, ...
    validationPerturbSpecies, speciesNames, validationSpeciesIdx, validationSpeciesSubset)

valid = find(~failedRuns & isfinite(validationScores));

fig1 = figure('Name','Validation agreement score distribution','Color','w');
histogram(validationScores(valid));
xlabel('Agreement with literature validation matrix (%)');
ylabel('Count');
title(sprintf('Validation score distribution; best run = %d, best score = %.2f%%', bestRun, bestScore));
set(fig1, 'Position', [100 100 800 500]);
saveas(fig1, fullfile(resultsDir, 'validation_score_histogram.png'));
savefig(fig1, fullfile(resultsDir, 'validation_score_histogram.fig'));

fig2 = figure('Name','Validation score by run','Color','w');
plot(valid, validationScores(valid), '.', 'MarkerSize', 10);
hold on;
plot(bestRun, bestScore, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Monte Carlo run');
ylabel('Agreement with literature validation matrix (%)');
title('Validation agreement by sampled constant parameter set');
set(fig2, 'Position', [100 100 900 450]);
saveas(fig2, fullfile(resultsDir, 'validation_score_by_run.png'));
savefig(fig2, fullfile(resultsDir, 'validation_score_by_run.fig'));

plotIdx = validationSpeciesIdx;
plotNames = speciesNames(plotIdx);
numPlotSpecies = numel(plotIdx);
numPert = numel(validationPerturbSpecies);

fig3 = figure('Name','Best validation direction matrix','Color','w');
imagesc(bestValidationDiff(:,plotIdx));
xticks(1:numPlotSpecies);
xticklabels(plotNames);
xtickangle(270);
yticks(1:numPert);
yticklabels(validationPerturbSpecies);
caxis([-0.5, 1]);
if exist('custom_cmap.m','file') == 2
    colormap custom_cmap;
else
    colormap parula;
end
colorbar;
set(gca,'TickLength',[0 0]);
set(gca,'Fontsize',8);
title(sprintf('Best constant-parameter run %d: %.2f%% agreement with validation matrix', bestRun, bestScore));
hold on;

for i = 1:numPlotSpecies
   plot([i-.5,i-.5],[.5,numPert+.5],'k-');
end
for i = 1:numPert
   plot([.5,numPlotSpecies+.5],[i+.5,i+.5],'k-');
end

for p = 1:numPert
    for jj = 1:numPlotSpecies
        speciesIdx = plotIdx(jj);
        if knownValidationMask(p,speciesIdx)
            if bestValidationQual(p,speciesIdx) == Exp_validation(p,speciesIdx)
                text(jj-0.25, p, char(hex2dec('2713')), 'FontWeight','bold');
            else
                text(jj-0.25, p, char(hex2dec('2717')), 'FontWeight','bold');
            end
        end
    end
end

set(fig3, 'Position', [100 100 1200 350]);
saveas(fig3, fullfile(resultsDir, 'best_validation_matrix_checks_x.png'));
savefig(fig3, fullfile(resultsDir, 'best_validation_matrix_checks_x.fig'));
end
