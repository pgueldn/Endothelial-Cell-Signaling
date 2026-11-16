%% run_EC_MonteCarlo_parameter_sweep_validation_fit.m

close all; clear; clc;

%% SECTION 1. USER SETTINGS
% Full search: 100000 samples. Baseline trajectories are optional and large.

version = 'AngII_v5';

numMC = 100000;

rngSeed = 1;

topNToSave = 50;

saveFullTimeCourses = false;

validationDirectionThreshold = 0.02;
outputFolder = fullfile('results','reaction_specific');

globalNominal.w = 0.81;
globalNominal.n = 1.68;
globalNominal.EC50 = 0.577;

defaultRange.wFactor = [0.30 1.50];      % multiplier around nominal w
defaultRange.nAbs    = [1.00 3.00];      % absolute Hill coefficient range
defaultRange.EC50Abs = [0.30 0.80];      % absolute EC50 range, clipped if needed

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

manualRangeSpec = {
    'Stress',  'Integrins', 'act',   0.90, 0.65, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'mechanosensitive integrin activation; tune to stress validation';
    'Stress',  'latentTGFb','act',   0.90, 0.60, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'stress-mediated latent TGFb activation; tune to stress validation';
    'Wss',     'IonChannels','act',  0.90, 0.60, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'rapid shear-sensitive ion channel activation; tune to Wss validation';
    'Wss',     'ETBR',      'act',   0.80, 0.50, 1.00, 1.60, 1.10, 2.80, 0.55, 0.40, 0.60, 'shear regulation of ETBR; tune to Wss validation';
    'Wss',     'Klf2',      'act',   0.95, 0.65, 1.00, 2.00, 1.10, 3.00, 0.55, 0.40, 0.60, 'flow/Klf2 response; tune to Wss validation';
    'TNFa',    'TNFaR',     'act',   0.95, 0.65, 1.00, 1.70, 1.10, 3.00, 0.55, 0.40, 0.60, 'TNFa receptor activation; tune to TNFa validation';
    'TNFaR',   'NFkB',      'act',   0.95, 0.65, 1.00, 2.00, 1.10, 3.00, 0.55, 0.40, 0.60, 'TNFaR/NFkB inflammatory signaling; tune to TNFa validation';
    'TNFaR',   'p38',       'act',   0.85, 0.50, 1.00, 1.60, 1.10, 3.00, 0.55, 0.40, 0.60, 'TNFaR/p38 inflammatory signaling; tune to TNFa validation';
    'TNFaR',   'JNK',       'act',   0.85, 0.50, 1.00, 1.60, 1.10, 3.00, 0.55, 0.40, 0.60, 'TNFaR/JNK inflammatory signaling; tune to TNFa validation';
    'TGFb',    'ALK1',      'act',   0.80, 0.40, 1.00, 1.60, 1.10, 3.00, 0.55, 0.40, 0.60, 'TGFb/ALK1 arm; tune to Smad1/5 behavior';
    'TGFb',    'ALK5',      'act',   0.90, 0.50, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'TGFb/ALK5 arm; tune to Smad2/3 behavior';
    'NO',      'ppET1',     'inhib', 0.95, 0.60, 1.00, 2.00, 1.10, 3.00, 0.55, 0.40, 0.60, 'NO inhibition of ppET1-related rules; tune NO/ET1 balance';
    'ET1',     'ETBR',      'act',   0.85, 0.50, 1.00, 1.60, 1.10, 3.00, 0.55, 0.40, 0.60, 'ET1/ETBR signaling';
    'ETBR',    'PI3K',      'act',   0.85, 0.50, 1.00, 1.60, 1.10, 3.00, 0.55, 0.40, 0.60, 'ETBR/PI3K signaling';
    'Hypoxia', 'HIF1a',     'act',   0.95, 0.65, 1.00, 1.90, 1.10, 3.00, 0.55, 0.40, 0.60, 'hypoxia stabilization of HIF1a; tune to hypoxia validation';
    'HIF1a',   'ppET1',     'act',   0.90, 0.50, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'HIF1a regulation of ppET1';
    'HIF1a',   'VEGF',      'act',   0.95, 0.60, 1.00, 2.00, 1.10, 3.00, 0.55, 0.40, 0.60, 'HIF1a/VEGF response';
    'Klf2',    'NFkB',      'inhib', 0.95, 0.60, 1.00, 2.00, 1.10, 3.00, 0.55, 0.40, 0.60, 'Klf2 inhibition of NFkB; tune Wss anti-inflammatory behavior';
    'Klf2',    'pSmad23',   'inhib', 0.90, 0.50, 1.00, 1.80, 1.10, 3.00, 0.55, 0.40, 0.60, 'Klf2 inhibition of pSmad23; tune Wss/EndMT behavior'
};

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
fprintf('Results will be saved to:\n  %s\n', resultsDir);

%% SECTION 3. INFER rpar EDGE LIST, PARAMETER RANGES, AND VALIDATION TARGETS

EdgeTable = inferEdgeTableFromODElist(odeFileMat, speciesNames, reactionRules, numEdges);
RangeTable = createDefaultRangeTable(EdgeTable, globalNominal, defaultRange);
RangeTable = applyManualRangeSpec(RangeTable, manualRangeSpec);

edgeTableFile = fullfile(resultsDir, 'inferred_rpar_edge_table.csv');
rangeTableFile = fullfile(resultsDir, 'rpar_sampling_ranges.csv');
writetable(EdgeTable, edgeTableFile);
writetable(RangeTable, rangeTableFile);

validationPerturbIdx = findSpeciesIndices(speciesNames, validationPerturbSpecies, true);
validationSpeciesIdx = findSpeciesIndices(speciesNames, validationSpeciesSubset, true);
outputIdx = findSpeciesIndices(speciesNames, outputSpeciesNames, true);
outputSpeciesNamesFound = speciesNames(outputIdx);

Exp_validation = CreateValidationMatrix(speciesNames, numel(validationPerturbSpecies));
knownValidationMask = abs(Exp_validation) ~= 0.5;
knownValidationN = nnz(knownValidationMask);

fprintf('Validation matrix contains %d known directional checks.\n', knownValidationN);

%% SECTION 4. BUILD INITIAL CONDITION

y0_base = zeros(1, numSpecies);
y0_base(ActiveInputs) = inputValue;
y0_base = applyInputOverrides(y0_base, speciesNames, inputOverrides);

%% SECTION 5. RUN MONTE CARLO SEARCH AND SCORE VALIDATION AGREEMENT

paramSamples = nan(3, numEdges, numMC);
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

if saveFullTimeCourses
    fullTimeCourses = [];
    tSim = [];
else
    fullTimeCourses = [];
    tSim = [];
end

fprintf('Running %d Monte Carlo validation simulations...\n', numMC);

for runIdx = 1:numMC
    rpar = sampleRparFromRanges(RangeTable);
    paramSamples(:,:,runIdx) = rpar;

    try
        [t, y] = ODE_master(rpar, y0_base, speciesNames, tau, ymax, numEdges, ODEfilename);

        if any(~isfinite(y(:)))
            error('Non-finite values were produced by baseline ODE_master.');
        end

        yBaselineFinal = y(end,:);
        finalSpecies(runIdx,:) = yBaselineFinal;
        finalOutputs(runIdx,:) = yBaselineFinal(outputIdx);

        if saveFullTimeCourses
            if isempty(tSim)
                tSim = t(:);
                fullTimeCourses = nan(numel(tSim), numSpecies, numMC);
            end
            fullTimeCourses(:,:,runIdx) = y;
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
fprintf('Best run: %d with %.2f%% agreement (%d of %d known checks).\n', ...
    bestRun, bestScore, validationMatches(bestRun), validationKnownCounts(bestRun));

%% SECTION 6. SAVE NUMERIC RESULTS AND BEST PARAMETER SETS

runInfo = struct();
runInfo.version = version;
runInfo.numMC = numMC;
runInfo.rngSeed = rngSeed;
runInfo.modelDir = modelDir;
runInfo.resultsDir = resultsDir;
runInfo.globalNominal = globalNominal;
runInfo.defaultRange = defaultRange;
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
runInfo.saveFullTimeCourses = saveFullTimeCourses;

bestRpar = paramSamples(:,:,bestRun);
bestValidationDiff = squeeze(validationDiffStore(bestRun,:,:));
bestValidationQual = squeeze(validationQualStore(bestRun,:,:));
bestValidationMatch = squeeze(validationMatchStore(bestRun,:,:));

save(fullfile(resultsDir, 'MonteCarlo_validation_fit_all_results.mat'), ...
    'runInfo', 'speciesNames', 'reactionRules', 'EdgeTable', 'RangeTable', ...
    'Exp_validation', 'knownValidationMask', 'paramSamples', 'finalSpecies', ...
    'finalOutputs', 'outputIdx', 'outputSpeciesNamesFound', 'failedRuns', ...
    'failureMessages', 'validationScores', 'validationMatches', ...
    'validationKnownCounts', 'validationDiffStore', 'validationQualStore', ...
    'validationMatchStore', 'bestRun', 'bestScore', 'bestRpar', ...
    'bestValidationDiff', 'bestValidationQual', 'bestValidationMatch', ...
    'tSim', 'fullTimeCourses', '-v7.3');

FinalSpeciesTable = [table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks'}), ...
    array2table(finalSpecies, 'VariableNames', matlab.lang.makeValidName(speciesNames))];
writetable(FinalSpeciesTable, fullfile(resultsDir, 'final_species_values_with_validation_scores.csv'));

FinalOutputTable = [table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks'}), ...
    array2table(finalOutputs, 'VariableNames', matlab.lang.makeValidName(outputSpeciesNamesFound))];
writetable(FinalOutputTable, fullfile(resultsDir, 'final_output_species_values_with_validation_scores.csv'));

SummaryStats = summarizeFinalOutputs(finalOutputs(validRuns,:), outputSpeciesNamesFound);
writetable(SummaryStats, fullfile(resultsDir, 'summary_statistics_outputs.csv'));

ParameterLongTable = buildLongParameterTable(paramSamples, EdgeTable);
writetable(ParameterLongTable, fullfile(resultsDir, 'sampled_parameters_long.csv'));

ValidationScoreTable = table((1:numMC)', failedRuns, validationScores, validationMatches, validationKnownCounts, ...
    'VariableNames', {'Run','Failed','ValidationAgreementPct','ValidationMatches','ValidationKnownChecks'});
writetable(ValidationScoreTable, fullfile(resultsDir, 'validation_scores_by_run.csv'));

BestParameterTable = buildSingleRunParameterTable(bestRpar, EdgeTable, bestRun, bestScore);
writetable(BestParameterTable, fullfile(resultsDir, 'best_rpar_parameter_set.csv'));
writematrix(bestRpar, fullfile(resultsDir, 'best_rpar_matrix_rows_w_n_EC50.csv'));

scoreForSort = validationScores;
scoreForSort(failedRuns | ~isfinite(scoreForSort)) = -Inf;
[~, sortedIdx] = sort(scoreForSort, 'descend');
sortedIdx = sortedIdx(scoreForSort(sortedIdx) > -Inf);
topRuns = sortedIdx(1:min(topNToSave, numel(sortedIdx)));
TopScoreTable = ValidationScoreTable(topRuns,:);
writetable(TopScoreTable, fullfile(resultsDir, 'top_validation_scores.csv'));

TopParameterLongTable = buildLongParameterTable(paramSamples(:,:,topRuns), EdgeTable, topRuns);
writetable(TopParameterLongTable, fullfile(resultsDir, 'top_parameter_sets_long.csv'));

BestValidationTable = buildValidationResultTable(bestValidationDiff, bestValidationQual, bestValidationMatch, ...
    Exp_validation, knownValidationMask, validationPerturbSpecies, speciesNames);
writetable(BestValidationTable, fullfile(resultsDir, 'best_validation_direction_matrix_long.csv'));

fprintf('Saved validation fit results and best parameter sets.\n');

%% SECTION 7. PLOTS

makeSummaryPlots(resultsDir, finalOutputs, outputSpeciesNamesFound, validRuns, tSim, fullTimeCourses, outputIdx);
makeValidationPlots(resultsDir, validationScores, failedRuns, bestRun, bestScore, ...
    bestValidationDiff, bestValidationQual, Exp_validation, knownValidationMask, ...
    validationPerturbSpecies, speciesNames, validationSpeciesIdx, validationSpeciesSubset);

fprintf('Saved plots in %s.\n', resultsDir);
fprintf('Done. Best sampled rpar is saved in best_rpar_parameter_set.csv and best_rpar_matrix_rows_w_n_EC50.csv.\n');

%% LOCAL FUNCTIONS

function EdgeTable = inferEdgeTableFromODElist(odeFileMat, speciesNames, reactionRules, numEdges)
S = load(odeFileMat, 'ODElist');
ODElist = S.ODElist(:);

EdgeIndex = (1:numEdges)';
Source = repmat({''}, numEdges, 1);
Target = repmat({''}, numEdges, 1);
Mode = repmat({''}, numEdges, 1);
ODERow = nan(numEdges,1);
ODEExpression = repmat({''}, numEdges, 1);
RuleGuess = repmat({''}, numEdges, 1);
UsedInODE = false(numEdges,1);

for lineIdx = 1:numel(ODElist)
    line = char(ODElist{lineIdx});

    targetTok = regexp(line, 'dydt\((\w+)\)', 'tokens', 'once');
    if isempty(targetTok)
        targetName = '';
    else
        targetName = targetTok{1};
    end

    actToks = regexp(line, 'act\(y\((\w+)\),rpar\(:,(\d+)\)\)', 'tokens');
    for j = 1:numel(actToks)
        src = actToks{j}{1};
        k = str2double(actToks{j}{2});
        [Source, Target, Mode, ODERow, ODEExpression, UsedInODE] = recordEdge( ...
            k, src, targetName, 'act', lineIdx, line, Source, Target, Mode, ODERow, ODEExpression, UsedInODE, numEdges);
    end

    inhibToks = regexp(line, 'inhib\(y\((\w+)\),rpar\(:,(\d+)\)\)', 'tokens');
    for j = 1:numel(inhibToks)
        src = inhibToks{j}{1};
        k = str2double(inhibToks{j}{2});
        [Source, Target, Mode, ODERow, ODEExpression, UsedInODE] = recordEdge( ...
            k, src, targetName, 'inhib', lineIdx, line, Source, Target, Mode, ODERow, ODEExpression, UsedInODE, numEdges);
    end

    andToks = regexp(line, 'AND\(rpar\(:,(\d+)\)', 'tokens');
    for j = 1:numel(andToks)
        k = str2double(andToks{j}{1});
        [Source, Target, Mode, ODERow, ODEExpression, UsedInODE] = recordEdge( ...
            k, 'logic_gate', targetName, 'AND_gate', lineIdx, line, Source, Target, Mode, ODERow, ODEExpression, UsedInODE, numEdges);
    end

    allRparToks = regexp(line, 'rpar\(:,(\d+)\)', 'tokens');
    for j = 1:numel(allRparToks)
        k = str2double(allRparToks{j}{1});
        if k >= 1 && k <= numEdges && ~UsedInODE(k)
            [Source, Target, Mode, ODERow, ODEExpression, UsedInODE] = recordEdge( ...
                k, 'unknown', targetName, 'unknown', lineIdx, line, Source, Target, Mode, ODERow, ODEExpression, UsedInODE, numEdges);
        end
    end
end

for k = 1:min(numEdges, numel(reactionRules))
    if isempty(Source{k}) && isempty(Target{k})
        rr = char(reactionRules{k});
        RuleGuess{k} = rr;
        parts = regexp(rr, '=>', 'split');
        if numel(parts) == 2
            lhs = strtrim(parts{1});
            rhs = strtrim(parts{2});
            if isempty(lhs)
                lhs = 'input';
            end
            Source{k} = lhs;
            Target{k} = rhs;
            Mode{k} = 'input_or_grouped_rule';
        end
    end
end

for k = 1:numEdges
    if isempty(RuleGuess{k})
        if ~isempty(Source{k}) && ~isempty(Target{k})
            if strcmpi(Mode{k}, 'inhib') || contains(lower(Mode{k}), 'inhib')
                RuleGuess{k} = ['!', Source{k}, '=>', Target{k}];
            elseif strcmpi(Source{k}, 'logic_gate')
                RuleGuess{k} = ['AND_gate=>', Target{k}];
            else
                RuleGuess{k} = [Source{k}, '=>', Target{k}];
            end
        else
            RuleGuess{k} = sprintf('rpar(:,%d)', k);
        end
    end
end

EdgeTable = table(EdgeIndex, Source, Target, Mode, UsedInODE, ODERow, RuleGuess, ODEExpression);
end

function [Source, Target, Mode, ODERow, ODEExpression, UsedInODE] = recordEdge(k, src, tgt, mode, lineIdx, line, Source, Target, Mode, ODERow, ODEExpression, UsedInODE, numEdges)
if isnan(k) || k < 1 || k > numEdges
    return;
end
Source{k} = appendUniqueToken(Source{k}, src);
Target{k} = appendUniqueToken(Target{k}, tgt);
Mode{k} = appendUniqueToken(Mode{k}, mode);
if isnan(ODERow(k))
    ODERow(k) = lineIdx;
end
if isempty(ODEExpression{k})
    ODEExpression{k} = line;
end
UsedInODE(k) = true;
end

function out = appendUniqueToken(in, token)
if isempty(token)
    out = in;
    return;
end
if isempty(in)
    out = token;
    return;
end
parts = regexp(in, '\s*;\s*', 'split');
if any(strcmp(parts, token))
    out = in;
else
    out = [in, '; ', token];
end
end

function RangeTable = createDefaultRangeTable(EdgeTable, globalNominal, defaultRange)
numEdges = height(EdgeTable);

NominalW = repmat(globalNominal.w, numEdges, 1);
LowW = globalNominal.w .* defaultRange.wFactor(1) .* ones(numEdges,1);
HighW = min(1.0, globalNominal.w .* defaultRange.wFactor(2)) .* ones(numEdges,1);

NominalN = repmat(globalNominal.n, numEdges, 1);
LowN = defaultRange.nAbs(1) .* ones(numEdges,1);
HighN = defaultRange.nAbs(2) .* ones(numEdges,1);

NominalEC50 = repmat(globalNominal.EC50, numEdges, 1);
LowEC50 = defaultRange.EC50Abs(1) .* ones(numEdges,1);
HighEC50 = defaultRange.EC50Abs(2) .* ones(numEdges,1);

ManualOverride = false(numEdges,1);
Notes = repmat({'default calibration range'}, numEdges, 1);

RangeTable = [EdgeTable, table(NominalW, LowW, HighW, NominalN, LowN, HighN, ...
    NominalEC50, LowEC50, HighEC50, ManualOverride, Notes)];
end

function RangeTable = applyManualRangeSpec(RangeTable, manualRangeSpec)
if isempty(manualRangeSpec)
    return;
end

for i = 1:size(manualRangeSpec,1)
    src = char(manualRangeSpec{i,1});
    tgt = char(manualRangeSpec{i,2});
    mode = char(manualRangeSpec{i,3});

    sourceMatch = contains(lower(RangeTable.Source), lower(src));
    targetMatch = strcmpi(RangeTable.Target, tgt) | contains(lower(RangeTable.Target), lower(tgt));

    if strcmpi(mode, 'any')
        modeMatch = true(height(RangeTable),1);
    else
        modeMatch = contains(lower(RangeTable.Mode), lower(mode));
    end

    matchIdx = find(sourceMatch & targetMatch & modeMatch);

    if isempty(matchIdx)
        warning('Manual range row %d did not match any inferred edge: %s -> %s (%s).', i, src, tgt, mode);
        continue;
    end

    RangeTable.NominalW(matchIdx) = manualRangeSpec{i,4};
    RangeTable.LowW(matchIdx) = manualRangeSpec{i,5};
    RangeTable.HighW(matchIdx) = manualRangeSpec{i,6};
    RangeTable.NominalN(matchIdx) = manualRangeSpec{i,7};
    RangeTable.LowN(matchIdx) = manualRangeSpec{i,8};
    RangeTable.HighN(matchIdx) = manualRangeSpec{i,9};
    RangeTable.NominalEC50(matchIdx) = manualRangeSpec{i,10};
    RangeTable.LowEC50(matchIdx) = manualRangeSpec{i,11};
    RangeTable.HighEC50(matchIdx) = manualRangeSpec{i,12};
    RangeTable.ManualOverride(matchIdx) = true;
    RangeTable.Notes(matchIdx) = manualRangeSpec(i,13);
end

RangeTable.LowW = max(0, min(RangeTable.LowW, RangeTable.HighW));
RangeTable.HighW = min(1, max(RangeTable.LowW, RangeTable.HighW));
RangeTable.LowN = max(0.1, min(RangeTable.LowN, RangeTable.HighN));
RangeTable.HighN = max(RangeTable.LowN, RangeTable.HighN);
RangeTable.LowEC50 = max(0.01, min(RangeTable.LowEC50, RangeTable.HighEC50));
RangeTable.HighEC50 = max(RangeTable.LowEC50, RangeTable.HighEC50);
end

function rpar = sampleRparFromRanges(RangeTable)
numEdges = height(RangeTable);
rpar = nan(3, numEdges);

for k = 1:numEdges
    w = uniformSample(RangeTable.LowW(k), RangeTable.HighW(k));
    n = uniformSample(RangeTable.LowN(k), RangeTable.HighN(k));

    validEC50Upper = 0.98 * (0.5 .^ (1 ./ n));
    ec50Low = RangeTable.LowEC50(k);
    ec50High = min(RangeTable.HighEC50(k), validEC50Upper);

    if ec50High <= ec50Low
        ec50High = validEC50Upper;
        ec50Low = min(ec50Low, 0.90 * validEC50Upper);
        ec50Low = max(0.01, ec50Low);
    end

    EC50 = uniformSample(ec50Low, ec50High);

    rpar(1,k) = w;
    rpar(2,k) = n;
    rpar(3,k) = EC50;
    rpar(1,6) = 0;
end
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
        idx(end+1) = thisIdx; 
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

function ParameterLongTable = buildLongParameterTable(paramSamples, EdgeTable, runIds)
[~, numEdges, numMC] = size(paramSamples);
if nargin < 3 || isempty(runIds)
    runIds = 1:numMC;
end
N = numEdges * numMC;
Run = zeros(N,1);
EdgeIndex = zeros(N,1);
Source = cell(N,1);
Target = cell(N,1);
Mode = cell(N,1);
W = zeros(N,1);
HillN = zeros(N,1);
EC50 = zeros(N,1);

row = 0;
for localRunIdx = 1:numMC
    for edgeIdx = 1:numEdges
        row = row + 1;
        Run(row) = runIds(localRunIdx);
        EdgeIndex(row) = edgeIdx;
        Source{row} = EdgeTable.Source{edgeIdx};
        Target{row} = EdgeTable.Target{edgeIdx};
        Mode{row} = EdgeTable.Mode{edgeIdx};
        W(row) = paramSamples(1,edgeIdx,localRunIdx);
        HillN(row) = paramSamples(2,edgeIdx,localRunIdx);
        EC50(row) = paramSamples(3,edgeIdx,localRunIdx);
    end
end

ParameterLongTable = table(Run, EdgeIndex, Source, Target, Mode, W, HillN, EC50);
end

function T = buildSingleRunParameterTable(rpar, EdgeTable, runId, scorePct)
numEdges = size(rpar,2);
Run = repmat(runId, numEdges, 1);
ValidationAgreementPct = repmat(scorePct, numEdges, 1);
W = rpar(1,:)';
HillN = rpar(2,:)';
EC50 = rpar(3,:)';
T = [table(Run, ValidationAgreementPct), EdgeTable, table(W, HillN, EC50)];
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

function makeSummaryPlots(resultsDir, finalOutputs, outputSpeciesNamesFound, validRuns, tSim, fullTimeCourses, outputIdx)
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
title('Baseline final output values across Monte Carlo runs');
set(fig2, 'Position', [100 100 1200 500]);
saveas(fig2, fullfile(resultsDir, 'final_output_heatmap.png'));
savefig(fig2, fullfile(resultsDir, 'final_output_heatmap.fig'));

if ~isempty(tSim) && ~isempty(fullTimeCourses)
    fig3 = figure('Name','Monte Carlo output time course envelopes','Color','w');
    nRows = ceil(sqrt(numOutputs));
    nCols = ceil(numOutputs / nRows);

    for j = 1:numOutputs
        subplot(nRows, nCols, j);
        y = squeeze(fullTimeCourses(:, outputIdx(j), validRuns));
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
title('Validation agreement by sampled parameter set');
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
title(sprintf('Best run %d: %.2f%% agreement with validation matrix', bestRun, bestScore));
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
