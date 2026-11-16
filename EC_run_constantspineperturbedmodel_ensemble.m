function ensembleResults = EC_run_constantspineperturbedmodel_ensemble(cfg)
% Optimize binary donor substitutions without altering the ODE model.

cfg = applyDefaults(cfg);
cfg = applyQuickTestSettings(cfg);
validateConfig(cfg);

enginePath = mfilename('fullpath');
modelDir = fileparts(enginePath);
originalDir = pwd;
cleanupDir = onCleanup(@() cd(originalDir));
cd(modelDir);

cfg.constantResultsDir = resolveFromModelDir(cfg.constantResultsDir, modelDir);
cfg.perturbedResultsDir = resolveFromModelDir(cfg.perturbedResultsDir, modelDir);
cfg.outputDir = resolveFromModelDir(cfg.outputDir, modelDir);

assert(exist(cfg.constantResultsDir, 'dir') == 7, ...
    'Constant result folder does not exist: %s', cfg.constantResultsDir);
assert(exist(cfg.perturbedResultsDir, 'dir') == 7, ...
    'Perturbed result folder does not exist: %s', cfg.perturbedResultsDir);
ensureDir(cfg.outputDir);

model = loadModelContext(cfg, modelDir);
EdgeTable = loadEdgeTable(cfg.perturbedResultsDir, model.numEdges);
model.EdgeTable = EdgeTable;

[constantSpines, ConstantCandidateTable] = loadConstantSpines(cfg, model);
[perturbedDonors, PerturbedCandidateTable] = loadPerturbedDonors(cfg, model);

assert(~isempty(constantSpines), 'No constant spine models passed the selection rules.');
assert(~isempty(perturbedDonors), 'No perturbed donor models passed the selection rules.');

writetable(ConstantCandidateTable, fullfile(cfg.outputDir, 'constant_spine_candidates.csv'));
writetable(PerturbedCandidateTable, fullfile(cfg.outputDir, 'perturbed_donor_candidates.csv'));
atomicSaveMat(fullfile(cfg.outputDir, 'ensemble_configuration.mat'), struct('cfg', cfg));

fprintf('\nConstant-spine / perturbed-donor ensemble\n');
fprintf('  Version: %s\n', cfg.version);
fprintf('  Selected constant spines: %d\n', numel(constantSpines));
fprintf('  Selected perturbed donors: %d\n', numel(perturbedDonors));
fprintf('  Pair combinations to screen: %d\n\n', numel(constantSpines) * numel(perturbedDonors));

screenSignature = buildScreenSignature(cfg, modelDir, constantSpines, perturbedDonors);
screenResults = runOrLoadPairScreening( ...
    constantSpines, perturbedDonors, model, cfg, screenSignature);
ScreeningTable = buildPairResultTable(screenResults);
ScreeningTable = sortPairResultTable(ScreeningTable);
writetable(ScreeningTable, fullfile(cfg.outputDir, 'screening_pair_results_ranked.csv'));

finalistRows = selectFinalistRows( ...
    ScreeningTable, cfg.numFinalistPairs, cfg.diversifyFinalistSources);
FinalistTable = ScreeningTable(finalistRows, :);
writetable(FinalistTable, fullfile(cfg.outputDir, 'selected_finalist_pairs.csv'));

fprintf('\nSelected %d finalist spine/donor pairs for the robust GA.\n', height(FinalistTable));
estimatedModelsPerRun = cfg.populationSize + ...
    cfg.numGenerations * (cfg.populationSize - cfg.eliteCount);
estimatedOdesPerRun = estimatedModelsPerRun * ...
    (1 + numel(model.validationPerturbIdx));
fprintf('Planned final runs: %d pairs x %d seeds = %d runs.\n', ...
    height(FinalistTable), numel(cfg.rngSeeds), ...
    height(FinalistTable) * numel(cfg.rngSeeds));
fprintf(['Approximate upper search scale per final run: %g unique masks, ', ...
    '%g ODE integrations before pruning and cache duplicates.\n'], ...
    estimatedModelsPerRun, estimatedOdesPerRun);
finalResults = runFinalistSearches( ...
    FinalistTable, screenResults, constantSpines, perturbedDonors, model, cfg, screenSignature);

FinalRunTable = buildPairResultTable(finalResults);
FinalRunTable = sortPairResultTable(FinalRunTable);
writetable(FinalRunTable, fullfile(cfg.outputDir, 'final_run_results_ranked.csv'));

representativeIdx = chooseRepresentativeResult(finalResults);
representative = finalResults{representativeIdx};

[EdgeEnsembleTable, ensembleStats] = buildEnsembleEdgeTable(finalResults, EdgeTable);
writetable(EdgeEnsembleTable, ...
    fullfile(cfg.outputDir, 'ensemble_edge_selection_and_parameter_ranges.csv'));

ConstantOutcomeTable = summarizeSourceOutcomes(finalResults, 'constant');
PerturbedOutcomeTable = summarizeSourceOutcomes(finalResults, 'perturbed');
writetable(ConstantOutcomeTable, fullfile(cfg.outputDir, 'constant_spine_outcomes.csv'));
writetable(PerturbedOutcomeTable, fullfile(cfg.outputDir, 'perturbed_donor_outcomes.csv'));

representativeDir = fullfile(cfg.outputDir, 'representative_best_model');
ensureDir(representativeDir);
writePairOutputs(representativeDir, representative, EdgeTable, model, cfg, ...
    cfg.makeRepresentativeValidationPlot);

overlayResults = struct();
if cfg.makeWssOverlay
    overlayResults = makeEnsembleWssOverlay(finalResults, representativeIdx, model, cfg);
end

ensembleResults = struct();
ensembleResults.configuration = cfg;
ensembleResults.constantSpines = constantSpines;
ensembleResults.perturbedDonors = perturbedDonors;
ensembleResults.screeningResults = screenResults;
ensembleResults.finalResults = finalResults;
ensembleResults.representativeIndex = representativeIdx;
ensembleResults.representative = representative;
ensembleResults.ensembleStats = ensembleStats;
ensembleResults.overlayResults = overlayResults;

atomicSaveMat(fullfile(cfg.outputDir, ...
    'constant_spine_perturbed_ensemble_results.mat'), ...
    struct('ensembleResults', ensembleResults));

fprintf('\nRepresentative hybrid:\n');
fprintf('  Constant run: %g\n', representative.constantRun);
fprintf('  Perturbed donor run: %g\n', representative.perturbedRun);
fprintf('  GA seed: %g\n', representative.rngSeed);
fprintf('  Validation: %.2f%% (%d/%d)\n', representative.bestEval.scorePct, ...
    representative.bestEval.nMatch, representative.bestEval.nKnown);
fprintf('  Selected perturbed edges: %d/%d\n', representative.bestEval.numSelected, ...
    representative.numCandidates);
fprintf('Results saved to:\n  %s\n', cfg.outputDir);

clear cleanupDir;
end

function cfg = applyDefaults(cfg)
defaults = struct();
defaults.constantResultsDir = fullfile('results','constant');
defaults.perturbedResultsDir = fullfile('results','reaction_specific');
defaults.version = 'AngII_v5';
defaults.outputDir = fullfile('results','hybrid');
defaults.constantScoreFile = 'sampled_global_constant_parameters.csv';
defaults.perturbedScoreFile = 'top_validation_scores.csv';
defaults.perturbedParameterFile = 'top_parameter_sets_long.csv';
defaults.rerankSourceModelsForCurrentValidation = true;
defaults.constantAllResultsMatFile = ...
    'MonteCarlo_constant_parameter_validation_fit_all_results.mat';
defaults.perturbedAllResultsMatFile = ...
    'MonteCarlo_validation_fit_all_results.mat';
defaults.sourceScoreChunkSize = 2000;
defaults.constantScoreToleranceMatches = 0;
defaults.keepAllMaxScoreConstant = true;
defaults.maxConstantSpines = Inf;
defaults.perturbedScoreToleranceMatches = 1;
defaults.keepAllMaxScorePerturbed = true;
defaults.maxPerturbedDonors = 10;
defaults.diversityUseModeledEdgesOnly = true;
defaults.runPairScreening = true;
defaults.reuseExistingScreening = false;
defaults.screenMaxPrunePasses = 2;
defaults.screenCheckpointEvery = 10;
defaults.useParallelPairScreening = true;
defaults.numFinalistPairs = 5;
defaults.diversifyFinalistSources = true;
defaults.populationSize = 250;
defaults.numGenerations = 200;
defaults.rngSeeds = 1:5;
defaults.eliteCount = 5;
defaults.tournamentSize = 3;
defaults.mutationRate = 0.03;
defaults.complexityPenalty = 1e-3;
defaults.finalMaxPrunePasses = 2;
defaults.runGreedyForward = false;
defaults.runGreedyPrune = true;
defaults.runGreedyForwardAfterPrune = false;
defaults.useParallelPopulationEvaluation = true;
defaults.parallelWorkers = [];
defaults.resumeCompletedFinalRuns = true;
defaults.validationDirectionThreshold = 0.02;
defaults.candidateMode = 'used';
defaults.activeInputs = 1:5;
defaults.inputValue = 0.2;
defaults.stressInputValue = 0.5;
defaults.perturbValue = 0.5;
defaults.validationPerturbSpecies = {'Stress','Wss','Hypoxia','TNFa'};
defaults.validationSpeciesSubset = {'VEGF','VEGFR','IL6','ppET1','ET1','eNOS', ...
    'FAK','TGFb','pSmad23','pSmad15','p38','ERK','JNK','PKC','DAG', ...
    'PI3K','Akt','mTOR','p70S6K','RhoA','ROCK','MLCK','Actin','NFkB', ...
    'HIF1a','MCP1','TSP1','Proliferation','Migration','Apoptosis','EndMT'};
defaults.saveIndividualFinalRunFiles = true;
defaults.saveIndividualValidationPlots = false;
defaults.makeRepresentativeValidationPlot = true;
defaults.makeWssOverlay = false;
defaults.shearMaxDyn = 36;
defaults.numShearPoints = 121;
defaults.maxOverlayModels = 25;
defaults.closeFiguresAfterSaving = true;
defaults.stopOnODEError = false;
defaults.verbose = true;
defaults.quickTestMode = false;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg, names{i})
        cfg.(names{i}) = defaults.(names{i});
    end
end
end

function cfg = applyQuickTestSettings(cfg)
if ~cfg.quickTestMode
    return;
end

cfg.outputDir = fullfile(cfg.outputDir, 'QuickTest');
cfg.keepAllMaxScoreConstant = false;
cfg.keepAllMaxScorePerturbed = false;
cfg.maxConstantSpines = 2;
cfg.maxPerturbedDonors = 2;
cfg.numFinalistPairs = 1;
cfg.populationSize = 8;
cfg.numGenerations = 1;
cfg.rngSeeds = 1;
cfg.eliteCount = 2;
cfg.screenMaxPrunePasses = 1;
cfg.finalMaxPrunePasses = 1;
cfg.screenCheckpointEvery = 1;
cfg.reuseExistingScreening = false;
cfg.resumeCompletedFinalRuns = false;
cfg.makeWssOverlay = false;
fprintf('QUICK TEST MODE: using reduced ensemble and GA settings.\n');
end

function validateConfig(cfg)
assert(cfg.constantScoreToleranceMatches >= 0, ...
    'constantScoreToleranceMatches must be nonnegative.');
assert(cfg.perturbedScoreToleranceMatches >= 0, ...
    'perturbedScoreToleranceMatches must be nonnegative.');
assert(cfg.numFinalistPairs >= 1 && isfinite(cfg.numFinalistPairs), ...
    'numFinalistPairs must be a positive finite integer.');
assert(cfg.populationSize >= 4 && isfinite(cfg.populationSize), ...
    'populationSize must be at least 4.');
assert(cfg.numGenerations >= 1 && isfinite(cfg.numGenerations), ...
    'numGenerations must be at least 1.');
assert(cfg.eliteCount >= 1 && cfg.eliteCount <= cfg.populationSize - 2, ...
    'eliteCount must be between 1 and populationSize - 2.');
assert(cfg.tournamentSize >= 2 && cfg.tournamentSize <= cfg.populationSize, ...
    'tournamentSize must be between 2 and populationSize.');
assert(cfg.mutationRate >= 0 && cfg.mutationRate <= 1, ...
    'mutationRate must be between 0 and 1.');
assert(isempty(cfg.parallelWorkers) || ...
    (isscalar(cfg.parallelWorkers) && cfg.parallelWorkers >= 1 && ...
    isfinite(cfg.parallelWorkers)), ...
    'parallelWorkers must be empty or a positive finite scalar.');
assert(~isempty(cfg.rngSeeds) && all(isfinite(cfg.rngSeeds)), ...
    'rngSeeds must contain at least one finite value.');
assert(any(strcmpi(cfg.candidateMode, {'used','all'})), ...
    'candidateMode must be ''used'' or ''all''.');
assert(cfg.screenCheckpointEvery >= 1, ...
    'screenCheckpointEvery must be at least 1.');
assert(cfg.sourceScoreChunkSize >= 1 && isfinite(cfg.sourceScoreChunkSize), ...
    'sourceScoreChunkSize must be a positive finite integer.');
end

function pathOut = resolveFromModelDir(pathIn, modelDir)
if isAbsolutePath(pathIn)
    pathOut = pathIn;
else
    pathOut = fullfile(modelDir, pathIn);
end
end

function tf = isAbsolutePath(pathIn)
tf = ~isempty(regexp(pathIn, '^[A-Za-z]:[\\/]|^[/\\]', 'once'));
end

function ensureDir(folder)
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end

function model = loadModelContext(cfg, modelDir)
refFile = fullfile(modelDir, ['refcase_EC_Network_', cfg.version, '.mat']);
ODEfilename = ['ODElist_EC_Network_', cfg.version];
odeFile = fullfile(modelDir, [ODEfilename, '.mat']);

assert(exist(refFile, 'file') == 2, 'Could not find model reference file: %s', refFile);
assert(exist(odeFile, 'file') == 2, 'Could not find ODE list file: %s', odeFile);
assert(exist(fullfile(modelDir, 'ODE_master.m'), 'file') == 2, ...
    'Could not find ODE_master.m in %s.', modelDir);
assert(exist(fullfile(modelDir, 'CreateValidationMatrix.m'), 'file') == 2, ...
    'Could not find CreateValidationMatrix.m in %s.', modelDir);

S = load(refFile, 'reactionRules', 'tau', 'ymax', 'speciesNames', 'numEdges');
model = S;
model.speciesNames = S.speciesNames(:);
model.ODEfilename = ODEfilename;
model.modelDir = modelDir;

y0 = zeros(1, numel(model.speciesNames));
activeInputs = cfg.activeInputs(cfg.activeInputs <= numel(y0));
y0(activeInputs) = cfg.inputValue;
stressIdx = find(strcmp(model.speciesNames, 'Stress'), 1);
if ~isempty(stressIdx)
    y0(stressIdx) = cfg.stressInputValue;
end
model.y0Base = y0;

model.validationPerturbIdx = findSpeciesIndices( ...
    model.speciesNames, cfg.validationPerturbSpecies, true);
assert(numel(model.validationPerturbIdx) == numel(cfg.validationPerturbSpecies), ...
    'Every validation perturbation species must exist in the model.');

model.validationSpeciesIdx = findSpeciesIndices( ...
    model.speciesNames, cfg.validationSpeciesSubset, true);
model.ExpValidation = CreateValidationMatrix( ...
    model.speciesNames, numel(cfg.validationPerturbSpecies));
model.knownValidationMask = abs(model.ExpValidation) ~= 0.5;
end

function EdgeTable = loadEdgeTable(resultsDir, numEdges)
edgeFile = fullfile(resultsDir, 'inferred_rpar_edge_table.csv');
bestFile = fullfile(resultsDir, 'best_rpar_parameter_set.csv');

if exist(edgeFile, 'file') == 2
    EdgeTable = readtable(edgeFile);
elseif exist(bestFile, 'file') == 2
    EdgeTable = readtable(bestFile);
else
    EdgeTable = table((1:numEdges)', repmat({''}, numEdges, 1), ...
        repmat({''}, numEdges, 1), repmat({''}, numEdges, 1), ...
        true(numEdges, 1), 'VariableNames', ...
        {'EdgeIndex','Source','Target','Mode','UsedInODE'});
end

if ~ismember('EdgeIndex', EdgeTable.Properties.VariableNames)
    EdgeTable.EdgeIndex = (1:height(EdgeTable))';
end
if ~ismember('UsedInODE', EdgeTable.Properties.VariableNames)
    EdgeTable.UsedInODE = true(height(EdgeTable), 1);
end

EdgeTable = sortrows(EdgeTable, 'EdgeIndex');
assert(height(EdgeTable) >= numEdges, ...
    'The edge table has %d rows, but the model expects %d.', height(EdgeTable), numEdges);
EdgeTable = EdgeTable(1:numEdges, :);
assert(isequal(EdgeTable.EdgeIndex(:), (1:numEdges)'), ...
    'The edge table must contain EdgeIndex 1:%d.', numEdges);
end

function [spines, CandidateTable] = loadConstantSpines(cfg, model)
numEdges = model.numEdges;
scoreFile = fullfile(cfg.constantResultsDir, cfg.constantScoreFile);
if exist(scoreFile, 'file') ~= 2
    fallback = fullfile(cfg.constantResultsDir, 'top_validation_scores.csv');
    assert(exist(fallback, 'file') == 2, ...
        'Could not find constant score file %s or fallback %s.', scoreFile, fallback);
    warning('Using fallback constant score file: %s', fallback);
    scoreFile = fallback;
end

T = readtable(scoreFile);
required = {'Run','W','HillN','EC50'};
assert(all(ismember(required, T.Properties.VariableNames)), ...
    'Constant score file must contain: %s.', strjoin(required, ', '));

if cfg.rerankSourceModelsForCurrentValidation
    resultFile = fullfile(cfg.constantResultsDir, cfg.constantAllResultsMatFile);
    [allMatches, allKnown, allScorePct, allValid] = ...
        reclassifySourceValidation(resultFile, model, cfg, 'constant');
    runIdx = round(T.Run);
    assert(all(runIdx >= 1 & runIdx <= numel(allMatches)), ...
        'Constant score table contains run indices outside the saved result workspace.');
    matches = allMatches(runIdx);
    known = allKnown(runIdx);
    scorePct = allScorePct(runIdx);
    valid = allValid(runIdx);
    scoreBasis = repmat("current validation + threshold", height(T), 1);
else
    [matches, known, scorePct, valid] = extractSavedScores(T);
    scoreBasis = repmat("saved source score", height(T), 1);
end
valid = valid & all(isfinite([T.W, T.HillN, T.EC50]), 2);
assert(any(valid), 'No valid constant models were found in %s.', scoreFile);

topMatches = max(matches(valid));
eligible = valid & matches >= topMatches - cfg.constantScoreToleranceMatches;
eligibleIdx = find(eligible);
featureMatrix = [T.W(eligibleIdx), T.HillN(eligibleIdx), T.EC50(eligibleIdx)];
requiredLocal = find(matches(eligibleIdx) == topMatches);

target = cfg.maxConstantSpines;
if cfg.keepAllMaxScoreConstant
    target = maxFiniteOrRequired(target, numel(requiredLocal));
else
    requiredLocal = requiredLocal(1:min(1, numel(requiredLocal)));
end

selectedLocal = selectDiverseRows(featureMatrix, requiredLocal, target);
selectedIdx = eligibleIdx(selectedLocal);

[~, order] = sortrows([matches(selectedIdx), -T.Run(selectedIdx)], [-1, -2]);
selectedIdx = selectedIdx(order);

spines = repmat(emptySourceModel(), numel(selectedIdx), 1);
for i = 1:numel(selectedIdx)
    row = selectedIdx(i);
    spines(i).run = T.Run(row);
    spines(i).savedMatches = matches(row);
    spines(i).savedKnown = known(row);
    spines(i).savedScorePct = scorePct(row);
    spines(i).W = T.W(row);
    spines(i).HillN = T.HillN(row);
    spines(i).EC50 = T.EC50(row);
    spines(i).rpar = repmat([T.W(row); T.HillN(row); T.EC50(row)], 1, numEdges);
end

Selected = ismember((1:height(T))', selectedIdx);
Eligible = eligible;
ExactTopScore = valid & matches == topMatches;
CandidateTable = table(T.Run, matches, known, scorePct, T.W, T.HillN, T.EC50, ...
    ExactTopScore, Eligible, Selected, scoreBasis, ...
    'VariableNames', {'Run','SavedValidationMatches','SavedValidationKnownChecks', ...
    'SavedValidationAgreementPct','W','HillN','EC50','ExactTopScore', ...
    'ScoreEligible','SelectedAsSpine','ScoreBasis'});
CandidateTable = CandidateTable(valid, :);
CandidateTable = sortrows(CandidateTable, ...
    {'SavedValidationMatches','Run'}, {'descend','ascend'});
end

function [donors, CandidateTable] = loadPerturbedDonors(cfg, model)
useCurrentScores = cfg.rerankSourceModelsForCurrentValidation;
if useCurrentScores
    resultFile = fullfile(cfg.perturbedResultsDir, cfg.perturbedAllResultsMatFile);
    [matches, known, scorePct, valid] = ...
        reclassifySourceValidation(resultFile, model, cfg, 'perturbed');
    ScoreTable = table((1:numel(matches))', 'VariableNames', {'Run'});
    scoreBasis = repmat("current validation + threshold", height(ScoreTable), 1);
else
    scoreFile = fullfile(cfg.perturbedResultsDir, cfg.perturbedScoreFile);
    assert(exist(scoreFile, 'file') == 2, ...
        'Could not find perturbed score file: %s', scoreFile);
    ScoreTable = readtable(scoreFile);
    assert(ismember('Run', ScoreTable.Properties.VariableNames), ...
        'Perturbed score file must contain a Run column.');
    [matches, known, scorePct, valid] = extractSavedScores(ScoreTable);
    scoreBasis = repmat("saved source score", height(ScoreTable), 1);
end
assert(any(valid), 'No valid perturbed models were found.');

topMatches = max(matches(valid));
eligible = valid & matches >= topMatches - cfg.perturbedScoreToleranceMatches;
eligibleRuns = ScoreTable.Run(eligible);

if useCurrentScores
    parameterFile = fullfile(cfg.perturbedResultsDir, cfg.perturbedAllResultsMatFile);
    ParameterMat = matfile(parameterFile);
    parameterSize = size(ParameterMat, 'paramSamples');
    assert(numel(parameterSize) == 3 && parameterSize(1) == 3 && ...
        parameterSize(2) == model.numEdges, ...
        'paramSamples in %s must have size 3 x %d x number-of-runs.', ...
        parameterFile, model.numEdges);
    availableRuns = (1:parameterSize(3))';
    ParameterTable = table();
else
    parameterFile = fullfile(cfg.perturbedResultsDir, cfg.perturbedParameterFile);
    assert(exist(parameterFile, 'file') == 2, ...
        'Could not find perturbed parameter file: %s', parameterFile);
    ParameterTable = readtable(parameterFile);
    required = {'Run','EdgeIndex','W','HillN','EC50'};
    assert(all(ismember(required, ParameterTable.Properties.VariableNames)), ...
        'Perturbed parameter file must contain: %s.', strjoin(required, ', '));
    availableRuns = unique(ParameterTable.Run);
end
missingRuns = setdiff(eligibleRuns, availableRuns);
if ~isempty(missingRuns)
    warning('%d score-qualified perturbed runs have no saved parameter matrix and will be skipped.', ...
        numel(missingRuns));
end
eligibleRuns = intersect(eligibleRuns, availableRuns, 'stable');
assert(~isempty(eligibleRuns), ...
    'None of the score-qualified perturbed runs have parameters in %s.', parameterFile);

allEligible = repmat(emptySourceModel(), numel(eligibleRuns), 1);
if useCurrentScores
    parameterStack = loadPerturbedParameterRuns( ...
        ParameterMat, eligibleRuns, model.numEdges);
end
for i = 1:numel(eligibleRuns)
    runId = eligibleRuns(i);
    if useCurrentScores
        rpar = parameterStack(:, :, i);
    else
        rows = ParameterTable.Run == runId;
        P = ParameterTable(rows, :);
        P = sortrows(P, 'EdgeIndex');
        assert(height(P) == model.numEdges, ...
            'Perturbed run %g has %d edges; expected %d.', ...
            runId, height(P), model.numEdges);
        assert(isequal(P.EdgeIndex(:), (1:model.numEdges)'), ...
            'Perturbed run %g does not contain EdgeIndex 1:%d.', ...
            runId, model.numEdges);
        rpar = [P.W'; P.HillN'; P.EC50'];
    end

    scoreRow = find(ScoreTable.Run == runId & valid, 1);
    allEligible(i).run = runId;
    allEligible(i).savedMatches = matches(scoreRow);
    allEligible(i).savedKnown = known(scoreRow);
    allEligible(i).savedScorePct = scorePct(scoreRow);
    allEligible(i).rpar = rpar;
end

if cfg.diversityUseModeledEdgesOnly
    diversityEdges = find(logical(model.EdgeTable.UsedInODE));
else
    diversityEdges = 1:model.numEdges;
end
features = flattenRparFeatures(allEligible, diversityEdges);
requiredLocal = find([allEligible.savedMatches] == topMatches);

target = cfg.maxPerturbedDonors;
if cfg.keepAllMaxScorePerturbed
    target = maxFiniteOrRequired(target, numel(requiredLocal));
else
    requiredLocal = requiredLocal(1:min(1, numel(requiredLocal)));
end

selectedLocal = selectDiverseRows(features, requiredLocal, target);
donors = allEligible(selectedLocal);

selectedRuns = [donors.run]';
ExactTopScore = matches == topMatches & valid;
ScoreEligible = eligible & ismember(ScoreTable.Run, availableRuns);
Selected = ismember(ScoreTable.Run, selectedRuns);
CandidateTable = table(ScoreTable.Run, matches, known, scorePct, ExactTopScore, ...
    ScoreEligible, Selected, scoreBasis, ...
    'VariableNames', {'Run','SavedValidationMatches','SavedValidationKnownChecks', ...
    'SavedValidationAgreementPct','ExactTopScore','ScoreEligible', ...
    'SelectedAsDonor','ScoreBasis'});
CandidateTable = CandidateTable(valid, :);
CandidateTable = sortrows(CandidateTable, ...
    {'SavedValidationMatches','Run'}, {'descend','ascend'});
end

function model = emptySourceModel()
model = struct('run', NaN, 'savedMatches', NaN, 'savedKnown', NaN, ...
    'savedScorePct', NaN, 'W', NaN, 'HillN', NaN, 'EC50', NaN, ...
    'rpar', []);
end

function [matches, known, scorePct, valid] = extractSavedScores(T)
n = height(T);
if ismember('ValidationKnownChecks', T.Properties.VariableNames)
    known = T.ValidationKnownChecks;
else
    known = nan(n, 1);
end
if ismember('ValidationMatches', T.Properties.VariableNames)
    matches = T.ValidationMatches;
elseif ismember('ValidationAgreementPct', T.Properties.VariableNames) && all(isfinite(known))
    matches = round(T.ValidationAgreementPct .* known ./ 100);
else
    error('Score table needs ValidationMatches, or agreement percentage plus known-check count.');
end
if ismember('ValidationAgreementPct', T.Properties.VariableNames)
    scorePct = T.ValidationAgreementPct;
else
    scorePct = 100 .* matches ./ known;
end
if ismember('Failed', T.Properties.VariableNames)
    valid = ~logical(T.Failed);
else
    valid = true(n, 1);
end
valid = valid & isfinite(matches) & isfinite(scorePct);
end

function [matches, known, scorePct, valid] = ...
        reclassifySourceValidation(resultFile, model, cfg, sourceLabel)
assert(exist(resultFile, 'file') == 2, ...
    'Could not find the %s all-results workspace: %s', sourceLabel, resultFile);

ResultMat = matfile(resultFile);
variableNames = who(ResultMat);
required = {'validationDiffStore','failedRuns'};
assert(all(ismember(required, variableNames)), ...
    '%s must contain validationDiffStore and failedRuns.', resultFile);

diffSize = size(ResultMat, 'validationDiffStore');
expectedSize = size(model.ExpValidation);
assert(numel(diffSize) == 3 && isequal(diffSize(2:3), expectedSize), ...
    ['validationDiffStore in %s has size %s; expected number-of-runs x ', ...
    '%d x %d.'], resultFile, mat2str(diffSize), expectedSize(1), expectedSize(2));

numRuns = diffSize(1);
knownCount = nnz(model.knownValidationMask);
matches = zeros(numRuns, 1);
known = repmat(knownCount, numRuns, 1);
expected = reshape(model.ExpValidation, [1, expectedSize]);
knownMask = reshape(model.knownValidationMask, [1, expectedSize]);
chunkSize = max(1, round(cfg.sourceScoreChunkSize));
numChunks = ceil(numRuns / chunkSize);

fprintf('Reclassifying %d %s source models using the current validation objective', ...
    numRuns, sourceLabel);
fprintf(' and a +/-%.4g no-change threshold.\n', ...
    cfg.validationDirectionThreshold);
for firstRun = 1:chunkSize:numRuns
    lastRun = min(numRuns, firstRun + chunkSize - 1);
    rows = firstRun:lastRun;
    diffBlock = ResultMat.validationDiffStore(rows, :, :);
    diffBlock(abs(diffBlock) < cfg.validationDirectionThreshold) = 0;
    qualitative = sign(diffBlock);
    matchBlock = (qualitative == expected) & knownMask;
    matches(rows) = reshape(sum(sum(matchBlock, 3), 2), [], 1);

    chunkIndex = ceil(lastRun / chunkSize);
    if cfg.verbose && (chunkIndex == 1 || chunkIndex == numChunks || ...
            mod(chunkIndex, max(1, floor(numChunks / 5))) == 0)
        fprintf('  %s reranking: %d/%d runs\n', sourceLabel, lastRun, numRuns);
    end
end

failedRuns = logical(ResultMat.failedRuns);
failedRuns = failedRuns(:);
assert(numel(failedRuns) == numRuns, ...
    'failedRuns in %s does not match validationDiffStore.', resultFile);
valid = ~failedRuns & isfinite(matches);
scorePct = 100 .* matches ./ known;
scorePct(~valid) = NaN;
end

function parameterStack = loadPerturbedParameterRuns(ParameterMat, runIds, numEdges)
runIds = round(runIds(:)');
parameterStack = nan(3, numEdges, numel(runIds));
try
    parameterStack = ParameterMat.paramSamples(:, :, runIds);
    parameterStack = reshape(parameterStack, 3, numEdges, numel(runIds));
catch
    for i = 1:numel(runIds)
        parameterStack(:, :, i) = ParameterMat.paramSamples(:, :, runIds(i));
    end
end
assert(all(isfinite(parameterStack(:))), ...
    'One or more selected perturbed parameter matrices contain non-finite values.');
end

function target = maxFiniteOrRequired(target, requiredCount)
if isinf(target)
    return;
end
target = max(target, requiredCount);
end

function selected = selectDiverseRows(X, requiredRows, target)
n = size(X, 1);
if n == 0
    selected = [];
    return;
end
if isinf(target) || target >= n
    selected = 1:n;
    return;
end

target = max(1, round(target));
requiredRows = unique(requiredRows(:)', 'stable');
requiredRows = requiredRows(requiredRows >= 1 & requiredRows <= n);
selected = requiredRows;
if isempty(selected)
    selected = 1;
end
if numel(selected) >= target
    return;
end

Z = standardizeColumns(X);
remaining = setdiff(1:n, selected, 'stable');
while numel(selected) < target && ~isempty(remaining)
    minDistance = inf(numel(remaining), 1);
    for s = selected
        delta = Z(remaining, :) - Z(s, :);
        distance = sum(delta .^ 2, 2);
        minDistance = min(minDistance, distance);
    end
    [~, localIdx] = max(minDistance);
    selected(end + 1) = remaining(localIdx); %#ok<AGROW>
    remaining(localIdx) = [];
end
end

function Z = standardizeColumns(X)
mu = mean(X, 1);
sigma = std(X, 0, 1);
sigma(~isfinite(sigma) | sigma < 1e-12) = 1;
Z = (X - mu) ./ sigma;
Z(~isfinite(Z)) = 0;
end

function X = flattenRparFeatures(models, edgeIdx)
numModels = numel(models);
numFeatures = 3 * numel(edgeIdx);
X = nan(numModels, numFeatures);
for i = 1:numModels
    X(i, :) = reshape(models(i).rpar(:, edgeIdx), 1, []);
end
end

function signature = buildScreenSignature(cfg, modelDir, spines, donors)
signature = struct();
signature.version = cfg.version;
signature.constantRuns = [spines.run];
signature.perturbedRuns = [donors.run];
signature.constantRpar = cat(3, spines.rpar);
signature.perturbedRpar = cat(3, donors.rpar);
signature.candidateMode = cfg.candidateMode;
signature.validationDirectionThreshold = cfg.validationDirectionThreshold;
signature.perturbValue = cfg.perturbValue;
signature.inputValue = cfg.inputValue;
signature.stressInputValue = cfg.stressInputValue;
signature.activeInputs = cfg.activeInputs;
signature.screenMaxPrunePasses = cfg.screenMaxPrunePasses;
signature.validationFile = fileSignature(fullfile(modelDir, 'CreateValidationMatrix.m'));
signature.odeFile = fileSignature(fullfile(modelDir, ['ODElist_EC_Network_', cfg.version, '.mat']));
signature.refFile = fileSignature(fullfile(modelDir, ['refcase_EC_Network_', cfg.version, '.mat']));
end

function signature = fileSignature(filename)
info = dir(filename);
signature = struct('name', info.name, 'bytes', info.bytes, 'datenum', info.datenum);
end

function screenResults = runOrLoadPairScreening(spines, donors, model, cfg, screenSignature)
checkpointFile = fullfile(cfg.outputDir, 'screening_checkpoint.mat');
finalFile = fullfile(cfg.outputDir, 'screening_results.mat');
numPairs = numel(spines) * numel(donors);
screenResults = cell(numPairs, 1);

if cfg.reuseExistingScreening
    loadFile = '';
    if exist(finalFile, 'file') == 2
        loadFile = finalFile;
    elseif exist(checkpointFile, 'file') == 2
        loadFile = checkpointFile;
    end
    if ~isempty(loadFile)
        S = load(loadFile, 'screenResults', 'savedScreenSignature');
        if isfield(S, 'savedScreenSignature') && isequaln(S.savedScreenSignature, screenSignature) ...
                && numel(S.screenResults) == numPairs
            screenResults = S.screenResults;
            fprintf('Loaded compatible screening results from %s.\n', loadFile);
        else
            warning('Existing screening results do not match current inputs and will be ignored.');
        end
    end
end

if ~cfg.runPairScreening
    assert(all(~cellfun(@isempty, screenResults)), ...
        'runPairScreening is false, but complete compatible screening results were not loaded.');
    return;
end

constantIndex = repelem((1:numel(spines))', numel(donors));
donorIndex = repmat((1:numel(donors))', numel(spines), 1);
parallelScreening = cfg.useParallelPairScreening && ...
    license('test', 'Distrib_Computing_Toolbox');
parallelPoolChecked = false;
batchSize = max(1, round(cfg.screenCheckpointEvery));

for batchStart = 1:batchSize:numPairs
    batchEnd = min(numPairs, batchStart + batchSize - 1);
    batchPairs = (batchStart:batchEnd)';
    pendingPairs = batchPairs(cellfun(@isempty, screenResults(batchPairs)));
    if isempty(pendingPairs)
        continue;
    end

    fprintf('Screening pair batch %d-%d of %d (%d pending pairs).\n', ...
        batchStart, batchEnd, numPairs, numel(pendingPairs));

    if parallelScreening && ~parallelPoolChecked
        parallelScreening = startParallelPool(cfg.parallelWorkers);
        parallelPoolChecked = true;
    end

    batchResults = cell(numel(pendingPairs), 1);
    batchSpines = cell(numel(pendingPairs), 1);
    batchDonors = cell(numel(pendingPairs), 1);
    for k = 1:numel(pendingPairs)
        pairIndex = pendingPairs(k);
        batchSpines{k} = spines(constantIndex(pairIndex));
        batchDonors{k} = donors(donorIndex(pairIndex));
    end

    if parallelScreening
        parfor k = 1:numel(pendingPairs)
            search = buildSearchSettings(cfg, 'screen', 0);
            batchResults{k} = optimizeHybridPair( ...
                batchSpines{k}, batchDonors{k}, model, cfg, search, []);
        end
    else
        for k = 1:numel(pendingPairs)
            pairIndex = pendingPairs(k);
            fprintf('  Pair %d: constant run %g, perturbed run %g\n', ...
                pairIndex, batchSpines{k}.run, batchDonors{k}.run);
            search = buildSearchSettings(cfg, 'screen', 0);
            batchResults{k} = optimizeHybridPair( ...
                batchSpines{k}, batchDonors{k}, model, cfg, search, []);
        end
    end

    for k = 1:numel(pendingPairs)
        pairIndex = pendingPairs(k);
        batchResults{k}.pairIndex = pairIndex;
        screenResults{pairIndex} = batchResults{k};
    end

    savedScreenSignature = screenSignature;
    atomicSaveMat(checkpointFile, struct( ...
        'screenResults', {screenResults}, ...
        'savedScreenSignature', savedScreenSignature));
end

savedScreenSignature = screenSignature;
atomicSaveMat(finalFile, struct( ...
    'screenResults', {screenResults}, ...
    'savedScreenSignature', savedScreenSignature));
end

function search = buildSearchSettings(cfg, mode, seed)
search = struct();
search.mode = mode;
search.rngSeed = seed;
search.populationSize = cfg.populationSize;
search.numGenerations = cfg.numGenerations;
search.eliteCount = cfg.eliteCount;
search.tournamentSize = cfg.tournamentSize;
search.mutationRate = cfg.mutationRate;
search.complexityPenalty = cfg.complexityPenalty;
search.runGreedyForward = cfg.runGreedyForward;
search.runGreedyForwardAfterPrune = cfg.runGreedyForwardAfterPrune;
search.useParallelPopulationEvaluation = cfg.useParallelPopulationEvaluation;
search.parallelWorkers = cfg.parallelWorkers;

switch lower(mode)
    case 'screen'
        search.runGeneticSearch = false;
        search.runGreedyPrune = true;
        search.maxPrunePasses = cfg.screenMaxPrunePasses;
    case 'final'
        search.runGeneticSearch = true;
        search.runGreedyPrune = cfg.runGreedyPrune;
        search.maxPrunePasses = cfg.finalMaxPrunePasses;
    otherwise
        error('Unknown search mode: %s', mode);
end
end

function result = optimizeHybridPair(constantModel, donorModel, model, cfg, search, seedFullMask)
rng(search.rngSeed, 'twister');

candidateEdges = findCandidateEdges( ...
    model.EdgeTable, constantModel.rpar, donorModel.rpar, cfg.candidateMode);
numCandidates = numel(candidateEdges);
scoreCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
numUniqueEvaluations = 0;
generationHistory = table();
parallelEnabled = search.useParallelPopulationEvaluation && ...
    license('test', 'Distrib_Computing_Toolbox');
parallelPoolChecked = false;

if numCandidates == 0
    evalResult = evaluateRparDirect(constantModel.rpar, model, cfg, search.complexityPenalty, 0);
    constantEval = evalResult;
    donorEval = evalResult;
    result = assembleResult(false(1, model.numEdges), constantModel.rpar, evalResult);
    return;
end

constantMask = false(1, numCandidates);
donorMask = true(1, numCandidates);
constantEval = evaluateMask(constantMask);
donorEval = evaluateMask(donorMask);

bestMask = chooseBetterMask(constantMask, donorMask);
if ~isempty(seedFullMask)
    seedMask = logical(seedFullMask(candidateEdges));
    bestMask = chooseBetterMask(bestMask, seedMask);
else
    seedMask = bestMask;
end

if strcmpi(search.mode, 'screen')
    prunedDonorMask = greedyPrune(donorMask);
    bestMask = chooseBetterMask(bestMask, prunedDonorMask);
else
    if search.runGreedyForward
        forwardMask = greedyForward(constantMask);
        bestMask = chooseBetterMask(bestMask, forwardMask);
    end
    if search.runGeneticSearch
        gaMask = runBinaryGA(seedMask);
        bestMask = chooseBetterMask(bestMask, gaMask);
    end
    if search.runGreedyPrune
        bestMask = greedyPrune(bestMask);
    end
    if search.runGreedyForwardAfterPrune
        bestMask = greedyForward(bestMask);
        bestMask = greedyPrune(bestMask);
    end
end

bestEval = evaluateMask(bestMask);
bestRpar = buildHybridRpar(bestMask);
selectedFull = false(1, model.numEdges);
selectedFull(candidateEdges(bestMask)) = true;
result = assembleResult(selectedFull, bestRpar, bestEval);

    function evalResult = evaluateMask(mask)
        key = maskToKey(mask);
        if isKey(scoreCache, key)
            evalResult = scoreCache(key);
            return;
        end
        thisRpar = buildHybridRpar(mask);
        evalResult = evaluateRparDirect( ...
            thisRpar, model, cfg, search.complexityPenalty, nnz(mask));
        scoreCache(key) = evalResult;
        numUniqueEvaluations = numUniqueEvaluations + 1;
    end

    function rpar = buildHybridRpar(mask)
        rpar = constantModel.rpar;
        edgesToSwap = candidateEdges(mask);
        rpar(:, edgesToSwap) = donorModel.rpar(:, edgesToSwap);
    end

    function betterMask = chooseBetterMask(maskA, maskB)
        evalA = evaluateMask(maskA);
        evalB = evaluateMask(maskB);
        if isBetterEval(evalB, evalA)
            betterMask = maskB;
        else
            betterMask = maskA;
        end
    end

    function outMask = greedyForward(startMask)
        outMask = startMask;
        improved = true;
        while improved
            improved = false;
            currentEval = evaluateMask(outMask);
            bestTrialMask = outMask;
            bestTrialEval = currentEval;
            for c = find(~outMask)
                trialMask = outMask;
                trialMask(c) = true;
                trialEval = evaluateMask(trialMask);
                if isBetterEval(trialEval, bestTrialEval)
                    bestTrialMask = trialMask;
                    bestTrialEval = trialEval;
                end
            end
            if isBetterEval(bestTrialEval, currentEval)
                outMask = bestTrialMask;
                improved = true;
            end
        end
    end

    function outMask = greedyPrune(startMask)
        outMask = startMask;
        changed = true;
        prunePass = 0;
        while changed && prunePass < search.maxPrunePasses
            prunePass = prunePass + 1;
            changed = false;
            currentEval = evaluateMask(outMask);
            for c = find(outMask)
                trialMask = outMask;
                trialMask(c) = false;
                trialEval = evaluateMask(trialMask);
                if ~trialEval.failed && trialEval.nMatch >= currentEval.nMatch
                    outMask = trialMask;
                    currentEval = trialEval;
                    changed = true;
                end
            end
        end
    end

    function bestGaMask = runBinaryGA(initialMask)
        pop = initializePopulation(initialMask);
        historyGeneration = (1:(search.numGenerations + 1))';
        historyMatches = nan(search.numGenerations + 1, 1);
        historyScore = nan(search.numGenerations + 1, 1);
        historySelected = nan(search.numGenerations + 1, 1);

        for gen = 1:search.numGenerations
            [pop, evals, ~] = evaluateAndSortPopulation(pop);
            historyMatches(gen) = evals{1}.nMatch;
            historyScore(gen) = evals{1}.scorePct;
            historySelected(gen) = evals{1}.numSelected;

            if cfg.verbose && (gen == 1 || mod(gen, max(1, round(search.numGenerations / 5))) == 0)
                fprintf('    GA %d/%d: %.2f%%, %d selected edges\n', ...
                    gen, search.numGenerations, evals{1}.scorePct, evals{1}.numSelected);
            end

            nextPop = pop(1:search.eliteCount, :);
            while size(nextPop, 1) < search.populationSize
                parentA = tournamentSelect(pop, evals, search.tournamentSize);
                parentB = tournamentSelect(pop, evals, search.tournamentSize);
                child = crossoverMasks(parentA, parentB);
                child = mutateMask(child, search.mutationRate);
                nextPop(end + 1, :) = child; %#ok<AGROW>
            end
            pop = nextPop;
        end

        [pop, evals, ~] = evaluateAndSortPopulation(pop);
        historyMatches(end) = evals{1}.nMatch;
        historyScore(end) = evals{1}.scorePct;
        historySelected(end) = evals{1}.numSelected;
        generationHistory = table(historyGeneration, historyMatches, historyScore, ...
            historySelected, 'VariableNames', ...
            {'EvaluatedPopulation','BestValidationMatches','BestValidationAgreementPct', ...
            'BestSelectedPerturbedEdges'});
        bestGaMask = pop(1, :);
    end

    function [pop, evals, fitness] = evaluateAndSortPopulation(pop)
        numPopulation = size(pop, 1);
        keys = cell(numPopulation, 1);
        for p = 1:numPopulation
            keys{p} = maskToKey(pop(p, :));
        end

        [uniqueKeys, uniqueRows] = unique(keys, 'stable');
        missing = false(numel(uniqueKeys), 1);
        for u = 1:numel(uniqueKeys)
            missing(u) = ~isKey(scoreCache, uniqueKeys{u});
        end

        missingRows = uniqueRows(missing);
        missingKeys = uniqueKeys(missing);
        numMissing = numel(missingRows);
        if numMissing > 0
            rparBatch = nan(3, model.numEdges, numMissing);
            selectedBatch = zeros(numMissing, 1);
            for u = 1:numMissing
                mask = pop(missingRows(u), :);
                rparBatch(:, :, u) = buildHybridRpar(mask);
                selectedBatch(u) = nnz(mask);
            end

            if parallelEnabled && ~parallelPoolChecked
                parallelEnabled = startParallelPool(search.parallelWorkers);
                parallelPoolChecked = true;
            end

            batchEvals = cell(numMissing, 1);
            complexityPenalty = search.complexityPenalty;
            if parallelEnabled
                parfor u = 1:numMissing
                    batchEvals{u} = evaluateRparDirect( ...
                        rparBatch(:, :, u), model, cfg, ...
                        complexityPenalty, selectedBatch(u));
                end
            else
                for u = 1:numMissing
                    batchEvals{u} = evaluateRparDirect( ...
                        rparBatch(:, :, u), model, cfg, ...
                        complexityPenalty, selectedBatch(u));
                end
            end

            for u = 1:numMissing
                scoreCache(missingKeys{u}) = batchEvals{u};
            end
            numUniqueEvaluations = numUniqueEvaluations + numMissing;
        end

        evals = cell(numPopulation, 1);
        fitness = -inf(numPopulation, 1);
        for p = 1:numPopulation
            evals{p} = scoreCache(keys{p});
            fitness(p) = evals{p}.fitness;
        end
        [~, order] = sort(fitness, 'descend');
        pop = pop(order, :);
        evals = evals(order);
        fitness = fitness(order);
    end

    function pop = initializePopulation(initialMask)
        pop = false(search.populationSize, numCandidates);
        row = 1;
        pop(row, :) = false(1, numCandidates);
        row = row + 1;
        if row <= search.populationSize
            pop(row, :) = true(1, numCandidates);
            row = row + 1;
        end
        if row <= search.populationSize
            pop(row, :) = initialMask;
            row = row + 1;
        end
        sparseRates = [0.02, 0.05, 0.10, 0.20, 0.35];
        while row <= search.populationSize
            rate = sparseRates(mod(row - 1, numel(sparseRates)) + 1);
            pop(row, :) = rand(1, numCandidates) < rate;
            row = row + 1;
        end
    end

    function resultOut = assembleResult(selectedFull, bestRpar, bestEval)
        resultOut = struct();
        resultOut.mode = search.mode;
        resultOut.constantRun = constantModel.run;
        resultOut.perturbedRun = donorModel.run;
        resultOut.constantSavedMatches = constantModel.savedMatches;
        resultOut.constantSavedScorePct = constantModel.savedScorePct;
        resultOut.perturbedSavedMatches = donorModel.savedMatches;
        resultOut.perturbedSavedScorePct = donorModel.savedScorePct;
        resultOut.rngSeed = search.rngSeed;
        resultOut.populationSize = search.populationSize;
        resultOut.numGenerations = search.numGenerations;
        resultOut.numCandidates = numCandidates;
        resultOut.candidateEdges = candidateEdges;
        resultOut.selectedFull = selectedFull;
        resultOut.bestRpar = bestRpar;
        resultOut.bestEval = bestEval;
        resultOut.constantEval = constantEval;
        resultOut.perturbedEval = donorEval;
        resultOut.generationHistory = generationHistory;
        resultOut.numUniqueModelEvaluations = numUniqueEvaluations;
        resultOut.numOdeIntegrationsApprox = numUniqueEvaluations * ...
            (1 + numel(model.validationPerturbIdx));
        resultOut.searchSettings = search;
        resultOut.pairIndex = NaN;
        resultOut.outputDir = '';
    end
end

function evalResult = evaluateRparDirect(rpar, model, cfg, complexityPenalty, numSelected)
try
    [scorePct, nMatch, nKnown, DiffStore, QualStore, MatchStore] = ...
        scoreExperimentalValidation(rpar, model, cfg);
    failed = false;
    errorMessage = '';
catch ME
    if cfg.stopOnODEError
        rethrow(ME);
    end
    failed = true;
    errorMessage = ME.message;
    scorePct = -Inf;
    nMatch = -Inf;
    nKnown = nnz(model.knownValidationMask);
    DiffStore = nan(size(model.ExpValidation));
    QualStore = nan(size(model.ExpValidation));
    MatchStore = false(size(model.ExpValidation));
end

evalResult = struct();
evalResult.failed = failed;
evalResult.errorMessage = errorMessage;
evalResult.scorePct = scorePct;
evalResult.nMatch = nMatch;
evalResult.nKnown = nKnown;
evalResult.DiffStore = DiffStore;
evalResult.QualStore = QualStore;
evalResult.MatchStore = MatchStore;
evalResult.numSelected = numSelected;
if failed
    evalResult.fitness = -Inf;
else
    evalResult.fitness = nMatch - complexityPenalty * numSelected;
end
end

function [scorePct, nMatch, nKnown, DiffStore, QualStore, MatchStore] = ...
        scoreExperimentalValidation(rpar, model, cfg)
[~, yBaseline] = ODE_master(rpar, model.y0Base, model.speciesNames, ...
    model.tau, model.ymax, model.numEdges, model.ODEfilename);
if any(~isfinite(yBaseline(:)))
    error('Non-finite values were produced by the baseline ODE simulation.');
end
yBaselineFinal = yBaseline(end, :);

numPert = numel(model.validationPerturbIdx);
numSpecies = numel(model.speciesNames);
DiffStore = nan(numPert, numSpecies);

for p = 1:numPert
    y0Perturbed = yBaselineFinal;
    idx = model.validationPerturbIdx(p);
    y0Perturbed(idx) = y0Perturbed(idx) + cfg.perturbValue;
    [~, yPerturbed] = ODE_master(rpar, y0Perturbed, model.speciesNames, ...
        model.tau, model.ymax, model.numEdges, model.ODEfilename);
    if any(~isfinite(yPerturbed(:)))
        error('Non-finite values were produced by perturbed ODE simulation %d.', p);
    end
    DiffStore(p, :) = yPerturbed(end, :) - yBaselineFinal;
end

DiffForSign = DiffStore;
DiffForSign(abs(DiffForSign) < cfg.validationDirectionThreshold) = 0;
QualStore = sign(DiffForSign);
MatchStore = (QualStore == model.ExpValidation) & model.knownValidationMask;
nMatch = nnz(MatchStore);
nKnown = nnz(model.knownValidationMask);
scorePct = 100 * nMatch / nKnown;
end

function tf = isBetterEval(evalA, evalB)
if evalA.failed ~= evalB.failed
    tf = ~evalA.failed;
elseif evalA.nMatch ~= evalB.nMatch
    tf = evalA.nMatch > evalB.nMatch;
elseif evalA.numSelected ~= evalB.numSelected
    tf = evalA.numSelected < evalB.numSelected;
else
    tf = evalA.scorePct > evalB.scorePct;
end
end

function parent = tournamentSelect(pop, evals, tournamentSize)
competitors = randi(size(pop, 1), [tournamentSize, 1]);
bestIdx = competitors(1);
for i = 2:numel(competitors)
    idx = competitors(i);
    if isBetterEval(evals{idx}, evals{bestIdx})
        bestIdx = idx;
    end
end
parent = pop(bestIdx, :);
end

function child = crossoverMasks(parentA, parentB)
useA = rand(size(parentA)) < 0.5;
child = parentA;
child(~useA) = parentB(~useA);
end

function mask = mutateMask(mask, mutationRate)
flip = rand(size(mask)) < mutationRate;
mask(flip) = ~mask(flip);
end

function key = maskToKey(mask)
key = char('0' + mask);
end

function candidateEdges = findCandidateEdges(EdgeTable, constantRpar, donorRpar, candidateMode)
changed = any(abs(constantRpar - donorRpar) > 1e-12, 1);
switch lower(candidateMode)
    case 'used'
        used = logical(EdgeTable.UsedInODE(:)');
        candidateEdges = find(used & changed);
    case 'all'
        candidateEdges = find(changed);
    otherwise
        error('Unknown candidateMode: %s.', candidateMode);
end
end

function rows = selectFinalistRows(ScreeningTable, numFinalists, diversifySources)
validRows = find(~ScreeningTable.Failed);
assert(~isempty(validRows), 'Every screened spine/donor pair failed.');
numToSelect = min(round(numFinalists), numel(validRows));
if ~diversifySources || numToSelect == 1
    rows = validRows(1:numToSelect);
    return;
end

bestMatches = ScreeningTable.BestValidationMatches(validRows(1));
topRows = validRows( ...
    ScreeningTable.BestValidationMatches(validRows) == bestMatches);
rows = topRows(1);

while numel(rows) < min(numToSelect, numel(topRows))
    candidates = setdiff(topRows, rows, 'stable');
    usedConstants = ScreeningTable.ConstantRun(rows);
    usedDonors = ScreeningTable.PerturbedRun(rows);
    addsConstant = ~ismember(ScreeningTable.ConstantRun(candidates), usedConstants);
    addsDonor = ~ismember(ScreeningTable.PerturbedRun(candidates), usedDonors);
    coverageGain = double(addsConstant) + double(addsDonor);
    candidateComplexity = ScreeningTable.SelectedPerturbedEdges(candidates);
    originalRank = (1:numel(candidates))';
    [~, order] = sortrows( ...
        [-coverageGain, candidateComplexity, originalRank], [1, 2, 3]);
    rows(end + 1, 1) = candidates(order(1)); %#ok<AGROW>
end

if numel(rows) < numToSelect
    remaining = setdiff(validRows, rows, 'stable');
    rows = [rows; remaining(1:(numToSelect - numel(rows)))];
end
end

function T = buildPairResultTable(results)
n = numel(results);
PairIndex = nan(n, 1);
ConstantRun = nan(n, 1);
PerturbedRun = nan(n, 1);
RngSeed = nan(n, 1);
Failed = true(n, 1);
SavedConstantMatches = nan(n, 1);
SavedPerturbedMatches = nan(n, 1);
RecomputedConstantMatches = nan(n, 1);
RecomputedPerturbedMatches = nan(n, 1);
BestValidationMatches = nan(n, 1);
ValidationKnownChecks = nan(n, 1);
BestValidationAgreementPct = nan(n, 1);
CandidateEdges = nan(n, 1);
SelectedPerturbedEdges = nan(n, 1);
UniqueModelEvaluations = nan(n, 1);
ApproxOdeIntegrations = nan(n, 1);
OutputDir = cell(n, 1);

for i = 1:n
    R = results{i};
    PairIndex(i) = R.pairIndex;
    ConstantRun(i) = R.constantRun;
    PerturbedRun(i) = R.perturbedRun;
    RngSeed(i) = R.rngSeed;
    Failed(i) = R.bestEval.failed;
    SavedConstantMatches(i) = R.constantSavedMatches;
    SavedPerturbedMatches(i) = R.perturbedSavedMatches;
    RecomputedConstantMatches(i) = R.constantEval.nMatch;
    RecomputedPerturbedMatches(i) = R.perturbedEval.nMatch;
    BestValidationMatches(i) = R.bestEval.nMatch;
    ValidationKnownChecks(i) = R.bestEval.nKnown;
    BestValidationAgreementPct(i) = R.bestEval.scorePct;
    CandidateEdges(i) = R.numCandidates;
    SelectedPerturbedEdges(i) = R.bestEval.numSelected;
    UniqueModelEvaluations(i) = R.numUniqueModelEvaluations;
    ApproxOdeIntegrations(i) = R.numOdeIntegrationsApprox;
    OutputDir{i} = R.outputDir;
end

T = table(PairIndex, ConstantRun, PerturbedRun, RngSeed, Failed, ...
    SavedConstantMatches, SavedPerturbedMatches, RecomputedConstantMatches, ...
    RecomputedPerturbedMatches, BestValidationMatches, ValidationKnownChecks, ...
    BestValidationAgreementPct, CandidateEdges, SelectedPerturbedEdges, ...
    UniqueModelEvaluations, ApproxOdeIntegrations, OutputDir);
end

function T = sortPairResultTable(T)
T = sortrows(T, {'Failed','BestValidationMatches','SelectedPerturbedEdges', ...
    'ConstantRun','PerturbedRun','RngSeed'}, ...
    {'ascend','descend','ascend','ascend','ascend','ascend'});
end

function finalResults = runFinalistSearches(FinalistTable, screenResults, ...
        spines, donors, model, cfg, screenSignature)
numRuns = height(FinalistTable) * numel(cfg.rngSeeds);
finalResults = cell(numRuns, 1);
finalRoot = fullfile(cfg.outputDir, 'final_runs');
ensureDir(finalRoot);

runIndex = 0;
for f = 1:height(FinalistTable)
    constantRun = FinalistTable.ConstantRun(f);
    perturbedRun = FinalistTable.PerturbedRun(f);
    spineIdx = find([spines.run] == constantRun, 1);
    donorIdx = find([donors.run] == perturbedRun, 1);
    screenIdx = FinalistTable.PairIndex(f);
    screenResult = screenResults{screenIdx};

    for s = 1:numel(cfg.rngSeeds)
        runIndex = runIndex + 1;
        seed = cfg.rngSeeds(s);
        runDir = fullfile(finalRoot, sprintf( ...
            'constant_%06g_perturbed_%06g_seed_%06g', constantRun, perturbedRun, seed));
        resultFile = fullfile(runDir, 'final_run_results.mat');
        finalSignature = buildFinalSignature( ...
            screenSignature, constantRun, perturbedRun, seed, cfg);

        loaded = false;
        if cfg.resumeCompletedFinalRuns && exist(resultFile, 'file') == 2
            try
                S = load(resultFile, 'pairResult', 'savedFinalSignature');
                if isfield(S, 'savedFinalSignature') && ...
                        isequaln(S.savedFinalSignature, finalSignature)
                    finalResults{runIndex} = S.pairResult;
                    loaded = true;
                    fprintf('Loaded completed final run %d/%d from %s.\n', ...
                        runIndex, numRuns, runDir);
                end
            catch ME
                quarantineFile = quarantineUnreadableMat(resultFile);
                warning('EC:UnreadableFinalCheckpoint', ...
                    ['Could not read final checkpoint %s (%s). ', ...
                    'It was moved to %s and this run will be recomputed.'], ...
                    resultFile, ME.message, quarantineFile);
            end
        end

        if loaded
            continue;
        end

        fprintf('\nFinal GA run %d/%d: constant %g, perturbed %g, seed %g\n', ...
            runIndex, numRuns, constantRun, perturbedRun, seed);
        search = buildSearchSettings(cfg, 'final', seed);
        pairResult = optimizeHybridPair(spines(spineIdx), donors(donorIdx), ...
            model, cfg, search, screenResult.selectedFull);
        pairResult.pairIndex = screenIdx;
        pairResult.outputDir = runDir;
        finalResults{runIndex} = pairResult;

        ensureDir(runDir);
        savedFinalSignature = finalSignature;
        atomicSaveMat(resultFile, struct( ...
            'pairResult', pairResult, ...
            'savedFinalSignature', savedFinalSignature));
        if cfg.saveIndividualFinalRunFiles
            writePairOutputs(runDir, pairResult, model.EdgeTable, model, cfg, ...
                cfg.saveIndividualValidationPlots);
        end
    end
end
end

function signature = buildFinalSignature(screenSignature, constantRun, perturbedRun, seed, cfg)
signature = struct();
signature.screenSignature = screenSignature;
signature.constantRun = constantRun;
signature.perturbedRun = perturbedRun;
signature.seed = seed;
signature.populationSize = cfg.populationSize;
signature.numGenerations = cfg.numGenerations;
signature.eliteCount = cfg.eliteCount;
signature.tournamentSize = cfg.tournamentSize;
signature.mutationRate = cfg.mutationRate;
signature.complexityPenalty = cfg.complexityPenalty;
signature.finalMaxPrunePasses = cfg.finalMaxPrunePasses;
signature.runGreedyForward = cfg.runGreedyForward;
signature.runGreedyPrune = cfg.runGreedyPrune;
signature.runGreedyForwardAfterPrune = cfg.runGreedyForwardAfterPrune;
signature.useParallelPopulationEvaluation = cfg.useParallelPopulationEvaluation;
signature.parallelWorkers = cfg.parallelWorkers;
end

function enabled = startParallelPool(workerCount)
enabled = false;
try
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(workerCount)
            pool = parpool('local');
        else
            pool = parpool('local', round(workerCount));
        end
    elseif ~isempty(workerCount) && pool.NumWorkers ~= round(workerCount)
        warning(['An existing parallel pool has %d workers; requested %d. ', ...
            'The existing pool will be used.'], pool.NumWorkers, round(workerCount));
    end
    fprintf('    Parallel evaluation enabled with %d workers.\n', pool.NumWorkers);
    enabled = true;
catch ME
    warning('EC:ParallelUnavailable', '%s', ...
        ['Parallel population scoring is unavailable; continuing serially: ', ME.message]);
end
end

function representativeIdx = chooseRepresentativeResult(results)
matches = cellfun(@(r) r.bestEval.nMatch, results);
selected = cellfun(@(r) r.bestEval.numSelected, results);
failed = cellfun(@(r) r.bestEval.failed, results);
maxMatches = max(matches(~failed));
candidates = find(~failed & matches == maxMatches);
minSelected = min(selected(candidates));
candidates = candidates(selected(candidates) == minSelected);

if isscalar(candidates)
    representativeIdx = candidates;
    return;
end

models = repmat(emptySourceModel(), numel(candidates), 1);
for i = 1:numel(candidates)
    models(i).rpar = results{candidates(i)}.bestRpar;
end
X = flattenRparFeatures(models, 1:size(models(1).rpar, 2));
Z = standardizeColumns(X);
distanceSums = zeros(numel(candidates), 1);
for i = 1:numel(candidates)
    delta = Z - Z(i, :);
    distanceSums(i) = sum(sqrt(sum(delta .^ 2, 2)));
end
[~, localIdx] = min(distanceSums);
representativeIdx = candidates(localIdx);
end

function [EdgeOut, stats] = buildEnsembleEdgeTable(results, EdgeTable)
numRuns = numel(results);
numEdges = height(EdgeTable);
selectedMatrix = false(numRuns, numEdges);
rparStack = nan(3, numEdges, numRuns);
matches = nan(numRuns, 1);
selectedCounts = nan(numRuns, 1);
failed = false(numRuns, 1);

for i = 1:numRuns
    selectedMatrix(i, :) = results{i}.selectedFull;
    rparStack(:, :, i) = results{i}.bestRpar;
    matches(i) = results{i}.bestEval.nMatch;
    selectedCounts(i) = results{i}.bestEval.numSelected;
    failed(i) = results{i}.bestEval.failed;
end

valid = ~failed;
maxMatches = max(matches(valid));
maxAccuracy = valid & matches == maxMatches;
minSelectedAtMax = min(selectedCounts(maxAccuracy));
bestParsimonious = maxAccuracy & selectedCounts == minSelectedAtMax;

EdgeOut = EdgeTable;
EdgeOut.SelectionCountAllFinalRuns = sum(selectedMatrix(valid, :), 1)';
EdgeOut.SelectionFrequencyAllFinalRuns = mean(selectedMatrix(valid, :), 1)';
EdgeOut.SelectionCountMaxAccuracy = sum(selectedMatrix(maxAccuracy, :), 1)';
EdgeOut.SelectionFrequencyMaxAccuracy = mean(selectedMatrix(maxAccuracy, :), 1)';
EdgeOut.SelectionCountBestParsimonious = sum(selectedMatrix(bestParsimonious, :), 1)';
EdgeOut.SelectionFrequencyBestParsimonious = mean(selectedMatrix(bestParsimonious, :), 1)';

labels = {'W','HillN','EC50'};
for p = 1:3
    values = reshape(permute(rparStack(p, :, maxAccuracy), [3 2 1]), ...
        nnz(maxAccuracy), numEdges);
    [minimum, med, maximum] = columnRange(values);
    EdgeOut.([labels{p}, 'MinMaxAccuracy']) = minimum';
    EdgeOut.([labels{p}, 'MedianMaxAccuracy']) = med';
    EdgeOut.([labels{p}, 'MaxMaxAccuracy']) = maximum';
end

stats = struct();
stats.numFinalRuns = numRuns;
stats.numValidFinalRuns = nnz(valid);
stats.maxValidationMatches = maxMatches;
stats.numMaxAccuracyRuns = nnz(maxAccuracy);
stats.minSelectedEdgesAtMaxAccuracy = minSelectedAtMax;
stats.numBestParsimoniousRuns = nnz(bestParsimonious);
stats.maxAccuracyMask = maxAccuracy;
stats.bestParsimoniousMask = bestParsimonious;
end

function [minimum, med, maximum] = columnRange(X)
minimum = min(X, [], 1);
med = median(X, 1);
maximum = max(X, [], 1);
end

function T = summarizeSourceOutcomes(results, sourceType)
switch lower(sourceType)
    case 'constant'
        ids = cellfun(@(r) r.constantRun, results);
        variableName = 'ConstantRun';
    case 'perturbed'
        ids = cellfun(@(r) r.perturbedRun, results);
        variableName = 'PerturbedRun';
    otherwise
        error('Unknown source type: %s', sourceType);
end

uniqueIds = unique(ids(:), 'stable');
n = numel(uniqueIds);
NumberFinalRuns = zeros(n, 1);
BestValidationMatches = nan(n, 1);
BestValidationAgreementPct = nan(n, 1);
MinimumSelectedEdgesAtBestScore = nan(n, 1);
MeanValidationAgreementPct = nan(n, 1);

for i = 1:n
    rows = find(ids == uniqueIds(i));
    scores = cellfun(@(r) r.bestEval.scorePct, results(rows));
    matches = cellfun(@(r) r.bestEval.nMatch, results(rows));
    selected = cellfun(@(r) r.bestEval.numSelected, results(rows));
    NumberFinalRuns(i) = numel(rows);
    BestValidationMatches(i) = max(matches);
    bestRows = matches == BestValidationMatches(i);
    BestValidationAgreementPct(i) = max(scores(bestRows));
    MinimumSelectedEdgesAtBestScore(i) = min(selected(bestRows));
    MeanValidationAgreementPct(i) = mean(scores);
end

T = table(uniqueIds, NumberFinalRuns, BestValidationMatches, ...
    BestValidationAgreementPct, MinimumSelectedEdgesAtBestScore, ...
    MeanValidationAgreementPct, 'VariableNames', ...
    {variableName,'NumberFinalRuns','BestValidationMatches', ...
    'BestValidationAgreementPct','MinimumSelectedEdgesAtBestScore', ...
    'MeanValidationAgreementPct'});
T = sortrows(T, {'BestValidationMatches','MinimumSelectedEdgesAtBestScore'}, ...
    {'descend','ascend'});
end

function writePairOutputs(outputDir, result, EdgeTable, model, cfg, makePlot)
ensureDir(outputDir);
writematrix(result.bestRpar, ...
    fullfile(outputDir, 'best_constant_spine_perturbed_rpar_matrix_rows_w_n_EC50.csv'));

Summary = table(result.constantRun, result.perturbedRun, result.rngSeed, ...
    result.bestEval.failed, result.bestEval.scorePct, result.bestEval.nMatch, ...
    result.bestEval.nKnown, result.numCandidates, result.bestEval.numSelected, ...
    result.numUniqueModelEvaluations, result.numOdeIntegrationsApprox, ...
    'VariableNames', {'ConstantRun','PerturbedRun','RngSeed','Failed', ...
    'ValidationAgreementPct','ValidationMatches','ValidationKnownChecks', ...
    'CandidateEdges','SelectedPerturbedEdges','UniqueModelEvaluations', ...
    'ApproxOdeIntegrations'});
writetable(Summary, fullfile(outputDir, 'constant_spine_perturbed_summary.csv'));

EdgeOut = EdgeTable;
EdgeOut.SelectedPerturbed = result.selectedFull(:);
EdgeOut.FinalW = result.bestRpar(1, :)';
EdgeOut.FinalHillN = result.bestRpar(2, :)';
EdgeOut.FinalEC50 = result.bestRpar(3, :)';
writetable(EdgeOut, fullfile(outputDir, 'selected_perturbed_edges.csv'));

ValidationTable = buildValidationResultTable(result.bestEval, model, cfg);
writetable(ValidationTable, ...
    fullfile(outputDir, 'best_validation_direction_matrix_long.csv'));

if ~isempty(result.generationHistory)
    writetable(result.generationHistory, fullfile(outputDir, 'ga_generation_history.csv'));
end
if makePlot
    makeValidationPlot(outputDir, result, model, cfg);
end
end

function T = buildValidationResultTable(evalResult, model, cfg)
numPert = size(evalResult.DiffStore, 1);
numSpecies = size(evalResult.DiffStore, 2);
N = numPert * numSpecies;
Perturbation = cell(N, 1);
Species = cell(N, 1);
ExpectedDirection = nan(N, 1);
ObservedDifference = nan(N, 1);
ObservedDirection = nan(N, 1);
KnownLiteratureCheck = false(N, 1);
Matched = false(N, 1);

row = 0;
for p = 1:numPert
    for j = 1:numSpecies
        row = row + 1;
        Perturbation{row} = cfg.validationPerturbSpecies{p};
        Species{row} = model.speciesNames{j};
        ExpectedDirection(row) = model.ExpValidation(p, j);
        ObservedDifference(row) = evalResult.DiffStore(p, j);
        ObservedDirection(row) = evalResult.QualStore(p, j);
        KnownLiteratureCheck(row) = model.knownValidationMask(p, j);
        Matched(row) = evalResult.MatchStore(p, j);
    end
end
T = table(Perturbation, Species, ExpectedDirection, ObservedDifference, ...
    ObservedDirection, KnownLiteratureCheck, Matched);
end

function makeValidationPlot(outputDir, result, model, cfg)
plotIdx = model.validationSpeciesIdx;
plotNames = model.speciesNames(plotIdx);
numPlotSpecies = numel(plotIdx);
numPert = numel(cfg.validationPerturbSpecies);

fig = figure('Name', 'Representative constant-spine hybrid', 'Color', 'w');
imagesc(result.bestEval.DiffStore(:, plotIdx));
xticks(1:numPlotSpecies);
xticklabels(plotNames);
xtickangle(270);
yticks(1:numPert);
yticklabels(cfg.validationPerturbSpecies);
clim([-0.5, 1]);
if exist('custom_cmap.m', 'file') == 2
    colormap custom_cmap;
else
    colormap parula;
end
colorbar;
set(gca, 'TickLength', [0 0], 'FontSize', 8);
title(sprintf('Hybrid %s: %.2f%% agreement, %d perturbed edges', ...
    cfg.version, result.bestEval.scorePct, result.bestEval.numSelected), ...
    'Interpreter', 'none');
hold on;

for i = 1:numPlotSpecies
    plot([i - 0.5, i - 0.5], [0.5, numPert + 0.5], 'k-');
end
for i = 1:numPert
    plot([0.5, numPlotSpecies + 0.5], [i + 0.5, i + 0.5], 'k-');
end
for p = 1:numPert
    for j = 1:numPlotSpecies
        speciesIdx = plotIdx(j);
        if model.knownValidationMask(p, speciesIdx)
            if result.bestEval.QualStore(p, speciesIdx) == model.ExpValidation(p, speciesIdx)
                marker = char(hex2dec('2713'));
            else
                marker = char(hex2dec('2717'));
            end
            text(j - 0.25, p, marker, 'FontWeight', 'bold');
        end
    end
end

set(fig, 'Position', [100 100 1200 350]);
saveas(fig, fullfile(outputDir, 'best_validation_matrix_checks_x.png'));
savefig(fig, fullfile(outputDir, 'best_validation_matrix_checks_x.fig'));
if cfg.closeFiguresAfterSaving
    close(fig);
end
end

function overlay = makeEnsembleWssOverlay(results, representativeIdx, model, cfg)
matches = cellfun(@(r) r.bestEval.nMatch, results);
failed = cellfun(@(r) r.bestEval.failed, results);
maxMatches = max(matches(~failed));
eligible = find(~failed & matches == maxMatches);

if numel(eligible) > cfg.maxOverlayModels
    sourceModels = repmat(emptySourceModel(), numel(eligible), 1);
    for i = 1:numel(eligible)
        sourceModels(i).rpar = results{eligible(i)}.bestRpar;
    end
    features = flattenRparFeatures(sourceModels, 1:model.numEdges);
    representativeLocal = find(eligible == representativeIdx, 1);
    selectedLocal = selectDiverseRows(features, representativeLocal, cfg.maxOverlayModels);
    eligible = eligible(selectedLocal);
end
if ~ismember(representativeIdx, eligible)
    eligible(end + 1) = representativeIdx;
end

modelShearDyn = linspace(0, cfg.shearMaxDyn, cfg.numShearPoints);
modelWss = min(max(modelShearDyn ./ cfg.shearMaxDyn, 0), 1);
speciesToPlot = {'eNOS','NO','ppET1','ET1'};
speciesIdx = findSpeciesIndices(model.speciesNames, speciesToPlot, true);
assert(numel(speciesIdx) == numel(speciesToPlot), ...
    'Every WSS overlay species must exist in the model.');
wssIdx = find(strcmp(model.speciesNames, 'Wss'), 1);
stressIdx = find(strcmp(model.speciesNames, 'Stress'), 1);

numModels = numel(eligible);
foldValues = nan(cfg.numShearPoints, numel(speciesIdx), numModels);
for m = 1:numModels
    rpar = results{eligible(m)}.bestRpar;
    finalValues = nan(cfg.numShearPoints, numel(speciesIdx));
    for i = 1:numel(modelWss)
        y0 = zeros(1, numel(model.speciesNames));
        active = cfg.activeInputs(cfg.activeInputs <= numel(y0));
        y0(active) = cfg.inputValue;
        y0(stressIdx) = cfg.stressInputValue;
        y0(wssIdx) = modelWss(i);
        [~, y] = ODE_master(rpar, y0, model.speciesNames, model.tau, ...
            model.ymax, model.numEdges, model.ODEfilename);
        finalValues(i, :) = y(end, speciesIdx);
    end
    baseline = finalValues(1, :);
    for j = 1:numel(speciesIdx)
        foldValues(:, j, m) = safeFold(finalValues(:, j), baseline(j));
    end
end

representativeLocal = find(eligible == representativeIdx, 1);
literature = literatureAnchors();
fig = figure('Name', ['WSS ensemble overlay ', cfg.version], 'Color', 'w');
set(fig, 'Position', [120, 80, 1100, 760]);

plotOverlayPanel(1, modelShearDyn, squeeze(foldValues(:, 1, :)), ...
    representativeLocal, literature.uematsuShear, literature.uematsuENOS, ...
    'eNOS fold vs static', 'Uematsu 1995 ecNOS mRNA', cfg);
plotOverlayPanel(2, modelShearDyn, squeeze(foldValues(:, 2, :)), ...
    representativeLocal, literature.uematsuNOShear, literature.uematsuNO, ...
    'NO fold vs static', 'Uematsu 1995 NO/nitrite', cfg);
plotOverlayPanel(3, modelShearDyn, squeeze(foldValues(:, 3, :)), ...
    representativeLocal, literature.malekShear, literature.malekET1, ...
    'ppET1 / Edn1 fold vs static', 'Malek 1995 ET-1 mRNA', cfg);
plotOverlayPanel(4, modelShearDyn, squeeze(foldValues(:, 4, :)), ...
    representativeLocal, literature.malekShear, literature.malekET1, ...
    'ET1 fold vs static', 'Malek 1995 ET-1 mRNA', cfg);

sgtitle(sprintf('WSS response ensemble, EC Network %s, literature not scored', ...
    cfg.version), 'Interpreter', 'none');
figureBase = fullfile(cfg.outputDir, ['wss_ensemble_literature_overlay_', cfg.version]);
saveas(fig, [figureBase, '.png']);
savefig(fig, [figureBase, '.fig']);
if cfg.closeFiguresAfterSaving
    close(fig);
end

ModelIndex = [];
ConstantRun = [];
PerturbedRun = [];
RngSeed = [];
ShearDynCm2 = [];
ModelWss = [];
Model_eNOS_Fold = [];
Model_NO_Fold = [];
Model_ppET1_Fold = [];
Model_ET1_Fold = [];
IsRepresentative = [];
for m = 1:numModels
    n = cfg.numShearPoints;
    result = results{eligible(m)};
    ModelIndex = [ModelIndex; repmat(eligible(m), n, 1)]; %#ok<AGROW>
    ConstantRun = [ConstantRun; repmat(result.constantRun, n, 1)]; %#ok<AGROW>
    PerturbedRun = [PerturbedRun; repmat(result.perturbedRun, n, 1)]; %#ok<AGROW>
    RngSeed = [RngSeed; repmat(result.rngSeed, n, 1)]; %#ok<AGROW>
    ShearDynCm2 = [ShearDynCm2; modelShearDyn(:)]; %#ok<AGROW>
    ModelWss = [ModelWss; modelWss(:)]; %#ok<AGROW>
    Model_eNOS_Fold = [Model_eNOS_Fold; foldValues(:, 1, m)]; %#ok<AGROW>
    Model_NO_Fold = [Model_NO_Fold; foldValues(:, 2, m)]; %#ok<AGROW>
    Model_ppET1_Fold = [Model_ppET1_Fold; foldValues(:, 3, m)]; %#ok<AGROW>
    Model_ET1_Fold = [Model_ET1_Fold; foldValues(:, 4, m)]; %#ok<AGROW>
    IsRepresentative = [IsRepresentative; repmat(eligible(m) == representativeIdx, n, 1)]; %#ok<AGROW>
end
ModelTable = table(ModelIndex, ConstantRun, PerturbedRun, RngSeed, ...
    ShearDynCm2, ModelWss, Model_eNOS_Fold, Model_NO_Fold, ...
    Model_ppET1_Fold, Model_ET1_Fold, IsRepresentative);
writetable(ModelTable, fullfile(cfg.outputDir, ...
    ['wss_ensemble_model_values_', cfg.version, '.csv']));

LiteratureTable = buildLiteratureTable(literature);
writetable(LiteratureTable, fullfile(cfg.outputDir, ...
    ['wss_ensemble_literature_values_', cfg.version, '.csv']));

overlay = struct();
overlay.resultIndices = eligible;
overlay.representativeLocalIndex = representativeLocal;
overlay.modelShearDyn = modelShearDyn;
overlay.modelWss = modelWss;
overlay.foldValues = foldValues;
overlay.literature = literature;
overlay.figurePng = [figureBase, '.png'];
overlay.figureFig = [figureBase, '.fig'];
end

function plotOverlayPanel(panelIndex, shear, curves, representativeLocal, ...
        literatureShear, literatureFold, yLabelText, literatureLabel, cfg)
subplot(2, 2, panelIndex);
hold on;
if isvector(curves)
    curves = curves(:);
end
hEnsemble = plot(shear, curves(:, 1), '-', 'Color', [0.65 0.75 0.85], ...
    'LineWidth', 0.8);
for i = 2:size(curves, 2)
    plot(shear, curves(:, i), '-', 'Color', [0.65 0.75 0.85], ...
        'LineWidth', 0.8);
end
hRepresentative = plot(shear, curves(:, representativeLocal), '-', ...
    'Color', [0 0.30 0.55], 'LineWidth', 2.4);
hLiterature = plot(literatureShear, literatureFold, '--o', ...
    'Color', [0.15 0.15 0.15], 'LineWidth', 1.8, 'MarkerSize', 5);
xlabel('Shear stress (dyn/cm^2)');
ylabel(yLabelText);
xlim([0 cfg.shearMaxDyn]);
allValues = [curves(:); literatureFold(:)];
upper = max(allValues(isfinite(allValues)));
if isempty(upper)
    upper = 1;
end
ylim([0, max(1.1, 1.08 * upper)]);
grid on;
box off;
legend([hEnsemble, hRepresentative, hLiterature], ...
    {'Max-accuracy finalist models','Representative hybrid',literatureLabel}, ...
    'Location', 'best');
end

function literature = literatureAnchors()
literature.uematsuShear = [0, 1.2, 3, 15];
literature.uematsuENOS = [1.00, 1.35, 1.75, 2.72];
literature.uematsuNOShear = [0, 15];
literature.uematsuNO = [1.00, 2.00];
literature.malekShear = [0, 3, 7, 15, 20, 36];
literature.malekET1 = [1.00, 0.95, 0.62, 0.35, 0.31, 0.20];
end

function T = buildLiteratureTable(literature)
ShearDynCm2 = [literature.uematsuShear(:); literature.uematsuNOShear(:); ...
    literature.malekShear(:)];
LiteratureFold = [literature.uematsuENOS(:); literature.uematsuNO(:); ...
    literature.malekET1(:)];
ReferenceCurve = [ ...
    repmat({'Uematsu 1995 ecNOS mRNA'}, numel(literature.uematsuShear), 1); ...
    repmat({'Uematsu 1995 NO/nitrite'}, numel(literature.uematsuNOShear), 1); ...
    repmat({'Malek 1995 ET-1 mRNA'}, numel(literature.malekShear), 1)];
T = table(ShearDynCm2, LiteratureFold, ReferenceCurve);
end

function fold = safeFold(values, baseline)
if abs(baseline) < eps
    fold = values ./ eps;
else
    fold = values ./ baseline;
end
end

function idx = findSpeciesIndices(speciesNames, requestedSpecies, warnIfMissing)
idx = [];
for i = 1:numel(requestedSpecies)
    thisIdx = find(strcmp(speciesNames, requestedSpecies{i}), 1);
    if isempty(thisIdx)
        if warnIfMissing
            warning('Requested species not found: %s', requestedSpecies{i});
        end
    else
        idx(end + 1) = thisIdx; %#ok<AGROW>
    end
end
end

function atomicSaveMat(targetFile, payload)
targetDir = fileparts(targetFile);
ensureDir(targetDir);

stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
stagingFile = [tempname(tempdir), '.mat'];
incomingFile = [targetFile, '.incoming_', stamp, '.mat'];
previousFile = [targetFile, '.previous_', stamp];
expectedVariables = fieldnames(payload);
previousMoved = false;

try
    save(stagingFile, '-struct', 'payload', '-v7');
    verifyMatVariables(stagingFile, expectedVariables);

    [copied, copyMessage] = copyfile(stagingFile, incomingFile, 'f');
    if ~copied
        error('EC:AtomicCopyFailed', ...
            'Could not copy the staged MAT file to OneDrive: %s', copyMessage);
    end
    verifyMatVariables(incomingFile, expectedVariables);

    if exist(targetFile, 'file') == 2
        [previousMoved, moveMessage] = movefile(targetFile, previousFile, 'f');
        if ~previousMoved
            error('EC:AtomicBackupFailed', ...
                'Could not move the existing target aside: %s', moveMessage);
        end
    end

    [installed, installMessage] = movefile(incomingFile, targetFile, 'f');
    if ~installed
        error('EC:AtomicInstallFailed', ...
            'Could not install the verified MAT file: %s', installMessage);
    end
    verifyMatVariables(targetFile, expectedVariables);

    safeDeleteFile(stagingFile);
    if previousMoved
        safeDeleteFile(previousFile);
    end
catch ME
    if previousMoved && exist(previousFile, 'file') == 2
        if exist(targetFile, 'file') == 2
            failedTarget = [targetFile, '.failed_', stamp];
            movefile(targetFile, failedTarget, 'f');
        end
        movefile(previousFile, targetFile, 'f');
    end

    recoveryNames = fieldnames(payload);
    for i = 1:numel(recoveryNames)
        assignin('base', ['EC_recovery_', recoveryNames{i}], ...
            payload.(recoveryNames{i}));
    end
    assignin('base', 'EC_recovery_stagingFile', stagingFile);
    assignin('base', 'EC_recovery_incomingFile', incomingFile);
    assignin('base', 'EC_recovery_targetFile', targetFile);

    error('EC:AtomicSaveFailed', ...
        ['Could not safely install %s. The completed data remain in %s ', ...
        'and in base-workspace EC_recovery_* variables. Cause: %s'], ...
        targetFile, stagingFile, ME.message);
end
end

function verifyMatVariables(filename, expectedVariables)
info = whos('-file', filename);
savedVariables = {info.name};
missing = setdiff(expectedVariables, savedVariables);
if ~isempty(missing)
    error('EC:MatVerificationFailed', ...
        'MAT file %s is missing variables: %s', filename, strjoin(missing, ', '));
end
end

function quarantineFile = quarantineUnreadableMat(filename)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
quarantineFile = [filename, '.unreadable_', stamp];
[moved, moveMessage] = movefile(filename, quarantineFile, 'f');
if ~moved
    warning('EC:CheckpointQuarantineFailed', ...
        'Could not quarantine unreadable checkpoint %s: %s', filename, moveMessage);
    quarantineFile = filename;
end
end

function safeDeleteFile(filename)
if exist(filename, 'file') ~= 2
    return;
end
try
    delete(filename);
catch ME
    warning('EC:TemporaryFileCleanupFailed', ...
        'Could not remove temporary file %s: %s', filename, ME.message);
end
end
