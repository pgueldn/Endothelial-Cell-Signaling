%% Constant-spine / perturbed-donor ensemble optimization

clearvars;
clc;

%% USER OPTIONS: folders and model

% Run the constant and reaction-specific searches first.
cfg = struct();
cfg.constantResultsDir = fullfile('results','constant');
cfg.perturbedResultsDir = fullfile('results','reaction_specific');
cfg.version = 'AngII_v5';
cfg.outputDir = fullfile('results','hybrid');

cfg.constantScoreFile = 'sampled_global_constant_parameters.csv';

cfg.perturbedScoreFile = 'top_validation_scores.csv';
cfg.perturbedParameterFile = 'top_parameter_sets_long.csv';

cfg.rerankSourceModelsForCurrentValidation = true;
cfg.constantAllResultsMatFile = ...
    'MonteCarlo_constant_parameter_validation_fit_all_results.mat';
cfg.perturbedAllResultsMatFile = ...
    'MonteCarlo_validation_fit_all_results.mat';
cfg.sourceScoreChunkSize = 2000;

%% USER OPTIONS: model ensemble selection

cfg.constantScoreToleranceMatches = 0;
cfg.keepAllMaxScoreConstant = true;
cfg.maxConstantSpines = Inf;

cfg.perturbedScoreToleranceMatches = 1;
cfg.keepAllMaxScorePerturbed = true;

cfg.maxPerturbedDonors = 10;
cfg.diversityUseModeledEdgesOnly = true;

%% USER OPTIONS: pair screening

cfg.runPairScreening = true;
cfg.reuseExistingScreening = false;
cfg.screenMaxPrunePasses = 2;
cfg.screenCheckpointEvery = 10;
cfg.useParallelPairScreening = true;

cfg.numFinalistPairs = 5;
cfg.diversifyFinalistSources = true;

%% USER OPTIONS: robust finalist genetic algorithm

cfg.populationSize = 250;
cfg.numGenerations = 200;
cfg.rngSeeds = 1:5;
cfg.eliteCount = 5;
cfg.tournamentSize = 3;
cfg.mutationRate = 0.03;
cfg.complexityPenalty = 1e-3;
cfg.finalMaxPrunePasses = 2;
cfg.runGreedyForward = false;
cfg.runGreedyPrune = true;
cfg.runGreedyForwardAfterPrune = false;

cfg.useParallelPopulationEvaluation = true;
cfg.parallelWorkers = [];

cfg.resumeCompletedFinalRuns = true;

%% USER OPTIONS: validation and model inputs

cfg.validationDirectionThreshold = 0.02;
cfg.candidateMode = 'used'; % 'used' or 'all'

cfg.activeInputs = 1:5;
cfg.inputValue = 0.2;
cfg.stressInputValue = 0.5;
cfg.perturbValue = 0.5;
cfg.validationPerturbSpecies = {'Stress','Wss','Hypoxia','TNFa'};
cfg.validationSpeciesSubset = {'VEGF','VEGFR','IL6','ppET1','ET1','eNOS', ...
    'FAK','TGFb','pSmad23','pSmad15','p38','ERK','JNK','PKC','DAG', ...
    'PI3K','Akt','mTOR','p70S6K','RhoA','ROCK','MLCK','Actin','NFkB', ...
    'HIF1a','MCP1','TSP1','Proliferation','Migration','Apoptosis','EndMT'};

%% USER OPTIONS: outputs and WSS overlays

cfg.saveIndividualFinalRunFiles = true;
cfg.saveIndividualValidationPlots = false;
cfg.makeRepresentativeValidationPlot = true;

cfg.makeWssOverlay = false; % Optional approximate literature guides, not Figure 5.
cfg.shearMaxDyn = 36;
cfg.numShearPoints = 121;
cfg.maxOverlayModels = 25;
cfg.closeFiguresAfterSaving = true;

%% USER OPTIONS: execution controls

cfg.stopOnODEError = false;
cfg.verbose = true;

cfg.quickTestMode = false;
% Full run: 5 finalist pairs x 5 seeds x 250 individuals x 200 generations.

%% RUN

ensembleResults = EC_run_constantspineperturbedmodel_ensemble(cfg);
