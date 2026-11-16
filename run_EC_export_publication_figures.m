%% Reproduce manuscript Figures 2-5 from the archived data

%% USER OPTIONS
cfg = struct();
cfg.outputDir = fullfile('results','figures');
cfg.constantResultsDir = fullfile('data','constant');
cfg.perturbedResultsDir = fullfile('data','reaction_specific');
cfg.hybridResultsDir = fullfile('data','hybrid');
cfg.figures = 2:5;

cfg.fontName = 'Arial';
cfg.fontSize = 11;
cfg.matrixTickFontSize = 10;
cfg.titleFontSize = 12;
cfg.lineWidth = 1.5;
cfg.matrixSizeInches = [7.4 2.55];
cfg.histogramSizeInches = [4.4 3.4];
cfg.curveSizeInches = [4.4 2.9];
cfg.resolution = 600;
cfg.formats = {'pdf','png'}; % Also supported: svg, tif, fig.
cfg.exportCombined = true;
cfg.showTitles = true;
cfg.closeFiguresAfterSaving = true;

cfg.exportValueLabeled = false;
cfg.valueOutputDir = fullfile(cfg.outputDir,'Value_Labeled');
cfg.valueMatrixSizeInches = [18.6 4.0];
cfg.valueFontSize = 8;
cfg.valueFormat = '%.4g';

cfg.figure2ReconstructionVersion = 'v6';
cfg.figure2ReconstructionParameters = [0.85 1.4 0.6]; % w, n, EC50
cfg.figure2ReconstructionInputValue = 0.2;
cfg.figure2ReconstructionStressValue = 0.5;
cfg.figure2ReconstructionPerturbation = 0.5;
cfg.figure2ReconstructionThreshold = 0.0001;

cfg.figure2NumericCsv = fullfile('data','figures','Figure_2_source.csv');
cfg.figure2UseReconstructionForChecks = false; % Set true to rerun the Figure 2 ODEs.

cfg.histogramBinEdges = 0:1:100;
cfg.histogramYLim = [0 5000];

cfg.wssVersion = 'AngII_v5';
cfg.wssParameters = [0.772 1.599 0.5702]; % w, n, EC50
cfg.wssValues = linspace(0,1,121);
cfg.wssRelTol = 1e-8;
cfg.wssAbsTol = 1e-12;
cfg.wssInputValue = 0.2;
cfg.wssStressValue = 0.5;
cfg.wssInputOverrides = {'BMP9_10',1};
cfg.reuseWssCache = true;
% Set empty to recompute all 121 WSS simulations using the settings above.
cfg.wssSourceCsv = fullfile('data','figures','Figure_5_wss_values.csv');

%% EXPORT
publicationResults = EC_export_publication_figures(cfg);
fprintf('\nPublication figures: %s\n',publicationResults.outputDir);
