function results = EC_export_publication_figures(cfg)
% Plot saved manuscript results; ODE recomputation is opt-in.
root = fileparts(mfilename('fullpath'));
cfg.outputDir = absolutePath(root,cfg.outputDir);
cfg.constantResultsDir = absolutePath(root,cfg.constantResultsDir);
cfg.perturbedResultsDir = absolutePath(root,cfg.perturbedResultsDir);
cfg.hybridResultsDir = absolutePath(root,cfg.hybridResultsDir);
if cfg.exportValueLabeled
    cfg.valueOutputDir = absolutePath(root,cfg.valueOutputDir);
    if ~isfolder(cfg.valueOutputDir), mkdir(cfg.valueOutputDir); end
    valueDataDir = fullfile(cfg.valueOutputDir,'source_data');
    if ~isfolder(valueDataDir), mkdir(valueDataDir); end
end
if ~isempty(cfg.figure2NumericCsv)
    cfg.figure2NumericCsv = absolutePath(root,cfg.figure2NumericCsv);
end
assert(any(strcmpi(listfonts,cfg.fontName)),...
    'EC:FontMissing','Requested font is not installed: %s',cfg.fontName);
assert(cfg.fontSize>=9 && cfg.matrixTickFontSize>=9,...
    'Use at least 9 point text at the exported physical size.');
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
dataDir = fullfile(cfg.outputDir,'source_data');
if ~isfolder(dataDir), mkdir(dataDir); end

species = {'VEGF','VEGFR','IL6','ppET1','ET1','eNOS','FAK','TGFb',...
    'pSmad23','pSmad15','p38','ERK','JNK','PKC','DAG','PI3K','Akt',...
    'mTOR','p70S6K','RhoA','ROCK','MLCK','Actin','NFkB','HIF1a',...
    'MCP1','TSP1','Proliferation','Migration','Apoptosis','EndMT'};
perturbations = {'Stress','Wss','Hypoxia','TNFa'};
audit = cell(0,6);
valueAudit = cell(0,6);
sources = struct();
sources.generatedAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
sources.figure2 = 'Not requested';
sources.figure3 = 'Not requested';
sources.figure4 = 'Not requested';
sources.figure5 = 'Not requested';

if ismember(2,cfg.figures)
    if cfg.figure2UseReconstructionForChecks
        [panel,reconstructedTable] = reconstructFigure2(root,species,perturbations,cfg);
        source = fullfile(dataDir,'Figure_2_reconstructed_values.csv');
        writetable(reconstructedTable,source);
        sources.figure2 = struct('file',source,'sha256',sha256(source),...
            'type','Numerical reconstruction from current model files',...
            'matches',panel.matches,'known',panel.known);
        figure2Label = sprintf('Investigator-selected parameters: current %s reconstruction',...
            cfg.figure2ReconstructionVersion);
    else
        source = cfg.figure2NumericCsv;
        panel = loadValidation(source,species,perturbations);
        sources.figure2 = struct('file',source,'sha256',sha256(source),...
            'type','Saved numerical validation table',...
            'matches',panel.matches,'known',panel.known);
        figure2Label = sprintf('Investigator-selected parameters: current %s reconstruction',...
            cfg.figure2ReconstructionVersion);
    end
    copyfile(source,fullfile(dataDir,'Figure_2_source.csv'));
    fig = newFigure('Figure 2',cfg.matrixSizeInches);
    drawValidation(fig,[0 0 1 1],panel,species,perturbations,...
        figure2Label,'',cfg);
    audit(end+1,:) = exportFigure(fig,'Figure_2_validation',cfg);
    if cfg.exportValueLabeled
        if isempty(panel.diff)
            [valuePanel,valueTable] = reconstructFigure2(root,species,perturbations,cfg);
            valueSource = fullfile(valueDataDir,'Figure_2_reconstructed_values.csv');
            writetable(valueTable,valueSource);
            valueLabel = sprintf('Investigator-selected parameters: current %s reconstruction',...
                cfg.figure2ReconstructionVersion);
        else
            valuePanel = panel;
            if cfg.figure2UseReconstructionForChecks
                valueSource = fullfile(valueDataDir,'Figure_2_reconstructed_values.csv');
            else
                valueSource = fullfile(valueDataDir,'Figure_2_values.csv');
            end
            copyfile(source,valueSource);
            valueLabel = 'Investigator-selected parameters';
        end
        fig = newFigure('Figure 2 values',cfg.valueMatrixSizeInches);
        drawValidationValues(fig,[0 0 1 1],valuePanel,species,perturbations,...
            valueLabel,'',cfg);
        valueCfg = cfg;
        valueCfg.outputDir = cfg.valueOutputDir;
        valueCfg.minimumFontSize = cfg.valueFontSize;
        valueAudit(end+1,:) = exportFigure(fig,'Figure_2_values',valueCfg);
    end
end

if ismember(3,cfg.figures)
    dirs = {cfg.constantResultsDir,cfg.perturbedResultsDir};
    labels = {'Constant parameters','Reaction-specific parameters'};
    panels = cell(1,2);
    stats = cell(2,6);
    sourceInfo = cell(1,2);
    for k = 1:2
        source = fullfile(dirs{k},'validation_scores_by_run.csv');
        T = readtable(source);
        valid = ~logical(T.Failed) & isfinite(T.ValidationAgreementPct);
        vals = T.ValidationAgreementPct(valid);
        assert(~isempty(vals),'No valid scores in %s',source);
        assert(all(vals>=0 & vals<=100),'Validation scores must be percentages.');
        panels{k} = vals;
        fig = newFigure(['Figure 3',char('A'+k-1)],cfg.histogramSizeInches);
        drawHistogram(fig,[0 0 1 1],vals,labels{k},'',cfg);
        audit(end+1,:) = exportFigure(fig,['Figure_3',char('A'+k-1),'_score_distribution'],cfg); %#ok<AGROW>
        stats(k,:) = {labels{k},numel(vals),nnz(~valid),mean(vals),std(vals),max(vals)};
        copyfile(source,fullfile(dataDir,['Figure_3',char('A'+k-1),'_scores.csv']));
        counts = histcounts(vals,cfg.histogramBinEdges);
        bins = table(cfg.histogramBinEdges(1:end-1)',cfg.histogramBinEdges(2:end)',counts',...
            'VariableNames',{'LeftEdgePct','RightEdgePct','Count'});
        assert(sum(counts)==numel(vals),'Histogram omitted valid scores.');
        writetable(bins,fullfile(dataDir,['Figure_3',char('A'+k-1),'_histogram_bins.csv']));
        sourceInfo{k} = struct('file',source,'sha256',sha256(source),'validRuns',numel(vals));
    end
    writetable(cell2table(stats,'VariableNames',...
        {'Model','ValidRuns','ExcludedRuns','MeanPct','StdPct','MaximumPct'}),...
        fullfile(dataDir,'Figure_3_summary.csv'));
    sources.figure3 = sourceInfo;
    if cfg.exportCombined
        fig = newFigure('Figure 3',[7.4 3.45]);
        for k = 1:2
            drawHistogram(fig,[(k-1)*0.5 0 0.5 1],panels{k},labels{k},char('A'+k-1),cfg);
        end
        audit(end+1,:) = exportFigure(fig,'Figure_3_combined',cfg);
    end
end

if ismember(4,cfg.figures)
    dirs = {cfg.constantResultsDir,cfg.perturbedResultsDir,...
        fullfile(cfg.hybridResultsDir,'representative_best_model')};
    labels = {'Constant parameters','Reaction-specific parameters',...
        'Constant spine with substitutions'};
    panels = cell(1,3);
    sourceInfo = cell(1,3);
    for k = 1:3
        source = fullfile(dirs{k},'best_validation_direction_matrix_long.csv');
        panels{k} = loadValidation(source,species,perturbations);
        fig = newFigure(['Figure 4',char('A'+k-1)],cfg.matrixSizeInches);
        drawValidation(fig,[0 0 1 1],panels{k},species,perturbations,labels{k},'',cfg);
        audit(end+1,:) = exportFigure(fig,['Figure_4',char('A'+k-1),'_validation'],cfg); %#ok<AGROW>
        copyfile(source,fullfile(dataDir,['Figure_4',char('A'+k-1),'_validation.csv']));
        sourceInfo{k} = struct('file',source,'sha256',sha256(source),...
            'matches',panels{k}.matches,'known',panels{k}.known,...
            'classification','Uses saved ObservedDirection and Matched, without rescoring');
        if cfg.exportValueLabeled
            fig = newFigure(['Figure 4',char('A'+k-1),' values'],cfg.valueMatrixSizeInches);
            drawValidationValues(fig,[0 0 1 1],panels{k},species,perturbations,...
                labels{k},'',cfg);
            valueCfg = cfg;
            valueCfg.outputDir = cfg.valueOutputDir;
            valueCfg.minimumFontSize = cfg.valueFontSize;
            valueAudit(end+1,:) = exportFigure(fig,...
                ['Figure_4',char('A'+k-1),'_values'],valueCfg); %#ok<AGROW>
            copyfile(source,fullfile(valueDataDir,...
                ['Figure_4',char('A'+k-1),'_values.csv']));
        end
    end
    sources.figure4 = sourceInfo;
    if cfg.exportCombined
        fig = newFigure('Figure 4',[cfg.matrixSizeInches(1) 3*cfg.matrixSizeInches(2)]);
        for k = 1:3
            drawValidation(fig,[0 (3-k)/3 1 1/3],panels{k},species,perturbations,...
                labels{k},char('A'+k-1),cfg);
        end
        audit(end+1,:) = exportFigure(fig,'Figure_4_combined',cfg);
    end
    if cfg.exportValueLabeled && cfg.exportCombined
        fig = newFigure('Figure 4 values',[cfg.valueMatrixSizeInches(1) 3*cfg.valueMatrixSizeInches(2)]);
        for k = 1:3
            drawValidationValues(fig,[0 (3-k)/3 1 1/3],panels{k},species,perturbations,...
                labels{k},char('A'+k-1),cfg);
        end
        valueCfg = cfg;
        valueCfg.outputDir = cfg.valueOutputDir;
        valueCfg.minimumFontSize = cfg.valueFontSize;
        valueAudit(end+1,:) = exportFigure(fig,'Figure_4_values_combined',valueCfg);
    end
end

if ismember(5,cfg.figures)
    if isfield(cfg,'wssSourceCsv') && ~isempty(cfg.wssSourceCsv)
        source = absolutePath(root,cfg.wssSourceCsv);
        W = readtable(source);
        signature = struct('file',source,'sha256',sha256(source),...
            'source','Archived Figure 5 numerical data; no ODEs recomputed');
        copyfile(source,fullfile(dataDir,'Figure_5_wss_values.csv'));
    else
        [W,signature] = simulateWss(root,dataDir,cfg);
    end
    names = {'NO','ET1','pSmad23','pSmad15','NFkB'};
    labels = {'Nitric oxide','Endothelin-1','pSmad2/3','pSmad1/5',['NF-',char(954),'B']};
    yLimits = {[0 10],[0 200],[0 1.25],[0.9 1.8],[0.5 2.5]};
    for k = 1:5
        fig = newFigure(['Figure 5',char('A'+k-1)],cfg.curveSizeInches);
        drawCurve(fig,[0 0 1 1],W.Wss,W.([names{k},'_Fold']),labels{k},yLimits{k},'',cfg,false);
        audit(end+1,:) = exportFigure(fig,['Figure_5',char('A'+k-1),'_',names{k}],cfg); %#ok<AGROW>
    end
    sources.figure5 = signature;
    if cfg.exportCombined
        fig = newFigure('Figure 5',[7.4 9.0]);
        for k = 1:5
            drawCurve(fig,[0 (5-k)/5 1 1/5],W.Wss,W.([names{k},'_Fold']),...
                labels{k},yLimits{k},char('A'+k-1),cfg,true);
        end
        audit(end+1,:) = exportFigure(fig,'Figure_5_combined',cfg);
    end
end

auditTable = cell2table(audit,'VariableNames',...
    {'Figure','WidthInches','HeightInches','MinimumFontPt','FontNames','AxesCount'});
writetable(auditTable,fullfile(cfg.outputDir,'export_audit.csv'));
if cfg.exportValueLabeled
    valueAuditTable = cell2table(valueAudit,'VariableNames',...
        {'Figure','WidthInches','HeightInches','MinimumFontPt','FontNames','AxesCount'});
    writetable(valueAuditTable,fullfile(cfg.valueOutputDir,'export_audit.csv'));
else
    valueAuditTable = table();
end
writeText(fullfile(dataDir,'provenance.json'),jsonencode(sources,PrettyPrint=true));
save(fullfile(dataDir,'export_configuration.mat'),'cfg','sources');
writeReadme(cfg,sources);
if cfg.exportValueLabeled, writeValueReadme(cfg); end
results = struct('outputDir',cfg.outputDir,'audit',auditTable,'sources',sources,...
    'valueOutputDir',cfg.valueOutputDir,'valueAudit',valueAuditTable);
end

function path = absolutePath(root,path)
if isempty(regexp(path,'^[A-Za-z]:[\\/]|^[/\\]','once'))
    path = fullfile(root,path);
end
end


function panel = loadValidation(file,species,perturbations)
T = readtable(file,'TextType','string');
required = {'Perturbation','Species','ExpectedDirection','ObservedDifference',...
    'ObservedDirection','KnownLiteratureCheck','Matched'};
assert(all(ismember(required,T.Properties.VariableNames)),'Invalid validation table: %s',file);
assert(all(isfinite(T.ObservedDifference)),'Nonfinite validation results: %s',file);
known = logical(T.KnownLiteratureCheck);
assert(isequal(known,abs(T.ExpectedDirection)~=0.5),'Known-check mask is inconsistent.');
assert(isequal(logical(T.Matched),known & T.ObservedDirection==T.ExpectedDirection),...
    'Stored validation matches are inconsistent in %s',file);
panel.matches = nnz(T.Matched);
panel.known = nnz(known);
panel.diff = nan(numel(perturbations),numel(species));
panel.markers = zeros(size(panel.diff));
for p = 1:numel(perturbations)
    for j = 1:numel(species)
        row = T.Perturbation==perturbations{p} & T.Species==species{j};
        assert(nnz(row)==1,'Missing or duplicate validation cell.');
        panel.diff(p,j) = T.ObservedDifference(row);
        if known(row)
            panel.markers(p,j) = 2*double(T.Matched(row))-1;
        end
    end
end
assert(nnz(panel.markers)==panel.known,'Species subset omits defined checks.');
panel.rgb = [];
end

function [panel,T] = reconstructFigure2(root,species,perturbations,cfg)
version = cfg.figure2ReconstructionVersion;
refFile = fullfile(root,['refcase_EC_Network_',version,'.mat']);
odeFile = fullfile(root,['ODElist_EC_Network_',version,'.mat']);
M = load(refFile,'speciesNames','tau','ymax','numEdges');
[found,speciesIdx] = ismember(species,M.speciesNames);
assert(all(found),'The Figure 2 reconstruction is missing a displayed species.');
[found,inputIdx] = ismember({'Stress','Wss','Hypoxia','Cav1','TNFa'},M.speciesNames);
assert(all(found),'The Figure 2 reconstruction is missing an input species.');
rpar = repmat(cfg.figure2ReconstructionParameters(:),1,M.numEdges);
rpar(1,6) = 0;
y0 = zeros(1,numel(M.speciesNames));
y0(inputIdx) = cfg.figure2ReconstructionInputValue;
y0(inputIdx(1)) = cfg.figure2ReconstructionStressValue;
[~,baseline] = ODE_master(rpar,y0,M.speciesNames,M.tau,M.ymax,M.numEdges,odeFile);
allDiff = zeros(4,numel(M.speciesNames));
perturbIdx = inputIdx([1 2 3 5]);
for p = 1:4
    perturbed = baseline(end,:);
    perturbed(perturbIdx(p)) = perturbed(perturbIdx(p))+...
        cfg.figure2ReconstructionPerturbation;
    [~,response] = ODE_master(rpar,perturbed,M.speciesNames,M.tau,M.ymax,M.numEdges,odeFile);
    allDiff(p,:) = response(end,:)-baseline(end,:);
end
expected = CreateValidationMatrix(M.speciesNames,4);
observed = allDiff;
observed(abs(observed)<cfg.figure2ReconstructionThreshold) = 0;
observed = sign(observed);
known = abs(expected)~=0.5;
matched = known & observed==expected;
panel.diff = allDiff(:,speciesIdx);
panel.rgb = [];
panel.matches = nnz(matched);
panel.known = nnz(known);
panel.markers = zeros(size(panel.diff));
subsetKnown = known(:,speciesIdx);
subsetMatched = matched(:,speciesIdx);
panel.markers(subsetKnown) = 2*double(subsetMatched(subsetKnown))-1;

rowCount = numel(perturbations)*numel(species);
Perturbation = strings(rowCount,1);
Species = strings(rowCount,1);
ExpectedDirection = zeros(rowCount,1);
ObservedDifference = zeros(rowCount,1);
ObservedDirection = zeros(rowCount,1);
KnownLiteratureCheck = false(rowCount,1);
Matched = false(rowCount,1);
row = 0;
for p = 1:numel(perturbations)
    for j = 1:numel(species)
        row = row+1;
        idx = speciesIdx(j);
        Perturbation(row) = perturbations{p};
        Species(row) = species{j};
        ExpectedDirection(row) = expected(p,idx);
        ObservedDifference(row) = allDiff(p,idx);
        ObservedDirection(row) = observed(p,idx);
        KnownLiteratureCheck(row) = known(p,idx);
        Matched(row) = matched(p,idx);
    end
end
T = table(Perturbation,Species,ExpectedDirection,ObservedDifference,...
    ObservedDirection,KnownLiteratureCheck,Matched);
end

function fig = newFigure(name,sizeInches)
fig = figure('Name',name,'Color','w','Visible','off',...
    'Units','inches','Position',[0.5 0.5 sizeInches],...
    'PaperUnits','inches','PaperSize',sizeInches,...
    'PaperPosition',[0 0 sizeInches],'InvertHardcopy','off');
end

function rect = within(region,local)
rect = [region(1)+local(1)*region(3),region(2)+local(2)*region(4),...
    local(3)*region(3),local(4)*region(4)];
end

function styleAxes(ax,cfg,sizePt)
set(ax,'FontName',cfg.fontName,'FontSize',sizePt,...
    'TickLabelInterpreter','none','LineWidth',0.7,...
    'XColor','k','YColor','k','Layer','top','Box','off',...
    'TickDir','out','TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
ax.Toolbar.Visible = 'off';
end

function drawValidation(fig,region,panel,species,perturbations,label,letter,cfg)
ax = axes(fig,'Position',within(region,[0.095 0.34 0.81 0.46]));
styleAxes(ax,cfg,cfg.matrixTickFontSize);
set(ax,'YDir','reverse','TickLength',[0 0],...
    'XLim',[0.5 numel(species)+0.5],'YLim',[0.5 4.5],...
    'XTick',1:numel(species),'XTickLabel',displaySpecies(species),...
    'XTickLabelRotation',270,'YTick',1:4,...
    'YTickLabel',{'Stress','WSS','Hypoxia',['TNF',char(945)]});
clim(ax,[-0.5 1]);
axes(ax);
map = custom_cmap(256);
colormap(ax,map);
hold(ax,'on');
for p = 1:numel(perturbations)
    for j = 1:numel(species)
        if isempty(panel.rgb)
            colorIndex = min(256,max(1,1+floor((panel.diff(p,j)+0.5)/1.5*256)));
            color = map(colorIndex,:);
        else
            color = reshape(panel.rgb(p,j,:),1,3);
        end
        patch(ax,j+[-.5 .5 .5 -.5],p+[-.5 -.5 .5 .5],color,...
            'EdgeColor','none','HandleVisibility','off');
    end
end
for j = 0.5:1:numel(species)+0.5
    line(ax,[j j],[0.5 4.5],'Color',[.2 .2 .2],'LineWidth',0.45);
end
for p = 0.5:1:4.5
    line(ax,[0.5 numel(species)+0.5],[p p],'Color',[.2 .2 .2],'LineWidth',0.45);
end
for p = 1:4
    for j = 1:numel(species)
        if panel.markers(p,j)==1
            line(ax,j+[-.16 -.035 .19],p+[-.005 .18 -.22],...
                'Color','k','LineWidth',1.0,'Clipping','on');
        elseif panel.markers(p,j)==-1
            line(ax,j+[-.14 .14 NaN -.14 .14],p+[-.18 .18 NaN .18 -.18],...
                'Color','k','LineWidth',1.0,'Clipping','on');
        end
    end
end
drawVectorColorbar(fig,within(region,[0.927 0.34 0.017 0.46]),map,cfg);
if cfg.showTitles
    labelText(fig,within(region,[0.095 0.82 0.77 0.14]),...
        label,cfg.titleFontSize,cfg,'left','bold');
end
labelText(fig,within(region,[0.70 0.82 0.25 0.14]),...
    sprintf('%d/%d (%.2f%%)',panel.matches,panel.known,100*panel.matches/panel.known),...
    cfg.matrixTickFontSize,cfg,'right','normal');
panelLetter(fig,region,letter,cfg);
end

function drawValidationValues(fig,region,panel,species,perturbations,label,letter,cfg)
assert(~isempty(panel.diff),'Numerical differences are required for value labels.');
ax = axes(fig,'Position',within(region,[0.06 0.31 0.86 0.51]));
styleAxes(ax,cfg,cfg.matrixTickFontSize);
set(ax,'YDir','reverse','TickLength',[0 0],...
    'XLim',[0.5 numel(species)+0.5],'YLim',[0.5 4.5],...
    'XTick',1:numel(species),'XTickLabel',displaySpecies(species),...
    'XTickLabelRotation',270,'YTick',1:4,...
    'YTickLabel',{'Stress','WSS','Hypoxia',['TNF',char(945)]});
clim(ax,[-0.5 1]);
axes(ax);
map = custom_cmap(256);
colormap(ax,map);
hold(ax,'on');
for p = 1:numel(perturbations)
    for j = 1:numel(species)
        value = panel.diff(p,j);
        colorIndex = min(256,max(1,1+floor((value+0.5)/1.5*256)));
        color = map(colorIndex,:);
        patch(ax,j+[-.5 .5 .5 -.5],p+[-.5 -.5 .5 .5],color,...
            'EdgeColor',[.25 .25 .25],'LineWidth',.45,'HandleVisibility','off');
        luminance = .2126*color(1)+.7152*color(2)+.0722*color(3);
        textColor = 'k';
        if luminance<0.42, textColor = 'w'; end
        text(ax,j,p,sprintf(cfg.valueFormat,value),...
            'FontName',cfg.fontName,'FontSize',cfg.valueFontSize,...
            'FontWeight','normal','HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',textColor,...
            'Interpreter','none','Clipping','on');
    end
end
drawVectorColorbar(fig,within(region,[0.94 0.31 0.012 0.51]),map,cfg);
if cfg.showTitles
    labelText(fig,within(region,[0.06 0.84 0.72 0.11]),...
        label,cfg.titleFontSize,cfg,'left','bold');
end
labelText(fig,within(region,[0.76 0.84 0.19 0.11]),...
    sprintf('%d/%d (%.2f%%)',panel.matches,panel.known,100*panel.matches/panel.known),...
    cfg.matrixTickFontSize,cfg,'right','normal');
panelLetter(fig,region,letter,cfg);
end

function drawVectorColorbar(fig,pos,map,cfg)
cb = axes(fig,'Position',pos);
styleAxes(cb,cfg,cfg.matrixTickFontSize);
set(cb,'YAxisLocation','right','XLim',[0 1],'XTick',[],...
    'YLim',[-.5 1],'YTick',[-.5 0 .5 1],'YTickLabel',{'-0.5','0','0.5','1'},...
    'TickLength',[0 0],'Box','on','LineWidth',.6);
hold(cb,'on');
edges = linspace(-.5,1,size(map,1)+1);
for k = 1:size(map,1)
    patch(cb,[0 1 1 0],[edges(k) edges(k) edges(k+1) edges(k+1)],...
        map(k,:),'EdgeColor','none','HandleVisibility','off');
end
end

function labels = displaySpecies(species)
labels = species;
for j = 1:numel(labels)
    switch labels{j}
        case 'pSmad23', labels{j} = 'pSmad2/3';
        case 'pSmad15', labels{j} = 'pSmad1/5';
        case 'TGFb', labels{j} = ['TGF',char(946)];
        case 'NFkB', labels{j} = ['NF-',char(954),'B'];
        case 'HIF1a', labels{j} = ['HIF1',char(945)];
    end
end
end

function drawHistogram(fig,region,values,label,letter,cfg)
ax = axes(fig,'Position',within(region,[0.18 0.2 0.77 0.65]));
styleAxes(ax,cfg,cfg.fontSize);
histogram(ax,values,cfg.histogramBinEdges,'FaceColor',[0.39 0.65 0.81],...
    'FaceAlpha',1,'EdgeColor',[0.2 0.2 0.2],'LineWidth',.4);
styleAxes(ax,cfg,cfg.fontSize);
xlim(ax,[0 100]);
ylim(ax,cfg.histogramYLim);
xticks(ax,0:20:100);
yticks(ax,0:1000:cfg.histogramYLim(2));
ax.YAxis.Exponent = 0;
xlabel(ax,'Validation agreement (%)','FontName',cfg.fontName,'FontSize',cfg.fontSize,'Interpreter','none');
ylabel(ax,'Number of models','FontName',cfg.fontName,'FontSize',cfg.fontSize,'Interpreter','none');
if cfg.showTitles
    title(ax,label,'FontName',cfg.fontName,'FontSize',cfg.titleFontSize,'Interpreter','none','FontWeight','normal');
end
panelLetter(fig,region,letter,cfg);
end

function drawCurve(fig,region,x,y,label,limits,letter,cfg,stacked)
if stacked
    pos = [0.16 0.34 0.80 0.48];
else
    pos = [0.18 0.23 0.77 0.57];
end
ax = axes(fig,'Position',within(region,pos));
styleAxes(ax,cfg,cfg.fontSize);
plot(ax,x,y,'Color',[0 0.447 0.741],'LineWidth',cfg.lineWidth);
styleAxes(ax,cfg,cfg.fontSize);
xlim(ax,[0 1]);
xticks(ax,0:.2:1);
if min(y)<limits(1) || max(y)>limits(2)
    span = max(y)-min(y);
    limits = [min(limits(1),min(y)-.04*span),max(limits(2),max(y)+.04*span)];
end
ylim(ax,limits);
xlabel(ax,'Normalized WSS','FontName',cfg.fontName,'FontSize',cfg.fontSize,'Interpreter','none');
ylabel(ax,'Fold vs. zero WSS','FontName',cfg.fontName,'FontSize',cfg.fontSize,'Interpreter','none');
if cfg.showTitles
    title(ax,label,'FontName',cfg.fontName,'FontSize',cfg.titleFontSize,'FontWeight','normal','Interpreter','none');
end
panelLetter(fig,region,letter,cfg);
end

function labelText(fig,pos,str,sizePt,cfg,alignment,weight)
annotation(fig,'textbox',pos,'String',str,'EdgeColor','none',...
    'FontName',cfg.fontName,'FontSize',sizePt,'Interpreter','none',...
    'HorizontalAlignment',alignment,'VerticalAlignment','middle',...
    'FontWeight',weight,'Margin',0,'FitBoxToText','off');
end

function panelLetter(fig,region,letter,cfg)
if ~isempty(letter)
    labelText(fig,within(region,[0.015 .88 .07 .10]),letter,13,cfg,'left','bold');
end
end

function [T,signature] = simulateWss(root,dataDir,cfg)
refFile = fullfile(root,['refcase_EC_Network_',cfg.wssVersion,'.mat']);
odeFile = fullfile(root,['ODElist_EC_Network_',cfg.wssVersion,'.mat']);
signature = struct('version',cfg.wssVersion,'parameters',cfg.wssParameters,...
    'wss',cfg.wssValues,'inputValue',cfg.wssInputValue,...
    'stressValue',cfg.wssStressValue,'inputOverrides',{cfg.wssInputOverrides},...
    'relativeTolerance',cfg.wssRelTol,'absoluteTolerance',cfg.wssAbsTol,...
    'refSha256',sha256(refFile),'odeSha256',sha256(odeFile),...
    'solverSha256',sha256(fullfile(root,'ODE_master.m')),...
    'source','SimulateODEs_Wss_Increasing.m settings, including BMP9_10 = 1',...
    'normalization','Each species divided by its final value at WSS=0');
cache = fullfile(dataDir,'Figure_5_wss_cache.mat');
if cfg.reuseWssCache && isfile(cache)
    C = load(cache,'T','savedSignature');
    if isequaln(C.savedSignature,signature)
        T = C.T;
        return
    end
end
M = load(refFile,'speciesNames','tau','ymax','numEdges');
names = {'NO','ET1','pSmad23','pSmad15','NFkB'};
[found,idx] = ismember(names,M.speciesNames);
assert(all(found),'A Figure 5 output species is absent.');
[found,inputIdx] = ismember({'Stress','Wss','Hypoxia','Cav1','TNFa'},M.speciesNames);
assert(all(found),'A required input species is absent.');
rpar = repmat(cfg.wssParameters(:),1,M.numEdges);
rpar(1,6) = 0;
y0 = zeros(1,numel(M.speciesNames));
y0(inputIdx) = cfg.wssInputValue;
y0(inputIdx(1)) = cfg.wssStressValue;
for k = 1:size(cfg.wssInputOverrides,1)
    j = find(strcmp(M.speciesNames,cfg.wssInputOverrides{k,1}));
    assert(isscalar(j),'Unknown or duplicate override species.');
    y0(j) = cfg.wssInputOverrides{k,2};
end
values = nan(numel(cfg.wssValues),numel(names));
odeOptions = odeset('RelTol',cfg.wssRelTol,'AbsTol',cfg.wssAbsTol);
for i = 1:numel(cfg.wssValues)
    yi = y0;
    yi(inputIdx(2)) = cfg.wssValues(i);
    [~,y] = ODE_master(rpar,yi,M.speciesNames,M.tau,M.ymax,M.numEdges,odeFile,odeOptions);
    assert(all(isfinite(y(:))),'Nonfinite WSS simulation.');
    values(i,:) = y(end,idx);
end
zero = find(cfg.wssValues==0,1);
assert(~isempty(zero) && all(values(zero,:)>0),'WSS fold normalization needs positive zero-WSS values.');
folds = values./values(zero,:);
T = table(cfg.wssValues(:),'VariableNames',{'Wss'});
for k = 1:numel(names)
    T.([names{k},'_Activity']) = values(:,k);
    T.([names{k},'_Fold']) = folds(:,k);
end
savedSignature = signature;
save(cache,'T','savedSignature');
writetable(T,fullfile(dataDir,'Figure_5_wss_values.csv'));
end

function row = exportFigure(fig,base,cfg)
drawnow;
objects = findall(fig,'-property','FontName');
set(objects,'FontName',cfg.fontName);
names = unique(string(get(objects,'FontName')));
assert(all(strcmpi(names,cfg.fontName)),'A figure contains a non-Arial font.');
fontObjects = findall(fig,'-property','FontSize');
sizes = cell2mat(get(fontObjects,'FontSize'));
minimumFontSize = 9;
if isfield(cfg,'minimumFontSize'), minimumFontSize = cfg.minimumFontSize; end
assert(all(sizes>=minimumFontSize),...
    'A figure contains text below the configured minimum font size.');
sizeInches = fig.Position(3:4);
basePath = fullfile(cfg.outputDir,base);
for k = 1:numel(cfg.formats)
    ext = cfg.formats{k};
    path = [basePath,'.',ext];
    stagedPath = [tempname,'.',ext];
    switch ext
        case 'fig'
            savefig(fig,stagedPath);
        case 'pdf'
            print(fig,stagedPath,'-dpdf','-vector');
        case 'svg'
            print(fig,stagedPath,'-dsvg','-vector');
        case {'png','tif'}
            driver = ['-d',ext];
            if strcmp(ext,'tif'), driver = '-dtiff'; end
            print(fig,stagedPath,driver,sprintf('-r%d',cfg.resolution));
        otherwise
            error('Unsupported figure format: %s',ext);
    end
    copyExport(stagedPath,path);
    delete(stagedPath);
    assert(isfile(path),'Missing export: %s',path);
end
row = {base,sizeInches(1),sizeInches(2),min(sizes),strjoin(names,', '),numel(findall(fig,'Type','axes'))};
fprintf('Exported %s (%g x %g in)\n',base,sizeInches);
if cfg.closeFiguresAfterSaving, close(fig); end
end

function copyExport(source,destination)
for attempt = 1:5
    try
        [ok,message] = copyfile(source,destination,'f');
    catch err
        ok = false;
        message = err.message;
    end
    if ok, return; end
    pause(0.5*attempt);
end
error('EC:PublicationWrite',...
    'Cannot write %s: %s\nThe completed export is retained at %s.',...
    destination,message,source);
end


function hash = sha256(file)
fid = fopen(file,'rb');
assert(fid~=-1,'Cannot open source: %s',file);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
hash = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end

function writeText(file,str)
fid = fopen(file,'w','n','UTF-8');
assert(fid~=-1,'Cannot write %s',file);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',str);
end

function writeReadme(cfg,~)
text = [
    "# Reproduced manuscript figures"
    ""
    "Generated by run_EC_export_publication_figures."
    "Saved data are plotted without rescoring; ODE recomputation is opt-in."
    "See the package README for model versions and historical scoring thresholds."
    "Figure 5 uses 121 normalized WSS values and folds relative to zero WSS."
    "source_data/provenance.json identifies the inputs to this export."
    "PDF/SVG exports are vector graphics. Raster resolution is set in cfg.resolution."
    ];
writeText(fullfile(cfg.outputDir,'README.md'),strjoin(text,newline));
end

function writeValueReadme(cfg)
text = [
    "# Numerical validation panels"
    ""
    "Labels show four significant digits of the saved absolute activity differences."
    "The source_data CSV files retain the saved numerical precision."
    "These panels use the same data and classifications as the check/cross panels."
    ];
writeText(fullfile(cfg.valueOutputDir,'README.md'),strjoin(text,newline));
end
