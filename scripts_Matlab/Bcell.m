%% ---------------- Relative Paths ----------------
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));

baseDir = fullfile(scriptDir, '..');

inputDir = fullfile(baseDir, "data_demo", "demo_Bcell_subset");
ctrlFile    = fullfile(inputDir, "Ctrl 1.xlsx");
vaccineFile = fullfile(inputDir, "Vaccine 1.xlsx");
negFile     = fullfile(inputDir, "negative.xlsx");


outDir = fullfile(baseDir, "output_example", "demo_Bcell_subset", "CD27_CD38_GlobalThreshold_withNegativeCheck");

if ~exist(outDir, "dir"); mkdir(outDir); end

outXls = fullfile(outDir, "Global_and_Results.xlsx");
condList  = ["Ctrl","Vaccine"];
sheetCD27 = "Global_CD27";
sheetCD38 = "Global_CD38";


if ~exist(ctrlFile, 'file') || ~exist(vaccineFile, 'file') || ~exist(negFile, 'file')
    error("no files", inputDir);
end









outDir = fullfile(baseDir, "output_example", "bcell_subset");
if ~exist(outDir, "dir"); mkdir(outDir); end

outXls = fullfile(outDir, "Global_and_Results.xlsx");

condList  = ["Ctrl","Vaccine"];
sheetCD27 = "Global_CD27";
sheetCD38 = "Global_CD38";

%% ---------------- READ + POOL ----------------
T27 = table();  % pooled long format for CD27
T38 = table();  % pooled long format for CD38

Tctrl = readLongExcel_AllDays_DaySpace(ctrlFile);
Tvac  = readLongExcel_AllDays_DaySpace(vaccineFile);

% Ensure numeric
Tctrl.Q = double(Tctrl.Q);  Tctrl.Z = double(Tctrl.Z);
Tvac.Q  = double(Tvac.Q);   Tvac.Z  = double(Tvac.Z);

% Build pooled long tables
for i = 1:numel(condList)
    cond = condList(i);

    if cond == "Ctrl"
        Tall = Tctrl;
    else
        Tall = Tvac;
    end

    days = sort(unique(Tall.Day));

    % ---- CD27 pooled ----
    for d = 1:numel(days)
        dayVal = days(d);

        x = Tall.Q(Tall.Day == dayVal);
        x = x(~isnan(x));
        if isempty(x), continue; end

        tmp = table( ...
            repmat(cond, numel(x), 1), ...
            repmat(dayVal, numel(x), 1), ...
            repmat("CD27", numel(x), 1), ...
            x, ...
            'VariableNames', {'Condition','Day','Marker','Intensity'} ...
        );
        T27 = [T27; tmp]; %#ok<AGROW>
    end

    % ---- CD38 pooled ----
    for d = 1:numel(days)
        dayVal = days(d);

        y = Tall.Z(Tall.Day == dayVal);
        y = y(~isnan(y));
        if isempty(y), continue; end

        tmp = table( ...
            repmat(cond, numel(y), 1), ...
            repmat(dayVal, numel(y), 1), ...
            repmat("CD38", numel(y), 1), ...
            y, ...
            'VariableNames', {'Condition','Day','Marker','Intensity'} ...
        );
        T38 = [T38; tmp]; %#ok<AGROW>
    end
end

%% ---------------- GLOBAL THRESHOLDS FROM REAL DATA ----------------
mu27 = mean(T27.Intensity);
sd27 = std(T27.Intensity);
thr27 = mu27;
thr27_2sd = mu27 + 2*sd27;   %#ok<NASGU>  % kept for record

mu38 = mean(T38.Intensity);
sd38 = std(T38.Intensity);
thr38 = mu38;
thr38_2sd = mu38 + 2*sd38;

summary27 = table(mu27, sd27, thr27, mu27 + 2*sd27, height(T27), ...
    'VariableNames', {'Mean','Std','Thr_Pos_Mean','Thr_MeanPlus2SD','N'});

summary38 = table(mu38, sd38, thr38, thr38_2sd, height(T38), ...
    'VariableNames', {'Mean','Std','Thr_Pos_Mean','Thr_MeanPlus2SD','N'});

%% ---------------- WRITE GLOBAL SHEETS ----------------
if exist(outXls, "file"); delete(outXls); end

writetable(summary27, outXls, 'Sheet', sheetCD27, 'Range', 'A1');
writetable(T27,       outXls, 'Sheet', sheetCD27, 'Range', 'A4');

writetable(summary38, outXls, 'Sheet', sheetCD38, 'Range', 'A1');
writetable(T38,       outXls, 'Sheet', sheetCD38, 'Range', 'A4');

fprintf("✅ Global pooled + thresholds written to: %s\n", outXls);

%% ---------------- READ NEGATIVE CONTROL ----------------
Tneg = readtable(negFile);

needNeg = ["CD27","CD38"];
if ~all(ismember(needNeg, string(Tneg.Properties.VariableNames)))
    error("negative.xlsx must contain columns: CD27 and CD38");
end

neg27 = double(Tneg.CD27);
neg38 = double(Tneg.CD38);

neg27 = neg27(~isnan(neg27));
neg38 = neg38(~isnan(neg38));

if isempty(neg27) || isempty(neg38)
    error("negative.xlsx has empty CD27 or CD38 values.");
end

%% ---------------- CHECK WHETHER NEGATIVE FITS GLOBAL THRESHOLDS ----------------
% Compare negative distribution against REAL-data thresholds
neg_mu27 = mean(neg27);
neg_sd27 = std(neg27);
neg_mu38 = mean(neg38);
neg_sd38 = std(neg38);

% Fraction of negative events above REAL thresholds
neg27_above_thr27      = sum(neg27 > thr27);
neg27_above_thr27_pct  = 100 * neg27_above_thr27 / numel(neg27);

neg38_above_thr38      = sum(neg38 > thr38);
neg38_above_thr38_pct  = 100 * neg38_above_thr38 / numel(neg38);

neg38_above_thr38pp     = sum(neg38 > thr38_2sd);
neg38_above_thr38pp_pct = 100 * neg38_above_thr38pp / numel(neg38);

% Optional interpretation
if neg27_above_thr27_pct < 5
    CD27_judgement = "Good";
elseif neg27_above_thr27_pct < 10
    CD27_judgement = "Acceptable";
else
    CD27_judgement = "Too high in negative";
end

if neg38_above_thr38_pct < 5
    CD38_judgement = "Good";
elseif neg38_above_thr38_pct < 10
    CD38_judgement = "Acceptable";
else
    CD38_judgement = "Too high in negative";
end

if neg38_above_thr38pp_pct < 1
    CD38pp_judgement = "Good";
elseif neg38_above_thr38pp_pct < 5
    CD38pp_judgement = "Acceptable";
else
    CD38pp_judgement = "Too high in negative";
end

negCheck = table( ...
    ["CD27"; "CD38"; "CD38++"], ...
    [neg_mu27; neg_mu38; neg_mu38], ...
    [neg_sd27; neg_sd38; neg_sd38], ...
    [thr27; thr38; thr38_2sd], ...
    [neg27_above_thr27; neg38_above_thr38; neg38_above_thr38pp], ...
    [neg27_above_thr27_pct; neg38_above_thr38_pct; neg38_above_thr38pp_pct], ...
    [CD27_judgement; CD38_judgement; CD38pp_judgement], ...
    'VariableNames', {'Marker','NegativeMean','NegativeStd','RealThreshold','NegAboveThreshold_n','NegAboveThreshold_pct','Judgement'});

negCheckFile = fullfile(outDir, "Negative_Threshold_Check.xlsx");
writetable(negCheck, negCheckFile, 'Sheet', 'Negative_Check', 'Range', 'A1');

%% ---------------- PLOT NEGATIVE vs REAL HISTOGRAMS ----------------
histFile = fullfile(outDir, "Histogram_Negative_vs_Real.png");

real27 = T27.Intensity;
real38 = T38.Intensity;

real27 = real27(~isnan(real27));
real38 = real38(~isnan(real38));

% unified x limits
xMax27 = max([real27; neg27]);
xMax38 = max([real38; neg38]);

xMax27 = ceil(xMax27 * 1.05);
xMax38 = ceil(xMax38 * 1.05);

figH = figure('Visible','off'); clf;

cNeg  = [0.45 0.65 0.85];   % blue
cReal = [0.95 0.55 0.60];   % pink/red

% -------- CD27 --------
subplot(1,2,1); hold on; box on;
histogram(neg27,  'Normalization','probability', ...
    'FaceColor', cNeg,  'EdgeColor','none', 'FaceAlpha',0.60);
histogram(real27, 'Normalization','probability', ...
    'FaceColor', cReal, 'EdgeColor','none', 'FaceAlpha',0.45);

% lines from REAL thresholds
xline(thr27, '--k', 'Real mean', ...
    'LineWidth',1.2, 'LabelOrientation','aligned', ...
    'LabelVerticalAlignment','middle');

xline(mu27 + 2*sd27, '--', 'Real mean+2SD', ...
    'Color',[0.85 0.20 0.20], 'LineWidth',1.2, ...
    'LabelOrientation','aligned', 'LabelVerticalAlignment','middle');

xlabel('CD27 Intensity (a.u.)');
ylabel('Normalized Frequency');
title(sprintf('CD27 | Neg vs Real | Neg>thr = %.1f%%', neg27_above_thr27_pct), 'Interpreter','none');
legend({'Negative Control','Real Data'}, 'Location','northeast');
xlim([0, max(1, xMax27)]);
ylim([0, 0.4]);

% -------- CD38 --------
subplot(1,2,2); hold on; box on;
histogram(neg38,  'Normalization','probability', ...
    'FaceColor', cNeg,  'EdgeColor','none', 'FaceAlpha',0.60);
histogram(real38, 'Normalization','probability', ...
    'FaceColor', cReal, 'EdgeColor','none', 'FaceAlpha',0.45);

xline(thr38, '--k', 'Real mean', ...
    'LineWidth',1.2, 'LabelOrientation','aligned', ...
    'LabelVerticalAlignment','middle');

xline(thr38_2sd, '--', 'Real mean+2SD', ...
    'Color',[0.85 0.20 0.20], 'LineWidth',1.2, ...
    'LabelOrientation','aligned', 'LabelVerticalAlignment','middle');

xlabel('CD38 Intensity (a.u.)');
ylabel('Normalized Frequency');
title(sprintf('CD38 | Neg vs Real | Neg>thr = %.1f%% | Neg>thr++ = %.1f%%', ...
    neg38_above_thr38_pct, neg38_above_thr38pp_pct), 'Interpreter','none');
legend({'Negative Control','Real Data'}, 'Location','northeast');
xlim([0, max(1, xMax38)]);
ylim([0, 0.4]);

set(figH, 'Position', [100 100 1200 360]);
saveas(figH, histFile);
close(figH);

disp("✅ Negative check table saved:");
disp(negCheckFile);
disp("✅ Negative vs Real histogram saved:");
disp(histFile);

%% ============================================================
% APPLY GLOBAL THRESHOLDS -> per Condition × Day scatter + percentages
%% ============================================================

plotDir = fullfile(outDir, "Scatter_byConditionDay");
if ~exist(plotDir, "dir"); mkdir(plotDir); end

maxPts = 20000;
results = table();

% Global scatter axis limits
allX = [Tctrl.Q; Tvac.Q; neg27];
allY = [Tctrl.Z; Tvac.Z; neg38];
allX = allX(~isnan(allX));
allY = allY(~isnan(allY));

xMax = ceil(max(allX) * 1.05);
yMax = ceil(max(allY) * 1.05);

xLimGlobal = [0, max(1, xMax)];
yLimGlobal = [0, max(1, yMax)];

for i = 1:numel(condList)
    cond = condList(i);

    if cond == "Ctrl"
        Tall = Tctrl;
    else
        Tall = Tvac;
    end

    days = sort(unique(Tall.Day));

    for d = 1:numel(days)
        dayVal = days(d);

        x = Tall.Q(Tall.Day == dayVal);
        y = Tall.Z(Tall.Day == dayVal);

        ok = ~isnan(x) & ~isnan(y);
        x = x(ok); y = y(ok);

        n = numel(x);
        if n == 0, continue; end

        % ---------- Gates ----------
        CD27pos = x > thr27;
        CD38pos = y > thr38;
        CD38pp  = y > thr38_2sd;

        is_naive  = ~CD27pos & ~CD38pos;
        is_preGC  = ~CD27pos &  CD38pos;
        is_memory =  CD27pos & ~CD38pos;
        is_plasma =  CD27pos &  CD38pp;
        is_GC     =  CD27pos &  CD38pos & ~CD38pp;

        pct = @(m) 100*sum(m)/n;
        pctNaive  = pct(is_naive);
        pctPreGC  = pct(is_preGC);
        pctGC     = pct(is_GC);
        pctMemory = pct(is_memory);
        pctPlasma = pct(is_plasma);

        results = [results; table(cond, dayVal, n, ...
            pctNaive, pctPreGC, pctGC, pctMemory, pctPlasma, ...
            'VariableNames', {'Condition','Day','Ncells','Naive','PreGC','GC','Memory','Plasma'})]; %#ok<AGROW>

        % ---------- Scatter ----------
        idx = 1:n;
        if n > maxPts
            idx = randperm(n, maxPts);
        end

        x_s = x(idx); y_s = y(idx);
        is_naive_s  = is_naive(idx);
        is_preGC_s  = is_preGC(idx);
        is_GC_s     = is_GC(idx);
        is_memory_s = is_memory(idx);
        is_plasma_s = is_plasma(idx);

        fig = figure('Visible','off'); clf; hold on;

        scatter(x_s(is_naive_s),  y_s(is_naive_s),  14, 'filled');
        scatter(x_s(is_preGC_s),  y_s(is_preGC_s),  14, 'filled');
        scatter(x_s(is_GC_s),     y_s(is_GC_s),     14, 'filled');
        scatter(x_s(is_memory_s), y_s(is_memory_s), 14, 'filled');
        scatter(x_s(is_plasma_s), y_s(is_plasma_s), 14, 'filled');

        xline(thr27, '--k', 'CD27+ (mean)');
        yline(thr38, '--b', 'CD38+ (mean)');
        yline(thr38_2sd, '-.r', 'CD38++ (mean+2SD)');

        xlabel('CD27');
        ylabel('CD38');
        title(sprintf('%s | Day %d | CD27 vs CD38 (n=%d)', cond, dayVal, n), 'Interpreter','none');
        grid on; box on;

        xlim(xLimGlobal);
        ylim(yLimGlobal);

        txt = sprintf(['Naive %.1f%%\npreGC %.1f%%\nGC %.1f%%\nMemory %.1f%%\nPlasma %.1f%%'], ...
            pctNaive, pctPreGC, pctGC, pctMemory, pctPlasma);

        ax = gca;
        xL = ax.XLim; 
        yL = ax.YLim;

        text(xL(1)+0.02*range(xL), yL(2)-0.05*range(yL), txt, ...
            'FontSize', 11, 'VerticalAlignment','top', 'BackgroundColor','w');

        legend({'Naive','preGC','GC','Memory','Plasma'}, 'Location','northwest');

        txt2 = sprintf([ ...
            'XLim = [%.0f, %.0f]\nYLim = [%.0f, %.0f]\n\n' ...
            'Global thr (used for gating)\nCD27 mean=%.1f\nCD27 mean+2SD=%.1f\nCD38 mean=%.1f\nCD38 mean+2SD=%.1f\n\n' ...
            'Negative check\nNeg CD27 > thr = %.1f%%\nNeg CD38 > thr = %.1f%%\nNeg CD38 > thr++ = %.1f%%' ...
            ], ...
            xLimGlobal(1), xLimGlobal(2), yLimGlobal(1), yLimGlobal(2), ...
            thr27, mu27 + 2*sd27, thr38, thr38_2sd, ...
            neg27_above_thr27_pct, neg38_above_thr38_pct, neg38_above_thr38pp_pct);

        annotation(fig, 'textbox', [0.72 0.12 0.26 0.78], ...
            'String', txt2, ...
            'FitBoxToText', 'off', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', [0.3 0.3 0.3], ...
            'FontSize', 10);

        saveas(fig, fullfile(plotDir, sprintf("Scatter_%s_Day%d.png", cond, dayVal)));
        close(fig);
    end
end

writetable(results, fullfile(outDir, "Bcell_SubsetPercent_byDay.csv"));
disp("✅ Scatter plots + percentages DONE.");
disp("Saved to: " + plotDir);

%% ============================================================
% EXPORT per-day subtype COUNTS + PERCENTAGES table
%% ============================================================

summaryRows = table();

for i = 1:numel(condList)
    cond = condList(i);

    if cond == "Ctrl"
        Tall = Tctrl;
    else
        Tall = Tvac;
    end

    days = sort(unique(Tall.Day));

    for d = 1:numel(days)
        dayVal = days(d);

        x = Tall.Q(Tall.Day == dayVal);
        y = Tall.Z(Tall.Day == dayVal);

        ok = ~isnan(x) & ~isnan(y);
        x = x(ok); y = y(ok);

        Ncells = numel(x);
        if Ncells == 0, continue; end

        CD27pos = x > thr27;
        CD38pos = y > thr38;
        CD38pp  = y > thr38_2sd;

        is_naive  = ~CD27pos & ~CD38pos;
        is_preGC  = ~CD27pos &  CD38pos;
        is_memory =  CD27pos & ~CD38pos;
        is_plasma =  CD27pos &  CD38pp;
        is_GC     =  CD27pos &  CD38pos & ~CD38pp;

        Naive_n  = sum(is_naive);
        PreGC_n  = sum(is_preGC);
        GC_n     = sum(is_GC);
        Memory_n = sum(is_memory);
        Plasma_n = sum(is_plasma);

        Naive_pct  = 100 * Naive_n  / Ncells;
        PreGC_pct  = 100 * PreGC_n  / Ncells;
        GC_pct     = 100 * GC_n     / Ncells;
        Memory_pct = 100 * Memory_n / Ncells;
        Plasma_pct = 100 * Plasma_n / Ncells;

        summaryRows = [summaryRows; table(cond, dayVal, Ncells, ...
            Naive_n, PreGC_n, GC_n, Memory_n, Plasma_n, ...
            Naive_pct, PreGC_pct, GC_pct, Memory_pct, Plasma_pct, ...
            'VariableNames', {'Condition','Day','Ncells', ...
                              'Naive_n','PreGC_n','GC_n','Memory_n','Plasma_n', ...
                              'Naive_pct','PreGC_pct','GC_pct','Memory_pct','Plasma_pct'})]; %#ok<AGROW>
    end
end

summaryRows.Condition = categorical(summaryRows.Condition);
summaryRows = sortrows(summaryRows, {'Condition','Day'});

csvOut = fullfile(outDir, "Bcell_Subtype_Counts_and_Percent_byDay.csv");
writetable(summaryRows, csvOut);

outSheetName = "Subtype_Summary_ByDay";
writetable(summaryRows, outXls, 'Sheet', outSheetName, 'Range', 'A1');

disp("✅ Subtype counts + percentages table exported:");
disp("CSV: " + csvOut);
disp("Excel sheet: " + outSheetName);

%% ======================= FUNCTIONS ==========================

function T = readLongExcel_AllDays_DaySpace(file)
% Read long-format Excel with day sheets like:
%   "Day 4", "Day_4", "Day-4", "day4"
% Columns required: Day Group Replicate Sheet Q Z

    if ~exist(file, "file")
        error("File not found: %s", file);
    end

    sh = sheetnames(file);
    sh = string(sh);
    shTrim = strtrim(sh);
    shLow  = lower(shTrim);

    isDaySheet = startsWith(shLow, "day") & ...
        ~cellfun(@isempty, regexp(cellstr(shLow), '^day[\s_\-]*\d+', 'once'));

    shDay = shTrim(isDaySheet);

    if isempty(shDay)
        error("No Day sheets found (expect 'Day 4' or 'Day_4' etc.) in: %s", file);
    end

    T = table();
    for i = 1:numel(shDay)
        sheetName = shDay(i);

        Ti = readtable(file, "Sheet", sheetName);

        needVars = ["Day","Group","Replicate","Sheet","Q","Z"];
        if ~all(ismember(needVars, string(Ti.Properties.VariableNames)))
            error("Sheet %s in %s missing required columns: %s", ...
                sheetName, file, strjoin(needVars, ", "));
        end

        Ti.Day = double(Ti.Day);
        Ti.Q   = double(Ti.Q);
        Ti.Z   = double(Ti.Z);
        Ti.Group = string(Ti.Group);
        Ti.Replicate = string(Ti.Replicate);
        Ti.Sheet = string(Ti.Sheet);

        keep = ~(isnan(Ti.Q) & isnan(Ti.Z));
        Ti = Ti(keep, :);

        if all(isnan(Ti.Day)) || isempty(Ti.Day)
            tok = regexp(lower(string(sheetName)), '^day[\s_\-]*(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                dayVal = str2double(tok{1});
                Ti.Day = repmat(dayVal, height(Ti), 1);
            end
        end

        T = [T; Ti]; %#ok<AGROW>
    end
end
