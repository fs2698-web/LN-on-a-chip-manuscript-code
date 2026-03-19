%% ==================== Data: T vs B activation (D1–D12) ====================
donor = ["D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11"]';

T = [446.94019, 171.20830, 128.58540, 153.50650, 169.02610, ...
     304.07970,  75.25710,  45.88367, 182.87610, ...
      57.47656,  52.72571]';
B = [545.04902, 378.10980, 179.56870, 167.54860, 127.31720, ...
     239.41350, 125.72150,  70.11630, 193.51800, ...
      84.48749,  51.86272]';

group_raw = ["High","High","High","High","Mid","Mid","Mid","Mid","Low","Low","Low"]';
group = group_raw; group(group_raw=="Mid") = "Intermediate";

x = T;  y = B;

%% ==================== Colors ====================
cHigh = [0.35 0.60 0.85];   
cMid  = [0.95 0.70 0.40];   
cLow  = [0.85 0.30 0.35];   
cFit  = [0.25 0.25 0.25];
isHigh = group=="High";
isMid  = group=="Intermediate";
isLow  = group=="Low";

%% ==================== Figure ====================
figure('Color','w','Units','pixels','Position',[80 80 950 720]); hold on;


sH = scatter(x(isHigh), y(isHigh), 300, 'filled', ...
    'MarkerFaceColor', cHigh, 'MarkerEdgeColor','k', 'DisplayName','High');
sM = scatter(x(isMid),  y(isMid),  300, 'filled', ...
    'MarkerFaceColor', cMid,  'MarkerEdgeColor','k', 'DisplayName','Intermediate');
sL = scatter(x(isLow),  y(isLow),  300, 'filled', ...
    'MarkerFaceColor', cLow,  'MarkerEdgeColor','k', 'DisplayName','Low');

dx = 0.02*range(x); dy = 0.02*range(y);
text(x+dx, y+dy, donor, 'FontSize',9, 'Color','k', 'Clipping','on');

%% ==================== Linear regression + 95% CI ====================
mdl = fitlm(x, y);                                   
xg = linspace(min(x)-0.05*range(x), max(x)+0.05*range(x), 400)';
[yg, yCI] = predict(mdl, xg, 'Alpha', 0.05);          

ph = patch([xg; flipud(xg)], [yCI(:,1); flipud(yCI(:,2))], cFit, ...
    'FaceAlpha', 0.15, 'EdgeColor','none', 'DisplayName','95% CI');
pl = plot(xg, yg, 'LineWidth', 2.2, 'Color', cFit, 'DisplayName','Linear fit');

%% ==================== Stats on plot (Modified for Image Style) ====================
R2_val = mdl.Rsquared.Ordinary; 
p_val = mdl.Coefficients.pValue(2);

statTxt = sprintf('R^2 = %.4f\nP = %.4f', R2_val, p_val);

xl = xlim; yl = ylim;
text(xl(1)+0.70*(xl(2)-xl(1)), yl(1)+0.15*(yl(2)-yl(1)), statTxt, ...
     'FontSize',14, 'FontWeight','normal', 'Color','k', 'FontName', 'Arial');

%% ==================== Axes & style ====================
xlabel('T cell activation (a.u.)');
ylabel('B cell activation (a.u.)');
ax = gca;
ax.LineWidth = 1.1; ax.FontSize = 11; ax.TickDir = 'out'; ax.Box = 'off';


grid off; 

xlim([min(x)-0.05*range(x), max(x)+0.10*range(x)]);
ylim([min(0, min(y)-0.05*range(y)), max(y)+0.10*range(y)]);


legend([sH sM sL pl ph], {'High','Intermediate','Low','Linear fit','95% CI'}, ...
       'Location','northeastoutside','Box','off');

title('T cell activation vs. B cell activation');

%% ==================== Export ====================
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));

outputDir = fullfile(scriptDir, '..', 'output_example');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

saveName = 'T_vs_B_activation.png';
savePath = fullfile(outputDir, saveName);

print(gcf, '-dpng', '-r600', savePath);

fprintf('file：\n<a href="matlab:winopen(''%s'')">%s</a>\n', outputDir, savePath);
