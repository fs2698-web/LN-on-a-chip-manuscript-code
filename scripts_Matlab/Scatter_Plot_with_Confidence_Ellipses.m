%% ========= Data: Flu AB (X) vs Clic (Y) =========
donor = ["D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11"]';

FluAB = [3.1984667, 3.1433333, 1.7913333, 3.3890000, 2.2640000, ...
         3.6350000, 2.0833333, 1.1573333, 1.2520000, ...
         2.2380000, 0.2313333]';
Clic  = [1.2710420, 1.1900770, 1.1657230, 1.1544400, 1.1327800, ...
         1.1175800, 1.1064000, 1.1038890, 1.0306760, ...
         1.0218720, 0.9746670]';

group_raw = ["High","High","High","High","Mid","Mid","Mid","Mid","Low","Low","Low"]';
group = group_raw; 
group(group_raw=="Mid") = "Intermediate";

x = FluAB; y = Clic;

%% ========= Colors & Settings =========
cHigh = [0.35 0.60 0.85];   
cMid  = [0.95 0.70 0.40];   
cLow  = [0.85 0.30 0.35];   
confLevel = 0.9; 

%% ========= Figure =========
figure('Color','w','Units','pixels','Position',[80 80 1000 750]); hold on;


draw_conf_ellipse(x(group=="High"),         y(group=="High"),         cHigh, confLevel);
draw_conf_ellipse(x(group=="Intermediate"), y(group=="Intermediate"), cMid,  confLevel);
draw_conf_ellipse(x(group=="Low"),          y(group=="Low"),          cLow,  confLevel);


sH = scatter(x(group=="High"),         y(group=="High"),         200, 'filled', ...
             'MarkerFaceColor',cHigh,'MarkerEdgeColor','k','DisplayName','High');
sM = scatter(x(group=="Intermediate"), y(group=="Intermediate"), 200, 'filled', ...
             'MarkerFaceColor',cMid, 'MarkerEdgeColor','k','DisplayName','Intermediate');
sL = scatter(x(group=="Low"),          y(group=="Low"),          200, 'filled', ...
             'MarkerFaceColor',cLow, 'MarkerEdgeColor','k','DisplayName','Low');


text(x + 0.05, y + 0.002, donor, 'FontSize',10, 'FontWeight','bold');


xlabel('Flu specific antibody from chip A450 (OD)');
ylabel('Flu specific antibody in serum (fold change)');
ax = gca; ax.LineWidth=1.2; ax.FontSize=12; ax.TickDir='out'; ax.Box='off';
grid off;


xlim([min(x)-1.2, max(x)+1.2]);
ylim([min(y)-0.15, max(y)+0.15]);

legend([sH sM sL], 'Location','northeastoutside', 'Box','off');
title(['Scatter Plot with ', num2str(confLevel*100), '% Confidence Ellipses']);

%% ========= Export =========

[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));


outputDir = fullfile(scriptDir, '..', 'output_example');


if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end


saveName = 'FluAB_vs_Clic_90_Ellipses.png';
savePath = fullfile(outputDir, saveName);
print(gcf, '-dpng', '-r600', savePath);


fprintf('save：\n<a href="matlab:winopen(''%s'')">%s</a>\n', outputDir, savePath);



function draw_conf_ellipse(x_in, y_in, faceColor, alphaLevel)
    if length(x_in) < 3, return; end 
    
    data = [x_in(:), y_in(:)];
    mu = mean(data);
    sigma = cov(data);
    
  
    [V, D] = eig(sigma);
    
    
    k = chi2inv(alphaLevel, 2); 
    
    
    t = linspace(0, 2*pi, 100);
    xy = [cos(t); sin(t)];
    
   
    k_scaled_xy = V * sqrt(D * k) * xy;
    
    
    ellipse_x = k_scaled_xy(1,:) + mu(1);
    ellipse_y = k_scaled_xy(2,:) + mu(2);
    
    
    patch(ellipse_x, ellipse_y, faceColor, ...
          'EdgeColor', faceColor, 'LineWidth', 2, ...
          'FaceAlpha', 0.1, 'HandleVisibility', 'off');
end