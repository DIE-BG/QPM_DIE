function plot_shd_dsc_comparative(MODEL, HIST, SHD, varargin)

% plot_shd_dsc realiza las gráficas de descomposición de choques del modelo
% suavizado.
%{
% ## Syntax ##
%
%     plot_shd_dsc(MODEL, varargin)
%
% ## Input Arguments ##
%
% __`MODEL`__ [ struct ] -
% Debe contener al menos la estructura del modelo `MODEL.MF`,
% `MODEL.shd_dsc` y `MODEL.DATES`
%
%
% ## Options ##
%
% * Level = [ `true`|`false` ] - Gráfica a partir del estado estacionario.
%
% * SavePath = fullfile(userpath, 'temp') [ `String` ] - Directorio donde
% guarda la gráfica.
%
% * Variables = get(MODEL.MF, 'xlist') [ `cell` ] - Nombre de variable a
% graficar.
%
% ## Output Arguments ##
%
%
% ## Description ##
%
%
% ## Example ##
%
%}

% -DIE
% -Octubre 2021

% Parametros opcionales
 p = inputParser;
    addParameter(p, 'Level', true);
    addParameter(p, 'SavePath', fullfile(userpath, 'temp'));
    addParameter(p, 'Variables', get(MODEL.MF, 'xlist'));
    addParameter(p, 'CloseAll', true);
    addParameter(p, 'OnlyHist', false);
    addParameter(p, 'Rng', MODEL.DATES.hist_start:MODEL.DATES.pred_end);
    addParameter(p, 'Esc', {});
    addParameter(p, 'sub_title', {});
parse(p, varargin{:});
params = p.Results; 

if params.OnlyHist
   MODEL.shd_dsc = dbclip(MODEL.shd_dsc, MODEL.DATES.hist_start:MODEL.DATES.hist_end);
   MODEL.F_pred = dbclip(MODEL.F_pred, MODEL.DATES.hist_start:MODEL.DATES.hist_end);
elseif ~params.OnlyHist && ~isempty(params.Rng)
   MODEL.shd_dsc = dbclip(MODEL.shd_dsc, MODEL.DATES.hist_end - 20:MODEL.DATES.hist_end + 20);
   MODEL.F_pred = dbclip(MODEL.F_pred, MODEL.DATES.hist_end - 20:MODEL.DATES.hist_end + 20);
end


% Variables a descomponer
var_plot = params.Variables;

% Choques a las variables
var_shd = get(MODEL.MF, 'elist');

% Paleta de colores
col = distinguishable_colors(length(var_shd) + 1, ...
    'b', ...
    @(x) colorspace('RGB->Lab',x));

% Diferencias
SHD_esc = params.Esc{1};

% Recorte de fechas
SHD = databank.clip(SHD, MODEL.DATES.pred_start, MODEL.DATES.pred_end);
SHD_esc = databank.clip(SHD_esc, MODEL.DATES.pred_start, MODEL.DATES.pred_end);

for i = 1:length(var_plot)
   SHD.(var_plot{i}) = SHD.(var_plot{i}) - SHD_esc.(var_plot{i});
end

% Iteración para la descomposicón de shocks por shcok
for i = 1:length(var_plot)
    
    figure('Position', [1 42.0182 1.6756e+03 825.6000]);
    hold on
    
    % Barras
    barcon(SHD.(var_plot{i}){:, 1:end}, ...
        'dateFormat=','YYFP', ...
        'colorMap=', col, ...
        'evenlySpread=', false); 
    % Líneas
    %plot(real(HIST.(var_plot{i})),'w','LineWidth',2);
    %plot(real(HIST.(var_plot{i})),'k.-','LineWidth',1, 'MarkerSize', 10);
    % Línea vertical en el fin de historia
    grfun.vline(MODEL.DATES.hist_end,'timePosition','middle');
    % Título
    set(gca,'FontSize',12);       
    title({var_plot{i}, sprintf('Comparativo %s_%s', params.sub_title{1}, params.sub_title{2})}, 'Interpreter', 'none');
    % Leyendas
    legend(var_shd,'location','northeastoutside','FontSize',11, 'Interpreter', 'none')
    grid on;
    hold off
    
     
    if i < 10
        num = sprintf("0%i", i);
    else
        num = sprintf('%i', i);
    end
    saveas(gcf, ...
        fullfile(params.SavePath, ...
        sprintf("%s_shd_dsc_comparative_all.png", var_plot{i}))...
    )

    close all
end


% iteración para la descomposición de shocks internos y externos
for i = 1:length(var_plot)

    % Cálculo de contribuciones
    int = sum(SHD.(var_plot{i})(:,1:15),2);
    ext = sum(SHD.(var_plot{i})(:,16:22),2);
    
    SHD_group.(var_plot{i}) = tseries(SHD.(var_plot{i}).Range, [int ext]);
    
    figure('Position', [1 42.0182 1.6756e+03 825.6000]);
    hold on
    
    % Barras
    barcon(SHD_group.(var_plot{i}), ...
        'dateFormat=','YYFP', ...
        'colorMap=',[0.235 0.004 0.988; 0.004 0.812 0.3], ...
        'evenlySpread=', false); 
    grfun.vline(MODEL.DATES.hist_end,'timePosition','middle');
    % Título
    set(gca,'FontSize',12);       
    title({var_plot{i}, sprintf('Comparativo %s_%s', params.sub_title{1}, params.sub_title{2})}, 'Interpreter', 'none');
    % Leyendas
    legend({'Internos', 'Externos'},'location','northeastoutside','FontSize',11, 'Interpreter', 'none')
    grid on;
    hold off
    
     
    if i < 10
        num = sprintf("0%i", i);
    else
        num = sprintf('%i', i);
    end
    saveas(gcf, ...
        fullfile(params.SavePath, ...
        sprintf("%s_shd_dsc_comparative.png", var_plot{i}))...
    )
    close all
end



end