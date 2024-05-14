% isOpen  = exportToPPTX();
% if ~isempty(isOpen)
%     % If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
% 
% exportToPPTX('open',fullfile('presentacion','dieTemplate.pptx'));

%% Encabezado

exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Tasa de Interés','Position','title','fontsize',48);

%% Tasa líder
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RS.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[0 0 0.9/2.54 0.9/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RS_short.png'),...
             'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);
         
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 0 0.9/2.54 0.9/2.54]);

%% Contribuciones Tasa de interes lider
if strcmp(folder_name{i}, 'v0')
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Contribuciones a Tasa de Interés de Política'),...
                     'Position',[8.47/2.54, 0.07/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('RS_%s.png', MODEL.CORR_DATE)),... 
                    'Position',[0 1.18/2.54 16.93/2.54 16.69/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'contributions','short',sprintf('RS_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[16.93/2.54 1.18/2.54 16.93/2.54 16.69/2.54]);  
% Escenario
else
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Contribuciones a Costos Marginales Reales'),...
                     'Position',[8.47/2.54, 0.13/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('RS_%s.png', 'Alterno')),... 
                    'Position',[0 1.18/2.54 16.93/2.54 16.69/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short',sprintf('RS_%s_short.png', 'Alterno')),...
                    'Position',[16.93/2.54 1.18/2.54 16.93/2.54 16.69/2.54]);  

end
        
%% Descomposición de shocks
if strcmp(folder_name{i}, 'v0')
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','long','RS_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','RS_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','short','RS_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','RS_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

% Escenario
else
    % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', 'Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','long','RS_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','long','RS_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', 'Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','short','RS_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','short','RS_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
end

%% Descomposición de shocks (primera diferencia)
if strcmp(folder_name{i}, 'v0')
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','diff','RS_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','diff','RS_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

% Escenario
else
    % short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', 'Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','diff','RS_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','diff','RS_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
end

%% Tasa Real
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[0 0 0.9/2.54 0.9/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR_short.png'),...
             'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]); 
         
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 0 0.9/2.54 0.9/2.54]);
         
%% Descomposición de shocks
if strcmp(folder_name{i}, 'v0')
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,'v0',... 
                             'shock_dec','long','RR_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','RR_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','short','RR_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','RR_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
                 
% Escenario
else
% Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','long','RR_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','long','RR_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','short','RR_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','short','RR_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
end

%% Descomposición de shocks (primera diferencia)
if strcmp(folder_name{i}, 'v0')
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','diff','RR_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','diff','RR_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

% Escenario
else
    % short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', 'Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','diff','RR_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','diff','RR_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
end

%% componentes R
if strcmp(folder_name{i}, 'v0')
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_R_%s.png', corr_date)),...
                    'Position',[0 0 16.93/2.54 17.57/2.54]);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_R_%s.png', MODEL.CORR_DATE)),...
                     'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);

        % Componentes R short
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_R_short_%s.png', corr_date)),...
                    'Position',[0 0 16.93/2.54 17.57/2.54]);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_R_short_%s.png', MODEL.CORR_DATE)),...
                     'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);
                 
% Escenario
else
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared',sprintf('Comp_R_%s.png', 'Alterno')),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared',sprintf('Comp_R_%s.png', corr_date)),...
             'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);
         
% Componentes R short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared',sprintf('Comp_R_short_%s.png', 'Alterno')),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared',sprintf('Comp_R_short_%s.png', corr_date)),...
             'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);     
end

%% Tasa de interés neutral de Política
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RSNEUTRAL.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RSNEUTRAL_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);

%% Tasa de interés real de política (Tendencia)
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR_BAR.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR_BAR_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);

%% Tasa de interés real de política (Brecha)
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR_GAP.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','RR_GAP_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);
        
%%
% exportToPPTX( ...
%     'save', ...
%     'prueba');
% exportToPPTX('close');  