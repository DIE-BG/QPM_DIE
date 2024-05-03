% isOpen  = exportToPPTX();
% if ~isempty(isOpen)
%     % If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
% 
% exportToPPTX('open',fullfile('presentacion','dieTemplate.pptx'));

%% Encabezado

exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Crecimiento Interno','Position','title','fontsize',48);

%% Producto Interno
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PreProcessing','q_GDP.png'),... 
            'Position',[0 0 16.93/2.54 17.57/2.54]); 

% Producto Interno short
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PreProcessing','q_GDP_short.png'),... 
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]); 
        
%% Simulacion
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','L_GDP.png'),...
            'Position',[0 0 16.93/2.54 8.79/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','L_GDP_GAP.png'),...
            'Position',[0 8.78/2.54 16.93/2.54 8.79/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[0 8.78/2.54 0.9/2.54 0.9/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_GDP.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 8.79/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4_GDP_SM.png'),...
            'Position',[16.93/2.54 8.78/2.54 16.93/2.54 8.79/2.54]);
% simulación short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_GDP_short.png'),...
            'Position',[0 0 16.93/2.54 8.79/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_GDP_GAP_short.png'),...
            'Position',[0 8.78/2.54 16.93/2.54 8.79/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[0 8.78/2.54 0.9/2.54 0.9/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_GDP_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 8.79/2.54]);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4_GDP_SM_short.png'),...
            'Position',[16.93/2.54 8.78/2.54 16.93/2.54 8.79/2.54]);

%% Descomposición de shocks
if strcmp(folder_name{i}, 'v0')
% Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','long','L_GDP_GAP_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','L_GDP_GAP_shd_dsc.png'),...
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
                             'shock_dec','short','L_GDP_GAP_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','L_GDP_GAP_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
                 
% Escenario
else
 exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Escenario'),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','long','L_GDP_GAP_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','long','L_GDP_GAP_shd_dsc.png'),...
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
                             'shock_dec','short','L_GDP_GAP_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','short','L_GDP_GAP_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
end

%% Indice de Condiciones Monetaria
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','MCI.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','MCI_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);

%% Remesas
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','REM_GDP.png'),...
            'Position',[0 0 16.93/2.54 17.57/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','REM_GDP_short.png'),...
            'Position',[16.93/2.54 0 16.93/2.54 17.57/2.54]);
%%
% exportToPPTX( ...
%     'save', ...
%     'prueba');
% exportToPPTX('close');  