% isOpen  = exportToPPTX();
% if ~isempty(isOpen)
%     % If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
% 
% exportToPPTX('open',fullfile('presentacion','dieTemplate.pptx'));
%%
exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Base Monetaria','Position','title','fontsize',48);

%% 
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_MB.png'),...
            'Position', L2x2_UPPER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_MB_GAP.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_MB.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_MB.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 8.78/2.54 0.9/2.54 0.9/2.54]);
        
% simulación short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_MB_short.png'),...
            'Position', L2x2_UPPER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_MB_GAP_short.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_MB_short.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_MB_short.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 8.78/2.54 0.9/2.54 0.9/2.54]);
        
%% Contribuciones a Base Monetaria
if strcmp(folder_name{i}, 'v0')
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long',sprintf('D4L_MB_%s.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('D4L_MB_%s.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]); 
        
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short', sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short',sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        % Primera Diferencia Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short', sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short',sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
% Escenario
else
        % Long
        exportToPPTX('addslide');   
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long',sprintf('D4L_MB_%s.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('D4L_MB_%s.png', 'Alterno')),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  


        % short
        exportToPPTX('addslide');     
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short',sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE)),...
                    'Position', [0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short', sprintf('D4L_MB_%s_short.png', 'Alterno')),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);


        %Primera Diferencia Short
        exportToPPTX('addslide');    
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short',sprintf('D4L_MB_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short', sprintf('D4L_MB_%s_short.png', 'Alterno')),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  


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
                             'shock_dec','long','D4L_MB_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','D4L_MB_shd_dsc.png'),...
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
                             'shock_dec','short','D4L_MB_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','D4L_MB_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
                 
% Escenario
else
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','long','D4L_MB_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','long','D4L_MB_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]); 
                
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);        
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','short','D4L_MB_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
    
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','short','D4L_MB_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
                
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
                             'shock_dec','diff','D4L_MB_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','diff','D4L_MB_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

% Escenario
else
    % short
        exportToPPTX('addslide');

        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','diff','D4L_MB_diff_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','diff','D4L_MB_diff_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  


end

%% Base Monetaria en quetzales
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','MB.png'),...
            'Position', L1x2_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_MB.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_MB.png'),...
            'Position', L2x2_LOWER_RIGHT);
% simulación short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','MB_short.png'),...
            'Position', L1x2_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_MB_short.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_MB_short.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
%% Shock de base monetaria
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','PM_D4L_MB.png'),...
            'Position',[0 0 16.93/2.54 17.65/2.54]);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','PM_D4L_MB_short.png'),...
            'Position',[16.94/2.54 0 16.93/2.54 17.65/2.54]);   


%% Velocidad de Circulación de dinero
exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Velocidad de Circulación de la Base Monetaria','Position','title','fontsize',48);
%% Velocidad SUBPLOT
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'otras','Velocidad (componentes).png'),...
            'Position', L1_SINGLE); 

% Velocidad subplot short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'otras','Velocidad (componentes)_short.png'),...
            'Position', L1_SINGLE); 
        
        
%%
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_VEL.png'),...
            'Position', L2x2_UPPER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_VEL_GAP.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_VEL.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_VEL.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 8.78/2.54 0.9/2.54 0.9/2.54]);
        
% simulación short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_VEL_short.png'),...
            'Position', L2x2_UPPER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_VEL_GAP_short.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_VEL_short.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_VEL_short.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position',[16.93/2.54 8.78/2.54 0.9/2.54 0.9/2.54]);

%% Componentes de la velocidad de circulacion
if strcmp(folder_name{i}, 'v0')
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_%s.png', MODEL.CORR_DATE_ANT)),...
                    'Position',[0 0 16.93/2.54 17.65/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_%s.png', MODEL.CORR_DATE)),...
                    'Position',[16.94/2.54 0 16.93/2.54 17.65/2.54]);

        % short
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_short_%s.png', MODEL.CORR_DATE_ANT)),...
                    'Position',[0 0 16.93/2.54 17.65/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_short_%s.png', MODEL.CORR_DATE)),...
                    'Position',[16.94/2.54 0 16.93/2.54 17.65/2.54]);
                
else
    
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_%s.png',corr_date)),...
                    'Position',[0 0 16.93/2.54 17.65/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_%s.png', 'Alterno' )),...
                    'Position',[16.94/2.54 0 16.93/2.54 17.65/2.54]);

        % short
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_short_%s.png', corr_date)),...
                    'Position',[0 0 16.93/2.54 17.65/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared', sprintf('Comp_VEL_short_%s.png', 'Alterno' )),...
                    'Position',[16.94/2.54 0 16.93/2.54 17.65/2.54]);
    
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
                             'shock_dec','long','D4L_VEL_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','D4L_VEL_shd_dsc.png'),...
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
                             'shock_dec','short','D4L_VEL_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','D4L_VEL_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
                 
else
    
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','long','D4L_VEL_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]); 
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','long','D4L_VEL_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]); 

        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','short','D4L_VEL_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','short','D4L_VEL_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);   
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
                             'shock_dec','diff','D4L_VEL_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','diff','D4L_VEL_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

% Escenario
else
    % short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, 'v0',... 
                             'shock_dec','diff','D4L_VEL_diff_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
                
        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date, folder_name{i},... 
                             'shock_dec','diff','D4L_VEL_diff_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  


end

%%
% exportToPPTX( ...
%     'save', ...
%     'prueba');
% exportToPPTX('close');  