% isOpen  = exportToPPTX();
% if ~isempty(isOpen)
%     % If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
% 
% exportToPPTX('open',fullfile('presentacion','dieTemplate.pptx'));

%% Encabezado TCN

exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Tipo de Cambio Nominal','Position','title','fontsize',48);


%%

exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','L_S.png'),...
            'Position', L2x2_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position', STAR_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','S.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_S.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_S.png'),...
            'Position', L2x2_LOWER_RIGHT);
        
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','L_S_short.png'),...
            'Position', L2x2_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position', STAR_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','S_short.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_S_short.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_S_short.png'),...
            'Position', L2x2_LOWER_RIGHT);

%% Contribuciones a Tipo de Cambio Nominal
if strcmp(folder_name{i}, 'v0')
        % Long
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('L_S_%s.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long',sprintf('L_S_%s.png', MODEL.CORR_DATE)),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]); 
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short', sprintf('L_S_%s_short.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short',sprintf('L_S_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        % Primera Diferencia Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                  
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short', sprintf('L_S_%s_short.png', MODEL.CORR_DATE_ANT)),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short',sprintf('L_S_%s_short.png', MODEL.CORR_DATE)),...
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
                             'contributions','long',sprintf('L_S_%s.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);   

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','long', sprintf('L_S_%s.png', 'Alterno')),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);               
        % short
        exportToPPTX('addslide');       
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short',sprintf('L_S_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);   

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'contributions','short', sprintf('L_S_%s_short.png', 'Alterno')),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);   
                
        %Primera Diferencia Short
        exportToPPTX('addslide');       
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short',sprintf('L_S_%s_short.png', MODEL.CORR_DATE)),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);   

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE, folder_name{i},... 
                             'diff_contributions','short', sprintf('L_S_%s_short.png', 'Alterno')),... 
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
                             'shock_dec','long','D4L_S_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','long','D4L_S_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  

        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','short','D4L_S_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);                
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','short','D4L_S_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  

%escenario
else
        % Long                         
        exportToPPTX('addslide');
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);        
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,'v0',... 
                             'shock_dec','long','D4L_S_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','long','D4L_S_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  


        % Short
        exportToPPTX('addslide');
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);        
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,'v0',... 
                             'shock_dec','short','D4L_S_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','short','D4L_S_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
end

%% Descomposición de shocks (primeras diferencias)

if strcmp(folder_name{i}, 'v0')
        % Short
        exportToPPTX('addslide');
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_ant),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','diff','D4L_S_diff_shd_dsc.png'),... 
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'shock_dec','diff','D4L_S_diff_shd_dsc.png'),...
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);

%escenario
else
        % Short
        exportToPPTX('addslide');
        
        exportToPPTX('addtext',sprintf('Corrimiento %s', MODEL.leg_act),...
                     'Position',[0/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18);        
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,'v0',... 
                             'shock_dec','diff','D4L_S_diff_shd_dsc.png'),...
                    'Position',[0 2.51/2.54 16.93/2.54 14.98/2.54]);

        exportToPPTX('addtext',MODEL.Esc.(folder_name{i}).name,...
                     'Position',[17.04/2.54, 1.2/2.54, 16.93/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',18); 
        exportToPPTX('addpicture',...
                    fullfile('plots',corr_date,folder_name{i},... 
                             'shock_dec','diff','D4L_S_diff_shd_dsc.png'),... 
                    'Position',[16.93/2.54 2.51/2.54 16.93/2.54 14.98/2.54]);  
end

%% PREM
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','PREM.png'),...
            'Position', L1x2_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','PREM_short.png'),...
            'Position', L1x2_RIGHT);

exportToPPTX('addtext',{'*Generada a partir de su ecuación en el modelo, utilizando el filtro de Kalman.'},...
                     'Position',[0/2.54, 17.94/2.54, 18.21/2.54, 0.77/2.54],...
                     'HorizontalAlignment', 'center','fontsize',11);

%% Encabezado IPE_Q
exportToPPTX('addslide','Layout','Encabezado de sección');
exportToPPTX('addtext','Precio de transables -Quetzales-','Position','title','fontsize',48);

%% IPEI QUETZALIZADO
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_IPEI_Q.png'),...
            'Position', L2x2_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position', STAR_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_IPEI_Q_GAP.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_IPEI_Q.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_IPEI_Q.png'),...
            'Position', L2x2_LOWER_RIGHT);

%% IPEI_Q short
exportToPPTX('addslide');
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_IPEI_Q_short.png'),...
            'Position', L2x2_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('src','presentacion','star.png'),...
            'Position', STAR_UPPER_LEFT);
        
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'PostProcessing','L_IPEI_Q_GAP_short.png'),...
            'Position', L2x2_LOWER_LEFT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','DLA_IPEI_Q_short.png'),...
            'Position', L2x2_UPPER_RIGHT);
exportToPPTX('addpicture',...
            fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                     'prediction_compared','D4L_IPEI_Q_short.png'),...
            'Position', L2x2_LOWER_RIGHT); 

%% Componentes IPEI_Q
if strcmp(folder_name{i}, 'v0')
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_%s.png', corr_date)),...
                    'Position', L1x2_LEFT);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_%s.png', MODEL.CORR_DATE)),...
                     'Position', L1x2_RIGHT);

        % Componentes Inflación short
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_short_%s.png', corr_date)),...
                    'Position', L1x2_LEFT);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_short_%s.png', MODEL.CORR_DATE)),...
                     'Position', L1x2_RIGHT);
                 
% Escenario
else
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_%s.png', 'Alterno')),...
                    'Position', L1x2_RIGHT);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_%s.png', corr_date)),...
                     'Position', L1x2_LEFT);

        % Componentes Inflación short
        exportToPPTX('addslide');
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_short_%s.png', 'Alterno')),...
                    'Position', L1x2_RIGHT);
        exportToPPTX('addpicture',...
                    fullfile('plots',MODEL.CORR_DATE,folder_name{i},... 
                             'prediction_compared',sprintf('Comp_IPEI_Q_short_%s.png', corr_date)),...
                     'Position', L1x2_LEFT);
end
%%
% exportToPPTX( ...
%     'save', ...
%     'prueba');
% exportToPPTX('close');  