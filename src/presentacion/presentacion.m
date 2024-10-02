%% === Presentación larga de corrimientos === 

%% Disposiciones para las gráficas
% Single plot
L1_SINGLE           = [3.25/2.54 0 27.38/2.54 17.55/2.54]; 
% One-row layout
L1x2_LEFT           = [0 0 16.93/2.54 17.57/2.54]; 
L1x2_RIGHT          = [16.93/2.54 0 16.93/2.54 17.57/2.54];
% 2x2 layout
L2x2_UPPER_LEFT     = [0 0 16.93/2.54 8.79/2.54]; 
L2x2_UPPER_RIGHT    = [16.93/2.54 0 16.93/2.54 8.79/2.54]; 
L2x2_LOWER_LEFT     = [0 8.78/2.54 16.93/2.54 8.79/2.54]; 
L2x2_LOWER_RIGHT    = [16.93/2.54 8.78/2.54 16.93/2.54 8.79/2.54]; 
% Big left and 2x1 right
L3_LEFT             = L1x2_LEFT; 
L3_UPPER_RIGHT      = L2x2_UPPER_RIGHT; 
L3_LOWER_RIGHT      = L2x2_LOWER_RIGHT; 
% Star positions
STAR_UPPER_RIGHT    = [16.93/2.54 0 0.9/2.54 0.9/2.54]; 
STAR_UPPER_LEFT     = [0 0 0.9/2.54 0.9/2.54]; 
STAR_LOWER_LEFT     = [0 8.79/2.54 0.9/2.54 0.9/2.54]; 

%% Nombres de los escenarios
% Se generan las presentaciones de acuerdo con las carpetas de cada
% escenario en plots/

name = dir(fullfile('plots', MODEL.CORR_DATE));
folder_name = {};
for i = 1:length(name)
   if  ~strcmp(name(i).name, '.') && ~strcmp(name(i).name, '..')
       folder_name{end+1} = name(i).name;
   end
end

for i = 1:length(folder_name)
        %% Inicializacion de presentacion
        isOpen  = exportToPPTX();
        if ~isempty(isOpen)
            % If PowerPoint already started, then close first and then open a new one
            exportToPPTX('close');
        end

        disp(['Generando presentación: ', folder_name{i}]); 
        exportToPPTX('open',fullfile('presentacion','dieTemplate.pptx'));

        %% Diapositiva de título
        exportToPPTX('addslide','Master',1,'Layout','Diapositiva de título');
        exportToPPTX('addtext',{'Modelo QPM', sprintf('Corrimiento %s: %s',MODEL.leg_act, MODEL.esc_names{i})},...
                    'Position','Title','fontsize',45);
        exportToPPTX('addtext',{'DEPARTAMENTO DE INVESTIGACIONES ECONÓMICAS','BANCO DE GUATEMALA'},...
                     'Position','Subtitle','fontsize',20);
        exportToPPTX('addtext',sprintf('%s', dat2char(ddtoday)),...
                     'Position',[7.00/2.54, 14.45/2.54, 20.75/2.54, 1.05/2.54],...
                     'HorizontalAlignment', 'center','fontsize',15);
                 
        exportToPPTX('addtext','**NUMERAR DIPOSITIVAS**',...
                     'Position',[7.00/2.54, 1.05/2.54, 20.75/2.54, 2.11/2.54],...
                     'HorizontalAlignment', 'center','fontsize', 45);
                 
        %% CONTENIDO
        exportToPPTX('addslide','Layout','Título y objetos');
        exportToPPTX('addtext','CONTENIDO','Position',...
                     'Title','fontsize',36,...
                     'HorizontalAlignment','Left');
        exportToPPTX('addtext',{...
                    sprintf('\t Principales cambios en los datos'),...
                    sprintf('\t Resultados Principales'),...
                    sprintf('\t Variables Externas'),...
                    sprintf('\t Tipo de Cambio Real'),...
                    sprintf('\t Tipo de Cambio Nominal'),...
                    sprintf('\t Precio de Transables -Quetzales-'),...
                    sprintf('\t Inflación'),...            
                    sprintf('\t Tasa de Interés'),...
                    sprintf('\t Crecimiento Interno'),...
                    sprintf('\t Base Monetaria'),...
                    sprintf('\t Velocidad de Circulación de la Base Monetaria')},...
                    'Position','content','fontsize',22);%
                 
        %% Cambios
        exportToPPTX('addslide','Layout','Título y objetos');
        exportToPPTX('addtext',sprintf("%s: Cambios en datos", MODEL.esc_names{i}),'Position',...
                     'title','fontsize',36,...
                     'HorizontalAlignment','Left');
                 
        %% Resultados
        exportToPPTX('addslide','Layout','Título y objetos');
        exportToPPTX('addtext','Descripción', 'Position',...
                     'title','fontsize',36,...
                     'HorizontalAlignment','Left');
       
        %% secciones
        if strcmp(folder_name{i}, 'v0')
           corr_date = MODEL.CORR_DATE_ANT;
        else
           corr_date = MODEL.CORR_DATE;
        end
        
        fprintf('Presentación: Variables externas... '); 
        VarExt;         fprintf('ok.\n'); 
        fprintf('Presentación: Tipo de cambio real... '); 
        TCReal;         fprintf('ok.\n'); 
        fprintf('Presentación: Tipo de cambio nominal... '); 
        TCNominal;      fprintf('ok.\n'); 
        fprintf('Presentación: Inflación... '); 
        Inflaciones;    fprintf('ok.\n'); 
        fprintf('Presentación: Tasas de interés... '); 
        TasasInteres;   fprintf('ok.\n'); 
        fprintf('Presentación: Crecimiento económico... '); 
        Crecimiento;    fprintf('ok.\n'); 
        fprintf('Presentación: Base monetaria y velocidad de circulación... '); 
        BM_Vel;         fprintf('ok.\n'); 

        %% Cierre
        exportToPPTX('addslide','Layout','Encabezado de sección');
        exportToPPTX('addtext','Muchas Gracias','Position','title','fontsize',48);

        %% Guardar y cerrar
        
        save_path = fullfile('Resultados', MODEL.CORR_DATE);
        if ~isfolder(save_path)
            mkdir(save_path)
        end
        
        exportToPPTX( ...
            'save', ...
            fullfile(save_path, sprintf('QPM Corrimiento %s %s - %s',...
                                         MODEL.CORR_DATE,...
                                         folder_name{i},...
                                         MODEL.esc_names{i})));
        exportToPPTX('close');   
        disp(['Presentación escenario ', folder_name{i}, ' completa.']);  
end