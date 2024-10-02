%% === Presentación larga de corrimientos === 
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
        
        disp('Presentación: Variables externas'); 
        VarExt;         disp('ok'); 
        disp('Presentación: Tipo de cambio real'); 
        TCReal;         disp('ok'); 
        disp('Presentación: Tipo de cambio nominal'); 
        TCNominal;      disp('ok'); 
        disp('Presentación: Inflación'); 
        Inflaciones;    disp('ok'); 
        disp('Presentación: Tasas de interés'); 
        TasasInteres;   disp('ok'); 
        disp('Presentación: Crecimiento económico'); 
        Crecimiento;    disp('ok'); 
        disp('Presentación: Base monetaria y velocidad de circulación '); 
        BM_Vel;         disp('ok'); 

        %%
        exportToPPTX('addslide','Layout','Encabezado de sección');
        exportToPPTX('addtext','Muchas Gracias','Position','title','fontsize',48);

        %% Guardar y cerrar
        disp(['Generando presentación: ', folder_name{i}, ' ok']);  
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
end