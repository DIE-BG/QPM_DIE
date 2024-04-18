%% presentacion

%% Elementos auxiliares

% Nombre del origen de los choques
temp_var_name = MODEL.FullListNames - MODEL.NoShockNames;
temp_var = MODEL.FullList - MODEL.NoShock;

temp_sname = cellfun( ...
    @(x) sprintf('s_%s', x), ...
    temp_var, ...
    'UniformOutput', false ...
);

temp_sname_var = cellfun( ...
    @(x) sprintf('s_%s', x), ...
    temp_var_name, ...
    'UniformOutput', false ...
);

% Nombre de las variables afectadas por el choque
xlist = MODEL.FullList;
xlistname = MODEL.FullListNames;

%% Inicializacion de presentacion
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end

exportToPPTX('open',fullfile('Models',MODEL.rootpath,'Initial_setup','presentacion','dieTemplate.pptx'));

%% Diapositiva de título
exportToPPTX('addslide','Master',1,'Layout','Diapositiva de título');
exportToPPTX('addtext',MODEL.rootpath,'Position','Título','fontsize',54);
exportToPPTX('addtext',['Múltiples Muestras de Estimación y Evaluación: ',...
            dat2char(MODEL.general_config.HistStart(1)),...
            ' a ' ,...
            dat2char(MODEL.general_config.HistStart(end))],...
            'Position','Subtitle');

exportToPPTX('addtext',...
            {"BANCO DE GUATEMALA",...
            "Departamento de Investigaciones Económicas"},...
            'Position',[1.3 0.6 7 0.85],'Vert','bottom','fontsize',14);
exportToPPTX('addtext',...
            "**Evaluación de Modelos Macroecómicos de Pronóstico y Análisis de Política Monetaria**",...
            'Position',[1.3 1.6 7 0.85],'fontsize',14);        


%% Diapositiva de contenido

exportToPPTX('addslide','Layout','Título y objetos');
exportToPPTX('addtext','CONTENIDO','Position',[1.5 0.2 7.40 1],'fontsize',34);
exportToPPTX('addtext',...
            {'- Informe de evaluación resumido',...
            ['- Modelo Evaluado: ',MODEL.rootpath],...            
            sprintf('\t\t- Descripción'),...
            sprintf('\t\t- Variables'),...
            sprintf('\t\t- Estimación'),...
            sprintf('\t\t- Criterios de Estado Estacionario'),...            
            sprintf('\t\t- Estabilidad de muestras'),...
            '- Evaluación de Causalidad',...
            sprintf('\t\t- Funciones de Impulso - Respuesta'),...
            sprintf('\t\t- Coeficientes de variación promedio'),...
            sprintf('\t\t- Descomposición de choques'),...
            '- Evaluación de Capacidad Predictiva',...
            sprintf('\t\t- Pronósticos'),...
            sprintf('\t\t- RMSE'),...                      
            sprintf('\t\t- Distribución de los errores'),...           
            '- ANEXOS'},...
            'Position',[1.3 1.1 11 5.5],'fontsize',18);
        
%% Informe de Evaluación resumido
exportToPPTX('addslide');
exportToPPTX('addtext','INFORME DE EVALUACIÓN RESUMIDO','Position',[1.5 4 7.40 1],'fontsize',48);

%% Informe de Evaluación resumido: CAPACIDAD EXPLICATIVA
exportToPPTX('addslide');
exportToPPTX('addtext','CAPACIDAD EXPLICATIVA','Position',[1 0.2 7.40 1],'fontsize',28);
exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\ir\CV_Prom\CV_PROM_',MODEL.rootpath,'.png'],'Position',[0.5 1 12 6]);

%% Informe de Evaluación resumido: CAPACIDAD PREDICTIVA
exportToPPTX('addslide');
exportToPPTX('addtext','CAPACIDAD PREDICTIVA','Position',[1 0.2 7.40 1],'fontsize',28);
% exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\cuadros_python\Resumen\Graf_Resumen1.PNG'],'Position',[2 1 10 6]);


%% MODELO EVALUADO:
exportToPPTX('addslide');
exportToPPTX('addtext',['MODELO: ',MODEL.rootpath],'Position',[1.5 4 7.40 1],'fontsize',48);

%% Diapositiva de descripcion del modelo

exportToPPTX('addslide')
exportToPPTX('addtext','DESCRIPCIÓN','Position',[1.5 1 7.40 1],'fontsize',40);
exportToPPTX('addtext','Colocar descripción AQUÍ','Position',[1.5 2.5 10 1],'fontsize',32);

%% Diapositivas con definiciones
exportToPPTX('addslide');
exportToPPTX('addtext','VARIABLES','Position',[1.5 1 7.40 1],'fontsize',40);
for i= 1:length(MODEL.listVar)    
    exportToPPTX('addtext',...
            ['- ',MODEL.FullListNames{i},': ',MODEL.listVar{i}],...
            'Position',[1.3 1.5+(0.5*i) 10 0.5],'Vert','bottom','fontsize',16);
end

pos_y=2.5+(0.5*length(MODEL.listVar));
exportToPPTX('addtext','- VARIABLES OBJETIVO:','Position',[1.3 pos_y 7.40 1],'fontsize',16);
for i= 1:length(MODEL.var_obj)    
    exportToPPTX('addtext',...
            [MODEL.obj_name{i}],...
            'Position',[(3+i) pos_y 4 0.5],'fontsize',16);
end    

%% Data
exportToPPTX('addslide');
exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\dbplot_',MODEL.rootpath,'.png'],'Position',[1 1 11 6]);

%% Diapositiva de descripcion del modelo

exportToPPTX('addslide')
exportToPPTX('addtext','ESTIMACIÓN','Position',[1.5 1 7.40 1],'fontsize',40);
exportToPPTX('addtext',...
            {'- Método de estimación: Mínimos Cuadrados Generalizados iterados (Eviews).',...
             '- Muestras de estimación: ',...
             '- Ajuste de estabilidad con factor de 0.99 en matriz de coeficientes rezagados.',...
             '- Muestra más pequeña: 20 trimestres ',...
             sprintf('\t\t- Backcast:'),...
             sprintf('\t\t- Forecast:')...
            },'Position',[1 2.5 11 1],'fontsize',32);

%% Diapositiva periodos de estimacion y evaluacion

exportToPPTX('addslide');
exportToPPTX('addtext','PERIODOS DE ESTIMACIÓN Y EVALUACIÓN','Position',[2.5 0.1 7.40 1],'fontsize',28);
exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'Initial_setup','presentacion','periodos.png'),'Position',[0.2 0.75 13 5.7]);

%% Criterios de estado estacionario
% Diapositiva de sección
exportToPPTX('addslide');
exportToPPTX('addtext','CRITERIOS DE ESTADO ESTACIONARIO','Position',[1.5 1 10 1],'fontsize',40);

exportToPPTX('addtext',...
            {...
            '**DESCRIBIR**',...
            },...
             'Position',[1.5 2.5 11 1],'fontsize',32);

%% Estabilidad
exportToPPTX('addslide');
exportToPPTX('addtext','ESTABILIDAD','Position',[1.5 1 7.40 1],'fontsize',40);

    exportToPPTX( ...
        'addpicture', ...
        fullfile('Models', MODEL.rootpath, 'results', 'plots', sprintf('Eigenvalues_%s_multiple.png', MODEL.rootpath)), ...
        'Position',[6.93/2.54 4.4/2.54 20/2.54 13.5/2.54] ...
    );

%% EVALUACIÓN DE CAUSALIDAD
% Diapositiva de sección
exportToPPTX('addslide');
exportToPPTX('addtext','EVALUACIÓN DE CAUSALIDAD','Position',[1.5 4 10.25 1],'fontsize',48);


%% Funciones impulso respuesta  
% Diapositiva de sección
exportToPPTX('addslide');
exportToPPTX('addtext','FUNCIONES IMPULSO - RESPUESTA','Position',[1.5 4 10.25 1],'fontsize',48);

for i=1:length(temp_sname)-1
    exportToPPTX('addslide');
    exportToPPTX('addtext',MODEL.listVar{i},'Position',[2.5 0.1 7.40 1],'fontsize',28);
    if i < 10
        exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','all',sprintf('0%s_s_v%s.png', dat2char(i),dat2char(i))),'Position',[0.2 1.75 13 5.7]);
    else
        exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','all',sprintf('%s_s_v%s.png', dat2char(i),dat2char(i))),'Position',[0.2 1.75 13 5.7]);    
    end
    
    for j= 1:length(xlist)-1
    exportToPPTX('addslide');
    exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','Distributions',sprintf('s_v%s', dat2char(i)),sprintf('s_v%s-v%s.png', dat2char(i),dat2char(j))),'Position',[0.2 0.75 13 5.7]);
    exportToPPTX('addslide');
    exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','Coef_var',sprintf('s_v%s', dat2char(i)),sprintf('s_v%s-v%s.png', dat2char(i),dat2char(j))),'Position',[0.2 0.75 13 5.7]);
    end
end

%% Coeficientes de Variación Promedio
exportToPPTX('addslide');
exportToPPTX('addtext',...
            {'Coeficientes de Variación Promedio',...
             'FUNCIONES IMPULSO - RESPUESTA'...
            },'Position',[1.5 4 10.25 1],'fontsize',48);

exportToPPTX('addslide');
exportToPPTX('addtext','Varibles Objetivo','Position',[1.5 4 7.40 1],'fontsize',32);

for var_obj = MODEL.var_obj
    exportToPPTX('addslide');    
    exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','CV_Prom',var_obj{:},sprintf('CVProm-%s.png', var_obj{:})),'Position',[0.2 0.75 13 5.7]);
    exportToPPTX('addslide'); 
    exportToPPTX('addpicture',fullfile('Models',MODEL.rootpath,'results','plots','ir','CV_Prom',var_obj{:},sprintf('CVPromAcum-%s.png', var_obj{:})),'Position',[0.2 0.75 13 5.7]);
end

exportToPPTX('addslide');
exportToPPTX('addtext','Promedio de todas las Variables Objetivo','Position',[1.5 4 7.40 1],'fontsize',32);

exportToPPTX('addslide');
exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\ir\CV_Prom\CV_PROM_',MODEL.rootpath,'.png'],'Position',[0.5 1 12 6]);

%% Descomposicion de choques
exportToPPTX('addslide');
exportToPPTX('addtext','DESCOMPOSICIÓN DE CHOQUES','Position',[1.5 4 10.25 1],'fontsize',48);

temp_date_name = sprintf( ...
                '%s-%s', ...
                dat2char(MODEL.general_config.HistStart(1)), ...
                dat2char(MODEL.general_config.HistStart(end)) ...
                );

for i = 1:length(xlist)-3
    exportToPPTX('addslide');
    exportToPPTX('addtext',MODEL.listVar{i},'Position',[1.5 0.5 7.40 1],'fontsize',28);
    exportToPPTX( ...
        'addpicture', ...
        ['Models\',MODEL.rootpath,'\results\plots\shdsc\',temp_date_name,'\',xlist{i}, '.png'], ...
        'Position',[2.65/2.54 3.81/2.54 11.75 14/2.54]);
end

exportToPPTX('addslide');
exportToPPTX( ...
        'addpicture', ...
        ['Models\',MODEL.rootpath,'\results\plots\shd\',temp_date_name,'\shd.png'], ...
        'Position',[2.00/2.54 1 11.75 14/2.54]);

%% EVALUACIÓN DE CAPACIDAD PREDICTIVA
% Diapositiva de sección
exportToPPTX('addslide');
exportToPPTX('addtext','EVALUACIÓN DE CAPACIDAD PREDICTIVA','Position',[1.5 4 10.25 1],'fontsize',48);

exportToPPTX('addslide');
exportToPPTX('addtext','PRONÓSTICOS','Position',[1.5 4 10.25 1],'fontsize',48);
%Gráficas de spaghetti
for j = 1:length(ScnG)
    for i = 1:length(VarGraph)
        exportToPPTX('addslide');    
        exportToPPTX( ...
            'addpicture', ...            
            ['Models\',MODEL.rootpath,'\results\plots\spaghetti\OS-Multiple\E',dat2char(j),'\','E',dat2char(j),'_',VarGraph{i},'.png'], ...
            'Position',[0.2 0.75 13 6]);
    end
end

%rmse
exportToPPTX('addslide');
exportToPPTX('addtext','RMSE','Position',[1.5 4 10.25 1],'fontsize',48);

exportToPPTX('addslide');
exportToPPTX('addtext','No. de observaciones por horizonte de pronóstico','Position',[1.5 1 10 1],'fontsize',40);

exportToPPTX('addtext',...
            {...
            '**copiar tablita de periodos.xlsx**',...
            },...
             'Position',[1.5 2.5 11 1],'fontsize',32);

exportToPPTX('addslide');
exportToPPTX('addtext','Variables OBJETIVO NO ANCLADAS','Position',[1 0.2 7.40 1],'fontsize',28);
% exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\cuadros_python\Horizontes\h_resumen_1.png'],'Position',[5 1 3.6 6]);

exportToPPTX('addslide');
exportToPPTX('addtext','**Ejercicio 1: RMSE**','Position',[1.3 0.3 7 0.75],'fontsize',24);
exportToPPTX('addtext','*Anclaje: Ninguno*','Position',[1.3 0.8 9 0.75],'fontsize',18);        
% exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\cuadros_python\Horizontes\h_select_01_1.png'],'Position',[5 1 3.6 6]);

exportToPPTX('addslide');
exportToPPTX('addtext','**Ejercicio 1: RMSE (Acumulado)**','Position',[1.3 0.3 7 0.75],'fontsize',24);
exportToPPTX('addtext','*Anclaje: Ninguno*','Position',[1.3 0.8 9 0.75],'fontsize',18);        
% exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\cuadros_python\Horizontes\h_select_01_2.png'],'Position',[5 1 3.6 6]);

exportToPPTX('addslide');
exportToPPTX('addtext','**Ejercicio 1: RMSE (Promedio Acumulado)**','Position',[1.3 0.3 7 0.75],'fontsize',24);
exportToPPTX('addtext','*Anclaje: Ninguno*','Position',[1.3 0.8 9 0.75],'fontsize',18);        
% exportToPPTX('addpicture',['Models\',MODEL.rootpath,'\results\plots\cuadros_python\Horizontes\h_select_01_3.png'],'Position',[5 1 3.6 6]);

%% Gráficas de errores

% Directorio temporal de las gráficas a importar.
temp_img_dir = fullfile('Models', MODEL.rootpath, 'results', 'plots', 'errors', 'OS-Multiple');

% Lista temporal del contenido dentro de la carpeta.
temp_dir_content = string( ...
    ls(...
        temp_img_dir ...
    )...
);

% Patron que deben seguir los nombres en las gráficas a importar.
temp_regular_pattern = '(?<=^)\S*\.jpg(?=\s*|$)';

% Lista de imagenes que cumplen el patron regular.
list = regexp(temp_dir_content, temp_regular_pattern, 'match');
list = list(cellfun(@(x) ~isempty(x), list));
list = cellfun(@(x) x, list);


exportToPPTX('addslide');
exportToPPTX('addtext','Distribución de los Errores','Position',[1.5 4 7.40 1],'fontsize',48);

for i = 1:length(list)
    exportToPPTX('addslide');    
        exportToPPTX( ...
        'addpicture', ...
        fullfile(temp_img_dir, char(list(i))), ...
        'Position',[2.5 1 8.5 5.9]);
end

% %% ANEXOS
% exportToPPTX('addslide');
% exportToPPTX('addtext','ANEXOS','Position',[1.5 4 7.40 1],'fontsize',48);
% 
% %% Ecuaciones del Modelo
% % Directorio temporal de las gráficas a importar.
% temp_img_dir = fullfile('Models', MODEL.rootpath, 'Initial_setup', 'presentacion', 'ecuaciones');
% 
% % Lista temporal del contenido dentro de la carpeta.
% temp_dir_content = lower(...
%     string( ...
%         ls(...
%             temp_img_dir ...
%         )...
%     )...
% );
% 
% % Patron que deben seguir los nombres en las gráficas a importar.
% temp_regular_pattern = '(?<=^)ec\d*\.png(?=\s*|$)';
% 
% % Lista de imagenes que cumplen el patron regular.
% temp_dir_content = regexp(temp_dir_content, temp_regular_pattern, 'match');
% temp_dir_content = temp_dir_content(cellfun(@(x) ~isempty(x), temp_dir_content));
% temp_dir_content = cellfun(@(x) x, temp_dir_content);
% 
% for i = 1:length(temp_dir_content)
%     if i == 1
%     exportToPPTX('addslide');
%     exportToPPTX('addtext','ECUACIONES DEL MODELO','Position',[1 0.2 7.40 1],'fontsize',38);
%     exportToPPTX( ...
%         'addpicture', ...
%         fullfile(temp_img_dir, char(temp_dir_content(i))), ...
%         'Position',[1 0.8 9.5 5.75]);    
%     else
%     exportToPPTX('addslide');    
%     exportToPPTX( ...
%         'addpicture', ...
%         fullfile(temp_img_dir, char(temp_dir_content(i))), ...
%         'Position',[1 0.8 9.5 5.75]);
%     end
% end

%% Guardar y cerrar
exportToPPTX( ...
    'save', ...
    ['Models\',MODEL.rootpath,'\results\reports\Results_',MODEL.rootpath,'.pptx']);
exportToPPTX('close');        
        
