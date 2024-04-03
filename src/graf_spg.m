%% Retrieving Historical Data
% Variables a graficar
err_list = {'L_GDP_RW_GAP','D4L_CPI_RW','RS_RW',...
            'D4L_CPI_NOSUBY','D4L_GDP', 'D4L_CPIXFE',...
            'D4L_S','D4L_MB','RS',...
            'D4L_CPI', 'D4L_Z','D4L_VEL', 'RR','D4_GDP_SM'};
% Datos históricos
HistData = databank.fromCSV(fullfile('output','kalm_his.csv'),'Select',err_list);
% Rango de evaluación
rng = qq(2005,1):qq(2023,1);
%% Retrieving each Scenario data

ResPath = 'fulldata';

% Data retrieving and creation of structures for each variable 
GenStruct = struct();
SimStruc = struct();
% Lista de archivos dentro de carpeta ResPath
FullDataList = ls(fullfile(ResPath));
FullDataList = cellstr(FullDataList(3:end,:));
    
for iter = 1:length(FullDataList)
        SimStruc.(strcat('Sim',num2str(iter))) = databank.fromCSV(...
            fullfile(...
            ResPath, ...
            FullDataList{iter}),...
                'Select', err_list...
            );
        % Creación de estructuras por cada variable (Aqui deben agregarse
        % si fuera necesario según la cantidad de variables)
           
        for var = 1:length(err_list)
        GenStruct.(err_list{var}).((strcat('Sim',num2str(iter)))) = SimStruc.(strcat('Sim',num2str(iter))).(err_list{var});
        end

end

% Crop the Data structure to match the data used for evaluation purposes
% (if the data_corr.csv has longer history backwards).
GenStruct = databank.clip(GenStruct, rng(1), rng(end));


%% Graphs
% Rangos para ajustar si es necesario
% y_l = {[-10 15], [-30 20], [-1 6], [-5 4], [-10 15], [-10 15], ...
%        [-6 12], [-10 25], [-1 9], [-10 15], [-25 15], [-3 -2],...
%        [-10 15],[-8 10]};
   

for ig = 1:length(err_list)
        figure('Position', [1 42.0181818181818 1675.63636363636 825.6]);
       
        h = plot([struct2array(eval(strcat('GenStruct.',err_list{ig}))), HistData.(err_list{ig}){rng}],...
             'Color', [.7 .7 .7]);
              set(h(end),'LineStyle','-','LineWidth',2,...
                         'Marker','o','MarkerFaceColor','k','MarkerSize',3, ...
                         'Color','k');
%         ylim([y_l{ig}]);    
        title(sprintf("%s \n Comparación Entre Historia y Pronóstico", err_list{ig}), ...
          'Interpreter', 'none', ...
          'FontSize', 15);
          
      
        saveas(gcf, ...
            fullfile(...
            'graphs',...
            sprintf( "%s.png", err_list{ig}))...
            );
close all
end

