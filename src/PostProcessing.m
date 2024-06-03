function MODEL = PostProcessing(MODEL, varargin)
%{
PostProcessing realiza el post procesamiento de los datos generados en los
ejercicios de simulación del QPM. 

Ejecuta para las variables logaritmicas en 'list':
    - Ajuste por estacionalidad (X12)
    - Tendencia HP de las variables ajustadas por estacionalidad.
    - Brecha relacionada a las dos anteriores

Para las variables logaritmicas en list_niv:
    - Recuperación del nivel (exp(x/100))
    - Ajuste por estacionalidad (x12)
    - Tendencia HP de las variables ajustadas por estacionalidad.
{
## Syntax ##

    MODEL = PpostProcessing(MODEL, varargin)

## Input Arguments ##

__`MODEL`__ [ struct ] -
Debe contener al menos la estructura con los resultados del proceso de
simulación MODEL.F_pred.

* 'list' = {} [ `Cell` ] - Nombres de las variables logaritmicas para
procesar.

* 'list_niv' = {} [ `Cell` ] - Nombres de las variables logaritmicas para
reconstruir el nivel original y luego procesar.

- DIE
- Marzo 2024
- MJGM/JGOR
%}
p = inputParser;
    addParameter(p, 'list', {});
    addParameter(p, 'list_niv', {});
    addParameter(p, 'Esc', {MODEL.CORR_VER, MODEL.F_pred});
parse(p, varargin{:});
params = p.Results; 


MODEL.PostProc.(params.Esc{1}) = struct();

%% Datos
if strcmp(params.Esc{1}, 'v0')
   pred = MODEL.F_pred;
   
else
   pred = params.Esc{2};
    
end
    
%% Logaritmos

temp_db = databank.copy(params.Esc{2}, params.list); 

% Se agregan los logs desestacionalizados del modelo a la estructura de
% post-procesamiento (para facilitar graficado)
list_gaps_mod = get(MODEL.MF, 'xlist');
list_gaps_mod = params.list + list_gaps_mod;

for i = 1:length(list_gaps_mod)
    if startsWith(list_gaps_mod{i}, 'L_') && ~endsWith(list_gaps_mod{i}, '_BAR') && ~endsWith(list_gaps_mod{i}, '_GAP')
       MODEL.PostProc.(params.Esc{1}).l_sa.(strcat(list_gaps_mod{i},'_SA')) = pred.(list_gaps_mod{i});
    end
end  

% Variables que es necesario desestacionalizar (porque la variable endogena
% del modelo es la D4L)
if isanystri('L_MB',params.list) && isanystri('L_VEL',params.list)
    MODEL.PostProc.(params.Esc{1}).l_sa.L_MB_SA = x12(temp_db.L_MB);
    MODEL.PostProc.(params.Esc{1}).l_sa.L_VEL_SA = x12(temp_db.L_VEL);
end
                           
% Tendencias HP de los logs desestacionalizados           
MODEL.PostProc.(params.Esc{1}).l_bar = databank.apply(@(x) hpf(x),...
                                                     MODEL.PostProc.(params.Esc{1}).l_sa,...
                                                    'EndsWith=','_SA',...
                                                    'RemoveEnd=',true,...
                                                    'Append=', '_BAR',...
                                                    'RemoveSource=',true);
                                                
% Se agregan las tendencias del modelo a la estructura de
% post-procesamiento (para facilitar graficado)
list_gaps_mod = get(MODEL.MF, 'xlist');
% cosa = struct;
for i = 1:length(list_gaps_mod)
    if startsWith(list_gaps_mod{i}, 'L_') && endsWith(list_gaps_mod{i}, '_BAR')
       MODEL.PostProc.(params.Esc{1}).l_bar.(list_gaps_mod{i}) = pred.(list_gaps_mod{i});
    end
end  
 
% Se agregan las brechas del modelo a la estructura de post-procesamiento
% (para facilitar el graficado)
list_gaps_mod = get(MODEL.MF, 'xlist');
% cosa = struct;
for i = 1:length(list_gaps_mod)
    if endsWith(list_gaps_mod{i}, '_GAP')
       MODEL.PostProc.(params.Esc{1}).l_gap.(list_gaps_mod{i}) = pred.(list_gaps_mod{i});
    end
end  

% brechas relacionadas
for i = 1:length(params.list)
   MODEL.PostProc.(params.Esc{1}).l_gap.(strcat(params.list{i},'_GAP')) = ...
                        MODEL.PostProc.(params.Esc{1}).l_sa.(strcat(params.list{i},'_SA')) - ...
                        MODEL.PostProc.(params.Esc{1}).l_bar.(strcat(params.list{i},'_BAR'));   
end



%% Niveles
temp_db = databank.copy(params.Esc{2}, params.list_niv); 


if isanystri('L_MB',params.list) 
    temp_db.L_MB = x12(temp_db.L_MB);
end

% Niveles 
MODEL.PostProc.(params.Esc{1}).niv = databank.apply(@(x) exp(x/100),...
                    temp_db,...
                    'StartsWith=','L_',...
                    'RemoveStart=', true,...
                    'RemoveSource=',true);
 
MODEL.PostProc.(params.Esc{1}).niv_bar = databank.apply(@(x) hpf(x),...
                                                        MODEL.PostProc.(params.Esc{1}).niv,...
                                                        'Append=', '_BAR',...   
                                                        'RemoveSource=',true);                                                  
end