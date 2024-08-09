function MODEL = balance_riesgos(MODEL, fan_params, varargin)

p = inputParser;
    addParameter(p, 'iter', 1);
    addParameter(p, 'FMSE', {});
    addParameter(p, 'sesgo', []);
parse(p, varargin{:});
res = p.Results;

 %%
 
  % Factor para cola izquierda/derecha
 fact_izq = res.sesgo(1);
 fact_der = res.sesgo(end);
 
 med_var = MODEL.F_pred.(fan_params.PlotList{res.iter}){MODEL.DATES.Fanchart};
 std_dev = tseries(MODEL.DATES.Fanchart,res.FMSE.data); 
 

 
 % Fin año 1
 
 if med_var{fan_params.AnnoRange(1)}.data < 3 
     % Si Pronóstico es menor a 3
     
     % Probabilidad de que sea menor a 3
     p_med_3 = normcdf(3,...
         med_var{fan_params.AnnoRange(1)}.data,...
         std_dev{fan_params.AnnoRange(1)}.data*fact_der)-0.5;
     p_less_3 =  0.5 + p_med_3;
     
     % Probabilidad de estar por encima de 5
     p_above_5 = 1 - normcdf(5,...
         med_var{fan_params.AnnoRange(1)}.data,...
         std_dev{fan_params.AnnoRange(1)}.data*fact_der);
     
     % Probabilidad de estar entre 3 y 5
     p_3_5 = 0.5 - p_med_3 - p_above_5;
     
 elseif med_var{fan_params.AnnoRange(1)}.data > 3 && med_var{fan_params.AnnoRange(1)}.data < 5
     % Pronóstico entre 3 y 5
     
     % Probabilidad de que sea menor a 3
     p_less_3 = normcdf(3,...
                        med_var{fan_params.AnnoRange(1)}.data,...
                        std_dev{fan_params.AnnoRange(1)}.data*fact_izq);
    % Probabilidad de estar entre 3 y 5
     p_3_med = 0.5 - p_less_3;   
     p_med_5 = normcdf(5,...
                      med_var{fan_params.AnnoRange(1)}.data,...
                      std_dev{fan_params.AnnoRange(1)}.data*fact_der)-0.5;
    p_3_5 = p_3_med + p_med_5;
    
    % Probabilidad de estar por encima de 5
    p_above_5 = 1 - normcdf(5,...
                      med_var{fan_params.AnnoRange(1)}.data,...
                      std_dev{fan_params.AnnoRange(1)}.data*fact_der);
     
     
 elseif med_var{fan_params.AnnoRange(1)}.data > 5 
     % Si el pronóstico es Mayor a 5
     
     % Probabilidad de que sea menor a 3
     p_less_3 = normcdf(3,...
                      med_var{fan_params.AnnoRange(1)}.data,...
                      std_dev{fan_params.AnnoRange(1)}.data*fact_izq);
     
     % Probabilidad de estar entre 3 y 5             
     p_3_5 = normcdf(5,...
                      med_var{fan_params.AnnoRange(1)}.data,...
                      std_dev{fan_params.AnnoRange(1)}.data*fact_izq) - p_less_3;
     
     % Probabilidad de estar por encima de 5
     p_above_5 = 0.5 + 0.5 - normcdf(5,...
                      med_var{fan_params.AnnoRange(1)}.data,...
                      std_dev{fan_params.AnnoRange(1)}.data*fact_izq);  
 end
    
    if  isempty(fan_params.Esc)
        MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.p_less_3 = p_less_3;
        MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.p_3_5 = p_3_5;
        MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.p_above_5 = p_above_5;
    else
        MODEL.Esc.(fan_params.Esc).(fan_params.PlotList{res.iter}).bal_riesgos.p_less_3 = p_less_3;
        MODEL.Esc.(fan_params.Esc).(fan_params.PlotList{res.iter}).bal_riesgos.p_3_5 = p_3_5;
        MODEL.Esc.(fan_params.Esc).(fan_params.PlotList{res.iter}).bal_riesgos.p_above_5 = p_above_5;
    end
    
end