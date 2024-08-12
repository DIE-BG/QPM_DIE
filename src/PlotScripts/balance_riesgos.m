function MODEL = balance_riesgos(MODEL, MODEL_ANT, fan_params, varargin)

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
 for j = 1:length(fan_params.AnnoRange)
     if med_var{fan_params.AnnoRange(j)}.data < 3
         % Si Pronóstico es menor a 3
         
         % Probabilidad de que sea menor a 3
         p_less_3 = normcdf(3,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_der)-0.5;
         
         % Probabilidad de estar por encima de 5
         p_above_5 = 1 - normcdf(5,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_der);
         
         % Probabilidad de estar entre 3 y 5
         p_3_5 = 0.5 - p_less_3 - p_above_5;
         
     elseif med_var{fan_params.AnnoRange(j)}.data > 3 && med_var{fan_params.AnnoRange(j)}.data < 5
         % Pronóstico entre 3 y 5
         
         % Probabilidad de que sea menor a 3
         p_less_3 = normcdf(3,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_izq);
         % Probabilidad de estar entre 3 y 5
         p_3_med = 0.5 - p_less_3;
         p_med_5 = normcdf(5,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_der)-0.5;
         p_3_5 = p_3_med + p_med_5;
         
         % Probabilidad de estar por encima de 5
         p_above_5 = 1 - normcdf(5,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_der);
         
         
     elseif med_var{fan_params.AnnoRange(j)}.data > 5
         % Si el pronóstico es Mayor a 5
         
         % Probabilidad de que sea menor a 3
         p_less_3 = normcdf(3,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_izq);
         
         % Probabilidad de estar entre 3 y 5
         p_3_5 = normcdf(5,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_izq) - p_less_3;
         
         % Probabilidad de estar por encima de 5
         p_above_5 = 0.5 + 0.5 - normcdf(5,...
             med_var{fan_params.AnnoRange(j)}.data,...
             std_dev{fan_params.AnnoRange(j)}.data*fact_izq);
     end
     if j == 1
         % Almacenamiento de resultados
         if  isempty(fan_params.Esc)
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3 = p_less_3;
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5 = p_3_5;
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5 = p_above_5;
             % Diferencia con Mes previo
             % Menor a 3
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_less_3 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3;
             % Entre 3 y 5
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_3_5 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5;
             % Arriba de 5
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_above_5 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5;
             
         else
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3 = p_less_3;
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5 = p_3_5;
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5 = p_above_5;
             % Diferencia con Mes previo
             % Menor a 3
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_less_3 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_less_3;
             % Entre 3 y 5
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_3_5 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_3_5;
             % Arriba de 5
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y1.p_above_5 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y1.p_above_5;
         end
         
     elseif j == 2
         if  isempty(fan_params.Esc)
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3 = p_less_3;
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5 = p_3_5;
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5 = p_above_5;
             % Diferencia con Mes previo
             % Menor a 3
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_less_3 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3;
             % Entre 3 y 5
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_3_5 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5;
             % Arriba de 5
             MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_above_5 = ...
                 MODEL.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5 - ...
                 MODEL_ANT.Esc.v0.Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5;
             
         else
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3 = p_less_3;
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5 = p_3_5;
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5 = p_above_5;
             % Diferencia con Mes previo
             % Menor a 3
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_less_3 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_less_3;
             % Entre 3 y 5
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_3_5 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_3_5;
             % Arriba de 5
             MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_diff.y2.p_above_5 = ...
                 MODEL.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5 - ...
                 MODEL_ANT.Esc.(fan_params.Esc).Fanchart.(fan_params.PlotList{res.iter}).bal_riesgos.y2.p_above_5;
         end
     end
 end
%% Gráfico con tabla

if  isempty(fan_params.Esc)
    st_br = MODEL.Esc.v0.Fanchart.D4L_CPI;
    st_br_ant = MODEL_ANT.Esc.v0.Fanchart.D4L_CPI;
else
    st_br = MODEL.Esc.(fan_params.Esc).Fanchart.D4L_CPI;
    st_br_ant = MODEL_ANT.Esc.(fan_params.Esc).Fanchart.D4L_CPI;
end


figure('Position', [1 42.0182 1117.1 776.73]);
x_br = [1 1 10 10];
y_br = [0 9 9 0];
fill(x_br,y_br,'w','LineStyle','none');
set(gca,'XTick',[],'YTick',[])

text(5.5,8.5,'Balance de Riesgos de la Tasa de Inflación Interanual','HorizontalAlignment','Center','FontSize',14,'FontWeight','bold');
text(5.5,8.1,['para ',dat2char(fan_params.AnnoRange(1)),' y ',dat2char(fan_params.AnnoRange(2))],'HorizontalAlignment','Center','FontSize',14,'FontWeight','bold');
hold on;

text(6.35,7.0,'Pronóstico de','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

text(5.5,6.6,'Corrimiento','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(7.2,6.6,'Corrimiento','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

text(5.5,6.2,MODEL.CORR_DATE_ANT,'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(7.2,6.2,MODEL.CORR_DATE,'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(8.7,5.8,'Diferencia','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

text(3.0,5.8,'Pronóstico para ','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(3.0,5.4,dat2char(fan_params.AnnoRange(1)),'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

text(5.5,5.4,'(A)','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(7.2,5.4,'(B)','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(8.7,5.4,'(C) = (B) - (A)','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

line([1.30 9.7],[5.15 5.15],'LineStyle','-','LineWidth',0.5,'Color','k')

text(3.00,4.8,['P(\pi \in [3.0%, 5.0%])'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,4.8,[num2str(round(st_br_ant.bal_riesgos.y1.p_3_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,4.8,[num2str(round(st_br.bal_riesgos.y1.p_3_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,4.8,[num2str(round(st_br.bal_riesgos.y1.p_3_5,3) - round(st_br_ant.bal_riesgos.y1.p_3_5,3),'%4.3f'),'%'],...
               'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');

text(3.00,4.4,['P(\pi > 5%)'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,4.4,[num2str(round(st_br_ant.bal_riesgos.y1.p_above_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,4.4,[num2str(round(st_br.bal_riesgos.y1.p_above_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,4.4,[num2str(round(st_br.bal_riesgos.y1.p_above_5,3) - round(st_br_ant.bal_riesgos.y1.p_above_5,3),'%4.3f'),'%'],...
               'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');

text(3.00,4.0,['P(\pi < 3.0%)'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,4.0,[num2str(1- round(st_br_ant.bal_riesgos.y1.p_3_5,3)-round(st_br_ant.bal_riesgos.y1.p_above_5,3),'%4.3f'),'%'],...
    'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,4.0,[num2str(1- round(st_br.bal_riesgos.y1.p_3_5,3)-round(st_br.bal_riesgos.y1.p_above_5,3),'%4.3f'),'%'],...
    'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,4.0,[num2str(st_br.bal_diff.y1.p_less_3,'%4.3f'),'%'],...
    'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');

text(3.0,2.8,'Pronóstico para','HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(3.0,2.4,dat2char(fan_params.AnnoRange(2)),'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');

line([1.30 9.7],[2.15 2.15],'LineStyle','-','LineWidth',0.5,'Color','k')

text(3.00,1.8,['P(\pi \in [3.0%, 5.0%])'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,1.8,[num2str(round(st_br_ant.bal_riesgos.y2.p_3_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,1.8,[num2str(round(st_br.bal_riesgos.y2.p_3_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,1.8,[num2str(round(st_br.bal_riesgos.y2.p_3_5,3) - round(st_br_ant.bal_riesgos.y2.p_3_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');

text(3.00,1.4,['P(\pi > 5%)'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,1.4,[num2str(round(st_br_ant.bal_riesgos.y2.p_above_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,1.4,[num2str(round(st_br.bal_riesgos.y2.p_above_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,1.4,[num2str(round(st_br.bal_riesgos.y2.p_above_5,3) - round(st_br_ant.bal_riesgos.y2.p_above_5,3),'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');

text(3.00,1.0,['P(\pi < 3.0%)'],'HorizontalAlignment','Center','FontSize',12,'FontWeight','bold','LineStyle','none');
text(6.00,1.0,[num2str(1- round(st_br_ant.bal_riesgos.y2.p_3_5,3)-round(st_br_ant.bal_riesgos.y2.p_above_5,3),'%4.3f'),'%'],...
    'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(7.60,1.0,[num2str(1- round(st_br.bal_riesgos.y2.p_3_5,3)-round(st_br.bal_riesgos.y2.p_above_5,3),'%4.3f'),'%'],...
    'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
text(9.20,1.0,[num2str(st_br.bal_diff.y2.p_less_3,'%4.3f'),'%'],'HorizontalAlignment','Right','FontSize',12,'LineStyle','none');
      


SimTools.scripts.pausaGuarda(...
            fullfile(fan_params.SavePath, ...
            sprintf("Bal_Riesgos %s.png", fan_params.PlotList{res.iter})), ...
            'AutoSave', fan_params.AutoSave ...
            );


 
end