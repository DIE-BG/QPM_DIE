%% Transformaciones a variables trimestrales
% Transformaciones relaconadas a GDP_RW
MODEL.PreProc.quarterly.L_GDP_RW_SM = MODEL.PreProc.quarterly.GDP_RW{(MODEL.PreProc.quarterly.GDP_RW.Start+3):MODEL.PreProc.quarterly.GDP_RW.End};
s = movsum(MODEL.PreProc.quarterly.GDP_RW.log.data, 4);
MODEL.PreProc.quarterly.L_GDP_RW_SM.Data = s(3:end-1);
MODEL.PreProc.quarterly.L_GDP_RW_SM.Comment = 'Suma Móvil de 4T del PIB de EEUU';
MODEL.PreProc.quarterly.L_GDP_RW_SM.Caption = 'Logaritmo';

MODEL.PreProc.quarterly.L_GDP_RW = MODEL.PreProc.quarterly.GDP_RW.log;
MODEL.PreProc.quarterly.L_GDP_RW.Comment = 'Producto Interno Bruto Real EEUU (Logaritmo)';

%Tasas de variación
% Interanual
MODEL.PreProc.quarterly.D4L_GDP_RW = MODEL.PreProc.quarterly.L_GDP_RW.diff(-4);
MODEL.PreProc.quarterly.D4L_GDP_RW.Comment = 'Tasa de Variación Interanual del PIB de EEUU';
% Interanual GDP_RW_SM
MODEL.PreProc.quarterly.D4L_GDP_RW_SM = MODEL.PreProc.quarterly.L_GDP_RW_SM.diff(-4);
MODEL.PreProc.quarterly.D4L_GDP_RW_SM.Comment = 'Tasa de Variación Interanual de la suma de 4 Trimestres del PIB de EEUU';