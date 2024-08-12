function h = colorarea(rangex,rangey,color,varargin)

%COLORAREA pone color a un área delimitada por los vertices
%          x(1),y(1)...x(n),y(n)
%colorarea(rangex,rangey,color)
%rangex = vector que contiene los valores en x de los vértices
%         de la figura, en orden ascendente
%rangey = vector que contiene los valores en y de los vértices 
%         de la figura, en orden ascendente
%color = vector columna que contiene los valores del color deseado,
%        estos valaores deben ser >=0 y <=1. Los colores "puros" eb 
%        MATLAB son:
%        [1 1 0] y yellow
%        [1 0 1] m magenta 
%        [0 1 1] c cyan
%        [1 0 0] r red
%        [0 1 0] g green
%        [0 0 1] b blue
%        [1 1 1] w white
%        [0 0 0] k black
%        Si se desea incluir un color "puro" se puede incluir como 
%        argumento de colorarea el nombre corto entre comillas simples
%

GRADO = 1;

% -----function GRAYAREA body----- %

if isempty(rangex) || isempty(rangey)
  return
end

endDate = max(rangex);
startDate = min(rangex);

status = get(gca,'nextplot');
set(gca,'nextplot','add');
y = rangey;
%Muestra la grilla sobre el área sombreada
set(gca,'layer','top');
% Crea el área sombreada
h = fill(dat2dec([startDate*[1,1],endDate*[1,1]])+0.13,[y(1),y(2),y(2),y(1)],GRADO*color);
set(h,'linestyle','none');
%Envía el área sombreada al fondo de la gráfica
ch = get(gca,'children');  
set(gca,'children',[ch(2:end);ch(1)]);
% Regresa a las propiedades axis a como estaban antes de colocar el área
% sombreada.
set(gca,'nextplot',status);

end
