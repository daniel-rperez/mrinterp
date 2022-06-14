%funcion auxiliar para extraer los numeros de las columnas de interes en un CSV
% c=getHeaderNos(h,...)
% c: vector con los números de las columnas seleccionadas
% h: cadena con los nombres de las columnas separados por comas
% ...: cadenas con los nombres de las columnas que se quieren seleccionar
function c=getHeaderNos(h,varargin)
	if( ischar(h) )
		h=strsplit(h, ',');
	endif

	c=[];
	for n=1:length(varargin)
		f=find(strcmp(h,varargin{n}));
    if( any(f) )
      c(n)=f(1);
    endif
	endfor
endfunction
