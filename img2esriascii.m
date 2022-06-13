%fn: ESRI-ASCII filename
%I: image 2-D matrix (single band)
%lat_min,lat_max,lon_min,lon_max: image boundaries
function img2esriascii(fn, I, x_min,x_max,y_min,y_max)

 %image info
 height=size(I,1);
 width=size(I,2);

 %header file
 f=fopen(fn, 'wt');

 %cell size
 c_size=(x_max-x_min)/width;

 %metadata
 fprintf(f, "NCOLS %d\n", width);
 fprintf(f, "NROWS %d\n", height);
 fprintf(f, "XLLCORNER %0.6f\n", x_min);
 fprintf(f, "YLLCORNER %0.6f\n", y_min);
 fprintf(f, "CELLSIZE %0.6f\n", c_size);
 fprintf(f, "NODATA_VALUE %d\n", -9999);

 switch( class(I) )
  case 'logical'
   I=I; %do nothing
  otherwise
   I(isnan(I(:)))=-9999;
 endswitch
 
 %image data
 for n=1:height
  fprintf(f, " %g", I(n,:));
  fprintf(f, "\n");
 endfor
 
 fclose(f);
 
endfunction
