%
% External variables: sId, x_min, x_max, y_min, y_max
%

%Open XYZ (lon,lat,depth) data in CSV
f=fopen(['data/' sId '.csv' ], "r");
fh=fgetl(f);
a=dlmread(f, ",");
fclose(f);

%Convert lat,lon to UTM (get utmZ and utmH)
y=a(:,getHeaderNos(fh,"y"));
x=a(:,getHeaderNos(fh,"x"));
if( ~any( abs(y)>180 | abs(x)>180 ) )
  [x,y, utmZ, utmH]=latlon2utmxy(-1, y, x );
  utmH=median(utmH);
  % ... again with the proper UTM zone
  [x,y]=latlon2utmxy(utmZ, a(:,getHeaderNos(fh,"y")), a(:,getHeaderNos(fh,"x")) );
endif
z=a( :, getHeaderNos(fh,"z") );
if( getHeaderNos(fh,"w") )
  w=a( :, getHeaderNos(fh,"w") );
else
  w=ones(size(z));
endif

clear a

%Remove NaNs
m=~( isnan(x) | isnan(y) | isnan(z) | isnan(w) );
x=x(m);
y=y(m);
z=z(m);
w=w(m);

clear m

%Interpolation parameters
mrTile0=0; %full area
mrSmooth=0; %no smoothing
#dxy=90; %m
#x_min=min(x); x_max=max(x);
#y_min=min(y); y_max=max(y);

%Export function
function writeRaster(r, fbn, xmn,xmx,ymn,ymx,uZ,uH)
  img2esriascii([fbn '.asc'], r, xmn, xmx, ymn, ymx);

  if( uH < 0 )
    uHs='south';
  else
    uHs='north';
  endif
  f=fopen([fbn '.crs'], "w");
  fprintf(f, "+proj=utm +zone=%d +%s +datum=WGS84 +units=m +no_defs\n", uZ, uHs);
  fclose(f);

  #disp( ['gdal_translate -ot Float32 -a_srs "$(cat ' [fbn '.crs'] ')" ' [fbn '.asc'] ' ' [fbn '.tif'] ] )
  #system( ['gdal_translate -ot Float32 -a_srs "$(cat ' [fbn '.crs'] ')" ' [fbn '.asc'] ' ' [fbn '.tif'] ] )
  #delete( [fbn '.asc'] )

endfunction


%Full data interpolation
disp(['Simple (weighted) multiresolution interpolation: ' sId])
[f_r, x_min,x_max, y_min,y_max]=mrinterpW(x,y,single(z),w, dxy, [mrTile0,mrSmooth], x_min,x_max,y_min,y_max);
writeRaster(f_r, ['data/' sId '-fr'], x_min,x_max,y_min,y_max, utmZ, utmH);

%disp(['x_min=', num2str(x_min), ', x_max=', num2str(x_max), ', y_min=', num2str(y_min), ', y_max=', num2str(y_max)])

clear f_r w

%%%%%%%%%%%%%%%%%%%%%
if( true )
%Full data fractal interpolation
disp(['Fractal multiresolution interpolation: ' sId])
%[f_fr, ~,~,~,~, ~, f_fH]=mrinterpfract(x,y,single(z), dxy, [mrTile0,mrSmooth], x_min,x_max,y_min,y_max);
f_fr=mrinterpfract(x,y,single(z), dxy, [mrTile0,mrSmooth], x_min,x_max,y_min,y_max);
writeRaster(f_fr, ['data/' sId '-ffr'], x_min,x_max,y_min,y_max, utmZ, utmH);
%writeRaster(f_fH, ['data/' sId '-ffH'], x_min,x_max,y_min,y_max, utmZ, utmH);

%disp(['x_min=', num2str(x_min), ', x_max=', num2str(x_max), ', y_min=', num2str(y_min), ', y_max=', num2str(y_max)])

clear f_fr
endif %if( true/false )

