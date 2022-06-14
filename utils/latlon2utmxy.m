%[x,y,zn,hm]=latlon2utmxy(zn, lat, lon)
%converts lat-lon coordinates to Universal Transverse Mercator projection
%zn, hm: UTM-zone and hemisphere (UTM-zone is computed from lon if zn<=0)
%lat, lon: geographical coordinates (degrees)
%x,y: UTM coordinates
function [x,y,zn,hm]=latlon2utmxy(zn, lat, lon)
 UTMScaleFactor = 0.9996;
 sm_a=6378137.0;
 sm_b=6356752.314;
 sm_EccSquared=6.69437999013e-03;

 ep2=(sm_a+sm_b)*(sm_a-sm_b)/sm_b^2;
 %ep2=sm_a^2/sm_b^2-1;
 
 %estimate best UTM-zone
 if(length(zn)==1 && zn<=0)
  zn=round((mean(lon(~isnan(lon)))+183)/6);
 else
  zn(zn<=0)=round((lon(zn<=0)+183)/6);
  zn(isnan(zn))=round((mean(lon(~isnan(lon)))+183)/6);
 endif
 
 lon0=6*zn-183;
 cphi=cos((pi/180)*lat);
 nu2=ep2*power(cphi,2);
 N=(sm_a^2/sm_b)*power(1 + nu2, -0.5);
 t=tan((pi/180)*lat);
 t2=power(t,2);

 l=(pi/180)*(lon-lon0);

 l3coef = 1 - t2 + nu2;
 l4coef = 5 - t2 + (9+4*nu2).*nu2;
 l5coef = 5 + (t2-18).*t2 + (14-58*t2).*nu2;
 l6coef = 61 + (t2-58).*t2 + (270-330*t2).*nu2;
 l7coef = 61 + ((179-t2).*t2-479).*t2;
 l8coef = 1385 + ((543-t2).*t2-3111).*t2;

 % Calculate easting (x)
 x = N.*cphi.*l.*( 1 +
       + (1/6)*cphi.*cphi .* l.*l.*( l3coef +
       + (1/20)*cphi.*cphi .* l.*l.*( l5coef +
       + (1/42)*cphi.*cphi .* l.*l.*l7coef ) ) );


 % Calculate northing (y)
 y = ArcLengthOfMeridian((pi/180)*lat)+0.5*t.*N.*(cphi.*cphi.*l.*l).*( 1 +
         (1/12)*cphi.*cphi.*l.*l.*( l4coef +
         (1/30)*cphi.*cphi.*l.*l.*( l6coef +
         (1/56)*cphi.*cphi.*l.*l.*l8coef ) ) );

 x=x*UTMScaleFactor+500000;
 y=y*UTMScaleFactor;
  y=y+10000000*(y<0);

 %
 hm=(lat>0)-(lat<=0);

endfunction



function result=ArcLengthOfMeridian(phi)
 sm_a=6378137.0;
 sm_b=6356752.314;

 n=(sm_a-sm_b)/(sm_a+sm_b);
 alpha=0.5*(sm_a+sm_b)*(1.0 + 0.25*n^2 * (1+0.0625*n^2));
 beta = (-1.5*n)*(1+(0.375*n^2)*(1-n^2/6));
 gamma = 0.9375*n^2*(1-0.5*n^2);
 delta = -(35/48)*n^3*(1+0.5625*n^2);
 epsilon = (315/512)*n^4;
 
 result=alpha*(phi
        + beta*sin(2*phi)
        + gamma*sin(4*phi)
        + delta*sin(6*phi)
        + epsilon*sin(8*phi));

endfunction
