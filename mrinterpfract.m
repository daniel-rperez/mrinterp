% Multiresolution fractal (random) interpolation function
% [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterpfract(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)
% x, y, z: coordinates and value
% dxy: spatial resolution of the raster (same for X and Y)
% slev: first interpolation level (optional) and smoothing level (nº de niveles de resolución en el que para la interpolación)
% mapa: raster matrix
% xmin, xmax, ymin, ymax: extent
% mapa0: raster with the projection of the original data points
function [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterpfract(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)

  %compute map size
  if( ~exist('xmin') || ~exist('xmax') )
    xmin=min(x);
    xmax=max(x);
  endif
  dx=xmax-xmin;
  wa=ceil(dx/dxy);
  Nlx=nextpow2(wa);
  w=pow2(Nlx);
  xmin=xmin-(w-wa)*dxy/2;
  xmax=xmin+w*dxy;

  if( ~exist('ymin') || ~exist('ymax') )
    ymin=min(y);
    ymax=max(y);
  endif
  dy=ymax-ymin;
  ha=ceil(dy/dxy);
  Nly=nextpow2(ha);
  h=pow2(Nly);
  ymin=ymin-(h-ha)*dxy/2;
  ymax=ymin+h*dxy;

  m=(xmin<=x & x<=xmax & ymin<y & y<=ymax);
  x=x(m);
  y=y(m);
  z=z(m);
  clear m

  Nl=min([Nlx,Nly]);
  bl=min([w,h]);

  %raw projection
  i=1+floor((ymax-y)/dxy);
  j=1+floor((x-xmin)/dxy);
  ij=i+h*(j-1);

  map0N=zeros(h,w);
  map0S=zeros(h,w);
  for n=1:length(z)
    ijn=ij(n);
    map0N(ijn)=map0N(ijn)+1;
    map0S(ijn)=map0S(ijn)+z(n);
  endfor

  %prune scarcely populated pixels
  nmin=1;
  while( sum(map0N(:)>nmin) > 0.99*sum(map0N(:)>0) )
    nmin=nmin+1;
  endwhile
  map0S(map0N>=nmin)=map0S(map0N>=nmin)./map0N(map0N>=nmin);
  map0S(map0N<nmin)=0;
  map0N=(map0N>=nmin); %start afresh from equally weighted pixels (!)

  mapa0=nan(size(map0S));
    mapa0(map0N>0) = map0S(map0N>0) ./ map0N(map0N>0);

  %pixel roughness (vertical distance to neighbors)
  map0R=zeros(h,w);
    m=(map0N>0);
    mK=[0.5 1 0.5; 1 0 1; 0.5 1 0.5];
      mK=mK/sum(mK(:));
    mapa0n=conv2(map0N, mK, 'same');
     mm=m & (mapa0n > 0); %exclude isolated points
    mapa0s=conv2(map0S, mK, 'same');
      mapa0s(mm)=mapa0s(mm) ./ mapa0n(mm);
      mapa0s(~mm)=0;
  map0R(mm)=power(mapa0(mm) - mapa0s(mm), 2);
    mmR=median([0.002^2; map0R(mm)]);
    map0R(m & ~mm)=mmR; %estimate isolated points from global roughness
  %map0R(m)=median(map0R(m));

  %scaling levels
  if( length(slev) == 1 )
    slev0=0;
    slev1=slev;
  else
    slev0=slev(1);
    slev1=slev(2);
  endif

  %multiscale iteration
  b=bl/2^slev0;
  for l=0:Nl-slev0
    ni=h/b;
    nj=w/b;

    %data sums and number
    maplS=zeros(ni,nj);
    maplN=zeros(ni,nj);
    maplV=zeros(ni,nj);
    maplR=zeros(ni,nj);
    if( l <= Nl - slev0 - slev1 ) %smoothing level
      for i=1:ni
        for j=1:nj
          maplN(i,j)=sum2( map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
          maplS(i,j)=sum2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
          maplV(i,j)=nvar2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b), map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
          maplR(i,j)=sum2( map0R(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
        endfor
      endfor
    endif

    if( l>0 )

      %compute low resolution mean height
      mapaz=zeros(size(mapaS));
      m=(mapaN>0);
      mapaz(m)=mapaS(m)./mapaN(m);
        mK=[0.5 1 0.5; 1 1 1; 0.5 1 0.5];
          mK=mK/sum(mK(:));
        mapaw=conv2(power(mapaz - conv2(mapaz,mK, 'same'),2),mK, 'same');

      %pad low-resolution map
      m=(maplN>0);
      mapaN(m)=maplN(m);
      mapaS(m)=maplS(m);
      mapaR(m)=maplR(m);
      mapaV(m)=maplV(m) + mapaw(m); %complete variance
      mapaz(m)=maplS(m)./maplN(m);
      m=~m;

      %compute two-resolution estimation
      %mK=[0 1 0; 1 0 1; 0 1 0];
      mK=[0.5 1 0.5; 1 0 1; 0.5 1 0.5];
        mK=mK/sum(mK(:));
      cmapaM=conv2((mapaN>0),mK, 'same');
        %cmapaM=cmapaM(2:end-1,2:end-1);
      cmapaN=conv2(mapaN,mK, 'same');
        %cmapaN=cmapaN(2:end-1,2:end-1);
      cmapaS=conv2(mapaS,mK, 'same');
        %cmapaS=cmapaS(2:end-1,2:end-1);
      cmapaV=conv2(mapaV,mK, 'same');
        %cmapaV=cmapaV(2:end-1,2:end-1);
      cmapaR=conv2(mapaR,mK, 'same');
        %cmapaV=cmapaR(2:end-1,2:end-1);
      cmapaZ=conv2(mapaz,mK, 'same'); %interpolate heights instead of S?
        %cmapaZ=cmapaZ(2:end-1,2:end-1);

      %
      m=m & (cmapaM>0);

      %interpolated N
      maplN(m)=cmapaN(m) ./ cmapaM(m);

      %interpolated Z and simulated variance
      H=0.5/(Nl-slev0-l+1)*log2(cmapaV(m) ./ cmapaR(m));
      maplR(m)=cmapaR(m) ./ cmapaM(m);
      maplV(m)=pow2(-2*H) .* cmapaV(m) ./ cmapaM(m);

      %add simulated dZ
      cmapaZ(m)=cmapaZ(m) + 1*(rand(sum(m(:)),1)-0.5) .* sqrt(12*maplV(m)) .* cmapaM(m);

      %convert Z to S again
      mm=m & ~isnan(cmapaZ);
      maplS(mm)=(cmapaZ(mm) .* cmapaN(mm)) ./ power(cmapaM(mm),2); %!testing here

      mm=m & isnan(cmapaZ); %!hack in case some Z turns NaN
      maplS(mm)=cmapaS(mm) ./ cmapaM(mm);

    endif

    %next resolution
    if( l<Nl-slev0)
      b=b/2;
      mapaS=dupl2(maplS)/4;
      mapaN=dupl2(maplN)/4;
      mapaV=dupl2(maplV);
      mapaR=dupl2(maplR)/4;
    endif

  endfor

  mapa=nan(size(maplS));
    mapa(maplN>0) = maplS(maplN>0) ./ maplN(maplN>0);

  if( any(size(mapa0) ~= size(mapa)) )
    warning "mapa0 and mapa have different sizes"
  endif

endfunction

%utility functions
%sum2
function s=sum2(x)
  s=sum(x(:));
endfunction

%nvar2
function s=nvar2(x,n)
  m=(n(:)>0);
  if( sum(m)>1 )
    s=sum( n(m).*x(m) )/sum(n(m));
    s=sum( n(m).*power(x(m)-s,2) )/sum(n(m));
  else
    s=0;
  endif
endfunction

%duplicate
function b=dupl2(a)
  ha=size(a,1);
  wa=size(a,2);
  b=nan(2*size(a));
  b(1:2:2*ha,1:2:2*wa)=a;
  b(2:2:2*ha,1:2:2*wa)=a;
  b(1:2:2*ha,2:2:2*wa)=a;
  b(2:2:2*ha,2:2:2*wa)=a;
endfunction
