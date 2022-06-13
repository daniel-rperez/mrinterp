% Multiresolution interpolation function
% [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterpR(x,y,z, dxy, slev, nR, xmin,xmax,ymin,ymax)
% x, y, z: coordinates and value
% dxy: spatial resolution of the raster (same for X and Y)
% slev: first interpolation level (optional) and smoothing level (nº de niveles de resolución en el que para la interpolación)
% nR: number of recursion
% mapa: raster matrix
% xmin, xmax, ymin, ymax: extent
% mapa0: raster with the projection of the original data points
function [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterpR(x,y,z, dxy, slev, nR, xmin,xmax,ymin,ymax)

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

  mapa0=nan(size(map0S));
    mapa0(map0N>0) = map0S(map0N>0) ./ map0N(map0N>0);

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
    if( l <= Nl - slev0 - slev1 ) %smoothing level
      if( b > 1 )
        for i=1:ni
          for j=1:nj
            maplN(i,j)=sum2( map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
            maplS(i,j)=sum2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
          endfor
        endfor
%        % Mucho más lento
%        mK=ones(b,b);
%        maplN=conv2(map0N,mK,'same');
%          maplN=maplN(1:b:end,1:b:end);
%        maplS=conv2(map0S,mK,'same');
%          maplS=maplS(1:b:end,1:b:end);
      else
        maplN=1.0*map0N;
        maplS=map0S;
      endif
    endif

    if( l>0 )
      %fill empty holes

      m0=(maplN>0);
      mapaN(m0)=maplN(m0);
      mapaS(m0)=maplS(m0);

      %compute two-resolution estimation
      mK=[0.614 1 0.614; 1 0 1; 0.614 1 0.614];
        mK=mK/sum(mK(:));

      %recursively iterate to obtain smoother interpolation
      for nr=0:nR

        %compute low resolution mean height
        mapaz=zeros(size(mapaS));
        m=(mapaN>0);
        mapaz(m)=mapaS(m)./mapaN(m);

        %compute new values only where no data exists at current resolution
        m=~m0;

        cmapaM=conv2((mapaN>0),mK, 'same');
        cmapaN=conv2(mapaN,mK, 'same');
        %cmapaS=conv2(mapaS,mK, 'same');
        cmapaZ=conv2(mapaz,mK, 'same'); %!testing here

        %
        mm=m & (cmapaM>0);

        mapaN(mm)=cmapaN(mm) ./ cmapaM(mm);
        mapaS(mm)=(cmapaZ(mm) .* cmapaN(mm)) ./ power(cmapaM(mm),2);

      endfor

      maplN(mm)=mapaN(mm);
      maplS(mm)=mapaS(mm);

      clear m mm mmm cmapaM cmapaN cmapaS cmapaZ mapaz

    endif

    %next resolution
    if( l<Nl-slev0)
      b=b/2; %!disp( [ 'b=' num2str(b) ] );
      mapaS=dupl2(maplS)/4;
      mapaN=dupl2(maplN)/4;
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
function s=sum2(m)
  s=sum(m(:));
endfunction

%svar2
function s=svar2(m,n)
  if( sum(n(:)>0)>1 )
    s=sum(power(m(n>0),2))-(sum(m(n>0)))^2/sum(n(:)>0);
  else
    s=inf;
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
