% Multiresolution fractal (random) interpolation function
% [mapa, xmin,xmax, ymin,ymax, mapa0, mapaH]=mrinterpfract(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)
% x, y, z: coordinates and value
% dxy: spatial resolution of the raster (same for X and Y)
% slev: first interpolation level (optional) and smoothing level (number of levels without exact interpolation)
% mapa: interpolated raster matrix
% xmin, xmax, ymin, ymax: extent
% mapa0: raster with the projection of the original data points
% mapaH: raster with the last calculated Hurst exponent
function [mapa, xmin,xmax, ymin,ymax, mapa0, mapaH]=mrinterpfract(x,y,z, dxy, slev, xmin,xmax,ymin,ymax, Hexp0)

  %data type
  dcls=class(z);

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
  
  %default value of Hurst exponent
  if( ~exist('Hexp0') )
    Hexp0=2.0;
  endif
  
  Nl=min([Nlx,Nly]);
  bl=min([w,h]);

  %raw projection
  i=1+floor((ymax-y)/dxy);
  j=1+floor((x-xmin)/dxy);
  ij=i+h*(j-1);
  
  %map number of points and sum of values
  map0N=zeros(h,w, dcls);
  map0S=zeros(h,w, dcls);
  for n=1:length(z)
    ijn=ij(n);
    map0N(ijn)=map0N(ijn)+1;
    map0S(ijn)=map0S(ijn)+z(n);
  endfor

  %OPTIONAL: prune scarcely populated pixels
  if( false )
    nmin=1;
    while( sum(map0N(:)>nmin) > 0.99*sum(map0N(:)>0) )
      nmin=nmin+1;
    endwhile
    map0S(map0N>=nmin)=map0S(map0N>=nmin)./map0N(map0N>=nmin);
    map0S(map0N<nmin)=0;
    map0N=(map0N>=nmin); %start afresh from equally weithted pixels (!)
  endif
  
  mapa0=nan(size(map0S), dcls);
    mapa0(map0N>0) = map0S(map0N>0) ./ map0N(map0N>0);
  
  clear x y z i j ij
  
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
    maplS=zeros(ni,nj, dcls);
    maplN=zeros(ni,nj, dcls);
    maplW=zeros(ni,nj, dcls);
    maplH=Hexp0*ones(ni,nj, dcls); %default globally estimated Hurst for larger scale
    if( l <= Nl - slev0 - slev1 ) %smoothing level
      if( b > 1 )
        for i=1:ni
          for j=1:nj
            maplN(i,j)=sum2( map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) ); %number of points inside cell
            maplS(i,j)=sum2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) ); %sum of heights inside cell
            maplW(i,j)=nvar2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b), map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) ); %variance of heights inside cell
          endfor
        endfor
      else
        maplN=cast(map0N,dcls);
        maplS=map0S;
        maplW=0.0*maplN;
      endif
    endif
    
    if( l>0 )
      %mK=[0 1 0; 1 0 1; 0 1 0];
      mK=[0.614 1 0.614; 1 0 1; 0.614 1 0.614];
        mK=mK/sum(mK(:));

      %compute low resolution mean height
      ma=(mapaN>0); %where data exists at lower resolution
      mapaz=zeros(size(mapaS), dcls);
      mapaz(ma)=mapaS(ma)./mapaN(ma);

      %pad low-resolution map with current resolution data
      m=(maplN>0); %pixels with data at current resolution
      mapaN(m)=maplN(m);
      %mapaS(m)=maplS(m);
      mapaz(m)=maplS(m)./maplN(m);

      %compute smoothed version of low resolution variables
      cmapaM=conv2(ma,mK, 'same');
      cmapaN=conv2(mapaN,mK, 'same');
      %cmapaS=conv2(mapaS,mK, 'same');
      cmapaZ=conv2(mapaz,mK, 'same');

      %compute smoothed versions of low resolution roughness and Hurst
      mma=isnan(mapaW);
      cmapaMM=conv2(~mma,mK, 'same');
        mapaW(mma)=0;
      cmapaW=conv2(mapaW,mK, 'same');
        mapaH(mma)=0;
      cmapaH=conv2(mapaH,mK, 'same');
      mma=(cmapaMM>0); %mask pixels where roughness and Hurst can be recovered
      
      %compute local roughness and Hurst from current resolution (exact or interpolated) and low-resolution roughness
      mm=mma & ~(maplW>0); %pixels with low-resolution roughness and no current resolution
      maplH(mm)=cmapaH(mm) ./ cmapaMM(mm); %smooth roughness and scaling where no current resolution data exist
      maplW(mm)=cmapaW(mm) ./ cmapaMM(mm) .* power(2, -2*maplH(mm));
      mm=mma & (maplW>0); %pixels with low-resolution roughness and current resolution
      maplH(mm)=-0.5*log2( maplW(mm) ./ (cmapaW(mm) ./ cmapaMM(mm)) ); %compute local Hurst: Wa = 2^2H W => H = (1/2) log2( Wa / W )
        maplH(mm & maplH<1 )=1;
        maplH(mm & maplH>3 )=3;
      maplH(~mma)=NaN; %where no data exist from lower resolutions, nothing can be said
      maplW(~mma)=NaN;
      
      %simulate fractal heights
      mm=~m & mma; %pixels with no data at current resolution but with roughness data at lower
      cmapaZ(mm)=cmapaZ(mm) + 2*(rand(sum(mm(:)),1)-0.5)/sqrt(3) .* sqrt(cmapaW(mm) ./ cmapaMM(mm)) .* cmapaM(mm);
      
      %fill in missing data at pixels with no data at current resolution
      m=~m & (cmapaM>0); %cmapaM should include cmapaMM
      maplN(m)=cmapaN(m) ./ cmapaM(m);
      maplS(m)=(cmapaZ(m) .* cmapaN(m)) ./ power(cmapaM(m),2);
      
      clear m mm ma mma
      clear cmapaM cmapaMM cmapaN cmapaZ cmapaW cmapaH %cmapaS
      clear mapaz mapaN mapaS mapaW mapaH
      
    endif
    
    %next resolution
    if( l<Nl-slev0)
      b=b/2;
      %disp( [ 'PROGRESS: b=' num2str(b) ] );
      mapaS=dupl2(maplS)/4;
      mapaN=dupl2(maplN)/4;
      mapaW=dupl2(maplW);
      mapaH=dupl2(maplH);
    endif

  endfor

  %return mapa
  mapa=nan(size(maplS), dcls);
    mapa(maplN>0) = maplS(maplN>0) ./ maplN(maplN>0);
  
  %return last Hurst exponent
  mapaH=maplW;
  
  %one last check
  if( any(size(mapa0) ~= size(mapa)) )
    warning "mapa0 and mapa have different sizes"
  endif
  
endfunction

%utility functions
%sum2
function s=sum2(m)
  s=sum(m(:));
endfunction

%nvar2
function s=nvar2(x,n)
  m=(n(:)>0);
  df=sum(m)-1;
  if( df > 0 )
    s=sum( n(m).*x(m) )/(1+df);
    s=sum( n(m).*power(x(m)-s,2) )/df;
  else
    s=0;
  endif
endfunction

%duplicate
function b=dupl2(a)
  ha=size(a,1);
  wa=size(a,2);
  b=nan(2*size(a),class(a));
  b(1:2:2*ha,1:2:2*wa)=a;
  b(2:2:2*ha,1:2:2*wa)=a;
  b(1:2:2*ha,2:2:2*wa)=a;
  b(2:2:2*ha,2:2:2*wa)=a;
endfunction
