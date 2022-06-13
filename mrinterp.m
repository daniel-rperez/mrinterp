% Multiresolution interpolation function
% [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterp(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)
% x, y, z: coordinates and value
% dxy: spatial resolution of the raster (same for X and Y)
% slev: first interpolation level (optional) and smoothing level (number of levels without exact interpolation)
% mapa: interpolated raster matrix
% xmin, xmax, ymin, ymax: extent
% mapa0: raster with the projection of the original data points
function [mapa, xmin,xmax, ymin,ymax, mapa0]=mrinterp(x,y,z, dxy, slev, xmin,xmax,ymin,ymax)

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
    if( l <= Nl - slev0 - slev1 ) %smoothing level
      if( b > 1 )
        for i=1:ni
          for j=1:nj
            maplN(i,j)=sum2( map0N(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
            maplS(i,j)=sum2( map0S(1+(i-1)*b:i*b, 1+(j-1)*b:j*b) );
          endfor
        endfor
      else
        maplN=1.0*map0N;
        maplS=map0S;
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
      m=(maplN>0); %where data exists at current resolution
      mapaN(m)=maplN(m);
      %mapaS(m)=maplS(m);
      mapaz(m)=maplS(m)./maplN(m);
      
      %compute smoothed version of low resolution variables
      cmapaM=conv2(ma,mK, 'same');
      cmapaN=conv2(mapaN,mK, 'same');
      %cmapaS=conv2(mapaS,mK, 'same');
      cmapaZ=conv2(mapaz,mK, 'same');
      
      %fill in missing data at pixels with no data at current resolution
      m=~m & (cmapaM>0); %where data does not exist at current resolution but it does at lower 
      maplN(m)=cmapaN(m) ./ cmapaM(m);
      maplS(m)=(cmapaZ(m) .* cmapaN(m)) ./ power(cmapaM(m),2);

      clear m ma
      clear cmapaM cmapaN cmapaZ %cmapaS
      clear mapaz mapaN mapaS
      
    endif
    
    %next resolution
    if( l<Nl-slev0)
      b=b/2;
      %disp( [ 'PROGRESS: b=' num2str(b) ] );
      mapaS=dupl2(maplS)/4;
      mapaN=dupl2(maplN)/4;
    endif

  endfor

  %return mapa
  mapa=nan(size(maplS), dcls);
    mapa(maplN>0) = maplS(maplN>0) ./ maplN(maplN>0);
  
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
