function [xxx,zzz,U] = GridUy(filename,vel_field,yloc)
%Produces a velocity field on a y-plane.
  
% Location (-1 to 1)
%yloc=1.;

% U,V,W
%vel_field = 'U';
if vel_field == 'U'
    ymode = 1;
elseif vel_field == 'W'
    ymode = 8;
elseif vel_field == 'V'
    ymode = 5;
else
  disp("not valid vel_field, returning")
  return
end
    
%Load file
[Lx,Lz,mpt] = readmpt(filename);
%Convert to u,v,w
s=mpt2sp(mpt,Lx,Lz);

%Extract velocity field at y
if vel_field == 'U' || vel_field == 'W'
spy=squeeze(s(ymode,:,:)...
    +sin(pi/2*yloc)*s(ymode+1,:,:)...
    +cos(pi*yloc)*s(ymode+2,:,:)...
    +sin(3*pi/2*yloc)*s(ymode+3,:,:)...
        );
 else
   spy=squeeze(cos(pi/2*yloc)*s(ymode,:,:)...
    +sin(pi*yloc)*s(ymode+1,:,:)...
    +cos(3*pi/2*yloc)*s(ymode+2,:,:)...
        );
end
%Transform in x
nnx = ceil(size(spy,1)/2);
spy_pad = cat(1,spy(1:nnx,:),zeros(1,size(spy,2)),spy(end-nnx+2:end,:));
gxsz=ifft(spy_pad,[],1)*size(spy_pad,1);
nnz = size(spy,2);
%Tranform in z
U=ifft(gxsz,2*nnz,2,'symmetric')*2*nnz;
Nx = size(U,1)
Nz = size(U,2)
%Transpose U so x-direction is horizontal
U=U';
%Make x & z grid
xx=[0:Nx-1]*Lx/Nx;
zz=[0:Nz-1]*Lz/Nz;
[xxx,zzz]=meshgrid(xx,zz);

end
