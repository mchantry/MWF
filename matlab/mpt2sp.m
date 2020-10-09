function [s] = mpt2sp(mp,Lx,Lz)
  %Converts mean-poloidal-torroidal to u,v,w
  
mpRe=squeeze(mp(:,:,:,1));
mpIm=squeeze(mp(:,:,:,2));

K0=(size(mpRe,1)+1)/2
MT=size(mpRe,2);
if (mod(MT,2)==0)
    MM=MT/2
else
    MM=(MT+1)/2
end
NT=size(mpRe,3)

%First derivative operator x 
ad_m1=[0:MM-1,-MT+MM:1:-1]'*(2*pi/Lx);
tmp1=[0:MM-1];
tmp2=[MT-MM:-1:1];
%second derivative operator x 
ad_m2=-[tmp1.^2,tmp2.^2]'*(2*pi/Lx)^2;

%Same for y
ad_n1=[0:NT-1]'*(2*pi/Lz);
ad_n2=([0:NT-1]'*(2*pi/Lz)).^2;

sRe=zeros(3*K0-1,MT,NT);
sIm=zeros(3*K0-1,MT,NT);
d_beta = pi/2;
ad_k1=zeros(K0);
for kk = 1 : K0
  k=kk-1;
  ad_k1(kk)=(-1.0)^k * k * d_beta;
end
for n=1:NT
    nn=n;
    for m=1:MT
        for k=1:K0
            k1=k;
            k2=k+K0-1;
            k3=k+2*K0-1;
            sRe(k1,m,n)=-(-mpIm(k1,m,n)*ad_n1(nn) ...
                 +mpIm(k2,m,n)*ad_k1(k)*ad_m1(m));
            sIm(k1,m,n)= (-mpRe(k1,m,n)*ad_n1(nn) ...
                 +mpRe(k2,m,n)*ad_k1(k)*ad_m1(m));

            sRe(k2,m,n)=-mpRe(k2,m,n)*(ad_m2(m)+ad_n2(nn));
            sIm(k2,m,n)=-mpIm(k2,m,n)*(ad_m2(m)+ad_n2(nn));

            sRe(k3,m,n)=-(mpIm(k1,m,n)*ad_m1(m) ...
                 +mpIm(k2,m,n)*ad_k1(k)*ad_n1(nn));
            sIm(k3,m,n)= (mpRe(k1,m,n)*ad_m1(m) ...
                 +mpRe(k2,m,n)*ad_k1(k)*ad_n1(nn));
        end
    end
end
%Constrain mean fields to be zero
sRe(1,1,1)=0;
sIm(1,1,1)=0;
sRe(2*K0,1,1)=0;
sIm(2*K0,1,1)=0;
for k=1:K0-1
   k1=k+1;
   k2=k+K0;
   k3=k+2*K0;
   sRe(k1,1,1)= mpRe(k1,1,1);
   sIm(k1,1,1)= 0;
   sRe(k3,1,1)= mpRe(k2,1,1);
   sIm(k3,1,1)= 0;
end
s=sRe+1i*sIm;
end

