% Example-30ne
% 25-9-2021
clc;clear
addpath bie fmm strip
%%
diskr =  0.099;
Coefc = [ 
       -0.90+0.50i     diskr
       -0.70+0.50i     diskr
       -0.50+0.50i     diskr
       -0.30+0.50i     diskr
       -0.10+0.50i     diskr
        0.10+0.50i     diskr
        0.30+0.50i     diskr
        0.50+0.50i     diskr
        0.70+0.50i     diskr
        0.90+0.50i     diskr
       -0.90+0.75i     diskr
       -0.70+0.75i     diskr
       -0.50+0.75i     diskr
       -0.30+0.75i     diskr
       -0.10+0.75i     diskr
        0.10+0.75i     diskr
        0.30+0.75i     diskr
        0.50+0.75i     diskr
        0.70+0.75i     diskr
        0.90+0.75i     diskr
       -0.90+0.25i     diskr
       -0.70+0.25i     diskr
       -0.50+0.25i     diskr
       -0.30+0.25i     diskr
       -0.10+0.25i     diskr
        0.10+0.25i     diskr
        0.30+0.25i     diskr
        0.50+0.25i     diskr
        0.70+0.25i     diskr
        0.90+0.25i     diskr
        ];
%
alphas    =  0.00+0.375i;
%
n         =  2^11;
%
ell    =  0;
ccent  =  Coefc(:,1);  
crad   =  Coefc(:,2);
p      =  length(ccent);
%
m      =  ell+p;
t      = (0:2*pi/n:2*pi-2*pi/n).';
%
zet_e  = []; zetp_e = [];
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    zet_e(Jk,1)  =  ecent(k)+exp(i*ethet(k)).*(ae(k).*cos(t)-i*be(k).*sin(t));
    zetp_e(Jk,1) =  exp(i*ethet(k)).*(-ae(k).*sin(t)-i*be(k).*cos(t));
end
%
for k=1:p
    Jk = (k-1)*n+1:k*n;
    zet_c(Jk,1)  =  ccent(k)+crad(k).*exp(-i.*t);
    zetp_c(Jk,1) =  crad(k)*(-i).*exp(-i.*t);
end
%
zet  = [zet_e  ; zet_c];
zetp = [zetp_e ; zetp_c];
%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:ell
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'r','LineWidth',1.5)
end
for k=ell+1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',1.5)
end
% plot(real(alphas),imag(alphas),'dr','LineWidth',1.5)
plot([-2 2],[1 1],'k','LineWidth',1.2)
hold on
plot([-2 2],[0 0],'k','LineWidth',1.2)
% 
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',18)
% set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
%%
[x,y] =  meshgrid(linspace(-1.6,1.6,2001),linspace(0.0001,0.9999,1001));
z     =  x+i*y;
for k=1:ell
    ze  =  exp(-i*ethet(k)).*(z-ecent(k));
    z((real(ze)./ae(k)).^2+(imag(ze)./be(k)).^2<1-1e-6)=NaN+i*NaN;
end
for k=1:p
    zc  =  z-ccent(k);
    z(abs(zc./crad(k))<1-1e-6)=NaN+i*NaN;
end
%
figure(1)
plot(real(z),imag(z),'.r')
%%
tic
u  =  Tofz (zet,zetp,alphas,n,ell,p,z);
toc
T    =  u.T;
Fp   =  u.Fp;
mu0  =  u.mu0;
q    = -conj(Fp);
Ty   = -imag(Fp);
%%
t1     =  2*pi-2*atan(exp(pi));
t2     =  2*pi-2*atan(exp(-pi));
mu12   =  tripoly(mu0,[t1;t2]);
Lam_y  =  1+0.5*(mu12(2)-mu12(1))
%%
cvc = [0:0.1:1];
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
contourf(real(z),imag(z),T,cvc)
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0 1])
colorbar 
for k=1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',0.5)
end
plot([-2 2],[1 1],'k','LineWidth',1.2)
hold on
plot([-2 2],[0 0],'k','LineWidth',1.2)
% 
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',18)
% set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
set(gcf,'Renderer','zbuffer')
print -depsc -r1000 figT30CNTne
%%
cvc = [0.0001,0.2,0.5,0.9,1.1,1.5,2];
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
contourf(real(z),imag(z),abs(q),cvc)
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0 2])
colorbar 
for k=1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',0.5)
end
plot([-4 4],[1 1],'k','LineWidth',1.2)
hold on
plot([-4 4],[0 0],'k','LineWidth',1.2)
axis equal
box on
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',18)
% set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
set(gcf,'Renderer','zbuffer')
print -depsc -r1000 figQ30CNTne
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
Z=z(1:50:end,1:25:end);
Q=q(1:50:end,1:25:end);
qq=quiver(real(Z),imag(Z),real(Q)./abs(Q),imag(Q)./abs(Q),0.25,'LineWidth',0.5,'MaxHeadSize',0.33);
for k=1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',0.5)
end
plot([-4 4],[1 1],'k','LineWidth',1.2)
hold on
plot([-4 4],[0 0],'k','LineWidth',1.2)
axis equal
box on
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',18)
% set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
% print -depsc figQr30CNTne
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
PhasePlot(z,q,'m');
caxis([0 1])
colorbar
hold on
for k=1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'k','LineWidth',0.5)
end
plot([-4 4],[1 1],'k','LineWidth',1.2)
hold on
plot([-4 4],[0 0],'k','LineWidth',1.2)
axis equal
box on
set(gca,'FontSize',18)
% axis off
% set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.6 1.6 -0.2 1.2])
set(gcf,'Renderer','zbuffer')
%%