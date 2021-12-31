% Example_2000_both_rand
% 12-10-2021
clc;clear
addpath bie fmm strip
%%
n       =  2^11;
r       =  0.0075;
t       = (0:2*pi/n:2*pi-2*pi/n).';
alphas  =  0.5i;
%% 
for kk=1:20
    kk
    ind(kk,1) = kk;
for k=1:10000
    ze(k)=(-1+2*rand)+i*rand;
end
%
for k=1:length(ze)
    if abs(imag(ze(k))-0.5)>0.5-3*r
        ze(k) = NaN+i*NaN;
    end
    if abs(real(ze(k)))>1-2*r
        ze(k) = NaN+i*NaN;
    end
end
ze(abs(ze-alphas)<4*r) = NaN+i*NaN;
ze = ze(abs(ze)>=0);
for k=1:length(ze)-1
    for j=k+1:length(ze)
        if abs(ze(k)-ze(j))<=4.2*r
            ze(j)=NaN;
        end
    end
end
ze=ze(abs(ze)>=0);
% 
length(ze)
cente = ze(1:1000);
thete = 2*pi*rand(size(cente));
ell     =  length(cente);
zet_e =[]; zetp_e = [];
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    zet_e(Jk,1)  =  cente(k)+r*exp(i*thete(k)).*( 2.*cos(t)-0.5i.*sin(t));
    zetp_e(Jk,1) =           r*exp(i*thete(k)).*(-2.*sin(t)-0.5i.*cos(t));
end
%
%
%
for k=1:20000
    zc(k)=(-1+2*rand)+i*rand;
end
%
for k=1:length(zc)
    if abs(imag(zc(k))-0.5)>0.5-2*r
        zc(k) = NaN+i*NaN;
    end
    if abs(real(zc(k)))>1-r
        zc(k) = NaN+i*NaN;
    end
end
zc(abs(zc-alphas)<4*r) = NaN+i*NaN;
zc = zc(abs(zc)>=0);
for k=1:length(zc)-1
    for j=k+1:length(zc)
        if abs(zc(k)-zc(j))<=2.5*r
            zc(j)=NaN;
        end
    end
end
zc=zc(abs(zc)>=0);
mm1 = length(zc)
pause
for k=1:length(zc)-1
    k
    if min(abs(zc(k)-zet_e)) <=1.5*r
        '======================'
        zc(k)=NaN;
    end
end
zc=zc(abs(zc)>=0);
%
mm2 = length(zc)
%
centc = zc;%(1:1000);
% 
p       =  length(centc);
m       =  ell+p;
%
zet_c =[]; zetp_c = [];
for k=1:p
    Jk = (k-1)*n+1:k*n;
    zet_c(Jk,1)  =  centc(k)+r.*exp(-i.*t);
    zetp_c(Jk,1) =  r*(-i).*exp(-i.*t);
end
%
zet  = [zet_e  ; zet_c];
zetp = [zetp_e ; zetp_c];
%%
figure(1)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:ell
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'r','LineWidth',0.75)
end
for k=ell+1:m
    crv = zet((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',0.75)
end
plot([-2 2],[1 1],'k','LineWidth',1.2)
hold on
plot([-2 2],[0 0],'k','LineWidth',1.2)
% 
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
set(gcf,'Renderer','zbuffer')
print -depsc -r1000 fig2000EC
drawnow
%%
rcentc = real(centc);
icentc = imag(centc);
%
save('rcentcp.mat',  'rcentc',  '-ascii', '-double')
save('icentcp.mat',  'icentc',  '-ascii', '-double')
%
rcente = real(cente);
icente = imag(cente);
%
save('rcentep.mat',  'rcente',  '-ascii', '-double')
save('icentep.mat',  'icente',  '-ascii', '-double')
%%

%%