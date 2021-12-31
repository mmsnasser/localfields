% Example_200ncc_Lam
% 25-9-2021
clear
addpath bie fmm strip
%%
av =  [0.0002,0.001,0.002:0.002:0.046,0.047,0.0475,0.048,0.0484,0.0487,0.049,0.0493,0.0495,0.04965,0.0498].';
%%
for kk=1:length(av)
    %
a     =  av(kk);  b = a;
%
[x,y] = meshgrid(-0.95:0.1:0.95,0.05:0.1:0.95);
z     = x+i*y;
ecent = z(:);
%
alphas    =  0.5i;
%
n         =  2^11;
%
ethet  =  zeros(size(ecent));
ae     =  a+zeros(size(ecent));  
be     =  b+zeros(size(ecent));
ell    =  length(ae);
p      =  0;
%
m      =  ell+p;
t      = (0:2*pi/n:2*pi-2*pi/n).';
%
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    zet_e(Jk,1)  =  ecent(k)+exp(i*ethet(k)).*(ae(k).*cos(t)-i*be(k).*sin(t));
    zetp_e(Jk,1) =  exp(i*ethet(k)).*(-ae(k).*sin(t)-i*be(k).*cos(t));
end
%
zet_c =[]; zetp_c = [];
for k=1:p
    Jk = (k-1)*n+1:k*n;
    zet_c(Jk,1)  =  ccent(k)+crad(k).*exp(-i.*t);
    zetp_c(Jk,1) =  crad(k)*(-i).*exp(-i.*t);
end
%
zet  = [zet_e  ; zet_c];
zetp = [zetp_e ; zetp_c];
%
figure(1)
clf
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
drawnow
if a==0.039
    print -depsc  fig200CNTncc
end
%%
tic
u  =  Tofz(zet,zetp,alphas,n,ell,p);
toc
mu0  =  u.mu0;
%%
t1     =  2*pi-2*atan(exp(pi));
t2     =  2*pi-2*atan(exp(-pi));
mu12   =  tripoly(mu0,[t1;t2]);
area   =  m*a*b*pi
Lam_y  =  1+0.5*(mu12(2)-mu12(1))
% 
concv(kk,1) = area/2;
Lamyv(kk,1) = Lam_y;
%
end
%%
[concv Lamyv]
%%
lame =  (1+concv)./(1-concv);
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% 
plot(concv,Lamyv,'-k','LineWidth',1.5)
hold on; box on
plot(concv,lame,'-.b','LineWidth',1.5)
% plot(pi/4,0.8,'pk','LineWidth',1.5)
plot([pi/4 pi/4],[0 33],':k','LineWidth',1.5)
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
% 
legend({'$\lambda_y$','$\lambda_e$'},...
    'Location','northwest','Interpreter','latex')
% text(pi/4-0.015,1.8,'$\frac{\pi}{4}$','Interpreter','latex','FontSize',18)
xlabel('$c\,$: concentration')
% ylabel('$\lambda_y$','Interpreter','latex')
xticks([0:0.1:0.8])
yticks([0:8:32])
%
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0 0.8 0.0 33])
print -depsc  figLam200CNTec
%%