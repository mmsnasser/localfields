% Example-50ne
% 25-9-2021
clear
addpath bie fmm strip
%%
diskrv =  [0.001,0.01:0.01:0.08,0.085,0.09,0.093,0.096,0.097,0.098,0.099].';
%%
for kk=1:length(diskrv)
    diskr =  diskrv(kk);
    %
Coefc = [ 
       -0.90+0.10i     diskr
       -0.70+0.10i     diskr
       -0.50+0.10i     diskr
       -0.30+0.10i     diskr
       -0.10+0.10i     diskr
        0.10+0.10i     diskr
        0.30+0.10i     diskr
        0.50+0.10i     diskr
        0.70+0.10i     diskr
        0.90+0.10i     diskr
        
       -0.90+0.30i     diskr
       -0.70+0.30i     diskr
       -0.50+0.30i     diskr
       -0.30+0.30i     diskr
       -0.10+0.30i     diskr
        0.10+0.30i     diskr
        0.30+0.30i     diskr
        0.50+0.30i     diskr
        0.70+0.30i     diskr
        0.90+0.30i     diskr
        
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
        
       -0.90+0.70i     diskr
       -0.70+0.70i     diskr
       -0.50+0.70i     diskr
       -0.30+0.70i     diskr
       -0.10+0.70i     diskr
        0.10+0.70i     diskr
        0.30+0.70i     diskr
        0.50+0.70i     diskr
        0.70+0.70i     diskr
        0.90+0.70i     diskr     
        
       -0.90+0.90i     diskr
       -0.70+0.90i     diskr
       -0.50+0.90i     diskr
       -0.30+0.90i     diskr
       -0.10+0.90i     diskr
        0.10+0.90i     diskr
        0.30+0.90i     diskr
        0.50+0.90i     diskr
        0.70+0.90i     diskr
        0.90+0.90i     diskr
        ];
%
alphas    =  0.00+0.4i;
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
% if diskr ==0.09
%     print -depsc  fig50CNTne
% end
drawnow
%%
tic
u  =  Tofz (zet,zetp,alphas,n,ell,p);
toc
mu0  =  u.mu0;
%%
t1     =  2*pi-2*atan(exp(pi));
t2     =  2*pi-2*atan(exp(-pi));
mu12   =  tripoly(mu0,[t1;t2]);
area   =  p*pi*diskr^2
Lam_y  =  1+0.5*(mu12(2)-mu12(1))
% 
concv(kk,1) = area/2;
Lamyv(kk,1) = Lam_y;
%
end
%%
CMA = @(p,c)((1+p*c)./(1-p*c));
cmac = CMA(-1,concv);
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% 
plot(concv,Lamyv,'k','LineWidth',1.5)
hold on; box on
plot(concv,cmac,'-.b','LineWidth',1.5)
plot([pi/4 pi/4],[0 1],':k','LineWidth',1.5)
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
% 
legend({'$\lambda_y$','$\lambda_e$'},...
    'Location','southwest','Interpreter','latex')
%
xlabel('$c\,$: concentration')
% ylabel('$\lambda_y$','Interpreter','latex')
xticks([0:0.2:0.8])
%
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0 0.8 0 1])
print -depsc  figLam50CNTne
%%