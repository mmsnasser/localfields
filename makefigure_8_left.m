% Example_5_nch_Lam
% 25-9-2021
clear
addpath bie fmm strip
%%
av =  [0.0001,0.001,0.01:0.01:0.18,0.185,0.19,0.193,0.196,0.197,0.198,0.199,0.1995,0.1999].';
%%
for kk=1:length(av)
    %
a     =  av(kk);  b = 0.1*a;
Coefe = [
        -0.80+0.50i  0.00*pi  a  b
        -0.40+0.50i  0.00*pi  a  b
         0.00+0.50i  0.00*pi  a  b
         0.40+0.50i  0.00*pi  a  b
         0.80+0.50i  0.00*pi  a  b
        ];
%
alphas    =  0.2+0.2i;
%
n         =  2^11;
%
ecent  =  Coefe(:,1);
ethet  =  Coefe(:,2);
ae     =  Coefe(:,3);  
be     =  Coefe(:,4);
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
plot(real(alphas),imag(alphas),'dr','LineWidth',1.5)
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
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% 
plot(concv,Lamyv,'k','LineWidth',1.5)
hold on; box on
plot([pi/100 pi/100],[0 1.1],':k','LineWidth',1.5)
%
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
% 
% text(pi/100-0.0002,1.003,'$\frac{\pi}{100}$','Interpreter','latex','FontSize',18)
xlabel('$c\,$: concentration')
ylabel('$\lambda_y$','Interpreter','latex')
xticks([0:0.01:0.04])
%
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0 0.04 1 1.04])
print -depsc  figLam5CNTnch
%%