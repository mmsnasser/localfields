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
    z(k)=(-1+2*rand)+i*rand;
end
%
for k=1:length(z)
    if abs(imag(z(k))-0.5)>0.48
        z(k) = NaN+i*NaN;
    end
    if abs(real(z(k)))>0.99
        z(k) = NaN+i*NaN;
    end
end
z(abs(z-alphas)<0.03) = NaN+i*NaN;
z = z(abs(z)>=0);
for k=1:length(z)-1
    for j=k+1:length(z)
        if abs(z(k)-z(j))<=0.02
            z(j)=NaN;
        end
    end
end
z=z(abs(z)>=0);
% 
cente = z(1:1000);
centc = z(1001:2000);
% 
length(z)
% 
ell     =  length(cente);
p       =  length(centc);
m       =  ell+p;
%
zet_e =[]; zetp_e = [];
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    zet_e(Jk,1)  =  cente(k)+r.*exp(-i.*t);
    zetp_e(Jk,1) =  r*(-i).*exp(-i.*t);
end
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
% print -depsc -r1000 fig2000NV
drawnow
%%
tic
u  =  Tofz (zet,zetp,alphas,n,ell,p);
toc
mu0  =  u.mu0;
%
t1     =  2*pi-2*atan(exp(pi));
t2     =  2*pi-2*atan(exp(-pi));
mu12   =  tripoly(mu0,[t1;t2]);
Lam_y  =  1+0.5*(mu12(2)-mu12(1))
%
Lam_yv(kk,1) = Lam_y;
%
end
%%
[ind Lam_yv]
c1     =  (m/2)*pi*r^2/2;
c2     =  (m/2)*pi*r^2/2;
Lam_e  = (1+c1-c2)./(1-c1+c2)
Lam_ev =  Lam_e+zeros(size(Lam_yv));
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(ind,Lam_yv,'b','LineWidth',1.0)
hold on; box on
plot(ind,Lam_ev,'r','LineWidth',1.0)
% 
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
% 
xticks([0:2:20])
xlabel('Number of the experiment','interpret','latex')
%
legend({'$\lambda_y$','$\lambda_e$'},...
    'Location','northeast','Interpreter','latex')
%
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([1 20 0.975 1.025])
print -depsc fig2000Lam
%%