% Example-2ne
% 10-7-2020
clear
addpath bie fmm strip
%%
phiv =  [0.4:0.1:1.3]';
cv   =  [0.1:0.1:0.5].';
%%
n    =  2^11;
t    = (0:2*pi/n:2*pi-2*pi/n).';
%
x = [-0.95:0.05:0.95].';
z = [x+0.125i;x+0.375i;x+0.625i;x+0.875i];
x = [-0.975:0.05:0.975].';
z = [z;x+0.25i;x+0.5i;x+0.75i];
%
m = length(z);
%%
for jj=1:length(cv)
    jj
    c = cv(jj); 
    for kk=1:length(phiv)
        kk
        %
        phi  = phiv(kk); b = sqrt(2*phi/m); ep =  2*c/(m*pi*b^2);  a = ep*b;
        if a>0.024
            abcd;
        end
        %
        ecent = z(:);
        %
        alphas    =  0.5i;
        %
        %
        ethet  =  zeros(size(ecent));
        ae     =  a+zeros(size(ecent));  
        be     =  b+zeros(size(ecent));
        ell    =  length(ae);
        p      =  0;
        %
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
        %%
figure(1)
clf
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
% plot(real(alphas),imag(alphas),'pk','LineWidth',1.5)
plot([-2 2],[1 1],'k','LineWidth',1.2)
hold on
plot([-2 2],[0 0],'k','LineWidth',1.2)
% 
xticks([-1.5:0.5:1.5])
axis equal
box on
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.5 1.5 -0.2 1.2])
set(gcf,'Renderer','zbuffer')
print -depsc -r1000 fig276
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
        Lam_y  =  1+0.5*(mu12(2)-mu12(1))
        % 
        Lamyv(kk,jj) = Lam_y;
        %
    end
end
%%
[phiv Lamyv]
%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% 
plot(phiv,Lamyv(:,1),'k','LineWidth',1.5)
hold on
plot(phiv,Lamyv(:,2),'r','LineWidth',1.5)
plot(phiv,Lamyv(:,3),'color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
plot(phiv,Lamyv(:,4),'m','LineWidth',1.5)
plot(phiv,Lamyv(:,5),'b','LineWidth',1.5)
hold on; box on
%
legend({'$c=0.1$','$c=0.2$','$c=0.3$','$c=0.4$','$c=0.5$'},...
    'Location','northwest','Interpreter','latex')
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
% 
xlabel('$\phi$','Interpreter','latex')
ylabel('$\lambda_y$','Interpreter','latex')
xticks([0.4:0.1:1.8])
%
set(gca,'FontSize',18)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0.4 1.3 0 15])
% set(gcf,'Renderer','zbuffer')
% print -depsc -r1000 figLam276
print -depsc figLam276
%%