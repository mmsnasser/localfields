% Example: plot the figures in the introduction
% 13-6-2020
clc;clear
addpath bie fmm strip
%%
Coefe = [
         0.70+0.50i  0.25*pi  0.30  0.030
         0.55+0.80i  0.75*pi  0.25  0.025
        -0.45+0.25i  0.50*pi  0.20  0.02
        -0.50+0.55i  0.00*pi  0.40  0.04
        ];
%
Coefc = [ 
        0.12+0.68i     0.10
        0.40+0.15i     0.12
       -0.90+0.30i     0.15
        ];
%
alphas    =  0.05+0.3i;
%
n         =  2^11;
%
ecent  =  Coefe(:,1);
ethet  =  Coefe(:,2);
ae     =  Coefe(:,3)./2;  
be     =  Coefe(:,4)./2;
ell    =  length(ae);
ccent  =  Coefc(:,1);  
crad   =  Coefc(:,2);
p      =  length(ccent);
%
t         = (0:2*pi/n:2*pi-2*pi/n).';
%
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
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
annotation('arrow',[0.52 0.52],[0.23 0.83],'LineWidth',1.0)
hold on
annotation('line',[0.52 0.52],[0.29 0.682],'LineWidth',1.5,'color',[0.9 0.9 0.9])
fill([4 -4 -4 4],[1 1 0 0], [0.9 0.9 0.9]);
plot([-4 4],[1 1],'k','LineWidth',1.2)
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    fill(real(zet_e(Jk,1)),imag(zet_e(Jk,1)),[0.8 0.8 1])
end
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    plot(real(zet_e(Jk,1)),imag(zet_e(Jk,1)),'b','LineWidth',0.5)
end
for k=1:p
    Jk = (k-1)*n+1:k*n;
    fill(real(zet_c(Jk,1)),imag(zet_c(Jk,1)),'w')
end
for k=1:p
    Jk = (k-1)*n+1:k*n;
    plot(real(zet_c(Jk,1)),imag(zet_c(Jk,1)),'k','LineWidth',0.5)
end
% plot(real(alphas),imag(alphas),'pk','MarkerFaceColor','k')

axis equal
set(gca,'LooseInset',get(gca,'TightInset'))

annotation('arrow',[0.05 0.999],[0.29 0.29],'LineWidth',1.0)
text(+1.52,-0.07,'{$x$}','FontSize',14,'Interpreter','latex');
text(+0.03,+1.35,'{$y$}','FontSize',14,'Interpreter','latex');
% 
text(-0.75,-0.07,'{$T=1$}','FontSize',14,'Interpreter','latex');
text(-0.75,+1.05,'{$T=0$}','FontSize',14,'Interpreter','latex');
% 
% text(+0.05,+0.25,'{$\alpha$}','FontSize',14,'Interpreter','latex');
text(-0.21,+0.45,'{$\Delta T=0$}','FontSize',14,'Interpreter','latex');
% 
text( 0.53,+0.15,'{$\frac{\partial T}{\partial {\bf n}}=0$}','FontSize',14,'Interpreter','latex');
text( 0.00, 0.86,'{$\frac{\partial T}{\partial {\bf n}}=0$}','FontSize',14,'Interpreter','latex');
text(-1.40,+0.30,'{$\frac{\partial T}{\partial {\bf n}}=0$}','FontSize',14,'Interpreter','latex');
% 
text(+0.75,+0.45,'{$T=\delta_1$}','FontSize',14,'Interpreter','latex');
text(+0.56,+0.85,'{$T=\delta_2$}','FontSize',14,'Interpreter','latex');
text(-0.42,+0.25,'{$T=\delta_3$}','FontSize',14,'Interpreter','latex');
text(-0.60,+0.62,'{$T=\delta_4$}','FontSize',14,'Interpreter','latex');
% 
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.6 1.6 -0.2 1.4])
print(gcf,'-depsc'  ,'dom_Om.eps')
%%

%%
eto  = exp(i.*t);
etop = i.*exp(i.*t);
alpha = Phiv(alphas);
%
if ell>0
    et_e = Phiv(zet_e);
    etp_e = Phivp(zet_e).*zetp_e;
end
%
if p>0
    et_c = Phiv(zet_c);
    etp_c = Phivp(zet_c).*zetp_c;
end
%%
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%
fill(real(eto),imag(eto), [0.9 0.9 0.9]);
hold on
plot(real(eto),imag(eto),'k','LineWidth',1.2)
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    fill(real(et_e(Jk,1)),imag(et_e(Jk,1)),[0.8 0.8 1])
end
for k=1:ell
    Jk = (k-1)*n+1:k*n;
    plot(real(et_e(Jk,1)),imag(et_e(Jk,1)),'b','LineWidth',0.5)
end
for k=1:p
    Jk = (k-1)*n+1:k*n;
    fill(real(et_c(Jk,1)),imag(et_c(Jk,1)),'w')
end
for k=1:p
    Jk = (k-1)*n+1:k*n;
    plot(real(et_c(Jk,1)),imag(et_c(Jk,1)),'k','LineWidth',0.5)
end
% plot(real(alpha),imag(alpha),'pk','MarkerFaceColor','k')
axis equal
% 
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.02 1.02 -1.02 1.02])
print(gcf,'-depsc'  ,'dom_G.eps')
%%
