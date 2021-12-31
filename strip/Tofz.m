function u = Tofz(zet,zetp,alphas,n,ell,p,z)
%
%
%
%%
t     = (0:2*pi/n:2*pi-2*pi/n).';
m     =  ell+p;
%
eto   = exp(i.*t); etop  = i.*exp(i.*t);
alpha = Phiv(alphas);
%
et    = [eto];   etp   = [etop];
%
if ell>0
    Je     =  1:ell*n;
    zet_e  =  zet(Je);      zetp_e =  zetp(Je);
    et_e   =  Phiv(zet_e);  etp_e  =  Phivp(zet_e).*zetp_e;
    et     = [et ;et_e ];   etp    = [etp;etp_e];
end
%
if p>0
    Jp     =  ell*n+1:(ell+p)*n;
    zet_c  =  zet(Jp);      zetp_c =  zetp(Jp);
    et_c   =  Phiv(zet_c);  etp_c  =  Phivp(zet_c).*zetp_c;
    et     = [et ;et_c ];   etp    = [etp;etp_c];
end
%%
figure(2);
hold on
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    plot(real(et(Jk,1)),imag(et(Jk,1)),'b','LineWidth',1.2)
end
plot(real(alpha),imag(alpha),'pk','MarkerFaceColor','k')
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-1.02 1.02 -1.02 1.02])
%%
thetk =  zeros(m+1,1);
thetk(ell+2:m+1) = pi/2;
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    thet(Jk,1) = thetk(k);
end
%
A  =  exp(-i.*thet).*(et-alpha);
%
gam = zeros((m+1)*n,1);
%
for k=2:ell+1
    Jk      =  (k-1)*n+1:k*n;
    gam(Jk) = -(1/pi).*imag(clog((1-et(Jk))./(1+et(Jk))));
end
for k=ell+2:ell+p+1
    Jk      =  (k-1)*n+1:k*n;
    gam(Jk) = (1/pi).*log(abs((1-et(Jk))./(1+et(Jk))));
end
%
[mun,h] =  fbie(et,etp,A,gam,n,5,[],1e-13,200);
u.mu0   =  mun(1:n);
%%
if( nargin == 7 )
    %
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    hk(k,1)=mean(h(Jk));
end
%
c = -hk(1);
for k=2:ell+1
    delt(k-1,1) =  hk(k)+0.5+c;
end
for k=ell+2:ell+p+1
    delt(k-1,1) =  hk(k);
end
% 
g       = (gam+h+i*mun)./A;
f       = (et-alpha).*g+c;
for k=1:m+1
    Jk        = (k-1)*n+1:k*n;
    fpe(Jk,1) =  derfft(real(f(Jk)))+i*derfft(imag(f(Jk)));
end
fp  =  fpe./etp;
%
fop  = @(z)((i/pi).*(1./(1-z)+1./(1+z)));
%
w    = Phiv(z);
figure(2);
plot(real(w),imag(w),'.r')
%
ww    =  w(:).';
wv    =  ww(abs(ww)>=0);
% 
fw    =  fcau(et,etp,f,wv);
uw    =  real(fw);
uow   = (1/pi).*imag(clog((1-wv)./(1+wv)))+1/2;
Uwv   =  uow+uw;
Uww   =  NaN(size(ww));
Uww(abs(ww)>=0)=Uwv;
U     =  NaN(size(w));
U(:)  =  Uww;
% 
fpw   =  fcau(et,etp,fp,wv);
Fpw   = (fpw+fop(wv))./Phip(wv);
Fww   = (1+i)*NaN(size(ww));
Fww(abs(ww)>=0)=Fpw;
Fp    =  (1+i)*NaN(size(w));
Fp(:) =  Fww;
%
u.T   =  U;
u.Fp  =  Fp;
% 
end
%%
end