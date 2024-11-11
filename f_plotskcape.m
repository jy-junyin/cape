function [CAPE,CIN,PLCL,LFC,LNB] = f_plotskcape(Tsrd,qsrd,plev)
g=9.81;
R=287;
cp=1005;
P0=10^5;
Rd=287;
theta=Tsrd(1);
q=qsrd(1);

Psrd=plev;
[ zLCL, TLCL ] = simpleLCL( theta,q );
PLCL=P0*(TLCL./theta).^(cp/R); % good
Tsrdx=Tsrd;
qsrdx=qsrd;
[ TAD, PAD ] = Adiabat( TLCL,PLCL );
Tsrd=interp1(Psrd,Tsrd,PAD);
qsrd=interp1(Psrd,qsrd,PAD);
Psrd=PAD;
[ CAPE,CIN,LFC,LNB,TADV,TsrdV ] = getAllV (Tsrd,Psrd,TAD,PAD,qsrd,q,PLCL);

TsrdxV=Tsrdx.*(1+0.61*qsrdx);
skTsrdxV=TsrdxV-40.*log(plev./P0);
skTsrdV=interp1(plev,skTsrdxV,Psrd);
skTADV=TADV-40.*log(Psrd./P0);
semilogy(skTsrdV,Psrd./100,'k','linewidth',0.8)
hold on
semilogy(skTADV,Psrd./100,'b--','linewidth',0.8)
set(gca, 'Ydir', 'reverse')
ylim([25000,100000]./100)
A=skTADV-skTsrdV;
LFC=Psrd(find(A>0,1));
A(1:find(A>0,1))=0;
LNB=Psrd(find(A<0,1)-2);
idx=Psrd<LFC & Psrd>LNB;
ypoints=Psrd(idx);
xupper=skTADV(idx);
xlower=skTsrdV(idx);
patch([xupper,xlower(end:-1:1)],[ypoints,ypoints(end:-1:1)]./100,[250,159,181]./255)
% xlim([305,350])

plot([280,310],[PLCL,PLCL]./100,':','color',[28,144,153]./255)
plot([280,310],[LFC,LFC]./100,':','color',[28,144,153]./255)
plot([280,310],[LNB,LNB]./100,':','color',[28,144,153]./255)
xlabel('$T_v$ (K)','interpreter','latex')


tt=280:310;
pp=100000:-500:25000;
[PP,TT]=meshgrid(pp,tt);
A=TT+40.*log(PP./100000);
contour(TT,PP./100,A,300:-10:230,'k:');
set(gca, 'Ydir', 'reverse')
set(gca,"YScale",'log')
end


