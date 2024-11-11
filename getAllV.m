function [ CAPE,CIN, LFC,LNB,TADV,TsrdV ] = getAllV (Tsrd,Psrd,TAD,PAD,qsrd,q,PLCL)
%GETALL Summary of this function goes here
%   Detailed explanation goes here
R=287;
% synchronize the surrouding index
Tsrd=interp1(Psrd,Tsrd,PAD);
qsrd=interp1(Psrd,qsrd,PAD);
Psrd=PAD;

% virtual potential temperature
index1=(PAD>PLCL); % below LCL
index2=(PAD<=PLCL); % above LCL
es=ArdenBuck(TAD(index2)-273.15);
qs=0.622*es./PAD(index2);
TADV(index2)=TAD(index2).*(1+0.61*qs);
TADV(index1)=TAD(index1).*(1+0.61*q);

TsrdV=Tsrd.*(1+0.61*qsrd);

% find climate index
A=TADV-TsrdV;
index=( A>0 & PAD<PLCL);
if sum(index)>2
    X=A(index);
    PP=PAD(index);
    CAPE=-trapz(log(PP),X*R);
    LFC=PP(1);
    LNB=PP(end);
    
    index=(PAD>LFC);
    X=A(index);
    PP=PAD(index);
    if length(PP)<=1
        CIN=0;
    else
        CIN=trapz(log(PP),X*R);
    end
else
    CAPE=nan;
    LFC=nan;
    LNB=nan;
    CIN=nan;
end




end

