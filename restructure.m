function restructure(baseline)
% This function is used to replece the results of stamps(4,4) 
% usage
%      restructure('SBAS')\restructure('PS')
% baseline:1)'SBAS';2)'PS'

load DPS.mat
load ps1.mat
ij=ij(DPS',:);
lonlat=lonlat(DPS',:);
xy=xy(DPS',:);
n_ps=length(DPS);

if strcmpi(baseline,'SBAS') 
   save('ps2.mat','bperp','day','ifgday','ifgday_ix','ij','ll0','lonlat','master_day','master_ix','mean_incidence','mean_range','n_ifg','n_image','n_ps','xy'); 
end
if strcmpi(baseline,'PS') 
   save('ps2.mat','bperp','day','ij','ll0','lonlat','master_day','master_ix','mean_incidence','mean_range','n_ifg','n_image','n_ps','xy');
end

clear all
load DPS.mat
load bp1.mat
bperp_mat=bperp_mat(DPS',:);
save('bp2.mat','bperp_mat');

clear all
load DPS.mat
load la1.mat
la=la(DPS',:);
save('la2.mat','la')

clear all
load DPS.mat
load hgt1.mat
hgt=hgt(DPS',:);
save('hgt2.mat','hgt')

clear all
load DPS.mat
load ph1.mat
ph=ph(DPS',:);
save('ph2.mat','ph')

clear all
load DPS.mat
load pm1.mat
C_ps=C_ps(DPS',:);
coh_ps=coh_ps(DPS',:);
K_ps=K_ps(DPS',:);
ph_patch=ph_patch(DPS',:);
ph_res=ph_res(DPS',:);
save('pm2.mat','C_ps','coh_ps','K_ps','ph_patch','ph_res');