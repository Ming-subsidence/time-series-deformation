function [disp,v]=sortout()
% This function is used to export the results of the stamps' processing
% output
%     - disp : the subsidence of LOS (mm)
%     - v : the subsidence rate of LOS (mm/yr)


small_baseline_flag=getparm('small_baseline_flag',1);
ref_vel=getparm('ref_velocity',1);
lambda=getparm('lambda',1);

load psver
psname=['ps',num2str(psver)];
ps=load(psname);

if strcmpi(small_baseline_flag,'y')
    n_image=ps.n_image;
else
    n_image=ps.n_ifg;
end

%lon&lat
lonlat=ps.lonlat; 

%master image
master_ix=sum(ps.master_day>ps.day)+1;

%number of DPS points
ref_ps=ps_setref;

%Calculate phase
ps_plot('u-dmo',-1);
load ps_plot_u-dmo
ph_uw=ph_disp;
ph_uw=ph_uw-repmat(nanmean(ph_uw(ref_ps,:)),ps.n_ps,1);

[row,col]=size(ph_uw);
 
% Phase conversion to Los
ph=[];
 for i=1:col  
    phi=ph_uw(:,i);
    ph_sorti=sort(phi);
    min_phi=ph_sorti(ceil(length(ph_sorti)*0.001));
    max_phi=ph_sorti(floor(length(ph_sorti)*0.999));
    %ph(ph<min_ph)=min_ph;
    %ph(ph>max_ph)=max_ph;
    phi=-phi*lambda*1000/4/pi;
    ph(:,i)=phi;
 end

 ph_1=ph(:,1); 
 ph_2=ph(:,2:col);
 ph_3=ph_1-ph_1; 
 
 ph_4=[];
 for p=1:col-1
     ph_4(:,p)=ph_2(:,p)- ph_1;
 end
 ph=[ph_3,ph_4];
 disp=[lonlat,ph];%Los
 
 clear ph_disp
 ps_plot('V-do',-1);
 load ps_plot_V-do
 v=[lonlat,ph_disp];%rate
 save('disp.mat','disp')
 save('v.mat','v')
 clear all
 load v.mat
 load disp.mat
 %导出数据
 xlswrite('LOS.xlsx',disp);
 xlswrite('rate.xlsx',v);