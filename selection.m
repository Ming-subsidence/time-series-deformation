function [DS,DPS]=selection(SHP,pcoh,pcoh_th,BroNum_ph)
% This function is used to select DS point and DPS point
% usage 
%      [DS,DPS]=selection(SHP,pcoh,pcoh_th,BroNum_ph)
% input:
%      - SHP: Homogeneous point set
%      - pcoh: Posterior coherence
%      - pcoh_th: 0.7 (Experience value)
%      - BroNum_ph: 0.5 (Experience value)

%select DS point
[ids,jds]=find(pcoh>pcoh_th&SHP.BroNum>BroNum_ph);
% imagesc(pcoh>pcoh_th&SHP.BroNum>SHP.BroNum_ph)

DSC=[ids,jds];
DSC=DSC-ones(size(DSC));

%Row and column numbers of all points in StaMPS format
load ps1.mat
rowcol=ij(:,2:3);

%DS point position in StaMPS
[~, iad, ibd] = intersect(rowcol,DSC,'rows'); 
DS=sort(iad,'ascend');

%DS+PS
load ps2.mat
[DPS, iap, ibp]=union(DS,xy(:,1),'rows');
