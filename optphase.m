function [pcoh,optintf]=optphase(intfstack,SHP)
% This function is used to optimize the original phase
%   Usage:
%      [pcoh,optintf]=optphase(intfstack,SHP);
%   Inputs:
%   - intfstack:  The stack of Interferogram
%   - SHP: Homogeneous point set
%   Outputs:
%   - pcoh: Posterior coherence
%   - optintf:  Optimized phase

%make data list
intfstack=(intfstack.datastack);
[nlines,nwidths,npages]=size(intfstack);

% if isreal(intfstack)  
%     intfstack = exp(1i*intfstack);
% else
%     intfstack(intfstack~=0) = intfstack(intfstack~=0)./abs(intfstack(intfstack~=0));
% end
% tic;

%Add a null value interferogram named intf
npages=npages+1;
intf=zeros(nlines,nwidths);
intf=exp(1i*intf);
tempintf = zeros(nlines,nwidths,npages,'single');
for pp=1:npages-1
tempintf(:,:,pp)=intfstack(:,:,pp);
end
tempintf(:,:,npages)=intf;
intfstack=tempintf;
intfstack=conj(intfstack);
ref=zeros(npages,1);
ref(npages,1)=1;
ref_idx=logical(ref);
%构建遍历窗口
FiltWin=[15 15];
CalWin=SHP.CalWin;
if FiltWin(1)> CalWin(1)
   FiltWin(1) = CalWin(1);
end
if FiltWin(2)> CalWin(2)
   FiltWin(2) = CalWin(2);
end
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2; 
  
intfstack = padarray(intfstack,[RadiusRow RadiusCol],'symmetric');
temp_mask = true(npages,npages);
temp_mask = triu(temp_mask,1);
shp=SHP.PixelInd';
%**************************************************************************
pcoh = zeros(nlines,nwidths,'single');
optintf = zeros(nlines,nwidths,npages,'single');
coh_pl = optintf;
num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/10);
num_shp=1;
for jj = 1:nwidths
    for kk= 1:nlines
        x_global  = jj+RadiusCol;
        y_global  = kk+RadiusRow;
        Z = intfstack(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol,:);
       
                Z = reshape(Z,[CalWin(1)*CalWin(2),npages]).';
                Z = Z(:,shp(num_shp,:));
                CpxCoh = double(Z*Z'./size(Z,2));
                Coh = abs(CpxCoh);
                
                coh_pl(kk,jj,:) = Coh(:,ref_idx);
                %positive definite detection
                [~,r] = chol(Coh); 
                e=1e-6;
                while ~(r == 0 && rank(Coh) == npages)
                    Coh = Coh+eye(npages)*e;
                    [~,r] = chol(Coh);
                    e=2*e;
                end   
                
                %phase optimization
                [V,~]=eig(inv(Coh).*CpxCoh);
                V=V(:,1);
                optintf(kk,jj,:)=V/V(ref_idx);
                %posteriori coherence
                rho = exp(1j*(angle(CpxCoh)-angle(V*V'))); 
                rho = rho(temp_mask);
                rho = sum(rho); 
                pcoh(kk,jj) = rho;
                num_shp=num_shp+1;   
                num=num+1;     
            if num == all_step * p
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end
    end
           
end
pcoh = real(pcoh)*2/npages/(npages-1);
optintf = coh_pl.*exp(-1j*angle(optintf)); 
end