function [h1,h2,h3,h4,meanh1,stdh1,meanh2,stdh2,meanh3,stdh3,meanh4,stdh4]=Monte_Carlo(spnum,CalWin,sigma)
%This function is used to validate the accuracy of homogeneous pixel selection algorithm
%   Usage:
%       simulation(spnum,CalWin,sigma);
%   Inputs:
%   - spnum:    The stack size (page)
%   - CalWin:   The image size (height,width)
%   - sigma:    Noise-free amplitude

if nargin < 3
    sigma =200;
end

if nargin < 2
    CalWin =[15 15]; %img size
end

if nargin < 1
    spnum =30;
end

ratio=1:.2:3;
rep=10000; %replication
alpha =0.05;

%the true parameter of the Rayleigh distribution 
sigma=sqrt(sigma/2);
Matrix=sigma*ones(CalWin(1),CalWin(2)); 

InitRow=(CalWin(1)+1)/2; % InitRow is CenterRow
InitCol=(CalWin(2)+1)/2;
tic
%------------------------%
h1=zeros(length(ratio),rep); %BWS-DIE
h2=h1; %KS
h3=h1; %BWS
h4=h1; %FaSHP
kstemp=zeros(CalWin(1),CalWin(2));
CRSHP=finv([alpha/2,1-alpha/2],2*spnum,2*spnum);
num=1;
p=1;
all = length(ratio)*rep;
all_step = floor(all/10);
for ii=1:length(ratio)
    Matrix(:,InitCol+1:end)=ratio(ii)*sigma;
    NEWMatrix = repmat(Matrix,[1,1,spnum]);
    for jj=1:rep
        Noise = random('rayl', 1,[CalWin(1),CalWin(2),spnum]);    
        NoiseAdd =Noise.*NEWMatrix;
        Ref = NoiseAdd(InitRow,InitCol,:);
       %BWS-DIE
       htemp = BWSDIE(NoiseAdd.^2,CRSHP,alpha);
       h1(ii,jj)=sum(htemp(:));      
       %KS
        for ll=1:CalWin(1)
            for kk=1:CalWin(2)
                temp = NoiseAdd(ll,kk,:);
                kstemp(ll,kk)=kstest2(Ref(:),temp(:));
            end
        end
        h2(ii,jj)=sum(kstemp(:));
        %BWS
        Xarry = repmat(Ref(:),[1,CalWin(1)*CalWin(2)]);
        Yarray= reshape(NoiseAdd,[CalWin(1)*CalWin(2),spnum])';
        bwstemp  = BWStest(Xarry,Yarray,alpha);
        h3(ii,jj)=sum(bwstemp(:));       
        %FaSHP
        htemp = shp(NoiseAdd.^2,Ref.^2,CRSHP,alpha);
        h4(ii,jj)=sum(htemp(:));        
        
        num=num+1;
        if num == all_step * p;
            disp(['progress: ', num2str(10*p),'%']);
            p = p+1;
        end        
        
    end
end
h1=h1/CalWin(1)/CalWin(2);
h2=h2/CalWin(1)/CalWin(2);
h3=h3/CalWin(1)/CalWin(2);
h4=h4/CalWin(1)/CalWin(2);
toc

meanh1=mean(h1,2);
stdh1 =std(h1,0,2);
meanh2=mean(h2,2);
stdh2 =std(h2,0,2);
meanh3=mean(h3,2);
stdh3 =std(h3,0,2);
meanh4=mean(h4,2);
stdh4 =std(h4,0,2);
figure;plot(ratio,meanh1,'x--',ratio,meanh2,'+--',ratio,meanh3,'*--',ratio,meanh4,'o-.');grid on;legend('BWS-DIE','KS','BWS','FaSHP');ylabel('Mean rejection');xlabel('\sigma_1/\sigma_2');
figure;plot(ratio,stdh1,'x--',ratio,stdh2,'+--',ratio,stdh3,'*--',ratio,stdh4,'o-.');grid on;legend('BWS-DIE','KS','BWS','FaSHP');ylabel('Std. rejection');xlabel('\sigma_1/\sigma_2');

%%
function y1 = shp(data,ref,CRSHP,Alpha)
[L,P,npages] = size(data);
LRT_nl = 3; %7*7 window size
LRT_nw = 3; 

InitRow=(L+1)/2; % InitRow is CenterRow
InitCol=(P+1)/2;

Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);

%Initial estimation (LRT)
Matrix = data(InitRow-LRT_nl:InitRow+LRT_nl,InitCol-LRT_nw:InitCol+LRT_nw,:);
ref = mean(ref,3);
temp = mean(Matrix,3);
T = ref./temp;
T = T>CRSHP(1)&T<CRSHP(2);
SeedPoint = mean(temp(T));
%iteration (Gamma Confidence interval)
MeanMatrix = mean(data,3);
SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
y1 = ~SeedPoint;

function y1 = BWSDIE(data,CRSHP,Alpha)
[L,P,npages] = size(data);
LRT_nl = 3; %7*7 window size
LRT_nw = 3; 

InitRow=(L+1)/2; % InitRow is CenterRow
InitCol=(P+1)/2;
RadiusRow=(L-1)/2;
RadiusCol=(P-1)/2;

Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);
meandata=mean(data,3);
temp = meandata(InitRow-LRT_nl:InitRow+LRT_nl,InitCol-LRT_nw:InitCol+LRT_nw,:);
            T = meandata(InitRow,InitCol)./temp;
            T = T>CRSHP(1)&T<CRSHP(2);
            SeedPoint = mean(temp(T));
            if LRT_nl<RadiusRow && LRT_nw<RadiusCol
               MeanMatrix = meandata(InitRow-LRT_nl:InitRow+LRT_nl,InitCol-LRT_nw:InitCol+LRT_nw,:);
               tempPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages;
               SeedPoint=mean(MeanMatrix(tempPoint));
               LRT_nl=LRT_nl+1;
               LRT_nw=LRT_nw+1;
            end
            MeanMatrix = meandata(InitRow-RadiusRow:InitRow+RadiusRow,InitCol-RadiusCol:InitCol+RadiusCol);
            SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages;
            y1 = ~SeedPoint;
            
