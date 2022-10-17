% example
clc;clear
%set path
mkdir workpath
cd workpath
mkdir DIFF
mkdir MLI
mkdir StaMPSprocess

%data input
curpath  = pwd;
mlipath  = [curpath,filesep,'MLI'];
diffpath = [curpath,filesep,'DIFF'];
nlines = 600;% the row of Interferogram
mlistack =ImgRead(mlipath,'mli',nlines,'float32');
intfstack =ImgRead(diffpath,'diff',nlines,'cpxfloat32');
save('mlistack.mat','mlistack')
save('intfstack.mat','intfstack')

%SHP selection
CalWin = [15 15]; %[row col]
Alpha=0.05;
[SHP]=SHP_BWSDIE(mlistack,CalWin,Alpha);
save('SHP.mat','SHP')

%Phase optimization
load intfstack.mat
[pcoh,optintf]=optphase(intfstack,SHP);
save('optintf.mat','optintf')
save('pcoh.mat','pcoh')

%export optintf
load intfstack.mat
load optintf.mat

mkdir OPTINTF
cd OPTINTF
files = intfstack.filename;
len=length(files);

for i=1:len
    oldname1=num2str(files(i,1));
    oldname2=num2str(files(i,2));
    newname=strcat(oldname1,'_',oldname2, '.diff');
    count=fwritebk(optintf(:,:,i),[num2str(newname)],'cpxfloat32');
end

%select DS+PS point
cd ../
load SHP.mat
load pcoh.mat
%setparm
pcoh_th=0.7;
BroNum_ph=20;
[DS,DPS]=selection(SHP,pcoh,pcoh_th,BroNum_ph);
save('DS.mat','DS')
save('DPS.mat','DPS')

%Deformation:stamps(1,3);restructure stamps(4,4);stamps(5,8)
cd StaMPSprocess
stamps(1,3);
restructure('SBAS');
stamps(5,8)

%export subsidence and rate
[disp,v]=sortout;


% Monte Carlo simulation
CalWin = [15 15]; %[row col]
Alpha=0.05;
stacksize = 30;%samples
sigma  = 1000; %  Noise-free amplitude
[h1,h2,h3,h4,meanh1,stdh1,meanh2,stdh2,meanh3,stdh3,meanh4,stdh4,]=Monte_Carlo(stacksize,CalWin,sigma);

