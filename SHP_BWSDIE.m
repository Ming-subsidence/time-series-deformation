function [SHP]=SHP_BWSDIE(mlistack,CalWin,Alpha)
%   Inputs:
%   - mlistack: A height by width by page matrix
%   - CalWin:   Fixed boxcar window size
%   - Alpha:    A value between 0 and 1 specifying the
%               significance level. Default is 0.05 for 5% significance.
%   Outputs:
%   - SHP.PixelInd: A CalWin(1)*CalWin(2) by size(mlistack,1)*size(mlistack,2) array with elements of type logical, containing a SHPs set per pixel 
%   - SHP.BroNum:   The SHP number per pixel (reference pixel is not included) 
%   - SHP.CalWin:   Fixed boxcar window size

mlistack=mlistack.datastack;
[nlines,nwidths,npages] = size(mlistack);
mlistack=single(mlistack);

%Parameter prepare:
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;
InitRow=(CalWin(1)+1)/2; % InitRow is CenterRow
InitCol=(CalWin(2)+1)/2; % InitCol is CenterCol

%Statistical threshold:
CR_lo = finv(Alpha/2,2*npages,2*npages);
CR_up = finv(1-Alpha/2,2*npages,2*npages);
Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);

%Edeg mirror-image Edeg
mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
meanmli = mean(mlistack,3);
[nlines_EP,nwidths_EP]= size(meanmli);
SHP.PixelInd=false(CalWin(1)*CalWin(2),nlines*nwidths); 

%estimate SHPs
num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/10);


    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow       
            %Initial estimation (Likelihood-ratio test)
            LRT_nl =3;
            LRT_nw =3;
            Matrix = mlistack(ll-LRT_nl:ll+LRT_nl,kk-LRT_nw:kk+LRT_nw,:);
            Ref = Matrix(LRT_nl+1,LRT_nw+1,:);
            T = BWStest(repmat(Ref(:),[1,(LRT_nl+LRT_nw+1)^2])...
                ,reshape(Matrix,[(LRT_nl+LRT_nwX+1)^2,npages])',Alpha);   
            temp=reshape(~T,[LRT_nl+LRT_nw+1,LRT_nl+LRT_nw+1]);
            temp=double(temp);
            SeedPoint=mean(mean(Matrix,3).*temp);
            if SeedPoint==0
                SeedPoint=mean(Ref);
            end
            %iteration (Gamma Confidence interval)
            if LRT_nl<RadiusRow && LRT_nw<RadiusCol
            MeanMatrix = meanmli(ll-LRT_nl:ll+LRT_nl,kk-LRT_nw:kk+LRT_nw);
            tempPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
            SeedPoint=mean(MeanMatrix(tempPoint));
            LRT_nl=LRT_nl+1;
            LRT_nw=LRT_nw+1;
            end
            %iteration (Gamma Confidence interval)
            MeanMatrix = meanmli(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol);
            SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
            SeedPoint(InitRow,InitCol)=true;
            %connection
            LL = bwlabel(SeedPoint); 
            SHP.PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);  
            num=num+1;
            if num == all_step * p;
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end
        end
         
    end
%SHPs map            
SHP.BroNum = sum(SHP.PixelInd,1);
SHP.BroNum = reshape(SHP.BroNum(:),[nlines,nwidths]);
SHP.BroNum = single((SHP.BroNum-1));          
SHP.CalWin = CalWin;            
t=toc;

figure;imagesc(SHP.BroNum);
colormap jet
ti=title ('Homogeneous Pixel Number');

end

