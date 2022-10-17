function Data=readData(imgpath,suffixname,nline,bkformat,machinefmt)
% see freadbkj

if nargin < 5
    machinefmt='b'; %GAMMA software, for example
end

if nargin < 4
    bkformat='float32'; %for *mli,*cc file
end

if isempty(strmatch(imgpath(end),filesep))
    imgpath=[imgpath,filesep];
end

tag_files = dir([imgpath,'*',suffixname]);
img_num = length(tag_files);
disp(['The number of the ', suffixname,' images:' num2str(img_num)]);

for ii=1:img_num
    tic;
    Data.datastack(:,:,ii)=single(freadbkj([imgpath,tag_files(ii).name],nline,bkformat,machinefmt));
    temp=regexp(tag_files(ii).name,'\d+','match');
    if length(temp)==1     %mli
        Data.filename(ii,1)=str2double(temp{1});
    elseif length(temp)==2 %intf
        Data.filename(ii,1)=str2double(temp{1});
        Data.filename(ii,2)=str2double(temp{2});
    else
        error('The format of file name should be: <yyyymmdd> or <yyyymmdd_yyyymmdd>.')
    end
    time=toc;
    fprintf('Reading Img %d / %d, time = %.0f sec\n',ii,img_num,time);      
end



