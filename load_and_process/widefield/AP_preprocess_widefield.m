%% TESTING: replacing pipelineHere

im_path = 'G:\test_widefield';

im_files = dir(fullfile(im_path,'*.tif'));


curr_im_fn = fullfile(im_path,im_files(1).name);




%% (faster than imread)


tic; V = tiffreadVolume(curr_im_fn); toc;


%% Do datastore to keep all files accessible without loading together?
% https://uk.mathworks.com/matlabcentral/answers/350406-how-do-i-read-large-set-of-multi-page-tiff-files-into-tall-array-using-a-datastore

ds = fileDatastore({im_path},'ReadFcn', @tiffreadVolume);


ds = imageDatastore(im_path,'FileExtensions','.tif');

%% Load current lilrig options 

load('bluePurpleOps.mat');


% averages 7500 frames together??
ops.Nframes = size(V,3)*17;


%% Build up covariance matrix?

% Nicks loaded in ~1000 frames and averaged groups of ~15 frames (~500ms)

% grab one channel
V2 = single(V(:,:,1:2:6540));

% option 1: average groups of N frames? 
V3 = squeeze(nanmean(reshape(V2,size(V2,1),size(V2,2),15,[]),3));

% option 2: just get covariance of these pixels?
V2s = single(V2);
x = reshape(V2s,[],size(V2,3))'*reshape(V2s,[],size(V2,3));


%% Nick or someone's timestamp extractor
% PCO puts binary code into image in first 14 pixels, apparently
% (I can't find this online? just use his function for now)

% PCO timestamps in binary coded decimal (BCD) as:
% Pixels code: 
% 1:4          5 6 7 8 9 10 11 12    13    14
% framenum(x4) Y Y M D H M  S  10kus 100us us
% BCD converts to decimal as: a decimal number that codes for one byte, 
% with two nibbles that are one decimal each
% (e.g. 10dec = 16bcd: 16 -> 0001 0000 = "1,0" = 10)

pco_timestamps_bcd = squeeze(V(1,1:14,:));
pco_timestamps_decimal = cellfun(@(x) num2str(bin2dec(reshape(dec2bin(x,8),4,[])')'), ...
    num2cell(pco_timestamps_bcd),'uni',false);

pco_frame_num = cell2mat(pco_timestamps_decimal(1:4,:)')

% (join consecutive values for items, remove spaces, convert to numbers)
pco_frame_num = cellfun(@(x) str2num(strrep(x,' ','')), ...
    join(pco_timestamps_decimal(1:4,:)'));

pco_timestamp = cellfun(@(x) strrep(x,' ',''), ...
    join(pco_timestamps_decimal(5:13,:)'),'uni',false);

% THIS DOESN'T WORK IT REQUIRES 3 AND ONLY 3 MS WHICH MEANS
% I'LL HAVE TO CUT A DIGIT OFF OF THE SECOND TO LAST ONE
pco_timestamp = cellfun(@(x) ...
    datenum(strrep(x,' ',''),'yyyymmddHHMMSSFFF'), ...
    join(pco_timestamps_decimal(5:12,:)'));


pco_timestamp_ms = cellfun(@(x) strrep(x,' ',''), ...
    join(pco_timestamps_decimal(12:13,:)'),'uni',false);


[curr_frames, curr_timestamps_day] = timeFromPCOBinaryMulti(squeeze(V(1,1:14,:)));
curr_timestamps = curr_timestamps_day*24*3600; % convert from seconds to days


pco_frame_num = pco_timestamps_bcd(1:4,:)

% fixing here: binary coded decimal in digits, 9 = 9, 10 = 16 (0001 0000 =
% "1,0")
r = 32;
bin2dec(reshape(dec2bin(r,8),4,[])')


r = pco_timestamps_bcd(1:4,:)';
a = bin2dec(reshape(dec2bin(r,8),4,[])')

r = cellfun(@(x) bin2dec(reshape(dec2bin(x,8),4,[])'), ...
    num2cell(pco_timestamps_bcd(1:4,1:3)),'uni',false);




str2num(num2str(bin2dec(reshape(dec2bin(r,8),4,[])'))')

% from Nick's function timeFromPCOBinaryMulti

records = squeeze(V(1,1:14,:));



records = double(records);

nfr = size(records,2);

dic = zeros(14,nfr,2);

dic(:,:,2) = mod(records,16);
dic(:,:,1) = (records-dic(:,:,2))/16;

dicd = dic(:,:,1)*10+dic(:,:,2);

frameNums = dicd(1,:)*1000000 + dicd(2,:)*10000 + dicd(3,:)*100 + dicd(4,:);

year = dicd(5,:)*100+dicd(6,:);
mo = dicd(7,:);
day = dicd(8,:);
hr = dicd(9,:);
min = dicd(10,:);
sec = dicd(11,:) + dicd(12,:)/100 + dicd(13,:)/10000 + dicd(14,:)/1000000;
timestampsAsDatenum = datenum(year, mo, day, hr, min, sec);


% finally found the manual! https://www.pco.de/fileadmin/user_upload/pco-manuals/MA_PCOSDK_V127.pdf
% binary coded decimal (BCD): 
% 

% pixel 1-4: image counter
% pixel 5








