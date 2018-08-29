function cnByPosn = getCNfromSegs(chr,posn,cnfile,K,delBool)
%% FUNCTION getCNfromSegs(chr,posn,cnfile)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : March 7, 2011
% delBool = 1 if include deletion state (cn=1); 0 if do not include deletion

%% Copy number segment data
%1       51599   76205   2   0.3453
fid = fopen(cnfile);
cnData = textscan(fid,'%s%s%d%d%d%f%d%s','delimiter','\t','headerlines',1,'TreatAsEmpty','NA');
cnChr = cnData{2}; cnChr(strcmp(cnChr,'X')) = {'23'}; cnChr(strcmp(cnChr,'Y')) = {'24'}; cnChr(strcmp(cnChr,'MT')) = {'25'};
cnChr = str2double(cnChr);
cnStart = double(cnData{3});
cnStop = double(cnData{4});
cnCall = double(cnData{7});
cnLogR = double(cnData{6});

%% get copy number for each position
N = length(posn);
cnByPosn = ones(N,1) + 1; %initialize to default: 2 copy number

for i=1:N
   cnInd = find(cnChr==chr(i) & cnStart<=posn(i) & posn(i)<=cnStop);
   if ~(isempty(cnInd))
    cnByPosn(i) = getCN(cnCall(cnInd(1)),K,delBool);    
   end
end
end

function newCall = getCN(call,K,delBool)
    newCall = call;
    if (K==11)
        if (call==1 || call==7)
            newCall = 0;
        elseif ((call==2 || call==8) && delBool==0)
            newCall = 2;
        elseif ((call==2 || call==8) && delBool==1)
            newCall = 1;
        elseif (call==3)
            newCall = 2;
        elseif (call==4 || call==9)
            newCall = 3;
        elseif (call==5 || call==10)
            newCall = 4;
        elseif (call==6 || call==11)
            newCall = 5;
        end
    end
end
