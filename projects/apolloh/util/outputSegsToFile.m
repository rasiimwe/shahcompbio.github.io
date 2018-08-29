function outputSegsToFile(segs,chr,posn,ratio,mirror,segoutfile)
%% FUNCTION outputSegsToFile(segs,chr,posn,ratio,segoutfile)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 15, 2011

numsegs=0;
chrs = unique(chr);
for i=1:max(chrs)
    [ms,ns] = size(segs{i});
    numsegs = numsegs + ms;
end
flatsegs = zeros(numsegs,5);
index=1;
for i=1:22
    [ms,ns] = size(segs{i});
    I = index:(index+ms-1);
    cInd = find(chr==i);
    flatsegs(I,1) = i;
    flatsegs(I,2) = posn(cInd(segs{i}(:,1)));
    flatsegs(I,3) = posn(cInd(segs{i}(:,2)));
    flatsegs(I,5) = segs{i}(:,3);
    %flatsegs(I(1),2) = posn(cInd(segs{i}(1,1))); % begin of chr
    %flatsegs(I(end),3) = posn(cInd(segs{i}(end,2))); % end of chr
    for j=1:length(I)
        flatsegs(I(j),4) = nanmedian(ratio((cInd(segs{i}(j,1))):(cInd(segs{i}(j,2)))));        
    end
    index = index+ms;
end
segCalls = decodeLOH(flatsegs(:,5)+1,mirror);
% write out the segments
fid = fopen(segoutfile,'wt');
for i = 1:numsegs
    fprintf(fid,'%d\t%d\t%d\t%1.4f\t%d\t%s\n',flatsegs(i,:),segCalls{i});
end
fclose(fid);
