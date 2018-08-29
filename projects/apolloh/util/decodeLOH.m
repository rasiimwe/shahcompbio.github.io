function [Z,CN] = decodeLOH(G,mirror)
%% FUNCTION Z = decodeLOH(G,mirror)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 12, 2011

T = length(G);
Z = cell(T,1);
CN = zeros(T,1);

if (mirror==1)
    DLOH = G==1;
    NLOH = G==2;
    LOH = G==4 | G==6 | G==9;
    HET = G==3 | G==5 | G==8 | G==11;
    ASCNA = G==7 | G==10;
else
    DLOH = G==1;
    NLOH = G==2 | G==4;
    LOH = G==5 | G==8 | G==9 | G==13 | G==14 | G==19;
    HET = G==3 | G==6 | G==7;
    ASCNA = G==10 | G==12 | G==15 | G==18;
    BCNA = G==11 | G==16 | G==17;
end
HOMD = G==0;

Z(HOMD) = {'HOMD'};
Z(DLOH) = {'DLOH'};
Z(NLOH) = {'NLOH'};
Z(LOH) = {'ALOH'};
Z(HET) = {'HET'};
Z(ASCNA) = {'ASCNA'};
Z(BCNA) = {'BCNA'};

CN(HOMD) = 0;
CN(DLOH) = 1;
CN(G==2 | G==3 | G==4) = 2;
CN(G>=5 & G<=8) = 3;
CN(G>=9 & G<=13) = 4;
CN(G>=14 & G<=19) = 5;

