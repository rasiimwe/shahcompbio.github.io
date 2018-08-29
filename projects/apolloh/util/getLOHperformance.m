function [precision,recall,accuracy] = getLOHperformance(Zpred,Ztruth)
%% FUNCTION [precision,recall,accuracy] = getLOHperformance(Gpred,Gtruth)
%
% author: Gavin Ha <gha@bccrc.ca>
%         Dept of Molecular Oncolgy/Centre for Translational and Applied Genomics
%         British Columbia Cancer Agency
%         University of British Columbia
% date  : April 12, 2011

LOHtruth = strcmp(Ztruth,'LOH') | strcmp(Ztruth,'DLOH');
LOHpred = strcmp(Zpred,'LOH') | strcmp(Zpred,'DLOH');
HETtruth = strcmp(Ztruth,'HET');
HETpred = strcmp(Zpred,'HET');
ASCNAtruth = strcmp(Ztruth,'ASCNA');
ASCNApred = strcmp(Zpred,'ASCNA');
TP = sum(LOHtruth & LOHpred);
TN = sum((HETtruth | ASCNAtruth) & (HETpred | ASCNApred));
FP = sum(~LOHtruth & LOHpred);
FN = sum(~(HETtruth | ASCNAtruth) & (HETpred | ASCNApred));

%precision
precision = TP / (TP + FP);
%recall
recall = TP / (TP + FN);
%accuracy
accuracy = (TP + TN) / (TP + TN + FP + FN);