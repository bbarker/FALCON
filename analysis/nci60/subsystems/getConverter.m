%creates the converter input for subsystemsMakeHist
function convert = getConverter(rec2,EntrezIPI)

%INPUT
%rec2 is the human recon 2 model
%EntrezIPI is Entrez_and_IPI_unique.csv
%OUTPUT
%convert is a (1x2) cell array with Entrez ids in {1,1} and the
%corresponding IPI ids in {1,2}.
%Narayanan Sadagopan, September 2013

%get values
convert{1,1} = 0;
[ecs ENTtoIPI] = GenesByECOneFile(rec2, EntrezIPI);

%convert format
convert{1,1} = keys(ENTtoIPI);
convert{1,2} = values(ENTtoIPI);
