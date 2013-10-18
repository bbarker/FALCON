function testLMoMA(FDBG)
% Check to make sure linearMOMAOneShot is working as expected.

%INPUT
% FDBG    Whether or not to turn on extra debugging info while 
%         running FALCON. This will increase the runtime.

REG = 0.1
objWeight = 0;

m = makeTestModel_InTriangleOut();

priorFlux = zeros(length(m.rxns), 1);
priorFlux(1:5) = 1;
rvec = zeros(length(m.rxns), 1);
rxnList = 1:5;


[solDel, solWT] = linearMOMAOneShot(m, priorFlux, objWeight, REG, 1, rvec, rxnList);
solDel.x'

%Branching seems to work:
priorFlux(2) = 0;
priorFlux(3) = -0.5;
priorFlux(4) = -0.5;
[solDel, solWT] = linearMOMAOneShot(m, priorFlux, objWeight, REG, 1, rvec, rxnList);
solDel.x'


