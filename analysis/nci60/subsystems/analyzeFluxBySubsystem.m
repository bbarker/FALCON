function [dist] = analyzeFluxBySubsystem (rec2, vsolAllFluxes,subsys,norm)
%this function gives the absolute value of the flux running through
%each subsystem

%INPUT
%rec2 is the human recon 2 model
%vsolAllFluxes is the matrix with 7440 rows containing fluxes and each
%column represents a tissue
%subsys is cell array containing subsystems in each position
%norm is 1 if output should be normalized by number of reactions in the
%subystem. 0 for not normalized output. 
%OUTPUT
%dist is a matrix with the rows containing the sum of the absolute value 
%of each flux for a particular subsystem. The columns are for each tissue. 

%Narayanan Sadagopan, 10/14/2013

dist = zeros(length(subsys), length(vsolAllFluxes(1, :)));

if (norm)
    %selects tissue
    for a = 1:length(vsolAllFluxes(1,:)) 
        %selects subsystem
        for x = 1:length(subsys) 
            sum = 0;
            count = 0;
            %selects reaction
            for y = 1:length(rec2.rxns)
                if (strcmp(rec2.subSystems{y},subsys{x}))
                    sum = sum + abs(vsolAllFluxes(y,a));
                    count = count +1;
                end
            end
            dist(x, a) = sum/count;
        end
    end
else
    for a = 1:length(vsolAllFluxes(1,:))
        for x = 1:length(subsys)
            sum = 0;
            for y = 1:length(rec2.rxns)
                if (strcmp(rec2.subSystems{y},subsys{x}))
                    sum = sum + abs(vsolAllFluxes(y,a));
                end
            end
            dist(x, a) = sum;
        end
    end
end


