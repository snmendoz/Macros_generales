function [concentrationMatrix, excRxnNames, timeVec, biomassVec] = dynamicpFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
% Performs dynamic FBA simulation using the static optimization approach
%
% USAGE:
%
%    [concentrationMatrix, excRxnNames, timeVec, biomassVec] = dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
%
% INPUTS:
%    model:                  COBRA model structure
%    substrateRxns:          List of exchange reaction names for substrates
%                            initially in the media that may change (e.g. not
%                            h2o or co2)
%    initConcentrations:     Initial concentrations of substrates (in the same
%                            structure as `substrateRxns`)
%    initBiomass:            Initial biomass (must be non zero)
%    timeStep:               Time step size
%    nSteps:                 Maximum number of time steps
%
% OPTIONAL INPUTS:
%    plotRxns:               Reactions to be plotted (Default = {'EX_glc(e)', 'EX_ac(e)', 'EX_for(e)'})
%    exclUptakeRxns:         List of uptake reactions whose substrate concentrations do not change
%                            (Default = {'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)'})
%
% OUTPUTS:
%    concentrationMatrix:    Matrix of extracellular metabolite concentrations
%    excRxnNames:            Names of exchange reactions for the EC metabolites
%    timeVec:                Vector of time points
%    biomassVec:             Vector of biomass values
%
% If no initial concentration is given for a substrate that has an open
% uptake in the model (i.e. `model.lb < 0`) the concentration is assumed to
% be high enough to not be limiting. If the uptake rate for a nutrient is
% calculated to exceed the maximum uptake rate for that nutrient specified
% in the model and the max uptake rate specified is > 0, the maximum uptake
% rate specified in the model is used instead of the calculated uptake
% rate.
%
% NOTE:
%
%    The dynamic FBA method implemented in this function is essentially
%    the same as the method described in
%    [`Varma, A., and B. O. Palsson. Appl. Environ. Microbiol. 60:3724 (1994)`].
%    This function does not implement the dynamic FBA using dynamic optimization approach
%    described in [`Mahadevan, R. et al. Biophys J, 83:1331-1340 (2003)`].
%
% .. Author: - Markus Herrgard 8/22/06

global WAITBAR_TYPE

if (nargin < 7)
    plotRxns = {'EX_glc(e)','EX_ac(e)','EX_for(e)'};
end

% Uptake reactions whose substrate concentrations do not change
if (nargin < 8)
    exclUptakeRxns = {'EX_co2(e)','EX_h2o(e)','EX_h(e)'};
end

% Find exchange rxns
excInd = findExcRxnsWithIDs(model);
% excInd = excInd & ~ismember(model.rxns,exclUptakeRxns);
excPos = excInd;
excRxnNames = model.rxns(excInd);
length(excRxnNames)
% Figure out if substrate reactions are correct
missingInd = find(~ismember(substrateRxns,excRxnNames));
if (~isempty(missingInd))
    for i = 1:length(missingInd)
        fprintf('%s\n',substrateRxns{missingInd(i)});
    end
    error('Invalid substrate uptake reaction!');
end

% Initialize concentrations
[~, substrateMatchInd] = ismember(substrateRxns,excRxnNames);
concentrations = zeros(length(excRxnNames),1);
concentrations(substrateMatchInd) = initConcentrations;

% Deal with reactions for which there are no initial concentrations
originalBound = -model.lb(excInd);
noInitConcentration = (concentrations == 0 & originalBound > 0);
% concentrations(noInitConcentration) = 1000;

biomass = initBiomass;

% Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

for i = 1:length(excRxnNames)
    %if is a substrate
    if model.lb(excPos(i))<0
        % if original lower bound is higher than the maximum calculated
        if -model.lb(excPos(i)) > uptakeBound(i)
            model.lb(excPos(i)) = -uptakeBound(i);
            disp('')
        end
    end
end

% Make sure bounds are not higher than what are specified in the model
% aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
% uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
% uptakeBound(originalBound<0) = originalBound(originalBound<0);
% model.lb(excInd) = -uptakeBound;

concentrationMatrix = sparse(concentrations);
biomassVec = biomass;
timeVec(1) = 0;

fprintf('Step number\tBiomass\n');
showprogress(0,'Dynamic FBA analysis in progress ...');
for stepNo = 1:nSteps
    % Run FBA
    %     sol = optimizeCbModel(model,'max','one');
    [~, ~, ~, solution] = pFBA_own(model, 'skipclass', 1);
    if isempty(solution)
        disp('')
        break;
    end
    sol.x = solution;
    sol.f = solution(find(model.c));
    if isempty(solution); sol.stat=0; else; sol.stat=1;end
    mu = sol.f;
    if (sol.stat ~= 1 || mu == 0)
        fprintf('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass);
        break;
    end
    uptakeFlux = sol.x(excInd);
    biomass = biomass*exp(mu*timeStep);
    %biomass = biomass*(1+mu*timeStep);
    biomassVec(end+1) = biomass;

    % Update concentrations
    concentrations = concentrations - uptakeFlux/mu*biomassVec(end-1)*(1-exp(mu*timeStep));
    %concentrations = concentrations + uptakeFlux*biomass*timeStep;
    concentrations(concentrations <= 0) = 0;
    concentrationMatrix(:,end+1) = sparse(concentrations);

        % Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
        % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    
    for i = 1:length(excRxnNames)
        %if is a substrate
        if model.lb(excPos(i))<0
            % if original lower bound is higher than the maximum calculated
            if -model.lb(excPos(i)) > uptakeBound(i)
                model.lb(excPos(i)) = -uptakeBound(i);
            end
        end
    end

%     % Figure out if the computed bounds were above the original bounds
%     aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
%     % Revert to original bounds if the rate was too high
%     uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
% %     uptakeBound(originalBound<0) = originalBound(originalBound<0);
%     uptakeBound(abs(uptakeBound) < 1e-9) = 0;
%     model.lb(excInd) = -uptakeBound;

    if WAITBAR_TYPE ~= 1
        fprintf('%d\t%f\n',stepNo,biomass);
    end
    showprogress(stepNo/nSteps);
    timeVec(stepNo+1) = stepNo*timeStep;
end

selNonZero = any(concentrationMatrix>0,2);
concentrationMatrix = concentrationMatrix(selNonZero,:);
excRxnNames = excRxnNames(selNonZero);
selPlot = ismember(excRxnNames,plotRxns);

% Plot concentrations as a function of time
% figure
% clf
% subplot(1,2,1);
% plot(timeVec,biomassVec);
% axis tight
% title('Biomass');
% subplot(1,2,2);
% plot(timeVec,concentrationMatrix(selPlot,:));
% axis tight
% legend(strrep(excRxnNames(selPlot),'EX_',''));
