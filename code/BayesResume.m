maxiteration = 1e4;

load('Bayes_results_2023_08_10_MPRK43I.mat')
newresultsMPRK43I = resume(results, 'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',maxiteration, 'AcquisitionFunctionName','expected-improvement-per-second-plus');
save("Bayes_resultsMPRK43Inew.mat")

load('Bayes_results_2023_08_09_MPRK43II.mat')
newresultsMPRK43II = resume(results, 'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',maxiteration, 'AcquisitionFunctionName','expected-improvement-per-second-plus');
save("Bayes_resultsMPRK43IInew.mat")

load('Bayes_results_2023_08_11_MPRK22.mat')
newresultsMPRK22 = resume(results, 'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',maxiteration, 'AcquisitionFunctionName','expected-improvement-per-second-plus');
save("Bayes_resultsMPRK22new.mat")