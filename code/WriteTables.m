clear, clc, close all;
%% standard controller
[s1,s2,s3,s4,s5,s6,s7,s8,s9] = standard_controller(); %gives s1 - s9
%%
S = [s1, s2, s3, s4, s5, s6, s7, s8, s9];
T = struct2table(S);
T.Bayes = repmat("-", length(S), 1);
str_pair = {{'MPRK22adap',{'3000'}},... % choose method and insert the number of Bayes iterations 
    {'MPRK43Iadap',{'3000'}},...
    {'MPRK43IIadap',{'3000'}},...
    {'ROS2',{'3000'}}...
    };
for i = 1:length(str_pair)
    for j = 1:length(str_pair{i}{2})
        load(string(pwd) + "\Bayes_results" + str_pair{i}{1} + "_" + str_pair{i}{2}{j} + "_nov2025.mat"); % reading "results"
        %% Take the cheapest (pay no attention to evaluation time)
         res_aux = results.XAtMinObjective;
        %% Find the best candidate in terms of efficiency (Runtime vs Costs) from the top three (lowest costs)
         % res_time = results.ObjectiveEvaluationTimeTrace(unique(results.IndexOfMinimumTrace)); 
         % res_cost = results.ObjectiveTrace(unique(results.IndexOfMinimumTrace));
         % res_eff_ratio = res_time(end-2:end) .* (res_cost(end-2:end)).^2;
         % [~,idx] = min(res_eff_ratio);
         % id = find(results.ObjectiveEvaluationTimeTrace == res_time(end - 3 + idx)); 
         % res_aux = results.XTrace(id, :);
        %%
        res_aux.Bayes = str_pair{i}{1} + "\_" + str_pair{i}{2}{j};
        res((i-1)*length(str_pair{max(i-1,1)}{2}) +  j,:) = res_aux;
    end
end

update_tab = {'MPRK43Iadap','MPRK43IIadap', 'MPRK22adap','ROS2'}; % updates the table
for i = 1:length(update_tab)
    varName = update_tab{i};
    assignin('base',char(varName), [res; T])
    assignin('base',char(varName), costfun(update_tab{i}, eval(char(varName))) );
    save("Costs_" + update_tab{i} + "c1" + "_nov2025.mat", char(varName))
    tab = eval(char(varName));
    table2latex(tab(:,2:end), char(string(fileparts(pwd)) + '\' + varName + ".tex"));
end


% x.beta1=     1.3651 ; %MPRK22
% x.beta2=    0.12708;
% x.beta3=   -0.57219;
% x.alpha2=    -0.37702;
% x.kappa2=2;
%
% x2.beta1=     1.951 ; %MPRK22
% x2.beta2=   -0.086559;
% x2.beta3=  -0.37409;
% x2.alpha2=    -0.48842;
% x2.kappa2=2;
%
% x3.beta1 = 1.487387;
% x3.beta2 = -0.086559;
% x3.beta3 = -0.473185;
% x3.alpha2 = -0.265014 ;
% x3.kappa2 = 3;
%
% y.beta1=   1.2119  ; %MPRK43II
% y.beta2=             0.023048   ;
% y.beta3=      -0.12786 ;
% y.alpha2=      -0.56604  ;
% y.kappa2=2;
%
% y2.beta1=   2.2556  ; %MPRK43II
% y2.beta2=   -1.1991   ;
% y2.beta3=     -0.15024;
% y2.alpha2=      -2.2167  ;
% y2.kappa2=2;
%
%
% z.beta1= 1.1786; %MPRK43I
% z.beta2=  -0.47694  ;
% z.beta3= -0.027401 ;
% z.alpha2=  -0.71155  ;
% z.kappa2=4;
%
% z2.beta1= 1.7706; %MPRK43I
% z2.beta2=   -0.27744  ;
% z2.beta3=  -0.37701;
% z2.alpha2=  -0.95947   ;
% z2.kappa2=3;



%costfun(x2),costfun(z2),costfun(y2)

% costfun(MPRK22_1000)%, costfun(s2)