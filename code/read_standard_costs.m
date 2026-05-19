clear
str_pair = {'MPRK22adap','MPRK43Iadap','MPRK43IIadap','ROS2'};
    T = table();
for i = 1:length(str_pair)
    varName = str_pair{i};
    load(string(pwd)  + "\" + "Costs_" + str_pair{i} + "c1" + "_nov2025.mat"); % reading "results"
    tab = eval(char(varName));
    ids = strcmp(tab.Bayes,"-");

    T.Controller = "(" + tab.beta1(ids) + "," + tab.beta2(ids) + "," + tab.beta3(ids) + ","+ tab.alpha2(ids) + ","+ tab.kappa2(ids) + ")";
    T.Costs = tab.Costs(ids);
    [~,idx] = min(str2double(T.Costs));
    T.Costs(idx) = "\textbf{" + T.Costs(idx) + "}";
    % T = struct2table(T);
    switch varName
        case 'MPRK22adap'
            T.Properties.VariableNames(strcmp(T.Properties.VariableNames, 'Costs')) = {'MPRK22(1)'};
        case 'MPRK43Iadap'
            T.Properties.VariableNames(strcmp(T.Properties.VariableNames, 'Costs')) = {'MPRK43Iadap(0.5,075)'};
        case 'MPRK43IIadap'
            T.Properties.VariableNames(strcmp(T.Properties.VariableNames, 'Costs')) = {'MPRK43IIadap(0.563)'};
        case 'ROS2'
            T.Properties.VariableNames(strcmp(T.Properties.VariableNames, 'Costs')) = {'ROS2'};
    end

end
table2latex(T, char(string(fileparts(pwd)) + "\" + "standard_costs" + ".tex"));