% Data structure MPRK22/MPRK43II [Beta1, Beta2, Beta3, tol, para,solvpara1, succstep, cnt_rej, err2exsol];
% Data struture MPRK43I [Beta1, Beta2, Beta3, tol, para,solvpara1,solvpara2, succstep, cnt_rej, err2exsol];
clc, clear, close all;
teststr = 'strat';
run(['test_' teststr])
if or(strcmp(teststr,'PR3'),strcmp(teststr,'PR4'))
    %A= fscanf(fopen('DATA_MPRK43IIadap.txt','r'),'%f %f %f %e %f %f %d %d %f',[9 Inf]);
    B= fscanf(fopen(append('DATA_MPRK43IIadap_',teststr,'.txt'),'r'),'%f %f %f %e %f %f %f %f %f',[9 Inf]);
    paraint =[0.2]; %%% set interval of parameter p
else
    %A= fscanf(fopen('DATA_MPRK43IIadap.txt','r'),'%f %f %f %e %f %f %f %f %f',[9 Inf]);
    B= fscanf(fopen(append('DATA_MPRK43IIadap_',teststr,'.txt'),'r'),'%f %f %f %e %f %f %f %f %f',[9 Inf]);
    paraint = [NaN];
end
%B =[A,B];
B=B';

% Beta1int = 0.1:dint:1;
% Beta2int = -0.4:dint:-0.05;
% Beta3int = 0:dint:0.1;
gammaint = [0.5, 0.75, 0.563];
gammaint=[0.750];


beta3 = 0;
%para = 0.9;
tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]; 
%tolint = [1e-5];
% z = B(:,7) + B(:,8);
% quant = quantile(z,0.1);

for gamma = gammaint
    for tol=tolint
        for para = paraint
            if isnan(para)
                C = B(B(:,3)==beta3 & B(:,4)==tol  & isnan(B(:,5)) & B(:,6)==gamma,:); %& ~isnan(B(:,7)) & ~isnan(B(:,8)) & ~isnan(B(:,9))
            else
                C = B(B(:,3)==beta3 & B(:,4)==tol  & B(:,5)==para & B(:,6)==gamma,:);
            end
            x = C(:,1);
            y = C(:,2);
            z1 = 3*(C(:,7) + C(:,8));
            figure
            scatter(x,y,100,z1,'square','filled');         
            c1 = colorbar;
            c1.Label.String = 'RHS evaluations';
            hold on 
            findnan7=find(isnan(C(:,7)));
            scatter(x(findnan7),y(findnan7),100,z1(findnan7),'square','filled',"r");
            findnan8=find(isnan(C(:,8)));
            scatter(x(findnan8),y(findnan8),100,z1(findnan8),'square','filled',"k");
            findnan9=find(isnan(C(:,9)));
            scatter(x(findnan9),y(findnan9),100,z1(findnan9),'square','filled',"w");
            xlim([0 1.1]);
            ylim([-0.5 0.3]);
            xlabel("beta1");
            ylabel("beta2");
            title(['beta3=',num2str(beta3),' tol=',num2str(tol,'%.e'),' p=',num2str(para),' gamma=',num2str(gamma)]);

            f= gca;
            %f.OuterPosition =[get( groot, 'Screensize' )];
            exportgraphics(f,['RHS_MPRK43IIadap_',teststr,'_beta3_',num2str(beta3),'_tol_',num2str(tol,'%.e'),'_p_',num2str(para),'_gamma_',num2str(gamma),'.eps'],'ContentType','vector');

            hold off
            %z2 = C(:,9)/norm(y0);
            z2 = C(:,9);
            figure
            scatter(x,y,100,z2,'square','filled');
            hold on 
            scatter(x(findnan7),y(findnan7),100,z2(findnan7),'square','filled',"r");
            scatter(x(findnan8),y(findnan8),100,z2(findnan8),'square','filled',"k");
            xlim([0 1.1]);
            ylim([-0.5 0.3]);
            c1 = colorbar;
            c1.Label.String = 'Error';
            set(gca,'ColorScale','log')
            xlabel("beta1");
            ylabel("beta2");
            title(['beta3=',num2str(beta3),' tol=',num2str(tol,'%.e'),' p=',num2str(para),' gamma=',num2str(gamma)]);

            f= gca;
            exportgraphics(f,['Error_MPRK43IIadap_',teststr,'_beta3_',num2str(beta3),'_tol_',num2str(tol,'%.e'),'_p_',num2str(para),'_gamma_',num2str(gamma),'.eps'],'ContentType','vector');
            hold off
        end
    end
end
figure
Markers= {'d','o','*','x','+','v','^','s','>','<'};
Markeropt ={'Markersize',8};
triple = {[0.1,-0.1,0],[0.17,-0.17,0],[0.6,-0.2,0],[0.7,-0.4,0],[0.41,-0.4,0]};
lgnd = cell(1,length(triple));
parameter = paraint(1); %%% value of p for the loglog-error plot
solvpara = gammaint(1);
for i=1:length(triple)
trip = triple{i};
trip1 = trip(1);
trip2 = trip(2);
trip3 = trip(3);
lgnd{i} = sprintf('[beta1,beta2,beta3]=[%.3f, %.3f,%.3f]',trip1,trip2,trip3);
if isnan(parameter)
    D=B(B(:,1)==trip1 & B(:,2)==trip2 & B(:,3)==trip3 & isnan(B(:,5)) & B(:,6)==solvpara,:);
else
    D=B(B(:,1)==trip1 & B(:,2)==trip2 & B(:,3)==trip3 & B(:,5)==parameter & B(:,6)==solvpara,:);
end
loglog(3*(D(:,7) + D(:,8)),D(:,9),Markers{i}, 'Linewidth',2, Markeropt{:})
hold on
end
title(['test: ',teststr,' p=',num2str(parameter),' gamma=',num2str(solvpara)]);
xlabel("Number of RHS evaluations");
ylabel("Error");
l = legend(lgnd,'location','ne');
l.Visible = 'on';
grid on
f= gca;
exportgraphics(f,['RHSvsError_MPRK43IIadap_',teststr,'_p_',num2str(para),'_gamma_',num2str(gamma),'.eps'],'ContentType','vector');
hold off
