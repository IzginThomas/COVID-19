% Data structure MPRK22/MPRK43II [Beta1, Beta2, Beta3, tol, para,solvpara1, succstep, cnt_rej, err2exsol];
% Data struture MPRK43I [Beta1, Beta2, Beta3, tol, para,solvpara1,solvpara2, succstep, cnt_rej, err2exsol];
clc, clear, close all;
teststr = 'strat';
run(['test_' teststr])
if or(strcmp(teststr,'PR3'),strcmp(teststr,'PR4'))
    %A= fscanf(fopen('DATA_MPRK43IIadap.txt','r'),'%f %f %f %e %f %f %d %d %f',[9 Inf]);
    B= fscanf(fopen(append('DATA_MPRK43Iadap_',teststr,'.txt'),'r'),'%f %f %f %e %f %f %f %f %f %f',[10 Inf]);
    paraint = [0.2]; %%% set interval of parameter p
else
    %A= fscanf(fopen('DATA_MPRK43IIadap.txt','r'),'%f %f %f %e %f %f %f %f %f',[9 Inf]);
    B= fscanf(fopen(append('DATA_MPRK43Iadap_',teststr,'.txt'),'r'),'%f %f %f %e %f %f %f %f %f %f',[10 Inf]);
    paraint = [NaN];
end
%B =[A,B];
B=B';

% Beta1int = 0.1:dint:1;
% Beta2int = -0.4:dint:-0.05;
% Beta3int = 0:dint:0.1;
alphabetapairs = [0.5, 0.75; 
                1, 0.5];
alphabetapairs = [1,0.5];



beta3 = 0;
%para = 0.9;
tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]; 
%tolint = [1e-5];
% z = B(:,7) + B(:,8);
% quant = quantile(z,0.1);

for i = 1 : size(alphabetapairs,1)
    alpha = alphabetapairs(i,1);
    beta = alphabetapairs(i,2);
    for tol=tolint
        for para = paraint
            if isnan(para)
                C = B(B(:,3)==beta3 & B(:,4)==tol  & isnan(B(:,5)) & B(:,6)==alpha & B(:,7) == beta & B(:,1)+B(:,2)>0.05 ,:); %& ~isnan(B(:,7)) & ~isnan(B(:,8)) & ~isnan(B(:,9))
            else
                C = B(B(:,3)==beta3 & B(:,4)==tol  & B(:,5)==para & B(:,6)==alpha & B(:,7) == beta & B(:,1)+B(:,2)>0.05,:);
            end
            x = C(:,1);
            y = C(:,2);
            z1 = 3*(C(:,8) + C(:,9));
            figure
            scatter(x,y,100,z1,'square','filled');         
            c1 = colorbar;
            c1.Label.String = 'RHS evaluations';
            hold on 
            findnan8=find(isnan(C(:,8)));
            scatter(x(findnan8),y(findnan8),100,z1(findnan8),'square','filled',"r");
            findnan9=find(isnan(C(:,9)));
            scatter(x(findnan9),y(findnan9),100,z1(findnan9),'square','filled',"k");
            findnan10=find(isnan(C(:,10)));
            scatter(x(findnan10),y(findnan10),100,z1(findnan10),'square','filled',"w");
            xlim([0 1.1]);
            ylim([-0.5 0.3]);
            xlabel("beta1");
            ylabel("beta2");
            title(['beta3=',num2str(beta3),' tol=',num2str(tol,'%.e'),' p=',num2str(para),' alpha=',num2str(alpha),' beta=',num2str(beta)]);

            f= gca;
            %f.OuterPosition =[get( groot, 'Screensize' )];
            exportgraphics(f,['RHS_MPRK43IIadap_',teststr,'_beta3_',num2str(beta3),'_tol_',num2str(tol,'%.e'),'_p_',num2str(para),'_alpha_',num2str(alpha),'_beta_',num2str(beta),'.eps'],'ContentType','vector');

            hold off
            %z2 = C(:,9)/norm(y0);
            z2 = C(:,10);
            figure
            scatter(x,y,100,z2,'square','filled');
            hold on 
            scatter(x(findnan8),y(findnan8),100,z2(findnan8),'square','filled',"r");
            scatter(x(findnan9),y(findnan9),100,z2(findnan9),'square','filled',"k");
            xlim([0 1.1]);
            ylim([-0.5 0.3]);
            c1 = colorbar;
            c1.Label.String = 'Error';
            set(gca,'ColorScale','log')
            xlabel("beta1");
            ylabel("beta2");
            title(['beta3=',num2str(beta3),' tol=',num2str(tol,'%.e'),' p=',num2str(para),' alpha=',num2str(alpha),' beta=',num2str(beta)]);

            f= gca;
            exportgraphics(f,['Error_MPRK43IIadap_',teststr,'_beta3_',num2str(beta3),'_tol_',num2str(tol,'%.e'),'_p_',num2str(para),'_alpha_',num2str(alpha),'_beta_',num2str(beta),'.eps'],'ContentType','vector');
            hold off
        end
    end
end
figure
Markers= {'d','o','*','x','+','v','^','s','>','<'};
Markeropt ={'Markersize',8};
triple = {[0.1,-0.1,0],[0.1,0.15,0],[0.17,-0.17,0],[0.6,-0.2,0],[0.7,-0.4,0],[0.41,-0.4,0]};
lgnd = cell(1,length(triple));
parameter = paraint(1); %%% value of p for the loglog-error plot
solvpara = alphabetapairs(1,:);
for i=1:length(triple)
trip = triple{i};
trip1 = trip(1);
trip2 = trip(2);
trip3 = trip(3);
lgnd{i} = sprintf('[beta1,beta2,beta3]=[%.3f, %.3f,%.3f]',trip1,trip2,trip3);
if isnan(parameter)
    D=B(B(:,1)==trip1 & B(:,2)==trip2 & B(:,3)==trip3 & isnan(B(:,5)) & B(:,6)==solvpara(1) & B(:,7)==solvpara(2),:);
else
    D=B(B(:,1)==trip1 & B(:,2)==trip2 & B(:,3)==trip3 & B(:,5)==parameter & B(:,6)==solvpara(1) & B(:,7)==solvpara(2),:);
end
loglog(3*(D(:,8) + D(:,9)),D(:,10),Markers{i}, 'Linewidth',2, Markeropt{:})
hold on
end
title(['test: ',teststr,' p=',num2str(parameter),' alpha=',num2str(solvpara(1)),' beta=',num2str(solvpara(2))]);
xlabel("Number of RHS evaluations");
ylabel("Error");
l = legend(lgnd,'location','ne');
l.Visible = 'on';
grid on
f= gca;
exportgraphics(f,['RHSvsError_MPRK43Iadap_',teststr,'_p_',num2str(para),'_alpha_',num2str(solvpara(1)),'_beta_',num2str(solvpara(2)),'.eps'],'ContentType','vector');
hold off
