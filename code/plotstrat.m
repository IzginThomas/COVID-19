

for i=1:length(points)
 for j = 1:size(errors{i},1)
            test_errors = errors{i}(j, :);
            total_steps = steps{i}(j, :) + rejected{i}(j, :);
            subplot(1,size(errors{i},1), j);
            loglog(total_steps,test_errors,Markers{i}, 'Linewidth',2, Markeropt{:})
            if length(testset{j}) == 2
                teststr = testset{j}{1};
                testpara = testset{j}{2};
                title([teststr,' p=',testpara]);
            else
                title([tests{j}]);
            end
            xlabel("Number of total steps");
            ylabel("Error");
            grid on
            l = legend(lgnd,'location','ne','ItemHitFcn',@cb_legend);
            l.Visible = 'on';
            hold on
        end

end
%lgnd{4} = sprintf('(%.3f, %.3f, %.3f, %.3f, %.3f)',0.99,-0.068,0.054,0.176,2.000); %MPRK22c=5
