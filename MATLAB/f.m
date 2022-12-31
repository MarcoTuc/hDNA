classdef f
    methods (Static)
%%%%%%%%%%%%
        function SD = explore(model,variant,param1,param2,outs,runtime)
            param1_array = logspace(param1.min,param1.max,param1.num);
            param2_array = logspace(param2.min,param2.max,param2.num);
            trials = [];
            for i = 1:length(param1_array)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,param2_array']; %#ok<*AGROW>
            end
            model_function = createSimFunction(model,{param1.name, param2.name},outs,[],variant,'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(trials,runtime);
        end
%%%%%%%%%%%%%
        function trials = maketrialsmatrix(min1,max1,num1,min2,max2,num2)
            param1_array = logspace(min1,max1,num1);
            param2_array = logspace(min2,max2,num2);
            trials = [];
            for i = 1:length(param1_array)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,param2_array']; %#ok<*AGROW>
            end
        end
%%%%%%%%%%%%%
        function SD = explore_duplexequalized(model,variant,param1,param2,outs,runtime)
            param1_array = logspace(param1.min,param1.max,param1.num);
            param2_array = logspace(param2.min,param2.max,param2.num);
            trials = [];
            for i = 1:length(param1_array)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,exp_vector,param2_array']; %#ok<*AGROW>
            end
            model_function = createSimFunction(model,{param1.name,'kf_2d_duplex', param2.name},outs,[],variant,'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(trials,runtime);
        end
%%%%%%%%%%%%%
        function SD = explore_duplexbulkequalized(model,variant,param1,param2,outs,runtime)
            param1_array = logspace(param1.min,param1.max,param1.num);
            param2_array = logspace(param2.min,param2.max,param2.num);
            trials = [];
            for i = 1:length(param1_array)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,exp_vector,param2_array',param2_array']; %#ok<*AGROW>
            end
            model_function = createSimFunction(model,{param1.name,'kf_2d_duplex', param2.name, 'kf_3d'},outs,[],variant,'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(trials,runtime);
        end
%%%%%%%%%%%%%
function SD = explore_check(model,variant,param1,param2,outs,runtime)
            param1_array = logspace(param1.min,param1.max,param1.num);
            param2_array = logspace(param2.min,param2.max,param2.num);
            trials = [];
            for i = 1:length(param1_array)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,exp_vector,exp_vector,exp_vector,param2_array']; %#ok<*AGROW>
            end
            model_function = createSimFunction(model,{'kf_TR1','kf_TL1','kf_duplex','kf_O1','kb_duplex'},outs,[],variant,'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(trials,runtime);
        end
%%%%%%%%%%%%%
        function SD = explore_experimental(model,params,outs,runtime)
            param_array = [];
            for i = 1:length(params)
            new_array = logspace(params(i).min,params(i).max,params(i).num);
            param_array(i) = [param_array; new_array];
            end
            trials = [];
            for i = 1:size(param_array,1)
                exp_vector = param1_array(i)*ones(size(param2_array))';
                trials = [trials;exp_vector,param2_array']; %#ok<*AGROW>
            end
            model_function = createSimFunction(model,{param1.name, param2.name},outs,[],'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(trials,runtime);
        end
%%%%%%%%%%%%%
        function SD = pointexplore(model,params_names,params_values,outs,runtime)
            model_function = createSimFunction(model,params_names,outs,[],'UseParallel',true,'AutoAccelerate',true);
            SD = model_function(params_values,runtime);
            plot(SD.time,SD.data)
        end
%%%%%%%%%%%%
        function maxsurf(data, param1, param2)
            a = [];
            for i = 1:length(data)
                [~, t] = max(data(i).data(:,1));
                a = [a; data(i).time(t)];
            end
            b = reshape(a,param2.num,[]);
            X = logspace(param1.min,param1.max,param1.num);
            Y = logspace(param2.min,param2.max,param2.num);
            Z = b;
            surf(X,Y,Z)
            xlabel(strcat('X = ',param1.name),'Interpreter','none');
            ylabel(strcat('Y = ',param2.name),'Interpreter','none');  
            zlabel('Z = peak time')
        end
%%%%%%%%%%%%
        function logmaxsurf(data, param1, param2)
            a = [];
            for i = 1:length(data)
                [~, t] = max(data(i).data(:,1));
                a = [a; data(i).time(t)];
            end
            b = reshape(a,param2.num,[]);
            X = logspace(param1.min,param1.max,param1.num);
            Y = logspace(param2.min,param2.max,param2.num);
            Z = b;
            surf(X,Y,Z)
            h = gca;
            set(h,'ColorScale','log')
            set(h,'zscale','log');
            xlabel(strcat('X = ',param1.name),'Interpreter','none');
            ylabel(strcat('Y = ',param2.name),'Interpreter','none');  
            zlabel('Z = peak time')
        end
%%%%%%%%%%%%
        function logallmaxsurf(data, param1, param2)
            a = [];
            for i = 1:length(data)
                [~, t] = max(data(i).data(:,1));
                a = [a; data(i).time(t)];
            end
            b = reshape(a,param2.num,[]);
            X = logspace(param1.min,param1.max,param1.num);
            Y = logspace(param2.min,param2.max,param2.num);
            Z = b;
            surf(X,Y,Z)
            h = gca;
            set(h,'ColorScale','log')
            set(h,'xscale','log');
            set(h,'yscale','log');
            set(h,'zscale','log');
            xlabel(strcat('X = ',param1.name),'Interpreter','none');
            ylabel(strcat('Y = ',param2.name),'Interpreter','none');  
            zlabel('Z = peak time')
        end
%%%%%%%%%%%%
        function maxsurf_mod(data, index, param1, param2, log)
            a = [];
            for i = 1:length(data)
                [~, t] = max(data(i).data(:,index));
                a = [a; data(i).time(t)];
            end
            b = reshape(a,param2.num,[]);
            x = logspace(param1.min,param1.max,param1.num);
            y = logspace(param2.min,param2.max,param2.num);
            [X,Y] = meshgrid(x,y);
            Z = b;
            surf(X,Y,Z)
            h = gca;
            if log == true
            set(h,'ColorScale','log')
            set(h,'xscale','log');
            set(h,'yscale','log');
            set(h,'zscale','log');
            end
            xlabel(strcat('X = ',param1.name),'Interpreter','none');
            ylabel(strcat('Y = ',param2.name),'Interpreter','none');  
            zlabel('Z = peak time')
        end
%%%%%%%%%%%%
        function peaktime = peakplot(data,pos)
            [m, t] = max(data.data(:,pos));
            peaktime = data.time(t);
            hold on
            plot(data.time(t),m,'r+','LineWidth',2)
            hold off 
        end
%%%%%%%%%%%%
        function peaktime = sbiopeakdata(data,pos)
            [~, t] = max(data.data(:,pos));
            peaktime = data.time(t);
        end
%%%%%%%%%%%%
        function peaktime = peakdataode(data,timevector,pos)
            [~, t] = max(data(:,pos));
            peaktime = timevector(t);
        end
%%%%%%%%%%%%
        function variantlegend(variantcontent,fontsize)
        legend = [];
        for i = 1:length(variantcontent)
        str = cellstr(string(variantcontent{:,i}(2)));
        num = sprintf('%.3e',string(variantcontent{:,i}(4)));
        legend = [legend, join([str,num],' = ')];
        end
        annotation('textbox',[0 .9 .1 .1],'String',join(legend,newline,1),'Interpreter','none','FontSize',fontsize)
        end
%%%%%%%%%%%%
        function setmodel_ode15s(model,abstol,reltol,runtime,timeunits)
            conf = model.config;
            set(conf, 'solvertype','ode15s');
            set(conf.solveroptions, 'absolutetolerance',abstol);
            set(conf.solveroptions, 'relativetolerance',reltol);
            set(conf, 'timeunits',timeunits);
            set(conf, 'stoptime',runtime);
            sbioaccelerate(model);
        end
     end
 end