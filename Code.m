%% Differential Evolution Algorithm (different crossover factors)
clearvars; clc;
%% Determine constants
F = 0.6;                           %mutation factor
gen = 500;                         %number of generation required
kmin = 1; kmax = 50;               %Gain limits 
T1min = 0.05; T1max = 0.5;         %Time constant 1 limits
T2min = 0.1; T2max = 0.1;          %Time constant 2 limits
limits = [kmin, T1min, T2max; kmax, T1max, T2max]; %limits mat. of par.
npop = 50;                         %Number of population
run = 1;

%% DE main loop
% Initialization
k = kmin + (kmax-kmin)*rand(npop,1);
T1 = T1min + (T1max-T1min)*rand(npop,1);
T2 = T2min + (T2max-T2min)*rand(npop,1);
init_pop = [k,T1,T2];

for CR = 0.1:0.1:1.0
    old_pop = init_pop;
    for iter =1:1:gen
        V = zeros(npop,3);
        TV = zeros(npop,3);
        % Mutation
        for i = 1:1:npop
            k = randperm(npop,3);
            V(i,:) = old_pop(k(1),:)+F*(old_pop(k(2),:)-old_pop(k(3),:));
            for j = 1:1:3
                if V(i,j) < limits(1,j) || V(i,j) > limits(2,j)
                    V(i,j) = min(V(i,j),limits(2,j));
                    V(i,j) = max(V(i,j),limits(1,j));
                end
            end
        end
        % Crossover
        for i = 1:1:npop
            kk = rand;
            if kk > CR
                TV(i,:) = old_pop(i,:);
            else
                TV(i,:) = V(i,:);
            end
        end
        % Selection
        for i = 1:1:npop
            old_pop(i,4) = ff(old_pop(i,:));
            TV(i,4) = ff(TV(i,:));
            if old_pop(i,4)> TV(i,4)
                new_pop(i,:) = old_pop(i,:);
            else
                new_pop(i,:) = TV(i,:);
            end
        end
        % Evaluation
        YYY = sortrows(new_pop,-4);
        optimal_per_run(iter,:) = YYY(1,:);
        if iter == 1
            optimal_over(iter,:) = optimal_per_run(iter,:);
        elseif optimal_per_run(iter,4) >= optimal_over(iter-1,4)
            optimal_over(iter,:) = optimal_per_run(iter,:);
        else
            optimal_over(iter,:) = optimal_over(iter-1,:);
        end
        old_pop = new_pop(:,[1:3]);
    end
    final(:,run)=optimal_over(:,4);
    tabulated_results(run,:)=optimal_over(gen,:);
    run = run + 1;
end

%% Display final resutls
tabulated_results(11,:)=min(tabulated_results([1:10],:));
tabulated_results(12,:)=max(tabulated_results([1:10],:));
tabulated_results(13,:)=mean(tabulated_results([1:10],:));
tabulated_results(14,:)=std(tabulated_results([1:10],:));
fprintf('\n             K         T1        T2       OV\n');
fprintf('Run 1 '), disp(tabulated_results(1,:)); 
fprintf('Run 2 '), disp(tabulated_results(2,:));
fprintf('Run 3 '), disp(tabulated_results(3,:)); 
fprintf('Run 4 '), disp(tabulated_results(4,:));
fprintf('Run 5 '), disp(tabulated_results(5,:)); 
fprintf('Run 6 '), disp(tabulated_results(6,:));
fprintf('Run 7 '), disp(tabulated_results(7,:)); 
fprintf('Run 8 '), disp(tabulated_results(8,:));
fprintf('Run 9 '), disp(tabulated_results(9,:)); 
fprintf('Run 10'), disp(tabulated_results(10,:));
fprintf(' Min. '), disp(tabulated_results(11,:)); 
fprintf(' Max. '), disp(tabulated_results(12,:));
fprintf(' Mean '), disp(tabulated_results(13,:));
fprintf('  SD  '), disp(tabulated_results(14,:));

%plot final results
plot(final,'LineWidth',2), title('Differential Evolution Algorithm'), ...
    xlabel('Generation','FontSize',12), ylabel('Objective function'), ...
    legend('CR = 0.1','CR = 0.2','CR = 0.3','CR = 0.4','CR = 0.5','CR = 0.6','CR = 0.7','CR = 0.8', ...
    'CR = 0.9','CR = 1.0'), ylim([min(min(final))-0.01 0.49])...
    , grid on;