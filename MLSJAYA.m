%  Multi-learning JAYA algorithm (MLS-JAYA)                                                                 
%                                                                                                     
%  Developed in MATLAB R2018a                                                                 
%                                                                                                     
%  Author : Lingyun Deng, Sanyang Liu                                                         
%                                                                                                     
%         e-Mail: lingyundeng@stu.xidian.edu.cn                                                            
%                 syliu@xidian.edu.cn 
                                                                                                                                                                                                        
%  Main paper:                                                                                        
%  Lingyun Deng,Sanyang Liu. MLS-JAYA: A multi-learning JAYA algorithm for numerical optimization and efficient engineering design
%  Journal name= Cluster Computing, DOI: https://doi.org/10.1007/s10586-025-05873-1
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run MLSJAYA: [Best_pos,Best_score,Convergence_curve]=SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________

function [GBestX,GBestF,curve]=MLSJAYA(pop,Max_iter,lb,ub,dim,fobj)
% tic
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%Initialization
X0=initialization_SSA(pop,dim,ub,lb);
X = X0;

%Fitness evaluation
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end

[fitness, index]= sort(fitness);%Sorting
GBestF = fitness(1);%Current best fitness

%According to the fitness sorting,X(1,:) represents the current best
%position and X(end,:)is the current worst position
for i = 1:pop
    X(i,:) = X0(index(i),:);
end

GBestX = X(1,:);%Current best position
curve(1)=GBestF;
X_new = X;

% Initialize the probabilistic model
K=pop/2;  
alpha=.01; % Learning factor
mu=mean(X);
for j=1:dim
    std_devs(j)=std(X(:,j));
end


% Main loop
t=2;
while t<=Max_iter
    
%     GBestX = X(1,:);
    Psecond=X(2,:);
    Pworst = X(end,:);
    
%     % Centroid position of K promising solutions in the current population
%     centr_X=mean(X(1:K,:));
    
  %% Population-based incremeatal learning (PBIL)
    % Update the mean and standard deviation of Gaussian distribution
    for j=1:dim
        mu(j)=(1-alpha)*mu(j)+alpha*(GBestX(j)+Psecond(j)-Pworst(j));
        std_devs(j)=(1-alpha)*std_devs(j)+alpha*std(X(1:K,j));
    end
    
    % Generate pop candidate solutions according to the probabilistic model
    for j = 1:dim
        X_pro(:, j) = normrnd(mu(j), std_devs(j), [pop 1]);
    end
    
    %Bound control
    for j = 1:pop
       for a = 1: dim
           if(X_pro(j,a)>ub(a))
               X_pro(j,a) =ub(a);
           end
           if(X_pro(j,a)<lb(a))
               X_pro(j,a) =lb(a);
           end
       end
   end
    
    %Evaluate the above pop solutions
    for i=1:pop
        fitness_pro(i)=fobj(X_pro(i,:));
        if fitness_pro(i)<fitness(i)
           X(i,:) = X_pro(i,:);
           fitness(i) = fitness_pro(i);     
       end   
       if(fitness_pro(i) < GBestF)
           GBestF = fitness_pro(i);
           GBestX = X_pro(i,:);   
       end
    end
    
   %% A modified position updating stage of MLS-JAYA algorithm
   for i = 1:pop
       X_new(i,:) = X(i,:) + rand(1,dim).*(GBestX - X(i,:))- rand(1,dim).*(Pworst - X(i,:));
%      X_new(i,:) = X(i,:) + rand(1,dim).*(Pbest - abs(X(i,:))) - rand(1,dim).*(Pworst - abs(X(i,:)));
   end
    
   % Bound control
   for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub(a))
               X_new(j,a) =ub(a);
           end
           if(X_new(j,a)<lb(a))
               X_new(j,a) =lb(a);
           end
       end
   end
   
   %Greedy selection
   for j=1:pop
       fitness_new(j) = fobj(X_new(j,:));
       if fitness_new(j)<fitness(j)
           X(j,:) = X_new(j,:);
           fitness(j) = fitness_new(j);     
       end   
       if(fitness_new(j) < GBestF)
           GBestF = fitness_new(j);
           GBestX = X_new(j,:);   
       end
   end
   
   [fitness, index]= sort(fitness);%Sorting
   X = X(index,:);
   
   %% Orthogonal opposition based learning
   Centroid_X=mean(X);
   for i=pop/2
       X_temp1(i,:)=2*Centroid_X-X(i,:);
       % Bound control
       for j=1:dim
           if X_temp1(i,j)>ub(j)
               X_temp1(i,j)=Centroid_X(j)+rand*(ub(j)-Centroid_X(j));
           elseif X_temp1(i,j)<lb(j)
               X_temp1(i,j)=lb(j)+rand*(Centroid_X(j)-lb(j));
           end
       end
       
       value_space=[GBestX',X_temp1(i,:)'];
       X_temp(i,:)=OL_select_operation(dim,value_space,fobj);
       f_temp(i)=fobj(X_temp(i,:));
       if f_temp(i)<fitness(i)
           fitness(i)=f_temp(i);
           X(i,:)=X_temp(i,:);
       end
   end
   
   GBestF=fitness(1);
   GBestX=X(1,:);   
   curve(t)=GBestF;
   t=t+1;
   
end
% toc
% disp(['The runtime of MLS-JAYA: ',num2str(toc)]);

end



