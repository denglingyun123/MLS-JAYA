%_________________________________________________________________________%
%Jaya算法             %
%_________________________________________________________________________%
function [GBestX,GBestF,curve,Div]=JAYA(pop,Max_iter,lb,ub,dim,fobj)
% tic
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%种群初始化
X0=initialization_SSA(pop,dim,ub,lb);
X = X0;

%计算初始适应度值
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end

[fitness, index]= sort(fitness);%排序
GBestF = fitness(1);%全局最优适应度值

%按适应度排序,X(1,:)代表最优位置，X(end,:)代表最差位置
for i = 1:pop
    X(i,:) = X0(index(i),:);
end

GBestX = X(1,:);%全局最优位置
curve(1)=GBestF;

X_new = X;
t=2;
while t<=Max_iter
    
    Pbest = X(1,:);
    Pworst = X(end,:);
   for i = 1:pop
    %% Position updating stage
        X_new(i,:) = X(i,:) + rand(1,dim).*(Pbest - abs(X(i,:))) - rand(1,dim).*(Pworst - abs(X(i,:)));
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
   
   %更新位置
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
   
   curve(t)=GBestF;
   t=t+1;
   
    %排序更新
   [fitness, index]= sort(fitness);%排序
   X = X(index,:);
 
end

% Best_pos = GBestX;
% Best_score = curve(end);
% toc
% disp(['The runtime of JAYA: ',num2str(toc)]);
end



