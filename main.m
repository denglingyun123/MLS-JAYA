clear all
clc
close all
SearchAgents_no=40; % Number of search agents
Max_iteration=1000; % Maximum number of iterations  
Function_name=1; %设定测试函数，1-29.其中大于11的要维度大于10，参考cec2017文档指定的维度
dim=30; %维度设定，维度可供选择范围[2,10,20,30,50,100]，其中Function_name>=11的最低维度设置为10.
lb=-100;%下边界
ub=100;%上边界
fobj = @(x) cec17_func(x',Function_name);


Max_test=30;
for i=1:Max_test
    disp(['The ',num2str(i),'-th experiment']);
    [Best_pos1(i,:),Best_score1(i),JAYA_curve(i,:)]=JAYA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); 
    [Best_pos2(i,:),Best_score2(i),MLSJAYA_curve(i,:)]=MLSJAYA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); 
    
end
% %结果对比
figure(1)
semilogy(mean(JAYA_curve),'color','[1,0.5,0]','linewidth',2.0,'Marker','o','MarkerIndices',1:50:length(mean(JAYA_curve)))
hold on
semilogy(mean(MLSJAYA_curve),'color','[0.62745,0.12549,0.94118]','linewidth',2.0,'Marker','o','MarkerIndices',1:50:length(mean(MLSJAYA_curve)))

xlabel('Iteration')
ylabel('Fitness')
title('Convergence curve of F2017_{1}')
legend('JAYA','MLS-JAYA')


disp('-------------------------------------------------')
display(['JAYA 30次实验最优适应度值(Best) : ', num2str(min(Best_score1))]);
display(['JAYA 30次实验最优解对应的平均适应度值(mean) : ', num2str(mean(Best_score1))]);
display(['JAYA 30次实验最差适应度值(worst) : ', num2str(max(Best_score1))]);
display(['JAYA 30次实验标准差（std） : ', num2str(std(Best_score1))]);


disp('-------------------------------------------------')
display(['MLS-JAYA 30次实验最优适应度值(Best) : ', num2str(min(Best_score2))]);
display(['MLS-JAYA 30次实验最优解对应的平均适应度值(mean) : ', num2str(mean(Best_score2))]);
display(['MLS-JAYA 30次实验最差适应度值(worst) : ', num2str(max(Best_score2))]);
display(['MLS-JAYA 30次实验标准差（std） : ', num2str(std(Best_score2))]);


