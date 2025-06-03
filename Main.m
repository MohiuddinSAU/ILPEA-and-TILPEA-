clear all 
clc
close all

 
for i=1:30
    
func_num=i
sum=0;
D=30;
VRmin=-100;
VRmax=100;
Pop_Number=50;
 Max_Gen=7000;

% Max_Gen=floor(((10^4)*D)/Pop_Number);
fhd=str2func('cec14_func');

for run=1:51
            %% Search the best results using ILPEA (fhd,Pop_Number,Max_gen,VRmin,VRmax,D,func_num)
         [bestFitness1,bestFitness_gobal1,bestSolution_gobal1]=ILPEA(fhd,Pop_Number,Max_Gen,VRmin,VRmax,D,func_num);

%        
         result1(run)= bestFitness_gobal1-100*i;
         
%         %% Search the best results using TILPEA
        [bestFitness2,bestFitness_gobal2,bestSolution_gobal2]=TILPEA(fhd,Pop_Number,Max_Gen,VRmin,VRmax,D,func_num);
% 
         result2(run)= bestFitness_gobal2-100*i;

 end
%  %
       min_re11(i)=min(result1);
       max_re11(i)=max(result1);
       med_re11(i)=median(result1);
       mean_re11(i)=mean(result1);
       std_re11(i)=std(result1);
% %%
       min_re12(i)=min(result2);
       max_re12(i)=max(result2);
       med_re12(i)=median(result2);
       mean_re12(i)=mean(result2);
       std_re12(i)=std(result2);
       
       

               


 %% store Best Result
B1(i,:)=result1;
B2(i,:)=result2;

%% Store MEan REsult
Me1(i,:)=M1;
Me2(i,:)=M2;

% 
 end

%       
% 
B=[min_re11',mean_re11',std_re11',max_re11',med_re11'];
C=[min_re12',mean_re12',std_re12',max_re12',med_re12'];









