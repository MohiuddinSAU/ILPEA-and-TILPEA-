function [bestFitness1,bestFitness_gobal1,bestSolution_gobal1,Fe2]=TLPEA(fhd,Pop_Number,Max_gen,VRmin,VRmax,D,func_num)
global initial_flag
initial_flag=0;
%% Initialize the first generation population
Vrmin=ones(1,D)*(VRmin+abs(VRmin));                           %Initialize low bound of variable 
Vrmax=ones(1,D)*(VRmax+abs(VRmin));                           %Initialize up bound of variable
Pop=repmat(Vrmin,Pop_Number,1)+(repmat(Vrmax,Pop_Number,1)-repmat(Vrmin,Pop_Number,1)).*rand(Pop_Number,D);   %Initialize the first generation population 


     Fitness_Pop=feval(fhd,Pop'-abs(VRmin),func_num);     %Calculate the function values of the first generation population 
     originPop(1)={Pop};
     [bestFitness1(1),~]=min(Fitness_Pop);                  %Store the fitness value of the first generation population  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:3      
%% Initialize parameters of DE
F=0.5;                   %Scaling factor or mutation factor,the value range is£¨0,1.2]
CR=0.9;                   %Crossover factor
jr=0.4;
%% Mutation operator of DE
%Three individuals are randomly selected from the population for generating a mutation individual
    for i=1:Pop_Number
        Pop_list=randperm(Pop_Number);                                 
        Pop_list(find(Pop_list==i))=[];                             
       
        %%

%          [bestfit,ind]=min(Fitness_Pop);       % best solution and location of best fitness
%           best = Pop(ind,:);
          
%           %% DE/rand/1
           Mutant(i,:)=Pop(Pop_list(1),:)+F*(Pop(Pop_list(2),:)-Pop(Pop_list(3),:));  
                          Mutant(i,:)= max(Vrmin,Mutant(i,:));
    end                               
   for i=1:Pop_Number
        for j=1:D
            if Mutant(i,j)<Vrmin(j)||Mutant(i,j)>Vrmax(j)        %Make sure that individuals are in the setting bounds.
                 Mutant(i,j)=Vrmin(j)+(Vrmax(j)-Vrmin(j))*rand;
            end
        end
   end
%% Crossover operator of DE
%Intersect the target individual and its corresponding mutation individual to generate trial individual of the target individual.
    for i=1:Pop_Number
        for j=1:D
            r=randperm(D);
            if rand<=CR||j==r(1)
                 trialPop(i,j)=Mutant(i,j);                     
            else
                 trialPop(i,j)=Pop(i,j);
            end
        end
    end
%% Selection operator of DE
%Compare the value of the target individuals and trial individuals for population evaluation. Individuals with better fitness value will be selected to enter the next iteration

Fitness_trial=feval(fhd,trialPop'-abs(VRmin),func_num);   %Evaluate the trial population
    for i=1:Pop_Number                                         
        if Fitness_trial(i)<Fitness_Pop(i)      
           Pop(i,:)=trialPop(i,:);
           Fitness_Pop(i)=Fitness_trial(i);
        end
    end
     
    originPop(I)={Pop};                              %Store the second, third generation population
    [bestFitness1(I),~]=min(Fitness_Pop);              %Store the fitness value of the best individual in each generation population
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for I=4:Max_gen
   % [bestFitness(I),ind]=min(Fitness_Pop);
%    t=0.01-(0.05/(Max_gen))*(i-(Max_gen));
   
    t=0.01-(3.99/Max_gen)*(I-Max_gen);         %The initial disturbance coefficient
    th=0.6;                                    % Difference Thrushold 
    Pop_list1=randperm(Pop_Number);
    Pop_list2=randperm(Pop_Number);
    Pop_list3=randperm(Pop_Number);
    


for k=1:Pop_Number
    for j=1:D
        X0=[originPop{1,1}(Pop_list1(k),j),originPop{1,2}(Pop_list2(k),j),originPop{1,3}(Pop_list3(k),j)];
 


          if abs(max(X0)-min(X0))<th &&(abs(X0(1)-X0(2))<th||abs(X0(1)-X0(3))<th||abs(X0(2)-X0(3))<th)                                               %GPEA concept                   
                 trialPop(k,j)=X0(1)+t*abs(max(X0)-min(X0))*(rand-0.5);                %Random perturbation model
           else
%                trialPop(k,j)=(X0(1)-X0(2)+X0(3));  % pi/2
               
%                  trialPop(k,j)=(-(1/4)*X0(1)+(1/2)*X0(2)+(3/4)*X0(3));

                 trialPop(k,j)=(X0(1)+0.73205*X0(2)-0.73205*X0(3));     
          end

        %               trialPop(k,j)=((4/3)*(X0(3))^(-0.5)+(1/3)*(X0(2))^(-0.5)-(2/3)*(X0(1))^(-0.5))^(-2);  %y=(ax+b)^-2
                    if  trialPop(k,j)<Vrmin(j)||trialPop(k,j)>Vrmax(j)     % Make sure that individuals are in the setting bounds.
                    trialPop(k,j)=Vrmin(j)+(Vrmax(j)-Vrmin(j))*rand; 
                    end
    end
end
    Fitness_trial=feval(fhd,trialPop'-abs(VRmin),func_num);     %Evaluate the trial population
     
     %% Topological LPEA
         [~,index]=min( Fitness_trial);         % 
         globalx =trialPop(index,:);
 
    for i=1:Pop_Number
        for j=1:D
            Xb=Vrmin(j)+Vrmax(j)-trialPop(i,j);
            [~,index2]=min([abs(globalx(j)-trialPop(i,j)),abs(globalx(j)-Xb)]);
            if index2==2
               trialPop(i,j)=Xb;
            end
        end
    end

     
     Pop=originPop{1,3};                              
     for i=1:Pop_Number
        if Fitness_trial(i)<Fitness_Pop(i)
           Pop(i,:)=trialPop(i,:);
           Fitness_Pop(i)=Fitness_trial(i);
        end
     end

     originPop(1,4)={Pop};                          %Generate the true population
     originPop=originPop(2:4);                      %Update the population chain to produce offspring
     [bestFitness1(I),~]=min(Fitness_Pop);           %Store the fitness value of the best individual in each generation population

end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output bestFitness_gobal and bestSolution_gobal
[bestFitness_gobal1,best_index]=min(Fitness_Pop);     %The fitness value of the gobal best individual 
bestSolution_gobal1=Pop(best_index,:);                %The gobal best individual                                
% bestFitness_gobal1=min(bestFitness);

end