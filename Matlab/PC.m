
                                                   %% PC Algorithm %%
%% Mohammed Ombadi, 2020, (ombadi@lbl.gov), refer to Ombadi, M., Nguyen, P., Sorooshian, S., & Hsu, 
%% K. L. (2020). Evaluation of methods for causal discovery in hydrometeorological systems. Water Resources Research, 56(7), e2020WR027251. https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2020WR027251

% Inputs 
% Matrix of time series of system varaiables of size = v*n; (v= # of variables, n= # of observational time steps) 
X= importdata('F:\Moh\2019\Review of causality measures\Data\det_X.mat'); 
X= X(:,1:100);

tau_max= 3; % Maximum time lag 
a= 0.05; % Significance threshold 
qmax=1; % Maximum number of combinations 


v= size(X,1); % Number of variables in the system 
n= size(X,2); % Number of observations  
pmax= v*tau_max; % Maximum number of parents 

% Use the time series to create lagged variables 
Data= nan(size(X,1),(tau_max+1),(n-tau_max));

for i=1:v
    lag=0;
    for j=tau_max+1:-1:1
        Data(i,j,:)=X(i,(tau_max+1-lag):(end-lag));
        lag=lag+1;
    end 
end 

final_parents=nan(pmax,v);

for i=1:v %Outer loop to identify parents of each varaiable in the Data matrix except those in the 1st column
    
    column= size(Data,2) + 1 - ceil(i/size(X,1)); %Define the variable from the end of the matrix
    row= i - (ceil(i/size(X,1)) -1)*size(Data,1);
    variable= reshape(Data(row,column,:),numel(Data(row,column,:)),1); %time series of the variable 
    parents= (1:1:(column-1)*size(Data,1))'; % linear indices of parents
    I_min= inf(size(parents,1),size(parents,2));
    
    for p=0:pmax %1st Inner loop to define the number of parents to condition on
        
        if numel(parents)-1 >= p
            
            p_value=nan(size(parents,1),size(parents,2));
            
            for j=1:length(parents) %linear index of y-parent
                
                column= ceil(parents(j)/size(X,1)); 
                row= parents(j) - (ceil(parents(j)/size(X,1)) -1)*size(Data,1);
                y_parent= reshape(Data(row,column,:),numel(Data(row,column,:)),1); %time series of the hypothesized parent
                
                
                    parents_n= parents; 
                    parents_n(j)=[];  
                    all_sub_parent= combntns(parents_n,p);
                    
                    if length(all_sub_parent)==0
                        ss_count=1;
                    else
                        ss_count=size(all_sub_parent,1);
                    end
                   
                    for ss=1:ss_count
                    
                    sub_parent= all_sub_parent(ss,:);
                    [i j sub_parent]
                    z_parent= nan(size(y_parent,1),p);
                    for k=1:length(sub_parent)
                        column= ceil(sub_parent(k)/size(X,1)); 
                        row= sub_parent(k) - (ceil(sub_parent(k)/size(X,1)) -1)*size(Data,1);
                        z_parent(:,k)= reshape(Data(row,column,:),numel(Data(row,column,:)),1); %time series of the conditioned-on parent
                    end
                    
                   %%%%%%%%%%%%%%%%%
                    if numel(z_parent)>0
                        [I pv]=CMI(variable,y_parent,z_parent);
                        %[I pv]=partialcorr(variable,y_parent,z_parent);
                        I= abs(I);
                    else 
                        [I pv]=CMI(variable,y_parent);
                        %[I pv]=partialcorr(variable,y_parent);
                        I= abs(I);
                    end
                    %%%%%%%%%%%%%%%%%%%%
                    
                    
                    if pv <= a
                        p_value(j)=pv; % store the minimum value for the link (y_parent >> variable)
                        break
                    end
                    
                    
                    end
                    
                    p_value(j)=pv;
                    
                end
                
            
            
            % Remove insignificant links and thier I values 
            me= parents.*(p_value<=a);
            me(me==0)=[];
            parents=me;
                      
        end
    end
    
    % Return the parents of the variable i (defined from the end of matrix)
    final_parents(1:length(parents),i)= parents; %linear indices of parents defined from the beginning of matrix
    
end       