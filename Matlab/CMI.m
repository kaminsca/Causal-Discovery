                                           %% Conditional Mutual Information (CMI) %%
%% Mohammed Ombadi, 2020, (ombadi@lbl.gov), refer to Ombadi, M., Nguyen, P., Sorooshian, S., & Hsu, 
%% K. L. (2020). Evaluation of methods for causal discovery in hydrometeorological systems. Water Resources Research, 56(7), e2020WR027251. https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2020WR027251
                                           
% This function returns the value of CMI for the General case CMI(X,Y,Z)
% using K nearest neighbors algorithm 
% CMI(X,Y,Z) is to test the independence of X from Y given Z 
% The algorithm also returns p-value (p) and optionally a range of CMI
% values for a range of k (k_range)

function [I_data,p,k_range]= CMI(X,Y,Z)

% X,Y and Z could be a column vectors (n*1) or matrix (n*m)

n= size(X,1); % length of samples
mx= size(X,2);
my= size(Y,2);

%Data_nn=nan(n,3,5);

% Free parameters 
k= ceil(0.2*n); %K_CMI (this is to define the k nearest neighbor)
k_perm= 5;
B= 500; % Number of surrogates to estimate the null distribution
%k_max= ceil(0.5*n); %maximum value of k_CMI for the algorithm to use 

% Initialization 
I_perm= nan(B,1);
%k_range= nan(k_max,1);



switch nargin
    
    %%
    case 3
    
mz= size(Z,2); 
        
                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                                                      % Pre-processing %
    
                                                      
    % Add small-amplitude random noise to break ties 
    % Debatabe question: Do we apply this for zero values as well (i.e. in
    % variables like precipitation and overland flow)
     for i=1:n
         for j=1:mx
               X(i,j)= X(i,j)+ normrnd(0,std(X(:,j))/1000); 
         end
         
         for j=1:my
              Y(i,j)= Y(i,j)+ normrnd(0,std(Y(:,j))/1000);
         end
    
         for j=1:mz
             Z(i,j)= Z(i,j)+ normrnd(0,std(Z(:,j))/1000);
         end
     end 


     Data_orig=[X,Y,Z];

     parfor perm=1:B+1
         if perm~=1
             Data= cmipermute(Data_orig,k_perm,mx,my,mz);
             %Data_nn(:,:,perm)=Data; %To record the permuted data 
         else
             Data=Data_orig;
         end  
         
      % Rank-transfer the data (From largest to smallest)
     
     X_rank= nan(n,mx);
     for j=1:mx
     X_sorted= sort(Data(:,j),'descend');
     [a X_rank(:,j)]= ismember(Data(:,j),X_sorted);
     end

     Y_rank= nan(n,my);
     for j=mx+1:mx+my
     Y_sorted= sort(Data(:,j),'descend');
     [a Y_rank(:,j-mx)]= ismember(Data(:,j),Y_sorted);
     end

     Z_rank= nan(n,mz);
     for j=(mx+my+1):(mx+my+mz)
     Z_sorted= sort(Data(:,j),'descend');
     [a Z_rank(:,j-mx-my)]= ismember(Data(:,j),Z_sorted);
     end

     X= X_rank;
     Y= Y_rank;
     Z= Z_rank;
     Data= nan(n,mx+my+mz);    
     Data= [X,Y,Z];        
         
                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                                                      % CMI Calculations%
     
        %for k=1:k_max
        summation=0;

        for i=1:n
            dist= nan(n,1);
            dist= max(abs(Data-Data(i,:)),[],2); %maximum norm of each point from the ith point 
            dist(dist==0)=nan; %replace the distance of the point itself by nan 
            dist_sorted= sort(dist);
            epsilon= dist_sorted(k); %distance to the k nearest neighbor in the joint space 
    
            % Calculate Kz(i) 
            dist= nan(n,1);
            dist= max(abs(Data(:,(mx+my+1):end)-Data(i,(mx+my+1):end)),[],2); %maximum norm of each point from the ith point in the subspace Z 
            dist(dist==0)=nan; %replace the distance of the point itself by nan 
            kz= numel(dist(dist<epsilon));
    
            % Calculate Kxz(i) 
            dist= nan(n,1);
            dist= max(abs(Data(:,[1:mx (mx+my+1):end])-Data(i,[1:mx (mx+my+1):end])),[],2); %maximum norm of each point from the ith point in the subspace X*Z 
            dist(dist==0)=nan; %replace the distance of the point itself by nan 
            kxz= numel(dist(dist<epsilon));
    
            % Calculate Kyz(i) 
            dist= nan(n,1);
            dist= max(abs(Data(:,(mx+1):end)-Data(i,(mx+1):end)),[],2); %maximum norm of each point from the ith point in the subspace Y*Z 
            dist(dist==0)=nan; %replace the distance of the point itself by nan 
            kyz= numel(dist(dist<epsilon));
    
            sum_i = psi(kz)-psi(kxz)-psi(kyz);
            summation= summation + sum_i;

       end 

       I = psi(k) + (1/n)*summation;

       I_perm(perm)=I;


     %end
     end

I_data= I_perm(1);     
I_perm= I_perm(2:end);    
     
     
    case 2

                                       %% Mutual Information (Not conditional) 

                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                                                      % Pre-processing %                                       
mz=0;          
                                                      
for i=1:n
    for j=1:mx
           %if X(i)~=0
               X(i,j)= X(i,j)+ normrnd(0,var(X(:,j))/1000); 
           %end 
    end
    
    for j=1:my
           %if X(i)~=0
               Y(i,j)= Y(i,j)+ normrnd(0,var(Y(:,j))/1000); 
           %end 
    end
        
end 

Data_orig=[X,Y];

parfor perm=1:B+1

if perm~=1
   Data= cmipermute(Data_orig,k_perm,mx,my,mz);    
   Data_nn(:,:,perm)=Data;
else
    Data=Data_orig;
end


% Rank-transfer the data (From largest to smallest)

X_sorted= sort(Data(:,1),'descend');
[a X_rank]= ismember(Data(:,1),X_sorted);

Y_sorted= sort(Data(:,2),'descend');
[a Y_rank]= ismember(Data(:,2),Y_sorted);


X= X_rank;
Y= Y_rank;
Data= nan(n,size(X,2)+size(Y,2));
Data= [X,Y];

                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
                                                      % CMI Calculations%
%for k=1:k_max
summation=0;

for i=1:n
    dist= nan(n,1);
    dist= max(abs(Data-Data(i,:)),[],2); %maximum norm of each point from the ith point 
    dist(dist==0)=nan; %replace the distance of the point itself by nan 
    dist_sorted= sort(dist);
    epsilon= dist_sorted(k); %distance to the k nearest neighbor in the joint space 
    
    % Calculate Kx(i) 
    dist= nan(n,1);
    dist= abs(Data(:,1)-Data(i,1)); %maximum norm of each point from the ith point in the subspace Z 
    dist(dist==0)=nan; %replace the distance of the point itself by nan 
    kx= numel(dist(dist<epsilon));
    
    % Calculate Ky(i) 
    dist= nan(n,1);
    dist= abs(Data(:,2)-Data(i,2)); %maximum norm of each point from the ith point in the subspace Z 
    dist(dist==0)=nan; %replace the distance of the point itself by nan 
    ky= numel(dist(dist<epsilon));
        
    sum_i = psi(kx+1)+psi(ky+1);
    summation= summation + sum_i;

end 

I = psi(k) - (1/n)*summation + psi(n);
%k_range(k)=I;
I_perm(perm)=I;
%k_range(k)= I;


end

I_data= I_perm(1);
I_perm= I_perm(2:end);

end


%end

% Calculate the p_value 

p= (1/B)*sum(double(I_perm>=I_data)); %The null hypothesis (H0) is that X is independent of Y 


end 



                                                     %% Permutation  %%

function[Data_perm]=cmipermute(Data,k_perm,mx,my,mz)

if mz>0
    
n= size(Data,1); %length of data 
U= []; %used indices 
Data_perm= nan(size(Data,1),size(Data,2));

for i=1:n 
    
  dist= nan(n,1);
  dist= max(abs(Data(:,(mx+my+1):end)-Data(i,(mx+my+1):end)),[],2); %maximum norm of each point from the ith point in the subspace Z  
  dist_sorted= sort(dist);
  epsilon= dist_sorted(k_perm); %distance to the k_perm nearest neighbor in the z space 
    
  % Identify nearest neighbors in the subspace Z 
  N= find(dist<=epsilon); %neighbors 
  
  % Shuffle the neighbors 
  N_n = N(randperm(length(N)));
  N= N_n;
  
  j= N(1);  % ??????
  m=0;
  while ismember(j,U) && m<(k_perm-1)
      m=m+1;
      j=N(m);
  end 
  Data_perm(i,1)= Data(j,1);
  U= [U j];
  
end
  
Data_perm(:,2:end)= Data(:,2:end); 

else
    
 n= size(Data,1); %length of data 
 X= Data(:,1);
 Data_perm= nan(size(Data,1),size(Data,2));
 Data_perm(:,1)= X(randperm(length(X))); %permute X to destroy dependence with Y
 Data_perm(:,2)=Data(:,2); %keep Y the same  

end

end


3, 2, 1, 4