%%%%%%%%% SCRIPT ATUAL %%%%%%%%%

%close all; clear all; clc;

n=1e4;

pk=zeros(1,n); % probability vector
infok=zeros(1,n); % info vector
x = 1:n;

for i=1:n
	
	pk(i)=rand(1); % random probability vector
	
	infok(i)=-log2(pk(i)); % info vector

end


[maximum,index_max]=max(pk);
disp('Most probable Caractere:');
index_max
disp('probability:');
maximum


[minimum,index_min]=min(pk);
disp('Less probable Caractere:');
index_min
disp('probability:');
minimum

hold on

plot(pk,infok,'r')
xlabel('Probability');
ylabel('Information');
grid on;