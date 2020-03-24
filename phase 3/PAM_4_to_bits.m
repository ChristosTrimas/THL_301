% function:PAM_4_to_bits(X,A)
% Projec Name: Thl_1
% Engineer: Christos Trimas, Alexandros Michael


function [est_bit] = PAM_4_to_bits(X,A)
   Y = [-3*A, -1*A, A , 3*A];

counter=1;
%again create default space
est_bit = zeros(1,2*length(X));
for i=1:length(X)
    if(X(i) == Y(1))
        est_bit(counter) =0;
        est_bit(counter+1) = 0;
    elseif(X(i) == Y(2))
        est_bit(counter) =0;
        est_bit(counter+1) = 1;
    elseif(X(i) == Y(3))
        est_bit(counter) =1;
        est_bit(counter+1) = 1;
    elseif(X(i) == Y(4))
        est_bit(counter) =1;
        est_bit(counter+1) = 0;
    end
    %the counter increases by 2 for 
    counter = counter+2;
end
end