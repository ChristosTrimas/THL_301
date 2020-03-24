% function:bits_to_4_PAM(b,A) function
% Projec Name: Thl_1
% Engineer: Christos Trimas, Alexandros Michael


function [X] = bits_to_4_PAM(b,A)
    %a helping counter
    count=1;
    %possible outcomes
    Z = [-3*A, -1*A, A , 3*A];
    %zero filling from begin to b/2 since we have 4pam
    X=zeros(1,length(b)/2);
    %gray code means one bit distance(e.g.00->01->11->10)
    for i=1:2:length(b)
        if(b(i)==0 && b(i+1)==0)
            X(count) = Z(1);
    
        elseif(b(i)==0 && b(i+1)==1)
            X(count) = Z(2); 
        
        elseif(b(i)==1 && b(i+1)==1)
            X(count) = Z(3);
        
        elseif(b(i)==1 && b(i+1)==0)
            X(count) = Z(4);
        end
        count=count+1;
    end
end