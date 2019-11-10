%function bits_to_2PAM
% Project Name: Thl_1
% Engineer: Christos Trimas, Alexandros Michael


function X = bits_to_2PAM(b)
 
%creating vector of zeros
    X=zeros(size(b));
    
    %matching 0 --> +1, 1 --> -1
    for k = 1:size(b)
        if(b(k)==0)
        X(k) = 1;
        
        elseif(b(k)==1)
        X(k) = -1;
        end;
    end;
end