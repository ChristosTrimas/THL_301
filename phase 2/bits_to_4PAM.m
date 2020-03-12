%function bits_to_4PAM
% Project Name: Thl_1
% Engineer: Christos Trimas, Alexandros Michael

function X = bits_to_4PAM(b)
 
    X = zeros(size(b)); %creation of vector with zeros
    i = 1; %helping variable
    
    %matching 00 --> +3, 01 --> +1, 11 --> -1, 10 --> -3
    for k = 1:2:size(b)
        
        if(b(k)==0 && b(k+1)==0)
            X(i) = 3;
            
        elseif(b(k)==0 && b(k+1)==1)
            X(k) = 1;
        
        elseif(b(k)==1 && b(k+1)==1)
            X(k) = -1;
        
        elseif(b(k)==1 && b(k)==0)
            X(k) = -3;
            
        end;
    end;
end