% function:detect_4_PAM(Y,A)
% Projec Name: Thl_1
% Engineer: Christos Trimas, Alexandros Michael


function [est_X] = detect_4_PAM(Y,A)
    X = [-3*A,-1*A,A,3*A];
    
    est_X = zeros(1,length(Y));
    
    for i=1:length(Y)
        %checking every distance for each element
        dist1 = norm(X(1)-Y(1,i));
        dist2 = norm(X(2)-Y(1,i));
        dist3 = norm(X(3)-Y(1,i));
        dist4 = norm(X(4)-Y(1,i));
        %finding minimum value
        minimum_val = min([dist1,dist2,dist3,dist4]);
        %checking shortest distance
        if(dist1==minimum_val)
            est_X(1,i) = X(1);
        elseif(dist2==minimum_val)
            est_X(1,i) = X(2);
        elseif(dist3==minimum_val)
            est_X(1,i) = X(3);
        elseif(dist4==minimum_val)
            est_X(1,i) = X(4);
        end
    end
end