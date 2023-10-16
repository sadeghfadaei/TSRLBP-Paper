function [table,newMax] = getmapping(samples,mappMode)

%       mappMode:
%       'no'   for normal LBP
%       'u2'   for uniform LBP
%       'ri'   for rotation-invariant LBP
%       'riu2' for uniform rotation-invariant LBP.

table = 0:2^samples-1;
newMax  = 0; %number of patterns in the resulting LBP code
index = 0;

%Uniform 2 mapping
if strcmp(mappMode,'no')
    newMax = 2^samples;
    table = 0:2^samples-1;
end
if strcmp(mappMode,'u2') 
    newMax = samples*(samples-1)+3;
    for i = 0:2^samples-1

        i_bin = dec2bin(i,samples);
        j_bin = circshift(i_bin',-1)';
        numt = sum(i_bin~=j_bin);                   

        if numt <= 2
            table(i+1) = index;
            index = index+1;
        else
            table(i+1) = newMax-1;
        end
    end
end

%Rotation invariant mapping
if strcmp(mappMode,'ri')
    tmpMap = zeros(2^samples,1)-1;
    for i = 0:2^samples-1
        rm = i;

        r_bin = dec2bin(i,samples);

        for j = 1:samples-1

            r = bin2dec(circshift(r_bin',-1*j)');  
            if r < rm
                rm = r;
            end
        end
        if tmpMap(rm+1) < 0
            tmpMap(rm+1) = newMax;
            newMax = newMax+1;
        end
        table(i+1) = tmpMap(rm+1);
    end
end

%Uniform & Rotation invariant mapping
if strcmp(mappMode,'riu2')
    newMax = samples+2;
    for i = 0:2^samples-1
        
        i_bin =  dec2bin(i,samples);
        j_bin = circshift(i_bin',-1)';
        numt = sum(i_bin~=j_bin);
  
        if numt <= 2
            table(i+1) = sum(bitget(i,1:samples));
        else
            table(i+1) = samples+1;
        end
    end
end

end % end of getmapping function