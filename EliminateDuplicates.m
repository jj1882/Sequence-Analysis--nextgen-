function [out1,out2] = EliminateDuplicates (in1, in2)

if(length(in1)~=length(in2))
    error('in1 and in2 must have same number of elements');
end

n=length(in1);

[sortedin1, sortindex]=sort(in1);

sortedin2=in2(sortindex);

out1(1)=sortedin1(1);
out2(1)=sortedin2(1);

index=2;

for(a=2:1:n)
       
    if(strcmp(out1(index-1), sortedin1(a)))
        out2(index-1)=out2(index-1)+sortedin2(a);
    else
        out1(index)=sortedin1(a);
        out2(index)=sortedin2(a);
        
        index=index+1;
    end
     
end

