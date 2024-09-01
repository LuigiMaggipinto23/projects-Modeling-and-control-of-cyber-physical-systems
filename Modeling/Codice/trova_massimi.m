function massimi=trova_massimi(x, j)
    
massimi=zeros(length(x),1);

for i=1:j
    m=max(x);
    idx=0;
    for k=1:length(x)
        if x(k)==m
            idx=k;
            x(k)=0;
            break;
        end
    end
    massimi(idx)=1;
end