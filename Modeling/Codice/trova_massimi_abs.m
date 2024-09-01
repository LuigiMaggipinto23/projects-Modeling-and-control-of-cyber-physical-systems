function massimi=trova_massimi_abs(x, j)
    
massimi=zeros(length(x),1);

for i=1:j
    m=max(abs(x));
    idx=0;
    for k=1:length(x)
        if abs(x(k))==m
            idx=k;
            x(k)=0;
            break;
        end
    end
    massimi(idx)=1;
end