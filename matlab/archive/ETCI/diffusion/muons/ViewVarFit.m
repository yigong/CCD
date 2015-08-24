function ViewVarFit(M2)
%function ViewVarFit(M2)

for i=1:55; 
    if ~isempty(M2{i}); 
        if isfield(M2{i},'DC'); 
            plot(M2{i}.DT,M2{i}.var,'.k'); hold on; 
            if isfield(M2{i},'yint')
                plot([0;2.5e-8],[M2{i}.yint;M2{i}.yint+2.5e-8*M2{i}.DC*2],'-r');
            else
                plot([0;2.5e-8],[19.9;19.9+2.5e-8*M2{i}.DC*2],'-r'); 
            end
            title([num2str(i),', ',num2str(M2{i}.DCrmse)]); 
            pause; 
            hold off; 
        end; 
    end; 
end