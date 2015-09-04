function parseGEANTetracks_func(in_name, out_name)

fid=fopen(in_name);
L=1;
linno=0;
nsteps=1e4;
ns=0;
out=zeros(nsteps,12);
while L~=-1
    L=fgets(fid);
    linno=linno+1;
    if mod(linno, 10000) == 0
        print linno
    end
    if strfind(L,'G4Track')
        tmp=strfind(L,'=');
        tmp2=strfind(L,',');
        particle=L(tmp(1)+1:tmp2(1)-1);
        trackid=L(tmp(2)+1:tmp2(2)-1);
        parentid=L(tmp(3)+1:end-1);
        fgets(fid);
        fgets(fid);
        fgets(fid);
        S=fgets(fid);
        linno=linno+4;
        
        while length(S)>1
            ns=ns+1;
            if ~mod(ns,nsteps) 
                ns
                out=[out; zeros(nsteps,12)];
            end
            tmp=strfind(S,' ');
            tc=zeros(1,8);
            indx=1;
            entry=[];
            for i=1:84 %length(S)
                if isempty(find(i==tmp,1))
                    entry=[entry S(i)];
                elseif ~isempty(entry) %if there's something in entry
                    if strfind('0 1 2 3 4 5 6 7 8 9 - .',entry(end)) %the entry was numeric
                        for i = 1:length(entry)
                            if strfind('0 1 2 3 4 5 6 7 8 9 - .',entry(i))
                                entry = entry(i:end);
                                break
                            end
                        end
                       
                        tc(indx)=str2double(entry);
                        indx=indx+1;
                    elseif strcmp(entry,'fm')
                        tc(indx-1)=tc(indx-1)/1e9;
                    elseif strcmp(entry,'pm')
                        tc(indx-1)=tc(indx-1)/1e6;
                    elseif strcmp(entry,'Ang')
                        tc(indx-1)=tc(indx-1)/1e4;
                    elseif strcmp(entry,'nm')
                        tc(indx-1)=tc(indx-1)/1e3;
                    elseif strcmp(entry,'mm')
                        tc(indx-1)=tc(indx-1)*1e3;
                    elseif strcmp(entry,'eV')
                        tc(indx-1)=tc(indx-1)/1e3;
                    elseif strcmp(entry,'MeV')
                        tc(indx-1)=tc(indx-1)*1e3;
                    end
                    entry=[];
                end
            end
            out(ns,1:8)=tc;
            if strcmp(particle,' e-')
                out(ns,9)=3;
            elseif strcmp(particle,' gamma')
                out(ns,9)=2;
            else
                out(ns,9)=-1;
                linno
                particle
            end            
            out(ns,10)=str2num(trackid);
            out(ns,11)=str2num(parentid);
            out(ns,12)=linno;
            S=fgets(fid);
            linno=linno+1;
            if linno==61
                linno;
            end
        end
    end
%     if ns > 20000
%         break
%     end
end
ns
out=out(1:ns,:);
fclose(fid);
save(out_name, 'out')

end
