function progressfile(varargin)
%function progressfile(varargin)
%
%For running scripts remotely.
%Creates a file on the host computer with the task progress in the filename.
%
%Syntax:
% progressfile
% progressfile(filename)
% progressfile(progress)
% progressfile(progress, filename)
%
% progressfile: initialize (sets clock)
% progressfile(filename): initialize with custom filename
%
% progressfile(progress): updates the file with 'progress' from 0 to 1 measuring progress of task
%
% progressfile(progress, filename): update with custom filename
%
%Filenames:
% use %p for the progress indicator ('34' for 34 percent)
% use %r for time remaining ('4h' for 4 hours)
% use %e for time elapsed ('22m' for 22 minutes)
%Default filename is 'progress%ppercent' e.g. 'progress34percent'
%Filename automatically given .mat extension

filename = [];
default = 'progress%ppercent';
endflag = false;

try
    curtime = toc;
catch
    %only start stopwatch anew if it has not been started by any other progressfile
    tic
    curtime = 0;
end

%Check input arguments

if nargin==0
    %initialize with default filename
    t0 = curtime; %initial time of task
    progress = 0;
    
elseif nargin==1
    if ischar(varargin{1})
        %initialize with custom filename
        t0 = curtime;
        filename = varargin{1};
        progress = 0;
    else
        %update with default filename
        progress = varargin{1};
        if progress >= 1
            endflag = true;
        elseif progress < 0
        elseif progress == 0
            t0 = curtime;
        elseif progress < 1 && progress > 0
            
        end
    end
    
elseif nargin==2
    %update with custom filename
    progress = varargin{1};
    if progress >= 1
        endflag = true;
    elseif progress < 0
    elseif progress == 0
        t0 = curtime;
    end
    filename = varargin{2};
end

%Use default filename if necessary
if isempty(filename)
    filename = default;
end
writename = filename;

%Find previous file
oldfilename = [];
if progress > 0
    files = dir;
    %to ID file, look at text between escape characters
    s = strfind(filename,'%');
    if isempty(s)
        return
    end
    temp = filename;
    for i = 1:length(s)
        %clear escape characters inserted in filename
        if strcmp(temp(1),'%')
            temp = temp(3:length(temp));
            s(i:length(s)) = s(i:length(s)) - 2;
        end
        %grab text as identifying string
        st{i} = temp(1:s(i)-1);
        temp = temp(s(i):length(temp));
        if i<length(s)
            s(i+1:length(s)) = s(i+1:length(s)) - (s(i)-1);
        end
    end
    
    if strcmp(temp(1),'%')
        temp = temp(3:length(temp));
    end
    st{i+1} = temp; %last part of filename
    clear temp
    
    %look through files
    for k = 1:length(files)
        next = false;
        if isdir(files(k).name); continue; end
        temp = files(k).name;
        for m = 1:length(st)
            if isempty(st{m})
                continue
            else
                pos = strfind(temp,st{m});
                if ~isempty(pos)
                    %cut temp up til end of this string, then look for next
                    temp = temp(pos + length(st{m}):length(temp));
                else
                    %string not found
                    %go to next file
                    next = true;
                    break
                end
            end
        end
        if next; continue; end %next file
        
        %If next is (still) false, this is the file.
        oldfilename = files(k).name;
        break
    end
    
    if isempty(oldfilename)
        %no file matches
        t0 = curtime;   %use current time as t0, even though progress>0
    else
        %grab t0 and delete file
        load(oldfilename,'t0')
        delete(oldfilename)
        if endflag
            %we're done
            return
        end
    end
end

%Create actual name of file

%Progress value
p = strfind(writename,'%p');
if ~isempty(p)
    writename = [writename(1:p-1),num2str(floor(progress*100)),writename(p+2:length(writename))];
end

%Time elapsed
te = strfind(writename,'%e');
elapsedsec = curtime - t0;
if ~isempty(te)
    %convert to good units
    if elapsedsec > 172800  %day
        elapsedday = floor(elapsedsec / 86400);
        elapsed = [num2str(elapsedday),'d',num2str(floor((elapsedsec - 86400*elapsedday)/3600)),'h'];
    elseif elapsedsec > 7200    %hour
        elapsed = [num2str(floor(elapsedsec / 3600)),'h'];
    elseif elapsedsec > 120 %minute
        elapsed = [num2str(floor(elapsedsec / 60)),'m'];
    else
        elapsed = [num2str(floor(elapsedsec)),'s'];
    end
    writename = [writename(1:te-1),elapsed,writename(te+2:length(writename))];
end

%Time remaining
tr = strfind(writename,'%r');
if progress == 0
    remaining = '_';
else
    remainingsec = elapsedsec / progress * (1-progress);
    %convert to good units
    if remainingsec > 172800  %day
        remainingday = floor(remainingsec / 86400);
        remaining = [num2str(remainingday),'d',num2str(floor((remainingsec - 86400*remainingday)/3600)),'h'];
    elseif remainingsec > 7200    %hour
        remaining = [num2str(floor(remainingsec / 3600)),'h'];
    elseif remainingsec > 120 %minute
        remaining = [num2str(floor(remainingsec / 60)),'m'];
    else
        remaining = [num2str(floor(remainingsec)),'s'];
    end
end

if ~isempty(tr)
    writename = [writename(1:tr-1),remaining,writename(tr+2:length(writename))];
end

% if ~strcmp(writename(length(writename)-3:length(writename)),'.mat')
%     %add extension

%.mat extension added by save

save(writename,'t0')


