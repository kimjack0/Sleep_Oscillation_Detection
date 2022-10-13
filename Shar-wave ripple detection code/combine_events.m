function [startidx,finishidx,start,finish] = ...
    combine_events(upidx,dwnidx,up,dwn,th)
% recursively combine events that are in close proximity to each other
% (defined by th). th and time units should be consistent.
% 
%   Inputs:
%       upidx  - upward thresh crossings (index)
%       up     - upward thresh crossings (time)
%       dwnidx - downward thresh crossings (index)
%       dwn    - downward thresh crossings (time)
% 
%   Outputs:
%       startidx  - start of events (index)
%       finishidx - start of events (time)
%       start     - end of events (index)
%       finish    - end of events (time)

startidx = [];
finishidx = [];
start = [];
finish = [];

flag = 0;
i=1;
while i<=(length(up)-1),
    startidx(end+1,1)   = upidx(i);
    finishidx(end+1,1)  = dwnidx(i);
    start(end+1,1)      = up(i);
    finish(end+1,1)     = dwn(i);
    if (up(i+1)-finish(end))<th,
        finishidx(end)  = dwnidx(i+1);
        finish(end)     = dwn(i+1);
        flag = 1;
        i = i + 2;
    else
        i = i + 1;
    end
end

if flag,
    [startidx,finishidx,start,finish] = ...
        combine_events(startidx,finishidx,start,finish,th);
end

