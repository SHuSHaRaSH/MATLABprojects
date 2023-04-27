%% function for Subtask 10
function y = ftfunc1(w)
    y = @(w) fillmissing(2*(w.*cos(3/2).*sin(w/2)-3.*sin(3/2).*cos(w/2))./(w.*w-9),'constant',2*sin(3));
end