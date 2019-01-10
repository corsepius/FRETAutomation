function ProcessStateAve=FRET_ProcessStateAve_3(StateAve)

% This script will assume that the first state after the dye on is the off
% state

for n=1:size(StateAve,1)
    ProcessStateAve{n,1}=[1];
    for i=2:size(StateAve{n,1},1)
        if (StateAve{n,1}(i-1,1)==0 || StateAve{n,1}(i,1)==0)
            ProcessStateAve{n,1}=[ProcessStateAve{n,1};1];
        %This is not a bleach state nor is it the state following a bleach
        %state, which we assume is off
        else
            if (ProcessStateAve{n,1}(i-1,1)==1)
                % Can only go up
                ProcessStateAve{n,1}=[ProcessStateAve{n,1};2];
            else
                % Can only go down
                ProcessStateAve{n,1}=[ProcessStateAve{n,1};1];
            end
        end
    end
end

           