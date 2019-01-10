function StateVec = FRET_StateVector_2(trans,bleach,ProcessStateAve,X0)

for i=1:size(trans,1)
    ttrans=sort([trans{i,1};bleach{i,1}]);
    m=1;
    state = ProcessStateAve{i,1}(m,1);
    for j=1:size(X0,1)
        if (m<=size(ttrans,1))
            if(j==ttrans(m,1))
                m=m+1;
                state = ProcessStateAve{i,1}(m,1);
            end
        end
        StateVec(j,i)=state;
    end
end