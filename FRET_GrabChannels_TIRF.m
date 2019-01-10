function A=FRET_GrabChannels_TIRF(X)

lenT=(size(X,2)-1)/3;

for i=1:lenT
    A(:,2*i-1)=X(:,3*i-1);
    A(:,2*i)=X(:,3*i);
end


    
        