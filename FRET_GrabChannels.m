function A=FRET_GrabChannels(X,I,J)

lenT=(size(X,2)-1)/4;

for i=1:lenT
    A(:,2*i-1)=X(:,4*i+(I-3));
    A(:,2*i)=X(:,4*i+(J-3));
end


    
        