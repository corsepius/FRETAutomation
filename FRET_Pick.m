function Pick=FRET_Pick(X, m, n, CUT)

Pick = [];
%     for I=1:size(X,2)/2
%         I
%         XX(:,2*I-1)=smoothdata(X(:,2*I-1));
%         XX(:,2*I)=smoothdata(X(:,2*I));
%         
%         for i=1:size(XX,1)-2*n
%             y1(i,1) = mean(XX(i+n:i+2*n,2*I-1)) - mean(XX(i:i+n,2*I-1));
%             y2(i,1) = mean(XX(i+n:i+2*n,2*I)) - mean(XX(i:i+n,2*I));
%         end
%         
%         yy = y1.*y2;
%         
%         if (min(yy) < -1*CUT)
%            Pick = [Pick; I];
%         end
%         
%     end



    
Meas = FRET_Measure(X,m,n);

for i=1:size(Meas,2)
    if (min(Meas(:,i)) < -1*CUT)
        Pick=[Pick;i];
    end
end

fprintf('Picked: %d \r',size(Pick,1));
    
        