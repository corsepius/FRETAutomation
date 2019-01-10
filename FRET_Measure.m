function Measure=FRET_Measure(X, m, n)

% m is the size of the intial smoothing window
% n is the size of the window used for measuring the difference

Measure = [];
    fprintf('Trace: ');
    for I=1:size(X,2)/2
        
        %Print frame number on screen
        if (I>10000)
            fprintf('\b\b\b\b\b');
        elseif (I>1000)
            fprintf('\b\b\b\b');
        elseif (I>100)
            fprintf('\b\b\b');
        elseif (I>10)
            fprintf('\b\b');
        elseif (I>1)
            fprintf('\b');
        end
        fprintf('%d',I);
        
        % Perform initial smoothing of window size m
        X(:,2*I-1)=smoothdata(X(:,2*I-1),m);
        X(:,2*I)=smoothdata(X(:,2*I),m);
        
        
        %Calculate FRET measure of window size n
        for i=1:size(X,1)-2*n
            y1(i,1) = mean(X(i+n:i+2*n,2*I-1)) - mean(X(i:i+n,2*I-1));
            y2(i,1) = mean(X(i+n:i+2*n,2*I)) - mean(X(i:i+n,2*I));
        end
        yy = y1.*y2;

        Measure = [Measure,yy];
    end
    fprintf('\r');
    
        