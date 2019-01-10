function [Measure,Trans, Bleach]=FRET_Measure_Square_Trans_All(X, m, n, cutoff, pos_cut)


%Measure peak width window for determining if at minimum (indicating transition)
window = 7;


Measure = [];
    for I=1:size(X,2)/2      
        X(:,2*I-1)=smooth(X(:,2*I-1),m);
        X(:,2*I)=smooth(X(:,2*I),m);
        
        
        for i=1:size(X,1)-2*n
            y1(i,1) = mean(X(i+n:i+2*n,2*I-1)) - mean(X(i:i+n,2*I-1));
            y2(i,1) = mean(X(i+n:i+2*n,2*I)) - mean(X(i:i+n,2*I));
        end
        
        yy = y1.*y2;

        Measure = [Measure,yy];
    end
    
    Y=Measure;
    Trans = cell(size(Y,2),1);
    Bleach = cell(size(Y,2),1);
   for j=1:size(Y,2)
        for i=1:size(Y,1)
            %On a peak, finding minimum
            if (Y(i,j)<cutoff)
                if (i-window<1) 
                    if (Y(i,j)==min(Y(1:i+window,j)))
                        Trans{j,1}(end+1,1)=i+n;
                    end
                elseif (i+window>size(Y,1))
                    if (Y(i,j)==min(Y(i-window:size(Y,1),j)))
                        Trans{j,1}(end+1,1)=i+n;
                    end
                else
                    if (Y(i,j)==min(Y(i-window:i+window,j)))
                        Trans{j,1}(end+1,1)=i+n;
                    end
                end
            %Finding initial on and bleaching events 
            elseif (Y(i,j)>pos_cut)
                if (i-window<1) 
                    if (Y(i,j)==max(Y(1:i+window,j)))
                        Bleach{j,1}(end+1,1)=i+n;
                    end
                elseif (i+window>size(Y,1))
                    if (Y(i,j)==max(Y(i-window:size(Y,1),j)))
                        Bleach{j,1}(end+1,1)=i+n;
                    end
                else
                    if (Y(i,j)==max(Y(i-window:i+window,j)))
                        Bleach{j,1}(end+1,1)=i+n;
                    end
                end
            end
        end
   end



        