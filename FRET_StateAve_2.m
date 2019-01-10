function StateAve=FRET_StateAve_2(X,Trans,Bleach)

clear StateAve


    for n=1:size(Trans,1)
         Trans{n,1}=sort([Trans{n,1};1;size(X,1);Bleach{n,1}]);
    end

    for n=1:size(Trans,1)
        StateAve{n,1} = [];
        for i=1:size(Trans{n,1},1)-1
            StateAve{n,1}=[StateAve{n,1};[mean(X(Trans{n,1}(i,1):Trans{n,1}(i+1,1),2*n-1)),mean(X(Trans{n,1}(i,1):Trans{n,1}(i+1,1),2*n))]];
        end
    end
           

    