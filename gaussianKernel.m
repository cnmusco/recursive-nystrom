function Ksub = gaussianKernel(D,rowInd,colInd,gamma)
    %% Guassian kernel generator
    % Outputs a submatrix of the Gaussian kernel with variance paramater 
    % gamma for the data rows of D. 
    %
    % usage : 
    %
    % input:
    %
    %  * D : A matrix with n rows (data points) and d columns (features)
    %
    %  * rowInd, colInd : Lists of indices between 1 and n. 
    %
    %  NOTE: colInd can be an empty list, in which case the **diagonal** 
    %  entries of the kernel will be output for the indices in rowInd.
    %  
    %  * gamma : kernel variance parameter
    %
    % output:
    %
    %  * Ksub : Let K(i,j) = e^-(gamma*||D(i,:)-D(j,:)||^2). Then Ksub = 
    %  K(rowInd,colInd). Or if colInd = [] then Ksub = diag(K)(rowInd).
    
    if(isempty(colInd))
        Ksub = ones(length(rowInd),1);
    else
        nsqRows = sum(D(rowInd,:).^2,2);
        nsqCols = sum(D(colInd,:).^2,2);
        Ksub = bsxfun(@minus,nsqRows,D(rowInd,:)*(2*D(colInd,:))');
        Ksub = bsxfun(@plus,nsqCols',Ksub);
        Ksub = exp(-gamma*Ksub);         
    end
end 