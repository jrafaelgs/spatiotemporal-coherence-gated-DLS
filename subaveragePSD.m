function psdOutSub = subaveragePSD(psdOut,n)

[row col] = size(psdOut);
col_aux = floor(col/n);

if n <= 1
    psdOutSub = psdOut;
    display(['PSD Matrix not reduced, still ' num2str(col) ' colums.'])
else
    psdOutSub = zeros(row,col_aux);
    aux_ind = (n-1:-1:0);
    for i = 1:col_aux
        psdOutSub(:,i) = mean((psdOut(:,(n*i-aux_ind)))');
    %     display((n*i-aux_ind))
    end
    display(['PSD Matrix reduced from ' num2str(col) ' to ' num2str(col_aux)])
end