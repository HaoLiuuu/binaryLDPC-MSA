function [hRowIndex,hColIndex,index_prefix,find_out_llr_index,zero_index,mink_tail,m,n] = Get_All_Index(H)
    [hRowIndex,hColIndex] = FindPosition(H);
    hColIndex = hColIndex';
    n = size(hRowIndex,2);
    m = size(hRowIndex,1);
    index_prefix = n * ((1 : m) - 1)';


    %% 得到由alpha_gamma转换为out_llr的index矩阵
    num_cnt      = zeros(m,2);
    num_cnt(:,1) = (1 : m)';
    num_cnt(:,2) = 1;
    find_out_llr_index = zeros(size(hColIndex,1),size(hColIndex,2));
    for j = 1 : size(hColIndex,2)
        for i = 1 : size(hColIndex,1)
            for k = 1 : m
                if hColIndex(i,j) == k
                    find_out_llr_index(i,j) = hColIndex(i,j) + (m * num_cnt(k,2));
                    num_cnt(k,2) = num_cnt(k,2) + 1;
                end
            end
        end
    end

    for i = 1 : size(find_out_llr_index,1)
        for j = 1 : size(find_out_llr_index,2)
            if find_out_llr_index(i,j) == 0
                find_out_llr_index(i,j) = 1;
            end
        end
    end

    zero_index = find(hRowIndex == 0);
    hRowIndex = hRowIndex + 1;
    hColIndex = hColIndex + 1;

    n = size(hRowIndex,2);
    m = size(hRowIndex,1);
    mink_tail = zeros(m,n - 2);
end