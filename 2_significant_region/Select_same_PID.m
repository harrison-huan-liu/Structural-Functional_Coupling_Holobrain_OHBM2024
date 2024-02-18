function [transform_ttest_i, transform_ttest_j, PID_com] = Select_same_PID(PID_i, PID_j, transform_ttest_i, transform_ttest_j)
    temp_i = 1;
    temp_j = 1;
    count_i = 1;
    count_j = 1;
    PID_com = cell(2,1);

    while count_i < size(PID_i,1)
        if any(contains(PID_j, PID_i{count_i}))
            temp_i = temp_i + 1;
            PID_com{temp_i,1} = PID_i{count_i};
        else
            transform_ttest_i(:,:,temp_i) = [];
        end
        count_i = count_i + 1;
    end
    while count_j < size(PID_j,1)
        if any(contains(PID_com, PID_j{count_j}))
            temp_j = temp_j + 1;
        else
            transform_ttest_j(:,:,temp_j) = [];
        end
        count_j = count_j + 1;
    end
end