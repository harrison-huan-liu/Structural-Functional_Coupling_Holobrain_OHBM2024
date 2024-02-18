function [Age, Gender] = Align_Confunder(label_information, PID)
    Age = zeros(size(PID, 1), 1);
    Gender = zeros(size(PID, 1), 1);
    for PID_i = 1:size(PID, 1)
        for i = 1:size(label_information,1)
            if str2num(PID{PID_i})==label_information{i,1}
                age_range = strsplit(label_information{i,5},{'-','+'});
                Age(PID_i,1) = str2num(age_range{1,1});
                if strcmp(label_information{i,4},'M')
                    Gender(PID_i,1)=1;
                else
                    Gender(PID_i,1)=0;
                end
            end
        end
    end
end