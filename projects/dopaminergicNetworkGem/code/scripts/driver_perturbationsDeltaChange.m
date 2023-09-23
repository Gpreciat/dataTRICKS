measured_dv = abs(measured_control_data.mean - measured_test_data.mean);
measured_dvsd = sqrt(measured_control_data.SD.^2 + measured_test_data.SD.^2);

for i=1:length(measured_dv)
    if measured_control_data.mean(i) > 0 && measured_test_data.mean(i) > 0
        if measured_control_data.mean(i) > measured_test_data.mean(i)
            measured_label{i,1} = 'SS_down';
        elseif measured_control_data.mean(i) < measured_test_data.mean(i)
            measured_label{i,1} = 'SS_up';
        else
            measured_label{i,1} = 'SS';
        end
    elseif measured_control_data.mean(i) < 0 && measured_test_data.mean(i) < 0
        if measured_control_data.mean(i) > measured_test_data.mean(i)
            measured_label{i,1} = 'UU_up';
        elseif measured_control_data.mean(i) < measured_test_data.mean(i)
            measured_label{i,1} = 'UU_down';
        else
            measured_label{i,1} = 'UU';
        end
    elseif measured_control_data.mean(i) < 0 && measured_test_data.mean(i) >= 0
        if measured_test_data.mean(i) == 0
            measured_label{i,1} = 'U0';
        else
            measured_label{i,1} = 'US';
        end  
    elseif measured_control_data.mean(i) > 0 && measured_test_data.mean(i) <= 0
        if measured_test_data.mean(i) == 0
            measured_label{i,1} = 'S0';
        else
            measured_label{i,1} = 'SU';
        end
    elseif measured_control_data.mean(i) == 0
        if measured_test_data.mean(i) == 0
            measured_label{i,1} = '00';
        elseif measured_test_data.mean(i) < 0
            measured_label{i,1} = '0U';
        elseif measured_test_data.mean(i) > 0
            measured_label{i,1} = '0S';
        end  
    end
end

predicted_dv = abs(predicted_control_data.v - predicted_test_data.v);

for i=1:length(predicted_dv)
    if predicted_control_data.v(i) > 0 && predicted_test_data.v(i) > 0
        if predicted_control_data.v(i) > predicted_test_data.v(i)
            predicted_label{i,1} = 'SS_down';
        elseif predicted_control_data.v(i) < predicted_test_data.v(i)
            predicted_label{i,1} = 'SS_up';
        else
            predicted_label{i,1} = 'SS';
        end
    elseif predicted_control_data.v(i) < 0 && predicted_test_data.v(i) < 0
        if predicted_control_data.v(i) > predicted_test_data.v(i)
            predicted_label{i,1} = 'UU_up';
        elseif predicted_control_data.v(i) < predicted_test_data.v(i)
            predicted_label{i,1} = 'UU_down';
        else
            predicted_label{i,1} = 'UU';
        end
    elseif predicted_control_data.v(i) < 0 && predicted_test_data.v(i) >= 0
        if predicted_control_data.v(i) == 0
            predicted_label{i,1} = 'U0';
        else
            predicted_label{i,1} = 'US';
        end  
    elseif predicted_control_data.v(i) > 0 && predicted_test_data.v(i) <= 0
        if predicted_control_data.v(i) == 0
            predicted_label{i,1} = 'S0';
        else
            predicted_label{i,1} = 'SU';
        end
    elseif predicted_control_data.v(i) == 0
        if predicted_control_data.v(i) == 0
            predicted_label{i,1} = '00';
        elseif predicted_test_data.v(i) < 0
            predicted_label{i,1} = '0U';
        elseif predicted_test_data.v(i) > 0
            predicted_label{i,1} = '0S';
        end  
    end
end