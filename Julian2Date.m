function [day, month] = Julian2Date(year,julian)
% This calculates the month and day, given the julian day and year
%
% Written by Cristen Adams.  June 19, 2009
%
% INPUT:  year, julian day that you want to calculate the date for 
%           (INTEGER)
% OUTPUT:  day, month (INTEGERS)
    
    % Determine whether it's a leap year
    leapyr = [1980; 1984; 1988; 1992; 1996; 2000; 2004; 2008; 2012; ...
        2016; 2020; 2024];
    is_leapyr = false;
    if sum(year == leapyr)
        is_leapyr = true;
    end

    % Set default values
    month = 0;
    day = 0;
    iter = julian;
    feb = 28;
    if is_leapyr
        feb = 29;
    end
    
    % Jan
    if iter <= 31
        month = 1;
        day = iter;
        return
    else
        iter = iter - 31;
    end
    
    % Feb
    if iter <= feb
        month = 2;
        day = iter;
        return
    else
        iter = iter - feb;
    end

    % Mar
    if iter <= 31
        month = 3;
        day = iter;
        return
    else
        iter = iter - 31;
    end
    
    % April
    if iter <= 30
        month = 4;
        day = iter;
        return
    else
        iter = iter - 30;
    end

    % May
    if iter <= 31
        month = 5;
        day = iter;
        return
    else
        iter = iter - 31;
    end
    
    % June
    if iter <= 30
        month = 6;
        day = iter;
        return
    else
        iter = iter - 30;
    end

    % July
    if iter <= 31
        month = 7;
        day = iter;
        return
    else
        iter = iter - 31;
    end

    % August
    if iter <= 31
        month = 8;
        day = iter;
        return
    else
        iter = iter - 31;
    end

    % September
    if iter <= 30
        month = 9;
        day = iter;
        return
    else
        iter = iter - 30;
    end

    % October
    if iter <= 31
        month = 10;
        day = iter;
        return
    else
        iter = iter - 31;
    end

    % November
    if iter <= 30
        month = 11;
        day = iter;
        return
    else
        iter = iter - 30;
    end
    
    % December
    month = 12;
    day = iter;
    
end