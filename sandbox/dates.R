get.dates = function(){
    day.of.week = c("Monday", "Tuesday", "Wednesday",
        "Thursday", "Friday", "Saturday", "Sunday")
    month = c("January", "February", "March", "April",
        "May", "June", "July", "August", "September",
        "October", "November", "December")
    ndays = c(31, 29, 31, 30, 31, 30,
        31, 31, 30, 31, 30, 31)

    years = seq(1900, 2020, by=1)
    # 1 if the year is common, 0 if it is a leap year
    non.leaps = 1-(((1-(years %% 4)>0) + (1-(years %% 100)>0) +
        (1-(years %% 400)>0)) %% 2)

    rep.days = rep(ndays, length(years))
    # gets all februarys for years of interest and substracts out
    # redundant day
    rep.days[seq(2, length(rep.days), by=12)] = rep.days[seq(2,
        length(rep.days), by=12)] - non.leaps

    rep.months = rep(rep(month, length(years)), rep.days)
    rep.days.in.month = unlist(apply(as.matrix(rep.days), 1, seq))
    rep.years = rep(years, 366-non.leaps)
    rep.day.of.week = rep(day.of.week, ceiling(length(rep.years)/7))
    rep.day.of.week = rep.day.of.week[1:length(rep.years)]


    dates = data.frame(matrix(0, length(rep.years), 4))
    dates[,1] = rep.day.of.week
    dates[,2] = rep.days.in.month
    dates[,3] = rep.months
    dates[,4] = rep.years
    names(dates) = c("Weekday", "Day", "Month", "Year")
    return (dates)
    }

dates = get.dates()

# dates[dates$Day == 25 & dates$Month == "March" & dates$Year == 2014, ]
# dates[dates$Weekday == "Sunday" & dates$Year == 2013, ]
