paychecks =
    c(9784.44,  # 2012
    12379.79,   # 2013
    20363.00,   # 2014
    30772.94,   # 2015
    40691.25,   # 2016
    2227.89,    # 3 jan 2017
    2227.89,    # 1 feb 2017
    2227.89)    # 1 mar 2017

tithing=
    c(895,  # 2012
    1277,   # 2013
    2015,   # 2014
    3150,   # 2015
    4070,   # 2016
    230, 210, 230)

offerings=
    c(158.25,   # 2012
    53,         # 2013
    355,        # 2014
    560,        # 2015
    390,        # 2016
    20, 40, 20)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
