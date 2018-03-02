paychecks =
    c(9784.44,  # 2012
    12379.79,   # 2013
    20363.00,   # 2014
    30772.94,   # 2015
    40691.25,   # 2016
    31943.29,   # 2017
    2294.78,    # 2 jan 2018
    250,        # 26 jan
    2294.78,    # 1 feb
    2294.78)    # 1 mar

tithing =
    c(895,  # 2012
    1277,   # 2013
    2015,   # 2014
    3150,   # 2015
    4070,   # 2016
    3190,   # 2017
    230,
    255,
    225)

# Don't think this is accurate
offerings =
    c(158.25,   # 2012
    53,         # 2013
    355,        # 2014
    560,        # 2015
    390,        # 2016
    20, 40, 20, 20, 10, 10, 10, 20, 20, 10, 10, 10, # 2017
    10, 10, 15)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
