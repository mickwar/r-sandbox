paychecks =
    c(9784.44, # 2012
    12379.79,  # 2013
    20363.00,  # 2014
    687.50) # 32 jan 2015

tithing=
    c(895,  # 2012
    1277,   # 2013
    160, 200, 200, 20, 150,
    260, 340, 145, 200, 70,
    65, 70, 65, 70)
    # do 140 (26 jan)


offerings=
    c(158.25,   # 2012
    43,         # 2013
    10, 6, 60, 10, 40, 30,
    35, 29, 35, 70)
    # d0 60 (26 jan)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
