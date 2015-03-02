paychecks =
    c(9784.44, # 2012
    12379.79,  # 2013
    20363.00,  # 2014
    687.50, # 23 jan 2015
    638.40, # 6 feb 2015
    638.40) # 20 feb 2015

tithing=
    c(895,  # 2012
    1277,   # 2013
    160, 200, 200, 20, 150,
    260, 340, 145, 200, 70,
    65, 70, 65, 70,
    140, 60)
    # do 70 (8 mar 2015)

offerings=
    c(158.25,   # 2012
    43,         # 2013
    10, 6, 60, 10, 40, 30,
    35, 29, 35, 70,
    60)
    # do 30 (8 mar 2015)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
