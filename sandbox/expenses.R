paychecks =
    c(9784.44, # 2012
    12379.79,  # 2013
    20363.00,  # 2014
    687.50, # 23 jan 2015
    638.40, # 6 feb 2015
    638.40, # 20 feb 2015
    638.40, # 6 mar 2015
    638.40) # 20 mar 2015

tithing=
    c(895,  # 2012
    1277,   # 2013
    2015,   # 2014
    140, 60, 70)

offerings=
    c(158.25,   # 2012
    53,         # 2013
    355,        # 2014
    60, 30)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
