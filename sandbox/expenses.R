paychecks =
    c(9784.44, # 2012
    12379.79,  # 2013
    20363.00,  # 2014
    687.50, # 23 jan 2015
    638.40, # 6 feb 2015
    638.40, # 20 feb 2015
    638.40, # 6 mar 2015
    638.40, # 20 mar 2015
    638.40, # 3 apr 2015
    638.40, # 17 apr 2015
    638.35, # 1 may 2015
    343.75, # 15 may 2015
    20987.04, # 19 jun through 18 sep 2015 (LLNL)
    2142.17,    # 30 oct 2015
    2142.17)    # 1 dec 2015
    

tithing=
    c(895,  # 2012
    1277,   # 2013
    2015,   # 2014
    140, 60, 70, 80, 100,
    100, 70, 2310, 220)

offerings=
    c(158.25,   # 2012
    53,         # 2013
    355,        # 2014
    60, 30, 20, 200, 30, 190, 30)

sum(paychecks)/10 - sum(tithing)
# Negative means how much I overpaid
# Positive means how much I need to pay
