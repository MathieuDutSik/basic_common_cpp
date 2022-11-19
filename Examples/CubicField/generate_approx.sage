

val = 2 *cos(2*pi/7)

def get_error(the_rat):
    thediff = abs(val - the_rat)
    thepow = 1
    while(True):
        thepow_new = thepow * 2
        if thediff > 1/thepow_new:
            return 1/thepow
        thepow = thepow_new


def get_approx_order(n):
    elist = continued_fraction_list(val, nterms=n)
    cf = continued_fraction(elist)
    the_rat = cf.value()
    return the_rat


for i in range(1,10):
    the_rat = get_approx_order(5*i)
    the_err = get_error(the_rat)
    print("{", the_rat, ",", the_err, "}")



    
