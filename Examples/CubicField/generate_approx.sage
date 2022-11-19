

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


f = open("CubicFieldDisc_49", "w")
f.write("3\n")
f.write("-1 2 1 1\n")
f.write(str(float(val)) + "\n");

f.write("10\n")
for i in range(1,10):
    the_rat = get_approx_order(5*i)
    the_err = get_error(the_rat)
    f.write(str(the_rat) + " " + str(the_err) + "\n")

f.close()
