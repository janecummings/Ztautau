import math


def frexp10( val ):
    float(val)
    exp = int(math.log10(val))
    #print exp
    return val/10**exp, exp

def threedig( val ):
    if not 0.1 <= val < 1:
        val = frexp10( val ) [0]
    #print val
    return int ( round ( float (val) * 1000, 0 ) )

def _frexp10_ ( val ):
    m,e = frexp10( val )
    return m*10,e-1

def round_N ( val , n):
    return round( val, n)

def round_ ( val, n ):
    #print val
    #print frexp10(val)[0] - round(frexp10(val)[0])

    if frexp10(val)[0] - round(frexp10(val)[0]) >= 0.5:

        return round(val+1)
    else:
        return round(val)  

def pdg_round( val, err ):

    if val == 0:
        return 0,err

    ne = threedig( err )
    #print ne

    if abs ( ne ) <= 354: 
        pr,pe = 2,2
    elif abs (ne ) <= 949:
        pr, pe = 1,1
    else:
        pr, pe = 2,1

    ev, bv = frexp10(val)
    ee, be = frexp10(err)
    ##
    e = round_N ( err, pe )
    e = err
    ##
    #print be
    q  = be  - pe 
    r  = 10**q
    #print q
    #print val/r
    #print e/r
    #print val/r
    v  = round_ ( val / r , q ) * r
    e  = round_ ( e / r, q ) * r
    #print v
    #print e

    return v,e


def testpdg(val,err):
    print '%g +/- %g' % pdg_round(val,err)


if __name__ == '__main__':

    testpdg( 0.01151, 0.00261)


