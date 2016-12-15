def row(string,nums,col=40,fcol=None, fmt = None, strip = '', sty = None, makeint = False, flt = None):
    if fcol is None: fcol = col
    string = str(string)
    #if len(string) > .8*fcol and ',' in string and '\n' not in string:
    #    strings = string.split(',')
    #    string = '\n'.join(strings)
    #if fmt=='latex': rowstr = string.replace('_','\_')
    #if fmt=='latex': rowstr = rowstr.replace('+/-','\pm')
    if fmt == 'twiki': 
        if sty == 'b': rowstr = '| * ' + string  + ' * ' 
        else: rowstr = '| ' + string 
    elif fmt == 'csv':
        rowstr = string +','
    else: rowstr = string
    if '\n' in string: string = string.split('\n')[-1]
    if not fmt: rowstr += ' '*(fcol-len(string))
    for n in nums:
        if type(n) is float: 
            if not flt:
                n = '%.4f' % n
            else:
                n = flt % n
            if makeint: n = int(n)
        else: n = str(n).lstrip(strip)
        if sty == 'b': n = ' * ' + n + ' * ' 
        if fmt == 'latex': rowstr+= ' & '
        elif fmt == 'twiki': rowstr+=' | '
        elif fmt == 'csv': rowstr+=','
        rowstr+=n
        if not fmt:
            rowstr+= ' ' * (col - len(n)) #+ '|'
    if fmt == 'latex': rowstr += ' \\\\'
    if fmt=='latex': rowstr = rowstr.replace('_','\_')
    if fmt == 'latex': rowstr = rowstr.replace(' +/- ',' $\pm$ ')
    elif fmt == 'twiki': rowstr += ' |'

    return rowstr

def section(string,l=40):
    return '#'*l +'\n%s\n' %string + '#'*l
