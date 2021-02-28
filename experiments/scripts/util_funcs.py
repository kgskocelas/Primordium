
'''Convert a string containing integers to a Python list of ints'''
# Format: 1,3,5->7,10 = [1,3,5,6,7,10]
# Seperate values with commas. Values can be single ints or inclusive ranges
def str_to_int_list(s):
    L = []
    vals = s.split(',')
    for val in vals:
        if '->' in val: # Range in form A->B
            val_split = val.split('->')
            if len(val_split) != 2:
                print('Error! Invalid integer range:', val)
                exit()
            lower = int(val_split[0])
            upper = int(val_split[1])
            for i in range(lower, upper + 1):
                L.append(i)
        else:
            L.append(int(val))
    return L

'''Convert a string containing integers to a Python list of floats'''
# Format: 1.1,2.3,5.9 = [1.1, 2.3, 5.9]
# Seperate values with commas.
def str_to_float_list(s):
    L = []
    vals = s.split(',')
    for val in vals:
        L.append(float(val))
    return L

'''Make sure there is a trailing slash (/) at the end of the string. Do not add if already present'''
# Useful for standardizing paths when part is user-submitted
def ensure_trailing_slash(s):
    if s[-1] != '/':
        return s + '/'
    return s


if __name__ == '__main__':
    # Test ints
    L1 = str_to_int_list('1,2,3')
    if(1 not in L1 or 2 not in L1 or 3 not in L1 or len(L1) != 3):
        print('Error! Expected:', [1,2,3], 'Received:', L1)
        exit()
    L2 = str_to_int_list('10,12->14,16')
    if(10 not in L2 or 12 not in L2 or 13 not in L2 or 14 not in L2 or 16 not in L2 or len(L2) != 5):
        print('Error! Expected:', [10,12,13,14,16], 'Received:', L2)
        exit()
    L3 = str_to_int_list('-2->2')
    if(-2 not in L3 or -1 not in L3 or 0 not in L3 or 1 not in L3 or 2 not in L3 or len(L3) != 5):
        print('Error! Expected:', [-2,-1,0,1,2], 'Received:', L3)
        exit()
    # Test floats
    L4 = str_to_float_list('2.4, 5.5, 6.1')
    if(2.4 not in L4 or 5.5 not in L4 or 6.1 not in L4 or len(L4) != 3):
        print('Error! Expected:', [2.4,5.5,6.1], 'Received:', L4)
        exit()
    L5 = str_to_float_list('-2.4, 5.5, -6.1')
    if(-2.4 not in L5 or 5.5 not in L5 or -6.1 not in L5 or len(L5) != 3):
        print('Error! Expected:', [-2.4,5.5,-6.1], 'Received:', L5)
        exit()
    # Trailing slash tests
    s1 = ensure_trailing_slash('./')
    if(s1 != './'):
        print('Error! Expected: ./  Received:', s1)
        exit()
    s2 = ensure_trailing_slash('foo')
    if(s2 != 'foo/'):
        print('Error! Expected: foo/  Received:', s2)
        exit()
    
    print('All tests pass!')
