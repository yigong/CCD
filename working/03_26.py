line = 'asdf asdf; adfas. fasdf,asdf,    foo'
import re
a = re.split(r'[;,\s]\s*', line)
print a
print 10*'-'

b = re.split(r'(;|,|\s)\s*', line)
print b 


