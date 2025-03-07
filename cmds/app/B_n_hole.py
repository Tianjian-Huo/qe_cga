import sys
sys.path.append('../src')
from material import *
from math import *
from structure import *

c1 = ring(8)
c2 = ring(16)
c2.move(Point(0,0,1))
c3 = ring(24)
c3.move(Point(0,0,2))
c = c1+c2+c3
c.set_element('B')
c.zoom(1.4)
c.write_car('B48_8hole.car')