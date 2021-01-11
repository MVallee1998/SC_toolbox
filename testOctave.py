from oct2py import octave
import timeit
octave.addpath('./')
start = timeit.default_timer()
octave.test(100)
stop = timeit.default_timer()
print(stop-start)