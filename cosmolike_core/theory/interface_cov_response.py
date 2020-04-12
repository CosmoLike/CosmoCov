from covNG_resp import *
import sys

def call_cov_response(k1,k2,z):
#	sys.stderr.write('args =## (%e, %e, %e)\n'%(k1,k2,z))
	return covNG_resp(k1, k2, z)

