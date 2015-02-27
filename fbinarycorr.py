#!/usr/bin/env python
#
#
# https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/util/debinorb.pro
# https://lost-contact.mit.edu/afs/physics.wisc.edu/home/craigm/lib/idl/bary/bnrybt.pro
 
import os
import yaml
import numpy
import pyfits
from optparse import OptionParser
from operator import div, mul
from datetime import datetime 

__name__	= 'fbinarycorr'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__	= '2015 Feb. 20'

# ================================================================================

def debinorb(time, asini, porb, t90, ecc, omega_d):
	# de-binary_orbit
	# time is array
	# see also, http://astro.uni-tuebingen.de/software/idl/aitlib/astro/binarypos.pro

	if len(time) == 0:
		print "no time array"
		quit()

	# Compute time shifts due to orbit around center of mass
	twopi   = 2.0 * numpy.pi
	t       = time 
	asini_d = asini
	print "[debinorb] ... start ... "
	print "[debinorb] asini", asini
	print "[debinorb] porb", porb
	print "[debinorb] t90", t90
	print "[debinorb] ecc", ecc
	print "[debinorb] omega_d", omega_d

	if ecc > 0.0:
		omega = numpy.radians(omega_d)
		sinw  = numpy.sin(omega)
		cosw  = numpy.cos(omega)
		for i in range(0,6):
			m = twopi * (t-t90)/porb + numpy.pi/2.0 - omega
			eanom = m
			for j in range(0,5):
				eanom = eanom - (eanom-ecc*numpy.sin(eanom)-m)/(1.0-ecc*numpy.cos(eanom))		
			m = 0 		
			sin_e = numpy.sin(eanom)
			cos_e = numpy.cos(eanom)
			f  = (t-time)+asini_d*(sinw*(cos_e-ecc)+numpy.sqrt(1.0-ecc*ecc)*cosw*sin_e)
			df = 1.0 + (twopi*asini_d/(porb*(1.0-ecc*cos_e)))*(numpy.sqrt(1.0-ecc*ecc)*cosw*cos_e-sinw*sin_e)
			sin_e = 0 
			cos_e = 0

			t = t - f/df
	else:
		# Note: the value for small values of asini/porb, the difference
		# between corrected and uncorrected times is approximated by:
		#
		# (ta - time)/asini = (- (1D/1D)*!dpi^0*(asini/porb)^0*cos(1D*phi) 
		#                      - (1D/1D)*!dpi^1*(asini/porb)^1*sin(2D*phi)
		#                      + (3D/2D)*!dpi^2*(asini/porb)^2*cos(3D*phi)
		#                      + (8D/3D)*!dpi^3*(asini/porb)^3*sin(4D*phi))
		#           [ phi = 2D*!dpi*(time-t90)/porb ]
		#    with the error [= abs(t - ta)] of LT asini^4/porb^3		
		for i in range(0,6):
			L = twopi*(t - t90)/porb + numpy.pi/2.0
			f = (t - time) + asini_d*numpy.sin(L)
			df = 1.0 + twopi*asini_d*numpy.cos(L)/porb
			t = t - f/df
	print "[debinorb] ... finished ... "			
	return t

"""
;==============================================================================
;   DEBINORB.PRO
;   09 Feb 2000
;   Markwardt (original due to Deepto Chakrabarty, Caltech)
;
;   Take an input time array and correct the values to place
;   them in a reference frame inertial wrt the pulsar. Uses
;   an exact elliptical orbit model, and inverts. 
;
;   ARGUMENTS:
;	time:	event time (TDB)
;	asini:	Projected semi-major axis [lt-secs]  
;	porb:	Orbital period
;	t90:	Epoch for mean longitude of 90 degrees
;	ecc:	Eccentricity
;	omega_d:Longitude of periastron [degrees]
;
;   ALL times must be in the same unit and coordinate system.
;
;   RETURNS an array of same type, length, and units as the input 
;	time array. 
;	
;
;   REVISIONS:
;   11-23-92 Deepto: Bug fixed for circular orbits.
;   09-01-97 Deepto: Bug fixed for ecc orbits: loop indexing error
;   2000 Feb 09 Markwardt: changed to debinorb from removeorb; changed
;                          to all same time units; !dpi instead of !pi
;==============================================================================

function debinorb, time, asini, porb, t90, ecc, omega_d

if (n_elements(time) eq 0)then begin
   print,$
  'ttcorr = removeorb(tt,asini,porb,t90,ecc,omega(deg))'
   return,-1
ENDIF 
; Compute time shifts due to orbit around center of mass
twopi = 2.0D*!dpi
t = time
asini_d = asini
if (ecc gt 0.0) then begin	
        omega = omega_d * !dpi/180.0D
        sinw = sin(omega)
        cosw = cos(omega)
	for i=0, 5 do begin
		m = twopi*(t - t90)/porb + !dpi/2.0D - omega
		eanom = m
		for j=0, 4 do begin
		  eanom = eanom $
			- (eanom - ecc*sin(eanom) - m)/(1.0D -ecc*cos(eanom))
		endfor
                m = 0
		sin_e = sin(eanom)
		cos_e = cos(temporary(eanom))
		f = (t - time) + $
                  asini_d*(sinw*(cos_e-ecc)+sqrt(1.0D -ecc*ecc)*cosw*sin_e)
		df = 1.0D + (twopi*asini_d/(porb*(1.0D -ecc*cos_e)))* $
                  (sqrt(1.0D -ecc*ecc)*cosw*cos_e - sinw*sin_e)
                sin_e = 0 & cos_e = 0
		t = temporary(t) - temporary(f)/temporary(df)
	endfor
endif else begin
    ;; Note: the value for small values of asini/porb, the difference
    ;; between corrected and uncorrected times is approximated by:

    ;; (ta - time)/asini = (- (1D/1D)*!dpi^0*(asini/porb)^0*cos(1D*phi) 
    ;;                      - (1D/1D)*!dpi^1*(asini/porb)^1*sin(2D*phi)
    ;;                      + (3D/2D)*!dpi^2*(asini/porb)^2*cos(3D*phi)
    ;;                      + (8D/3D)*!dpi^3*(asini/porb)^3*sin(4D*phi))
    ;;           [ phi = 2D*!dpi*(time-t90)/porb ]
    ;;    with the error [= abs(t - ta)] of LT asini^4/porb^3
	for i=0, 5 do begin
		L = twopi*(t - t90)/porb + !dpi/2.0D
		f = (t - time) + asini_d*sin(L)
		df = 1.0D + twopi*asini_d*cos(temporary(L))/porb
		t = temporary(t) - temporary(f)/temporary(df)
	endfor
endelse
return, t
end
"""

# ================================================================================

def bnrybt(ct, a1, pb, t0, e, omz,
	pbdot=0.0, xpbdot=0.0, edot=0.0, xdot=0.0,
	omdot=0.0, gamma=0.0):

	print "[bnrybt] ... start ... "
	print "[bnrybt] a1", a1
	print "[bnrybt] pb", pb
	print "[bnrybt] t0", t0
	print "[bnrybt] e", e
	print "[bnrybt] omz", omz

	twopi = 2.0 * numpy.pi
	rad   = 180.0 / numpy.pi 
	n     = len(ct)

	tt0   = (ct-t0)
	orbits = tt0 / pb - 0.5 * (pbdot+xpbdot)*(tt0/pb)**2.0
	ecc = e + edot*tt0
	asini = a1 + xdot * tt0

	norbits = numpy.floor(orbits)
	phase = twopi * (orbits-norbits)
	omega=(omz+omdot*tt0/(86400.0*365.25))/rad

	ep=phase + ecc*numpy.sin(phase)*(1.0+ecc*numpy.cos(phase))
	dep = ep*0.0
	denom=1.0-ecc*numpy.cos(ep)

	wh = (numpy.arange(n))
	while True:
		dep[wh]=map(mul,(phase[wh]-(ep[wh]-map(mul,ecc[wh],numpy.sin(ep[wh])))),1./denom[wh])		
		ep[wh]=ep[wh]+dep[wh]
		wh = numpy.where(numpy.fabs(dep)>1e-12)[0]
		if len(wh) == 0:
			break

	bige = ep
	tt   = 1.0 - ecc**2
	som  = numpy.sin(omega)
	com  = numpy.cos(omega)
	alpha= asini*som
	beta = asini*com*numpy.sqrt(tt)
	sbe  = numpy.sin(bige)
  	cbe  = numpy.cos(bige)
  	q    = alpha*(cbe-ecc) + (beta+gamma)*sbe
  	r    =-alpha*sbe + beta*cbe
  	s    = 1.0/(1.0-ecc*cbe)
	torb = -q + (twopi/pb)*q*r*s 

	return torb


"""
;
; BNRYBT - binary time delay model of Blandford and Teukolsky 1976
;
; This model computes the pulse time of arrival delay due to a binary
; pulsar system.  It follows the BNRYBT.F derivation of the TEMPO
; pulse timing software.
;
;  This version uses the method of Blandford and Teukolsky (ApJ 205,
;  580,1976)--including a solution for gamma.  For notation, see Taylor
;  et al, ApJ Lett 206, L53, 1976.
;
; CT - Input time
; A1 - a1 sin(i) - projected semimajor axis (lt-s)
; PB - binary period (s)
; T0 - time of reference periastron passage (and epochs of Pspin &
;      OMZ) (s)
; E - eccentricity
; OMZ - omega_z - longitude of periastron (deg)
;
; PBDOT - binary period derivative (s/s)
; XPBDOT - additional orbital period derivative (s/s)
; EDOT - eccentricity derivative (/s)
; XDOT - derivative of A1*SIN(I)
;
; RETURNS - time delay [s] due to binary orbit.  
;
; $Id: bnrybt.pro,v 1.8 2010/05/04 21:02:35 craigm Exp $
;
function bnrybt, ct, a1, pb, t0, e, omz, $
                 pbdot=pbdot0, xpbdot=xpbdot0, edot=edot0, $
                 xdot=xdot0, omdot=omdot0, gamma=gamma0

  if n_elements(pbdot0) GT 0 then pbdot = pbdot0(0) $ ;; binary orbital pdot
  else pbdot = 0   
  if n_elements(xpbdot0) GT 0 then xpbdot = xpbdot0(0) $ ;; addl binary pdot
  else xpbdot = 0  
  if n_elements(edot0) GT 0 then edot = edot0(0) $ ;; eccentricity derivative
  else edot = 0    
  if n_elements(xdot0) GT 0 then xdot = xdot0(0) $ ;; asini derivative
  else xdot = 0    
  if n_elements(omdot0) GT 0 then omdot = omdot0(0) $ ;; omega derivative
  else omdot = 0   
  if n_elements(gamma0) GT 0 then gamma = gamma0(0) $ ;; PPN gamma value
  else gamma = 0   

  twopi = 2D * !dpi
  rad = 180d/!dpi
  n = n_elements(ct)
  

  ;; NOTE: original BNRYBT code accounts for multiple planets, but
  ;; this code does not.
  ;; Apply omdot, pbdot, xdot and edot to first orbit.

  tt0=(ct-t0)
  orbits=tt0/pb - 0.5d0*(pbdot+xpbdot)*(tt0/pb)^2
  ecc=e+edot*tt0
  asini=a1+xdot*tt0

  norbits=floor(orbits)
  phase=twopi*(orbits-norbits)
  omega=(omz+omdot*tt0/(86400.d0*365.25d0))/rad

  ;;  Use Pat Wallace's method of solving Kepler's equation
  ;; ep = E' of BT equation 2.32
  ;; phase = r.h.s. of eqn 2.32
  ep=phase + ecc*sin(phase)*(1.d0+ecc*cos(phase))
  dep = ep*0
  denom=1.d0-ecc*cos(ep)

  wh = lindgen(n)
  repeat begin
      dep(wh)=(phase(wh) - (ep(wh)-ecc*sin(ep(wh))) )/denom(wh)
      ep(wh)=ep(wh)+dep(wh)
      wh = where(abs(dep) GT 1d-12, ctx)
  endrep until ctx EQ 0
  
  bige=ep
  tt=1.d0-ecc*ecc
  som=sin(omega)
  com=cos(omega)
  alpha=asini*som
  beta=asini*com*sqrt(tt)
  sbe=sin(bige)
  cbe=cos(bige)
  q=alpha*(cbe-ecc) + (beta+gamma)*sbe
  r=-alpha*sbe + beta*cbe
  s=1.d0/(1.d0-ecc*cbe)
  ;; BT equation 2.33
  ;; -q is first two terms, and remainder is last term of 2.33
  torb= -q + (twopi/pb)*q*r*s 

  return, torb
end
"""

usage  = "usage: %s in.fits out.fits [options] \n" % __name__
usage += "   - binary correction with an input yaml file ephemeris.\n"
usage += "\n"
usage += "propt> cat ephemeris.yaml\n"
usage += "orb_porb_day:  1.70016723662    # Orbital Period (day)\n"
usage += "orb_t0_sec:    181827904.948 # Time of periastron passage (sec), Mission Time \n"
usage += "orb_eccen:     4.2e-4           # Eccetricity (no-dimension)\n"
usage += "orb_asini_sec: 13.1831          # Asini (light second) a * sin(i) / c (s)\n"
usage += "orb_w_deg:     96.0             # Longitude Periastron (degree)\n"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--filldelay",
	help="fill a new column of DELAY TIME (bool)?", 
	action="store_true", dest="filldelay", default=False)	
parser.add_option("-c", "--clobber",
	help="overwrite the output file (bool)?", 
	action="store_true", dest="clobber", default=False)	
parser.add_option("-t", "--functype", type='str',
	help="used method (debinorb, bnrybt)", 
	action="store", dest="functype", default='debinorb')	
(options, args) = parser.parse_args()

if len(args) != 3:
	print "> %s.py inputfits outputfits ephem.yaml " % __name__
	quit()
inputfits  = args[0]	
outputfits = args[1]
yamlfile   = args[2]

if not os.path.exists(yamlfile):
	print "%s does not exists." % yamlfile
	quit()
if not os.path.exists(inputfits):
	print "input file does not exists: %s" % inputfits
	quit()
if options.clobber:
	cmd = 'rm -f %s' % outputfits
	print cmd; os.system(cmd)
elif os.path.exists(outputfits):
	print "output file has already existed: %s " % outputfits
	quit()

param    = yaml.load(open(yamlfile))
porb = float(param['orb_porb_day'])*24.0*60.0**2

headerfile = 'temp_header.txt'
f = open(headerfile,'w')
out = """HISTORY -----------------------------------------------------
HISTORY  %s version %s at %s
HISTORY -----------------------------------------------------
HISTORY   inputfits='%s'
HISTORY   outputfits='%s'
HISTORY   orb_porb_day  = %.10f / Orbital Period (day)
HISTORY   orb_t0_sec    = %.6f  / Time of periastron passage (sec), Mission Time 
HISTORY   orb_eccen     = %.5e  / Eccetricity (no-dimension)
HISTORY   orb_asini_sec = %.5f  / Asini (light second) a * sin(i) / c (s)
HISTORY   orb_w_deg     = %.5f  / Longitude Periastron (degree)
HISTORY   function_type = %s    / Function Type
""" % (__name__, __version__, datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
	inputfits, outputfits, 
	param['orb_porb_day'], param['orb_t0_sec'],
	param['orb_eccen'], param['orb_asini_sec'],
	param['orb_w_deg'], options.functype)
f.write(out)
f.close()
print out

hdu      = pyfits.open(inputfits)
time_org = hdu['EVENTS'].data['TIME']
start_org = hdu['GTI'].data['START']
stop_org = hdu['GTI'].data['STOP']
if options.functype == 'debinorb':
	time_new = debinorb(time=time_org,
		asini=float(param['orb_asini_sec']), 
		porb=porb, 
		t90=float(param['orb_t0_sec']), 
		ecc=float(param['orb_eccen']), 
		omega_d=float(param['orb_w_deg']))
	start_new = debinorb(time=start_org,
		asini=float(param['orb_asini_sec']), 
		porb=porb, 
		t90=float(param['orb_t0_sec']), 
		ecc=float(param['orb_eccen']), 
		omega_d=float(param['orb_w_deg']))
	stop_new = debinorb(time=stop_org,
		asini=float(param['orb_asini_sec']), 
		porb=porb, 
		t90=float(param['orb_t0_sec']), 
		ecc=float(param['orb_eccen']), 
		omega_d=float(param['orb_w_deg']))
elif options.functype == 'bnrybt':
	time_new = time_org + bnrybt(
		ct=time_org, 
		a1=float(param['orb_asini_sec']), 
		pb=porb, 
		t0=float(param['orb_t0_sec']), 
		e=float(param['orb_eccen']), 
		omz=float(param['orb_w_deg']))
	start_new = start_org + bnrybt(
		ct=start_org, 
		a1=float(param['orb_asini_sec']), 
		pb=porb, 
		t0=float(param['orb_t0_sec']), 
		e=float(param['orb_eccen']), 
		omz=float(param['orb_w_deg']))	
	stop_new = stop_org + bnrybt(
		ct=stop_org, 
		a1=float(param['orb_asini_sec']), 
		pb=porb, 
		t0=float(param['orb_t0_sec']), 
		e=float(param['orb_eccen']), 
		omz=float(param['orb_w_deg']))		
else:
	print "invalid binary correction function."
	quit()
print "... Corrected time is calculated."


print "... Start recording events."
if options.filldelay:
	delay_time = time_new - time_org
	delay_start = start_new - start_org
	delay_stop  = stop_new  - stop_org
Ndisp = 100000
print "...[events]..."
for i in range(len(time_org)):
	if i % Ndisp == 0: 
		print "...%d/%d (%.1f%%)" % (i,len(time_org),float(i)/float(len(time_org))*100.0)
	time_org[i] = time_new[i]	
print "...[gti]..."		
for i in range(len(start_new)):	
	start_org[i] = start_new[i]
	stop_org[i]  = stop_new[i]
hdu.writeto(outputfits)

if options.filldelay:
	print "... Prepare blank fits column for the delta-time"
	cmd  = 'fcalc infile=%s+1 ' % outputfits
	cmd += 'outfile=%s ' % outputfits
	cmd += 'clobber=yes '
	cmd += 'clname=\"DELTA_TIME\" expr=\"0.0\" rowrange=\"-\"' 
	print cmd; os.system(cmd)
	hdu       = pyfits.open(outputfits)
	new_array = hdu['EVENTS'].data['DELTA_TIME']
	for i in range(len(new_array)):
		if i % Ndisp == 0: 
			print "...%d/%d (%.1f%%)" % (i,len(time_org),float(i)/float(len(time_org))*100.0)
		new_array[i] = delay_time[i]
	hdu.writeto(outputfits,clobber=True)	

cmd = ''
for i in range(0,3):
	cmd += 'fthedit %s[%d] \@%s \n' % (outputfits,i,headerfile)
	cmd += 'fparkey %.10f "%s[%d]" ORBPDAY comm="Orbital Period (day)" add=yes\n' % (
		param['orb_porb_day'], outputfits, i)
	cmd += 'fparkey %.6f "%s[%d]" ORBT0 comm="Periastron passage (Mission Time s)" add=yes\n' % (
		param['orb_t0_sec'], outputfits, i)			
	cmd += 'fparkey %.5e "%s[%d]" ORBECCEN comm="Eccetricity (no-dimension)" add=yes\n' % (
		param['orb_eccen'], outputfits, i)						
	cmd += 'fparkey %.5f "%s[%d]" ORBASINI comm="Asini (light second) a * sin(i) / c (s)" add=yes\n' % (
		param['orb_asini_sec'], outputfits, i)			
	cmd += 'fparkey %.5f "%s[%d]" ORBWDEG comm="Longitude Periastron (degree)" add=yes\n' % (
		param['orb_w_deg'], outputfits, i)						
cmd += 'rm -f %s' % headerfile
print cmd; os.system(cmd)

quit()

