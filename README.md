# ftime

FITS (HEASOFT) based spectral tools

faddphase.py
- adds a new "PHASE" column based on the PERIOD.
- example:
	../faddphase.py -i ae102019010hxd_0_pin_rep_cl_bary_energy.evt 
	-o ae102019010hxd_0_pin_rep_cl_bary_energy_phase.evt 
	-p 0.033600912903 
	-d 4.20507696e-13 
	-e 227701951.819787 
