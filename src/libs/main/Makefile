# MAKEFILE: src/main/makefile
# Author:   Lee Margetts
# Date:     14/05/20

FC       = gfortran
FFLAGS   = -O2 -Wall -mtune=intel

SOURCE = bandred.f03 bandwidth.f03 banmul.f03 bantmul.f03 beam_ge.f03         \
	 beam_km.f03 beam_mm.f03 beamdis.f03 bee8.f03 beemat.f03              \
  	 bisect.f03 bmat_nonaxi.f03 checon.f03 chobk1.f03 chobk2.f03          \
	 cholin.f03 comred.f03 comsub.f03 contour.f03 cross_product.f03       \
	 deemat.f03 determinant.f03 dismsh.f03 dismsh_ensi.f03 ecmat.f03      \
	 elmat.f03 exc_nods.f03 fkdiag.f03 fmacat.f03 fmdsig.f03 fmkdke.f03   \
	 fmplat.f03 fmrmat.f03 form_s.f03 formaa.f03 formkb.f03 formkc.f03    \
	 formke.f03 formku.f03 formlump.f03 formm.f03 formtb.f03 formupv.f03  \
	 fsparv.f03 gauss_band.f03 getname.f03 glob_to_axial.f03              \
	 glob_to_loc.f03 hinge.f03 interp.f03 invar.f03 invert.f03            \
	 linmul_sky.f03 lnblnk.f03 load_function.f03 loc_to_glob.f03          \
	 main_int.f03 mcdpl.f03 mesh.f03 mesh_ensi.f03 mocouf.f03 mocouq.f03  \
	 norm.f03 num_to_g.f03 pin_jointed.f03 pmsh_ensi.f03 rect_km.f03      \
	 rigid_jointed.f03 rod_km.f03 rod_mm.f03 sample.f03 seep4.f03         \
	 shape_der.f03 shape_fun.f03 solve_band.f03 spabac.f03                \
	 spabac_gauss.f03 sparin.f03 sparin_gauss.f03 stability.f03           \
         stiff10.f03 stiff15.f03 stiff3.f03 stiff4.f03 stiff6.f03 vecmsh.f03  \
         vmdpl.f03 vmflow.f03
         
OBJS = $(SOURCE:%.f03=%.o)

lib: $(OBJS)
	ar -r main.a $(OBJS)
	
%.o: %.f03
	$(FC) -c  $(FFLAGS) $<
	
clean:
	rm -f $(OBJS)
	rm main.a
	rm main.modmake