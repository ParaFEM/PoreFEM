# MAKEFILE: src/geom/makefile
# SYSTEM:   Windows 10 with MSYS2
# Author:   Lee Margetts
# Date:     14/05/20

FC       = gfortran
FFLAGS   = -O2 -Wall -mtune=intel

SOURCE = beta.f95 chkne3.f95 chknxe.f95 chol2.f95 corn2d.f95 corn3d.f95      \
	daxpy.f95 dchol2.f95 dcvaa1.f95 dcvaa2.f95 dcvaa3.f95 dcvab1.f95     \
 	dcvab2.f95 dcvab3.f95 dcvit1.f95 dcvit2.f95 dcvit3.f95 dcvmt2.f95    \
 	dcvmt3.f95 ddot.f95 derf.f95 dismsh2f.f95 dismsh3f.f95 dismsh4f.f95  \
 	dlace1.f95 dlafr1.f95 dlafr2.f95 dlafs3.f95 dlavx1.f95 dlavx2.f95    \
 	dlavx22.f95 dlavx3.f95 dlsep2.f95 dlsep3.f95 dlsfr2.f95 dlsfr3.f95   \
 	dlsmp1.f95 dlspx1.f95 dlspx2.f95 dlspx3.f95 dot3d.f95 dot3db.f95     \
 	dphi.f95 dphinv.f95 dsifa.f95 dsisl.f95 dswap.f95 edge3d.f95         \
 	exca.f95 exca2.f95 fem2rf.f95 fem2rf2.f95 fem3rf.f95 gaf95_int.f95   \
 	gamln.f95 gammaf.f95 gausv.f95 genfld.f95 genfld2.f95 genprp.f95     \
 	getfsp.f95 gmdst.f95 idamax.f95 intr2d.f95 intr3d.f95 iseed.f95      \
 	iseed2.f95 las1g.f95 las1i.f95 las2g.f95 las2g2.f95 las2i.f95        \
 	las2i2.f95 las3g.f95 las3i.f95 ldigit.f95 lnblnk.f95 meshWT.f95      \
 	msla.f95 ntwom.f95 permprop.f95 phi.f95 plan3d.f95 pltfld.f95        \
 	pltflds.f95 prfmte.f95 prfmtf.f95 prfmti.f95 print1.f95 print2.f95   \
 	print3.f95 print4.f95 printv.f95 propst.f95 randu.f95 randu2.f95     \
 	rect.f95 rect2.f95 rect2mod.f95 rect2mod1.f95 second.f95             \
 	side2d.f95 side3d.f95 sim1dm.f95 sim1dm1.f95 sim1dm2.f95 sim2bc1.f95 \
 	sim2bc2.f95 sim2dk.f95 sim2sd_2.f95 sim2sd1.f95 sim2sd2.f95          \
  	simpl3.f95 simpl3cfld.f95 thin1d.f95 uservf.f95 vnorm.f95            \
 	vnorm2.f95
         
OBJS = $(SOURCE:%.f95=%.o)

lib: $(OBJS)
	ar -r gaf95.a $(OBJS)
	
%.o: %.f95
	$(FC) -c  $(FFLAGS) $<
	
clean:
	rm -f $(OBJS)
	rm gaf95.a
	rm gaf95.mod