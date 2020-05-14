# MAKEFILE: src/geom/makefile
# Author:   Lee Margetts
# Date:     14/05/20

FC       = gfortran
FFLAGS   = -O2 -Wall -mtune=intel

SOURCE = bc_rect.f03 emb_2d_bc.f03 emb_2d_geom.f03 emb_3d_bc.f03              \
         fmcoem.f03 fmglem.f03 formnf.f03 geom_freesurf.f03                   \
         geom_int.f03 geom_rect.f03 hexahedron_xz.f03 mesh_size.f03
         
OBJS = $(SOURCE:%.f03=%.o)

lib: $(OBJS)
	ar -r geom.a $(OBJS)
	
%.o: %.f03
	$(FC) -c  $(FFLAGS) $<
	
clean:
	rm -f $(OBJS)
	rm geom.a
	rm geom.mod