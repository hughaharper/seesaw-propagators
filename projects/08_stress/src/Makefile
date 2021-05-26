#	Makefile for topo_stress 

GMT     = -L$(GMTHOME)/lib
CDF     = -L$(NETCDFHOME)/lib
FFT     = ../lib/libfftfault.a
FFLAGS =  -O3 -m64 -Bstatic
FC = gfortran

all: topo_stress principal_stress eqrockimg spheregrav eqrockgrd flexure2d predictgrav

topo_stress:	topo_stress.o imgfftio.o readsrtmmask.o
	$(FC) $(FFLAGS) topo_stress.o imgfftio.o readsrtmmask.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/topo_stress

principal_stress:	principal_stress.o imgfftio.o 
	$(FC) $(FFLAGS) principal_stress.o imgfftio.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/principal_stress

eqrockimg:	eqrockimg.o imgfftio.o readsrtmmask.o
	$(FC) $(FFLAGS) eqrockimg.o imgfftio.o readsrtmmask.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/eqrockimg

spheregrav:	spheregrav.o imgfftio.o
	$(FC) $(FFLAGS) spheregrav.o imgfftio.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/spheregrav

eqrockgrd:	eqrockgrd.o readsrtmmask.o
	$(FC) $(FFLAGS) eqrockgrd.o readsrtmmask.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/eqrockgrd

flexure2d:	flexure2d.o imgfftio.o readsrtmmask.o
	$(FC) $(FFLAGS) flexure2d.o imgfftio.o readsrtmmask.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/flexure2d

predictgrav:	predictgrav.o
	$(FC) $(FFLAGS) predictgrav.o \
		$(FFT) $(CDF) -lnetcdf $(GMT) -lgmt -o ../bin/predictgrav

principal: principal.o

