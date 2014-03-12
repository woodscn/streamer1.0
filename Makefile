execname = 2_Db
compiler = gfortran

objects = $(execname).o \
	physical_data.o \
	window_data.o \
	global_data.o \
	primary_variables.f95 \
	vars.o \
	init.o \
	write_files_mod.o \
	geom_update_mod.o \
	h_update_mod.o \
	muscl_mod.o \
	riemann.o \
	flow_update_mod.o \
	physical_fluxes.o \
	boundary_conditions_mod.o
	
o_compile = $(compiler) -O0 -fdefault-real-8 -c -fbounds-check -g3
	
all : $(execname)
$(execname) : $(objects)
	$(compiler) -O3 -fdefault-real-8 -fbounds-check $(objects)
	mv a.out 2_Db
	
boundary_conditions_mod.o : boundary_conditions_mod.f95 vars.o
	$(o_compile) boundary_conditions_mod.f95

physical_data.o : physical_data.f95
	$(o_compile) physical_data.f95
	
physical_fluxes.o : physical_fluxes.f95
	$(o_compile) physical_fluxes.f95

window_data.o : window_data.f95
	$(o_compile) window_data.f95

global_data.o : global_data.f95
	$(o_compile) global_data.f95

primary_variables.o : primary_variables.f95
	$(o_compile) primary_variables.f95

vars.o : vars.f95 physical_data.o
	$(o_compile) vars.f95

init.o : init.f95 global_data.o primary_variables.o window_data.o boundary_conditions_mod.o
	$(o_compile) init.f95

write_files_mod.o : write_files_mod.f95 vars.o
	$(o_compile) write_files_mod.f95

geom_update_mod.o : geom_update_mod.f95 global_data.o primary_variables.o window_data.o boundary_conditions_mod.o
	$(o_compile) geom_update_mod.f95

h_update_mod.o : h_update_mod.f95 #boundary_conditions_mod.o
	$(o_compile) h_update_mod.f95

muscl_mod.o : muscl_mod.f95
	$(o_compile) muscl_mod.f95

#riemann.o : riemann2.f95
#	$(o_compile) riemann2.f95
#	mv riemann2.o   riemann.o
riemann.o : riemann.f95
	$(o_compile) riemann.f95

#flow_update_mod.o : flow_update_mod.f95  vars.o riemann.o muscl_mod.o physical_fluxes.o boundary_conditions_mod.o
#	$(o_compile) flow_update_mod.f95


flow_update_mod.o : flow_update_mod2.f95  vars.o riemann.o muscl_mod.o physical_fluxes.o boundary_conditions_mod.o
	$(o_compile) flow_update_mod2.f95
	mv flow_update_mod2.o   flow_update_mod.o

$(execname).o : $(execname).f95 global_data.o primary_variables.o init.o write_files_mod.o geom_update_mod.o flow_update_mod.o h_update_mod.o vars.o physical_data.o #boundary_data.o
	$(o_compile) 2_Db.f95


clean : 
	rm *.o
	rm *.mod
	rm $(execname)