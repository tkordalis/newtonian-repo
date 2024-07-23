#---------------------------------------------------------------------------

NEXE="exe"
echo $NEXE

#----------------------------------------------------------------------------


OPT="-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -gen-interfaces -debug all -implicitnone -fstack-protector -warn none"


echo $OPT

#----------------------------------------------------------------------------
# gia na doume
MKL="-qmkl=parallel"

#----------------------------------------------------------------------------

rm -f $NEXE
rm -f nohup.out

								# ./src/Remesh/RemeshProcedure.f90\
 ifort $OPT -o $NEXE                           \
								./src/utilities/check_for_floating_point_exceptions.f90 \
								./src/utilities/system_tools.f90 \
								./src/utilities/formats.f90    \
				 				./src/utilities/arraytools.f90 \
								./src/utilities/FileModule.f90 \
								./src/utilities/geometry.f90   \
								./src/Remesh/MeshGeneration.f90\
				        ./src/export/unv.f90           \
								./src/export/WriteTecplot.f90  \
								./src/export/ReadTecplot.f90   \
								./src/Interpolation/TecplotInterpolation.f90\
								./src/export/tecplot.f90       \
	              ./src/storage.f90            \
	              ./src/Fem2D_mod.f90            \
								./src/Bulk_Equations.f90       \
							  ./src/extraequations.f90       \
								./src/NumericalExtraJacobian.f90\
		            ./src/Boundary_Equations.f90     	 \
								./src/BoundaryEquations/FixWallBoundary.f90 \
								./src/BoundaryEquations/MovingWallBoundary.f90 \
								./src/FieldFunctions.f90\
								./src/io_module.f90\
							./src/InitializeTypes.f90\
							./src/former_external_subroutines.f90\
                     Fem2D_prg.f90            \
	             $MKL -qopenmp -L./src/export/TECLIB/lib/ -ltecio -lstdc++
 
#----------------------------------------------------------------------------
ctags -R .
./Clear.sh

gdb ./$NEXE
#----------------------------------------------------------------------------
