!****************************************************************
! Author : Dionisis Pettas
! Date   : 03/12/2020
!****************************************************************

Module MeshGeneration
    Use Formats
    Use FileModule

    private 
    public  :: SalomeMeshGeneration, python


    Type SalomeMeshGeneration
        Character(len=8)                       :: randomId
        Character(len=:), Allocatable          :: output_pfile
        Character(len=:), Allocatable          :: unvfilename
        Type(fFile)                            :: pfile ! python file
        Type(fFile), Dimension(:), Allocatable :: bnd
        contains
        procedure :: writeBoundaryNodes
        procedure :: setUnvFilename
        procedure :: execute
        procedure :: convertToBinaryTecplot
        procedure :: close
    End Type SalomeMeshGeneration


    interface SalomeMeshGeneration
        module procedure constructor
    end interface SalomeMeshGeneration

    Character(len=*), Parameter :: python = "python3.8"
    Character(len=*), Parameter :: salome = "salome -t"

    contains

        Function constructor(filename) Result(this)
            Implicit None 
            Character(len=*), Intent(In) :: filename
            Type(SalomeMeshGeneration)   :: this

            this%pfile = ReadFile(filename)

            ! Check If basic lines exist in file
            if (.not. this%pfile%existInFile("import read_datfile")) Then 
                print*, "[Error] SalomeMeshGeneration : the python file does not contain the 'import read_datfile'"
                print*, "please check"
                stop    
            end if
        End Function constructor

        Subroutine execute(this) 
            Implicit None 
            Class(SalomeMeshGeneration), Intent(InOut) :: this

            Character(len=:), Allocatable              :: pfilename
            Real(8)                                    :: r ! random number from [0,1]
            Integer                                    :: ibnd  
            Character(len=:), Allocatable              :: line
            Character(len=:), Allocatable              :: newline
            Character(len=:), Allocatable              :: bndfilename
            Character(len=:), Allocatable              :: newbndfilename
            integer                                    :: lineid

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! For the uniqueness of the case concatenate a random number at the end of the generated files
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            call random_number(r)       
            r = 1.d0 + r
            this%randomId = toStr(int(1e+7*r))


            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! Create a unique boundary numbers 
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            do ibnd = 1, size(this%bnd)
                bndfilename     = this%bnd(ibnd)%filename
                lineid          = this%pfile%findLineIdThatContains( bndfilename )
                newbndfilename  = concatenate(bndfilename,"_",this%randomId)

                line            = this%pfile%getLineThatContains(bndfilename )
                newline         = replace( line, bndfilename,  newbndfilename)

                ! replace the name of the boundary on the python file 
                call this%pfile%replaceLinewith(lineid, newline)
                
                ! dump the boundary file to the new name
                call this%bnd(ibnd)%dumpToFile(newbndfilename)
            end do

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! dump the salome script to the new name
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            this%output_pfile = concatenate( this%pfile%filename,"_", this%randomId )
            call this%pfile%dumpToFile     ( this%output_pfile )

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! Call Salome Via Terminal
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
            call execute_command_line( concatenate(salome," ", this%output_pfile," ","> /dev/null") )   





            deallocate ( line          )
            deallocate ( newline       )
            deallocate ( bndfilename   )
            deallocate ( newbndfilename)
        End Subroutine execute

        Subroutine convertToBinaryTecplot(this,filename, renameCoords)
            Implicit None
            Class(SalomeMeshGeneration), Intent(In) :: this
            Character(len=*),            Intent(In) :: filename
            Character(len=*),            Intent(In) :: renameCoords

            Integer                                 :: punit

            open(newunit = punit, file = concatenate("tmp.py","_",this%randomId) )
            write(punit,'(a)') 'import unv'
            write(punit,'(a)') replace('f = unv.Reader("*")','*', this%unvfilename)
            write(punit,'(a)') replace('f.toAsciiTecplot("*","remesh",renameCoords = "*")',"*",[filename, renameCoords] )

            call execute_command_line( concatenate(python," ", concatenate("tmp.py","_",this%randomId)) )
            close(punit,status = 'delete')
            call execute_command_line( replace("preplot * * > /dev/null","*", [filename                 , concatenate(filename,"_")] ) )
            call execute_command_line( replace("mv      * *"            ,"*", [concatenate(filename,"_"),                  filename] ) )
        End Subroutine convertToBinaryTecplot



        Subroutine close(this)
            Implicit None 
            Class(SalomeMeshGeneration) :: this
            Integer :: ibnd

            ! remove the generated files
            do ibnd = 1, size(this%bnd)
            call execute_command_line(concatenate("rm -f"," ",concatenate(this%bnd(ibnd)%filename,"_",this%randomId) ) )
            end do
            call execute_command_line(concatenate("rm -f"," ",concatenate(this%pfile%filename    ,"_",this%randomId) ) )

            do ibnd = 1, size(this%bnd)
                deallocate( this%bnd(ibnd)%filename )
                deallocate( this%bnd(ibnd)%line     )
                this%bnd(ibnd)%nrows = 0
            end do

            deallocate( this%pfile%filename )
            deallocate( this%pfile%line     )
            this%pfile%nrows = 0

            deallocate( this%bnd          )
            deallocate( this%output_pfile )
            deallocate( this%unvfilename  )

            call execute_command_line("killall omniNames")
            call execute_command_line("killall SALOME_ConnectionManagerServer")
            call execute_command_line("killall SALOME_Session_Server")
            call execute_command_line("killall SALOME_LauncherServer")
            call execute_command_line("killall SALOME_Container")
            call execute_command_line("killall SALOME_ModuleCatalog_Server")
            call execute_command_line("killall SALOME_Registry_Server") 
            call execute_command_line("killall SALOMEDS_Server")

            call execute_command_line("rm -rf /tmp/.salome/.PortManager.lock")
            call execute_command_line("rm -rf /tmp/.salome_PortManager.cfg")
            call execute_command_line("rm -rf /tmp/.omni*")

        End Subroutine close



        Subroutine setUnvFilename(this, unvname) 
            Implicit None 
            Class(SalomeMeshGeneration)  :: this
            Character(len=*), Intent(In) :: unvname

            Integer                      :: exportUnvLineId
            Character(len=:),Allocatable :: OldLine
            Character(len=:),Allocatable :: NewLine 
            Integer                      :: left_par
            Integer                      :: right_par

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            exportUnvLineId = this%pfile%findLineIdThatContains("ExportUNV")

            if (exportUnvLineId == -1) Then
                print*, "[Error] SalomeMeshGeneration: the line does not exist in file"
                stop
            end if 
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! Store the unvfilename in the type
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            this%unvfilename = unvname
            

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! Replace the name of the unvfilename in the pfile
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            OldLine     = this%pfile%getLine(exportUnvLineId)
            ! Create the new replacement
            ! the pattern is like 
            ! Mesh_1.ExportUNV("Test.unv")
            ! so the unvname should be placed between the parentheses 

            left_par = index(OldLine,"(")
            right_par= index(OldLine,")")

            NewLine  = concatenate( OldLine(1:left_par), concatenate('"',unvname,'"') , OldLine(right_par:) )

        call this%pfile%replaceLinewith( exportUnvLineId, NewLine)

            deallocate(OldLine)
            deallocate(NewLine)
        End Subroutine setUnvFilename



        Subroutine writeBoundaryNodes(this, filename, x, y)
            Implicit None 
            Class(SalomeMeshGeneration)             :: this
            Character(len=*)           , Intent(In) :: filename
            Real(8), Dimension(:)      , Intent(In) :: x
            Real(8), Dimension(:)      , Intent(In) :: y

            Type(fFile), Dimension(:), Allocatable  :: tmp
            Integer                                 :: nnodes
            Integer                                 :: i
            Type(fFile)                             :: bnd

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! Check if the input arrays have the same size
            if (size(x) /= size(y) ) Then 
                print*, "[Error] SalomeMeshGeneration : the two entities have different size."
                stop
            end if
            
            ! Check If filename exist in file
            ! if (.not. this%pfile%existInFile(filename)) Then 
                ! print*, "[Error] SalomeMeshGeneration : the python file does not contain the 'boundary nodes' file"
                ! print*, "please check"
                ! stop    
            ! end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
            nnodes = size(x)
            bnd    = OpenFile(filename) 

            ! Write the values in the fFile Object
            do i = 1, nnodes
                call bnd%push_back( replace("{} {}","{}", [toStr(x(i),'f20.12'), toStr(y(i),'f20.12')]) )
            end do


            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            ! push back the values on boundaries on bnd
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            if ( .not. allocated(this%bnd) ) Then 
                allocate(this%bnd(1) )
                this%bnd(1) = bnd
            else 
                ! tmp array 
                allocate( tmp( size(this%bnd)    ) )
                tmp(1:size(tmp)) = this%bnd

                deallocate( this%bnd               )
                allocate  ( this%bnd(size(tmp)+1)  )
                this%bnd(1:size(tmp)  ) = tmp
                this%bnd(  size(tmp)+1) = bnd
                deallocate(tmp)
            end if  
        End Subroutine writeBoundaryNodes




End Module MeshGeneration
