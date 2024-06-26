
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!   Author : Dionisis Pettas
! Last Modification : 3/12/2020
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    Module Formats

        private 
        public  :: pathJoin, getFilenameFromPath, concatenate, replace, replace_abstract, split, remove, &
                     toStr, toDouble, toInt, dataFormatColumns , &
                             createFolder, containFolders, reverse

        ! replaces text that is contained in ... (elipsis)
        ! e.g. str = replace_abstract("writeFile('/home/user/file')","writeFile(...)",'"file.dat"')
        ! print
        ! writeFile("file.dat")
        interface replace_abstract
            module procedure replace_abstract_simple
        end interface replace_abstract



        interface replace
            module procedure replace_str_simple
            module procedure replace_str_multiple
        end interface replace

        interface remove 
            module procedure remove_simple
            module procedure remove_multiple
        end interface remove

        ! Converts an int or double to String 
        ! toStr(15.d0) or toStr(15.d0,'f12.4')
        interface toStr
            module procedure toStr_int
            module procedure toStr_int_fmt
            module procedure toStr_double
            module procedure toStr_double_fmt
        end interface toStr

        interface toDouble
            module procedure toDouble_simple
        end interface toDouble

        interface toInt
            module procedure toInt_simple
        end interface toInt

        interface concatenate
            module procedure concatenate_str_2
            module procedure concatenate_str_3
            module procedure concatenate_str_4
            module procedure concatenate_str_5
            module procedure concatenate_str_6
            module procedure concatenate_str_7
        end interface concatenate

        interface pathJoin
            module procedure pathJoin_simple
        end interface pathJoin

        contains 

            Function getFilenameFromPath(path) Result(filename)
                Implicit None 
                Character(len=*), Intent(In) :: path 
                Character(len=:), Allocatable:: filename
            
                Character(len=500), Dimension(:), Allocatable :: path_
                Character(len=:), Allocatable:: str
                
                path_ = split(path, sep = "/")
                filename = trim(adjustl(path_(size(path_))))
            End Function getFilenameFromPath

            !------------------------------------------------------
            ! Reverse a string
            !------------------------------------------------------
            Function reverse(str) Result(output)
                Implicit None 
                Character(len=*)             , Intent(In) :: str
                Character(len=:), Allocatable             :: output
                Integer                                   :: i
                Integer                                   :: n

                n = len(str)
                allocate(output, mold = str)

                do i = 0, n - 1
                    output( n-i:n-i ) = str( i+1:i+1)
                end do 
            End Function reverse

            !------------------------------------------------------
            ! Formats the results in a specific format
            !------------------------------------------------------
            Function dataFormatColumns(fmt) Result(fmt_)
                Implicit None 
                Character(len=*), Intent(In)  :: fmt
                Character(len=:), Allocatable :: fmt_
            
                fmt_ = replace('(*(__,1x))','__', fmt)
            End Function dataFormatColumns


            !---------------------------------------------------------------------
            ! Return a copy with all occurrences of substring old replaced by new. 
            !----------------------------------------------------------------------
            
            Function replace_str_simple(str, old, new) Result(output)
                ! Replaces a substring in a str by the substring new
                Implicit None 
                Character(len=*), Intent(In)  :: str
                Character(len=*), Intent(In)  :: old 
                Character(len=*), Intent(In)  :: new 

                Character(len=:), Allocatable :: output

                Integer                       :: idx
                Integer                       :: len_

                idx  = index(str, old)
                len_ = len(old)

                if ( idx == 0 ) then ; output = str
                else                 ; output = str(1:idx-1)//new//str(idx+len_:) 
                end if  
            End Function replace_str_simple


            Function replace_str_multiple(str, old, new) Result(output)
                ! Replaces a substring in a str by the substring new
                Implicit None 
                Character(len=*),               Intent(In)  :: str
                Character(len=*),               Intent(In)  :: old 
                Character(len=*), Dimension(:), Intent(In)  :: new 
            
                Integer                                     :: i    
                Character(len=:), Allocatable               :: str_new
                Character(len=:), Allocatable               :: str_old
                Character(len=:), Allocatable               :: output

                str_old = str
                do i = 1, size(new)
                    
                    str_new  = replace( str_old, old, trim(adjustl(new(i))) )
                    if (str_new == str_old ) exit
                    str_old  = str_new
                end do  
                output = str_new
            End Function replace_str_multiple


            Function replace_abstract_simple(str, abc, new) Result(output)
                ! abc must contains the ellipsis symbol "..."
                Implicit None 
                Character(len=*), Intent(In) :: str
                Character(len=*), Intent(In) :: abc 
                Character(len=*), Intent(In) :: new 
                Character(len=:), Allocatable:: output

                Character(len=500), Dimension(:), Allocatable :: words  
                Character(len=:)                , Allocatable :: tmp

                words = split(abc, sep = "...") ! this is always 2 (left , right)

                ! Check if substrings are the same
                ! Left substring
                tmp = trim(words(1))
                if ( str(1:len(tmp)) /= tmp  ) Then 
                    print*, "[Error] replace_abstract: The abstract str is not compatible with the str."
                    print*, replace("str          : *", "*", str)
                    print*, replace("abstract str : *", "*", abc)
                    stop
                end if 

                ! Right Substring
                tmp = trim(words(2))
                if ( str(len(str) - len(tmp)+1:len(str)) /= tmp  ) Then 
                    print*, "[Error] replace_abstract: The abstract str is not compatible with the str."
                    print*, replace("str          : *", "*", str)
                    print*, replace("abstract str : *", "*", abc)
                    stop
                end if 



                output = concatenate(trim(words(1)), new, trim(words(2)) )

                deallocate(tmp)
            End Function replace_abstract_simple





            !----------------------------------------------------------------------
            ! Split String according to separator
            !----------------------------------------------------------------------
            Function split(str, sep) Result(output)
                Implicit None 
                Character(len=*)  , Intent(In)                :: str
                Character(len=*)  , Intent(In)                :: sep
                Character(len=500), Dimension(:), Allocatable :: output
    
                Character(len=500), Dimension(:), Allocatable :: store  
                Character(len=:), Allocatable                 :: str_   
                Character(len=:), Allocatable                 :: itm
                Integer                                       :: isub
                Integer                                       :: fsub
                Integer                                       :: itm_counter    
                str_   = str    

                ! If starts with the sep / remove sep
                if ( str_(1:len(sep)) == sep) str_ = replace(str,sep,"")


                allocate( store ( 100 ) )

                itm_counter = 0
                do while (.True.)
                    isub = 1
                    fsub = index(str_,sep)
                    if (fsub == 0) exit
                    itm  = str_(isub:fsub-1) ! exclude the sep
                    str_ = replace(str_, itm,"")
                    str_ = trim(adjustl(str_(len(sep)+1:len(str_))))
                    ! Store 
                    itm_counter = itm_counter + 1
                    if (itm_counter > 100) then; print*, "[Error] : split. The items exceed the limit of 100 please update the routine."; stop
                    end if

                    store(itm_counter) = itm
                end do
                itm_counter = itm_counter + 1
                store(itm_counter) = str_
                ! return
                output = store(1:itm_counter)
                deallocate(store)
            End Function split

            
            Function toStr_int_fmt( val, fmt ) Result(output) 
                Implicit None 
                Integer         , Intent(In)  :: val 
                Character(len=*), Intent(In)  :: fmt
                Character(len=:), Allocatable :: output
                

                Integer, Parameter            :: nbuf = 100
                Character(len=nbuf)           :: buffer 

                write(buffer, replace("(*)","*",fmt) ) val
                output = trim(adjustl(buffer) )

            End Function toStr_int_fmt

            Function toStr_int( val ) Result(output) 
                Implicit None 
                Integer         , Intent(In)  :: val 
                Character(len=:), Allocatable :: output
                

                Integer, Parameter            :: nbuf = 100
                Character(len=nbuf)           :: buffer 

                write(buffer, '(i)' ) val

                output = trim(adjustl(buffer) )

            End Function toStr_int

            Function toStr_double_fmt( val, fmt ) Result(output) 
                Implicit None 
                Real(8)         , Intent(In)  :: val 
                Character(len=*), Intent(In)  :: fmt
                Character(len=:), Allocatable :: output
                

                Integer, Parameter            :: nbuf = 100
                Character(len=nbuf)           :: buffer 

                write(buffer, replace("(*)","*",fmt) ) val
                output = trim(adjustl(buffer) )

            End Function toStr_double_fmt

            Function toStr_double( val ) Result(output) 
                Implicit None 
                Real(8)         , Intent(In)  :: val 
                Character(len=:), Allocatable :: output
                

                Integer, Parameter            :: nbuf = 100
                Character(len=nbuf)           :: buffer 

                write(buffer, '(f12.4)' ) val

                output = trim(adjustl(buffer) )

            End Function toStr_double



            Function toDouble_simple(str) Result(output)
                Implicit None 
                Character(len=*), Intent(In) :: str
                Real(8)                      :: output
                Character(len=:), Allocatable:: str_

                str_ = trim(adjustl(str))

                Read(str_,*) output
            End Function toDouble_simple

            Function toInt_simple(str) Result(output)
                Implicit None 
                Character(len=*), Intent(In) :: str
                Integer                      :: output
                Character(len=:), Allocatable:: str_

                str_ = trim(adjustl(str))
                Read(str_,*) output
            End Function toInt_simple

            Function remove_simple(str, sub) Result(output)
                Implicit None 
                Character(len=*), Intent(In) :: str
                Character(len=*), Intent(In) :: sub 
                Character(len=:), Allocatable:: output

                output = replace( str, sub, "")
            End Function remove_simple


            Function remove_multiple(str, sub) Result(output)
                Implicit None 
                Character(len=*),               Intent(In) :: str
                Character(len=*), Dimension(:), Intent(In) :: sub 
                Character(len=:), Allocatable              :: output

                Integer                                    :: i
                Character(len=:), Allocatable              :: str_

                str_ = str
                do i = 1, size(sub)
                    str_ = replace(str_, trim(adjustl(sub(i))),"")
                end do
                output = str_
            End Function remove_multiple



            Function containFolders(filename) Result(output)
                Implicit None 
                Character(len=*), Intent(In) :: filename 
                Logical                      :: output

                output = index(filename, '/') /= 0

            End Function containFolders

            Subroutine createFolder(filename)
                ! The path is the total string minus the name of the file. 
                ! e.g. /home/user/folder/data/filename.dat
                !      ---------------------- ------------
                !         folder                file
                Implicit None 
                Character(len=*), Intent(In)  :: filename

                Character(len=:), Allocatable :: rfilename
                Integer                       :: idx

                rfilename = reverse(filename)   
                idx       = index( rfilename, '/')
                idx       = len(filename) -  idx

                call execute_command_line("mkdir -p "// filename(1:idx) )
            End Subroutine createFolder


            Function concatenate_str_2 (str1, str2) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                                
                output = str1//str2

            End Function concatenate_str_2

            Function concatenate_str_3 (str1, str2, str3) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                Character(len=*), Intent(In) :: str3
                                
                output = str1//str2//str3

            End Function concatenate_str_3

            Function concatenate_str_4 (str1, str2, str3, str4) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                Character(len=*), Intent(In) :: str3
                Character(len=*), Intent(In) :: str4
                                
                output = str1//str2//str3//str4

            End Function concatenate_str_4

            Function concatenate_str_5 (str1, str2, str3, str4, str5) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                Character(len=*), Intent(In) :: str3
                Character(len=*), Intent(In) :: str4
                Character(len=*), Intent(In) :: str5
                                
                output = str1//str2//str3//str4//str5

            End Function concatenate_str_5

            Function concatenate_str_6 (str1, str2, str3, str4, str5, str6) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                Character(len=*), Intent(In) :: str3
                Character(len=*), Intent(In) :: str4
                Character(len=*), Intent(In) :: str5
                Character(len=*), Intent(In) :: str6
                                
                output = str1//str2//str3//str4//str5//str6

            End Function concatenate_str_6


            Function concatenate_str_7 (str1, str2, str3, str4, str5, str6, str7) Result(output)
                Implicit None 
                Character(len=:), Allocatable:: output
                Character(len=*), Intent(In) :: str1
                Character(len=*), Intent(In) :: str2
                Character(len=*), Intent(In) :: str3
                Character(len=*), Intent(In) :: str4
                Character(len=*), Intent(In) :: str5
                Character(len=*), Intent(In) :: str6
                Character(len=*), Intent(In) :: str7
                                
                output = str1//str2//str3//str4//str5//str6//str7

            End Function concatenate_str_7

            Function pathJoin_simple(str1, str2) Result(output)
                Implicit None 
                Character(len=*), Intent(In)  :: str1
                Character(len=*), Intent(In)  :: str2
                Character(len=:), Allocatable :: output

                Character(len=1)              :: lastCharacter
                
                lastCharacter = str1(len(str1):len(str1))

                if (lastCharacter == "/") then ; output = str1//str2
                else                           ; output = str1//"/"//str2
                end if 
            End Function pathJoin_simple

    End Module Formats
