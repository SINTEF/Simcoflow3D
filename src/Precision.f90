Module PrecisionVar
    Implicit none
    integer,parameter:: sp = selected_real_kind(6,35)   ! 6 digit after ","  and form 10^-35 to 10^35
    integer,parameter:: dp = selected_real_kind(15,307) ! 15 digit after "," and from 10^-307 to 10^307
    integer,parameter:: it1b = selected_int_kind(2) ! from -10^2 to 10^2
    integer,parameter:: it2b = selected_int_kind(3) ! from -10^3 to 10^3
    integer,parameter:: it4b = selected_int_kind(6) ! from -10^6 to 10^6
    integer,parameter:: it8b = selected_int_kind(12) ! from -10^12 to 10^12
End module PrecisionVar
