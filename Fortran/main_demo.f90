program main
    use func_Fluspect_KM, only : nw,nwlf,nwle,Fluspect_KM
    implicit none

    real(8) :: Cab,Cca,V2Z,Cw,Cdm,Cs,Cant,Cbc,Cp,N,fqeI,fqeII
    real(8) :: LRT(nw,10),FLUO_I(nwlf,nwle,2),FLUO_II(nwlf,nwle,2)

    Cab = 50.0
    Cca = 10.0
    V2Z = 0.50
    Cw = 0.010
    Cdm = 0.015
    Cs = 0.0
    Cant = 0.0
    N = 1.6
    fqeI = 0.002
    fqeII = 0.01
    Cbc = 0.0
    Cp = 0.0

    ! Be careful with the usage of Cdm, Cbc+Cp
    call Fluspect_KM(Cab,Cca,V2Z,Cw,Cdm,Cs,Cant,Cbc,Cp,N,fqeI,fqeII,LRT,FLUO_I,FLUO_II)

end program
