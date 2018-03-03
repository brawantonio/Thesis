    program molecular_dynamics
    
    integer nparticles,particles_per_side
    real dt,time,Tmax
    real LJpotential,LJforce,periodic
    
    character*20 :: filename1,filename2,filename3
    character*62 :: fileplace
    
    real position_x(1001),position_y(1001),position_z(1001)
    common /block1/ position_x,position_y,position_z
    real position0_x(1001),position0_y(1001),position0_z(1001)
    common /block2/ position0_x,position0_y,position0_z
    real v_x(1001),v_y(1001),v_z(1001)
    common /block3/ v_x,v_y,v_z
    real f_x(1001),f_y(1001),f_z(1001)
    common /block4/ f_x,f_y,f_z
    real Uen,Ten,Ken,Temp  
    common /block5/ Uen,Ten,Ken,Temp
    real box
    common /block6/ box
    real g(1000),bins(1000)
    integer nhis
    common /block7/ g,bins,nhis
    
    
    fileplace="/home/brawantonio/Documents/Current semester/SMSP/MDP/P8/data/"
    call random_seed()
    
    !Everything is in normalized units
    Temp=0.78 
    particles_per_side=10
    nparticles=particles_per_side**3 
    box=10.73
    
    dt=1e-2
    Tmax=25
    time=0.0
    nhis=70
    
    i=0
    k=1
    call init(particles_per_side,nparticles,dt) 
    open(1,file='energy.dat')
    do while(time .lt. Tmax)
        call force(nparticles)
        call integrate(nparticles,dt)
        write(1,*) time,Uen,Ten,Ken,Temp
        
            if(mod(i,10)==0 .and. time>0.5) then
                call gr(nparticles)
                write(filename1,'("positions",I4,".dat")') k
                write(filename2,'("velocities",I4,".dat")') k
                write(filename3,'("rdf",I4,".dat")') k
                open(2,file=fileplace//filename1)
                open(3,file=fileplace//filename2)
                open(4,file=fileplace//filename3)
                    do j=1,nparticles
                        write(2,*) v_x(j),v_y(j),v_z(j)
                        write(3,*) position_x(j),position_y(j),position_z(j)
                    enddo
                    do j=1,nhis
                        write(4,*) bins(j),g(j)
                    enddo
                close(2)
                close(3)
                close(4)
                k=k+1
            end if
            
        time=time+dt
        i=i+1
    enddo
    close(1)
    
    write(*,*) k-1

    
    
    
    end
    
    
    
    
    
    

    subroutine init(particles_per_side,nparticles,dt)
    real position_x(1001),position_y(1001),position_z(1001)
    common /block1/ position_x,position_y,position_z
    real position0_x(1001),position0_y(1001),position0_z(1001)
    common /block2/ position0_x,position0_y,position0_z
    real v_x(1001),v_y(1001),v_z(1001)
    common /block3/ v_x,v_y,v_z
    real f_x(1001),f_y(1001),f_z(1001)
    common /block4/ f_x,f_y,f_z
    real Uen,Ten,Ken,Temp  
    common /block5/ Uen,Ten,Ken,Temp
    real box
    common /block6/ box
        
    integer counter,particles_per_side,nparticles
    real fs,svx,svy,svz,sv2
        
    svx=0.0
    svy=0.0
    svz=0.0
    sv2=0.0
    
    counter=0
    do i=1,particles_per_side
        do j=1,particles_per_side
            do k=1,particles_per_side
                position_x(k+counter)=((i-1.0)/(particles_per_side-1.0)*box-box/2.0)*0.9 !(rand()-0.5)*box!
                position_y(k+counter)=((j-1.0)/(particles_per_side-1.0)*box-box/2.0)*0.9
                position_z(k+counter)=((k-1.0)/(particles_per_side-1.0)*box-box/2.0)*0.9
                
                v_x(k+counter)=rand()-0.5
                v_y(k+counter)=rand()-0.5
                v_z(k+counter)=rand()-0.5
                
                svx=svx+v_x(k+counter)
                svy=svy+v_y(k+counter)
                svz=svz+v_z(k+counter)
                
                sv2=sv2+v_x(k+counter)**2+v_y(k+counter)**2+v_z(k+counter)**2
            enddo
            counter=counter+particles_per_side
        enddo
    enddo
    
    svx=svx/nparticles
    svy=svy/nparticles
    svz=svz/nparticles
    sv2=sv2/nparticles
    fs=sqrt(3.0*Temp/sv2)
    do i=1,nparticles
        v_x(i)=(v_x(i)-svx)*fs
        v_y(i)=(v_y(i)-svy)*fs
        v_z(i)=(v_z(i)-svz)*fs
        
        position0_x(i)=position_x(i)-v_x(i)*dt
        position0_y(i)=position_y(i)-v_y(i)*dt
        position0_z(i)=position_z(i)-v_z(i)*dt
            
    enddo
    return
    end
    
    subroutine force(nparticles)
    real LJpotential,LJforce
    integer nparticles
        
    real position_x(1001),position_y(1001),position_z(1001)
    common /block1/ position_x,position_y,position_z
    real position0_x(1001),position0_y(1001),position0_z(1001)
    common /block2/ position0_x,position0_y,position0_z
    real v_x(1001),v_y(1001),v_z(1001)
    common /block3/ v_x,v_y,v_z
    real f_x(1001),f_y(1001),f_z(1001)
    common /block4/ f_x,f_y,f_z
    real Uen,Ten,Ken,Temp  
    common /block5/ Uen,Ten,Ken,Temp
    real box
    common /block6/ box
    
    Uen=0.0
    do i=1,nparticles
        f_x(i)=0.0
        f_y(i)=0.0
        f_z(i)=0.0
    enddo

    do i=1,nparticles-1
        do j=i+1,nparticles
            xr=position_x(i)-position_x(j)
            yr=position_y(i)-position_y(j)
            zr=position_z(i)-position_z(j)
                
            xr=xr-box*nint(xr/box)
            yr=yr-box*nint(yr/box)
            zr=zr-box*nint(zr/box)
                
            r2=xr**2.+yr**2.+zr**2.
                
            f_x(i)=LJforce(r2)*xr+f_x(i)
            f_y(i)=LJforce(r2)*yr+f_y(i)
            f_z(i)=LJforce(r2)*zr+f_z(i)
            
            f_x(j)=-LJforce(r2)*xr+f_x(j)
            f_y(j)=-LJforce(r2)*yr+f_y(j)
            f_z(j)=-LJforce(r2)*zr+f_z(j)
            
            Uen=LJpotential(r2)+Uen
        enddo
    enddo
    
    return
    end
    
    subroutine integrate(nparticles,dt)
    integer nparticles
     
    real position_x(1001),position_y(1001),position_z(1001)
    common /block1/ position_x,position_y,position_z
    real position0_x(1001),position0_y(1001),position0_z(1001)
    common /block2/ position0_x,position0_y,position0_z
    real v_x(1001),v_y(1001),v_z(1001)
    common /block3/ v_x,v_y,v_z
    real f_x(1001),f_y(1001),f_z(1001)
    common /block4/ f_x,f_y,f_z
    real Uen,Ten,Ken,Temp  
    common /block5/ Uen,Ten,Ken,Temp
    real box
    common /block6/ box
       
    sv2=0.0
    Ken=0.0
    Ten=0.0
    do i=1,nparticles
        xx=2.0*position_x(i)-position0_x(i)+f_x(i)*dt**2
        yy=2.0*position_y(i)-position0_y(i)+f_y(i)*dt**2
        zz=2.0*position_z(i)-position0_z(i)+f_z(i)*dt**2
        
        v_x(i)=(xx-position0_x(i))/(2.0*dt)
        v_y(i)=(yy-position0_y(i))/(2.0*dt)
        v_z(i)=(zz-position0_z(i))/(2.0*dt)
        
        dx=xx-position_x(i)
        dy=yy-position_y(i)
        dz=zz-position_z(i)
        
        position_x(i)=periodic(xx,box)
        position_y(i)=periodic(yy,box)
        position_z(i)=periodic(zz,box)
        
        position0_x(i)=position_x(i)-dx
        position0_y(i)=position_y(i)-dy
        position0_z(i)=position_z(i)-dz
        
        sv2=sv2+v_x(i)**2+v_y(i)**2+v_z(i)**2
    enddo
    
    Temp=sv2/(3.0*nparticles)
    Ken=0.5*sv2/nparticles
    Uen=Uen/nparticles
    Ten=Uen+Ken
    return
    end
    
    subroutine gr(nparticles)
    real position_x(1001),position_y(1001),position_z(1001)
    common /block1/ position_x,position_y,position_z
    real box
    common /block6/ box
    real g(1000),bins(1000)
    integer nhis
    common /block7/ g,bins,nhis
    
    integer ig
    
    pi=4*atan(1.0)
    rho=nparticles/box**3

    delg=box/(2.*nhis)
    do i=0,nhis
        g(i)=0
    enddo
        
    do i=1,nparticles-1
        do j=i+1,nparticles
            xr=position_x(i)-position_x(j)
            yr=position_y(i)-position_y(j)
            zr=position_z(i)-position_z(j)
                
            xr=xr-box*nint(xr/box)
            yr=yr-box*nint(yr/box)
            zr=zr-box*nint(zr/box)
                
            r=sqrt(xr**2+yr**2+zr**2)
            if(r .le. 0.5*box) then
                ig=int(r/delg)
                g(ig)=g(ig)+2
            end if
        enddo
    enddo
        
    do i=1,nhis
        r=delg*i
        vb=4*pi*r**2*delg !(4./3.)*pi*((r+delg)**3-r**3)
        nid=vb*rho
        g(i)=g(i)/(nid*nparticles)
        bins(i)=r
    enddo
    
    return
    end
    
   

   
    real function LJpotential(r2)
        LJpotential=4.0*(1.0/r2**6-1.0/r2**3)
    end
    
    real function LJforce(r2)
        LJforce=24.0*(2.0/r2**7-1.0/r2**4)
    end
    
    real function periodic(arg,box)
        if(abs(arg)<0.5*box) then
            periodic=arg
        else if (arg .ge. 0.5*box) then
            periodic=-0.5*box+mod(arg,0.5*box)
        else if(arg .le. -0.5*box) then
            periodic=0.5*box+mod(arg,0.5*box)
        end if
    end
    
