module nhmg_debug

  use mg_mpi
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_in(nx,ny,nz,Hzba,Hza,Hzha,uba,vba,wba,ua,va,wa,rua,rva,rwa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hza
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzha
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: uba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: vba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: rua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: rva
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: rwa

    real(kind=rp), dimension(:,:,:), pointer :: Hzb,Hz,Hzh,ub,vb,wb,u,v,w,ru,rv,rw

    integer(kind=ip), save :: iter_write_pred_in=0
    iter_write_pred_in = iter_write_pred_in + 1

    Hzb => Hzba
    Hz => Hza
    Hzh => Hzha
    ub => uba
    vb => vba
    wb => wba
    u => ua
    v => va
    w => wa
    ru => rua
    rv => rva
    rw => rwa

    call write_netcdf(Hzb,vname='Hzb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hz,vname='Hz',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hzh,vname='Hzh',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(ub,vname='ub',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(vb,vname='vb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(wb,vname='wb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(ru,vname='ru',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(rv,vname='rv',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(rw,vname='rw',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)

  end subroutine nhmg_write_pred_in

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_inter1(nx,ny,nz,ua,va,wa,rufrca,rvfrca)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rufrca
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rvfrca

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:),   pointer :: rufrc,rvfrc

    integer(kind=ip), save :: iter_write_pred_inter1=0
    iter_write_pred_inter1 = iter_write_pred_inter1 + 1

    u => ua
    v => va
    w => wa
    rufrc => rufrca
    rvfrc => rvfrca

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_inter1.nc',rank=myrank,iter=iter_write_pred_inter1)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_inter1.nc',rank=myrank,iter=iter_write_pred_inter1)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_inter1.nc',rank=myrank,iter=iter_write_pred_inter1)
    call write_netcdf(rufrc,vname='rufrc',netcdf_file_name='wrt_pred_inter1.nc',rank=myrank,iter=iter_write_pred_inter1)
    call write_netcdf(rvfrc,vname='rvfrc',netcdf_file_name='wrt_pred_inter1.nc',rank=myrank,iter=iter_write_pred_inter1)

  end subroutine nhmg_write_pred_inter1

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_inter2(nx,ny,nz,ua,va,wa,rufrca,rvfrca)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rufrca
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rvfrca

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:),   pointer :: rufrc,rvfrc

    integer(kind=ip), save :: iter_write_pred_inter2=0
    iter_write_pred_inter2 = iter_write_pred_inter2 + 1

    u => ua
    v => va
    w => wa
    rufrc => rufrca
    rvfrc => rvfrca

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_inter2.nc',rank=myrank,iter=iter_write_pred_inter2)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_inter2.nc',rank=myrank,iter=iter_write_pred_inter2)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_inter2.nc',rank=myrank,iter=iter_write_pred_inter2)
    call write_netcdf(rufrc,vname='rufrc',netcdf_file_name='wrt_pred_inter2.nc',rank=myrank,iter=iter_write_pred_inter2)
    call write_netcdf(rvfrc,vname='rvfrc',netcdf_file_name='wrt_pred_inter2.nc',rank=myrank,iter=iter_write_pred_inter2)

  end subroutine nhmg_write_pred_inter2

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_2d(nx,ny,zetaa,du1a,dv1a,du2a,dv2a)

    integer(kind=ip), intent(in) :: nx, ny

    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: zetaa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: du1a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: dv1a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: du2a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: dv2a

    real(kind=rp), dimension(:,:), pointer :: zeta,du1,dv1,du2,dv2

    integer(kind=ip), save :: iter_write_pred_2d=0
    iter_write_pred_2d = iter_write_pred_2d + 1

    zeta => zetaa
    du1 => du1a
    dv1 => dv1a
    du2 => du2a
    dv2 => dv2a

    call write_netcdf(zeta,vname='zeta',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(du1,vname='du1',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(dv1,vname='dv1',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(du2,vname='du2',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(dv2,vname='dv2',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)

  end subroutine nhmg_write_pred_2d

  !--------------------------------------------------------------
  subroutine nhmg_write_pred_out(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_pred_out=0
    iter_write_pred_out = iter_write_pred_out + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)

  end subroutine nhmg_write_pred_out

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_in(nx,ny,nz,ua,va,wa,rua,rva,rwa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: rua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: rva
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: rwa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w,ru,rv,rw

    integer(kind=ip), save :: iter_write_corr_in=0
    iter_write_corr_in = iter_write_corr_in + 1

    u => ua
    v => va
    w => wa
    ru => rua
    rv => rva
    rw => rwa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(ru,vname='ru',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(rv,vname='rv',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(rw,vname='rw',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)

  end subroutine nhmg_write_corr_in

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_inter1(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_corr_inter1=0
    iter_write_corr_inter1 = iter_write_corr_inter1 + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_inter1.nc',rank=myrank,iter=iter_write_corr_inter1)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_inter1.nc',rank=myrank,iter=iter_write_corr_inter1)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_inter1.nc',rank=myrank,iter=iter_write_corr_inter1)

  end subroutine nhmg_write_corr_inter1

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_inter2(nx,ny,nz,ua,va,wa,hza,du_avga,dv_avga)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: hza
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: du_avga
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: dv_avga

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w
    real(kind=rp), dimension(:,:,:), pointer :: hz
    real(kind=rp), dimension(:,:),   pointer :: du_avg,dv_avg

    integer(kind=ip), save :: iter_write_corr_inter2=0
    iter_write_corr_inter2 = iter_write_corr_inter2 + 1

    u => ua
    v => va
    w => wa
    hz => hza
    du_avg => du_avga
    dv_avg => dv_avga

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)
    call write_netcdf(hz,vname='hz',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)
    call write_netcdf(du_avg,vname='du',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)
    call write_netcdf(dv_avg,vname='dv',netcdf_file_name='wrt_corr_inter2.nc',rank=myrank,iter=iter_write_corr_inter2)

  end subroutine nhmg_write_corr_inter2

  !--------------------------------------------------------------
  subroutine nhmg_write_corr_out(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz+1), target, intent(in) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w

    integer(kind=ip), save :: iter_write_corr_out=0
    iter_write_corr_out = iter_write_corr_out + 1

    u => ua
    v => va
    w => wa

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(w,vname='w',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)

  end subroutine nhmg_write_corr_out


  !--------------------------------------------------------------
  subroutine h_write_pred_in(nx,ny,nz,Hzba,Hza,Hzha,uba,vba,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hza
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: Hzha
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: uba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: vba
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: Hzb,Hz,Hzh,ub,vb,u,v

    integer(kind=ip), save :: iter_write_pred_in=0
    iter_write_pred_in = iter_write_pred_in + 1

    Hzb => Hzba
    Hz => Hza
    Hzh => Hzha
    ub => uba
    vb => vba
    u => ua
    v => va

    call write_netcdf(Hzb,vname='Hzb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hz,vname='Hz',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(Hzh,vname='Hzh',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(ub,vname='ub',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(vb,vname='vb',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_in.nc',rank=myrank,iter=iter_write_pred_in)

  end subroutine h_write_pred_in

  !--------------------------------------------------------------
  subroutine h_write_pred_inter(nx,ny,nz,ua,va,rua,rva)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2),        target, intent(in) :: rva

    real(kind=rp), dimension(:,:,:), pointer :: u,v
    real(kind=rp), dimension(:,:),   pointer :: ru,rv

    integer(kind=ip), save :: iter_write_pred_inter=0
    iter_write_pred_inter = iter_write_pred_inter + 1

    u => ua
    v => va
    ru => rua
    rv => rva

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(ru,vname='ru',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)
    call write_netcdf(rv,vname='rv',netcdf_file_name='wrt_pred_inter.nc',rank=myrank,iter=iter_write_pred_inter)

  end subroutine h_write_pred_inter

  !--------------------------------------------------------------
  subroutine h_write_pred_2d(nx,ny,zetaa,du1a,dv1a,du2a,dv2a)

    integer(kind=ip), intent(in) :: nx, ny

    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: zetaa
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: du1a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: dv1a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: du2a
    real(kind=rp), dimension(-1:nx+2,-1:ny+2), target, intent(in) :: dv2a

    real(kind=rp), dimension(:,:), pointer :: zeta,du1,dv1,du2,dv2

    integer(kind=ip), save :: iter_write_pred_2d=0
    iter_write_pred_2d = iter_write_pred_2d + 1

    zeta => zetaa
    du1 => du1a
    dv1 => dv1a
    du2 => du2a
    dv2 => dv2a

    call write_netcdf(zeta,vname='zeta',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(du1,vname='du1',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(dv1,vname='dv1',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(du2,vname='du2',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)
    call write_netcdf(dv2,vname='dv2',netcdf_file_name='wrt_pred_2d.nc',rank=myrank,iter=iter_write_pred_2d)

  end subroutine h_write_pred_2d

  !--------------------------------------------------------------
  subroutine h_write_pred_out(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_pred_out=0
    iter_write_pred_out = iter_write_pred_out + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_pred_out.nc',rank=myrank,iter=iter_write_pred_out)

  end subroutine h_write_pred_out

  !--------------------------------------------------------------
  subroutine h_write_corr_in(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_in=0
    iter_write_corr_in = iter_write_corr_in + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_in.nc',rank=myrank,iter=iter_write_corr_in)

  end subroutine h_write_corr_in

  !--------------------------------------------------------------
  subroutine h_write_corr_inter(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_inter=0
    iter_write_corr_inter = iter_write_corr_inter + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_inter.nc',rank=myrank,iter=iter_write_corr_inter)

  end subroutine h_write_corr_inter

  !--------------------------------------------------------------
  subroutine h_write_corr_out(nx,ny,nz,ua,va)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: ua
    real(kind=rp), dimension(-1:nx+2,-1:ny+2,1:nz),   target, intent(in) :: va

    real(kind=rp), dimension(:,:,:), pointer :: u,v

    integer(kind=ip), save :: iter_write_corr_out=0
    iter_write_corr_out = iter_write_corr_out + 1

    u => ua
    v => va

    call write_netcdf(u,vname='u',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)
    call write_netcdf(v,vname='v',netcdf_file_name='wrt_corr_out.nc',rank=myrank,iter=iter_write_corr_out)

  end subroutine h_write_corr_out

end module nhmg_debug

