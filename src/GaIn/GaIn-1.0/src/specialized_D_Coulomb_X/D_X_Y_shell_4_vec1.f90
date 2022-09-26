module mod_D_X_Y_shell_optv1_4
  
  implicit none
  
  private
  public :: D_X_Y_shell_optv1_4
  integer, parameter :: NVEC=1
  
contains
  
  

!> Compute D X Y/R coulomb integral for the right hand shell 
!!
recursive function D_X_Y_shell_optv1_4(a12,a3,r3,c,lmax,value)
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12          !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3(NVEC)     !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(NVEC,3)   !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)         !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax         !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value(*)     !< result value
  
  ! return value
  logical :: D_X_Y_shell_optv1_4
  
  ! local variables
  integer      :: i
  real(kind=8) :: d(NVEC)
  real(kind=8) :: p(NVEC)
  real(kind=8) :: pd(NVEC)
  real(kind=8) :: pd2(NVEC)
  real(kind=8) :: erfpd(NVEC)
  real(kind=8) :: exppd2(NVEC)
  real(kind=8) :: x(NVEC)
  real(kind=8) :: y(NVEC)
  real(kind=8) :: z(NVEC)
  
  ! temporary variables
  real(kind=8) :: u2885(NVEC)
  real(kind=8) :: u2886(NVEC)
  real(kind=8) :: u2934(NVEC)
  real(kind=8) :: u2936(NVEC)
  real(kind=8) :: u2937(NVEC)
  real(kind=8) :: u2980(NVEC)
  real(kind=8) :: u3770(NVEC)
  real(kind=8) :: u3772(NVEC)
  real(kind=8) :: u3773(NVEC)
  real(kind=8) :: u3774(NVEC)
  real(kind=8) :: u3775(NVEC)
  real(kind=8) :: u3776(NVEC)
  real(kind=8) :: u3777(NVEC)
  real(kind=8) :: u3778(NVEC)
  real(kind=8) :: u3779(NVEC)
  real(kind=8) :: u3780(NVEC)
  real(kind=8) :: u3781(NVEC)
  real(kind=8) :: u3782(NVEC)
  real(kind=8) :: u3783(NVEC)
  real(kind=8) :: u3784(NVEC)
  real(kind=8) :: u3785(NVEC)
  real(kind=8) :: u3786(NVEC)
  real(kind=8) :: u3788(NVEC)
  real(kind=8) :: u3789(NVEC)
  real(kind=8) :: u3791(NVEC)
  real(kind=8) :: u3792(NVEC)
  real(kind=8) :: u3793(NVEC)
  real(kind=8) :: u3794(NVEC)
  real(kind=8) :: u3795(NVEC)
  real(kind=8) :: u3796(NVEC)
  real(kind=8) :: u3798(NVEC)
  real(kind=8) :: u3799(NVEC)
  real(kind=8) :: u3800(NVEC)
  real(kind=8) :: u3801(NVEC)
  real(kind=8) :: u3805(NVEC)
  real(kind=8) :: u3806(NVEC)
  real(kind=8) :: u3807(NVEC)
  real(kind=8) :: u3808(NVEC)
  real(kind=8) :: u3809(NVEC)
  real(kind=8) :: u3810(NVEC)
  real(kind=8) :: u3811(NVEC)
  real(kind=8) :: u3812(NVEC)
  real(kind=8) :: u3813(NVEC)
  real(kind=8) :: u3814(NVEC)
  real(kind=8) :: u3816(NVEC)
  real(kind=8) :: u3817(NVEC)
  real(kind=8) :: u3818(NVEC)
  real(kind=8) :: u3819(NVEC)
  real(kind=8) :: u382(NVEC)
  real(kind=8) :: u3820(NVEC)
  real(kind=8) :: u3821(NVEC)
  real(kind=8) :: u3822(NVEC)
  real(kind=8) :: u3823(NVEC)
  real(kind=8) :: u3824(NVEC)
  real(kind=8) :: u3825(NVEC)
  real(kind=8) :: u3826(NVEC)
  real(kind=8) :: u3827(NVEC)
  real(kind=8) :: u3828(NVEC)
  real(kind=8) :: u3829(NVEC)
  real(kind=8) :: u383(NVEC)
  real(kind=8) :: u3830(NVEC)
  real(kind=8) :: u3831(NVEC)
  real(kind=8) :: u3832(NVEC)
  real(kind=8) :: u3833(NVEC)
  real(kind=8) :: u3834(NVEC)
  real(kind=8) :: u3835(NVEC)
  real(kind=8) :: u3836(NVEC)
  real(kind=8) :: u3837(NVEC)
  real(kind=8) :: u3838(NVEC)
  real(kind=8) :: u3839(NVEC)
  real(kind=8) :: u384(NVEC)
  real(kind=8) :: u3840(NVEC)
  real(kind=8) :: u3841(NVEC)
  real(kind=8) :: u3842(NVEC)
  real(kind=8) :: u3843(NVEC)
  real(kind=8) :: u3844(NVEC)
  real(kind=8) :: u3845(NVEC)
  real(kind=8) :: u3846(NVEC)
  real(kind=8) :: u3847(NVEC)
  real(kind=8) :: u3848(NVEC)
  real(kind=8) :: u3849(NVEC)
  real(kind=8) :: u385(NVEC)
  real(kind=8) :: u3850(NVEC)
  real(kind=8) :: u3851(NVEC)
  real(kind=8) :: u3852(NVEC)
  real(kind=8) :: u3853(NVEC)
  real(kind=8) :: u3854(NVEC)
  real(kind=8) :: u3855(NVEC)
  real(kind=8) :: u3856(NVEC)
  real(kind=8) :: u3857(NVEC)
  real(kind=8) :: u3858(NVEC)
  real(kind=8) :: u3859(NVEC)
  real(kind=8) :: u386(NVEC)
  real(kind=8) :: u3860(NVEC)
  real(kind=8) :: u3861(NVEC)
  real(kind=8) :: u3862(NVEC)
  real(kind=8) :: u3863(NVEC)
  real(kind=8) :: u3864(NVEC)
  real(kind=8) :: u3865(NVEC)
  real(kind=8) :: u3866(NVEC)
  real(kind=8) :: u3867(NVEC)
  real(kind=8) :: u3868(NVEC)
  real(kind=8) :: u3869(NVEC)
  real(kind=8) :: u387(NVEC)
  real(kind=8) :: u3870(NVEC)
  real(kind=8) :: u3871(NVEC)
  real(kind=8) :: u3872(NVEC)
  real(kind=8) :: u3873(NVEC)
  real(kind=8) :: u3874(NVEC)
  real(kind=8) :: u3875(NVEC)
  real(kind=8) :: u3876(NVEC)
  real(kind=8) :: u3877(NVEC)
  real(kind=8) :: u3878(NVEC)
  real(kind=8) :: u3879(NVEC)
  real(kind=8) :: u388(NVEC)
  real(kind=8) :: u3880(NVEC)
  real(kind=8) :: u3881(NVEC)
  real(kind=8) :: u3882(NVEC)
  real(kind=8) :: u3883(NVEC)
  real(kind=8) :: u3884(NVEC)
  real(kind=8) :: u3885(NVEC)
  real(kind=8) :: u3886(NVEC)
  real(kind=8) :: u3887(NVEC)
  real(kind=8) :: u3888(NVEC)
  real(kind=8) :: u3889(NVEC)
  real(kind=8) :: u389(NVEC)
  real(kind=8) :: u3890(NVEC)
  real(kind=8) :: u3891(NVEC)
  real(kind=8) :: u3892(NVEC)
  real(kind=8) :: u3893(NVEC)
  real(kind=8) :: u3894(NVEC)
  real(kind=8) :: u3895(NVEC)
  real(kind=8) :: u3896(NVEC)
  real(kind=8) :: u3897(NVEC)
  real(kind=8) :: u3898(NVEC)
  real(kind=8) :: u3899(NVEC)
  real(kind=8) :: u39(NVEC)
  real(kind=8) :: u390(NVEC)
  real(kind=8) :: u3900(NVEC)
  real(kind=8) :: u3901(NVEC)
  real(kind=8) :: u3902(NVEC)
  real(kind=8) :: u3903(NVEC)
  real(kind=8) :: u3904(NVEC)
  real(kind=8) :: u3905(NVEC)
  real(kind=8) :: u3906(NVEC)
  real(kind=8) :: u3907(NVEC)
  real(kind=8) :: u3908(NVEC)
  real(kind=8) :: u3909(NVEC)
  real(kind=8) :: u391(NVEC)
  real(kind=8) :: u3910(NVEC)
  real(kind=8) :: u3911(NVEC)
  real(kind=8) :: u3912(NVEC)
  real(kind=8) :: u3913(NVEC)
  real(kind=8) :: u3914(NVEC)
  real(kind=8) :: u3915(NVEC)
  real(kind=8) :: u3916(NVEC)
  real(kind=8) :: u3917(NVEC)
  real(kind=8) :: u3918(NVEC)
  real(kind=8) :: u3919(NVEC)
  real(kind=8) :: u392(NVEC)
  real(kind=8) :: u3920(NVEC)
  real(kind=8) :: u3921(NVEC)
  real(kind=8) :: u3922(NVEC)
  real(kind=8) :: u3923(NVEC)
  real(kind=8) :: u3924(NVEC)
  real(kind=8) :: u3925(NVEC)
  real(kind=8) :: u3926(NVEC)
  real(kind=8) :: u3927(NVEC)
  real(kind=8) :: u3928(NVEC)
  real(kind=8) :: u3929(NVEC)
  real(kind=8) :: u393(NVEC)
  real(kind=8) :: u3930(NVEC)
  real(kind=8) :: u3931(NVEC)
  real(kind=8) :: u3932(NVEC)
  real(kind=8) :: u3933(NVEC)
  real(kind=8) :: u3934(NVEC)
  real(kind=8) :: u3935(NVEC)
  real(kind=8) :: u3936(NVEC)
  real(kind=8) :: u3937(NVEC)
  real(kind=8) :: u3938(NVEC)
  real(kind=8) :: u3939(NVEC)
  real(kind=8) :: u394(NVEC)
  real(kind=8) :: u3940(NVEC)
  real(kind=8) :: u3941(NVEC)
  real(kind=8) :: u3942(NVEC)
  real(kind=8) :: u3943(NVEC)
  real(kind=8) :: u3944(NVEC)
  real(kind=8) :: u3945(NVEC)
  real(kind=8) :: u3946(NVEC)
  real(kind=8) :: u3947(NVEC)
  real(kind=8) :: u3948(NVEC)
  real(kind=8) :: u3949(NVEC)
  real(kind=8) :: u395(NVEC)
  real(kind=8) :: u3950(NVEC)
  real(kind=8) :: u3952(NVEC)
  real(kind=8) :: u3953(NVEC)
  real(kind=8) :: u3954(NVEC)
  real(kind=8) :: u3955(NVEC)
  real(kind=8) :: u3956(NVEC)
  real(kind=8) :: u3957(NVEC)
  real(kind=8) :: u3958(NVEC)
  real(kind=8) :: u3959(NVEC)
  real(kind=8) :: u396(NVEC)
  real(kind=8) :: u3960(NVEC)
  real(kind=8) :: u3961(NVEC)
  real(kind=8) :: u3962(NVEC)
  real(kind=8) :: u3963(NVEC)
  real(kind=8) :: u3964(NVEC)
  real(kind=8) :: u3965(NVEC)
  real(kind=8) :: u3966(NVEC)
  real(kind=8) :: u3967(NVEC)
  real(kind=8) :: u3968(NVEC)
  real(kind=8) :: u3969(NVEC)
  real(kind=8) :: u397(NVEC)
  real(kind=8) :: u3970(NVEC)
  real(kind=8) :: u3971(NVEC)
  real(kind=8) :: u3972(NVEC)
  real(kind=8) :: u3973(NVEC)
  real(kind=8) :: u3974(NVEC)
  real(kind=8) :: u3975(NVEC)
  real(kind=8) :: u3976(NVEC)
  real(kind=8) :: u3977(NVEC)
  real(kind=8) :: u3978(NVEC)
  real(kind=8) :: u3979(NVEC)
  real(kind=8) :: u398(NVEC)
  real(kind=8) :: u3980(NVEC)
  real(kind=8) :: u3981(NVEC)
  real(kind=8) :: u3982(NVEC)
  real(kind=8) :: u3983(NVEC)
  real(kind=8) :: u3984(NVEC)
  real(kind=8) :: u3985(NVEC)
  real(kind=8) :: u3986(NVEC)
  real(kind=8) :: u3987(NVEC)
  real(kind=8) :: u3988(NVEC)
  real(kind=8) :: u3989(NVEC)
  real(kind=8) :: u399(NVEC)
  real(kind=8) :: u3990(NVEC)
  real(kind=8) :: u3991(NVEC)
  real(kind=8) :: u3992(NVEC)
  real(kind=8) :: u3993(NVEC)
  real(kind=8) :: u3994(NVEC)
  real(kind=8) :: u3995(NVEC)
  real(kind=8) :: u3996(NVEC)
  real(kind=8) :: u3997(NVEC)
  real(kind=8) :: u3998(NVEC)
  real(kind=8) :: u3999(NVEC)
  real(kind=8) :: u4(NVEC)
  real(kind=8) :: u40(NVEC)
  real(kind=8) :: u400(NVEC)
  real(kind=8) :: u4000(NVEC)
  real(kind=8) :: u4001(NVEC)
  real(kind=8) :: u4002(NVEC)
  real(kind=8) :: u4003(NVEC)
  real(kind=8) :: u4004(NVEC)
  real(kind=8) :: u4005(NVEC)
  real(kind=8) :: u4006(NVEC)
  real(kind=8) :: u4007(NVEC)
  real(kind=8) :: u4008(NVEC)
  real(kind=8) :: u4009(NVEC)
  real(kind=8) :: u401(NVEC)
  real(kind=8) :: u4010(NVEC)
  real(kind=8) :: u4011(NVEC)
  real(kind=8) :: u4012(NVEC)
  real(kind=8) :: u4013(NVEC)
  real(kind=8) :: u4014(NVEC)
  real(kind=8) :: u4015(NVEC)
  real(kind=8) :: u4016(NVEC)
  real(kind=8) :: u4017(NVEC)
  real(kind=8) :: u4018(NVEC)
  real(kind=8) :: u4019(NVEC)
  real(kind=8) :: u402(NVEC)
  real(kind=8) :: u4020(NVEC)
  real(kind=8) :: u4021(NVEC)
  real(kind=8) :: u4022(NVEC)
  real(kind=8) :: u4023(NVEC)
  real(kind=8) :: u4024(NVEC)
  real(kind=8) :: u4025(NVEC)
  real(kind=8) :: u4026(NVEC)
  real(kind=8) :: u4027(NVEC)
  real(kind=8) :: u4028(NVEC)
  real(kind=8) :: u4029(NVEC)
  real(kind=8) :: u403(NVEC)
  real(kind=8) :: u4030(NVEC)
  real(kind=8) :: u4031(NVEC)
  real(kind=8) :: u4032(NVEC)
  real(kind=8) :: u4033(NVEC)
  real(kind=8) :: u4034(NVEC)
  real(kind=8) :: u4035(NVEC)
  real(kind=8) :: u4036(NVEC)
  real(kind=8) :: u4037(NVEC)
  real(kind=8) :: u4038(NVEC)
  real(kind=8) :: u4039(NVEC)
  real(kind=8) :: u404(NVEC)
  real(kind=8) :: u4040(NVEC)
  real(kind=8) :: u4041(NVEC)
  real(kind=8) :: u4042(NVEC)
  real(kind=8) :: u4043(NVEC)
  real(kind=8) :: u4044(NVEC)
  real(kind=8) :: u4045(NVEC)
  real(kind=8) :: u4046(NVEC)
  real(kind=8) :: u4047(NVEC)
  real(kind=8) :: u4048(NVEC)
  real(kind=8) :: u4049(NVEC)
  real(kind=8) :: u405(NVEC)
  real(kind=8) :: u4050(NVEC)
  real(kind=8) :: u4051(NVEC)
  real(kind=8) :: u4052(NVEC)
  real(kind=8) :: u4053(NVEC)
  real(kind=8) :: u4054(NVEC)
  real(kind=8) :: u4055(NVEC)
  real(kind=8) :: u4056(NVEC)
  real(kind=8) :: u4057(NVEC)
  real(kind=8) :: u4058(NVEC)
  real(kind=8) :: u4059(NVEC)
  real(kind=8) :: u406(NVEC)
  real(kind=8) :: u4060(NVEC)
  real(kind=8) :: u4061(NVEC)
  real(kind=8) :: u4062(NVEC)
  real(kind=8) :: u4063(NVEC)
  real(kind=8) :: u4064(NVEC)
  real(kind=8) :: u4065(NVEC)
  real(kind=8) :: u4066(NVEC)
  real(kind=8) :: u4067(NVEC)
  real(kind=8) :: u4068(NVEC)
  real(kind=8) :: u4069(NVEC)
  real(kind=8) :: u407(NVEC)
  real(kind=8) :: u4070(NVEC)
  real(kind=8) :: u4071(NVEC)
  real(kind=8) :: u4072(NVEC)
  real(kind=8) :: u4073(NVEC)
  real(kind=8) :: u4074(NVEC)
  real(kind=8) :: u4075(NVEC)
  real(kind=8) :: u4076(NVEC)
  real(kind=8) :: u4077(NVEC)
  real(kind=8) :: u4078(NVEC)
  real(kind=8) :: u4079(NVEC)
  real(kind=8) :: u408(NVEC)
  real(kind=8) :: u4080(NVEC)
  real(kind=8) :: u4081(NVEC)
  real(kind=8) :: u4082(NVEC)
  real(kind=8) :: u4083(NVEC)
  real(kind=8) :: u4084(NVEC)
  real(kind=8) :: u4085(NVEC)
  real(kind=8) :: u4086(NVEC)
  real(kind=8) :: u4087(NVEC)
  real(kind=8) :: u4088(NVEC)
  real(kind=8) :: u4089(NVEC)
  real(kind=8) :: u409(NVEC)
  real(kind=8) :: u4090(NVEC)
  real(kind=8) :: u4091(NVEC)
  real(kind=8) :: u4092(NVEC)
  real(kind=8) :: u4093(NVEC)
  real(kind=8) :: u4094(NVEC)
  real(kind=8) :: u4095(NVEC)
  real(kind=8) :: u4096(NVEC)
  real(kind=8) :: u4097(NVEC)
  real(kind=8) :: u4098(NVEC)
  real(kind=8) :: u4099(NVEC)
  real(kind=8) :: u41(NVEC)
  real(kind=8) :: u410(NVEC)
  real(kind=8) :: u4100(NVEC)
  real(kind=8) :: u4101(NVEC)
  real(kind=8) :: u4102(NVEC)
  real(kind=8) :: u4103(NVEC)
  real(kind=8) :: u4104(NVEC)
  real(kind=8) :: u4105(NVEC)
  real(kind=8) :: u4106(NVEC)
  real(kind=8) :: u4107(NVEC)
  real(kind=8) :: u4108(NVEC)
  real(kind=8) :: u4109(NVEC)
  real(kind=8) :: u411(NVEC)
  real(kind=8) :: u4110(NVEC)
  real(kind=8) :: u4111(NVEC)
  real(kind=8) :: u4112(NVEC)
  real(kind=8) :: u4113(NVEC)
  real(kind=8) :: u4114(NVEC)
  real(kind=8) :: u4115(NVEC)
  real(kind=8) :: u4116(NVEC)
  real(kind=8) :: u4117(NVEC)
  real(kind=8) :: u4118(NVEC)
  real(kind=8) :: u4119(NVEC)
  real(kind=8) :: u412(NVEC)
  real(kind=8) :: u4120(NVEC)
  real(kind=8) :: u4121(NVEC)
  real(kind=8) :: u4122(NVEC)
  real(kind=8) :: u4123(NVEC)
  real(kind=8) :: u4124(NVEC)
  real(kind=8) :: u4125(NVEC)
  real(kind=8) :: u4126(NVEC)
  real(kind=8) :: u4127(NVEC)
  real(kind=8) :: u4128(NVEC)
  real(kind=8) :: u4129(NVEC)
  real(kind=8) :: u413(NVEC)
  real(kind=8) :: u4130(NVEC)
  real(kind=8) :: u4131(NVEC)
  real(kind=8) :: u4132(NVEC)
  real(kind=8) :: u4133(NVEC)
  real(kind=8) :: u4134(NVEC)
  real(kind=8) :: u4135(NVEC)
  real(kind=8) :: u4136(NVEC)
  real(kind=8) :: u4137(NVEC)
  real(kind=8) :: u4138(NVEC)
  real(kind=8) :: u4139(NVEC)
  real(kind=8) :: u414(NVEC)
  real(kind=8) :: u4140(NVEC)
  real(kind=8) :: u4141(NVEC)
  real(kind=8) :: u4142(NVEC)
  real(kind=8) :: u4143(NVEC)
  real(kind=8) :: u4144(NVEC)
  real(kind=8) :: u4145(NVEC)
  real(kind=8) :: u4146(NVEC)
  real(kind=8) :: u4147(NVEC)
  real(kind=8) :: u4148(NVEC)
  real(kind=8) :: u4149(NVEC)
  real(kind=8) :: u415(NVEC)
  real(kind=8) :: u4150(NVEC)
  real(kind=8) :: u4151(NVEC)
  real(kind=8) :: u4152(NVEC)
  real(kind=8) :: u4153(NVEC)
  real(kind=8) :: u4154(NVEC)
  real(kind=8) :: u4155(NVEC)
  real(kind=8) :: u4156(NVEC)
  real(kind=8) :: u4157(NVEC)
  real(kind=8) :: u4158(NVEC)
  real(kind=8) :: u4159(NVEC)
  real(kind=8) :: u416(NVEC)
  real(kind=8) :: u4160(NVEC)
  real(kind=8) :: u4161(NVEC)
  real(kind=8) :: u4162(NVEC)
  real(kind=8) :: u4163(NVEC)
  real(kind=8) :: u4164(NVEC)
  real(kind=8) :: u4165(NVEC)
  real(kind=8) :: u4166(NVEC)
  real(kind=8) :: u4167(NVEC)
  real(kind=8) :: u4168(NVEC)
  real(kind=8) :: u4169(NVEC)
  real(kind=8) :: u417(NVEC)
  real(kind=8) :: u4170(NVEC)
  real(kind=8) :: u4171(NVEC)
  real(kind=8) :: u4172(NVEC)
  real(kind=8) :: u4173(NVEC)
  real(kind=8) :: u4174(NVEC)
  real(kind=8) :: u4175(NVEC)
  real(kind=8) :: u4176(NVEC)
  real(kind=8) :: u4177(NVEC)
  real(kind=8) :: u4178(NVEC)
  real(kind=8) :: u4179(NVEC)
  real(kind=8) :: u418(NVEC)
  real(kind=8) :: u4180(NVEC)
  real(kind=8) :: u4181(NVEC)
  real(kind=8) :: u4182(NVEC)
  real(kind=8) :: u4183(NVEC)
  real(kind=8) :: u4184(NVEC)
  real(kind=8) :: u4185(NVEC)
  real(kind=8) :: u4186(NVEC)
  real(kind=8) :: u4187(NVEC)
  real(kind=8) :: u4188(NVEC)
  real(kind=8) :: u4189(NVEC)
  real(kind=8) :: u419(NVEC)
  real(kind=8) :: u4190(NVEC)
  real(kind=8) :: u4191(NVEC)
  real(kind=8) :: u4192(NVEC)
  real(kind=8) :: u4193(NVEC)
  real(kind=8) :: u4194(NVEC)
  real(kind=8) :: u4195(NVEC)
  real(kind=8) :: u4196(NVEC)
  real(kind=8) :: u4197(NVEC)
  real(kind=8) :: u4198(NVEC)
  real(kind=8) :: u4199(NVEC)
  real(kind=8) :: u42(NVEC)
  real(kind=8) :: u420(NVEC)
  real(kind=8) :: u4200(NVEC)
  real(kind=8) :: u4201(NVEC)
  real(kind=8) :: u4202(NVEC)
  real(kind=8) :: u4203(NVEC)
  real(kind=8) :: u4204(NVEC)
  real(kind=8) :: u4205(NVEC)
  real(kind=8) :: u4206(NVEC)
  real(kind=8) :: u4207(NVEC)
  real(kind=8) :: u4208(NVEC)
  real(kind=8) :: u4209(NVEC)
  real(kind=8) :: u421(NVEC)
  real(kind=8) :: u4210(NVEC)
  real(kind=8) :: u4211(NVEC)
  real(kind=8) :: u4212(NVEC)
  real(kind=8) :: u422(NVEC)
  real(kind=8) :: u423(NVEC)
  real(kind=8) :: u424(NVEC)
  real(kind=8) :: u425(NVEC)
  real(kind=8) :: u426(NVEC)
  real(kind=8) :: u427(NVEC)
  real(kind=8) :: u428(NVEC)
  real(kind=8) :: u429(NVEC)
  real(kind=8) :: u43(NVEC)
  real(kind=8) :: u430(NVEC)
  real(kind=8) :: u431(NVEC)
  real(kind=8) :: u432(NVEC)
  real(kind=8) :: u433(NVEC)
  real(kind=8) :: u434(NVEC)
  real(kind=8) :: u435(NVEC)
  real(kind=8) :: u436(NVEC)
  real(kind=8) :: u437(NVEC)
  real(kind=8) :: u438(NVEC)
  real(kind=8) :: u439(NVEC)
  real(kind=8) :: u44(NVEC)
  real(kind=8) :: u440(NVEC)
  real(kind=8) :: u441(NVEC)
  real(kind=8) :: u442(NVEC)
  real(kind=8) :: u443(NVEC)
  real(kind=8) :: u444(NVEC)
  real(kind=8) :: u445(NVEC)
  real(kind=8) :: u446(NVEC)
  real(kind=8) :: u447(NVEC)
  real(kind=8) :: u448(NVEC)
  real(kind=8) :: u449(NVEC)
  real(kind=8) :: u45(NVEC)
  real(kind=8) :: u450(NVEC)
  real(kind=8) :: u451(NVEC)
  real(kind=8) :: u452(NVEC)
  real(kind=8) :: u453(NVEC)
  real(kind=8) :: u454(NVEC)
  real(kind=8) :: u455(NVEC)
  real(kind=8) :: u456(NVEC)
  real(kind=8) :: u457(NVEC)
  real(kind=8) :: u458(NVEC)
  real(kind=8) :: u459(NVEC)
  real(kind=8) :: u46(NVEC)
  real(kind=8) :: u460(NVEC)
  real(kind=8) :: u461(NVEC)
  real(kind=8) :: u462(NVEC)
  real(kind=8) :: u463(NVEC)
  real(kind=8) :: u464(NVEC)
  real(kind=8) :: u465(NVEC)
  real(kind=8) :: u466(NVEC)
  real(kind=8) :: u467(NVEC)
  real(kind=8) :: u468(NVEC)
  real(kind=8) :: u469(NVEC)
  real(kind=8) :: u47(NVEC)
  real(kind=8) :: u470(NVEC)
  real(kind=8) :: u471(NVEC)
  real(kind=8) :: u472(NVEC)
  real(kind=8) :: u473(NVEC)
  real(kind=8) :: u474(NVEC)
  real(kind=8) :: u475(NVEC)
  real(kind=8) :: u476(NVEC)
  real(kind=8) :: u477(NVEC)
  real(kind=8) :: u478(NVEC)
  real(kind=8) :: u479(NVEC)
  real(kind=8) :: u48(NVEC)
  real(kind=8) :: u480(NVEC)
  real(kind=8) :: u481(NVEC)
  real(kind=8) :: u482(NVEC)
  real(kind=8) :: u483(NVEC)
  real(kind=8) :: u484(NVEC)
  real(kind=8) :: u485(NVEC)
  real(kind=8) :: u486(NVEC)
  real(kind=8) :: u487(NVEC)
  real(kind=8) :: u488(NVEC)
  real(kind=8) :: u489(NVEC)
  real(kind=8) :: u49(NVEC)
  real(kind=8) :: u490(NVEC)
  real(kind=8) :: u491(NVEC)
  real(kind=8) :: u492(NVEC)
  real(kind=8) :: u493(NVEC)
  real(kind=8) :: u494(NVEC)
  real(kind=8) :: u495(NVEC)
  real(kind=8) :: u496(NVEC)
  real(kind=8) :: u497(NVEC)
  real(kind=8) :: u498(NVEC)
  real(kind=8) :: u499(NVEC)
  real(kind=8) :: u5(NVEC)
  real(kind=8) :: u50(NVEC)
  real(kind=8) :: u500(NVEC)
  real(kind=8) :: u501(NVEC)
  real(kind=8) :: u502(NVEC)
  real(kind=8) :: u503(NVEC)
  real(kind=8) :: u504(NVEC)
  real(kind=8) :: u505(NVEC)
  real(kind=8) :: u506(NVEC)
  real(kind=8) :: u507(NVEC)
  real(kind=8) :: u508(NVEC)
  real(kind=8) :: u509(NVEC)
  real(kind=8) :: u51(NVEC)
  real(kind=8) :: u510(NVEC)
  real(kind=8) :: u511(NVEC)
  real(kind=8) :: u512(NVEC)
  real(kind=8) :: u513(NVEC)
  real(kind=8) :: u514(NVEC)
  real(kind=8) :: u515(NVEC)
  real(kind=8) :: u516(NVEC)
  real(kind=8) :: u517(NVEC)
  real(kind=8) :: u518(NVEC)
  real(kind=8) :: u519(NVEC)
  real(kind=8) :: u52(NVEC)
  real(kind=8) :: u520(NVEC)
  real(kind=8) :: u521(NVEC)
  real(kind=8) :: u522(NVEC)
  real(kind=8) :: u523(NVEC)
  real(kind=8) :: u524(NVEC)
  real(kind=8) :: u525(NVEC)
  real(kind=8) :: u526(NVEC)
  real(kind=8) :: u527(NVEC)
  real(kind=8) :: u528(NVEC)
  real(kind=8) :: u529(NVEC)
  real(kind=8) :: u53(NVEC)
  real(kind=8) :: u530(NVEC)
  real(kind=8) :: u531(NVEC)
  real(kind=8) :: u532(NVEC)
  real(kind=8) :: u533(NVEC)
  real(kind=8) :: u534(NVEC)
  real(kind=8) :: u535(NVEC)
  real(kind=8) :: u536(NVEC)
  real(kind=8) :: u537(NVEC)
  real(kind=8) :: u538(NVEC)
  real(kind=8) :: u539(NVEC)
  real(kind=8) :: u54(NVEC)
  real(kind=8) :: u540(NVEC)
  real(kind=8) :: u541(NVEC)
  real(kind=8) :: u542(NVEC)
  real(kind=8) :: u543(NVEC)
  real(kind=8) :: u544(NVEC)
  real(kind=8) :: u545(NVEC)
  real(kind=8) :: u546(NVEC)
  real(kind=8) :: u547(NVEC)
  real(kind=8) :: u549(NVEC)
  real(kind=8) :: u55(NVEC)
  real(kind=8) :: u550(NVEC)
  real(kind=8) :: u551(NVEC)
  real(kind=8) :: u552(NVEC)
  real(kind=8) :: u553(NVEC)
  real(kind=8) :: u554(NVEC)
  real(kind=8) :: u555(NVEC)
  real(kind=8) :: u556(NVEC)
  real(kind=8) :: u557(NVEC)
  real(kind=8) :: u558(NVEC)
  real(kind=8) :: u559(NVEC)
  real(kind=8) :: u56(NVEC)
  real(kind=8) :: u560(NVEC)
  real(kind=8) :: u562(NVEC)
  real(kind=8) :: u563(NVEC)
  real(kind=8) :: u564(NVEC)
  real(kind=8) :: u565(NVEC)
  real(kind=8) :: u566(NVEC)
  real(kind=8) :: u567(NVEC)
  real(kind=8) :: u568(NVEC)
  real(kind=8) :: u569(NVEC)
  real(kind=8) :: u57(NVEC)
  real(kind=8) :: u570(NVEC)
  real(kind=8) :: u571(NVEC)
  real(kind=8) :: u572(NVEC)
  real(kind=8) :: u573(NVEC)
  real(kind=8) :: u574(NVEC)
  real(kind=8) :: u575(NVEC)
  real(kind=8) :: u576(NVEC)
  real(kind=8) :: u577(NVEC)
  real(kind=8) :: u578(NVEC)
  real(kind=8) :: u579(NVEC)
  real(kind=8) :: u58(NVEC)
  real(kind=8) :: u580(NVEC)
  real(kind=8) :: u581(NVEC)
  real(kind=8) :: u582(NVEC)
  real(kind=8) :: u583(NVEC)
  real(kind=8) :: u584(NVEC)
  real(kind=8) :: u585(NVEC)
  real(kind=8) :: u586(NVEC)
  real(kind=8) :: u587(NVEC)
  real(kind=8) :: u588(NVEC)
  real(kind=8) :: u589(NVEC)
  real(kind=8) :: u59(NVEC)
  real(kind=8) :: u590(NVEC)
  real(kind=8) :: u591(NVEC)
  real(kind=8) :: u592(NVEC)
  real(kind=8) :: u593(NVEC)
  real(kind=8) :: u594(NVEC)
  real(kind=8) :: u595(NVEC)
  real(kind=8) :: u596(NVEC)
  real(kind=8) :: u597(NVEC)
  real(kind=8) :: u598(NVEC)
  real(kind=8) :: u599(NVEC)
  real(kind=8) :: u6(NVEC)
  real(kind=8) :: u60(NVEC)
  real(kind=8) :: u600(NVEC)
  real(kind=8) :: u601(NVEC)
  real(kind=8) :: u602(NVEC)
  real(kind=8) :: u604(NVEC)
  real(kind=8) :: u605(NVEC)
  real(kind=8) :: u606(NVEC)
  real(kind=8) :: u607(NVEC)
  real(kind=8) :: u608(NVEC)
  real(kind=8) :: u609(NVEC)
  real(kind=8) :: u61(NVEC)
  real(kind=8) :: u611(NVEC)
  real(kind=8) :: u612(NVEC)
  real(kind=8) :: u613(NVEC)
  real(kind=8) :: u614(NVEC)
  real(kind=8) :: u615(NVEC)
  real(kind=8) :: u616(NVEC)
  real(kind=8) :: u617(NVEC)
  real(kind=8) :: u618(NVEC)
  real(kind=8) :: u619(NVEC)
  real(kind=8) :: u62(NVEC)
  real(kind=8) :: u620(NVEC)
  real(kind=8) :: u621(NVEC)
  real(kind=8) :: u622(NVEC)
  real(kind=8) :: u623(NVEC)
  real(kind=8) :: u624(NVEC)
  real(kind=8) :: u625(NVEC)
  real(kind=8) :: u627(NVEC)
  real(kind=8) :: u628(NVEC)
  real(kind=8) :: u629(NVEC)
  real(kind=8) :: u63(NVEC)
  real(kind=8) :: u630(NVEC)
  real(kind=8) :: u631(NVEC)
  real(kind=8) :: u632(NVEC)
  real(kind=8) :: u633(NVEC)
  real(kind=8) :: u634(NVEC)
  real(kind=8) :: u635(NVEC)
  real(kind=8) :: u636(NVEC)
  real(kind=8) :: u637(NVEC)
  real(kind=8) :: u638(NVEC)
  real(kind=8) :: u639(NVEC)
  real(kind=8) :: u64(NVEC)
  real(kind=8) :: u640(NVEC)
  real(kind=8) :: u641(NVEC)
  real(kind=8) :: u642(NVEC)
  real(kind=8) :: u643(NVEC)
  real(kind=8) :: u644(NVEC)
  real(kind=8) :: u645(NVEC)
  real(kind=8) :: u646(NVEC)
  real(kind=8) :: u647(NVEC)
  real(kind=8) :: u648(NVEC)
  real(kind=8) :: u649(NVEC)
  real(kind=8) :: u65(NVEC)
  real(kind=8) :: u650(NVEC)
  real(kind=8) :: u651(NVEC)
  real(kind=8) :: u652(NVEC)
  real(kind=8) :: u653(NVEC)
  real(kind=8) :: u654(NVEC)
  real(kind=8) :: u655(NVEC)
  real(kind=8) :: u656(NVEC)
  real(kind=8) :: u657(NVEC)
  real(kind=8) :: u658(NVEC)
  real(kind=8) :: u659(NVEC)
  real(kind=8) :: u66(NVEC)
  real(kind=8) :: u660(NVEC)
  real(kind=8) :: u661(NVEC)
  real(kind=8) :: u662(NVEC)
  real(kind=8) :: u663(NVEC)
  real(kind=8) :: u664(NVEC)
  real(kind=8) :: u665(NVEC)
  real(kind=8) :: u666(NVEC)
  real(kind=8) :: u667(NVEC)
  real(kind=8) :: u668(NVEC)
  real(kind=8) :: u669(NVEC)
  real(kind=8) :: u67(NVEC)
  real(kind=8) :: u670(NVEC)
  real(kind=8) :: u671(NVEC)
  real(kind=8) :: u672(NVEC)
  real(kind=8) :: u673(NVEC)
  real(kind=8) :: u674(NVEC)
  real(kind=8) :: u675(NVEC)
  real(kind=8) :: u676(NVEC)
  real(kind=8) :: u677(NVEC)
  real(kind=8) :: u678(NVEC)
  real(kind=8) :: u679(NVEC)
  real(kind=8) :: u68(NVEC)
  real(kind=8) :: u680(NVEC)
  real(kind=8) :: u681(NVEC)
  real(kind=8) :: u682(NVEC)
  real(kind=8) :: u683(NVEC)
  real(kind=8) :: u684(NVEC)
  real(kind=8) :: u685(NVEC)
  real(kind=8) :: u686(NVEC)
  real(kind=8) :: u687(NVEC)
  real(kind=8) :: u688(NVEC)
  real(kind=8) :: u689(NVEC)
  real(kind=8) :: u69(NVEC)
  real(kind=8) :: u690(NVEC)
  real(kind=8) :: u691(NVEC)
  real(kind=8) :: u692(NVEC)
  real(kind=8) :: u693(NVEC)
  real(kind=8) :: u694(NVEC)
  real(kind=8) :: u695(NVEC)
  real(kind=8) :: u696(NVEC)
  real(kind=8) :: u697(NVEC)
  real(kind=8) :: u698(NVEC)
  real(kind=8) :: u699(NVEC)
  real(kind=8) :: u7(NVEC)
  real(kind=8) :: u70(NVEC)
  real(kind=8) :: u700(NVEC)
  real(kind=8) :: u701(NVEC)
  real(kind=8) :: u702(NVEC)
  real(kind=8) :: u703(NVEC)
  real(kind=8) :: u704(NVEC)
  real(kind=8) :: u705(NVEC)
  real(kind=8) :: u706(NVEC)
  real(kind=8) :: u707(NVEC)
  real(kind=8) :: u708(NVEC)
  real(kind=8) :: u709(NVEC)
  real(kind=8) :: u71(NVEC)
  real(kind=8) :: u710(NVEC)
  real(kind=8) :: u711(NVEC)
  real(kind=8) :: u712(NVEC)
  real(kind=8) :: u713(NVEC)
  real(kind=8) :: u714(NVEC)
  real(kind=8) :: u715(NVEC)
  real(kind=8) :: u716(NVEC)
  real(kind=8) :: u717(NVEC)
  real(kind=8) :: u718(NVEC)
  real(kind=8) :: u719(NVEC)
  real(kind=8) :: u72(NVEC)
  real(kind=8) :: u720(NVEC)
  real(kind=8) :: u721(NVEC)
  real(kind=8) :: u722(NVEC)
  real(kind=8) :: u723(NVEC)
  real(kind=8) :: u724(NVEC)
  real(kind=8) :: u725(NVEC)
  real(kind=8) :: u726(NVEC)
  real(kind=8) :: u727(NVEC)
  real(kind=8) :: u728(NVEC)
  real(kind=8) :: u729(NVEC)
  real(kind=8) :: u73(NVEC)
  real(kind=8) :: u730(NVEC)
  real(kind=8) :: u731(NVEC)
  real(kind=8) :: u732(NVEC)
  real(kind=8) :: u733(NVEC)
  real(kind=8) :: u734(NVEC)
  real(kind=8) :: u735(NVEC)
  real(kind=8) :: u736(NVEC)
  real(kind=8) :: u737(NVEC)
  real(kind=8) :: u738(NVEC)
  real(kind=8) :: u739(NVEC)
  real(kind=8) :: u74(NVEC)
  real(kind=8) :: u740(NVEC)
  real(kind=8) :: u741(NVEC)
  real(kind=8) :: u742(NVEC)
  real(kind=8) :: u743(NVEC)
  real(kind=8) :: u744(NVEC)
  real(kind=8) :: u745(NVEC)
  real(kind=8) :: u746(NVEC)
  real(kind=8) :: u747(NVEC)
  real(kind=8) :: u748(NVEC)
  real(kind=8) :: u749(NVEC)
  real(kind=8) :: u75(NVEC)
  real(kind=8) :: u750(NVEC)
  real(kind=8) :: u751(NVEC)
  real(kind=8) :: u752(NVEC)
  real(kind=8) :: u753(NVEC)
  real(kind=8) :: u754(NVEC)
  real(kind=8) :: u755(NVEC)
  real(kind=8) :: u756(NVEC)
  real(kind=8) :: u757(NVEC)
  real(kind=8) :: u758(NVEC)
  real(kind=8) :: u759(NVEC)
  real(kind=8) :: u76(NVEC)
  real(kind=8) :: u760(NVEC)
  real(kind=8) :: u761(NVEC)
  real(kind=8) :: u762(NVEC)
  real(kind=8) :: u763(NVEC)
  real(kind=8) :: u764(NVEC)
  real(kind=8) :: u765(NVEC)
  real(kind=8) :: u766(NVEC)
  real(kind=8) :: u767(NVEC)
  real(kind=8) :: u768(NVEC)
  real(kind=8) :: u769(NVEC)
  real(kind=8) :: u77(NVEC)
  real(kind=8) :: u770(NVEC)
  real(kind=8) :: u771(NVEC)
  real(kind=8) :: u772(NVEC)
  real(kind=8) :: u773(NVEC)
  real(kind=8) :: u774(NVEC)
  real(kind=8) :: u775(NVEC)
  real(kind=8) :: u776(NVEC)
  real(kind=8) :: u777(NVEC)
  real(kind=8) :: u778(NVEC)
  real(kind=8) :: u779(NVEC)
  real(kind=8) :: u78(NVEC)
  real(kind=8) :: u780(NVEC)
  real(kind=8) :: u781(NVEC)
  real(kind=8) :: u782(NVEC)
  real(kind=8) :: u783(NVEC)
  real(kind=8) :: u784(NVEC)
  real(kind=8) :: u785(NVEC)
  real(kind=8) :: u786(NVEC)
  real(kind=8) :: u787(NVEC)
  real(kind=8) :: u788(NVEC)
  real(kind=8) :: u789(NVEC)
  real(kind=8) :: u79(NVEC)
  real(kind=8) :: u790(NVEC)
  real(kind=8) :: u791(NVEC)
  real(kind=8) :: u792(NVEC)
  real(kind=8) :: u793(NVEC)
  real(kind=8) :: u794(NVEC)
  real(kind=8) :: u795(NVEC)
  real(kind=8) :: u796(NVEC)
  real(kind=8) :: u797(NVEC)
  real(kind=8) :: u798(NVEC)
  real(kind=8) :: u799(NVEC)
  real(kind=8) :: u8(NVEC)
  real(kind=8) :: u80(NVEC)
  real(kind=8) :: u800(NVEC)
  real(kind=8) :: u801(NVEC)
  real(kind=8) :: u802(NVEC)
  real(kind=8) :: u803(NVEC)
  real(kind=8) :: u804(NVEC)
  real(kind=8) :: u805(NVEC)
  real(kind=8) :: u806(NVEC)
  real(kind=8) :: u807(NVEC)
  real(kind=8) :: u808(NVEC)
  real(kind=8) :: u809(NVEC)
  real(kind=8) :: u81(NVEC)
  real(kind=8) :: u810(NVEC)
  real(kind=8) :: u811(NVEC)
  real(kind=8) :: u812(NVEC)
  real(kind=8) :: u813(NVEC)
  real(kind=8) :: u814(NVEC)
  real(kind=8) :: u815(NVEC)
  real(kind=8) :: u816(NVEC)
  real(kind=8) :: u817(NVEC)
  real(kind=8) :: u818(NVEC)
  real(kind=8) :: u819(NVEC)
  real(kind=8) :: u82(NVEC)
  real(kind=8) :: u820(NVEC)
  real(kind=8) :: u821(NVEC)
  real(kind=8) :: u822(NVEC)
  real(kind=8) :: u823(NVEC)
  real(kind=8) :: u824(NVEC)
  real(kind=8) :: u825(NVEC)
  real(kind=8) :: u826(NVEC)
  real(kind=8) :: u827(NVEC)
  real(kind=8) :: u828(NVEC)
  real(kind=8) :: u829(NVEC)
  real(kind=8) :: u83(NVEC)
  real(kind=8) :: u830(NVEC)
  real(kind=8) :: u831(NVEC)
  real(kind=8) :: u832(NVEC)
  real(kind=8) :: u833(NVEC)
  real(kind=8) :: u834(NVEC)
  real(kind=8) :: u835(NVEC)
  real(kind=8) :: u836(NVEC)
  real(kind=8) :: u837(NVEC)
  real(kind=8) :: u838(NVEC)
  real(kind=8) :: u839(NVEC)
  real(kind=8) :: u84(NVEC)
  real(kind=8) :: u840(NVEC)
  real(kind=8) :: u841(NVEC)
  real(kind=8) :: u842(NVEC)
  real(kind=8) :: u843(NVEC)
  real(kind=8) :: u844(NVEC)
  real(kind=8) :: u845(NVEC)
  real(kind=8) :: u846(NVEC)
  real(kind=8) :: u847(NVEC)
  real(kind=8) :: u848(NVEC)
  real(kind=8) :: u849(NVEC)
  real(kind=8) :: u85(NVEC)
  real(kind=8) :: u850(NVEC)
  real(kind=8) :: u851(NVEC)
  real(kind=8) :: u852(NVEC)
  real(kind=8) :: u853(NVEC)
  real(kind=8) :: u854(NVEC)
  real(kind=8) :: u855(NVEC)
  real(kind=8) :: u856(NVEC)
  real(kind=8) :: u857(NVEC)
  real(kind=8) :: u858(NVEC)
  real(kind=8) :: u859(NVEC)
  real(kind=8) :: u86(NVEC)
  real(kind=8) :: u860(NVEC)
  real(kind=8) :: u861(NVEC)
  real(kind=8) :: u862(NVEC)
  real(kind=8) :: u863(NVEC)
  real(kind=8) :: u864(NVEC)
  real(kind=8) :: u865(NVEC)
  real(kind=8) :: u866(NVEC)
  real(kind=8) :: u867(NVEC)
  real(kind=8) :: u868(NVEC)
  real(kind=8) :: u869(NVEC)
  real(kind=8) :: u87(NVEC)
  real(kind=8) :: u870(NVEC)
  real(kind=8) :: u871(NVEC)
  real(kind=8) :: u872(NVEC)
  real(kind=8) :: u873(NVEC)
  real(kind=8) :: u874(NVEC)
  real(kind=8) :: u875(NVEC)
  real(kind=8) :: u876(NVEC)
  real(kind=8) :: u877(NVEC)
  real(kind=8) :: u878(NVEC)
  real(kind=8) :: u879(NVEC)
  real(kind=8) :: u88(NVEC)
  real(kind=8) :: u880(NVEC)
  real(kind=8) :: u881(NVEC)
  real(kind=8) :: u882(NVEC)
  real(kind=8) :: u883(NVEC)
  real(kind=8) :: u884(NVEC)
  real(kind=8) :: u885(NVEC)
  real(kind=8) :: u886(NVEC)
  real(kind=8) :: u887(NVEC)
  real(kind=8) :: u888(NVEC)
  real(kind=8) :: u889(NVEC)
  real(kind=8) :: u89(NVEC)
  real(kind=8) :: u890(NVEC)
  real(kind=8) :: u891(NVEC)
  real(kind=8) :: u892(NVEC)
  real(kind=8) :: u893(NVEC)
  real(kind=8) :: u894(NVEC)
  real(kind=8) :: u895(NVEC)
  real(kind=8) :: u896(NVEC)
  real(kind=8) :: u897(NVEC)
  real(kind=8) :: u898(NVEC)
  real(kind=8) :: u899(NVEC)
  real(kind=8) :: u9(NVEC)
  real(kind=8) :: u90(NVEC)
  real(kind=8) :: u900(NVEC)
  real(kind=8) :: u901(NVEC)
  real(kind=8) :: u902(NVEC)
  real(kind=8) :: u903(NVEC)
  real(kind=8) :: u904(NVEC)
  real(kind=8) :: u905(NVEC)
  real(kind=8) :: u906(NVEC)
  real(kind=8) :: u907(NVEC)
  real(kind=8) :: u908(NVEC)
  real(kind=8) :: u909(NVEC)
  real(kind=8) :: u91(NVEC)
  real(kind=8) :: u910(NVEC)
  real(kind=8) :: u911(NVEC)
  real(kind=8) :: u912(NVEC)
  real(kind=8) :: u913(NVEC)
  real(kind=8) :: u914(NVEC)
  real(kind=8) :: u915(NVEC)
  real(kind=8) :: u916(NVEC)
  real(kind=8) :: u917(NVEC)
  real(kind=8) :: u918(NVEC)
  real(kind=8) :: u919(NVEC)
  real(kind=8) :: u92(NVEC)
  real(kind=8) :: u920(NVEC)
  real(kind=8) :: u921(NVEC)
  real(kind=8) :: u922(NVEC)
  real(kind=8) :: u923(NVEC)
  real(kind=8) :: u924(NVEC)
  real(kind=8) :: u925(NVEC)
  real(kind=8) :: u926(NVEC)
  real(kind=8) :: u927(NVEC)
  real(kind=8) :: u928(NVEC)
  real(kind=8) :: u929(NVEC)
  real(kind=8) :: u93(NVEC)
  real(kind=8) :: u930(NVEC)
  real(kind=8) :: u931(NVEC)
  real(kind=8) :: u932(NVEC)
  real(kind=8) :: u933(NVEC)
  real(kind=8) :: u934(NVEC)
  real(kind=8) :: u935(NVEC)
  real(kind=8) :: u936(NVEC)
  real(kind=8) :: u937(NVEC)
  real(kind=8) :: u938(NVEC)
  real(kind=8) :: u939(NVEC)
  real(kind=8) :: u94(NVEC)
  real(kind=8) :: u940(NVEC)
  real(kind=8) :: u941(NVEC)
  real(kind=8) :: u942(NVEC)
  real(kind=8) :: u943(NVEC)
  real(kind=8) :: u944(NVEC)
  real(kind=8) :: u945(NVEC)
  real(kind=8) :: u946(NVEC)
  real(kind=8) :: u947(NVEC)
  real(kind=8) :: u948(NVEC)
  real(kind=8) :: u949(NVEC)
  real(kind=8) :: u95(NVEC)
  real(kind=8) :: u950(NVEC)
  real(kind=8) :: u951(NVEC)
  real(kind=8) :: u952(NVEC)
  real(kind=8) :: u953(NVEC)
  real(kind=8) :: u954(NVEC)
  real(kind=8) :: u955(NVEC)
  real(kind=8) :: u956(NVEC)
  real(kind=8) :: u957(NVEC)
  real(kind=8) :: u958(NVEC)
  real(kind=8) :: u959(NVEC)
  real(kind=8) :: u96(NVEC)
  real(kind=8) :: u960(NVEC)
  real(kind=8) :: u961(NVEC)
  real(kind=8) :: u962(NVEC)
  real(kind=8) :: u963(NVEC)
  real(kind=8) :: u964(NVEC)
  real(kind=8) :: u965(NVEC)
  real(kind=8) :: u966(NVEC)
  real(kind=8) :: u967(NVEC)
  real(kind=8) :: u968(NVEC)
  real(kind=8) :: u969(NVEC)
  real(kind=8) :: u97(NVEC)
  real(kind=8) :: u970(NVEC)
  real(kind=8) :: u971(NVEC)
  real(kind=8) :: u972(NVEC)
  real(kind=8) :: u973(NVEC)
  real(kind=8) :: u974(NVEC)
  real(kind=8) :: u975(NVEC)
  real(kind=8) :: u976(NVEC)
  real(kind=8) :: u977(NVEC)
  real(kind=8) :: u978(NVEC)
  real(kind=8) :: u979(NVEC)
  real(kind=8) :: u98(NVEC)
  real(kind=8) :: u980(NVEC)
  real(kind=8) :: u981(NVEC)
  real(kind=8) :: u982(NVEC)
  real(kind=8) :: u983(NVEC)
  real(kind=8) :: u984(NVEC)
  real(kind=8) :: u985(NVEC)
  real(kind=8) :: u986(NVEC)
  real(kind=8) :: u987(NVEC)
  real(kind=8) :: u988(NVEC)
  real(kind=8) :: u989(NVEC)
  real(kind=8) :: u99(NVEC)
  real(kind=8) :: u990(NVEC)
  real(kind=8) :: u991(NVEC)
  real(kind=8) :: u992(NVEC)
  real(kind=8) :: u993(NVEC)
  real(kind=8) :: u994(NVEC)
  real(kind=8) :: u995(NVEC)
  real(kind=8) :: u996(NVEC)
  real(kind=8) :: u997(NVEC)
  real(kind=8) :: u998(NVEC)
  real(kind=8) :: u999(NVEC)
  real(kind=8) :: value_1_(NVEC)
  real(kind=8) :: value_2_(NVEC)
  real(kind=8) :: value_3_(NVEC)
  real(kind=8) :: value_4_(NVEC)
  real(kind=8) :: value_5_(NVEC)
  real(kind=8) :: value_6_(NVEC)
  real(kind=8) :: value_7_(NVEC)
  real(kind=8) :: value_8_(NVEC)
  real(kind=8) :: value_9_(NVEC)

  
  ! init values
  value_1_ = 0.0d0
  value_2_ = 0.0d0
  value_3_ = 0.0d0
  value_4_ = 0.0d0
  value_5_ = 0.0d0
  value_6_ = 0.0d0
  value_7_ = 0.0d0
  value_8_ = 0.0d0
  value_9_ = 0.0d0
  
  ! init computation flag
  D_X_Y_shell_optv1_4=.true.
  
  ! compute local quantities
  d=sqrt(sum(r3**2,dim=2))+1.0d-32
  p=sqrt((a12*a3)/(a12+a3))
  pd=p*d
  pd2=pd*pd
  erfpd  = erf (pd)
  exppd2 = exp(-pd2)
  
  ! renormalize xyz
  x=0.0d0
  y=0.0d0
  z=0.0d0
  x=r3(:,1)/d
  y=r3(:,2)/d
  z=r3(:,3)/d
  
  ! computation section
  u3812=A7_v(pd,exppd2,erfpd)
  u954=a3**(-4)
  u765=u954*p**5*pd2*u3812
  u2980=-1.64281880555408736d1*u765
  u3773=+u2980
  u3964=-u2980
  u2980=x**2
  u2885=y**2
  u2934=x*y
  value_1_=value_1_+(c(1)*(u2934*(u3773*u2885+u3964*u2980)))
  u3773=-1.16164831766807941d1*u765
  u3964=+u3773
  u4025=3.48494495300423824d1*u765
  value_2_=value_2_+(c(1)*(y*(u3964*u2885+u4025*u2980)*z))
  u3964=3.72556286454539258d1*u765
  u3819=-6.20927144090898762d0*u765
  u2886=z**2
  value_3_=value_3_+(c(1)*(u2934*(u3964*u2886+u3819*u2885+u3819*u2980))&
 &)
  u3964=1.75624717683788407d1*u765
  u3819=-1.31718538262841305d1*u765
  u2936=y*z
  u3975=(u3964*u2886+u3819*u2885+u3819*u2980)
  value_4_=value_4_+(c(1)*(u2936*u3975))
  u2937=5.55374121304822595d0*u765
  u4191=-1.66612236391446779d1*u765
  u61=2.08265295489308473d0*u765
  u717=4.16530590978616946d0*u765
  u800=y**4
  value_5_=value_5_+(c(1)*(u2886*(u4191*u2980+u4191*u2885+u2937*u2886)+&
 &u61*u800+u2980*(u61*u2980+u717*u2885)))
  u2937=x*z
  value_6_=value_6_+(c(1)*(u2937*u3975))
  u61=-1.86278143227269629d1*u765
  u3819=+u61
  u717=3.10463572045449381d0*u765
  u3964=-u61
  u61=-u717
  u4191=x**4
  value_7_=value_7_+(c(1)*((u3964*u2980+u3819*u2885)*u2886+u717*u800+u6&
 &1*u4191))
  u3819=-u4025
  u3964=-u3773
  value_8_=value_8_+(c(1)*(x*(u3819*u2885+u3964*u2980)*z))
  u3819=4.10704701388521839d0*u765
  u3964=-2.46422820833113103d1*u765
  value_9_=value_9_+(c(1)*(u3819*u800+u2980*(u3819*u2980+u3964*u2885)))
  if ( lmax .eq. 0 ) go to 100
  u3964=u954*p**6*pd
  u3819=u3964*u3812
  u61=1.64281880555408736d1*u3819
  u765=-4.92845641666226206d1*u3819
  u3773=+u765
  u3975=-5.67185232289765129d1*exppd2
  u717=+u3975*pd2
  u4025=u3964*(2.96880505764235461d3*u3812+u717)
  u3812=-4.98024254301440461d-2*u4025
  u3964=4.98024254301440461d-2*u4025
  u3859=u3964*u2980
  value_1_=value_1_+(c(2)*(y*((u61+u3812*u2980)*u2885+u2980*(u3773+u385&
 &9))))
  u3773=-6.96988990600847647d1*u3819
  u3989=+u3773
  u805=-3.5215632741192216d-2*u4025
  u4002=1.05646898223576648d-1*u4025
  u3977=u805*u2885
  u3978=u4002*u2980
  value_2_=value_2_+(c(2)*(u2934*(u3977+u3978+u3989)*z))
  u3989=-3.72556286454539258d1*u3819
  u422=+u3989
  u3888=6.20927144090898762d0*u3819
  u3986=1.86278143227269629d1*u3819
  u407=1.12941284893714386d-1*u4025
  u571=-1.8823547482285731d-2*u4025
  u3860=u571*u2980
  u4130=u3986+u3860
  value_3_=value_3_+(c(2)*(y*((u422+u407*u2980)*u2886+(u3888+u3860)*u28&
 &85+u2980*(u4130))))
  u3987=2.6343707652568261d1*u3819
  u4024=5.32410322828448158d-2*u4025
  u629=-3.99307742121336118d-2*u4025
  u934=u2934*z
  u3811=u629*u2980
  u3979=u4024*u2886
  u3980=u629*u2885
  u870=(u934*(u3979+u3980+u3811+u3987))
  value_4_=value_4_+(c(2)*u870)
  u3820=3.33224472782893557d1*u3819
  u419=1.68362926992343652d-2*u4025
  u3886=-8.33061181957233893d0*u3819
  u4011=-5.05088780977030955d-2*u4025
  u806=6.31360976221288694d-3*u4025
  u4118=1.26272195244257739d-2*u4025
  u639=u4011*u2980+u4011*u2885+u419*u2886
  u3902=u4118*u2980+u806*u2885
  u3981=u806*u2980
  u642=(u2886*(u3820+u639)+u2885*(u3886+u3902)+u2980*(u3886+u3981))
  value_5_=value_5_+(c(2)*(x*u642))
  u3821=-1.75624717683788407d1*u3819
  u527=1.31718538262841305d1*u3819
  u3988=3.95155614788523915d1*u3819
  value_6_=value_6_+(c(2)*(z*((u3821+u4024*u2980)*u2886+(u527+u3811)*u2&
 &885+u2980*(u3988+u3811))))
  u859=-5.6470642446857193d-2*u4025
  u3863=9.41177374114286549d-3*u4025
  u4003=1.24185428818179753d1*u3819
  u3864=5.6470642446857193d-2*u4025
  u3939=-9.41177374114286549d-3*u4025
  u3865=u3864*u2980+u859*u2885
  u3982=u3939*u2980
  value_7_=value_7_+(c(2)*(x*((u422+u3865)*u2886+u3863*u800+u2980*(u400&
 &3+u3982))))
  u3817=3.48494495300423824d1*u3819
  u474=-u3817
  u4179=-1.05646898223576648d-1*u4025
  u72=3.5215632741192216d-2*u4025
  u3861=u72*u2980
  value_8_=value_8_+(c(2)*(((u3817+u4179*u2980)*u2885+u2980*(u474+u3861&
 &))*z))
  u4000=-u765
  u765=1.24506063575360115d-2*u4025
  u3784=-u61
  u61=-7.47036381452160691d-2*u4025
  u4025=u61*u2980+u765*u2885
  u3862=u765*u2980
  value_9_=value_9_+(c(2)*(x*(u2885*(u4000+u4025)+u2980*(u3784+u3862)))&
 &)
  u3983=u3812*u2885
  value_1_=value_1_+(c(3)*(x*(u2885*(u4000+u3859+u3983)+u3784*u2980)))
  u824=u3978+u3977
  u3984=u474*u2980
  value_2_=value_2_+(c(3)*((u2885*(u3817+u824)+u3984)*z))
  u3985=u571*u2885
  value_3_=value_3_+(c(3)*(x*((u422+u407*u2885)*u2886+u2885*(u4130+u398&
 &5)+u3888*u2980)))
  u3888=u3811+u3980
  u3986=u527*u2980
  value_4_=value_4_+(c(3)*(z*((u3821+u4024*u2885)*u2886+u2885*(u3988+u3&
 &888)+u3986)))
  value_5_=value_5_+(c(3)*(y*u642))
  value_6_=value_6_+(c(3)*u870)
  u3886=-u3989
  u3989=-u4003
  u3987=u3863*u2885
  value_7_=value_7_+(c(3)*(y*((u3886+u3865)*u2886+u2885*(u3989+u3987)+u&
 &3939*u4191)))
  u3989=-u3773
  u3988=u4179*u2885
  value_8_=value_8_+(c(3)*(u2934*(u3988+u3861+u3989)*z))
  value_9_=value_9_+(c(3)*(y*(u2885*(u3784+u4025)+u2980*(u4000+u3862)))&
 &)
  value_1_=value_1_+(c(4)*(u2934*(u3983+u3859)*z))
  u3784=1.16164831766807941d1*u3819
  value_2_=value_2_+(c(4)*(y*((u824)*u2886+u3784*u2885+u3984)))
  u805=-7.45112572909078515d1*u3819
  value_3_=value_3_+(c(4)*(u934*(u407*u2886+u3985+u3860+u805)))
  u407=-5.2687415305136522d1*u3819
  u474=(u2886*(u407+u3888+u3979)+u527*u2885+u3986)
  value_4_=value_4_+(c(4)*(y*u474))
  u805=-2.22149648521929038d1*u3819
  value_5_=value_5_+(c(4)*(z*(u2886*(u805+u639)+u2885*(u3820+u3902)+u29&
 &80*(u3820+u3981))))
  value_6_=value_6_+(c(4)*(x*u474))
  value_7_=value_7_+(c(4)*(z*((u3865)*u2886+u2885*(u3886+u3987)+u2980*(&
 &u422+u3982))))
  u805=-u3784
  value_8_=value_8_+(c(4)*(x*((u3861+u3988)*u2886+u3817*u2885+u805*u298&
 &0)))
  value_9_=value_9_+(c(4)*((u765*u800+u2980*(u3862+u61*u2885))*z))
  if ( lmax .eq. 1 ) go to 100
  u805=A9_v(pd,exppd2,erfpd)
  u629=2.96880505764235461d3*u805
  u4011=u954*p**7*pd2
  u859=u4011*(u629-u3975)
  u3939=3.32016169534293641d-2*u859
  u4179=u4011*u805
  u61=4.43561077499603586d2*u4179
  u4024=-1.03497584749907503d3*u4179
  u806=+u4024
  u765=3.26568556340659007d4*u805
  u3988=-1.13437046457953026d2*exppd2
  u3982=+u3988*pd2
  u407=u4011*(u765+u3982)
  u419=-4.98024254301440461d-2*u407
  u4118=4.98024254301440461d-2*u407
  u3863=u419*u2980
  u3864=u4118*u2980
  value_1_=value_1_+(c(5)*(u2934*((u61+u3863)*u2885+u2980*(u806+u3864)+&
 &u3939)))
  u3939=2.3477088494128144d-2*u859
  u806=1.04548348590127147d2*u4179
  u72=-1.56822522885190721d3*u4179
  u422=+u72
  u474=-3.5215632741192216d-2*u407
  u3784=1.05646898223576648d-1*u407
  u3812=u3784*u2980
  value_2_=value_2_+(c(5)*(y*((u806+u474*u2980)*u2885+u2980*(u422+u3812&
 &)+u3939)*z))
  u422=u4011*(-u3975+u629)
  u527=-1.2549031654857154d-2*u422
  u3820=-1.005901973427256d3*u4179
  u4000=+u3820
  u3886=1.67650328904542666d2*u4179
  u3819=3.9118410077726622d2*u4179
  u4002=1.12941284893714386d-1*u407
  u3964=-1.8823547482285731d-2*u407
  u3773=u3964*u2980
  u3989=u4002*u2980
  value_3_=value_3_+(c(5)*(u2934*((u4000+u3989)*u2886+(u3886+u3773)*u28&
 &85+u2980*(u3819+u3773)+u527)))
  u4003=-8.87350538047413596d-3*u422
  u3817=-1.58062245915409566d2*u4179
  u3821=1.18546684436557175d2*u4179
  u642=5.92733422182785873d2*u4179
  u571=5.32410322828448158d-2*u407
  u870=-3.99307742121336118d-2*u407
  u639=u870*u2980
  u3865=u571*u2980
  u3902=(u2936*((u3817+u3865)*u2886+(u3821+u639)*u2885+u2980*(u642+u639&
 &)+u4003))
  value_4_=value_4_+(c(5)*u3902)
  u4025=-5.61209756641145505d-3*u422
  u824=-4.99836709174340336d1*u4179
  u3888=2.07816354034964823d4*u805
  u4130=8.41814634961718258d-3*u4011*(u3888-u3975)
  u4050=-1.68694889346339863d2*u4179
  u4001=1.40302439160286376d-2*u4011*(1.36565032651548312d4*u805-u3975)
  u452=5.99804051009208403d2*u4179
  u4069=1.68362926992343652d-2*u407
  u3861=-4.53748185831812103d2*exppd2
  u3862=+u3861*pd2
  u3889=-1.26272195244257739d-2*u4011*(1.69221888285614213d5*u805+u3862&
 &)
  u3904=5.68224878599159824d-2*u407
  u416=-3.18645902098641964d2*u4179
  u82=-5.05088780977030955d-2*u407
  u3786=6.31360976221288694d-2*u407
  u547=6.31360976221288694d-3*u407
  u3813=u547*u2980
  u3866=u82*u2980
  u3990=u3786*u2980
  u3991=u3904*u2980
  value_5_=value_5_+(c(5)*(u2886*(u4025+u2980*(u3866+u452)+(u4069*u2980&
 &+u824)*u2886)+u2885*(u4130+u2980*(u3990+u3889)+(u3991+u4050)*u2885)+u&
 &2980*(u4001+u2980*(u3813+u416))+u4025))
  u4195=-2.66205161414224079d-2*u422
  u472=-4.74186737746228698d2*u4179
  u3963=3.55640053309671524d2*u4179
  u504=8.29826791055900222d2*u4179
  value_6_=value_6_+(c(5)*(u2937*((u472+u3865)*u2886+(u3963+u639)*u2885&
 &+u2980*(u504+u639)+u4195)))
  u854=-6.274515827428577d-3*u422
  u461=1.8823547482285731d-2*u859
  u3936=1.48440252882117731d4*u805
  u3987=2.83592616144882564d1*exppd2
  u4021=u4011*(u3936+u3987)
  u3811=1.2549031654857154d-2*u4021
  u3874=-1.9559205038863311d2*u4179
  u402=+u3874
  u3981=-7.08981540362206411d0*exppd2
  u3829=u4011*(u629+u3981)
  u3903=5.0196126619428616d-2*u3829
  u966=u4011*(3.85944657493506099d4*u805+u3982)
  u418=-5.6470642446857193d-2*u966
  u3816=6.58824161880000585d-2*u407
  u4004=8.38251644522713329d1*u4179
  u426=5.6470642446857193d-2*u407
  u509=-9.41177374114286549d-3*u407
  u3774=u509*u2980
  u3814=u426*u2980
  value_7_=value_7_+(c(5)*((u461+u2980*(u3814+u4000))*u2886+u2885*(u381&
 &1+u2980*(u3814+u418)+(u3816*u2980+u402)*u2885)+u2980*(u3903+u2980*(u3&
 &774+u4004))+u854))
  u854=9.40935137311144324d2*u4179
  u461=-7.31838440130890029d2*u4179
  u3811=+u461
  u402=-1.05646898223576648d-1*u407
  u3903=3.5215632741192216d-2*u407
  u418=u402*u2980
  u3816=u3903*u2980
  value_8_=value_8_+(c(5)*(x*((u854+u418)*u2885+u2980*(u3811+u3816)+u39&
 &39)*z))
  u3811=-1.6600808476714682d-2*u422
  u895=-3.69634231249669655d1*u4179
  u4065=1.6600808476714682d-2*u859
  u4068=1.10890269374900896d3*u4179
  u3850=1.24506063575360115d-2*u407
  u861=-3.32670808124702689d2*u4179
  u897=-7.47036381452160691d-2*u407
  u56=u3850*u2980
  u3867=u897*u2980
  value_9_=value_9_+(c(5)*(u2885*(u3811+u2980*(u3867+u4068)+(u56+u895)*&
 &u2885)+u2980*(u4065+u2980*(u56+u861))))
  u4192=1.47853692499867862d2*u4179
  u4023=-u4192
  u3868=u4023*u2980
  value_1_=value_1_+(c(6)*(u2885*(u3811+u4118*u4191+(u3863+u4192)*u2885&
 &)+u2980*(u4065+u3868)))
  u3998=-3.13645045770381441d2*u4179
  u853=+u3998
  u3992=u474*u2885
  u553=u3812+u3992
  u3993=u853*u2980
  u470=u3993+u3939
  value_2_=value_2_+(c(6)*(x*(u2885*(u853+u553)+u470)*z))
  u670=1.2549031654857154d-2*u859
  u774=-1.8823547482285731d-2*u4011*(-u3975+u3888)
  u3995=2.25882569787428772d-1*u4011*(u3888+u717)
  u943=-1.31764832376000117d-1*u407
  u3994=u943*u2980
  value_3_=value_3_+(c(6)*(u2885*(u774+u2980*(u3994+u3995)+(u3994+u3819&
 &)*u2885)+u2980*(u774+u3819*u2980)+u670))
  u3995=u870*u2885
  u943=u639+u3995
  u774=u3821*u2980+u4003
  u3869=u571*u2885
  u715=(u2937*((u3817+u3869)*u2886+u2885*(u642+u943)+u774))
  value_4_=value_4_+(c(6)*u715)
  u3872=5.61209756641145505d-3*u859
  u3971=-1.49951012752302101d2*u4179
  u523=1.26272195244257739d-2*u407
  u3870=u82*u2885
  u3996=u4069*u2886
  u844=u3866+u3870+u3996
  u3871=u547*u2885
  u646=u523*u2980+u3871
  value_5_=value_5_+(c(6)*(u2934*(u2886*(u452+u844)+u2885*(u3971+u646)+&
 &u2980*(u3971+u3813)+u3872)))
  value_6_=value_6_+(c(6)*u3902)
  u3872=-1.11766885936361777d2*u4179
  u3971=+u3872
  u3902=-5.6470642446857193d-2*u407
  u392=9.41177374114286549d-3*u407
  u394=-u3872
  u3872=u3902*u2885
  u423=u3814+u3872
  u3775=u392*u2885
  value_7_=value_7_+(c(6)*(u2934*((u423)*u2886+u2885*(u3971+u3775)+u298&
 &0*(u394+u3774))))
  u3971=-2.3477088494128144d-2*u422
  u394=-u3998
  value_8_=value_8_+(c(6)*(y*((u394+u418)*u2885+u2980*(u394+u3816)+u397&
 &1)*z))
  u3998=-3.32016169534293641d-2*u422
  u3999=2.95707384999735724d2*u4179
  u3873=u3850*u2885
  u4039=u3867+u3873
  value_9_=value_9_+(c(6)*(u2934*(u2885*(u3999+u4039)+u2980*(u3999+u56)&
 &+u3998)))
  u3999=-u61
  value_1_=value_1_+(c(7)*(y*((u4192+u3863)*u2885+u2980*(u3999+u3864))*&
 &z))
  u899=-6.27290091540762882d2*u4179
  u541=+u899
  value_2_=value_2_+(c(7)*(u2934*((u541+u553)*u2886+u806*u2885+u470)))
  u3939=2.5098063309714308d-2*u859
  u470=-3.35300657809085332d2*u4179
  u660=+u470
  u820=5.58834429681808886d1*u4179
  u560=-5.02950986713627998d2*u4179
  value_3_=value_3_+(c(7)*(u2936*((u660+u3989)*u2886+(u820+u3773)*u2885&
 &+u2980*(u560+u3773)+u3939)))
  u4022=-2.37093368873114349d2*u4179
  u453=u943+u571*u2886
  u473=(u2934*(u2886*(u4022+u453)+u3821*u2885+u774))
  value_4_=value_4_+(c(7)*u473)
  u774=-2.24483902656458202d-2*u422
  u651=9.99673418348680671d1*u4179
  u962=2.24926519128453151d2*u4179
  u446=(u2886*(u651+u844)+u2885*(u962+u646)+u2980*(u962+u3813)+u774)
  value_5_=value_5_+(c(7)*(u2937*u446))
  u844=-4.43675269023706798d-3*u422
  u646=2.21837634511853399d-2*u859
  u3818=-u3821
  u522=4.43675269023706798d-3*u4011*(u3888+u3988)
  u3989=exppd2*pd2
  u3863=-2.83592616144882564d1*u3989
  u3888=u4011*(8.90641517292706383d3*u805+u3863)
  u805=-1.59723096848534447d-1*u3888
  u604=3.99307742121336118d-2*u407
  u3997=u604*u2980
  u974=u3997+u805
  value_6_=value_6_+(c(7)*(u2886*(u646+u2980*(u639+u4022)+(u3865+u3817)&
 &*u2886)+u2885*(u3821+u2980*(u974)+(u3997+u3818)*u2885)+u522*u2980+u84&
 &4))
  u639=-u470
  u470=-2.23533771872723554d2*u4179
  u3883=+u470
  value_7_=value_7_+(c(7)*(u2937*((u660+u423)*u2886+u2885*(u639+u3775)+&
 &u2980*(u3883+u3774)+u3939)))
  u3883=-1.1738544247064072d-2*u422
  u3890=1.1738544247064072d-2*u859
  u900=1.1738544247064072d-2*u4011*(u765-u3988)
  u938=-4.22587592894306592d-1*u3888
  u4202=-4.18193394360508588d2*u4179
  u834=+u4202
  value_8_=value_8_+(c(7)*((u3890+u2980*(u3816+u541))*u2886+u2885*(u394&
 &+u2980*(u3812+u938)+(u3812+u853)*u2885)+u2980*(u900+u834*u2980)+u3883&
 &))
  value_9_=value_9_+(c(7)*(x*(u2885*(u61+u4039)+u2980*(u4023+u56))*z))
  u541=-u4024
  u900=u3864+u419*u2885
  value_1_=value_1_+(c(8)*(u2934*(u2885*(u541+u900)+u3999*u2980+u3998))&
 &)
  u4024=-u461
  u541=-u854
  u3998=u541*u2980
  value_2_=value_2_+(c(8)*(y*(u2885*(u4024+u553)+u3998+u3971)*z))
  u4024=u3773+u3964*u2885
  u3999=u4002*u2885
  value_3_=value_3_+(c(8)*(u2934*((u4000+u3999)*u2886+u2885*(u3819+u402&
 &4)+u3886*u2980+u527)))
  u4000=u3963*u2980
  value_4_=value_4_+(c(8)*(u2936*((u472+u3869)*u2886+u2885*(u504+u943)+&
 &u4000+u4195)))
  value_5_=value_5_+(c(8)*(u2886*(u4025+u2885*(u3870+u452)+(u4069*u2885&
 &+u824)*u2886)+u2885*(u4001+u2980*(u3991+u3889)+u2885*(u3871+u3990+u41&
 &6))+u2980*(u4130+u4050*u2980)+u4025))
  value_6_=value_6_+(c(8)*u715)
  u3889=6.274515827428577d-3*u859
  u824=-1.8823547482285731d-2*u422
  u4025=-5.0196126619428616d-2*u3829
  u4195=-u3820
  u4050=-u4004
  u416=-1.2549031654857154d-2*u4011*(u3987+u3936)
  u4001=5.6470642446857193d-2*u966
  u472=-u3874
  u3819=-6.58824161880000585d-2*u407
  value_7_=value_7_+(c(8)*((u824+u2885*(u3872+u4195))*u2886+u2885*(u402&
 &5+u2980*(u3819*u2980+u4001)+u2885*(u3775+u3902*u2980+u4050))+u2980*(u&
 &416+u472*u2980)+u3889))
  u416=-u72
  u3889=-u806
  u824=u3816+u402*u2885
  u4050=u3889*u2980+u3971
  value_8_=value_8_+(c(8)*(x*(u2885*(u416+u824)+u4050)*z))
  value_9_=value_9_+(c(8)*(u2885*(u4065+u2980*(u56+u4068)+u2885*(u3873+&
 &u3867+u861))+u2980*(u3811+u895*u2980)))
  value_1_=value_1_+(c(9)*(x*(u2885*(u61+u900)+u3868)*z))
  u416=-1.1738544247064072d-2*u4011*(-u3988+u765)
  u4065=-u899
  u4001=-u4202
  u3811=4.22587592894306592d-1*u3888
  value_2_=value_2_+(c(9)*((u3883+u2885*(u3992+u4065))*u2886+u2885*(u41&
 &6+u2980*(u418+u3811)+(u418+u4001)*u2885)+u2980*(u853+u394*u2980)+u389&
 &0))
  u3811=u820*u2980+u3939
  value_3_=value_3_+(c(9)*(u2937*((u660+u3999)*u2886+u2885*(u560+u4024)&
 &+u3811)))
  value_4_=value_4_+(c(9)*(u2886*(u646+u2885*(u3995+u4022)+(u3869+u3817&
 &)*u2886)+u2885*(u522+u2980*(u604*u2885+u974))+u2980*(u3821+u3818*u298&
 &0)+u844))
  value_5_=value_5_+(c(9)*(u2936*u446))
  value_6_=value_6_+(c(9)*u473)
  u4001=-2.5098063309714308d-2*u422
  u604=-u470
  value_7_=value_7_+(c(9)*(u2936*((u639+u423)*u2886+u2885*(u604+u3775)+&
 &u2980*(u660+u3774)+u4001)))
  u4001=u394*u2885
  value_8_=value_8_+(c(9)*(u2934*((u4065+u824)*u2886+u4001+u4050)))
  value_9_=value_9_+(c(9)*(y*(u2885*(u4023+u4039)+u2980*(u61+u56))*z))
  value_1_=value_1_+(c(10)*(u2934*((u900)*u2886+u4192*u2885+u3868)))
  value_2_=value_2_+(c(10)*(u2936*((u553)*u2886+u4001+u3998)))
  u604=-1.67650328904542666d3*u4179
  value_3_=value_3_+(c(10)*(u2934*(u2886*(u604+u4024+u4002*u2886)+u820*&
 &u2885+u3811)))
  u604=3.54940215218965439d-2*u859
  u3811=-1.10643572140786696d3*u4179
  u897=(u2886*(u3811+u453)+u3963*u2885+u4000+u604)
  value_4_=value_4_+(c(10)*(u2936*u897))
  u844=1.40302439160286376d-3*u859
  u3883=2.10453658740429565d-2*u859
  u474=-4.49853038256906302d2*u4179
  u3971=-2.52544390488515477d-2*u4011*(u3987+u629)
  u523=7.87242816949586029d2*u4179
  u774=1.87438765940377626d1*u4179
  u4118=5.05088780977030955d-2*u3888
  u4003=-1.26272195244257739d-2*u407
  u4002=u4003*u2980
  value_5_=value_5_+(c(10)*(u2886*(u3883+u2980*(u3813+u523)+u2885*(u387&
 &1+u523)+u2886*(u3996+u3870+u3866+u474))+u2885*(u3971+u2980*(u4002+u41&
 &18)+(u4002+u774)*u2885)+u2980*(u3971+u774*u2980)+u844))
  value_6_=value_6_+(c(10)*(u2937*u897))
  u4069=8.38251644522713329d2*u4179
  u4003=-2.79417214840904443d1*u4179
  u474=+u4003
  u4118=-u4069
  u547=-u4003
  value_7_=value_7_+(c(10)*(u2886*(u2980*(u3774+u4118)+u2885*(u3775+u40&
 &69)+(u3872+u3814)*u2886)+u2885*(u527+u474*u2885)+u2980*(u670+u547*u29&
 &80)))
  value_8_=value_8_+(c(10)*(u2937*((u824)*u2886+u854*u2885+u3993)))
  u474=-8.30040423835734101d-3*u422
  u4118=8.30040423835734101d-3*u859
  u547=1.6600808476714682d-2*u4021
  u4069=-2.21780538749801793d2*u4179
  u426=-2.58743961874768758d2*u4179
  u3902=-2.98814552580864276d-1*u3888
  u3774=7.47036381452160691d-2*u407
  u4003=u3774*u2980
  value_9_=value_9_+(c(10)*((u4118+u2980*(u56+u4069)+u2885*(u3873+u4069&
 &))*u2886+u2885*(u547+u2980*(u4003+u3902)+(u4003+u426)*u2885)+u2980*(u&
 &547+u426*u2980)+u474))
  if ( lmax .eq. 2 ) go to 100
  u474=A11_v(pd,exppd2,erfpd)
  u4118=pd2*u474
  u3902=3.26568556340659007d4*u4118
  u3774=2.0d0*pd2
  u3850=(-u3975*(u3774+1.1d1)+u3902)
  u547=u954*p**8*pd
  u4069=u547*u3850
  u426=-3.01832881394812401d-3*u4069
  u3903=3.26568556340659007d4*u474
  u392=pd2*(-u3988+u3903)
  u3971=u547*u392
  u604=-1.3582479662766558d-2*u3971
  u527=pd2*(u3903-u3988)
  u3811=u547*u527
  u523=1.22242316964899022d-1*u3811
  u774=u547*u4118
  u844=9.75834370499127888d3*u774
  u3883=-1.95166874099825578d4*u774
  u897=+u3883
  u509=4.24539123242856709d5*u474
  u4003=-2.26874092915906051d2*exppd2
  u3873=+u4003*pd2
  u3939=pd2*(u509+u3873)
  u4065=u547*u3939
  u670=-4.98024254301440461d-2*u4065
  u646=4.98024254301440461d-2*u4065
  u522=u646*u2980
  u3890=u670*u2980
  value_1_=value_1_+(c(11)*(y*((u604+u2980*(u3890+u844))*u2885+u2980*(u&
 &523+u2980*(u522+u897))+u426)))
  u523=1.15251161698447252d-1*u3811
  u897=3.45009550347419585d3*u774
  u3889=-3.10508595312677627d4*u774
  u870=+u3889
  u859=-3.5215632741192216d-2*u4065
  u422=1.05646898223576648d-1*u4065
  u4021=u422*u2980
  u3817=u859*u2980
  value_2_=value_2_+(c(11)*(u2934*((u897+u3817)*u2885+u2980*(u870+u4021&
 &)+u523)*z))
  u523=(u3902-u3975*(1.1d1+u3774))
  u870=u547*u523
  u853=1.14082105953246854d-3*u870
  u4023=3.08021686073766507d-2*u3811
  u660=-5.13369476789610845d-3*u3971
  u4022=-4.62032529110649761d-2*u3971
  u824=-2.21298434153996319d4*u774
  u4050=3.68830723589993865d3*u774
  u416=7.3766144717998773d3*u774
  u472=1.12941284893714386d-1*u4065
  u895=-1.8823547482285731d-2*u4065
  u861=u895*u2980
  u541=u861+u416
  u560=u2980*(u541)
  u3818=u472*u2980
  value_3_=value_3_+(c(11)*(y*((u4023+u2980*(u3818+u824))*u2886+(u660+u&
 &2980*(u861+u4050))*u2885+u2980*(u4022+u560)+u853)))
  u900=-4.35608445950548493d-2*u3971
  u4025=-5.21605411520851568d3*u774
  u4195=3.91204058640638676d3*u774
  u82=1.17361217592191603d4*u774
  u402=5.32410322828448158d-2*u4065
  u3963=-3.99307742121336118d-2*u4065
  u394=u3963*u2980
  u3819=u402*u2980
  u642=(u934*((u4025+u3819)*u2886+(u4195+u394)*u2885+u2980*(u82+u394)+u&
 &900))
  value_4_=value_4_+(c(11)*u642)
  u820=1.63284278170329504d5*u4118
  u962=2.0d1*pd2
  u452=-3.06114412713352094d-3*u547*(u3987*(u962+1.1d1)+u820)
  u3886=-4.13254457163025327d-2*u3971
  u504=-1.64946114027532311d3*u774
  u4068=+u504
  u651=pd2*(u509-u4003)
  u639=u547*u651
  u4130=1.37751485721008442d-2*u639
  u805=-5.56693134842921549d3*u774
  u3784=1.8142697574481056d4*u474
  u4179=u547*pd2
  u629=1.65301782865210131d-1*u4179*(u3784+u3987)
  u765=1.31956891222025849d4*u774
  u3936=1.68362926992343652d-2*u4065
  u4011=-1.26272195244257739d-2*u4179
  u3872=-9.07496371663624206d2*u3989
  u3888=+u4011*(2.51457788382307436d6*u474+u3872)
  u407=5.68224878599159824d-2*u4065
  u571=-4.74220077829155393d3*u774
  u61=-5.05088780977030955d-2*u4065
  u4024=6.31360976221288694d-2*u4065
  u806=6.31360976221288694d-3*u4065
  u72=u806*u2980
  u3820=u61*u2980
  u3821=u407*u2980
  u3874=u4024*u2980
  u4004=u3936*u2980
  value_5_=value_5_+(c(11)*(x*(u2886*(u3886+u2980*(u3820+u765)+(u4004+u&
 &4068)*u2886)+u2885*(u4130+u2980*(u3874+u3888)+(u3821+u805)*u2885)+u29&
 &80*(u629+u2980*(u72+u571))+u452)))
  u854=2.42004692194749163d-3*u870
  u461=1.45202815316849498d-2*u3811
  u4192=-1.08902111487637123d-2*u3971
  u899=-9.80119003388734109d-2*u3971
  u3992=-1.04321082304170314d4*u774
  u4202=7.82408117281277352d3*u774
  u3964=1.5648162345625547d4*u774
  u419=u2980*(u394+u3964)
  value_6_=value_6_+(c(11)*(z*((u461+u2980*(u3819+u3992))*u2886+(u4192+&
 &u2980*(u394+u4202))*u2885+u2980*(u899+u419)+u854)))
  u446=7.51107679583515716d5*u4118
  u3814=2.3d1*pd2
  u715=-u3988*(u3814-2.2d1)
  u473=(u715+u446)
  u966=u547*u473
  u453=-5.70410529766234272d-4*u966
  u553=7.70054215184416268d-2*u3811
  u943=9.79705669021977021d4*u474
  u423=pd2*(u943+u3987)
  u4039=u547*u423
  u938=6.16043372147533014d-2*u4039
  u834=-6.45453766282489263d3*u774
  u3904=+u834
  u974=u4179*(1.95941133804395404d4*u474+u3975)
  u3786=5.13369476789610845d-2*u974
  u3829=-1.84415361794996932d4*u774
  u56=+u3829
  u470=5.55166545779120312d5*u474
  u418=pd2*(u470+u3873)
  u494=u547*u418
  u569=-5.6470642446857193d-2*u494
  u464=6.58824161880000585d-2*u4065
  u4012=2.76623042692495399d3*u774
  u776=5.6470642446857193d-2*u4065
  u395=-9.41177374114286549d-3*u4065
  u66=u395*u2980
  u4018=u776*u2980
  u4005=u464*u2980
  value_7_=value_7_+(c(11)*(x*((u553+u2980*(u4018+u56))*u2886+u2885*(u9&
 &38+u2980*(u4018+u569)+(u4005+u3904)*u2885)+u2980*(u3786+u2980*(u66+u4&
 &012))+u453)))
  u453=-2.13428077219346764d-3*u4069
  u553=-2.88127904246118131d-2*u3971
  u938=8.64383712738354392d-2*u3811
  u3786=2.07005730208451751d4*u774
  u569=-1.38003820138967834d4*u774
  u3782=+u569
  u468=-1.05646898223576648d-1*u4065
  u3905=3.5215632741192216d-2*u4065
  u667=u3905*u2980
  u3995=u468*u2980
  value_8_=value_8_+(c(11)*(((u553+u2980*(u3995+u3786))*u2885+u2980*(u9&
 &38+u2980*(u667+u3782))+u453)*z))
  u938=-8.14948779765993481d-2*u3971
  u4188=-1.21979296312390986d3*u774
  u654=+u4188
  u554=5.43299186510662321d-2*u3811
  u505=2.19562733362303775d4*u774
  u4167=1.24506063575360115d-2*u4065
  u3783=-6.09896481561954931d3*u774
  u887=-7.47036381452160691d-2*u4065
  u797=u4167*u2980
  u3776=u887*u2980
  value_9_=value_9_+(c(11)*(x*(u2885*(u938+u2980*(u3776+u505)+(u797+u65&
 &4)*u2885)+u2980*(u554+u2980*(u797+u3783))+u426)))
  u4121=4.87917185249563944d3*u774
  u873=3.16924525464553021d-2*u3811
  u818=-6.50556246999418592d3*u774
  u4066=+u818
  u809=-1.62639061749854648d3*u774
  u3837=+u809
  u4006=u3837*u2980
  value_1_=value_1_+(c(12)*(x*(u2885*(u604+u2980*(u522+u4066)+(u3890+u4&
 &121)*u2885)+u2980*(u873+u4006)+u426)))
  u873=9.60426347487060435d-3*u3811
  u787=1.15003183449139862d3*u774
  u3778=4.80213173743530218d-2*u3811
  u54=-u897
  u3875=u54*u2980
  value_2_=value_2_+(c(12)*((u2885*(u873+u2980*(u4021+u3782)+(u3817+u78&
 &7)*u2885)+u2980*(u3778+u3875)+u453)*z))
  u41=2.28597989438461305d5*u4118
  u832=1.41796308072441282d1*exppd2
  u3993=5.6d1*pd2
  u855=4.56328423812987418d-3*u547*(u41+u832*(1.1d1+u3993))
  u62=2.28597989438461305d5*u474
  u3779=-5.6470642446857193d-2*u4179*(-u3988+u62)
  u3788=1.29090753256497853d4*u774
  u453=5.67185232289765129d2*exppd2
  u4010=-8.55615794649351409d-3*u4179*(u453+5.94354772539999393d5*u474)
  u3833=7.6199329812820435d4*u474
  u903=pd2*(u3833+u3863)
  u780=9.03530279149715087d-1*u547*u903
  u858=-1.31764832376000117d-1*u4065
  u3879=4.30302510854992842d3*u774
  u3876=u858*u2980
  u3780=u2980*(u3876+u780)
  value_3_=value_3_+(c(12)*(x*(u2885*(u3779+u3780+(u3876+u3788)*u2885)+&
 &u2980*(u4010+u3879*u2980)+u855)))
  u856=8.06682307315830542d-4*u547*(6.20480257047252114d5*u4118-u3975*(&
 &1.1d1+3.8d1*pd2))
  u499=pd2*(-u4003+u509)
  u509=u547*u499
  u3793=-4.84009384389498325d-3*u509
  u804=1.73868470506950523d3*u774
  u857=-3.2670633446291137d-2*u3971
  u807=1.30401352880212892d3*u774
  u4064=-3.2670633446291137d-2*u4179*(-u3988+1.92312594289499193d5*u474&
 &)
  u540=pd2*(u3903-1.41796308072441282d1*u3989)
  u998=u547*u540
  u406=8.51856516525517053d-1*u998
  u575=-5.32410322828448158d-2*u4065
  u4034=6.5200676440106446d3*u774
  u4007=u575*u2980
  u4070=(z*(u2886*(u3793+u2980*(u4007+u406)+(u4007+u804)*u2886)+u2885*(&
 &u857+u419+(u394+u807)*u2885)+u2980*(u4064+u4034*u2980)+u856))
  value_4_=value_4_+(c(12)*u4070)
  u419=(u832*(u3993+1.1d1)+u41)
  u454=-2.04076275142234729d-3*u547*u419
  u3892=-4.59171619070028141d-3*u3971
  u3880=-5.49820380091774369d2*u774
  u3991=4.89852834510988511d5*u474
  u4183=pd2*(u3991-u3861)
  u530=u547*u4183
  u3823=4.59171619070028141d-3*u530
  u627=-1.8556437828097385d3*u774
  u965=u4179*(2.33263254529042148d4*u474+u832)
  u412=2.57136106679215759d-1*u965
  u3773=6.59784456110129243d3*u774
  u95=+u4011*(2.44926417255494255d6*u474+u3872)
  u4011=-7.62875777377336937d3*u774
  value_5_=value_5_+(c(12)*(y*(u2886*(u3892+u2980*(u3820+u3773)+(u4004+&
 &u3880)*u2886)+u2885*(u3823+u2980*(u3874+u95)+(u3821+u627)*u2885)+u298&
 &0*(u412+u2980*(u72+u4011))+u454)))
  value_6_=value_6_+(c(12)*u642)
  u642=-1.54010843036883254d-2*u3971
  u3885=1.54010843036883254d-2*u3811
  u408=u4179*(u3903+u3987)
  u4117=7.52941899291429239d-2*u408
  u3940=-2.15151255427496421d3*u774
  u833=+u3940
  u661=pd2*(u470-u3975)
  u999=u547*u661
  u3824=1.02673895357922169d-2*u999
  u3897=-1.10649217076998159d4*u774
  u813=+u3897
  u826=pd2*(5.76937782868497579d5*u474+u3873)
  u4013=u547*u826
  u823=-5.6470642446857193d-2*u4013
  u3826=-u4012
  value_7_=value_7_+(c(12)*(y*((u3885+u2980*(u4018+u813))*u2886+u2885*(&
 &u4117+u2980*(u4018+u823)+(u4005+u833)*u2885)+u2980*(u3824+u2980*(u66+&
 &u3826))+u642)))
  u823=-3.84170538994824174d-2*u3971
  u833=1.03502865104225876d4*u774
  u3824=-u787
  value_8_=value_8_+(c(12)*(u2934*((u833+u3995)*u2885+u2980*(u3824+u667&
 &)+u823)*z))
  u823=3.01832881394812401d-3*u870
  u4117=-9.05498644184437202d-3*u3971
  u464=-4.0659765437463662d2*u774
  u585=-5.43299186510662321d-2*u3971
  u91=1.05715390137405521d4*u774
  u961=-u4188
  u4188=(u797+u464)*u2885
  value_9_=value_9_+(c(12)*(y*(u2885*(u4117+u2980*(u3776+u91)+u4188)+u2&
 &980*(u585+u2980*(u797+u961))+u823)))
  u4004=2.7164959325533116d-2*u3811
  u3967=-1.13847343224898254d4*u774
  u4026=+u3967
  value_1_=value_1_+(c(13)*(u2934*((u4121+u3890)*u2885+u2980*(u4026+u52&
 &2)+u4004)*z))
  u4004=exppd2*(pd2+1.0d0)
  u4026=-1.1738544247064072d-2*u547*(1.13437046457953026d2*u4004+u3902)
  u403=2.88127904246118131d-2*u3811
  u860=3.20142115829020145d-3*u639
  u871=1.85055515259706771d5*u474
  u3887=2.88127904246118131d-2*u4179*(u871-u4003)
  u524=-u3786
  u3834=-5.63450123859075456d-1*u998
  u386=-6.9001910069483917d3*u774
  u4061=+u386
  u4008=u4061*u2980
  value_2_=value_2_+(c(13)*(y*((u403+u2980*(u4021+u524))*u2886+u2885*(u&
 &860+u2980*(u667+u3834)+(u667+u3824)*u2885)+u2980*(u3887+u4008)+u4026)&
 &))
  u4026=5.13369476789610845d-2*u3811
  u860=1.84415361794996932d3*u774
  u3887=-3.07358936324994888d3*u774
  value_3_=value_3_+(c(13)*(u934*((u813+u3818)*u2886+(u860+u861)*u2885+&
 &u2980*(u3887+u861)+u4026)))
  u3834=8.16421390851647518d5*u4118
  u4007=2.5d1*pd2
  u808=-u3988*(u4007-1.1d1)
  u457=u547*(u808+u3834)
  u427=-4.03341153657915271d-4*u457
  u937=1.81503519146061872d-2*u3811
  u4067=-u804
  u3871=3.63007038292123744d-3*u639
  u628=-u807
  u3901=1.01236252465604292d6*u474
  u3900=pd2*(u3901+u4003)
  u997=3.63007038292123744d-3*u547
  u668=u997*u3900
  u404=-2.60802705760425784d3*u774
  u4112=+u404
  u944=-6.38892387394137789d-1*u998
  u968=3.99307742121336118d-2*u4065
  u3777=u968*u2980
  u4205=u3777+u944
  u613=u2980*(u4205)
  u784=(y*(u2886*(u937+u2980*(u394+u4112)+(u3819+u4067)*u2886)+u2885*(u&
 &3871+u613+(u3777+u628)*u2885)+u2980*(u668+u4112*u2980)+u427))
  value_4_=value_4_+(c(13)*u784)
  u455=-5.10190687855586823d-4*u966
  u966=8.57120355597385863d-2*u965
  u965=-2.19928152036709748d3*u774
  u4005=6.88757428605042211d-3*u3811
  u465=-2.06182642534415388d2*u774
  u985=+u465
  u794=1.63284278170329504d5*u474
  u456=2.52544390488515477d-2*u4179*(u794+u3988)
  u875=pd2*(5.44280927234431679d4*u474+u3863)
  u451=-4.04071024781624764d-1*u547*u875
  u3838=6.73451707969374607d-2*u4065
  u440=-3.71128756561947699d3*u774
  u498=-3.5051049230850616d3*u774
  u459=1.26272195244257739d-2*u4065
  u3877=u459*u2980
  value_5_=value_5_+(c(13)*(z*(u2886*(u966+u451*u2980+(u3838*u2980+u965&
 &)*u2886)+u2885*(u4005+u2980*(u3877+u440)+(u72+u985)*u2885)+u2980*(u45&
 &6+u2980*(u72+u498))+u455)))
  u4005=7.0d0*pd2
  u985=-u3988*(u4005-1.1d1)
  u456=-1.21002346097374581d-3*u547*(u985+u41)
  u451=3.2670633446291137d-2*u3811
  u3838=-u4195
  u498=pd2*(u62+u3861)
  u4162=u997*u498
  u997=-u404
  value_6_=value_6_+(c(13)*(x*(u2886*(u451+u2980*(u394+u997)+(u3819+u40&
 &25)*u2886)+u2885*(u4195+u613+(u3777+u3838)*u2885)+u4162*u2980+u456)))
  u404=1.01236252465604292d6*u4118
  u546=-5.70410529766234272d-4*u547
  u3775=3.1d1*pd2
  u803=+u546*(-u3988*(u3775+2.2d1)+u404)
  u415=6.53137112681318014d3*u474
  u963=4.10695581431688676d-1*u4179*(u415-u3981)
  u3819=-u860
  u4126=5.13369476789610845d-3*u3811
  u754=-3.07358936324994887d2*u774
  u460=+u754
  u812=2.25882569787428772d-1*u408
  u614=-1.12941284893714386d-1*u4179*(3.59225411974724908d5*u474+u3982)
  u3899=9.41177374114286549d-3*u4065
  value_7_=value_7_+(c(13)*(z*(u2886*(u963+u2980*(u3818+u614)+(u4018+u3&
 &819)*u2886)+u2885*(u4126+u3819*u2980+(u3899*u2980+u460)*u2885)+u2980*&
 &(u812+u2980*(u66+u3904))+u803)))
  u803=9.47048813387911121d5*u4118
  u963=2.9d1*pd2
  u4126=(-u3988*(u963+1.1d1)+u803)
  u812=-1.06714038609673382d-3*u547*u4126
  u614=-u833
  u3904=1.60018592606922914d6*u474
  u3819=1.81499274332724841d3*exppd2
  u400=3.20142115829020145d-3*u4179*(u3904+u3819)
  u839=-1.15003183449139862d4*u774
  u4020=+u839
  u4136=-1.69035037157722637d0*u998
  u4015=+u4136
  u607=-4.60012733796559447d3*u774
  u398=+u607
  u3881=u2980*(u4021+u4015)
  value_8_=value_8_+(c(13)*(x*((u3778+u2980*(u667+u4020))*u2886+u2885*(&
 &u833+u3881+(u4021+u614)*u2885)+u2980*(u400+u398*u2980)+u812)))
  u812=1.3582479662766558d-2*u3811
  u400=1.21979296312390986d4*u774
  u398=-3.65937888937172958d3*u774
  value_9_=value_9_+(c(13)*((u2885*(u604+u2980*(u3776+u400)+u4188)+u298&
 &0*(u812+u2980*(u797+u398)))*z))
  u3778=-3.16924525464553021d-2*u3971
  u4188=-u809
  u809=-u818
  u818=-u4121
  u4014=(u3890+u4188)*u2885
  u3878=u818*u2980
  value_1_=value_1_+(c(14)*(y*(u2885*(u3778+u2980*(u522+u809)+u4014)+u2&
 &980*(u812+u3878)+u823)))
  u3778=3.84170538994824174d-2*u3811
  u3822=u859*u2885
  u827=u4021+u3822
  u4009=u614*u2980
  value_2_=value_2_+(c(14)*(u2934*(u2885*(u787+u827)+u4009+u3778)*z))
  value_3_=value_3_+(c(14)*(y*(u2885*(u4010+u3780+(u3876+u3879)*u2885)+&
 &u2980*(u3779+u3788*u2980)+u855)))
  u3778=u3963*u2885
  u858=u394+u3778
  u3779=u402*u2885
  u4010=u4195*u2980
  u780=(u934*((u4025+u3779)*u2886+u2885*(u82+u858)+u4010+u900))
  value_4_=value_4_+(c(14)*u780)
  u3788=u806*u2885
  u855=u3788+u3874
  u3780=u61*u2885
  u3879=u3936*u2885
  value_5_=value_5_+(c(14)*(x*(u2886*(u3892+u2885*(u3780+u3773)+(u3879+&
 &u3880)*u2886)+u2885*(u412+u2980*(u3821+u95)+u2885*(u855+u4011))+u2980&
 &*(u3823+u627*u2980)+u454)))
  value_6_=value_6_+(c(14)*u4070)
  u575=-1.02673895357922169d-2*u4179*(-u3975+u470)
  u95=-u3897
  u3793=-5.6470642446857193d-2*u4065
  u857=-7.52941899291429239d-2*u4179*(u3987+u3903)
  u4064=5.6470642446857193d-2*u4013
  u406=-u3940
  u3892=-6.58824161880000585d-2*u4065
  u3823=u3793*u2885
  u412=u2885*(u3823+u95)
  u454=u3899*u2885
  u3880=u3793*u2980
  u627=u454+u3880
  u4011=u3892*u2980
  value_7_=value_7_+(c(14)*(x*((u642+u412)*u2886+u2885*(u575+u2980*(u40&
 &11+u4064)+u2885*(u627+u4012))+u2980*(u857+u406*u2980)+u3885)))
  u857=2.13428077219346764d-3*u870
  u4064=-4.80213173743530218d-2*u3971
  u406=-9.60426347487060435d-3*u3971
  u3885=-u569
  u4034=(u3995+u897)
  value_8_=value_8_+(c(14)*((u2885*(u4064+u2980*(u667+u3885)+u4034*u288&
 &5)+u2980*(u406+u3824*u2980)+u857)*z))
  u3824=u4167*u2885
  u3773=u3824+u3776
  u4012=u464*u2980
  value_9_=value_9_+(c(14)*(x*(u2885*(u585+u2980*(u797+u91)+u2885*(u377&
 &3+u961))+u2980*(u4117+u4012)+u823)))
  value_1_=value_1_+(c(15)*((u2885*(u604+u646*u4191+u4014)+u2980*(u812+&
 &u4006))*z))
  u4117=u547*(u3834+u808)
  u585=1.06714038609673382d-3*u4117
  u91=1.14298994719230652d6*u474
  u961=pd2*(-u4003+u91)
  u3834=-9.60426347487060435d-3*u547*u961
  u808=-u386
  u569=-u839
  u804=-9.60426347487060435d-3*u509
  u807=-u4136
  u3940=u2980*(u3995+u807)
  u3897=u2885*(u3822+u808)
  value_2_=value_2_+(c(15)*(x*((u406+u3897)*u2886+u2885*(u3834+u3940+(u&
 &3995+u569)*u2885)+u2980*(u804+u897*u2980)+u585)))
  u585=1.14082105953246854d-3*u4117
  u3834=-1.02673895357922169d-2*u509
  u804=6.14717872649989775d2*u774
  u4117=7.51107679583515716d5*u474
  u386=pd2*(-u3988+u4117)
  u839=-1.54010843036883254d-2*u547*u386
  u4136=1.80706055829943018d0*u998
  u4070=-1.12941284893714386d-1*u4065
  u4014=1.16796395803498057d4*u774
  u4013=u4070*u2980
  value_3_=value_3_+(c(15)*(z*(u2886*(u3834+u2980*(u4013+u4136)+(u4013+&
 &u4050)*u2886)+u2885*(u642+u560+(u861+u804)*u2885)+u2980*(u839+u4014*u&
 &2980)+u585)))
  u585=(x*(u2886*(u937+u2885*(u3778+u4112)+(u3779+u4067)*u2886)+u2885*(&
 &u668+u613+(u3777+u4112)*u2885)+u2980*(u3871+u628*u2980)+u427))
  value_4_=value_4_+(c(15)*u585)
  u642=-3.21420133349019699d-2*u3971
  u839=4.39856304073419495d3*u774
  u4070=-u504
  u4014=u3936*u2886
  value_5_=value_5_+(c(15)*(u934*(u2886*(u839+u3820+u3780+u4014)+u2885*&
 &(u4070+u3877+u3788)+u2980*(u4070+u72)+u642)))
  value_6_=value_6_+(c(15)*u784)
  u642=2.4588714905999591d3*u774
  u4070=-u642
  value_7_=value_7_+(c(15)*(u934*((u4018+u3823)*u2886+u2885*(u642+u454)&
 &+u2980*(u4070+u66))))
  u4070=-1.06714038609673382d-3*u457
  u457=9.60426347487060435d-3*u639
  u642=pd2*(u91-u4003)
  u4136=u547*u642
  u504=9.60426347487060435d-3*u4136
  u3834=(u4021+u54)
  value_8_=value_8_+(c(15)*(y*((u873+u2980*(u667+u4061))*u2886+u2885*(u&
 &457+u3881+u3834*u2885)+u2980*(u504+u4020*u2980)+u4070)))
  u4070=-2.7164959325533116d-2*u3971
  u504=3.25278123499709296d3*u774
  value_9_=value_9_+(c(15)*(u2934*(u2885*(u504+u3776+u3824)+u2980*(u504&
 &+u797)+u4070)*z))
  u504=4.52749322092218601d-3*u639
  u4015=1.3582479662766558d-2*u639
  u4020=-u844
  u457=-7.96838806882304737d-1*u998
  value_1_=value_1_+(c(16)*(y*((u812+u2980*(u522+u4020))*u2886+u2885*(u&
 &504+u2980*(u522+u457)+(u522+u3837)*u2885)+u2980*(u4015+u4066*u2980)+u&
 &604)))
  u457=5.76255808492236261d-2*u3811
  u4015=u897*u2885
  value_2_=value_2_+(c(16)*(u934*((u4061+u827)*u2886+u4015+u4009+u457))&
 &)
  u457=+u546*(-u3988*(1.3d1*pd2+2.2d1)+4.24539123242856709d5*u4118)
  u504=5.6470642446857193d-2*u3811
  u3837=-u4050
  u4066=1.71123158929870282d-3*u639
  u546=-u804
  u804=5.13369476789610845d-3*u530
  u530=-3.01176759716571696d-1*u998
  u784=1.8823547482285731d-2*u4065
  u560=-1.22943574529997955d3*u774
  u3881=u784*u2980
  u827=u2980*(u3881+u530)
  value_3_=value_3_+(c(16)*(y*(u2886*(u504+u2980*(u861+u56)+(u3818+u383&
 &7)*u2886)+u2885*(u4066+u827+(u3881+u546)*u2885)+u2980*(u804+u560*u298&
 &0)+u457)))
  u898=7.26014076584247488d-3*u3811
  u411=-9.56276587788227874d3*u774
  u4016=u402*u2886
  u4017=u4195*u2885
  u3898=(u934*(u2886*(u411+u858+u4016)+u4017+u4010+u898))
  value_4_=value_4_+(c(16)*u3898)
  u4010=4.3d1*pd2
  u858=1.27547671963896706d-4*u547*(1.40424479226483373d6*u4118-u3988*(&
 &8.8d1+u4010))
  u877=-3.09940842872268995d-2*u3971
  u793=u4179*(-u3975+u3833)
  u874=-2.06627228581512663d-2*u793
  u896=8.65967098644544631d3*u774
  u810=1.03091321267207694d3*u774
  u4175=-4.13254457163025327d-2*u4179*(u3987+u3784)
  u3784=7.83494041630778476d3*u774
  u828=2.02035512390812382d-1*u998
  u967=-1.26272195244257739d-2*u4065
  u409=-u465
  u465=u4014+u3780
  u840=u2886*(u465+u3820+u4068)
  u3882=u967*u2980
  u475=u2980*(u3882+u828)
  value_5_=value_5_+(c(16)*(x*(u2886*(u877+u2980*(u72+u3784)+u2885*(u37&
 &88+u896)+u840)+u2885*(u874+u475+(u3882+u810)*u2885)+u2980*(u4175+u409&
 &*u2980)+u858)))
  u4060=1.14298994719230652d6*u4118
  u3921=3.5d1*pd2
  u458=-4.03341153657915271d-4*u547*(-u3988*(u3921+4.4d1)+u4060)
  u393=2.42004692194749163d-3*u4179*(u91+2.09858535947213098d3*exppd2)
  u4114=-3.04269823387163415d3*u774
  u4113=3.63007038292123744d-3*u4136
  u4136=3.80996649064102175d5*u474
  u4042=pd2*(u4136+u3982)
  u3781=-7.98615484242672237d-2*u547*u4042
  u4013=9.31718064949784276d-2*u4065
  value_6_=value_6_+(c(16)*(z*(u2886*(u393+u3781*u2980+(u4013*u2980+u41&
 &14)*u2886)+u4113*u2980+u458)))
  u4113=-2.28164211906493709d-3*u4069
  u3781=-1.02673895357922169d-2*u3971
  u4013=9.22076808974984662d3*u774
  u87=6.84492635719481127d-3*u3811
  u47=-7.99133234444986707d3*u774
  u3796=+u47
  u542=-u754
  u754=u3823+u4018
  u4018=u460*u2885
  u4019=u542*u2980
  value_7_=value_7_+(c(16)*(x*(u2886*(u4026+u2980*(u66+u3796)+u2885*(u4&
 &54+u4013)+(u754+u3837)*u2886)+u2885*(u3781+u4018)+u2980*(u87+u4019)+u&
 &4113)))
  u4113=u4179*(u3833-u3975)
  u3781=5.76255808492236261d-2*u4113
  u87=2.88127904246118131d-2*u639
  u3796=2.93911700706593106d5*u474
  u639=pd2*(u3796+u3982)
  u3833=-2.11293796447153296d-1*u547*u639
  u421=1.40862530964768864d-1*u4065
  value_8_=value_8_+(c(16)*(z*(u2886*(u3781+u2980*(u421*u2980+u3833)+u3&
 &834*u2886)+u2980*(u87+u3782*u2980)+u553)))
  u3781=-2.0373719494149837d-2*u3971
  u87=2.0373719494149837d-2*u3811
  u3833=1.3582479662766558d-2*u999
  u421=-2.43958592624781972d3*u774
  u999=-7.72535543311809578d3*u774
  u3782=9.96048508602880921d-2*u408
  u408=-4.0659765437463662d3*u774
  u3834=-1.19525821032345711d0*u998
  u578=7.47036381452160691d-2*u4065
  u973=-2.84618358062245634d3*u774
  u3827=u2885*(u3824+u421)
  u3825=u578*u2980
  u501=u2980*(u3825+u3834)
  value_9_=value_9_+(c(16)*(x*((u87+u2980*(u797+u408)+u3827)*u2886+u288&
 &5*(u3833+u501+(u3825+u999)*u2885)+u2980*(u3782+u973*u2980)+u3781)))
  u3884=-1.22242316964899022d-1*u3971
  u3965=-u3883
  u3883=u670*u2885
  value_1_=value_1_+(c(17)*(x*(u2885*(u3884+u4020*u2980+u2885*(u3883+u5&
 &22+u3965))+u812*u2980+u823)))
  u3884=-8.64383712738354392d-2*u3971
  u3965=u3822+u4021
  u4020=u524*u2980
  u4021=u403*u2980
  value_2_=value_2_+(c(17)*((u2885*(u3884+u4020+u2885*(u3965+u3885))+u4&
 &021+u857)*z))
  u3884=u472*u2885
  u3885=u895*u2885
  value_3_=value_3_+(c(17)*(x*((u4023+u2885*(u3884+u824))*u2886+u2885*(&
 &u4022+u4050*u2980+u2885*(u3885+u541))+u660*u2980+u853)))
  u416=u3778+u394
  u4022=u4202*u2980
  u4023=u4192*u2980
  value_4_=value_4_+(c(17)*(z*((u461+u2885*(u3779+u3992))*u2886+u2885*(&
 &u899+u4022+u2885*(u416+u3964))+u4023+u854)))
  value_5_=value_5_+(c(17)*(y*(u2886*(u3886+u2885*(u3780+u765)+(u3879+u&
 &4068)*u2886)+u2885*(u629+u2980*(u3821+u3888)+u2885*(u855+u571))+u2980&
 &*(u4130+u805*u2980)+u452)))
  value_6_=value_6_+(c(17)*u780)
  u82=(u446+u715)
  u824=5.70410529766234272d-4*u547
  u571=u824*u82
  u4024=-7.70054215184416268d-2*u3971
  u3886=-5.13369476789610845d-2*u974
  u805=-u3829
  u3876=pd2*(u3987+u943)
  u3874=-6.16043372147533014d-2*u547*u3876
  u765=5.6470642446857193d-2*u494
  u974=-u834
  value_7_=value_7_+(c(17)*(y*((u4024+u2885*(u3823+u805))*u2886+u2885*(&
 &u3886+u2980*(u4011+u765)+u2885*(u627+u3826))+u2980*(u3874+u974*u2980)&
 &+u571)))
  u974=-1.15251161698447252d-1*u3971
  u765=-u3889
  u3886=u468*u2885
  u3826=u667+u3886
  value_8_=value_8_+(c(17)*(u2934*(u2885*(u765+u3826)+u3875+u974)*z))
  value_9_=value_9_+(c(17)*(y*(u2885*(u554+u2980*(u797+u505)+u2885*(u37&
 &73+u3783))+u2980*(u938+u654*u2980)+u426)))
  u974=-u3967
  u765=u522+u3883
  value_1_=value_1_+(c(18)*(u2934*(u2885*(u974+u765)+u3878+u4070)*z))
  u974=(u803-u3988*(1.1d1+u963))
  u505=1.06714038609673382d-3*u547*u974
  u4024=-3.20142115829020145d-3*u4179*(u3819+u3904)
  u3783=-u607
  u3874=u2980*(u614+u833*u2980)
  value_2_=value_2_+(c(18)*(y*((u4064+u2885*(u3822+u569))*u2886+u2885*(&
 &u4024+u3940+(u3995+u3783)*u2885)+u3874+u505)))
  u505=u861+u3885
  u4024=u860*u2980
  value_3_=value_3_+(c(18)*(u934*((u813+u3884)*u2886+u2885*(u3887+u505)&
 &+u4024+u4026)))
  u3783=u2980*(u4195+u3838*u2980)
  value_4_=value_4_+(c(18)*(y*(u2886*(u451+u2885*(u3778+u997)+(u3779+u4&
 &025)*u2886)+u2885*(u4162+u2980*(u968*u2885+u4205))+u3783+u456)))
  u813=-2.52544390488515477d-2*u3971
  u4064=1.44327849774090772d3*u774
  u4070=pd2*(u4117-u3988)
  u805=6.88757428605042211d-3*u547*u4070
  u571=-8.08142049563249528d-1*u998
  u451=5.05088780977030955d-2*u4065
  u4011=-5.15456606336038471d3*u774
  u4025=u451*u2980
  value_5_=value_5_+(c(18)*(z*(u2886*(u966+u2980*(u4025+u571)+u2885*(u3&
 &780+u839)+(u3879+u4025+u965)*u2886)+u2885*(u813+u2980*(u72+u440)+u288&
 &5*(u3788+u3877+u4064))+u2980*(u805+u4011*u2980)+u455)))
  value_6_=value_6_+(c(18)*u585)
  u4064=u824*(u404-u3988*(2.2d1+u3775))
  u451=-4.10695581431688676d-1*u4179*(-u3981+u415)
  u452=-5.6470642446857193d-2*u3971
  u937=9.22076808974984662d2*u774
  u4011=9.03530279149715087d-1*u998
  u459=5.83981979017490286d3*u774
  value_7_=value_7_+(c(18)*(z*(u2886*(u451+u2980*(u3880+u4011)+u412+(u3&
 &880+u860)*u2886)+u2885*(u452+u2980*(u66+u860)+u2885*(u454+u937))+u298&
 &0*(u575+u459*u2980)+u4064)))
  u3880=exppd2*(1.0d0+pd2)
  u575=1.1738544247064072d-2*u547*(u3902+1.13437046457953026d2*u3880)
  u451=-2.88127904246118131d-2*u4179*(-u4003+u871)
  u4011=-3.20142115829020145d-3*u509
  u459=5.63450123859075456d-1*u998
  value_8_=value_8_+(c(18)*(x*((u553+u2885*(u3886+u3786))*u2886+u2885*(&
 &u451+u2980*(u3817+u459)+(u3817+u808)*u2885)+u2980*(u4011+u787*u2980)+&
 &u575)))
  value_9_=value_9_+(c(18)*((u2885*(u812+u2980*(u797+u400)+u2885*(u3773&
 &+u398))+u2980*(u604+u4012))*z))
  u4011=-1.3582479662766558d-2*u509
  u459=-4.52749322092218601d-3*u509
  u575=7.96838806882304737d-1*u998
  value_1_=value_1_+(c(19)*(x*((u604+u2885*(u3883+u844))*u2886+u2885*(u&
 &4011+u2980*(u3890+u575)+(u3890+u809)*u2885)+u2980*(u459+u4188*u2980)+&
 &u812)))
  u4011=-5.76255808492236261d-2*u793
  u604=-5.76255808492236261d-2*u3971
  value_2_=value_2_+(c(19)*(z*(u2886*(u4011+u3940+u3897+u4034*u2886)+u2&
 &885*(u604+u4015)+u3874+u403)))
  value_3_=value_3_+(c(19)*(x*(u2886*(u504+u2885*(u3885+u56)+(u3884+u38&
 &37)*u2886)+u2885*(u804+u827+(u3881+u560)*u2885)+u2980*(u4066+u546*u29&
 &80)+u457)))
  value_4_=value_4_+(c(19)*(z*(u2886*(u393+u613+u2885*(u3778+u411)+(u37&
 &79+u3777+u4114)*u2886)+u2885*(u898+u4017)+u3783+u458)))
  value_5_=value_5_+(c(19)*(y*(u2886*(u877+u2980*(u72+u896)+u2885*(u378&
 &8+u3784)+u840)+u2885*(u4175+u475+(u3882+u409)*u2885)+u2980*(u874+u810&
 &*u2980)+u858)))
  value_6_=value_6_+(c(19)*u3898)
  u967=2.28164211906493709d-3*u870
  u4011=-5.13369476789610845d-2*u3971
  u459=-6.84492635719481127d-3*u3971
  u784=-u47
  u807=1.02673895357922169d-2*u3811
  u575=-u4013
  value_7_=value_7_+(c(19)*(y*(u2886*(u4011+u2980*(u66+u575)+u2885*(u45&
 &4+u784)+(u754+u4050)*u2886)+u2885*(u459+u4018)+u2980*(u807+u4019)+u96&
 &7)))
  value_8_=value_8_+(c(19)*(u934*((u808+u3826)*u2886+u833*u2885+u3875+u&
 &604)))
  value_9_=value_9_+(c(19)*(y*((u87+u2980*(u797+u421)+u2885*(u3824+u408&
 &))*u2886+u2885*(u3782+u501+(u3825+u973)*u2885)+u2980*(u3833+u999*u298&
 &0)+u3781)))
  value_1_=value_1_+(c(20)*(u934*((u765)*u2886+u4121*u2885+u3878)))
  value_2_=value_2_+(c(20)*(y*(u2886*(u4020+u808*u2885+(u3965)*u2886)+u&
 &406*u2885+u4021)))
  u459=1.23208674429506603d-1*u3811
  u406=-3.31947651230994478d4*u774
  value_3_=value_3_+(c(20)*(u934*(u2886*(u406+u505+u472*u2886)+u860*u28&
 &85+u4024+u459)))
  u459=-3.22672922926332217d-3*u4069
  u406=1.30682533785164548d-1*u3811
  u807=-2.08642164608340627d4*u774
  u575=(u2886*(u406+u4022+u4202*u2885+u2886*(u4016+u416+u807))+u4192*u2&
 &885+u4023+u459)
  value_4_=value_4_+(c(20)*(y*u575))
  u4021=5.0d0*pd2
  u4011=-1.27547671963896706d-4*u547*(-u3988*(u4021+1.76d2)+u820)
  u4069=8.16421390851647518d5*u474
  u4022=2.2958580953501407d-3*u4179*(u4069+3.45982991696756728d3*exppd2&
 &)
  u670=-u3784
  u859=-5.8544381431428588d-2*u3971
  u967=1.52575155475467387d4*u774
  u784=-6.18547927603246165d2*u774
  u3834=+u784
  u3874=4.82107447446300359d2*exppd2
  u422=-1.37751485721008442d-2*u4179*(u3874+u62)
  u646=1.26272195244257739d-2*u4179*(1.73081334860549274d6*u474+u3873)
  u472=-6.31360976221288694d-2*u4065
  u604=-u784
  u784=-6.31360976221288694d-3*u4065
  value_5_=value_5_+(c(20)*(z*(u2886*(u4022+u2980*(u784*u2980+u646)+u28&
 &85*(u3788+u967)+u2886*(u465+u472*u2980+u670))+u2885*(u859+u3834*u2885&
 &)+u2980*(u422+u604*u2980)+u4011)))
  value_6_=value_6_+(c(20)*(x*u575))
  u784=-6.16043372147533014d-2*u3971
  u967=1.65973825615497239d4*u774
  u4011=-u937
  u670=6.16043372147533014d-2*u3811
  u859=-u967
  value_7_=value_7_+(c(20)*(z*(u2886*(u2980*(u66+u859)+u2885*(u454+u967&
 &)+(u754)*u2886)+u2885*(u784+u4011*u2885)+u2980*(u670+u937*u2980))))
  value_8_=value_8_+(c(20)*(x*(u2886*(u4008+u3786*u2885+(u3886+u667)*u2&
 &886)+u553*u2885+u873*u2980)))
  u784=4.07474389882996741d-2*u4113
  u4011=8.14948779765993481d-2*u4039
  u670=-7.47036381452160691d-2*u494
  u859=-8.53855074186736902d3*u774
  u967=8.71542445027520806d-2*u4065
  value_9_=value_9_+(c(20)*(z*(u2886*(u784+u2980*(u967*u2980+u670)+u382&
 &7+(u3825+u421)*u2886)+u2885*(u87+u654*u2885)+u2980*(u4011+u859*u2980)&
 &+u3781)))
  if ( lmax .eq. 3 ) go to 100
  u784=p**9
  u4011=u954*u784
  u670=u4011*u392
  u859=-2.71649593255331161d-1*u670
  u967=u4011*u4118
  u3834=-2.43958592624781972d4*u967
  u422=+u3834
  u646=1.21979296312390986d5*u967
  u472=u4011*u3939
  u604=4.98024254301440461d-1*u472
  u3963=-8.9644365774259283d-1*u472
  u468=6.36808684864285064d6*u474
  u61=u4011*pd2
  u3886=(1.5d1+u3774)
  u4022=u3989*u3886
  u4167=u61*(u468-2.26874092915906051d2*u4022)
  u402=-4.98024254301440461d-2*u4167
  u806=4.98024254301440461d-2*u4167
  u3905=u806*u2980
  u3883=u402*u2980
  value_1_=value_1_+(c(21)*(u2934*((u422+u2980*(u3883+u604))*u2885+u298&
 &0*(u646+u2980*(u3905+u3963))+u859)))
  u859=-1.15251161698447252d-1*u670
  u3963=-3.45009550347419585d3*u967
  u776=+u3963
  u452=1.34553724635493638d5*u967
  u4192=2.11293796447153296d-1*u472
  u578=-1.47905657513007307d0*u472
  u553=+u578
  u406=-3.5215632741192216d-2*u4167
  u807=1.05646898223576648d-1*u4167
  u3781=u406*u2980
  u575=(u3781+u4192)
  u877=u807*u2980
  value_2_=value_2_+(c(21)*(y*((u776+u2980*u575)*u2885+u2980*(u452+u298&
 &0*(u877+u553))+u859)*z))
  u553=u4011*u527
  u874=1.02673895357922169d-1*u553
  u4175=5.53246085384990798d4*u967
  u828=-9.22076808974984662d3*u967
  u451=-4.61038404487492331d4*u967
  u4064=-1.12941284893714386d0*u472
  u4065=1.8823547482285731d-1*u472
  u998=3.38823854681143158d-1*u472
  u812=1.12941284893714386d-1*u4167
  u403=-1.8823547482285731d-2*u4167
  u937=u403*u2980
  u3871=u812*u2980
  value_3_=value_3_+(c(21)*(u2934*((u4175+u2980*(u3871+u4064))*u2886+(u&
 &828+u2980*(u937+u4065))*u2885+u2980*(u451+u2980*(u937+u998))+u874)))
  u668=4.35608445950548493d-2*u553
  u494=5.21605411520851568d3*u967
  u873=-3.91204058640638676d3*u967
  u898=+u873
  u87=-5.08565276232830279d4*u967
  u4179=-3.19446193697068895d-1*u472
  u504=2.39584645272801671d-1*u472
  u4066=5.59030838969870566d-1*u472
  u804=5.32410322828448158d-2*u4167
  u393=-3.99307742121336118d-2*u4167
  u3833=u393*u2980
  u3782=u804*u2980
  u870=(u2936*((u494+u2980*(u3782+u4179))*u2886+(u898+u2980*(u3833+u504&
 &))*u2885+u2980*(u87+u2980*(u3833+u4066))+u668))
  value_4_=value_4_+(c(21)*u870)
  u3971=3.06114412713352094d-3*u4011*(u820+u3987*(1.1d1+u962))
  u3811=4.13254457163025327d-2*u553
  u509=1.64946114027532311d3*u967
  u4039=u4011*u499
  u793=-1.37751485721008442d-2*u4039
  u4113=5.56693134842921549d3*u967
  u4112=u61*(-u3861+u3991)
  u54=-2.75502971442016884d-2*u4112
  u614=-5.44322176290856626d4*u967
  u4068=-1.01017756195406191d-1*u472
  u464=3.78816585732773216d-2*u61
  u4008=-1.81499274332724841d3*u3989
  u4067=u464*(4.21273437679450119d6*u474+u4008)
  u628=-3.40934927159495895d-1*u472
  u421=6.31360976221288694d-3*u61
  u56=u421*(8.98063529936812269d6*u474+u4008)
  u654=6.56615415270140241d-1*u472
  u813=1.68362926992343652d-2*u4167
  u3992=4.0d0*pd2
  u4013=u3989*(5.1d1+u3992)
  u460=-2.52544390488515477d-2*u61*(2.16514952853856922d7*u474-2.268740&
 &92915906051d2*u4013)
  u818=5.68224878599159824d-2*u4167
  u3837=-1.89408292866386608d-1*u472
  u411=-5.05088780977030955d-2*u4167
  u824=6.31360976221288694d-2*u4167
  u805=6.31360976221288694d-3*u4167
  u571=u805*u2980
  u607=u818*u2980
  u3783=u411*u2980
  u3826=u824*u2980
  u3887=u813*u2980
  value_5_=value_5_+(c(21)*(u2886*(u3811+u2980*(u2980*(u654+u3783)+u614&
 &)+(u2980*(u4068+u3887)+u509)*u2886)+u2885*(u793+u2980*(u2980*(u460+u3&
 &826)+u4067)+(u2980*(u628+u607)+u4113)*u2885)+u2980*(u54+u2980*(u2980*&
 &(u3837+u571)+u56))+u3971))
  u965=2.17804222975274246d-1*u553
  u440=2.60802705760425784d4*u967
  u398=-1.95602029320319338d4*u967
  u546=-9.78010146601596691d4*u967
  u560=-5.32410322828448158d-1*u472
  u4114=3.99307742121336118d-1*u472
  u999=7.18753935818405013d-1*u472
  value_6_=value_6_+(c(21)*(u2937*((u440+u2980*(u3782+u560))*u2886+(u39&
 &8+u2980*(u3833+u4114))*u2885+u2980*(u546+u2980*(u3833+u999))+u965)))
  u3875=5.70410529766234272d-4*u4011*u82
  u973=-7.70054215184416268d-2*u670
  u887=u4011*u3876
  u4202=-6.16043372147533014d-2*u887
  u808=6.45453766282489263d3*u967
  u82=u61*(2.67786216199340386d5*u474-u453)
  u3876=-2.56684738394805423d-2*u82
  u416=8.29869128077486196d4*u967
  u4188=9.47048813387911121d5*u474
  u765=u61*(u4188+u3862)
  u505=1.69411927340571579d-1*u765
  u997=-3.95294497128000351d-1*u472
  u400=-2.82353212234285965d-2*u61*(-u3862+u943)
  u809=-8.47059636702857895d-1*u472
  u95=u61*(4.6699303556714238d6*u474+u3873*(1.1d1+pd2))
  u569=-1.12941284893714386d-1*u95
  u839=6.58824161880000585d-2*u4167
  u896=1.50588379858285848d-1*u472
  u810=5.6470642446857193d-2*u4167
  u409=-9.41177374114286549d-3*u4167
  u3884=u409*u2980
  u774=u810*u2980
  u3888=u839*u2980
  value_7_=value_7_+(c(21)*((u973+u2980*(u2980*(u809+u774)+u416))*u2886&
 &+u2885*(u4202+u2980*(u2980*(u569+u774)+u505)+(u2980*(u997+u3888)+u808&
 &)*u2885)+u2980*(u3876+u2980*(u2980*(u896+u3884)+u400))+u3875))
  u3875=-1.92085269497412087d-1*u670
  u973=-5.17514325521129378d4*u967
  u4202=+u973
  u3876=8.62523875868548964d4*u967
  u505=1.05646898223576648d0*u472
  u997=-6.33881389341459887d-1*u472
  u400=-1.05646898223576648d-1*u4167
  u809=3.5215632741192216d-2*u4167
  u569=u400*u2980
  u896=u809*u2980
  value_8_=value_8_+(c(21)*(x*((u4202+u2980*(u569+u505))*u2885+u2980*(u&
 &3876+u2980*(u896+u997))+u3875)*z))
  u3875=3.01832881394812401d-3*u4011*u523
  u523=8.14948779765993481d-2*u553
  u3939=1.21979296312390986d3*u967
  u392=-1.90154715278731812d-1*u670
  u527=-9.51438511236649692d4*u967
  u499=-7.47036381452160691d-2*u472
  u4130=5.00115114880803043d4*u967
  u3778=1.04585093403302497d0*u472
  u803=1.24506063575360115d-2*u4167
  u871=-2.73913339865792253d-1*u472
  u415=-7.47036381452160691d-2*u4167
  u3904=u803*u2980
  u454=u415*u2980
  value_9_=value_9_+(c(21)*(u2885*(u523+u2980*(u2980*(u3778+u454)+u527)&
 &+(u2980*(u499+u3904)+u3939)*u2885)+u2980*(u392+u2980*(u2980*(u871+u39&
 &04)+u4130))+u3875))
  u900=1.3582479662766558d-2*u553
  u547=-4.87917185249563944d3*u967
  u3995=+u547
  u715=-1.22242316964899022d-1*u670
  u3838=1.46375155574869183d4*u967
  u938=2.98814552580864276d-1*u472
  u3793=1.95166874099825578d4*u967
  u4061=-4.48221828871296415d-1*u472
  u524=-4.98024254301440461d-2*u472
  u4026=u524*u2980
  value_1_=value_1_+(c(22)*(u2885*(u900+u2980*(u2980*(u4061+u3905)+u383&
 &8)+(u2980*(u938+u3883)+u3995)*u2885)+u2980*(u715+u2980*(u4026+u3793))&
 &+u3875))
  u900=3.10508595312677627d4*u967
  u715=1.05646898223576648d-1*u472
  u4061=-8.45175185788613183d-1*u472
  u844=-1.05646898223576648d-1*u472
  u897=(u3781+u715)
  u3889=u844*u2980
  value_2_=value_2_+(c(22)*(x*(u2885*(u900+u2980*(u877+u4061)+u897*u288&
 &5)+u2980*(u900+u3889)+u859)*z))
  u859=-4.56328423812987418d-3*u4011*u419
  u419=5.6470642446857193d-2*u61
  u4050=u419*(u62-u3988)
  u4195=-1.29090753256497853d4*u967
  u834=7.70054215184416268d-2*u61*(3.20037185213845827d5*u474-u3988)
  u3829=-1.69411927340571579d-1*u61
  u3786=+u3829*(2.05738190494615174d6*u474+u3872)
  u4121=7.90588994256000702d-1*u472
  u787=6.85793968315383915d5*u474
  u833=u61*(u787+u3873)
  u3967=-1.12941284893714386d-1*u833
  u3825=u3989*(4.9d1+u3992)
  u860=u419*(2.08024170388999788d7*u474-2.26874092915906051d2*u3825)
  u419=-1.31764832376000117d-1*u4167
  u3965=1.31764832376000117d-1*u472
  u3784=u419*u2980
  u3773=u2980*(u860+u3784)
  value_3_=value_3_+(c(22)*(u2885*(u4050+u2980*(u3773+u3786)+(u2980*(u4&
 &121+u3784)+u4195)*u2885)+u2980*(u834+u2980*(u3965*u2980+u3967))+u859)&
 &)
  u47=u4011*u423
  u423=1.74243378380219397d-1*u47
  u554=u61*(u3991+u3873)
  u3877=-1.59723096848534447d-1*u554
  u895=1.59723096848534447d-1*u472
  u585=-4.30324464504702544d4*u967
  u944=1.19792322636400836d-1*u472
  u780=u61*(2.38395046128681075d6*u474+u3872)
  u3898=-3.99307742121336118d-2*u780
  u840=2.12269561621428355d6*u474
  u4017=u3989*(1.0d1+pd2)
  u613=u61*(u840-1.13437046457953026d2*u4017)
  u3940=2.12964129131379263d-1*u613
  u412=-5.32410322828448158d-2*u4167
  u3897=6.38892387394137789d-1*u472
  u827=1.99653871060668059d-1*u472
  u475=u2980*(u3833+u3897)
  u3827=u412*u2980
  u501=(u2937*(u2886*(u3877+u2980*(u3827+u3940)+(u3827+u895)*u2886)+u28&
 &85*(u585+u475+(u3833+u944)*u2885)+u2980*(u3898+u827*u2980)+u423))
  value_4_=value_4_+(c(22)*u501)
  u855=-1.65301782865210131d-1*u887
  u629=-1.4845150262477908d4*u967
  u72=-5.05088780977030955d-2*u472
  u754=u61*(u91+u3862)
  u3964=7.57633171465546432d-2*u754
  u3817=-1.70467463579747947d-1*u472
  u465=1.27361736972857013d6*u474
  u656=7.57633171465546432d-2*u61*(u465+u3862)
  u660=4.04071024781624764d-1*u472
  u4016=8.0d0*pd2
  u899=u3989*(9.5d1+u4016)
  u461=-1.26272195244257739d-2*u61*(4.03312167080713874d7*u474-2.268740&
 &92915906051d2*u899)
  u407=-2.71485219775154138d-1*u472
  value_5_=value_5_+(c(22)*(u2934*(u2886*(u629+u2980*(u3783+u660)+(u388&
 &7+u72)*u2886)+u2885*(u3964+u2980*(u3826+u461)+(u607+u3817)*u2885)+u29&
 &80*(u656+u2980*(u571+u407))+u855)))
  value_6_=value_6_+(c(22)*u870)
  u870=u61*(3.10240128523626057d6*u474+u3861)
  u4034=-5.13369476789610845d-3*u870
  u4014=2.76623042692495399d4*u967
  u4205=u61*(2.72140463617215839d5*u474+u3982)
  u408=3.38823854681143158d-1*u4205
  u404=-1.97647248564000175d-1*u472
  u542=u4011*u875
  u875=1.35529541872457263d0*u542
  u797=-5.6470642446857193d-1*u472
  u3960=6.0d0*pd2
  u627=u3989*(6.5d1+u3960)
  u395=u61*(2.75950430107856861d7*u474-2.26874092915906051d2*u627)
  u3936=-1.8823547482285731d-2*u395
  u966=-2.82353212234285965d-2*u472
  u530=u2980*(u3884+u966)
  value_7_=value_7_+(c(22)*(u2934*((u4014+u2980*(u774+u797))*u2886+u288&
 &5*(u408+u2980*(u774+u3936)+(u3888+u404)*u2885)+u2980*(u875+u530)+u403&
 &4)))
  u4034=3.84170538994824174d-2*u553
  u408=-1.03502865104225876d4*u967
  u404=+u408
  u797=6.33881389341459887d-1*u472
  u3936=-2.11293796447153296d-1*u472
  u4162=u2980*(u896+u3936)
  value_8_=value_8_+(c(22)*(y*((u404+u2980*(u569+u797))*u2885+u2980*(u4&
 &04+u4162)+u4034)*z))
  u3899=1.3582479662766558d-1*u553
  u522=-3.73518190726080346d-2*u472
  u394=6.22530317876800577d-1*u472
  u667=(u3904+u522)
  u3892=u667*u2885
  value_9_=value_9_+(c(22)*(u2934*(u2885*(u422+u2980*(u454+u394)+u3892)&
 &+u2980*(u422+u2980*u667)+u3899)))
  u667=-2.7164959325533116d-2*u670
  u541=4.3912546672460755d4*u967
  u994=-5.97629105161728553d-1*u472
  value_1_=value_1_+(c(23)*(y*((u3995+u2980*(u3883+u938))*u2885+u2980*(&
 &u541+u2980*(u3905+u994))+u667)*z))
  u994=-2.88127904246118131d-2*u4112
  u4112=-u973
  u973=1.05646898223576648d-1*u554
  u66=u61*(u4069+u3873)
  u4069=1.05646898223576648d-1*u66
  u968=-u505
  u861=-1.40862530964768864d-1*u613
  u3890=u3936*u2980
  value_2_=value_2_+(c(23)*(u2934*((u4112+u2980*(u877+u968))*u2886+u288&
 &5*(u973+u2980*(u896+u861)+(u896+u844)*u2885)+u2980*(u4069+u3890)+u994&
 &)))
  u994=-5.13369476789610845d-2*u670
  u4069=1.10649217076998159d4*u967
  u3788=-1.84415361794996932d3*u967
  u471=-6.77647709362286316d-1*u472
  u4115=1.12941284893714386d-1*u472
  value_3_=value_3_+(c(23)*(u2936*((u4069+u2980*(u3871+u471))*u2886+(u3&
 &788+u2980*(u937+u4115))*u2885+u2980*(u4014+u403*u4191)+u994)))
  u4006=-1.08902111487637123d-2*u61
  u3985=+u4006*(u4188+u3861)
  u4188=1.17361217592191603d4*u967
  u65=-1.59723096848534447d-1*u472
  u439=1.19792322636400836d-1*u554
  u4002=-1.19792322636400836d-1*u472
  u383=7.98615484242672237d-2*u472
  u815=-1.59723096848534447d-1*u613
  u802=3.99307742121336118d-2*u4167
  u531=-7.98615484242672237d-2*u472
  u631=u802*u2980
  u891=u2980*(u631+u815)
  u572=(u2934*(u2886*(u4188+u2980*(u3833+u383)+(u3782+u65)*u2886)+u2885&
 &*(u439+u891+(u631+u4002)*u2885)+u2980*(u944+u531*u2980)+u3985))
  value_4_=value_4_+(c(23)*u572)
  u586=-1.92852080009411819d-1*u61*(6.06484461775509585d4*u474+u3975)
  u909=1.21221307434487429d0*u542
  u658=-2.02035512390812382d-1*u472
  u4186=9.89676684165193864d3*u967
  u947=-1.89408292866386608d-2*u472
  u842=u61*(u943+u717)
  u663=6.06106537172437146d-1*u842
  u4024=3.0d0*pd2
  u4020=u3989*(2.0d1+u4024)
  u462=-1.34690341593874921d-1*u61*(u840-5.67185232289765129d1*u4020)
  u840=6.73451707969374607d-2*u4167
  u92=-1.6415385381753506d-1*u472
  u872=-1.452130245308964d-1*u472
  u825=1.26272195244257739d-2*u4167
  u3891=u825*u2980
  u4027=u840*u2980
  value_5_=value_5_+(c(23)*(u2937*(u2886*(u909+u462*u2980+(u4027+u658)*&
 &u2886)+u2885*(u4186+u2980*(u3891+u92)+(u571+u947)*u2885)+u2980*(u663+&
 &u2980*(u571+u872))+u586)))
  u767=1.21002346097374581d-3*u4011*(u41+u985)
  u41=-3.2670633446291137d-2*u670
  u985=-u873
  u873=-2.17804222975274246d-2*u4011*u498
  u498=u61*(8.81735102119779319d5*u474+u3862)
  u538=1.19792322636400836d-1*u498
  u4134=-2.39584645272801671d-1*u472
  u3847=3.99307742121336118d-2*u61
  u680=u3847*(u62+u3873)
  u3937=2.79515419484935283d-1*u472
  u459=u3989*(2.1d1+u3774)
  u959=u61*(8.9153215880999909d6*u474-2.26874092915906051d2*u459)
  u444=-3.99307742121336118d-2*u959
  u3961=u444+u631
  u845=u2980*(u3961)
  value_6_=value_6_+(c(23)*(u2886*(u41+u2980*(u2980*(u3937+u3833)+u985)&
 &+(u2980*(u4179+u3782)+u494)*u2886)+u2885*(u898+u2980*(u845+u538)+(u29&
 &80*(u4134+u631)+u985)*u2885)+u2980*(u873+u680*u2980)+u767))
  u876=7.77233164090768437d5*u474
  u570=-2.56684738394805423d-2*u61
  u835=+u570*(-u3861+u876)
  u3980=u61*(u794+u717)
  u3779=6.77647709362286316d-1*u3980
  u747=-1.69411927340571579d-1*u472
  u4085=5.53246085384990797d3*u967
  u576=u61*(3.15682937795970374d5*u474+u3982)
  u852=3.38823854681143158d-1*u576
  u4018=u3989*(1.5d1+pd2)
  u941=u61*(u468-2.26874092915906051d2*u4018)
  u785=-1.12941284893714386d-1*u941
  u753=-5.6470642446857193d-2*u472
  u811=9.41177374114286549d-3*u4167
  u789=-1.41176606117142983d-1*u472
  u3828=u811*u2980
  u4028=u753*u2980
  value_7_=value_7_+(c(23)*(u2937*(u2886*(u3779+u2980*(u3871+u785)+(u77&
 &4+u747)*u2886)+u2885*(u4085+u4028+(u3828+u966)*u2885)+u2980*(u852+u29&
 &80*(u3884+u789))+u835)))
  u835=1.06714038609673382d-3*u4011*u974
  u3779=-4.80213173743530218d-2*u670
  u852=-u408
  u785=-5.76255808492236261d-2*u4039
  u789=3.16940694670729944d-1*u498
  u974=u61*(u4117+u3873)
  u408=1.05646898223576648d-1*u974
  u3836=-5.28234491117883239d-1*u472
  u51=-1.05646898223576648d-1*u959
  u711=-1.40862530964768864d-1*u472
  value_8_=value_8_+(c(23)*((u3779+u2980*(u2980*(u3836+u896)+u4112))*u2&
 &886+u2885*(u404+u2980*(u2980*(u51+u877)+u789)+(u2980*(u997+u877)+u852&
 &)*u2885)+u2980*(u785+u2980*(u711*u2980+u408))+u835))
  u835=-2.92750311149738367d4*u967
  u3779=+u835
  u785=6.72332743306944622d-1*u472
  u789=-1.86759095363040173d-1*u472
  value_9_=value_9_+(c(23)*(x*(u2885*(u3779+u2980*(u454+u785)+u3892)+u2&
 &980*(u3793+u2980*(u3904+u789))+u667)*z))
  u408=1.49407276290432138d-1*u472
  u3836=-u3834
  u51=-1.49407276290432138d-1*u472
  u711=(u3883+u408)
  u3834=u711*u2885
  u3892=u51*u2980
  value_1_=value_1_+(c(24)*(u2934*(u2885*(u422+u806*u4191+u3834)+u2980*&
 &(u3836+u3892))))
  u4049=-3.84170538994824174d-2*u670
  u4201=-1.15003183449139862d3*u967
  u671=+u4201
  u3870=3.5215632741192216d-2*u472
  u4059=4.48512415451645461d4*u967
  u4031=-2.81725061929537728d-1*u472
  u4052=-3.16940694670729944d-1*u472
  u4029=u4052*u2980
  value_2_=value_2_+(c(24)*(y*(u2885*(u671+u2980*(u877+u4031)+(u3781+u3&
 &870)*u2885)+u2980*(u4059+u4029)+u4049)*z))
  u671=u4011*u642
  u4031=3.08021686073766507d-2*u671
  u642=-1.69411927340571579d-1*u754
  u637=3.95294497128000351d-1*u472
  u4062=1.48588693134999848d7*u474
  u3878=u3989*(3.5d1+u4024)
  u766=u61*(u4062-2.26874092915906051d2*u3878)
  u4032=7.52941899291429239d-2*u766
  u4030=u637*u2980
  value_3_=value_3_+(c(24)*(u2934*(u2885*(u642+u2980*(u3784+u4032)+(u37&
 &84+u637)*u2885)+u2980*(u642+u4030)+u4031)))
  u4031=u412*u2885
  u4032=u944*u2980
  u642=(u2936*(u2886*(u3877+u2885*(u4031+u3940)+(u4031+u895)*u2886)+u28&
 &85*(u3898+u475+(u3833+u827)*u2885)+u2980*(u585+u4032)+u423))
  value_4_=value_4_+(c(24)*u642)
  u475=1.02038137571117365d-3*u4011*(u446+u3987*(1.1d1+9.2d1*pd2))
  u446=u61*(u787+u3988)
  u787=4.59171619070028141d-3*u446
  u566=u4011*u639
  u639=-3.36725853984687303d-2*u566
  u850=1.68362926992343652d-2*u472
  u500=u421*(6.17214571483845523d6*u474+u4008)
  u3926=-5.68224878599159824d-2*u472
  u399=-4.13254457163025327d-2*u61*(-u3988+3.51968332944932485d5*u474)
  u3893=-5.05088780977030955d-2*u498
  u577=1.68362926992343652d-2*u959
  u463=-1.68362926992343652d-2*u4167
  u4023=-3.62998548665449682d3*u3989
  u3772=u464*(7.93561591907801387d6*u474+u4023)
  u4031=1.6d1*pd2
  u464=-6.31360976221288694d-3*u61*(1.0316300694801418d8*u474+u3873*(2.&
 &43d2+u4031))
  u862=u421*(5.91089086976592803d6*u474+u4008)
  u430=1.01017756195406191d-1*u472
  u445=1.13644975719831965d-1*u4167
  u4033=u463*u2980
  value_5_=value_5_+(c(24)*(u2886*(u787+u2980*(u430*u2980+u3893)+u2886*&
 &((u850+u4033)*u2886+u2980*(u577+u4033)+u639))+u2885*(u855+u2980*(u298&
 &0*(u464+u607)+u3772)+u2885*((u3926+u607)*u2885+u2980*(u464+u445*u2980&
 &)+u500))+u2980*(u399+u2980*(u3926*u2980+u862))+u475))
  value_6_=value_6_+(c(24)*u501)
  u475=u4011*u961
  u787=-5.13369476789610845d-3*u475
  u639=2.82353212234285965d-2*u754
  u850=-6.58824161880000585d-2*u472
  u500=5.13369476789610845d-3*u671
  u4033=u3989*(3.5d1+u3992)
  u430=u61*(u4062-2.26874092915906051d2*u4033)
  u577=-2.82353212234285965d-2*u430
  u3772=-2.82353212234285965d-2*u754
  u862=2.82353212234285965d-2*u430
  u430=6.58824161880000585d-2*u472
  u445=-6.58824161880000585d-2*u4167
  u3893=u445*u2980
  value_7_=value_7_+(c(24)*(u2885*(u787+u4191*(u3893+u862)+u2885*((u850&
 &+u3888)*u2885+u577*u2980+u639))+u2980*(u500+u2980*(u430*u2980+u3772))&
 &))
  u862=-u4059
  u639=3.16940694670729944d-1*u472
  u430=-u4201
  u787=2.81725061929537728d-1*u472
  u577=-3.5215632741192216d-2*u472
  u3772=(u569+u639)
  value_8_=value_8_+(c(24)*(x*(u2885*(u862+u2980*(u896+u787)+u3772*u288&
 &5)+u2980*(u430+u577*u2980)+u4034)*z))
  u862=-3.01832881394812401d-3*u4011*u3850
  u430=5.43299186510662321d-2*u553
  u787=-u3939
  u850=-1.24506063575360115d-2*u472
  u500=-5.12313044512042141d4*u967
  u3926=2.61462733508256242d-1*u472
  u4034=u850*u2980
  value_9_=value_9_+(c(24)*(u2885*(u430+u2980*(u2980*(u3926+u3904)+u500&
 &)+u2885*((u850+u3904)*u2885+u2980*(u3926+u454)+u787))+u2980*(u430+u29&
 &80*(u4034+u787))+u862))
  u500=1.13847343224898254d4*u967
  u3926=-1.99209701720576184d-1*u472
  value_1_=value_1_+(c(25)*(x*(u2885*(u3995+u2980*(u3905+u3926)+u3834)+&
 &u2980*(u500+u4026)+u667)*z))
  u430=-u3981*(u4031-1.1d1)
  u787=u4011*(u430+u3902)
  u3850=-1.70742461775477411d-2*u787
  u961=3.84170538994824174d-2*u47
  u399=u61*(u3796+u3873)
  u458=3.5215632741192216d-2*u399
  u4201=1.41513041080952236d5*u474
  u4059=u61*(u4201-u3987)
  u501=1.15251161698447252d-1*u4059
  u3834=-2.11293796447153296d-1*u498
  u3796=u61*(u465+u3873*(3.0d0+u3774))
  u79=-3.5215632741192216d-2*u3796
  u63=u4011*u418
  u418=-1.05646898223576648d-1*u63
  u3857=u3989*(1.9d1+u3774)
  u645=u61*(8.06624334161427748d6*u474-2.26874092915906051d2*u3857)
  u971=1.05646898223576648d-1*u645
  u4105=-7.04312654823844319d-2*u4167
  u4035=u715*u2980
  value_2_=value_2_+(c(25)*(u2885*(u961+u2980*(u2980*(u971+u569)+u3834)&
 &+u2885*((u577+u896)*u2885+u2980*(u79+u4105*u2980)+u458))+u2980*(u501+&
 &u2980*(u4035+u418))+u3850))
  u3850=1.02673895357922169d-1*u61*(3.06974442960219467d5*u474+u3975)
  u961=-3.38823854681143158d-1*u554
  u458=-2.02856897974496626d4*u967
  u501=5.6470642446857193d-2*u472
  u3834=+u3829*(1.02687668271562777d6*u474+u3862)
  u79=4.51765139574857544d-1*u613
  u418=-1.12941284893714386d-1*u4167
  u971=3.01176759716571696d-1*u472
  u4105=3.57647402163428889d-1*u472
  u577=u2980*(u937+u971)
  u3829=u418*u2980
  value_3_=value_3_+(c(25)*(u2937*(u2886*(u961+u2980*(u3829+u79)+(u3829&
 &+u998)*u2886)+u2885*(u458+u577+(u937+u501)*u2885)+u2980*(u3834+u4105*&
 &u2980)+u3850)))
  u993=u4011*(u3902+u430)
  u3902=6.45345845852664434d-3*u993
  u430=1.45202815316849498d-2*u446
  u414=-1.06482064565689632d-1*u566
  u4054=5.32410322828448158d-2*u472
  u915=u61*(u3903-3.54490770181103205d0*exppd2)
  u3903=-3.48486756760438794d-1*u915
  u420=3.99307742121336118d-2*u833
  u4099=-3.99307742121336118d-2*u472
  u50=-4.35608445950548493d-2*u4059
  u908=-1.59723096848534447d-1*u498
  u863=5.32410322828448158d-2*u959
  u643=2.39584645272801671d-1*u498
  u77=u3989*(2.7d1+u3774)
  u723=u61*(1.14625563275571312d7*u474-2.26874092915906051d2*u77)
  u428=-3.99307742121336118d-2*u723
  u4055=3.99307742121336118d-2*u63
  u4056=3.19446193697068895d-1*u472
  u864=7.98615484242672237d-2*u4167
  u3879=(u2886*(u430+u2980*(u4056*u2980+u908)+u2886*((u4054+u3827)*u288&
 &6+u2980*(u863+u3827)+u414))+u2885*(u3903+u2980*(u2980*(u428+u631)+u64&
 &3)+u2885*((u4099+u631)*u2885+u2980*(u428+u864*u2980)+u420))+u2980*(u5&
 &0+u2980*(u4099*u2980+u4055))+u3902)
  value_4_=value_4_+(c(25)*u3879)
  u397=-1.28568053339607879d-1*u61*(1.07301097083359388d5*u474-u3987)
  u539=1.21221307434487429d0*u61*(6.16851717532355902d4*u474+u3863)
  u4109=-1.68362926992343652d-1*u472
  u801=1.51526634293109286d-1*u61*(5.11624071600365778d5*u474+u3873)
  u841=-2.02035512390812382d-1*u613
  u819=5.05088780977030955d-2*u4167
  u4097=-1.57840244055322173d-1*u472
  u4025=-u509
  u4161=1.34690341593874921d-1*u472
  u917=-1.38899414768683513d-1*u472
  u4019=1.89408292866386608d-2*u472
  u4036=u819*u2885
  value_5_=value_5_+(c(25)*(u2936*(u2886*(u539+u2980*(u3783+u4161)+u288&
 &5*(u4036+u841)+(u4036+u3887+u4109)*u2886)+u2885*(u801+u2980*(u3891+u9&
 &17)+(u571+u4097)*u2885)+u2980*(u4025+u2980*(u571+u4019))+u397)))
  value_6_=value_6_+(c(25)*u572)
  u572=-1.65973825615497239d4*u967
  u3962=+u572
  u615=u4011*u540
  u540=2.71059083744914526d0*u615
  u679=u61*(1.26998883021367392d5*u474+u717)
  u3966=6.77647709362286316d-1*u679
  u4057=-2.25882569787428772d-1*u613
  u4051=-1.78823701081714444d-1*u472
  u3925=2.39739970333496012d4*u967
  u520=-3.38823854681143158d-1*u472
  u3791=-9.4117737411428655d-2*u472
  u4037=u810*u2885
  value_7_=value_7_+(c(25)*(u2936*(u2886*(u540+u2980*(u774+u520)+u2885*&
 &(u4037+u4057)+(u4037+u747)*u2886)+u2885*(u3966+u3791*u2980+(u3828+u40&
 &51)*u2885)+u2980*(u3925+u530)+u3962)))
  u3962=-9.60426347487060435d-3*u870
  u3966=1.72504775173709793d4*u967
  u4057=3.16940694670729944d-1*u554
  u4051=u61*(5.26138229659950623d5*u474+u3873)
  u3791=3.16940694670729944d-1*u4051
  u530=-3.5215632741192216d-1*u472
  u865=-4.22587592894306592d-1*u613
  u40=(u877+u4052)
  u4038=u530*u2980
  value_8_=value_8_+(c(25)*(u2934*((u3966+u2980*(u896+u530))*u2886+u288&
 &5*(u4057+u2980*(u877+u865)+u40*u2885)+u2980*(u3791+u4038)+u3962)))
  u3962=2.7164959325533116d-2*u553
  u3791=-3.25278123499709296d3*u967
  u904=-u3793
  u4058=3.237157652959363d-1*u472
  u758=3.73518190726080346d-2*u472
  value_9_=value_9_+(c(25)*(y*(u2885*(u3791+u2980*(u454+u4058)+(u3904+u&
 &850)*u2885)+u2980*(u904+u2980*(u3904+u758))+u3962)*z))
  u3835=-1.3582479662766558d-2*u475
  u3957=1.49407276290432138d-1*u554
  u3996=u61*(5.98709019957874846d5*u474+u3873)
  u3795=1.49407276290432138d-1*u3996
  u4063=-4.98024254301440461d-1*u472
  u866=-1.99209701720576184d-1*u613
  value_1_=value_1_+(c(26)*(u2934*((u3836+u2980*(u3905+u4063))*u2886+u2&
 &885*(u3957+u2980*(u3905+u866)+(u3905+u51)*u2885)+u2980*(u3795+u3926*u&
 &2980)+u3835)))
  u3835=-2.88127904246118131d-2*u4039
  u3795=2.11293796447153296d-1*u566
  u3926=6.21017190625355254d4*u967
  u4039=u809*u2885
  value_2_=value_2_+(c(26)*(u2936*(u2886*(u3795+u2980*(u877+u997)+u2885&
 &*(u4039+u861)+(u4039+u844)*u2886)+u2885*(u973+u844*u2885)+u2980*(u392&
 &6+u4029)+u3835)))
  u3795=+u570*(-u4003+2.80848958452966746d5*u474)
  u973=5.7168762156449049d4*u967
  u570=5.6470642446857193d-2*u554
  u4119=5.6470642446857193d-2*u4011*u826
  u3827=-4.89412234539429006d-1*u472
  u826=-7.52941899291429239d-2*u613
  u814=1.8823547482285731d-2*u4167
  u447=-3.7647094964571462d-2*u472
  u3816=u814*u2980
  u431=u2980*(u3816+u826)
  value_3_=value_3_+(c(26)*(u2934*(u2886*(u973+u2980*(u937+u3827)+(u387&
 &1+u520)*u2886)+u2885*(u570+u431+(u3816+u753)*u2885)+u2980*(u4119+u447&
 &*u2980)+u3795)))
  u742=-3.63007038292123744d-3*u61*(-u4003+3.29834241904065597d6*u474)
  u568=2.39584645272801671d-1*u61*(3.01168779736385529d5*u474+u3982)
  u949=-1.73033354919245651d-1*u472
  u880=2.60802705760425784d3*u967
  u52=-2.92825677555646487d-1*u472
  u3894=u802*u2885
  u630=(u2936*(u2886*(u568+u2980*(u3833+u52)+u2885*(u3894+u815)+(u3894+&
 &u3782+u949)*u2886)+u2885*(u439+u4002*u2885)+u2980*(u880+u4032)+u742))
  value_4_=value_4_+(c(26)*u630)
  u3959=6.37738359819483529d-5*u4011*(3.36365613030878777d6*u4118-u3988&
 &*(1.03d2*pd2-1.76d2))
  u598=-1.89408292866386608d-2*u61*(u3991+u3988)
  u508=2.02035512390812382d-1*u3980
  u995=-6.73451707969374607d-2*u472
  u4120=1.72189357151260553d-3*u61*(1.40424479226483373d6*u474+u3988)
  u567=u61*(u943+u3863)
  u746=-5.05088780977030955d-2*u567
  u556=6.31360976221288694d-3*u472
  u4100=-1.14792904767507035d-3*u61
  u486=-3.00608173113575518d3*exppd2
  u388=+u4100*(u4117+u486)
  u4117=1.89408292866386608d-2*u61*(5.71494973596153262d6*u474+u4023)
  u396=-5.68224878599159824d-2*u498
  u3782=u3989*(3.3d1+u3774)
  u867=u421*(1.40097910670142714d7*u474-2.26874092915906051d2*u3782)
  u421=-6.31360976221288694d-3*u4167
  u466=-5.05088780977030955d-2*u842
  u842=-9.4704146433193304d-2*u472
  u4116=1.01017756195406191d-1*u61*(u465+u717*(1.2d1+pd2))
  u465=-1.89408292866386608d-2*u4167
  u429=-1.26272195244257739d-2*u4167
  u3830=u429*u2980
  u3831=u421*u2980
  u4040=u465*u2980
  u4041=u556*u2980
  value_5_=value_5_+(c(26)*(u2886*(u598+u2980*(u2980*(u842+u571)+u4117)&
 &+u2886*((u995+u4027)*u2886+u841*u2980+u508))+u2885*(u4120+u2980*(u298&
 &0*(u4116+u3830)+u396)+u2885*((u556+u3831)*u2885+u2980*(u867+u4040)+u7&
 &46))+u2980*(u388+u2980*(u4041+u466))+u3959))
  u4120=-1.08902111487637123d-2*u475
  u746=2.39584645272801671d-1*u4011*u4042
  u388=-2.79515419484935283d-1*u472
  u4117=1.19792322636400836d-1*u61*(u4136+u3873)
  u466=-2.66205161414224079d-2*u766
  u842=9.31718064949784276d-2*u4167
  u4042=u842*u2980
  value_6_=value_6_+(c(26)*(u2937*(u2886*(u746+u466*u2980+(u4042+u388)*&
 &u2886)+u4117*u2980+u4120)))
  u4136=1.73081334860549274d6*u4118
  u4029=5.3d1*pd2
  u4118=2.85205264883117136d-4*u4011*(u4136-u3988*(4.4d1+u4029))
  u777=1.02542526690966928d6*u474
  u662=-1.28342369197402711d-2*u61*(u453+u777)
  u4027=5.6470642446857193d-2*u833
  u4098=2.56684738394805423d-3*u446
  u4039=-1.8823547482285731d-2*u566
  u4053=9.41177374114286549d-3*u472
  u868=-2.15615180251636555d-1*u61*(u3987+2.95466789070120054d4*u474)
  u4044=8.47059636702857895d-2*u780
  u780=-5.6470642446857193d-2*u723
  u4172=-2.82353212234285965d-2*u498
  u955=9.41177374114286549d-3*u959
  u853=-5.36471103245143333d-1*u472
  u3832=u501*u2980
  u4043=u4053*u2980
  value_7_=value_7_+(c(26)*(u2886*(u662+u2980*(u2980*(u853+u3884)+u4044&
 &)+u2886*((u753+u774)*u2886+u2980*(u780+u3871)+u4027))+u2885*(u4098+u2&
 &980*(u3832+u4172)+u2885*((u4053+u3884)*u2885+u2980*(u955+u3884)+u4039&
 &))+u2980*(u868+u2980*(u4043+u4085))+u4118))
  u780=-2.88127904246118131d-2*u475
  u4172=6.33881389341459887d-1*u4205
  u955=3.16940694670729944d-1*u3996
  u853=-7.04312654823844319d-2*u766
  u4118=-4.22587592894306592d-1*u472
  u868=1.40862530964768864d-1*u4167
  u4044=u868*u2980
  value_8_=value_8_+(c(26)*(u2937*(u2886*(u4172+u2980*(u4044+u853)+u40*&
 &u2886)+u2980*(u955+u4118*u2980)+u780)))
  u780=1.69780995784581975d-2*u553
  u4172=-1.69780995784581975d-2*u670
  u955=1.46955850353296553d6*u474
  u853=u61*(-u3988+u955)
  u4118=-3.39561991569163951d-3*u853
  u4039=u954*exppd2
  u4027=pd2**2
  u4098=u4039*u784*u4027
  u929=2.82472002361899568d0*u4098
  u662=1.24506063575360115d-2*u472
  u784=-6.11211584824495111d-2*u61*(-u3975+2.35855068468253727d5*u474)
  u40=1.82968944468586479d4*u967
  u4122=1.86759095363040173d-1*u498
  u4032=u3989*(u3774-1.5d1)
  u467=-1.24506063575360115d-2*u61*(2.26874092915906051d2*u4032+u468)
  u468=-1.24506063575360115d-2*u4167
  u410=2.98814552580864276d-1*u3980
  u843=-2.98814552580864276d-1*u613
  u869=6.22530317876800576d-2*u4167
  u463=-8.71542445027520806d-2*u472
  u817=7.47036381452160691d-2*u4167
  u587=u817*u2980
  u4045=u468*u2980
  u4046=u869*u2980
  value_9_=value_9_+(c(26)*((u4172+u2980*(u2980*(u789+u3904)+u40))*u288&
 &6+u2885*(u4118+u2980*(u2980*(u843+u587)+u4122)+u2885*((u662+u4045)*u2&
 &885+u2980*(u467+u4046)+u929))+u2980*(u784+u2980*(u463*u2980+u410))+u7&
 &80))
  u3785=1.22242316964899022d-1*u553
  u469=4.98024254301440461d-2*u472
  u952=-1.3582479662766558d-2*u670
  u426=-u3838
  u4047=4.48221828871296415d-1*u472
  u928=-u547
  u547=-2.98814552580864276d-1*u472
  u3895=u547*u2980
  u3896=u928*u2980
  value_1_=value_1_+(c(27)*(u2885*(u3785+u2980*(u3895+u426)+u2885*((u46&
 &9+u3883)*u2885+u2980*(u4047+u3905)+u904))+u2980*(u952+u3896)+u862))
  u3785=u406*u2885
  u952=u3785+u877
  u4047=u997*u2980
  u4048=u852*u2980
  value_2_=value_2_+(c(27)*(x*(u2885*(u852+u4047+u2885*(u952+u4192))+u4&
 &048+u4049)*z))
  value_3_=value_3_+(c(27)*(u2885*(u834+u2980*(u4121*u2980+u3786)+u2885&
 &*((u3965+u3784)*u2885+u3773+u3967))+u2980*(u4050+u4195*u2980)+u859))
  u3786=u393*u2885
  u3967=u3786+u3833
  u4195=u804*u2885
  u4049=u504*u2980
  u4050=u898*u2980
  u834=(u2937*((u494+u2885*(u4195+u4179))*u2886+u2885*(u87+u4049+u2885*&
 &(u3967+u4066))+u4050+u668))
  value_4_=value_4_+(c(27)*u834)
  u3887=u805*u2885
  u3773=u411*u2885
  u3833=u813*u2885
  value_5_=value_5_+(c(27)*(u2934*(u2886*(u629+u2885*(u3773+u660)+(u383&
 &3+u72)*u2886)+u2885*(u656+u2980*(u607+u461)+u2885*(u3887+u3826+u407))&
 &+u2980*(u3964+u3817*u2980)+u855)))
  value_6_=value_6_+(c(27)*u642)
  u642=5.13369476789610845d-3*u870
  u660=-u4014
  u3898=-u875
  u895=5.6470642446857193d-1*u472
  u3817=2.82353212234285965d-2*u472
  u407=-5.6470642446857193d-2*u4167
  u629=-3.38823854681143158d-1*u4205
  u3897=1.8823547482285731d-2*u395
  u827=1.97647248564000175d-1*u472
  u855=u811*u2885
  u585=u407*u2885
  u3940=u407*u2980
  value_7_=value_7_+(c(27)*(u2934*((u660+u2885*(u585+u895))*u2886+u2885&
 &*(u3898+u2980*(u3893+u3897)+u2885*(u855+u3940+u3817))+u2980*(u629+u82&
 &7*u2980)+u642)))
  u642=1.15251161698447252d-1*u553
  u629=-u900
  u895=8.45175185788613183d-1*u472
  value_8_=value_8_+(c(27)*(y*(u2885*(u629+u2980*(u896+u895)+(u569+u715&
 &)*u2885)+u2980*(u629+u3889)+u642)*z))
  u895=u803*u2885
  u629=u895+u454
  u3897=u522*u2980
  value_9_=value_9_+(c(27)*(u2934*(u2885*(u422+u2980*(u3904+u394)+u2885&
 &*(u629+u522))+u2980*(u422+u3897)+u3899)))
  u3898=-u500
  u660=1.99209701720576184d-1*u472
  value_1_=value_1_+(c(28)*(y*(u2885*(u3898+u2980*(u3905+u660)+(u3883+u&
 &469)*u2885)+u2980*(u928+u3892)+u3962)*z))
  u3898=9.60426347487060435d-3*u870
  u394=-u3966
  u900=-3.16940694670729944d-1*u4051
  u419=3.5215632741192216d-1*u472
  u3899=-3.16940694670729944d-1*u554
  u4051=4.22587592894306592d-1*u613
  u500=u2980*(u569+u4051)
  u3877=u2980*(u3899+u639*u2980)
  value_2_=value_2_+(c(28)*(u2934*((u394+u2885*(u3785+u419))*u2886+u288&
 &5*(u900+u500+(u569+u419)*u2885)+u3877+u3898)))
  u3898=u418*u2885
  value_3_=value_3_+(c(28)*(u2936*(u2886*(u961+u2885*(u3898+u79)+(u3898&
 &+u998)*u2886)+u2885*(u3834+u577+(u937+u4105)*u2885)+u2980*(u458+u3832&
 &)+u3850)))
  u3850=u2980*(u439+u4002*u2980)
  u900=(u2934*(u2886*(u4188+u2885*(u3786+u383)+(u4195+u65)*u2886)+u2885&
 &*(u944+u891+(u631+u531)*u2885)+u3850+u3985))
  value_4_=value_4_+(c(28)*u900)
  u4105=u3887+u3891
  u3834=u819*u2980
  value_5_=value_5_+(c(28)*(u2937*(u2886*(u539+u2980*(u3834+u841)+u2885&
 &*(u3773+u4161)+(u3833+u3834+u4109)*u2886)+u2885*(u4025+u2980*(u571+u9&
 &17)+u2885*(u4105+u4019))+u2980*(u801+u4097*u2980)+u397)))
  value_6_=value_6_+(c(28)*u3879)
  u4054=-u572
  u4055=-u540
  u414=1.69411927340571579d-1*u472
  u420=-u3925
  u875=-6.77647709362286316d-1*u679
  u4099=2.25882569787428772d-1*u613
  u643=9.4117737411428655d-2*u472
  u3925=1.78823701081714444d-1*u472
  value_7_=value_7_+(c(28)*(u2937*(u2886*(u4055+u2980*(u3940+u4099)+u28&
 &85*(u585+u998)+(u3940+u414)*u2886)+u2885*(u420+u2980*(u3884+u643)+u28&
 &85*(u855+u3817))+u2980*(u875+u3925*u2980)+u4054)))
  u4054=1.70742461775477411d-2*u993
  u4055=-1.15251161698447252d-1*u4059
  u420=1.05646898223576648d-1*u63
  u875=-3.84170538994824174d-2*u887
  u4099=2.11293796447153296d-1*u498
  u643=-1.05646898223576648d-1*u645
  u3925=-3.5215632741192216d-2*u399
  u645=3.5215632741192216d-2*u3796
  u572=7.04312654823844319d-2*u4167
  u79=(u844+u877)
  value_8_=value_8_+(c(28)*(u2885*(u4055+u2980*(u2980*(u645+u3781)+u409&
 &9)+u2885*(u79*u2885+u2980*(u643+u572*u2980)+u420))+u2980*(u875+u2980*&
 &(u3870*u2980+u3925))+u4054))
  value_9_=value_9_+(c(28)*(x*(u2885*(u904+u2980*(u3904+u4058)+u2885*(u&
 &629+u758))+u2980*(u3791+u4034)+u3962)*z))
  u4055=u61*(u794+u3975)
  u643=-2.7164959325533116d-2*u4055
  u875=4.98024254301440461d-2*u554
  u645=2.7164959325533116d-2*u4055
  u3925=-4.98024254301440461d-2*u554
  value_1_=value_1_+(c(29)*(u2885*(u643+u4191*(u3883+u806)+u2885*((u524&
 &+u3905)*u2885+u3883+u875))+u2980*(u645+u2980*(u469*u2980+u3925))))
  u643=u4011*u3900
  u875=2.88127904246118131d-2*u643
  u645=u61*(2.50369226527838572d5*u474+u3982)
  u3925=-6.33881389341459887d-1*u645
  u4055=-2.07005730208451751d4*u967
  u860=+u4055
  value_2_=value_2_+(c(29)*(u2937*(u2886*(u3925+u500+u2885*(u3785+u4192&
 &)+u3772*u2886)+u2885*(u860+u715*u2885)+u3877+u875)))
  u875=-1.82531369525194967d-2*u787
  u3925=3.08021686073766507d-2*u446
  u860=-2.25882569787428772d-1*u566
  u572=-1.6427823257267547d-1*u915
  u540=1.8823547482285731d-2*u833
  u4051=-1.8823547482285731d-2*u472
  u4054=1.02673895357922169d-2*u4011*u661
  u4099=-3.38823854681143158d-1*u498
  u420=1.12941284893714386d-1*u959
  u3899=1.12941284893714386d-1*u498
  u469=-1.8823547482285731d-2*u723
  u915=1.8823547482285731d-2*u61*(u943+u3873)
  u4056=6.77647709362286316d-1*u472
  u430=3.7647094964571462d-2*u4167
  value_3_=value_3_+(c(29)*(u2886*(u3925+u2980*(u4056*u2980+u4099)+u288&
 &6*((u4115+u3829)*u2886+u2980*(u420+u3829)+u860))+u2885*(u572+u2980*(u&
 &2980*(u469+u3816)+u3899)+u2885*((u4051+u3816)*u2885+u2980*(u469+u430*&
 &u2980)+u540))+u2980*(u4054+u2980*(u4051*u2980+u915))+u875))
  u875=(u2937*(u2886*(u568+u891+u2885*(u3786+u52)+(u4195+u631+u949)*u28&
 &86)+u2885*(u880+u944*u2885)+u3850+u742))
  value_4_=value_4_+(c(29)*u875)
  u469=8.03550333372549246d-3*u61*(5.92488666503767056d5*u474-u3861)
  u430=-2.84532046697493236d4*u967
  u572=5.05088780977030955d-2*u472
  u860=-1.51526634293109286d-1*u61*(u4201+u717)
  u4054=2.39917170964089704d-1*u472
  u4051=3.15680488110644347d-2*u472
  u3925=5.05088780977030955d-2*u613
  u3899=u813*u2886
  u915=u3899+u3773
  value_5_=value_5_+(c(29)*(u2934*(u2886*(u430+u2980*(u571+u4054)+u2885&
 &*(u3887+u4054)+u2886*(u915+u3783+u572))+u2885*(u860+u2980*(u3830+u392&
 &5)+(u3830+u4051)*u2885)+u2980*(u860+u4051*u2980)+u469)))
  value_6_=value_6_+(c(29)*u630)
  u4051=-2.4588714905999591d3*u967
  u430=+u4051
  u3925=2.44706117269714503d-1*u472
  u860=-9.41177374114286549d-3*u472
  u4054=-u4051
  u4051=-2.44706117269714503d-1*u472
  u469=u585+u774
  value_7_=value_7_+(c(29)*(u2934*(u2886*(u2980*(u3884+u4051)+u2885*(u8&
 &55+u3925)+(u469)*u2886)+u2885*(u430+u860*u2885)+u2980*(u4054+u4043)))&
 &)
  u430=-2.88127904246118131d-2*u643
  u3925=6.33881389341459887d-1*u645
  u4054=-u4055
  u4051=u807*u2885
  value_8_=value_8_+(c(29)*(u2936*(u2886*(u3925+u4162+u2885*(u4051+u865&
 &)+(u4051+u4052)*u2886)+u2885*(u4057+u4052*u2885)+u2980*(u4054+u3889)+&
 &u430)))
  u430=-6.79123983138327901d-3*u870
  u3925=1.21979296312390986d4*u967
  u4054=8.9644365774259283d-1*u679
  u4052=-1.24506063575360115d-1*u472
  u4053=-2.36561520793184219d-1*u472
  value_9_=value_9_+(c(29)*(u2934*((u3925+u2980*(u3904+u4052)+u2885*(u8&
 &95+u4052))*u2886+u2885*(u4054+u2980*(u587+u843)+(u587+u4053)*u2885)+u&
 &2980*(u4054+u4053*u2980)+u430)))
  u430=2.39051642064691421d0*u615
  u3925=-u835
  u4052=u806*u2885
  value_1_=value_1_+(c(30)*(u2936*(u2886*(u430+u2980*(u3905+u547)+u2885&
 &*(u4052+u866)+(u4052+u51)*u2886)+u2885*(u3957+u51*u2885)+u2980*(u3925&
 &+u3892)+u426)))
  u3925=-5.76255808492236261d-2*u670
  u645=4.14011460416903502d4*u967
  u4053=u4192*u2885
  u4054=u776*u2885
  value_2_=value_2_+(c(30)*(u2934*(u2886*(u645+u4047+u4053+(u952+u3936)&
 &*u2886)+u4054+u4048+u3925)))
  u3925=6.20480257047252114d5*u474
  u4047=9.07496371663624206d2*exppd2
  u4099=-1.54010843036883254d-2*u61*(u4047+u3925)
  u420=1.12941284893714386d-1*u61*(u470+u3982)
  u540=4.42596868307992638d4*u967
  u643=-1.01647156404342947d0*u472
  u4055=u814*u2885
  u4056=u753*u2885
  value_3_=value_3_+(c(30)*(u2936*(u2886*(u420+u2980*(u937+u643)+u2885*&
 &(u4055+u826)+(u4055+u3871+u747)*u2886)+u2885*(u570+u4056)+u2980*(u540&
 &+u3832)+u4099)))
  u397=-7.26014076584247488d-3*u670
  u539=3.12963246912510941d4*u967
  u801=-5.59030838969870566d-1*u472
  u500=u804*u2886+u3967
  u4057=u504*u2885
  u4058=u898*u2885
  u835=(u2934*(u2886*(u539+u4049+u4057+u2886*(u500+u801))+u4058+u4050+u&
 &397))
  value_4_=value_4_+(c(30)*u835)
  u4097=1.14792904767507035d-3*u61
  u4050=7.25997097330899365d3*exppd2
  u850=u4097*(5.32306746835274182d6*u474+u4050)
  u524=-3.03053268586218573d-1*u4011*u903
  u63=-1.13644975719831965d-1*u472
  u679=-2.10306295385103696d4*u967
  u399=4.67207122403753633d-1*u472
  u3850=-7.57633171465546432d-2*u61*(u3991+u3982)
  u787=1.26272195244257739d-2*u61
  u4049=u3989*(5.5d1+u3774)
  u870=u787*(2.3349651778357119d7*u474-2.26874092915906051d2*u4049)
  u446=-6.31360976221288694d-2*u4167
  u4059=u446*u2980
  value_5_=value_5_+(c(30)*(u2937*(u2886*(u524+u2980*(u3831+u870)+u2885&
 &*(u3887+u399)+u2886*(u915+u4059+u63))+u2885*(u679+u947*u2885)+u2980*(&
 &u3850+u4019*u2980)+u850)))
  u915=4.03341153657915271d-4*u4011*(u4060-u3988*(4.4d1+u3921))
  u3902=+u4006*(1.47468160395338933d3*exppd2+u91)
  u841=u3847*(u91+u3873)
  u904=-9.31718064949784276d-2*u472
  u917=-3.63007038292123744d-3*u475
  u3847=1.19792322636400836d-1*u754
  u754=u3989*(3.5d1+u3774)
  u458=-3.99307742121336118d-2*u61*(u4062-2.26874092915906051d2*u754)
  value_6_=value_6_+(c(30)*(u2886*(u3902+u3847*u2980+u2886*((u904+u4042&
 &)*u2886+u458*u2980+u841))+u917*u2980+u915))
  u917=-1.23208674429506603d-1*u670
  u3847=3.31947651230994478d4*u967
  u458=-1.12941284893714386d-1*u472
  u4025=-2.21298434153996319d4*u967
  u3889=+u4025
  u3877=5.08235782021714737d-1*u472
  u997=1.84415361794996932d4*u967
  u661=-4.70588687057143275d-1*u472
  u3900=u3817*u2980
  u4060=u966*u2885
  value_7_=value_7_+(c(30)*(u2937*(u2886*(u3847+u2980*(u3884+u661)+u288&
 &5*(u855+u3877)+(u469+u458)*u2886)+u2885*(u3889+u4060)+u2980*(u997+u39&
 &00)+u917)))
  u458=2.88127904246118131d-2*u553
  u3889=u4011*u386
  u3877=-2.88127904246118131d-2*u3889
  u661=u61*(u3925+u3873)
  u386=1.05646898223576648d-1*u661
  u3925=3.16940694670729944d-1*u61*(u3901+u3862)
  u3991=-1.05646898223576648d-1*u61*(9.76439983458570431d6*u474+u3873*(&
 &2.3d1+u3774))
  u470=1.38003820138967834d4*u967
  value_8_=value_8_+(c(30)*(u2886*(u3877+u2980*(u4061*u2980+u3925)+u288&
 &6*(u79*u2886+u2980*(u3991+u4044)+u386))+u2980*(u3835+u470*u2980)+u458&
 &))
  u458=-2.19562733362303775d4*u967
  u3877=3.58577463097037132d0*u615
  u386=-2.24110914435648207d-1*u472
  u3925=7.31875777874345916d3*u967
  u3991=4.48221828871296415d-1*u4205
  u470=-2.4901212715072023d-2*u395
  u395=-2.61462733508256242d-1*u472
  u844=8.71542445027520806d-2*u4167
  u4061=u844*u2980
  value_9_=value_9_+(c(30)*(u2937*(u2886*(u3877+u2980*(u4061+u470)+u288&
 &5*(u895+u499)+(u587+u386)*u2886)+u2885*(u3925+u522*u2885)+u2980*(u399&
 &1+u395*u2980)+u458)))
  u3835=2.71649593255331161d-1*u553
  u51=-u646
  u4062=8.9644365774259283d-1*u472
  u3901=u402*u2885
  u646=u3901+u3905
  value_1_=value_1_+(c(31)*(u2934*(u2885*(u51+u4063*u2980+u2885*(u646+u&
 &4062))+u3836*u2980+u3835)))
  u3835=1.92085269497412087d-1*u553
  u51=-u3876
  u4062=u968*u2980
  u4063=u4112*u2980
  value_2_=value_2_+(c(31)*(y*(u2885*(u51+u4062+u2885*(u952+u797))+u406&
 &3+u3835)*z))
  u3835=u403*u2885
  u51=u3835+u937
  u3836=u812*u2885
  value_3_=value_3_+(c(31)*(u2934*((u4175+u2885*(u3836+u4064))*u2886+u2&
 &885*(u451+u4065*u2980+u2885*(u51+u998))+u828*u2980+u874)))
  u4064=u4114*u2980
  u4065=u398*u2980
  value_4_=value_4_+(c(31)*(u2936*((u440+u2885*(u4195+u560))*u2886+u288&
 &5*(u546+u4064+u2885*(u3967+u999))+u4065+u965)))
  value_5_=value_5_+(c(31)*(u2886*(u3811+u2885*(u2885*(u654+u3773)+u614&
 &)+(u2885*(u4068+u3833)+u509)*u2886)+u2885*(u54+u2980*(u628*u2980+u406&
 &7)+u2885*(u2885*(u3837+u3826+u3887)+u2980*(u460+u607)+u56))+u2980*(u7&
 &93+u4113*u2980)+u3971))
  value_6_=value_6_+(c(31)*u834)
  u965=-5.70410529766234272d-4*u4011*u473
  u4068=7.70054215184416268d-2*u553
  u87=2.56684738394805423d-2*u82
  u3837=-u416
  u654=2.82353212234285965d-2*u61*(u943-u3862)
  u4066=8.47059636702857895d-1*u472
  u874=-1.50588379858285848d-1*u472
  u868=6.16043372147533014d-2*u47
  u4067=-1.69411927340571579d-1*u765
  u903=1.12941284893714386d-1*u95
  u82=-u808
  value_7_=value_7_+(c(31)*((u4068+u2885*(u2885*(u4066+u585)+u3837))*u2&
 &886+u2885*(u87+u2980*(u4030+u4067)+u2885*(u2885*(u874+u3940+u855)+u29&
 &80*(u903+u3893)+u654))+u2980*(u868+u82*u2980)+u965))
  u4067=-u452
  u903=-u578
  u82=-u3963
  u3837=u400*u2885
  u47=u3837+u896
  u4066=u82*u2980
  value_8_=value_8_+(c(31)*(x*(u2885*(u4067+u3890+u2885*(u47+u903))+u40&
 &66+u642)*z))
  value_9_=value_9_+(c(31)*(u2885*(u392+u2980*(u499*u2980+u527)+u2885*(&
 &u2885*(u871+u454+u895)+u2980*(u3778+u3904)+u4130))+u2980*(u523+u3939*&
 &u2980)+u3875))
  u4067=-u541
  u903=5.97629105161728553d-1*u472
  value_1_=value_1_+(c(32)*(x*(u2885*(u4067+u3895+u2885*(u646+u903))+u3&
 &896+u3962)*z))
  u4067=-1.06714038609673382d-3*u4011*u4126
  u903=4.80213173743530218d-2*u553
  u3778=u4011*u651
  u868=5.76255808492236261d-2*u3778
  u4068=-1.05646898223576648d-1*u974
  u3962=5.28234491117883239d-1*u472
  u974=1.40862530964768864d-1*u472
  u3875=-3.16940694670729944d-1*u498
  u871=1.05646898223576648d-1*u959
  u874=u2980*(u797*u2980+u3875)
  u3811=u2980*(u871+u569)
  u965=u2980*(u852+u404*u2980)
  value_2_=value_2_+(c(32)*((u903+u2885*(u2885*(u3962+u3785)+u4202))*u2&
 &886+u2885*(u868+u874+u2885*((u974+u569)*u2885+u3811+u4068))+u965+u406&
 &7))
  u4067=u3788*u2980
  value_3_=value_3_+(c(32)*(u2937*((u4069+u2885*(u3836+u471))*u2886+u28&
 &85*(u4014+u2980*(u3835+u4115)+u403*u800)+u4067+u994)))
  u3962=u2980*(u4134*u2980+u538)
  u868=u2980*(u898+u985*u2980)
  value_4_=value_4_+(c(32)*(u2886*(u41+u2885*(u2885*(u3937+u3786)+u985)&
 &+(u2885*(u4179+u4195)+u494)*u2886)+u2885*(u873+u3962+u2885*(u2980*(u3&
 &961+u3894)+u680))+u868+u767))
  u4068=u947*u2980
  value_5_=value_5_+(c(32)*(u2936*(u2886*(u909+u462*u2885+(u840*u2885+u&
 &658)*u2886)+u2885*(u663+u2980*(u571+u92)+u2885*(u4105+u872))+u2980*(u&
 &4186+u4068)+u586)))
  value_6_=value_6_+(c(32)*u900)
  u903=2.56684738394805423d-2*u61*(u876-u3861)
  u872=-6.77647709362286316d-1*u3980
  u680=-3.38823854681143158d-1*u576
  u4179=1.12941284893714386d-1*u941
  u909=1.41176606117142983d-1*u472
  u3937=-u4085
  value_7_=value_7_+(c(32)*(u2936*(u2886*(u872+u2885*(u3898+u4179)+(u58&
 &5+u414)*u2886)+u2885*(u680+u2980*(u3884+u501)+u2885*(u855+u909))+u298&
 &0*(u3937+u3900)+u903)))
  u903=2.88127904246118131d-2*u4011
  u872=u903*u4183
  u680=-1.05646898223576648d-1*u66
  u4179=-1.05646898223576648d-1*u554
  u909=1.40862530964768864d-1*u613
  u3937=u2980*(u3781+u909)
  u414=u2980*(u4179+u4035)
  value_8_=value_8_+(c(32)*(u2934*((u4202+u2885*(u3837+u505))*u2886+u28&
 &85*(u680+u3937+u575*u2885)+u414+u872)))
  value_9_=value_9_+(c(32)*(y*(u2885*(u3793+u2980*(u3904+u785)+u2885*(u&
 &629+u789))+u2980*(u3779+u3897)+u667)*z))
  u872=1.3582479662766558d-2*u671
  u680=-1.49407276290432138d-1*u3996
  u471=-1.49407276290432138d-1*u554
  u658=1.99209701720576184d-1*u613
  u785=u2980*(u3883+u658)
  u554=u2980*(u471+u408*u2980)
  value_1_=value_1_+(c(33)*(u2934*((u422+u2885*(u3901+u604))*u2886+u288&
 &5*(u680+u785+(u3883+u660)*u2885)+u554+u872)))
  u660=2.88127904246118131d-2*u671
  u680=-6.33881389341459887d-1*u4205
  u872=-3.16940694670729944d-1*u3996
  u41=7.04312654823844319d-2*u766
  u994=4.22587592894306592d-1*u472
  u873=-1.40862530964768864d-1*u4167
  value_2_=value_2_+(c(33)*(u2936*(u2886*(u680+u2885*(u873*u2885+u41)+(&
 &u3837+u639)*u2886)+u2885*(u872+u994*u2885)+u660)))
  u660=u2980*(u570+u4028)
  value_3_=value_3_+(c(33)*(u2934*(u2886*(u973+u2885*(u3835+u3827)+(u38&
 &36+u520)*u2886)+u2885*(u4119+u431+(u3816+u447)*u2885)+u660+u3795)))
  value_4_=value_4_+(c(33)*(u2936*(u2886*(u746+u466*u2885+(u842*u2885+u&
 &388)*u2886)+u4117*u2885+u4120)))
  u680=u4097*(3.55959726411318318d6*u474-u486)
  u872=-2.53604650317330928d4*u967
  u41=-1.01017756195406191d-1*u567
  u994=2.08349122153025269d-1*u472
  u873=-1.72189357151260553d-3*u853
  u4167=1.51526634293109286d-1*u498
  u388=-5.05088780977030955d-2*u959
  u383=1.43239448782705802d0*u4098
  u604=-3.03053268586218573d-1*u472
  value_5_=value_5_+(c(33)*(u2886*(u598+u2980*(u604*u2980+u4167)+u2885*&
 &(u2885*(u994+u3887)+u872)+u2886*((u995+u3834+u3833)*u2886+u2885*(u572&
 &+u3773)+u2980*(u388+u3834)+u508))+u2885*(u680+u2980*(u2980*(u867+u383&
 &1)+u396)+u2885*((u556+u3830)*u2885+u2980*(u4116+u4040)+u41))+u2980*(u&
 &873+u2980*(u4041+u383))+u3959))
  value_6_=value_6_+(c(33)*u875)
  u439=-2.85205264883117136d-4*u4011*(-u3988*(u4029+4.4d1)+u4136)
  u872=1.28342369197402711d-2*u61*(u777+u453)
  u4116=-5.6470642446857193d-2*u833
  u388=1.07807590125818278d-1*u61*(7.77544181763473826d3*u474-u3975)
  u944=-5.25583781115741257d4*u967
  u4167=2.56684738394805423d-3*u61*(u955-u3988)
  u949=-1.69411927340571579d-1*u498
  u52=5.6470642446857193d-2*u959
  u995=2.82353212234285965d-2*u498
  u396=-2.13528763025153118d0*u4098
  u4002=-9.41177374114286549d-3*u959
  value_7_=value_7_+(c(33)*(u2886*(u872+u2980*(u998*u2980+u949)+u2885*(&
 &u2885*(u827+u855)+u944)+u2886*((u501+u3940)*u2886+u2885*(u998+u585)+u&
 &2980*(u52+u3940)+u4116))+u2885*(u388+u2980*((u753+u3828)*u2885+u2980*&
 &(u4002+u3828)+u995)+u860*u800)+u2980*(u4167+u2980*(u860*u2980+u396))+&
 &u439))
  u949=2.88127904246118131d-2*u3778
  u872=-2.11293796447153296d-1*u566
  u860=-u3926
  value_8_=value_8_+(c(33)*(u2937*(u2886*(u872+u3937+u2885*(u3837+u797)&
 &+u897*u2886)+u2885*(u860+u639*u2885)+u414+u949)))
  value_9_=value_9_+(c(33)*((u4172+u2885*(u2885*(u789+u895)+u40))*u2886&
 &+u2885*(u784+u2980*(u2980*(u467+u4045)+u4122)+u2885*((u463+u587)*u288&
 &5+u2980*(u843+u4046)+u410))+u2980*(u4118+u2980*(u662*u2980+u929))+u78&
 &0))
  u872=-u430
  value_1_=value_1_+(c(34)*(u2937*(u2886*(u872+u785+u2885*(u3901+u938)+&
 &u711*u2886)+u2885*(u3779+u408*u2885)+u554+u3838)))
  u872=-2.88127904246118131d-2*u670
  u860=u903*u4070
  u949=-1.05646898223576648d-1*u661
  u52=5.76255808492236261d-2*u553
  u959=-u645
  value_2_=value_2_+(c(34)*(u2886*(u860+u874+u2885*(u4053+u959)+u2886*(&
 &(u715+u569)*u2886+u2885*(u4192+u3785)+u3811+u949))+u2885*(u52+u4054)+&
 &u965+u872))
  value_3_=value_3_+(c(34)*(u2937*(u2886*(u420+u431+u2885*(u3835+u643)+&
 &(u3836+u3816+u747)*u2886)+u2885*(u540+u501*u2885)+u660+u4099)))
  value_4_=value_4_+(c(34)*(u2886*(u3902+u3962+u2885*(u4057+u539)+u2886&
 &*((u904+u631+u4195)*u2886+u2885*(u801+u3786)+u845+u841))+u2885*(u397+&
 &u4058)+u868+u915))
  value_5_=value_5_+(c(34)*(u2936*(u2886*(u524+u2980*(u571+u399)+u2885*&
 &(u421*u2885+u870)+u2886*(u3899+u446*u2885+u3783+u63))+u2885*(u3850+u4&
 &019*u2885)+u2980*(u679+u4068)+u850)))
  value_6_=value_6_+(c(34)*u835)
  u872=1.23208674429506603d-1*u553
  u4134=-u3847
  u539=-u997
  u949=4.70588687057143275d-1*u472
  u747=-u4025
  u860=-5.08235782021714737d-1*u472
  value_7_=value_7_+(c(34)*(u2936*(u2886*(u4134+u2980*(u3884+u860)+u288&
 &5*(u855+u949)+(u469+u4115)*u2886)+u2885*(u539+u4060)+u2980*(u747+u390&
 &0)+u872)))
  value_8_=value_8_+(c(34)*(u2934*(u2886*(u959+u3890+u797*u2885+(u47+u4&
 &192)*u2886)+u404*u2885+u4066+u52)))
  value_9_=value_9_+(c(34)*(u2936*(u2886*(u3877+u2980*(u3904+u499)+u288&
 &5*(u844*u2885+u470)+(u817*u2885+u386)*u2886)+u2885*(u3991+u395*u2885)&
 &+u2980*(u3925+u3897)+u458)))
  value_1_=value_1_+(c(35)*(u2934*(u2886*(u3895+u938*u2885+(u646)*u2886&
 &)+u3995*u2885+u3896)))
  value_2_=value_2_+(c(35)*(u2936*(u2886*(u4062+u419*u2885+(u952)*u2886&
 &)+u394*u2885+u4063)))
  u539=1.43843982200097607d5*u967
  u949=-1.5811779885120014d0*u472
  value_3_=value_3_+(c(35)*(u2934*(u2886*(u539+u4115*u2980+u4115*u2885+&
 &u2886*(u812*u2886+u51+u949))+u3788*u2885+u4067+u917)))
  u539=-2.90405630633698995d-1*u670
  u949=1.30401352880212892d5*u967
  u747=-9.58338581091206684d-1*u472
  u860=(u2886*(u949+u4064+u4114*u2885+u2886*(u500+u747))+u398*u2885+u40&
 &65+u539)
  value_4_=value_4_+(c(35)*(u2936*u860))
  u872=1.27547671963896706d-4*u4011*(u820-u3988*(1.76d2+u4021))
  u968=+u4100*(2.50695872672076187d4*exppd2+5.06181262328021461d6*u474)
  u52=u787*(4.7352440669395556d6*u474-u3873)
  u4134=-3.5777121985873026d-1*u472
  u4167=5.8544381431428588d-2*u553
  u959=-6.68031761811505858d4*u967
  u915=7.19751512892269111d-1*u472
  u547=6.18547927603246165d2*u967
  u658=-3.78816585732773216d-2*u472
  u522=1.37751485721008442d-2*u61*(u62+u3874)
  u4002=-3.78816585732773216d-2*u61*(2.64520530635933796d6*u474+u3862)
  u3936=2.52544390488515477d-2*u61*(1.65570258064714117d7*u474+u3873*(3&
 &.9d1+pd2))
  u966=-u547
  u947=3.78816585732773216d-2*u472
  value_5_=value_5_+(c(35)*(u2886*(u968+u2980*(u947*u2980+u4002)+u2885*&
 &(u658*u2885+u959)+u2886*(u2886*(u4134+u4059+u3773+u3899)+u2885*(u915+&
 &u3887)+u2980*(u3936+u3831)+u52))+u2885*(u4167+u547*u2885)+u2980*(u522&
 &+u966*u2980)+u872))
  value_6_=value_6_+(c(35)*(u2937*u860))
  u915=6.16043372147533014d-2*u553
  u3936=-7.19219911000488036d4*u967
  u658=+u3936
  u947=9.22076808974984662d2*u967
  u522=-6.16043372147533014d-2*u670
  u4002=-u3936
  u3936=-7.90588994256000702d-1*u472
  u966=-u947
  value_7_=value_7_+(c(35)*(u2886*(u2980*(u3832+u4002)+u2885*(u4056+u65&
 &8)+u2886*((u774+u585)*u2886+u2885*(u4121+u855)+u2980*(u3936+u3884)))+&
 &u2885*(u915+u947*u2885)+u2980*(u522+u966*u2980)))
  value_8_=value_8_+(c(35)*(u2937*(u2886*(u4038+u505*u2885+(u47)*u2886)&
 &+u4202*u2885+u3966*u2980)))
  u915=2.0373719494149837d-2*u553
  u658=-2.0373719494149837d-2*u3889
  u522=7.47036381452160691d-2*u661
  u4002=-2.0373719494149837d-2*u670
  u3936=-8.14948779765993481d-2*u887
  u966=2.24110914435648207d-1*u765
  u947=-1.49407276290432138d-1*u95
  u872=8.53855074186736902d3*u967
  u4167=-5.22925467016512484d-1*u472
  value_9_=value_9_+(c(35)*(u2886*(u658+u2980*(u4167*u2980+u966)+u2885*&
 &(u499*u2885+u3838)+u2886*((u499+u587)*u2886+u2885*(u499+u895)+u2980*(&
 &u947+u4061)+u522))+u2885*(u4002+u3939*u2885)+u2980*(u3936+u872*u2980)&
 &+u915))
  if ( lmax .eq. 4 ) go to 100
  u915=A13_v(pd,exppd2,erfpd)
  u658=pd2*u915
  u522=4.24539123242856709d5*u658
  u4002=p**10*pd
  u3936=u954*u4002
  u966=(1.3d1+u3774)
  u947=u3936*(u522-u3988*u966)
  u872=2.0896122558102397d-2*u947
  u4061=4.24539123242856709d5*u915
  u4167=u3936*pd2
  u959=u4167*(u4061-u4003)
  u95=5.74643370347815917d-2*u959
  u753=u4167*(-u4003+u4061)
  u4061=-1.09182240366085024d0*u753
  u747=u3936*u658
  u539=-9.51438511236649692d5*u747
  u949=+u539
  u52=3.4886078745343822d6*u747
  u547=6.36808684864285064d6*u915
  u968=u4167*(u547+u3862)
  u4134=7.47036381452160692d-1*u968
  u530=-1.24506063575360115d0*u968
  u860=+u530
  u801=u4167*(1.08257476426928461d8*u915+u3862*(1.7d1+u3774))
  u995=-4.98024254301440461d-2*u801
  u396=4.98024254301440461d-2*u801
  u388=u396*u2980
  u463=u995*u2980
  value_1_=value_1_+(c(36)*(y*((u95+u2980*(u2980*(u4134+u463)+u949))*u2&
 &885+u2980*(u4061+u2980*(u2980*(u860+u388)+u52))+u872)))
  u4061=-7.3140160308629987d-1*u753
  u860=-2.24256207725822731d5*u747
  u643=+u860
  u524=4.26086794679063188d6*u747
  u63=3.5215632741192216d-1*u968
  u3850=-2.11293796447153296d0*u968
  u904=+u3850
  u386=-3.5215632741192216d-2*u801
  u395=1.05646898223576648d-1*u801
  u3875=u395*u2980
  u4179=u386*u2980
  value_2_=value_2_+(c(36)*(u2934*((u643+u2980*(u4179+u63))*u2885+u2980&
 &*(u524+u2980*(u3875+u904))+u4061)*z))
  u4061=(u3774+1.3d1)
  u904=u3936*(-u3988*u4061+u522)
  u471=-7.89799195060939762d-3*u904
  u4116=-1.30316867185055061d-1*u753
  u944=2.17194778641758435d-2*u959
  u504=4.12670079419341026d-1*u959
  u4192=2.15765973300146411d6*u747
  u715=-3.59609955500244018d5*u747
  u439=+u715
  u501=-1.31856983683422807d6*u747
  u998=+u501
  u4114=-1.69411927340571579d0*u968
  u819=2.82353212234285965d-1*u968
  u4115=4.70588687057143275d-1*u968
  u408=1.12941284893714386d-1*u801
  u639=-1.8823547482285731d-2*u801
  u3817=u639*u2980
  u841=u4115+u3817
  u383=u408*u2980
  value_3_=value_3_+(c(36)*(y*((u4116+u2980*(u2980*(u4114+u383)+u4192))&
 &*u2886+(u944+u2980*(u2980*(u819+u3817)+u439))*u2885+u2980*(u504+u2980&
 &*(u2980*(u841)+u998))+u471)))
  u4019=2.76443821468617313d-1*u959
  u570=3.39043517488553519d5*u747
  u742=-2.5428263811641514d5*u747
  u568=+u742
  u556=-1.61045670807062922d6*u747
  u4121=-5.32410322828448158d-1*u968
  u538=3.99307742121336118d-1*u968
  u4046=7.98615484242672237d-1*u968
  u397=5.32410322828448158d-2*u801
  u604=-3.99307742121336118d-2*u801
  u3778=u604*u2980
  u994=u397*u2980
  u903=(u934*((u570+u2980*(u994+u4121))*u2886+(u568+u2980*(u3778+u538))&
 &*u2885+u2980*(u556+u2980*(u3778+u4046))+u4019))
  value_4_=value_4_+(c(36)*u903)
  u909=2.97177386269999696d6*u658
  u407=2.8d1*pd2
  u41=(3.9d1+u407)
  u873=1.05962681323852648d-2*u3936*(u909-u3975*u41)
  u680=2.91397373640594782d-1*u959
  u3937=1.07214974117896002d5*u747
  u785=2.6321425641057116d6*u915
  u3795=u4167*(-u4003+u785)
  u4119=-1.45698686820297391d-1*u3795
  u3786=3.61850537647899007d5*u747
  u508=9.90591287566665655d5*u915
  u4120=-3.78816585732773216d-1*u4167*(-u4003+u508)
  u468=-1.82265456000423203d6*u747
  u4117=-1.68362926992343652d-1*u968
  u4172=6.31360976221288694d-2*u4167*(6.75017205956142168d7*u915+u4023)
  u4118=-5.68224878599159824d-1*u968
  u662=6.31360976221288694d-3*u4167
  u784=u662*(1.93165301075499803d8*u915+u4023)
  u868=9.59668683856358814d-1*u968
  u410=1.68362926992343652d-2*u801
  u827=u4167*(2.10146866005214071d8*u915-4.53748185831812103d2*u3782)
  u660=-5.05088780977030955d-2*u827
  u414=5.68224878599159824d-2*u801
  u3898=-2.39917170964089704d-1*u968
  u814=-5.05088780977030955d-2*u801
  u420=6.31360976221288694d-2*u801
  u399=6.31360976221288694d-3*u801
  u3902=u410*u2980
  u917=u4117+u3902
  u458=u399*u2980
  u3877=u414*u2980
  u3991=u814*u2980
  u520=u420*u2980
  value_5_=value_5_+(c(36)*(x*(u2886*(u680+u2980*(u2980*(u868+u3991)+u4&
 &68)+(u2980*(u917)+u3937)*u2886)+u2885*(u4119+u2980*(u2980*(u660+u520)&
 &+u4172)+(u2980*(u4118+u3877)+u3786)*u2885)+u2980*(u4120+u2980*(u2980*&
 &(u3898+u458)+u784))+u873)))
  u472=-1.6754170998098019d-2*u904
  u498=-6.14319603263594028d-2*u753
  u554=4.60739702447695521d-2*u959
  u566=8.7540543465062149d-1*u959
  u833=1.01713055246566056d6*u747
  u4205=-7.62847914349245418d5*u747
  u615=-2.79710901928056653d6*u747
  u4058=-7.98615484242672237d-1*u968
  u3980=5.98961613182004178d-1*u968
  u3962=9.98269355303340296d-1*u968
  value_6_=value_6_+(c(36)*(z*((u498+u2980*(u2980*(u4058+u994)+u833))*u&
 &2886+(u554+u2980*(u2980*(u3980+u3778)+u4205))*u2885+u2980*(u566+u2980&
 &*(u2980*(u3962+u3778)+u615))+u472)))
  u3996=(u4005-2.6d1)
  u765=-u4003*u3996
  u661=u3936*(u909+u765)
  u874=5.92349396295704821d-3*u661
  u3811=-4.56109035147692713d-1*u753
  u965=2.97177386269999696d6*u915
  u523=u4167*(-u3988+u965)
  u66=-1.30316867185055061d-1*u523
  u576=4.19544948083618021d5*u747
  u974=-4.31060776540221498d3*exppd2
  u567=u4167*(u965+u974)
  u780=-2.17194778641758435d-2*u567
  u867=2.51726968850170813d6*u747
  u850=1.48588693134999848d7*u915
  u3827=u4167*(u850+u3872)
  u553=2.82353212234285965d-1*u3827
  u670=-6.58824161880000585d-1*u968
  u887=-2.82353212234285965d-2*u4167
  u671=+u887*(-u3872+u850)
  u47=-1.18588349138400105d0*u968
  u82=+u47
  u853=4.45766079404999545d7*u915
  u802=u3989*(1.4d1+pd2)
  u654=u4167*(u853-2.26874092915906051d2*u802)
  u898=-2.25882569787428772d-1*u654
  u422=6.58824161880000585d-2*u801
  u398=2.25882569787428772d-1*u968
  u404=5.6470642446857193d-2*u801
  u87=-9.41177374114286549d-3*u801
  u4202=u87*u2980
  u3995=u404*u2980
  u3788=u422*u2980
  value_7_=value_7_+(c(36)*(x*((u3811+u2980*(u2980*(u82+u3995)+u867))*u&
 &2886+u2885*(u66+u2980*(u2980*(u898+u3995)+u553)+(u2980*(u670+u3788)+u&
 &576)*u2885)+u2980*(u780+u2980*(u2980*(u398+u4202)+u671))+u874)))
  u874=1.47757899613393913d-2*u947
  u3811=1.21900267181049978d-1*u959
  u66=-7.72035025479983196d-1*u753
  u780=-2.01830586953240457d6*u747
  u553=+u780
  u670=2.46681828498405004d6*u747
  u671=1.58470347335364972d0*u968
  u82=-8.803908185298054d-1*u968
  u898=-1.05646898223576648d-1*u801
  u398=3.5215632741192216d-2*u801
  u776=u398*u2980
  u3779=u898*u2980
  value_8_=value_8_+(c(36)*(((u3811+u2980*(u2980*(u671+u3779)+u553))*u2&
 &885+u2980*(u66+u2980*(u2980*(u82+u776)+u670))+u874)*z))
  u874=3.13441838371535955d-2*u947
  u66=5.17179033313034325d-1*u959
  u82=7.9286542603054141d4*u747
  u4040=-6.32107707382597508d-1*u753
  u451=-3.01288861891605736d6*u747
  u614=-1.24506063575360115d-1*u968
  u546=1.3478712242519204d6*u747
  u3828=1.49407276290432138d0*u968
  u394=1.24506063575360115d-2*u801
  u679=-3.73518190726080346d-1*u968
  u4025=-7.47036381452160691d-2*u801
  u3876=u394*u2980
  u56=u2980*(u614+u3876)
  u840=u4025*u2980
  value_9_=value_9_+(c(36)*(x*(u2885*(u66+u2980*(u2980*(u3828+u840)+u45&
 &1)+(u56+u82)*u2885)+u2980*(u4040+u2980*(u2980*(u679+u3876)+u546))+u87&
 &4)))
  u3884=-5.74643370347815917d-2*u753
  u4112=-3.17146170412216564d5*u747
  u4188=+u4112
  u880=-2.87321685173907958d-1*u753
  u928=-u539
  u539=4.98024254301440461d-1*u968
  u4175=3.80575404494659877d5*u747
  u4113=-7.47036381452160692d-1*u968
  u440=-4.98024254301440461d-2*u968
  u4130=u3884+u2980*(u2980*(u4113+u388)+u928)
  u4069=u440*u2980
  value_1_=value_1_+(c(37)*(x*(u2885*(u4130+(u2980*(u539+u463)+u4188)*u&
 &2885)+u2980*(u880+u2980*(u4069+u4175))+u872)))
  u880=8.86547397680363479d-3*u947
  u4186=-7.3140160308629987d-2*u753
  u973=-4.48512415451645461d4*u747
  u40=+u973
  u855=-3.16940694670729944d-1*u753
  u3773=1.48009097099043002d6*u747
  u628=2.11293796447153296d-1*u968
  u3959=6.27917381632303645d5*u747
  u967=-1.37340967690649642d0*u968
  u473=+u967
  u4126=-1.05646898223576648d-1*u968
  u3889=(u628+u4179)
  u4183=u2980*u3889
  u4070=u4126*u2980
  value_2_=value_2_+(c(37)*((u2885*(u4186+u2980*(u2980*(u473+u3875)+u37&
 &73)+(u4183+u40)*u2885)+u2980*(u855+u2980*(u4070+u3959))+u880)*z))
  u880=1.48588693134999848d7*u658
  u4186=(u3921+1.3d1)
  u855=u3936*(-u4003*u4186+u880)
  u473=-3.94899597530469881d-3*u855
  u3874=2.82353212234285965d-1*u4167*(u965-u4003)
  u4136=-8.39089896167236043d5*u747
  u943=+u4136
  u62=2.43145134220908843d6*u915
  u3826=2.38914256505934278d-1*u4167*(u62-u4003)
  u777=3.26895124896999666d7*u915
  u3816=-2.82353212234285965d-1*u4167*(u777+u4008)
  u4097=1.31764832376000117d0*u968
  u766=-4.51765139574857544d-1*u4167*(u965+u3982)
  u474=5.6470642446857193d-2*u4167*(4.0118947146449959d8*u915+u3862*(6.&
 &3d1+u3992))
  u929=-1.31764832376000117d-1*u801
  u4098=1.31764832376000117d-1*u968
  u4006=u929*u2980
  u3847=u2980*(u474+u4006)
  value_3_=value_3_+(c(37)*(x*(u2885*(u3874+u2980*(u3847+u3816)+(u2980*&
 &(u4097+u4006)+u943)*u2885)+u2980*(u3826+u2980*(u4098*u2980+u766))+u47&
 &3)))
  u4100=1.27361736972857013d6*u658
  u541=2.4d1*pd2
  u430=-1.34033367984784152d-2*u3936*(u3987*(u541+1.3d1)+u4100)
  u3833=u4167*(u547-u3861)
  u793=1.22863920652718806d-2*u3833
  u4011=-6.78087034977107038d4*u747
  u645=+u4011
  u426=1.01362734538493015d-1*u959
  u499=-5.08565276232830279d4*u747
  u997=+u499
  u803=3.85685957957040288d3*exppd2
  u787=9.21479404895391042d-3*u4167*(5.13692339123856618d7*u915+u803)
  u663=-1.59723096848534447d-1*u3827
  u454=3.19446193697068895d-1*u968
  u646=-1.37312624582864175d6*u747
  u3963=2.39584645272801671d-1*u968
  u452=2.08024170388999788d7*u915
  u3892=u4167*(u452+u3872)
  u509=-7.98615484242672237d-2*u3892
  u872=u3989*(2.5d1+u3774)
  u808=u4167*(1.59202171216071266d8*u915-4.53748185831812103d2*u872)
  u416=5.32410322828448158d-2*u808
  u4056=-5.32410322828448158d-2*u801
  u505=8.38546258454805848d-1*u968
  u4051=1.99653871060668059d-1*u968
  u3838=u4056*u2980
  u3793=u2980*(u416+u3838)
  u900=u2980*(u505+u3778)
  u4014=(z*(u2886*(u793+u2980*(u3793+u663)+(u2980*(u454+u3838)+u645)*u2&
 &886)+u2885*(u426+u2980*(u900+u646)+(u2980*(u3963+u3778)+u997)*u2885)+&
 &u2980*(u787+u2980*(u4051*u2980+u509))+u430))
  value_4_=value_4_+(c(37)*u4014)
  u875=1.27155217588623177d-2*u3936*(u4100+u3987*(1.3d1+u541))
  u4100=3.49676848368713738d-2*u959
  u541=2.14429948235792004d4*u747
  u4085=1.65098547927777609d6*u915
  u835=u4167*(-u4003+u4085)
  u500=-5.24515272553070607d-2*u835
  u572=7.23701075295798013d4*u747
  u540=-1.92322266602792556d-1*u4167*(-u4003+u62)
  u62=-7.07618829178113613d5*u747
  u3925=-1.01017756195406191d-1*u968
  u469=3.78816585732773216d-2*u4167*(7.0048955335071357d7*u915+u4023)
  u3926=-3.40934927159495895d-1*u968
  u486=1.89408292866386608d-2*u4167
  u4099=u486*(9.66534070582903775d7*u915+u4023)
  u3870=6.56615415270140241d-1*u968
  u4105=9.55213027296427596d7*u915
  u51=u4167*(u4105-4.53748185831812103d2*u4018)
  u475=-1.01017756195406191d-1*u51
  u3791=-3.15680488110644347d-1*u968
  value_5_=value_5_+(c(37)*(y*(u2886*(u4100+u2980*(u2980*(u3870+u3991)+&
 &u62)+(u2980*(u3925+u3902)+u541)*u2886)+u2885*(u500+u2980*(u2980*(u475&
 &+u520)+u469)+(u2980*(u3926+u3877)+u572)*u2885)+u2980*(u540+u2980*(u29&
 &80*(u3791+u458)+u4099))+u875)))
  value_6_=value_6_+(c(37)*u903)
  u903=4.03312167080713874d7*u658
  u941=-u4003*(9.5d1*pd2-2.6d1)
  u667=u3936*(u903+u941)
  u3985=3.94899597530469881d-4*u667
  u789=-6.51584335925275304d-2*u753
  u3957=9.64861643733765248d5*u915
  u3834=u4167*(-u3988+u3957)
  u642=-9.55657026023737112d-2*u3834
  u79=8.39089896167236043d4*u747
  u630=9.12759114972141925d7*u915
  u4043=u4167*(u630+u974)
  u586=-4.34389557283516869d-3*u4043
  u834=1.07882986650073206d6*u747
  u758=u4167*(1.62739997243095072d7*u915+u3872)
  u767=1.69411927340571579d-1*u758
  u4161=-3.95294497128000351d-1*u968
  u3871=u4167*(1.20286084918809401d7*u915+u3872)
  u4062=8.47059636702857895d-2*u3871
  u431=-8.47059636702857895d-1*u968
  u4162=3.18404342432142532d7*u915
  u4038=u3989*(4.0d1+u4024)
  u845=u4167*(u4162-5.67185232289765129d1*u4038)
  u577=-3.01176759716571696d-1*u845
  u3895=3.7647094964571462d-2*u968
  value_7_=value_7_+(c(37)*(y*((u789+u2980*(u2980*(u431+u3995)+u834))*u&
 &2886+u2885*(u642+u2980*(u2980*(u577+u3995)+u767)+(u2980*(u4161+u3788)&
 &+u79)*u2885)+u2980*(u586+u2980*(u2980*(u3895+u4202)+u4062))+u3985)))
  u431=8.12668447873666523d-2*u959
  u789=-6.72768623177468191d5*u747
  u642=+u789
  u586=-u860
  u767=1.05646898223576648d0*u968
  u4062=-4.22587592894306592d-1*u968
  value_8_=value_8_+(c(37)*(u2934*((u642+u2980*(u3779+u767))*u2885+u298&
 &0*(u586+u2980*(u776+u4062))+u431)*z))
  u431=-1.04480612790511985d-2*u904
  u577=1.58573085206108282d4*u747
  u3895=2.87321685173907958d-1*u959
  u3985=-1.11001159644275797d6*u747
  u860=-7.47036381452160691d-2*u968
  u629=-2.37859627809162423d5*u747
  u419=+u629
  u637=9.96048508602880922d-1*u968
  u446=u2980*(u860+u3876)
  u825=(u446+u577)*u2885
  value_9_=value_9_+(c(37)*(y*(u2885*(u95+u2980*(u2980*(u637+u840)+u398&
 &5)+u825)+u2980*(u3895+u2980*(u56+u419))+u431)))
  u56=-2.29857348139126367d-1*u753
  u3971=1.58573085206108282d6*u747
  u65=-8.9644365774259283d-1*u968
  value_1_=value_1_+(c(38)*(u2934*((u4188+u2980*(u463+u539))*u2885+u298&
 &0*(u3971+u2980*(u388+u65))+u56)*z))
  u56=6.36808684864285064d6*u658
  u65=2.2163684942009087d-3*u3936
  u531=1.5d1*pd2
  u613=u65*(u56-u4003*(2.6d1+u531))
  u92=-1.21900267181049978d-1*u753
  u54=u4167*(-u3861+u547)
  u494=-8.12668447873666522d-3*u54
  u999=-u973
  u973=3.03242230887754792d5*u915
  u585=u4167*(-u3975+u973)
  u3796=-1.36528299242775976d0*u585
  u50=-u780
  u780=1.05646898223576648d-1*u3827
  u891=-2.11293796447153296d-1*u968
  u91=u4167*(u850+u3862)
  u598=1.05646898223576648d-1*u91
  u4045=-u671
  u3939=-3.5215632741192216d-2*u808
  u411=u2980*(u891+u776)
  u3903=u891*u2980
  value_2_=value_2_+(c(38)*(y*((u92+u2980*(u2980*(u4045+u3875)+u50))*u2&
 &886+u2885*(u494+u2980*(u2980*(u3939+u776)+u780)+(u411+u999)*u2885)+u2&
 &980*(u3796+u2980*(u3903+u598))+u613)))
  u613=-1.73755822913406748d-1*u753
  u3796=7.19219911000488037d5*u747
  u598=-1.19869985166748006d5*u747
  u746=-u715
  u715=-1.12941284893714386d0*u968
  u527=+u715
  u61=1.8823547482285731d-1*u968
  u3879=1.12941284893714386d-1*u968
  value_3_=value_3_+(c(38)*(u934*((u3796+u2980*(u383+u527))*u2886+(u598&
 &+u2980*(u3817+u61))*u2885+u2980*(u746+u2980*(u3817+u3879))+u613)))
  u828=8.37708549904900948d-4*u3936*(1.23116345740428446d7*u658-u4003*(&
 &u963-2.6d1))
  u406=-2.76443821468617313d-2*u753
  u963=-u4011
  u4011=-9.21479404895391042d-3*u54
  u797=-u499
  u499=-1.84295880979078208d-1*u4167*(1.44343301902571281d6*u915+u3988)
  u3940=-3.19446193697068895d-1*u968
  u952=1.19792322636400836d-1*u3827
  u3967=-2.39584645272801671d-1*u968
  u895=1.19792322636400836d-1*u4167
  u908=u895*(6.93413901296665958d6*u915+u3862)
  u4201=2.79515419484935283d-1*u968
  u560=-3.99307742121336118d-2*u808
  u392=3.99307742121336118d-2*u801
  u961=-7.98615484242672237d-2*u968
  u723=u392*u2980
  u794=u2980*(u560+u723)
  u955=(y*(u2886*(u406+u2980*(u2980*(u4201+u3778)+u797)+(u2980*(u3940+u&
 &994)+u963)*u2886)+u2885*(u4011+u2980*(u794+u952)+(u2980*(u3967+u723)+&
 &u797)*u2885)+u2980*(u499+u2980*(u961*u2980+u908))+u828))
  value_4_=value_4_+(c(38)*u955)
  u876=2.75502971442016884d-2*u3936*(u522-u3975*(u3992-7.0d0))
  u985=2.12269561621428355d6*u915
  u4109=u4167*(u3987+u985)
  u3783=-3.10823865216634434d-2*u4109
  u542=8.57719792943168016d4*u747
  u4122=-2.33117898912475825d-2*u753
  u578=8.04112305884220015d3*u747
  u575=1.3160712820528558d7*u915
  u897=-2.33117898912475825d-2*u4167*(u575-u3819)
  u711=7.07565205404761182d5*u915
  u3772=u4167*(u711+u717)
  u993=2.42442614868974858d0*u3772
  u3900=-4.04071024781624764d-1*u968
  u651=3.37727168471372406d5*u747
  u3966=-3.78816585732773216d-2*u968
  u3961=u486*(5.70297555556237513d7*u915+u4023)
  u938=-1.34690341593874921d-1*u4167*(u4162-1.13437046457953026d2*u4020&
 &)
  u423=6.73451707969374607d-2*u801
  u3883=-2.2728995143966393d-1*u968
  u971=-1.89408292866386608d-1*u968
  u447=1.26272195244257739d-2*u801
  u3904=u423*u2980
  u3905=u447*u2980
  value_5_=value_5_+(c(38)*(z*(u2886*(u3783+u2980*(u938*u2980+u993)+(u2&
 &980*(u3900+u3904)+u542)*u2886)+u2885*(u4122+u2980*(u2980*(u3883+u3905&
 &)+u651)+(u2980*(u3966+u458)+u578)*u2885)+u2980*(u897+u2980*(u2980*(u9&
 &71+u458)+u3961))+u876)))
  u877=4.18854274952450474d-3*u661
  u661=-4.60739702447695521d-2*u753
  u571=-u742
  u742=-9.21479404895391042d-2*u4167*(u508+u3861)
  u508=-4.23804396860691899d5*u747
  u569=+u508
  u3887=u4167*(1.40097910670142714d7*u915+u3872)
  u896=1.99653871060668059d-1*u3887
  u937=-3.99307742121336118d-1*u968
  u774=3.99307742121336118d-2*u4167*(u965+u3862)
  u631=5.19100064757736954d-1*u968
  u852=1.71938344913356967d8*u915
  u820=u4167*(u852-4.53748185831812103d2*u77)
  u77=-3.99307742121336118d-2*u820
  u4195=u77+u723
  u3964=u2980*(u4195)
  value_6_=value_6_+(c(38)*(x*(u2886*(u661+u2980*(u2980*(u631+u3778)+u5&
 &69)+(u2980*(u4121+u994)+u570)*u2886)+u2885*(u568+u2980*(u3964+u896)+(&
 &u2980*(u937+u723)+u571)*u2885)+u2980*(u742+u774*u2980)+u877)))
  u607=5.05201556658999484d7*u658
  u587=3.94899597530469881d-4*u3936
  u3965=1.19d2*pd2
  u656=u587*(u607-u4003*(1.3d2+u3965))
  u817=-3.6488722811815417d-1*u585
  u72=7.19219911000488036d4*u747
  u668=-1.30316867185055061d-2*u753
  u581=1.19869985166748006d4*u747
  u3855=1.24389963110157016d8*u915
  u945=-4.34389557283516869d-3*u4167
  u46=1.3385571482038457d4*exppd2
  u846=+u945*(u46+u3855)
  u3845=1.06134780810714177d7*u915
  u60=u4167*(u3845+u3862)
  u760=3.38823854681143158d-1*u60
  u4160=-3.38823854681143158d-1*u968
  u582=1.43843982200097607d5*u747
  u976=-5.6470642446857193d-2*u968
  u558=u4167*(1.99533387924142653d7*u915+u3872)
  u816=8.47059636702857895d-2*u558
  u809=u3989*(2.0d1+pd2)
  u870=u4167*(u4162-1.13437046457953026d2*u809)
  u550=-4.51765139574857544d-1*u870
  u405=9.41177374114286549d-3*u801
  u4177=-7.52941899291429239d-2*u968
  u888=u2980*(u4160+u3995)
  u3994=u405*u2980
  u4064=u2980*(u976+u3994)
  u3839=u976*u2980
  value_7_=value_7_+(c(38)*(z*(u2886*(u817+u2980*(u2980*(u550+u383)+u76&
 &0)+(u888+u72)*u2886)+u2885*(u668+u2980*(u3839+u582)+(u4064+u581)*u288&
 &5)+u2980*(u846+u2980*(u2980*(u4177+u4202)+u816))+u656)))
  u656=1.57079475599856982d7*u658
  u817=3.7d1*pd2
  u668=(2.6d1+u817)
  u846=3.69394749033484783d-3*u3936*(u656-u4003*u668)
  u760=-2.84433956755783283d-1*u753
  u816=-u789
  u550=7.21716509512856406d6*u915
  u4177=-8.12668447873666523d-2*u4167*(u4047+u550)
  u789=1.56979345408075911d6*u747
  u950=5.28234491117883239d-1*u3887
  u641=-u767
  u3972=u4167*(u575+u3862)
  u575=1.05646898223576648d-1*u3972
  u4156=-7.39528287565036535d-1*u968
  u4123=-1.05646898223576648d-1*u820
  u507=-1.40862530964768864d-1*u968
  value_8_=value_8_+(c(38)*(x*((u760+u2980*(u2980*(u4156+u776)+u789))*u&
 &2886+u2885*(u642+u2980*(u2980*(u4123+u3875)+u950)+(u2980*(u641+u3875)&
 &+u816)*u2885)+u2980*(u4177+u2980*(u507*u2980+u575))+u846)))
  u846=2.0896122558102397d-3*u947
  u760=6.895720444173791d-2*u959
  u4177=-1.60900143697388457d-1*u753
  u950=-1.2368700646076446d6*u747
  u575=6.50149649345043956d5*u747
  u4156=1.04585093403302497d0*u968
  u4123=-2.73913339865792253d-1*u968
  value_9_=value_9_+(c(38)*((u2885*(u760+u2980*(u2980*(u4156+u840)+u950&
 &)+u825)+u2980*(u4177+u2980*(u2980*(u4123+u3876)+u575))+u846)*z))
  u507=-6.34292340824433128d4*u747
  u825=+u507
  u4026=-1.72393011104344775d-1*u753
  u437=2.98814552580864276d-1*u968
  u942=6.34292340824433128d5*u747
  u579=-2.49012127150720231d-1*u968
  u4129=-1.49407276290432138d-1*u968
  u465=u2980*(u437+u463)
  u481=(u465+u825)*u2885
  u4071=u4129*u2980
  value_1_=value_1_+(c(39)*(y*(u2885*(u95+u2980*(u2980*(u579+u388)+u418&
 &8)+u481)+u2980*(u4026+u2980*(u4071+u942)))))
  u4026=-2.43800534362099957d-1*u753
  u579=1.05646898223576648d-1*u968
  u3860=1.12128103862911365d6*u747
  u476=-7.0431265482384432d-1*u968
  u779=-3.16940694670729944d-1*u968
  value_2_=value_2_+(c(39)*(u2934*(u2885*(u586+u2980*(u3875+u476)+(u417&
 &9+u579)*u2885)+u2980*(u3860+u779*u2980)+u4026)*z))
  u476=-2.36939758518281929d-3*u855
  u492=u4167*(u4085-u4003)
  u4085=1.17285180466549555d-1*u492
  u3915=-1.67817979233447209d5*u747
  u931=1.30316867185055061d-2*u4167
  u849=u931*(7.42943465674999241d7*u915-u974)
  u974=-1.69411927340571579d-1*u4167*(3.46706950648332979d7*u915+u4008)
  u659=7.90588994256000702d-1*u968
  u4125=4.95295643783332827d6*u915
  u5=u4167*(u4125+u3873)
  u4084=-6.77647709362286316d-1*u5
  u4155=1.8823547482285731d-2*u4167
  u844=1.2d1*pd2
  u409=u3989*(1.75d2+u844)
  u477=u4155*(1.11441519851249886d9*u915-4.53748185831812103d2*u409)
  u402=3.95294497128000351d-1*u968
  u503=u2980*(u477+u4006)
  u3906=u402*u2980
  value_3_=value_3_+(c(39)*(y*(u2885*(u4085+u2980*(u503+u974)+(u2980*(u&
 &659+u4006)+u3915)*u2885)+u2980*(u849+u2980*(u3906+u4084))+u476)))
  u425=1.65866292881170388d-1*u492
  u940=-3.99307742121336118d-2*u4167
  u609=+u940*(u853+u4008)
  u564=2.22883039702499772d8*u915
  u821=u4167*(u564-4.53748185831812103d2*u754)
  u478=5.32410322828448158d-2*u821
  u487=-9.31718064949784276d-2*u801
  u4072=u487*u2980
  u488=(u2934*(u2885*(u609+u2980*(u4072+u478)+(u4072+u4201)*u2885)+u298&
 &0*(u609+u4201*u2980)+u425)*z)
  value_4_=value_4_+(c(39)*u488)
  u633=2.12269561621428355d6*u658
  u3952=(3.9d1+1.6d2*pd2)
  u878=1.69540290118164237d-2*u3936*(u633-u3981*u3952)
  u4166=u4167*(u3845+u4003)
  u53=1.16558949456237913d-2*u4166
  u3851=u4167*(u985+u3982)
  u4=-2.02035512390812382d-1*u3851
  u4102=5.05088780977030955d-2*u968
  u429=1.33729823821499863d8*u915
  u4087=-5.82794747281189563d-3*u4167*(6.57934869456127549d3*exppd2+u42&
 &9)
  u4108=u486*(u4105+u4023)
  u81=-1.70467463579747947d-1*u968
  u4088=-1.22386896929049808d-1*u4167*(-u4003+2.52701859073128994d6*u91&
 &5)
  u611=-5.05088780977030955d-2*u3827
  u479=1.68362926992343652d-2*u808
  u3931=-1.68362926992343652d-2*u801
  u851=5.15815034740070902d8*u915
  u743=1.26272195244257739d-2*u4167
  u890=u743*(u851-2.90398838932359746d4*u3989)
  u905=-6.31360976221288694d-3*u4167
  u818=+u905*(2.00594735732249795d9*u915+u3862*(3.15d2+u4031))
  u664=u662*(9.97666939620713267d7*u915+u4023)
  u4086=1.01017756195406191d-1*u968
  u3916=+u905*(1.89132179404692664d9*u915+u3862*(2.97d2+u4031))
  u480=1.13644975719831965d-1*u801
  u636=-5.68224878599159824d-2*u968
  u3907=u3931*u2980
  u4073=u480*u2980
  value_5_=value_5_+(c(39)*(x*(u2886*(u53+u2980*(u4086*u2980+u611)+u288&
 &6*((u4102+u3907)*u2886+u2980*(u479+u3907)+u4))+u2885*(u4087+u2980*(u2&
 &980*(u3916+u3877)+u890)+u2885*((u81+u3877)*u2885+u2980*(u818+u4073)+u&
 &4108))+u2980*(u4088+u2980*(u636*u2980+u664))+u878)))
  value_6_=value_6_+(c(39)*u4014)
  u4014=-7.89799195060939762d-4*u855
  u3932=u4167*(u850-u4047)
  u552=-4.34389557283516869d-3*u3932
  u4000=8.47059636702857895d-2*u3827
  u4063=-1.97647248564000175d-1*u968
  u69=3.04072690098461808d-2*u3833
  u866=-1.12941284893714386d-1*u3827
  u847=u4167*(u564-4.53748185831812103d2*u4033)
  u3867=-2.82353212234285965d-2*u847
  u3792=-2.82353212234285965d-2*u3892
  u795=3.12036255583499681d8*u915
  u879=u4167*(u795-4.53748185831812103d2*u3825)
  u956=2.82353212234285965d-2*u879
  u4135=6.58824161880000585d-2*u968
  u640=-6.58824161880000585d-2*u801
  u3789=u640*u2980
  value_7_=value_7_+(c(39)*(x*(u2885*(u552+u2980*(u2980*(u956+u3789)+u8&
 &66)+u2885*((u4063+u3788)*u2885+u3867*u2980+u4000))+u2980*(u69+u2980*(&
 &u4135*u2980+u3792))+u4014)))
  u4014=-2.95515799226787826d-3*u904
  u552=1.05646898223576648d-1*u959
  u4000=-1.34553724635493638d5*u747
  u69=+u4000
  u866=2.43800534362099957d-2*u959
  u3867=-9.41876072448455468d5*u747
  u3792=+u3867
  u956=6.33881389341459887d-1*u968
  u4135=8.97024830903290922d4*u747
  u3920=-3.5215632741192216d-2*u968
  u490=u2980*(u956+u3779)
  value_8_=value_8_+(c(39)*((u2885*(u552+u2980*(u2980*(u579+u776)+u3792&
 &)+(u490+u69)*u2885)+u2980*(u866+u2980*(u3920*u2980+u4135))+u4014)*z))
  u4014=-3.73518190726080346d-2*u968
  u552=5.60277286089120519d-1*u968
  u866=1.86759095363040173d-1*u968
  u3792=-1.24506063575360115d-2*u968
  u4074=u3792*u2980
  value_9_=value_9_+(c(39)*(x*(u2885*(u3895+u2980*(u2980*(u866+u3876)+u&
 &3985)+u2885*((u4014+u3876)*u2885+u2980*(u552+u840)+u419))+u2980*(u95+&
 &u2980*(u4074+u577))+u431)))
  u563=1.14928674069563183d-2*u959
  u829=-1.03435806662606865d-1*u753
  u4157=1.90287702247329938d5*u747
  u906=2.53716936329773251d5*u747
  u4030=-4.48221828871296415d-1*u968
  value_1_=value_1_+(c(40)*((u2885*(u563+u2980*(u2980*(u4030+u388)+u415&
 &7)+u481)+u2980*(u829+u2980*(u4069+u906))+u846)*z))
  u563=(u962-1.3d1)
  u829=-u3975*u563
  u962=u3936*(u829+u633)
  u481=-1.77309479536072696d-2*u962
  u470=u4167*(u985+u3987)
  u495=1.95040427489679965d-1*u470
  u677=u4167*(u985+u3862)
  u893=1.05646898223576648d-1*u677
  u4132=2.75950430107856861d7*u915
  u669=u4167*(u4132+u4008)
  u964=-2.11293796447153296d-1*u669
  u467=-u3862*(u3774-5.0d0)
  u782=3.5215632741192216d-2*u4167*(u4162+u467)
  u892=u4167*(9.76439983458570431d6*u915+u3862)
  u89=-1.05646898223576648d-1*u892
  u432=1.05646898223576648d-1*u808
  u3780=-7.04312654823844319d-2*u801
  u4042=u2980*(u432+u3779)
  u3840=u579*u2980
  u4075=u3780*u2980
  value_2_=value_2_+(c(40)*(x*(u2885*(u495+u2980*(u4042+u964)+u2885*((u&
 &4126+u776)*u2885+u2980*(u782+u4075)+u893))+u2980*(u495+u2980*(u3840+u&
 &89))+u481)))
  u481=-1.57959839012187952d-3*u3936*(-u3975*(1.88d2*pd2-6.5d1)+1.99533&
 &387924142653d7*u658)
  u893=2.60633734370110121d-2*u3833
  u964=-u582
  u782=4.77828513011868556d-2*u959
  u495=-2.39739970333496012d4*u747
  u506=+u495
  u856=u931*(6.66526423491285034d7*u915-u4003)
  u931=-3.38823854681143158d-1*u3827
  u4111=6.77647709362286316d-1*u968
  u858=-6.47297919900439233d5*u747
  u464=u4167*(8.9153215880999909d6*u915+u3862)
  u815=-3.38823854681143158d-1*u464
  u433=1.12941284893714386d-1*u808
  u4158=-1.12941284893714386d-1*u801
  u4072=3.57647402163428889d-1*u968
  u3841=u4158*u2980
  u491=u2980*(u433+u3841)
  u3812=u2980*(u402+u3817)
  value_3_=value_3_+(c(40)*(z*(u2886*(u893+u2980*(u491+u931)+(u2980*(u4&
 &111+u3841)+u964)*u2886)+u2885*(u782+u2980*(u3812+u858)+(u2980*(u3879+&
 &u3817)+u506)*u2885)+u2980*(u856+u2980*(u4072*u2980+u815))+u481)))
  u3928=u3936*(u633+u829)
  u829=6.70166839923920758d-3*u3928
  u960=3.68591761958156417d-2*u4166
  u574=-6.38892387394137789d-1*u3851
  u806=1.59723096848534447d-1*u968
  u4089=u4167*(u4125+u3975)
  u80=-1.10577528587446925d-1*u4089
  u608=1.19792322636400836d-1*u60
  u555=-1.19792322636400836d-1*u968
  u3927=-7.37183523916312834d-2*u4109
  u951=5.73127816377856558d7*u915
  u466=u4167*(u951+u4023)
  u400=7.98615484242672237d-2*u466
  u387=-3.99307742121336118d-2*u821
  u992=3.99307742121336118d-2*u892
  u972=-3.99307742121336118d-2*u827
  u434=7.98615484242672237d-2*u801
  u412=-3.99307742121336118d-2*u968
  u4076=u434*u2980
  u4176=(x*(u2886*(u960+u2980*(u454*u2980+u663)+u2886*((u806+u3838)*u28&
 &86+u3793+u574))+u2885*(u80+u2980*(u2980*(u972+u723)+u400)+u2885*((u55&
 &5+u723)*u2885+u2980*(u387+u4076)+u608))+u2980*(u3927+u2980*(u412*u298&
 &0+u992))+u829))
  value_4_=value_4_+(c(40)*u4176)
  u3793=u4167*(u711+u3988)
  u711=-2.09806109021228243d-1*u3793
  u652=-u3937
  u4103=-5.05088780977030955d-2*u968
  u4104=7.57633171465546432d-2*u3827
  u4200=7.57633171465546432d-2*u3871
  u958=3.36725853984687303d-1*u968
  u562=-1.26272195244257739d-2*u4167
  u526=+u562*(6.04968250621070811d8*u915-4.53748185831812103d2*u899)
  value_5_=value_5_+(c(40)*(u934*(u2886*(u652+u2980*(u3991+u958)+(u3902&
 &+u4103)*u2886)+u2885*(u4104+u2980*(u520+u526)+(u3877+u81)*u2885)+u298&
 &0*(u4200+u2980*(u458+u81))+u711)))
  value_6_=value_6_+(c(40)*u955)
  u3901=5.30673904053570887d7*u915
  u955=+u945*(u3819+u3901)
  u945=3.53782602702380591d6*u915
  u899=u4167*(u945+u3873)
  u616=3.38823854681143158d-1*u899
  u88=6.77647709362286316d-1*u3851
  u648=-5.6470642446857193d-1*u968
  u822=u4167*(4.13925645161785292d8*u915-4.53748185831812103d2*u627)
  u627=-1.8823547482285731d-2*u822
  u489=-1.41176606117142983d-1*u968
  value_7_=value_7_+(c(40)*(u934*((u746+u2980*(u3995+u648))*u2886+u2885&
 &*(u616+u2980*(u3995+u627)+(u3788+u4063)*u2885)+u2980*(u88+u2980*(u420&
 &2+u489))+u955)))
  u955=7.38789498066969565d-4*u667
  u616=-4.06334223936833261d-2*u753
  u88=-2.43800534362099957d-2*u54
  u648=-u4000
  u627=-3.90080854979359931d-1*u4109
  u667=3.16940694670729944d-1*u3827
  u4000=-6.33881389341459887d-1*u968
  u4101=u4167*(9.19834767026189537d6*u915+u3862)
  u678=3.16940694670729944d-1*u4101
  u4036=-5.28234491117883239d-1*u968
  u462=-1.05646898223576648d-1*u808
  u4133=-3.5215632741192216d-1*u968
  u502=u2980*(u4000+u3875)
  u655=(u502+u648)
  u3908=u4133*u2980
  value_8_=value_8_+(c(40)*(y*((u616+u2980*(u2980*(u4036+u776)+u816))*u&
 &2886+u2885*(u88+u2980*(u2980*(u462+u3875)+u667)+u655*u2885)+u2980*(u6&
 &27+u2980*(u3908+u678))+u955)))
  u955=1.14928674069563183d-1*u959
  u627=6.22530317876800577d-1*u968
  u4131=(u3876+u4014)
  value_9_=value_9_+(c(40)*(u2934*(u2885*(u4188+u2980*(u840+u627)+u4131&
 &*u2885)+u2980*(u4188+u2980*u4131)+u955)*z))
  u3835=(1.3d1+u3921)
  u3921=u3936*(u880-u4003*u3835)
  u880=1.04480612790511985d-3*u3921
  u4131=-1.14928674069563183d-2*u54
  u606=-u507
  u507=-6.895720444173791d-2*u54
  u957=1.49407276290432138d-1*u3827
  u902=-2.98814552580864276d-1*u968
  u3849=1.49407276290432138d-1*u60
  u393=-4.98024254301440461d-2*u808
  u674=-1.99209701720576184d-1*u968
  u461=u2980*(u902+u388)
  value_1_=value_1_+(c(41)*(y*((u4130)*u2886+u2885*(u4131+u2980*(u2980*&
 &(u393+u388)+u957)+(u461+u606)*u2885)+u2980*(u507+u2980*(u674*u2980+u3&
 &849))+u880)))
  u880=-7.3140160308629987d-2*u4167*(u4047+u945)
  u507=1.91042605459285519d7*u915
  u4130=u4167*(u507+u3862)
  u632=1.05646898223576648d-1*u4130
  u781=u4167*(u4162-2.26874092915906051d2*u4017)
  u3805=-1.40862530964768864d-1*u781
  u3909=u4062*u2980
  value_2_=value_2_+(c(41)*(u934*((u816+u2980*(u3875+u641))*u2886+u2885&
 &*(u579+u2980*(u776+u3805)+(u776+u4126)*u2885)+u2980*(u632+u3909)+u880&
 &)))
  u880=u587*(1.82551822994428385d7*u658-u4003*(6.5d1+u4010))
  u632=-1.34660762757890229d-1*u753
  u3805=-4.34389557283516869d-3*u54
  u587=-u495
  u495=u4167*(u832+1.41513041080952236d5*u915)
  u922=-1.25104192497652858d0*u495
  u573=1.36651783090092727d6*u747
  u600=-6.77647709362286316d-1*u968
  u535=5.6470642446857193d-2*u3827
  u657=-1.12941284893714386d-1*u968
  u3938=5.6470642446857193d-2*u464
  u4001=-1.8823547482285731d-2*u808
  u401=1.8823547482285731d-2*u801
  u441=-3.7647094964571462d-2*u968
  u427=u401*u2980
  u712=u2980*(u4001+u427)
  value_3_=value_3_+(c(41)*(y*(u2886*(u632+u2980*(u2980*(u4161+u3817)+u&
 &573)+(u2980*(u600+u383)+u582)*u2886)+u2885*(u3805+u2980*(u712+u535)+(&
 &u2980*(u657+u427)+u587)*u2885)+u2980*(u922+u2980*(u441*u2980+u3938))+&
 &u880)))
  u4159=-9.21479404895391042d-3*u3932
  u3934=-u508
  u508=-1.59723096848534447d-1*u968
  u771=1.19792322636400836d-1*u968
  u421=u895*(u4125+u3862)
  u4125=-1.3310258070711204d-1*u968
  u435=-1.59723096848534447d-1*u781
  u3784=u723+u435
  u4077=u421*u2980
  u86=(u934*(u2886*(u3934+u2980*(u3778+u4125)+(u994+u508)*u2886)+u2885*&
 &(u771+u2980*(u3784)+(u723+u555)*u2885)+u4077+u4159))
  value_4_=value_4_+(c(41)*u86)
  u4069=(u4021-2.8d1)
  u881=1.72189357151260553d-3*u3936*(u633-u4003*u4069)
  u633=-1.45698686820297391d-3*u4167*(2.18637648470071205d8*u915-2.4729&
 &2761278337596d4*exppd2)
  u673=6.06106537172437146d-1*u3851
  u612=-2.02035512390812382d-1*u968
  u545=4.37096060460892173d-3*u4167*(2.90101734215952085d7*u915+u4003)
  u4169=u4167*(u547+u3873)
  u4127=-3.78816585732773216d-2*u4169
  u551=1.89408292866386608d-2*u968
  u681=6.11934484645249042d-2*u4167*(1.68467906048752662d5*u915-u3988)
  u485=8.27851290323570583d7*u915
  u3902=-7.25997097330899365d3*u3989
  u3935=u486*(u485+u3902)
  u482=-2.02035512390812382d-1*u781
  u4128=-8.20769269087675302d-2*u4167*(1.42057322008186668d7*u915+u3872&
 &)
  u3808=2.86563908188928279d8*u915
  u714=u3989*(4.5d1+u3774)
  u483=u662*(u3808-4.53748185831812103d2*u714)
  u4037=-6.31360976221288694d-3*u801
  u894=-1.26272195244257739d-2*u968
  u4124=-1.32585805006470626d-1*u968
  u830=2.52544390488515477d-2*u51
  u438=-1.89408292866386608d-2*u801
  u792=6.31360976221288694d-3*u968
  u4181=-1.26272195244257739d-2*u801
  u3869=u4181*u2980
  u3868=u2980*(u830+u3869)
  u3888=u4037*u2980
  u3910=u438*u2980
  u4078=u792*u2980
  value_5_=value_5_+(c(41)*(x*(u2886*(u633+u2980*(u2980*(u4124+u458)+u3&
 &935)+u2886*((u612+u3904)*u2886+u482*u2980+u673))+u2885*(u545+u2980*(u&
 &3868+u4128)+u2885*((u551+u3888)*u2885+u2980*(u483+u3910)+u4127))+u298&
 &0*(u681+u2980*(u4078+u894))+u881)))
  u831=8.37708549904900948d-4*u3921
  u634=-6.14319603263594028d-3*u4167*(3.5165484401965438d3*exppd2+u850)
  u39=1.18665231120993732d5*u747
  u635=-1.84295880979078208d-2*u3932
  u98=4.79169290545603342d-1*u5
  u5=-5.59030838969870566d-1*u968
  u848=u4167*(u564-4.53748185831812103d2*u3878)
  u644=-2.66205161414224079d-2*u848
  u424=9.31718064949784276d-2*u801
  u3911=u424*u2980
  value_6_=value_6_+(c(41)*(z*(u2886*(u634+u2980*(u644*u2980+u98)+(u298&
 &0*(u5+u3911)+u39)*u2886)+u2980*(u635+u4077)+u831)))
  u3954=1.74061040529571251d7*u658
  u823=4.1d1*pd2
  u4077=(5.2d1+u823)
  u882=9.87248993826174703d-4*u3936*(u3954-u4003*u4077)
  u528=2.52600778329499742d8*u915
  u4184=2.88130098003200685d4*exppd2
  u676=-2.17194778641758435d-3*u4167*(u4184+u528)
  u99=1.69411927340571579d-1*u60
  u3858=-1.69411927340571579d-1*u968
  u946=6.51584335925275303d-3*u4166
  u559=-1.12941284893714386d-1*u3851
  u442=2.82353212234285965d-2*u968
  u672=-3.04072690098461808d-1*u585
  u585=u4167*(4.17463471188809097d7*u915+u4008)
  u883=8.47059636702857895d-2*u585
  u59=-5.6470642446857193d-2*u821
  u773=-2.82353212234285965d-2*u3827
  u484=9.41177374114286549d-3*u808
  u623=4.79479940666992024d4*u747
  u448=-4.8000046079828614d-1*u968
  u96=5.6470642446857193d-2*u968
  u3929=9.41177374114286549d-3*u968
  u3842=u96*u2980
  u4079=u3929*u2980
  value_7_=value_7_+(c(41)*(x*(u2886*(u676+u2980*(u2980*(u448+u4202)+u8&
 &83)+u2886*((u3858+u3995)*u2886+u2980*(u59+u383)+u99))+u2885*(u946+u29&
 &80*(u3842+u773)+u2885*((u442+u4202)*u2885+u2980*(u484+u4202)+u559))+u&
 &2980*(u672+u2980*(u4079+u623))+u882)))
  u882=2.2163684942009087d-3*u3921
  u676=-1.78787058532206635d-1*u3834
  u99=-1.46280320617259974d-1*u54
  u672=2.53552555736583955d0*u3851
  u883=3.16940694670729944d-1*u60
  u59=-2.11293796447153296d-1*u51
  u448=1.40862530964768864d-1*u801
  u3912=u448*u2980
  value_8_=value_8_+(c(41)*(z*(u2886*(u676+u2980*(u2980*(u59+u3912)+u67&
 &2)+u655*u2886)+u2980*(u99+u2980*(u3909+u883))+u882)))
  u676=2.61201531976279962d-4*u3936
  u3832=(5.2d1+3.05d2*pd2)
  u883=u676*(1.29484432589071296d8*u658-u4003*u3832)
  u59=-1.00562589810867785d-1*u753
  u99=-8.61965055521723874d-3*u4167*(-u4003+u853)
  u655=7.47036381452160691d-2*u4167
  u975=u655*(u985-u3873)
  u653=3.73518190726080346d-2*u968
  u3864=-5.74643370347815916d-3*u4167*(4.87779299769198011d3*exppd2+6.1&
 &5581728702142228d7*u915)
  u650=5.55005798221378987d5*u747
  u638=u4167*(1.39154490396269699d7*u915+u3872)
  u94=3.36166371653472311d-1*u638
  u3856=-1.24506063575360115d-2*u4167
  u3893=-u3862*(u3774-3.5d1)
  u457=+u3856*(u3893+u564)
  u584=-1.24506063575360115d-2*u801
  u772=1.14625563275571312d7*u915
  u3881=u655*(u772+u3862)
  u655=-2.61462733508256242d-1*u968
  u4165=u4167*(u485+u3862*(1.3d1+pd2))
  u485=-1.49407276290432138d-1*u4165
  u436=6.22530317876800576d-2*u801
  u3824=-8.71542445027520806d-2*u968
  u413=7.47036381452160691d-2*u801
  u45=u413*u2980
  u3998=u584*u2980
  u3843=u436*u2980
  value_9_=value_9_+(c(41)*(x*((u59+u2980*(u2980*(u655+u3876)+u650))*u2&
 &886+u2885*(u99+u2980*(u2980*(u485+u45)+u94)+u2885*((u653+u3998)*u2885&
 &+u2980*(u457+u3843)+u975))+u2980*(u3864+u2980*(u3824*u2980+u3881))+u8&
 &83)))
  u3914=1.72393011104344775d-1*u959
  u580=-u942
  u942=1.49407276290432138d-1*u968
  u666=-u4112
  u4112=2.49012127150720231d-1*u968
  u3969=(u942+u463)
  u3844=u902*u2980
  u3913=u606*u2980
  value_1_=value_1_+(c(42)*(x*(u2885*(u3914+u2980*(u3844+u666)+u2885*(u&
 &3969*u2885+u2980*(u4112+u388)+u580))+u2980*(u3884+u3913))))
  u3914=2.95515799226787826d-3*u947
  u580=-2.43800534362099957d-2*u753
  u4112=-u4135
  u884=3.5215632741192216d-2*u968
  u4082=-1.05646898223576648d-1*u753
  u4083=-u3867
  u3867=(u4126+u3875)
  u4080=u4000*u2980
  u4081=u648*u2980
  value_2_=value_2_+(c(42)*((u2885*(u580+u2980*(u4080+u4083)+u2885*((u8&
 &84+u4179)*u2885+u2980*u3867+u4112))+u2980*(u4082+u4081)+u3914)*z))
  value_3_=value_3_+(c(42)*(x*(u2885*(u849+u2980*(u659*u2980+u974)+u288&
 &5*((u402+u4006)*u2885+u503+u4084))+u2980*(u4085+u3915*u2980)+u476)))
  u3914=u4056*u2885
  u974=u2885*(u416+u3914)
  u4082=u3963*u2980
  u4083=u997*u2980
  u476=(z*(u2886*(u793+u2885*(u974+u663)+(u2885*(u454+u3914)+u645)*u288&
 &6)+u2885*(u787+u2980*(u4082+u646)+u2885*((u4051+u3778)*u2885+u900+u50&
 &9))+u2980*(u426+u4083)+u430))
  value_4_=value_4_+(c(42)*u476)
  u4084=u3931*u2885
  u4085=u81*u2980
  value_5_=value_5_+(c(42)*(y*(u2886*(u53+u2885*(u4086*u2885+u611)+u288&
 &6*((u4102+u4084)*u2886+u2885*(u479+u4084)+u4))+u2885*(u4088+u2980*(u2&
 &980*(u818+u3877)+u890)+u2885*((u636+u3877)*u2885+u2980*(u3916+u4073)+&
 &u664))+u2980*(u4087+u2980*(u4085+u4108))+u878)))
  value_6_=value_6_+(c(42)*u488)
  u609=7.89799195060939762d-4*u3921
  u3921=-3.04072690098461808d-2*u54
  u4=2.82353212234285965d-2*u3892
  u611=-6.58824161880000585d-2*u968
  u636=4.34389557283516869d-3*u3932
  u4087=1.12941284893714386d-1*u3827
  u477=-2.82353212234285965d-2*u879
  u4088=-8.47059636702857895d-2*u3827
  u480=2.82353212234285965d-2*u847
  u664=1.97647248564000175d-1*u968
  u4086=u664*u2980
  value_7_=value_7_+(c(42)*(y*(u2885*(u3921+u2980*(u2980*(u480+u3789)+u&
 &4087)+u2885*((u611+u3788)*u2885+u477*u2980+u4))+u2980*(u636+u2980*(u4&
 &086+u4088))+u609)))
  u609=2.43800534362099957d-1*u959
  u636=-u3860
  u4=3.16940694670729944d-1*u968
  u4087=7.0431265482384432d-1*u968
  value_8_=value_8_+(c(42)*(u2934*(u2885*(u636+u2980*(u776+u4087)+(u377&
 &9+u4)*u2885)+u2980*(u643+u4070)+u609)*z))
  u636=(u3792+u3876)*u2885
  u4087=u4014*u2980
  value_9_=value_9_+(c(42)*(y*(u2885*(u95+u2980*(u2980*(u552+u3876)+u39&
 &85)+u2885*(u636+u2980*(u866+u840)+u577))+u2980*(u3895+u2980*(u4087+u4&
 &19))+u431)))
  value_1_=value_1_+(c(43)*(u2934*(u2885*(u4188+u396*u4191+(u463+u942)*&
 &u2885)+u2980*(u666+u4071))*z))
  u477=-5.91031598453575652d-3*u962
  u4088=u4039*u4002*pd2
  u480=7.98951473461766613d0*u4088
  u3921=u4167*(2.35855068468253727d5*u915+u3981)
  u611=2.34048512987615958d0*u3921
  u552=u4167*(u4162+u4008)
  u866=-7.04312654823844319d-2*u552
  u849=u4167*(u951+u3862*(9.0d0+u3774))
  u425=-3.5215632741192216d-2*u849
  u53=7.783217259452373d6*u915
  u487=u4167*(u53+u3862)
  u818=-3.16940694670729944d-1*u487
  u3916=3.5215632741192216d-2*u822
  u3915=u4*u2980
  value_2_=value_2_+(c(43)*(y*(u2885*(u480+u2980*(u2980*(u3916+u3779)+u&
 &866)+u2885*((u3920+u776)*u2885+u2980*(u425+u4075)+u884))+u2980*(u611+&
 &u2980*(u3915+u818))+u477)))
  u3916=2.60633734370110121d-2*u3932
  u866=-1.69411927340571579d-1*u3827
  u477=7.52941899291429239d-2*u848
  value_3_=value_3_+(c(43)*(u2934*(u2885*(u866+u2980*(u4006+u477)+(u400&
 &6+u402)*u2885)+u2980*(u866+u3906)+u3916)*z))
  u3916=u555*u2980
  u477=(y*(u2886*(u960+u2885*(u454*u2885+u663)+u2886*((u806+u3914)*u288&
 &6+u974+u574))+u2885*(u3927+u2980*(u2980*(u387+u723)+u400)+u2885*((u41&
 &2+u723)*u2885+u2980*(u972+u4076)+u992))+u2980*(u80+u2980*(u3916+u608)&
 &)+u829))
  value_4_=value_4_+(c(43)*u477)
  u425=2.82567150196940394d-3*u3936*(7.21716509512856406d6*u658+u832*(2&
 &.72d2*pd2-9.1d1))
  u3780=7.0330968803930876d3*exppd2
  u4002=u4167*(u3780+u951)
  u3860=-1.29509943840264347d-3*u4002
  u4108=1.27361736972857013d6*u915
  u974=u4167*(u4108+u3873)
  u488=-3.36725853984687303d-2*u974
  u900=1.68362926992343652d-2*u968
  u503=-1.16558949456237913d-2*u4167*(-u4003+3.353859073618568d7*u915)
  u3788=1.51526634293109286d-1*u3827
  u580=-3.03053268586218573d-1*u968
  u890=u486*(7.08980335815570705d7*u915+u4023)
  u3794=-5.05088780977030955d-2*u808
  u417=5.05088780977030955d-2*u801
  u996=-1.57840244055322173d-1*u968
  u991=-4.66235797824951651d-2*u4167*(-u3975+6.08506076648094617d6*u915&
 &)
  u3897=3.36725853984687303d-2*u552
  u885=1.68362926992343652d-2*u849
  u583=4.34220645177478808d5*u747
  u3914=u486*(6.5237511938318981d7*u915+u4023)
  u486=-1.68362926992343652d-2*u822
  u907=3.36725853984687303d-2*u801
  u3917=u417*u2885
  u497=u2885*(u3794+u3917)
  u3918=u417*u2980
  value_5_=value_5_+(c(43)*(z*(u2886*(u3860+u2980*(u2980*(u486+u3918)+u&
 &3897)+u2885*(u497+u3788)+u2886*((u900+u3907)*u2886+u2885*(u580+u3917)&
 &+u2980*(u885+u907*u2980)+u488))+u2885*(u503+u2980*(u2980*(u81+u458)+u&
 &583)+u2885*((u996+u458)*u2885+u2980*(u81+u3905)+u890))+u2980*(u991+u2&
 &980*(u996*u2980+u3914))+u425)))
  value_6_=value_6_+(c(43)*u4176)
  u425=4.88219991729285216d7*u915
  u583=u4167*(u425+u3861)
  u488=-8.68779114567033738d-3*u583
  u900=1.69411927340571579d-1*u3827
  u503=u4167*(1.76891301351190296d7*u915+u3872)
  u890=8.47059636702857895d-2*u503
  u991=-5.6470642446857193d-2*u808
  u3897=-1.78823701081714444d-1*u968
  u885=8.68779114567033738d-3*u583
  u583=3.38823854681143158d-1*u968
  u3914=-8.47059636702857895d-2*u503
  u486=5.6470642446857193d-2*u808
  u907=-5.6470642446857193d-2*u801
  u996=1.41176606117142983d-1*u968
  u3931=1.78823701081714444d-1*u968
  u3860=u907*u2980
  u4176=u2980*(u486+u3860)
  u3919=u404*u2885
  u813=u2885*(u991+u3919)
  value_7_=value_7_+(c(43)*(z*(u2886*(u2980*(u4176+u866)+u2885*(u813+u9&
 &00)+(u2885*(u4160+u3919)+u2980*(u583+u3860))*u2886)+u2885*(u488+u4191&
 &*(u4202+u996)+u2885*((u3897+u3994)*u2885+u489*u2980+u890))+u2980*(u88&
 &5+u2980*(u3931*u2980+u3914)))))
  u488=5.91031598453575652d-3*u3928
  u890=-u611
  u3897=3.16940694670729944d-1*u487
  u885=-u480
  u3914=7.04312654823844319d-2*u552
  u3931=-3.5215632741192216d-2*u822
  u489=3.5215632741192216d-2*u849
  u487=7.04312654823844319d-2*u801
  u3921=(u779+u3875)
  u4088=u487*u2980
  value_8_=value_8_+(c(43)*(x*(u2885*(u890+u2980*(u2980*(u489+u4179)+u3&
 &914)+u2885*(u3921*u2885+u2980*(u3931+u4088)+u3897))+u2980*(u885+u2980&
 &*(u884*u2980+u3920))+u488)))
  u488=-2.0896122558102397d-3*u904
  u480=4.59714696278252733d-2*u959
  u884=-u577
  u489=-6.66006957865654784d5*u747
  u3914=2.61462733508256242d-1*u968
  value_9_=value_9_+(c(43)*((u2885*(u480+u2980*(u2980*(u3914+u3876)+u48&
 &9)+u2885*(u636+u2980*(u3914+u840)+u884))+u2980*(u480+u2980*(u4074+u88&
 &4))+u488)*z))
  u489=-4.17922451162047939d-3*u962
  u3914=-6.895720444173791d-2*u3793
  u3792=u4167*(u507+u3988)
  u480=7.66191160463754555d-3*u3792
  u884=u4167*(u4105-4.53748185831812103d2*u4022)
  u3931=-4.98024254301440461d-2*u884
  u890=-4.98024254301440461d-2*u464
  u885=u4167*(u429-4.53748185831812103d2*u459)
  u4105=4.98024254301440461d-2*u885
  u807=4.98024254301440461d-2*u968
  value_1_=value_1_+(c(44)*(x*(u2885*(u3914+u2980*(u2980*(u4105+u463)+u&
 &674)+u2885*((u4129+u388)*u2885+u3931*u2980+u942))+u2980*(u480+u2980*(&
 &u807*u2980+u890))+u489)))
  u489=-8.86547397680363479d-3*u962
  u3914=1.62533689574733305d-2*u3833
  u480=-2.43800534362099957d-2*u4166
  u3931=4.22587592894306592d-1*u3851
  u890=7.3140160308629987d-2*u4166
  u4105=-3.16940694670729944d-1*u3827
  u674=-1.26776277868291978d0*u3851
  u611=+u674
  u3920=u398*u2885
  u636=u2885*(u3939+u3920)
  value_2_=value_2_+(c(44)*(z*(u2886*(u3914+u2980*(u4042+u4105)+u2885*(&
 &u636+u780)+(u2885*(u891+u3920)+u490+u4112)*u2886)+u2885*(u480+u2885*(&
 &u4126*u2885+u3931))+u2980*(u890+u2980*(u3915+u611))+u489)))
  u489=-7.89799195060939762d-3*u962
  u3914=7.81901203110330364d-2*u4166
  u480=-1.35529541872457263d0*u3851
  u3931=-5.21267468740220243d-2*u4089
  u890=5.6470642446857193d-2*u60
  u611=7.23982595472528116d-2*u4167*(7.64170421837142077d5*u915-u3988)
  u4076=3.7647094964571462d-2*u466
  u466=-1.8823547482285731d-2*u821
  u4112=u4155*(3.82085210918571038d6*u915+u3862)
  u4155=-1.8823547482285731d-2*u827
  u490=3.7647094964571462d-2*u801
  u849=-1.8823547482285731d-2*u968
  u4089=u490*u2980
  value_3_=value_3_+(c(44)*(x*(u2886*(u3914+u2980*(u4111*u2980+u931)+u2&
 &886*((u583+u3841)*u2886+u491+u480))+u2885*(u3931+u2980*(u2980*(u4155+&
 &u427)+u4076)+u2885*((u976+u427)*u2885+u2980*(u466+u4089)+u890))+u2980&
 &*(u611+u2980*(u849*u2980+u4112))+u489)))
  u491=-1.1169447332065346d-3*u962
  u588=2.04773201087864676d-3*u4167*(u951-u3819)
  u450=-2.66205161414224079d-2*u4167*(1.65570258064714117d7*u915+u3872)
  u3848=5.32410322828448158d-2*u968
  u3907=-2.76443821468617313d-2*u4166
  u3930=4.79169290545603342d-1*u3851
  u953=2.76443821468617313d-2*u492
  u4067=-1.46412838777823243d-1*u4167*(1.33150906835259604d7*u915+u3872&
 &)
  u403=u3989*(4.5d1+u3992)
  u492=2.66205161414224079d-2*u4167*(u3808-4.53748185831812103d2*u403)
  u3808=9.58338581091206684d-1*u3772
  u3772=-1.3310258070711204d-2*u4167*(u564+u3862*(3.5d1+u3960))
  u4074=-1.3310258070711204d-2*u801
  u4070=u392*u2885
  u948=u2885*(u560+u4070)
  u3885=(z*(u2886*(u588+u2980*(u2980*(u3772+u723)+u4067)+u2885*(u948+u9&
 &52)+u2886*((u3848+u3838)*u2886+u2885*(u3967+u4070)+u2980*(u492+u4074*&
 &u2980)+u450))+u2885*(u3907+u2885*(u555*u2885+u3930))+u2980*(u953+u298&
 &0*(u3916+u3808))+u491))
  value_4_=value_4_+(c(44)*u3885)
  u3933=4.41511172182719366d-5*u3936
  u886=u3933*(1.58353092969585553d8*u658-u4003*(3.73d2*pd2-3.64d2))
  u543=-3.93386454414802955d-2*u4167*(9.38703172503649835d6*u915+u4003)
  u4106=3.03053268586218573d-1*u4167*(4.38690427350951933d6*u915+u3873)
  u4107=-4.85662289400991303d-4*u4002
  u449=-1.26272195244257739d-2*u974
  u3838=3.10823865216634434d-2*u4167*(4.6699303556714238d6*u915+1.45341&
 &215774252314d2*exppd2)
  u647=-3.29686045412530206d5*u747
  u811=-6.94497073843417563d-2*u4167*(1.44729246560064787d7*u915+u3872)
  u493=u662*(2.48355387097071175d8*u915+u3862*(3.9d1+u3774))
  u810=-3.03053268586218573d-1*u4167*(u4108+u717)
  u4108=2.08349122153025269d-1*u968
  u879=3.15680488110644347d-2*u968
  value_5_=value_5_+(c(44)*(y*(u2886*(u543+u2980*(u2980*(u4108+u458)+u6&
 &47)+u2885*(u580*u2885+u3788)+u2886*((u917+u3917)*u2886+u497+u2980*(u4&
 &102+u3991)+u4106))+u2885*(u4107+u2980*(u3868+u811)+u2885*((u792+u3888&
 &)*u2885+u2980*(u493+u3910)+u449))+u2980*(u3838+u2980*(u879*u2980+u810&
 &))+u886)))
  value_6_=value_6_+(c(44)*u86)
  u86=2.82353212234285965d-2*u959
  u3868=1.57079475599856982d7*u915
  u3882=+u887*(-u4003+u3868)
  u887=1.69411927340571579d-1*u464
  u497=-7.23982595472528115d-4*u4002
  u917=-1.8823547482285731d-2*u974
  u4110=u4167*(2.22377635984353514d5*u915+u3975)
  u649=1.82443614059077085d-1*u4110
  u675=6.83258915450463634d5*u747
  u665=-4.70588687057143275d-2*u3887
  u839=9.41177374114286549d-3*u820
  u93=-u623
  u4015=9.4117737411428655d-2*u968
  u4090=u4015*u2980
  value_7_=value_7_+(c(44)*(y*(u2886*(u3882+u2980*(u2980*(u4063+u4202)+&
 &u675)+u2885*(u4160*u2885+u900)+u2886*((u3858+u3919)*u2886+u813+u888+u&
 &887))+u2885*(u497+u2980*(u4090+u665)+u2885*((u3929+u4202)*u2885+u2980&
 &*(u839+u4202)+u917))+u2980*(u649+u2980*(u4079+u93))+u86)))
  u887=-2.43800534362099957d-2*u3932
  u497=-4.22587592894306592d-1*u781
  value_8_=value_8_+(c(44)*(u934*((u586+u2980*(u776+u4133))*u2886+u2885&
 &*(u4+u2980*(u3875+u497)+(u3875+u779)*u2885)+u2980*(u3897+u3909)+u887)&
 &))
  u887=u676*(5.73127816377856558d7*u658-u4003*(1.35d2*pd2-5.2d1))
  u497=-1.43660842586953979d-2*u753
  u917=-9.57738950579693194d-4*u4002
  u813=-2.4901212715072023d-2*u974
  u665=1.24506063575360115d-2*u968
  u839=-1.3791440888347582d-1*u4167*(u945-u832)
  u86=-u629
  u3882=1.6185788264796815d-1*u4167*(1.51854378698406438d7*u915+u3872)
  u3919=+u3856*(-u3862*(u3774-9.0d0)+u951)
  u4079=8.9644365774259283d-1*u3851
  u3929=-1.86759095363040173d-1*u968
  u900=-4.98024254301440461d-2*u848
  u3897=-2.36561520793184219d-1*u968
  value_9_=value_9_+(c(44)*(y*((u497+u2980*(u2980*(u3929+u3876)+u86))*u&
 &2886+u2885*(u917+u2980*(u2980*(u900+u45)+u3882)+u2885*((u665+u3998)*u&
 &2885+u2980*(u3919+u3843)+u813))+u2980*(u839+u2980*(u3897*u2980+u4079)&
 &)+u887)))
  u4002=-3.4478602220868955d-2*u54
  u991=-4.98024254301440461d-1*u968
  u3856=-1.99209701720576184d-1*u781
  value_1_=value_1_+(c(45)*(u934*((u666+u2980*(u388+u991))*u2886+u2885*&
 &(u942+u2980*(u388+u3856)+(u388+u4129)*u2885)+u2980*(u3849+u3844)+u400&
 &2)))
  u4002=2.88127904246118131d-2*u3936*(u522+2.26874092915906051d2*u3880)
  u3856=-1.21900267181049978d-1*u3795
  u3909=1.05646898223576648d-1*u892
  u3795=-1.46280320617259974d-1*u753
  u676=1.61464469562592366d6*u747
  value_2_=value_2_+(c(45)*(y*(u2886*(u3856+u2980*(u4080+u676)+u2885*(u&
 &891*u2885+u780)+u2886*((u4126+u3920)*u2886+u636+u502+u3909))+u2885*(u&
 &494+u999*u2885)+u2980*(u3795+u4081)+u4002)))
  u4002=-1.43348553903560567d-1*u4167*(-u3861+1.35080630122727135d6*u91&
 &5)
  u3856=-u501
  u3795=5.6470642446857193d-2*u91
  u780=-9.4117737411428655d-1*u968
  u494=-7.52941899291429239d-2*u781
  u4075=u427+u494
  value_3_=value_3_+(c(45)*(u934*(u2886*(u3856+u2980*(u3817+u780)+(u383&
 &+u4160)*u2886)+u2885*(u96+u2980*(u4075)+(u427+u976)*u2885)+u3795*u298&
 &0+u4002)))
  u832=2.79236183301633649d-4*u3936*(4.28784514475285276d7*u658-u4003*(&
 &1.3d1+1.01d2*pd2))
  u501=-1.01362734538493015d-1*u4167*(-u4003+3.43490745169220428d6*u915&
 &)
  u629=u895*(1.08965041632333222d7*u915+u3862)
  u895=-1.73033354919245651d-1*u968
  u3849=-6.14319603263594028d-3*u753
  u649=4.06852220986264223d5*u747
  u4062=(y*(u2886*(u501+u2980*(u4082+u649)+u2885*(u3967*u2885+u952)+u28&
 &86*((u895+u994+u4070)*u2886+u948+u2980*(u5+u3778)+u629))+u2885*(u4011&
 &+u797*u2885)+u2980*(u3849+u4083)+u832))
  value_4_=value_4_+(c(45)*u4062)
  u4063=(8.95d2*pd2-8.32d2)
  u888=u3933*(3.79962515302356755d8*u658-u4003*u4063)
  u3933=-2.99167970271010643d-1*u4167*(u3987+6.89186888381260892d5*u915&
 &)
  u502=u662*(u630+u4023)
  u3920=-6.73451707969374607d-2*u968
  u636=4.37096060460892173d-3*u4166
  u948=-1.89408292866386608d-2*u3827
  u630=3.78816585732773216d-2*u968
  u3939=-7.57633171465546432d-2*u3851
  u994=6.31360976221288694d-3*u808
  u4083=-9.30183780955214811d3*exppd2
  u939=-3.39963602580693912d-3*u4167*(u425+u4083)
  u848=u662*(4.05434862696928157d8*u915-3.0854876636563223d4*u3989)
  u425=-2.52544390488515477d-2*u4167*(3.37508602978071084d8*u915+u3862*&
 &(5.3d1+u3992))
  u759=-1.89408292866386608d-2*u4167
  u4082=+u759*(u3845+u3872)
  u3997=6.31360976221288694d-3*u847
  u3845=u4037*u2885
  u3846=u551*u2980
  u4091=u551*u2885
  value_5_=value_5_+(c(45)*(z*(u2886*(u3933+u2980*(u2980*(u3997+u3888)+&
 &u848)+u2885*(u2885*(u994+u3845)+u948)+u2886*((u3920+u3904)*u2886+u288&
 &5*(u630+u3845)+u2980*(u425+u3869)+u502))+u2885*(u636+u2885*(u4091+u39&
 &39))+u2980*(u939+u2980*(u3846+u4082))+u888)))
  u636=-2.48799439321755581d-1*u835
  u948=1.19792322636400836d-1*u91
  u3939=-2.79515419484935283d-1*u968
  u994=-3.07159801631797014d-3*u3932
  value_6_=value_6_+(c(45)*(x*(u2886*(u636+u952*u2980+u2886*((u3939+u39&
 &11)*u2886+u387*u2980+u948))+u994*u2980+u831)))
  u939=4.88219991729285216d7*u658
  u848=1.15d2*pd2
  u3904=(1.04d2+u848)
  u4082=5.92349396295704821d-4*u3936*(u939-u4003*u3904)
  u3997=1.04450101750226651d6*u915
  u425=-2.73665421088615628d-1*u4167*(-u3988+u3997)
  u901=u4167*(2.25005735318714056d7*u915+u3872)
  u4080=2.82353212234285965d-2*u901
  u605=-2.93212951166373887d-1*u835
  u835=1.41176606117142983d-1*u4167*(3.69349037221285337d7*u915+u4008)
  u64=-5.6470642446857193d-2*u827
  u90=5.39414933250366028d5*u747
  u889=-9.88236242820000877d-1*u968
  u4092=u87*u2885
  u4093=u442*u2980
  value_7_=value_7_+(c(45)*(z*(u2886*(u425+u2980*(u2980*(u889+u4202)+u8&
 &35)+u2885*(u2885*(u484+u4092)+u773)+u2886*((u976+u3995)*u2886+u2885*(&
 &u96+u4092)+u2980*(u64+u383)+u4080))+u2885*(u946+u2885*(u442*u2885+u55&
 &9))+u2980*(u605+u2980*(u4093+u90))+u4082)))
  u4082=5.15511792509183147d6*u915
  u90=-1.7066037405346997d-1*u4167*(-u4003+u4082)
  u4080=2.3349651778357119d7*u915
  u484=3.62998548665449682d3*exppd2
  u835=-8.12668447873666522d-3*u4167*(u484+u4080)
  u64=3.16940694670729944d-1*u503
  u4092=u3989*(8.5d1+u3960)
  u889=u4167*(5.41287382134642304d8*u915-4.53748185831812103d2*u4092)
  u559=-3.5215632741192216d-2*u889
  u773=1.79404966180658184d5*u747
  u946=-8.45175185788613183d-1*u968
  value_8_=value_8_+(c(45)*(x*(u2886*(u90+u2980*(u946*u2980+u64)+u2886*&
 &(u3921*u2886+u2980*(u559+u3912)+u678))+u2980*(u835+u773*u2980)+u882))&
 &)
  u882=4.30982527760861937d-2*u959
  u835=-6.895720444173791d-1*u495
  u64=8.61965055521723874d-3*u4166
  u559=-3.73518190726080346d-2*u3827
  u946=7.47036381452160691d-2*u968
  u678=-1.49407276290432138d-1*u3851
  u495=1.24506063575360115d-2*u808
  u90=-7.75768549969551487d-2*u4167*(-u4003+u53)
  u53=1.12055457217824104d-1*u552
  u552=1.12055457217824104d-1*u4167*(u507+u3872)
  u507=3.50244776675356785d8*u915
  u4197=u3989*(5.5d1+u3992)
  u605=-3.73518190726080346d-2*u4167*(u507-4.53748185831812103d2*u4197)
  u425=8.71542445027520806d-2*u801
  u3921=u425*u2980
  u4094=u584*u2885
  u4095=u655*u2980
  value_9_=value_9_+(c(45)*(z*(u2886*(u835+u2980*(u2980*(u605+u3921)+u5&
 &3)+u2885*(u2885*(u495+u4094)+u559)+(u2885*(u946+u4094)+u2980*(u4030+u&
 &45)+u82)*u2886)+u2885*(u64+u2885*(u653*u2885+u678))+u2980*(u90+u2980*&
 &(u4095+u552))+u882)))
  u496=-2.0896122558102397d-2*u904
  u827=-u4175
  u4175=(u807+u463)*u2885
  u3922=u991*u2980
  u3923=u666*u2980
  value_1_=value_1_+(c(46)*(y*(u2885*(u3895+u2980*(u3922+u949)+u2885*(u&
 &4175+u2980*(u4134+u388)+u827))+u2980*(u95+u3923)+u496)))
  u827=-8.12668447873666523d-2*u753
  u4094=4.22587592894306592d-1*u968
  u4071=u386*u2885
  u847=u4071+u3875
  u3924=u816*u2980
  u4096=u641*u2980
  value_2_=value_2_+(c(46)*(u2934*(u2885*(u643+u4096+u2885*(u847+u4094)&
 &)+u3924+u827)*z))
  value_3_=value_3_+(c(46)*(y*(u2885*(u3826+u2980*(u4097*u2980+u3816)+u&
 &2885*((u4098+u4006)*u2885+u3847+u766))+u2980*(u3874+u943*u2980)+u473)&
 &))
  u3874=u604*u2885
  u3826=u3874+u3778
  u473=u397*u2885
  u4097=u538*u2980
  u4098=u568*u2980
  u827=(u934*((u570+u2885*(u473+u4121))*u2886+u2885*(u556+u4097+u2885*(&
 &u3826+u4046))+u4098+u4019))
  value_4_=value_4_+(c(46)*u827)
  u4006=u399*u2885
  u3917=u520+u4006
  u929=u814*u2885
  u3847=u410*u2885
  value_5_=value_5_+(c(46)*(x*(u2886*(u4100+u2885*(u2885*(u3870+u929)+u&
 &62)+(u2885*(u3925+u3847)+u541)*u2886)+u2885*(u540+u2980*(u3926*u2980+&
 &u469)+u2885*(u2885*(u3791+u3917)+u2980*(u475+u3877)+u4099))+u2980*(u5&
 &00+u572*u2980)+u875)))
  value_6_=value_6_+(c(46)*u476)
  u645=u3936*(u941+u903)
  u3925=-3.94899597530469881d-4*u645
  u3926=6.51584335925275304d-2*u959
  u3791=4.34389557283516869d-3*u4043
  u505=-u834
  u475=-8.47059636702857895d-2*u3871
  u500=8.47059636702857895d-1*u968
  u540=u4167*(u3957-u3988)
  u469=9.55657026023737112d-2*u540
  u4099=-1.69411927340571579d-1*u758
  u4043=3.01176759716571696d-1*u845
  u4051=-u79
  u903=u405*u2885
  u4081=u3860+u903
  u3871=u907*u2885
  value_7_=value_7_+(c(46)*(x*((u3926+u2885*(u2885*(u500+u3871)+u505))*&
 &u2886+u2885*(u3791+u2980*(u3906+u4099)+u2885*(u2885*(u441+u4081)+u298&
 &0*(u4043+u3789)+u475))+u2980*(u469+u4051*u2980)+u3925)))
  u4043=-8.86547397680363479d-3*u904
  u4051=3.16940694670729944d-1*u959
  u3925=-u3959
  u3926=7.3140160308629987d-2*u959
  u3791=-u3773
  u505=-u967
  u469=(u579+u3779)
  u4099=u999*u2980
  value_8_=value_8_+(c(46)*((u2885*(u4051+u2980*(u3903+u3791)+u2885*(u4&
 &69*u2885+u2980*(u505+u776)+u3925))+u2980*(u3926+u4099)+u4043)*z))
  u4043=u394*u2885
  u4051=u840+u4043
  u3925=u860*u2980
  u3926=u577*u2980
  value_9_=value_9_+(c(46)*(x*(u2885*(u3895+u2980*(u3925+u3985)+u2885*(&
 &u2885*(u614+u4051)+u2980*(u637+u3876)+u419))+u2980*(u95+u3926)+u431))&
 &)
  u3791=1.03435806662606865d-1*u959
  u505=-u906
  u475=-1.14928674069563183d-2*u753
  u500=-u4157
  u637=4.48221828871296415d-1*u968
  value_1_=value_1_+(c(47)*((u2885*(u3791+u2980*(u3844+u500)+u2885*(u41&
 &75+u2980*(u637+u388)+u505))+u2980*(u475+u3913)+u488)*z))
  u3791=-7.38789498066969565d-4*u645
  u505=4.06334223936833261d-2*u959
  u807=3.90080854979359931d-1*u470
  u500=-3.16940694670729944d-1*u4101
  u637=5.28234491117883239d-1*u968
  u645=2.43800534362099957d-2*u3833
  u475=u2980*(u956*u2980+u4105)
  u875=u2980*(u645+u69*u2980)
  value_2_=value_2_+(c(47)*(x*((u505+u2885*(u2885*(u637+u4071)+u642))*u&
 &2886+u2885*(u807+u475+u2885*((u63+u3779)*u2885+u4042+u500))+u875+u379&
 &1)))
  u3791=u4158*u2885
  u807=u2885*(u433+u3791)
  u4100=u3879*u2980
  u4101=u506*u2980
  value_3_=value_3_+(c(47)*(z*(u2886*(u893+u2885*(u807+u931)+(u2885*(u4&
 &111+u3791)+u964)*u2886)+u2885*(u856+u2980*(u4100+u858)+u2885*((u4072+&
 &u3817)*u2885+u3812+u815))+u2980*(u782+u4101)+u481)))
  u62=u2980*(u3967*u2980+u952)
  u646=u2980*(u4011+u797*u2980)
  u430=(x*(u2886*(u406+u2885*(u2885*(u4201+u3874)+u797)+(u2885*(u3940+u&
 &473)+u963)*u2886)+u2885*(u499+u62+u2885*((u961+u723)*u2885+u794+u908)&
 &)+u646+u828))
  value_4_=value_4_+(c(47)*u430)
  value_5_=value_5_+(c(47)*(u934*(u2886*(u652+u2885*(u929+u958)+(u3847+&
 &u4103)*u2886)+u2885*(u4200+u2980*(u3877+u526)+u2885*(u4006+u520+u81))&
 &+u2980*(u4104+u4085)+u711)))
  value_6_=value_6_+(c(47)*u477)
  u541=4.34389557283516869d-3*u4167
  u572=u541*(u3901+u3819)
  u806=-6.77647709362286316d-1*u3851
  u400=5.6470642446857193d-1*u968
  u652=-3.38823854681143158d-1*u899
  u402=1.8823547482285731d-2*u822
  value_7_=value_7_+(c(47)*(u934*((u439+u2885*(u3871+u400))*u2886+u2885&
 &*(u806+u2980*(u3789+u402)+u2885*(u903+u3860+u996))+u2980*(u652+u4086)&
 &+u572)))
  u572=1.77309479536072696d-2*u3928
  u806=-1.95040427489679965d-1*u4109
  u400=2.11293796447153296d-1*u669
  u652=-1.05646898223576648d-1*u677
  u402=-3.5215632741192216d-2*u4167*(u467+u4162)
  value_8_=value_8_+(c(47)*(y*(u2885*(u806+u2980*(u2980*(u402+u4179)+u4&
 &00)+u2885*(u3867*u2885+u2980*(u462+u4088)+u3909))+u2980*(u806+u2980*(&
 &u3840+u652))+u572)))
  value_9_=value_9_+(c(47)*(u2934*(u2885*(u4188+u2980*(u3876+u627)+u288&
 &5*(u4043+u840+u4014))+u2980*(u4188+u4087)+u955)*z))
  u572=4.17922451162047939d-3*u3928
  u400=-7.66191160463754555d-3*u3792
  u4014=4.98024254301440461d-2*u464
  u402=6.895720444173791d-2*u3793
  u652=1.99209701720576184d-1*u968
  u806=-4.98024254301440461d-2*u885
  u955=4.98024254301440461d-2*u884
  u3927=u942*u2980
  value_1_=value_1_+(c(48)*(y*(u2885*(u400+u2980*(u2980*(u955+u463)+u65&
 &2)+u2885*((u440+u388)*u2885+u806*u2980+u4014))+u2980*(u402+u2980*(u39&
 &27+u4129))+u572)))
  u4014=2.43800534362099957d-2*u3932
  u402=4.22587592894306592d-1*u781
  value_2_=value_2_+(c(48)*(u934*((u643+u2885*(u4071+u63))*u2886+u2885*&
 &(u818+u2980*(u3779+u402)+(u3779+u4094)*u2885)+u2980*(u779+u3915)+u401&
 &4)))
  value_3_=value_3_+(c(48)*(y*(u2886*(u3914+u2885*(u4111*u2885+u931)+u2&
 &886*((u583+u3791)*u2886+u807+u480))+u2885*(u611+u2980*(u2980*(u466+u4&
 &27)+u4076)+u2885*((u849+u427)*u2885+u2980*(u4155+u4089)+u4112))+u2980&
 &*(u3931+u2980*(u3839+u890))+u489)))
  u4014=(u934*(u2886*(u3934+u2885*(u3874+u4125)+(u473+u508)*u2886)+u288&
 &5*(u421+u2980*(u4070+u3784))+u2980*(u771+u3916)+u4159))
  value_4_=value_4_+(c(48)*u4014)
  u402=u2980*(u830+u3910)
  u806=u3918+u3847
  value_5_=value_5_+(c(48)*(x*(u2886*(u543+u2980*(u580*u2980+u3788)+u28&
 &85*(u2885*(u4108+u4006)+u647)+u2886*((u4117+u806)*u2886+u2885*(u4102+&
 &u929)+u2980*(u3794+u3918)+u4106))+u2885*(u3838+u2980*(u2980*(u493+u38&
 &88)+u811)+u2885*((u879+u3869)*u2885+u402+u810))+u2980*(u4107+u2980*(u&
 &4078+u449))+u886)))
  value_6_=value_6_+(c(48)*u3885)
  u4078=-2.82353212234285965d-2*u753
  u572=2.82353212234285965d-2*u4167
  u543=u572*(u3868-u4003)
  u400=-1.69411927340571579d-1*u464
  u3772=1.69411927340571579d-1*u968
  u955=-1.82443614059077085d-1*u4110
  u434=-u675
  u807=-9.41177374114286549d-3*u968
  u4067=7.23982595472528115d-4*u4167*(u951+u3780)
  u450=4.70588687057143275d-2*u3887
  u3918=-9.4117737411428655d-2*u968
  u449=1.8823547482285731d-2*u974
  u3915=-9.41177374114286549d-3*u820
  u433=u2885*(u583+u3871)
  u4102=u807*u2885
  value_7_=value_7_+(c(48)*(x*(u2886*(u543+u2980*(u583*u2980+u866)+u288&
 &5*(u2885*(u664+u903)+u434)+u2886*((u3772+u3860)*u2886+u433+u4176+u400&
 &))+u2885*(u955+u2980*(u2980*(u3915+u3994)+u450)+u2885*(u4102+u2980*(u&
 &3918+u3994)+u623))+u2980*(u4067+u2980*(u807*u2980+u449))+u4078)))
  u4078=8.86547397680363479d-3*u3928
  u543=-1.62533689574733305d-2*u54
  u400=-7.3140160308629987d-2*u4166
  u955=-u674
  u434=2.43800534362099957d-2*u4166
  u4067=-1.05646898223576648d-1*u3827
  u450=-4.22587592894306592d-1*u3851
  u449=3.5215632741192216d-2*u808
  u3915=u2980*(u449+u4179)
  u3928=u395*u2885
  u958=u2885*(u462+u3928)
  value_8_=value_8_+(c(48)*(z*(u2886*(u543+u2980*(u3915+u4067)+u2885*(u&
 &958+u667)+(u2885*(u4000+u3928)+u4183+u4135)*u2886)+u2885*(u400+u2885*&
 &(u779*u2885+u955))+u2980*(u434+u2980*(u3840+u450))+u4078)))
  value_9_=value_9_+(c(48)*(x*((u497+u2885*(u2885*(u3929+u4043)+u86))*u&
 &2886+u2885*(u839+u2980*(u2980*(u3919+u3998)+u3882)+u2885*((u3897+u45)&
 &*u2885+u2980*(u900+u3843)+u4079))+u2980*(u917+u2980*(u665*u2980+u813)&
 &)+u887)))
  u4078=-3.4478602220868955d-2*u4166
  u543=5.97629105161728553d-1*u3851
  u400=3.4478602220868955d-2*u4166
  u86=-1.49407276290432138d-1*u3827
  u955=-5.97629105161728553d-1*u3851
  u450=4.98024254301440461d-2*u808
  u434=u2980*(u450+u463)
  u3929=u396*u2885
  u497=u2885*(u393+u3929)
  value_1_=value_1_+(c(49)*(z*(u2886*(u2980*(u434+u86)+u2885*(u497+u957&
 &)+(u2885*(u902+u3929)+u465)*u2886)+u2885*(u4078+u2885*(u4129*u2885+u5&
 &43))+u2980*(u400+u2980*(u3927+u955)))))
  u4078=1.3160712820528558d7*u658
  u543=-2.2163684942009087d-3*u3936
  u400=-u4003*(u3775-1.3d1)
  u465=+u543*(u400+u4078)
  u955=u4167*(u777+u4003)
  u492=2.43800534362099957d-2*u955
  u4106=u4167*(8.63229550593808642d6*u915+u3862)
  u493=-3.16940694670729944d-1*u4106
  u4087=4.87601068724199913d-2*u959
  u4107=-5.38214898541974553d5*u747
  u4104=+u4107
  u419=u2885*(u628+u4071)
  value_2_=value_2_+(c(49)*(x*(u2886*(u492+u475+u2885*(u628*u2885+u4104&
 &)+u2886*((u4+u3779)*u2886+u419+u4042+u493))+u2885*(u4087+u40*u2885)+u&
 &875+u465)))
  u465=-1.89551806814625543d-2*u962
  u492=2.08506987496088097d-1*u470
  u470=-1.12941284893714386d-1*u892
  u4087=-1.30316867185055061d-2*u4166
  u4104=2.25882569787428772d-1*u3851
  u962=u4167*(u853-u4003)
  u4105=1.30316867185055061d-2*u962
  u493=-5.08235782021714737d-1*u638
  u972=2.25882569787428772d-1*u4165
  u996=u4167*(-u3873+u985)
  u3931=-1.12941284893714386d-1*u996
  u890=u4167*(u564+u3893)
  u4108=1.8823547482285731d-2*u890
  u475=-9.4117737411428655d-2*u801
  u3848=u401*u2885
  u856=u2885*(u4001+u3848)
  u3930=u976*u2885
  value_3_=value_3_+(c(49)*(z*(u2886*(u492+u2980*(u2980*(u4108+u427)+u4&
 &93)+u2885*(u856+u535)+u2886*((u3879+u3841)*u2886+u2885*(u657+u3848)+u&
 &2980*(u972+u475*u2980)+u470))+u2885*(u4087+u2885*(u3930+u4104))+u2980&
 &*(u4105+u2980*(u3839+u3931))+u465)))
  u465=u723+u473
  u492=(x*(u2886*(u501+u62+u2885*(u3963*u2885+u649)+u2886*((u895+u465)*&
 &u2886+u2885*(u5+u3874)+u794+u629))+u2885*(u3849+u997*u2885)+u646+u832&
 &))
  value_4_=value_4_+(c(49)*u492)
  u470=9.71324578801982605d-4*u4167
  u493=u470*(1.12502867659357028d8*u915+3.35773657515540956d4*exppd2)
  u475=-4.55663640001058008d5*u747
  u4104=-7.57633171465546432d-2*u4169
  u4105=4.41952683354902085d-1*u968
  u4087=5.05088780977030955d-2*u781
  u3931=u410*u2886
  value_5_=value_5_+(c(49)*(u934*(u2886*(u475+u2980*(u458+u4105)+u2885*&
 &(u4006+u4105)+u2886*(u3931+u929+u3991+u4103))+u2885*(u4104+u2980*(u38&
 &69+u4087)+(u3869+u551)*u2885)+u2980*(u4104+u3846)+u493)))
  value_6_=value_6_+(c(49)*u4062)
  u4105=-2.39739970333496012d5*u747
  u475=+u4105
  u4087=-2.82353212234285965d-2*u968
  u4104=-u4105
  u4105=-4.70588687057143275d-1*u968
  u4103=u4087*u2885
  value_7_=value_7_+(c(49)*(u934*(u2886*(u2980*(u4202+u4105)+u2885*(u90&
 &3+u4115)+(u3871+u3995)*u2886)+u2885*(u475+u4103)+u2980*(u4104+u4093))&
 &))
  u475=u65*(u4078+u400)
  u4104=-2.43800534362099957d-2*u955
  u955=3.16940694670729944d-1*u4106
  u4105=-4.87601068724199913d-2*u753
  u493=-u4107
  value_8_=value_8_+(c(49)*(y*(u2886*(u4104+u2980*(u3903+u493)+u2885*(u&
 &4000*u2885+u667)+u2886*((u779+u3928)*u2886+u958+u411+u955))+u2885*(u8&
 &8+u648*u2885)+u2980*(u4105+u4099)+u475)))
  u493=-1.72393011104344775d-2*u3932
  u475=1.58573085206108282d5*u747
  u4105=4.48221828871296415d-1*u899
  u4104=-2.98814552580864276d-1*u781
  value_9_=value_9_+(c(49)*(u934*((u475+u2980*(u3876+u614)+u2885*(u4043&
 &+u614))*u2886+u2885*(u4105+u2980*(u45+u4104)+(u45+u655)*u2885)+u2980*&
 &(u4105+u4095)+u493)))
  u493=3.4478602220868955d-2*u959
  u475=u4167*(-u4003+u772)
  u4104=-3.4478602220868955d-2*u475
  u4105=1.49407276290432138d-1*u464
  u400=-6.895720444173791d-2*u753
  u4078=7.61150808989319753d5*u747
  value_1_=value_1_+(c(50)*(y*(u2886*(u4104+u2980*(u3844+u4078)+u2885*(&
 &u902*u2885+u957)+u2886*((u4129+u3929)*u2886+u497+u461+u4105))+u2885*(&
 &u4131+u606*u2885)+u2980*(u400+u3913)+u493)))
  u493=8.97024830903290922d5*u747
  u4104=u63*u2885
  u4105=u643*u2885
  value_2_=value_2_+(c(50)*(u934*(u2886*(u493+u4096+u4104+(u847+u891)*u&
 &2886)+u4105+u3924+u4026)))
  u400=1.9d1*pd2
  u891=1.18469879259140964d-3*u3936*(8.06624334161427748d6*u658-u4003*(&
 &5.2d1+u400))
  u902=-3.51855541399648664d-1*u4167*(-u4003+8.33354575254496503d5*u915&
 &)
  u955=u4167*(u452+u3862)
  u999=5.6470642446857193d-2*u955
  u957=-1.04253493748044049d-1*u753
  u411=1.86997176860126889d6*u747
  u88=-1.5811779885120014d0*u968
  value_3_=value_3_+(c(50)*(y*(u2886*(u902+u2980*(u4100+u411)+u2885*(u6&
 &57*u2885+u535)+u2886*((u3858+u383+u3848)*u2886+u856+u2980*(u88+u3817)&
 &+u999))+u2885*(u3805+u587*u2885)+u2980*(u957+u4101)+u891)))
  u856=-1.53579900815898507d-1*u753
  u875=1.35617406995421408d6*u747
  u972=-8.78477032666939461d-1*u968
  u4106=u538*u2885
  u4107=u397*u2886
  u4108=u568*u2885
  u4062=(u934*(u2886*(u875+u4097+u4106+u2886*(u4107+u3826+u972))+u4108+&
 &u4098+u856))
  value_4_=value_4_+(c(50)*u4062)
  u497=-8.83022344365438732d-5*u3936*(-u4003*(1.63d2*pd2+4.16d2)+6.9199&
 &8770885856436d7*u658)
  u917=u470*(1.40522449793385571d8*u915+2.2006787012842887d4*exppd2)
  u470=-3.78816585732773216d-2*u4167
  u839=+u470*(1.5566434518904746d6*u915+u3862)
  u3882=-2.31499024614472521d-1*u968
  u648=4.95375535189011129d-2*u959
  u645=-8.68441290354957616d5*u747
  u62=7.19751512892269111d-1*u968
  u958=5.82794747281189563d-3*u4167*(u547+1.92842978978520144d3*exppd2)
  u627=+u470*(3.60858254756428203d7*u915+u3872)
  u470=u3989*(4.0d1+pd2)
  u892=2.02035512390812382d-1*u4167*(u4162-5.67185232289765129d1*u470)
  u960=-6.31360976221288694d-2*u801
  u646=-u578
  u810=u960*u2980+u929+u3931
  value_5_=value_5_+(c(50)*(x*(u2886*(u917+u2980*(u630*u2980+u627)+u288&
 &5*(u3966*u2885+u645)+u2886*(u2886*(u3882+u810)+u2885*(u62+u4006)+u298&
 &0*(u892+u3888)+u839))+u2885*(u648+u578*u2885)+u2980*(u958+u646*u2980)&
 &+u497)))
  u777=4.9d1*pd2
  u3929=(1.04d2+u777)
  u893=1.39618091650816825d-3*u3936*(2.08024170388999788d7*u658-u4003*u&
 &3929)
  u782=-1.53579900815898507d-2*u4167*(5.21810413706583918d3*exppd2+u452&
 &)
  u953=3.99307742121336118d-2*u955
  u955=-9.31718064949784276d-2*u968
  u900=-7.67899504079492536d-2*u4167*(-u3861+4.16048340777999575d6*u915&
 &)
  u899=1.99653871060668059d-1*u3892
  u3826=u3989*(4.9d1+u3774)
  u829=+u940*(u795-4.53748185831812103d2*u3826)
  value_6_=value_6_+(c(50)*(z*(u2886*(u782+u899*u2980+u2886*((u955+u391&
 &1)*u2886+u829*u2980+u953))+u900*u2980+u893)))
  u900=9.47759034073127714d-3*u947
  u899=-3.38823854681143158d-1*u753
  u829=6.71271916933788834d5*u747
  u588=5.21267468740220243d-2*u959
  u4097=-9.34985884300634447d5*u747
  u647=+u4097
  u3784=-4.34389557283516869d-2*u753
  u992=-u4136
  u641=-7.5294189929142924d-1*u968
  u3895=-u581
  u4096=u3995+u3871
  u3932=u3895*u2980
  u4109=u581*u2885
  value_7_=value_7_+(c(50)*(x*(u2886*(u899+u2980*(u3842+u992)+u2885*(u3&
 &930+u647)+u2886*((u657+u4096)*u2886+u2885*(u659+u903)+u2980*(u641+u42&
 &02)+u829))+u2885*(u588+u4109)+u2980*(u3784+u3932)+u900)))
  u900=4.10387819134761486d6*u915
  u899=u4167*(-u4003+u900)
  u588=-1.21900267181049978d-1*u899
  u3780=1.05646898223576648d-1*u60
  u3784=-1.21900267181049978d-1*u54
  u992=5.28234491117883239d-1*u4167*(u3868+u3872)
  u641=-1.05646898223576648d-1*u4167*(1.84674518610642669d8*u915+u3862*&
 &(2.9d1+u3774))
  u647=-1.40862530964768864d0*u968
  value_8_=value_8_+(c(50)*(z*(u2886*(u588+u2980*(u647*u2980+u992)+u288&
 &6*(u3867*u2886+u2980*(u641+u3912)+u3780))+u2980*(u3784+u493*u2980)+u3&
 &811)))
  u588=5.17179033313034325d-2*u959
  u3780=-5.17179033313034325d-2*u475
  u3784=2.24110914435648207d-1*u464
  u992=-2.24110914435648207d-1*u968
  u641=-1.72393011104344775d-2*u753
  u647=-1.26421541476519502d-1*u3834
  u475=2.24110914435648207d-1*u758
  u606=-3.98419403441152369d-1*u845
  u3892=1.11001159644275797d5*u747
  u431=-5.22925467016512484d-1*u968
  u3816=u2885*(u860+u4043)
  value_9_=value_9_+(c(50)*(x*(u2886*(u3780+u2980*(u431*u2980+u475)+u28&
 &85*(u860*u2885+u4157)+u2886*((u992+u45)*u2886+u3816+u2980*(u606+u3921&
 &)+u3784))+u2885*(u641+u577*u2885)+u2980*(u647+u3892*u2980)+u588)))
  u3793=1.09182240366085024d0*u959
  u3794=-u52
  u818=-u530
  u3792=u995*u2885
  value_1_=value_1_+(c(51)*(x*(u2885*(u3793+u928*u2980+u2885*(u2885*(u8&
 &18+u388+u3792)+u4113*u2980+u3794))+u3884*u2980+u496)))
  u3793=-1.47757899613393913d-2*u904
  u3794=7.72035025479983196d-1*u959
  u818=-u670
  u4113=8.803908185298054d-1*u968
  u845=u3875+u4071
  u4110=u50*u2980
  u4111=u4045*u2980
  u4112=u92*u2980
  value_2_=value_2_+(c(51)*((u2885*(u3794+u4110+u2885*(u2885*(u4113+u84&
 &5)+u4111+u818))+u4112+u3793)*z))
  u3793=u408*u2885
  u3794=u639*u2885
  value_3_=value_3_+(c(51)*(x*((u4116+u2885*(u2885*(u4114+u3793)+u4192)&
 &)*u2886+u2885*(u504+u439*u2980+u2885*(u2885*(u841+u3794)+u819*u2980+u&
 &998))+u944*u2980+u471)))
  u818=u3778+u3874
  u4113=u4205*u2980
  u4114=u3980*u2980
  u4115=u554*u2980
  value_4_=value_4_+(c(51)*(z*((u498+u2885*(u2885*(u4058+u473)+u833))*u&
 &2886+u2885*(u566+u4113+u2885*(u2885*(u3962+u818)+u4114+u615))+u4115+u&
 &472)))
  value_5_=value_5_+(c(51)*(y*(u2886*(u680+u2885*(u2885*(u868+u929)+u46&
 &8)+(u2885*(u4117+u3847)+u3937)*u2886)+u2885*(u4120+u2980*(u4118*u2980&
 &+u4172)+u2885*(u2885*(u3898+u3917)+u2980*(u660+u3877)+u784))+u2980*(u&
 &4119+u3786*u2980)+u873)))
  value_6_=value_6_+(c(51)*u827)
  u4117=-5.92349396295704821d-3*u3936*(u765+u909)
  u4118=4.56109035147692713d-1*u959
  u3898=2.17194778641758435d-2*u567
  u3875=-u867
  u426=u572*(u850-u3872)
  u4116=-u47
  u4072=-2.25882569787428772d-1*u968
  u4119=1.30316867185055061d-1*u4167*(u965-u3988)
  u4120=-2.82353212234285965d-1*u3827
  u827=2.25882569787428772d-1*u654
  u4095=-u576
  u3906=6.58824161880000585d-1*u968
  value_7_=value_7_+(c(51)*(y*((u4118+u2885*(u2885*(u4116+u3871)+u3875)&
 &)*u2886+u2885*(u3898+u2980*(u3906*u2980+u4120)+u2885*(u2885*(u4072+u4&
 &081)+u2980*(u827+u3789)+u426))+u2980*(u4119+u4095*u2980)+u4117)))
  u4095=7.3140160308629987d-1*u959
  u4118=-u524
  u4117=-u3850
  u3875=u898*u2885
  u426=u3875+u776
  u4116=u586*u2980
  value_8_=value_8_+(c(51)*(u2934*(u2885*(u4118+u3908+u2885*(u426+u4117&
 &))+u4116+u4095)*z))
  value_9_=value_9_+(c(51)*(y*(u2885*(u4040+u2980*(u614*u2980+u451)+u28&
 &85*(u2885*(u679+u4051)+u2980*(u3828+u3876)+u546))+u2980*(u66+u82*u298&
 &0)+u874)))
  u4095=2.29857348139126367d-1*u959
  u4118=-u3971
  u4117=8.9644365774259283d-1*u968
  u827=u3792+u388
  value_1_=value_1_+(c(52)*(u2934*(u2885*(u4118+u3922+u2885*(u827+u4117&
 &))+u3923+u4095)*z))
  u4095=(u817+2.6d1)
  u4026=-3.69394749033484783d-3*u3936*(-u4003*u4095+u656)
  u4117=2.84433956755783283d-1*u959
  u4118=8.12668447873666523d-2*u4167*(u550+u4047)
  u3898=-u789
  u4072=-1.05646898223576648d-1*u3972
  u4119=7.39528287565036535d-1*u968
  u4120=1.40862530964768864d-1*u968
  u3906=-5.28234491117883239d-1*u3887
  u498=1.05646898223576648d-1*u820
  u874=u2980*(u767*u2980+u3906)
  u3885=u2980*(u498+u3779)
  u4056=u2980*(u816+u642*u2980)
  value_2_=value_2_+(c(52)*(y*((u4117+u2885*(u2885*(u4119+u4071)+u3898)&
 &)*u2886+u2885*(u4118+u874+u2885*((u4120+u3779)*u2885+u3885+u4072))+u4&
 &056+u4026)))
  u3887=u3794+u3817
  u4117=u61*u2980
  u4118=u598*u2980
  value_3_=value_3_+(c(52)*(u934*((u3796+u2885*(u3793+u527))*u2886+u288&
 &5*(u746+u4117+u2885*(u3887+u3879))+u4118+u613)))
  u3898=u2980*(u937*u2980+u896)
  u4072=u2980*(u568+u571*u2980)
  value_4_=value_4_+(c(52)*(y*(u2886*(u661+u2885*(u2885*(u631+u3874)+u5&
 &69)+(u2885*(u4121+u473)+u570)*u2886)+u2885*(u742+u3898+u2885*(u2980*(&
 &u4195+u4070)+u774))+u4072+u877)))
  u4119=u423*u2885
  u4120=u3966*u2980
  u4121=u578*u2980
  value_5_=value_5_+(c(52)*(z*(u2886*(u3783+u2885*(u938*u2885+u993)+(u2&
 &885*(u3900+u4119)+u542)*u2886)+u2885*(u897+u2980*(u4120+u651)+u2885*(&
 &u2885*(u971+u3905+u4006)+u2980*(u3883+u458)+u3961))+u2980*(u4122+u412&
 &1)+u876)))
  value_6_=value_6_+(c(52)*u430)
  u897=-3.94899597530469881d-4*u3936*(-u4003*(u3965+1.3d2)+u607)
  u4026=u4167*(u973-u3975)
  u813=3.6488722811815417d-1*u4026
  u631=-u72
  u430=u541*(u3855+u46)
  u4019=-3.38823854681143158d-1*u60
  u3883=-8.47059636702857895d-2*u558
  u499=4.51765139574857544d-1*u870
  u870=7.52941899291429239d-2*u968
  u993=1.30316867185055061d-2*u959
  value_7_=value_7_+(c(52)*(z*(u2886*(u813+u2885*(u2885*(u499+u3791)+u4&
 &019)+(u433+u631)*u2886)+u2885*(u430+u2980*(u3842+u964)+u2885*(u2885*(&
 &u870+u903)+u2980*(u96+u4202)+u3883))+u2980*(u993+u3932)+u897)))
  u499=+u543*(-u4003*(u531+2.6d1)+u56)
  u870=1.36528299242775976d0*u4026
  u993=-1.05646898223576648d-1*u91
  u430=8.12668447873666522d-3*u3833
  u4019=u2980*(u628*u2980+u4067)
  u3883=u2980*(u430+u40*u2980)
  value_8_=value_8_+(c(52)*(x*((u3811+u2885*(u2885*(u671+u3875)+u553))*&
 &u2886+u2885*(u870+u4019+u2885*(u3889*u2885+u3915+u993))+u3883+u499)))
  value_9_=value_9_+(c(52)*((u2885*(u4177+u2980*(u3925+u950)+u2885*(u28&
 &85*(u4123+u4051)+u2980*(u4156+u3876)+u575))+u2980*(u760+u3926)+u846)*&
 &z))
  u499=-1.04480612790511985d-3*u855
  u870=6.895720444173791d-2*u3833
  u993=-1.49407276290432138d-1*u60
  u4156=1.14928674069563183d-2*u3833
  u846=u2980*(u437*u2980+u86)
  u462=u2980*(u4156+u825*u2980)
  value_1_=value_1_+(c(53)*(x*((u95+u2885*(u2885*(u4134+u3792)+u949))*u&
 &2886+u2885*(u870+u846+u2885*((u652+u463)*u2885+u434+u993))+u462+u499)&
 &))
  u499=-2.2163684942009087d-3*u855
  u870=1.78787058532206635d-1*u540
  u3779=1.46280320617259974d-1*u3833
  u4123=-u672
  u813=-3.16940694670729944d-1*u60
  u897=2.11293796447153296d-1*u51
  u631=-1.40862530964768864d-1*u801
  u3783=u2885*(u956+u3875)
  u4122=u631*u2885
  value_2_=value_2_+(c(53)*(z*(u2886*(u870+u2885*(u2885*(u897+u4122)+u4&
 &123)+(u3783+u69)*u2886)+u2885*(u3779+u2885*(u4094*u2885+u813))+u499))&
 &)
  u870=u2980*(u657*u2980+u535)
  u897=u2980*(u3805+u587*u2980)
  value_3_=value_3_+(c(53)*(x*(u2886*(u632+u2885*(u2885*(u4161+u3794)+u&
 &573)+(u2885*(u600+u3793)+u582)*u2886)+u2885*(u922+u870+u2885*((u441+u&
 &427)*u2885+u712+u3938))+u897+u880)))
  u4123=u424*u2885
  value_4_=value_4_+(c(53)*(z*(u2886*(u634+u2885*(u644*u2885+u98)+(u288&
 &5*(u5+u4123)+u39)*u2886)+u2885*(u635+u421*u2885)+u831)))
  value_5_=value_5_+(c(53)*(y*(u2886*(u633+u2885*(u2885*(u4124+u4006)+u&
 &3935)+u2886*((u612+u4119)*u2886+u482*u2885+u673))+u2885*(u681+u2980*(&
 &u2980*(u483+u3888)+u4128)+u2885*((u792+u3869)*u2885+u402+u894))+u2980&
 &*(u545+u2980*(u3846+u4127))+u881)))
  value_6_=value_6_+(c(53)*u4014)
  u813=(u823+5.2d1)
  u3779=-9.87248993826174703d-4*u3936*(-u4003*u813+u3954)
  u3897=2.17194778641758435d-3*u4167*(u528+u4184)
  u4014=-1.69411927340571579d-1*u60
  u4125=3.04072690098461808d-1*u4026
  u4026=-8.47059636702857895d-2*u585
  u3903=5.6470642446857193d-2*u821
  u830=4.8000046079828614d-1*u968
  u873=-6.51584335925275303d-3*u4166
  u4127=2.82353212234285965d-2*u3827
  u4128=1.12941284893714386d-1*u3851
  u894=-9.41177374114286549d-3*u808
  u4124=u4087*u2980
  value_7_=value_7_+(c(53)*(y*(u2886*(u3897+u2885*(u2885*(u830+u903)+u4&
 &026)+u2886*((u3772+u3871)*u2886+u2885*(u3903+u3791)+u4014))+u2885*(u4&
 &125+u2980*(u2980*(u894+u3994)+u4127)+u2885*(u4102+u4064+u93))+u2980*(&
 &u873+u2980*(u4124+u4128))+u3779)))
  u4125=7.3140160308629987d-2*u4167*(u945+u4047)
  u3897=-1.05646898223576648d-1*u4130
  u807=1.40862530964768864d-1*u781
  value_8_=value_8_+(c(53)*(u934*((u642+u2885*(u3875+u767))*u2886+u2885&
 &*(u3897+u2980*(u4179+u807)+(u4179+u4094)*u2885)+u2980*(u4126+u3840)+u&
 &4125)))
  u4125=u653*u2980
  value_9_=value_9_+(c(53)*(y*((u59+u2885*(u2885*(u655+u4043)+u650))*u2&
 &886+u2885*(u3864+u2980*(u2980*(u457+u3998)+u94)+u2885*((u3824+u45)*u2&
 &885+u2980*(u485+u3843)+u3881))+u2980*(u99+u2980*(u4125+u975))+u883)))
  u3897=3.4478602220868955d-2*u3833
  u807=1.99209701720576184d-1*u781
  value_1_=value_1_+(c(54)*(u934*((u4188+u2885*(u3792+u539))*u2886+u288&
 &5*(u993+u2980*(u463+u807)+(u463+u437)*u2885)+u2980*(u4129+u3927)+u389&
 &7)))
  u3897=1.7066037405346997d-1*u4167*(u4082-u4003)
  u4129=8.12668447873666522d-3*u4167*(u4080+u484)
  u807=-3.16940694670729944d-1*u503
  u3903=3.5215632741192216d-2*u889
  u4026=-u773
  u873=8.45175185788613183d-1*u968
  value_2_=value_2_+(c(54)*(y*(u2886*(u3897+u2885*(u873*u2885+u807)+u28&
 &86*((u4+u3875)*u2886+u2885*(u3903+u4122)+u500))+u2885*(u4129+u4026*u2&
 &885)+u499)))
  value_3_=value_3_+(c(54)*(u934*(u2886*(u3856+u2885*(u3794+u780)+(u379&
 &3+u4160)*u2886)+u2885*(u3795+u2980*(u3848+u4075))+u2980*(u96+u3839)+u&
 &4002)))
  value_4_=value_4_+(c(54)*(y*(u2886*(u636+u952*u2885+u2886*((u3939+u41&
 &23)*u2886+u387*u2885+u948))+u994*u2885+u831)))
  u3897=3.39963602580693912d-3*u4167*(4.57895768640509736d7*u915-u4083)
  u4129=+u759*(5.16522599945475663d7*u915+u4008)
  u807=2.52544390488515477d-2*u968
  u3903=+u759*(u4132+u3872)
  u4026=u662*(7.32329987593927823d8*u915+u3862*(1.15d2+u3992))
  u873=-4.37096060460892173d-3*u4167*(-u4003+6.29733032810237452d7*u915&
 &)
  u4127=u662*(u851-3.35773657515540956d4*u3989)
  u4128=+u562*(6.68649119107499317d8*u915+u3862*(1.05d2+u4016))
  u894=7.57633171465546432d-2*u4167*(u985-u3982)
  u830=+u905*(-u3862*(u3774-5.5d1)+u507)
  u3779=4.41952683354902085d-2*u801
  u4126=u960*u2885
  value_5_=value_5_+(c(54)*(z*(u2886*(u3933+u2980*(u2980*(u830+u3888)+u&
 &4127)+u2885*(u2885*(u4026+u3845)+u4129)+u2886*((u3920+u806)*u2886+u28&
 &85*(u807+u4126)+u2980*(u4128+u3779*u2980)+u502))+u2885*(u3897+u2885*(&
 &u4091+u3903))+u2980*(u873+u2980*(u3846+u894))+u888)))
  value_6_=value_6_+(c(54)*u492)
  u3897=(u848+1.04d2)
  u895=-5.92349396295704821d-4*u3936*(-u4003*u3897+u939)
  u807=2.73665421088615628d-1*u4167*(u3997-u3988)
  u3903=-2.82353212234285965d-2*u901
  u4026=2.93212951166373887d-1*u959
  u873=-1.25863484425085406d6*u747
  u4127=-1.79804977750122009d5*u747
  u4128=4.23529818351428947d-1*u968
  u894=6.51584335925275303d-3*u962
  u3779=-2.54117891010857368d-1*u638
  u4064=1.12941284893714386d-1*u4165
  u830=-5.6470642446857193d-2*u996
  u4123=9.41177374114286549d-3*u890
  u4129=-4.70588687057143275d-2*u801
  value_7_=value_7_+(c(54)*(z*(u2886*(u807+u2980*(u2980*(u4123+u3994)+u&
 &3779)+u2885*(u2885*(u4128+u903)+u873)+u2886*((u96+u3860)*u2886+u433+u&
 &2980*(u4064+u4129*u2980)+u3903))+u2885*(u4026+u2885*(u4103+u4127))+u2&
 &980*(u894+u2980*(u4124+u830))+u895)))
  u4127=-2.88127904246118131d-2*u3936*(2.26874092915906051d2*u4004+u522&
 &)
  u894=1.21900267181049978d-1*u4167
  u4128=u894*(u785-u4003)
  u830=1.46280320617259974d-1*u959
  u3903=-u676
  value_8_=value_8_+(c(54)*(x*(u2886*(u4128+u4019+u2885*(u956*u2885+u39&
 &03)+u2886*((u579+u4179)*u2886+u3783+u3915+u89))+u2885*(u830+u69*u2885&
 &)+u3883+u4127)))
  u4127=u425*u2885
  u4128=u413*u2885
  value_9_=value_9_+(c(54)*(z*(u2886*(u835+u2980*(u2980*(u495+u3998)+u5&
 &59)+u2885*(u2885*(u605+u4127)+u53)+(u2885*(u4030+u4128)+u2980*(u946+u&
 &3998)+u82)*u2886)+u2885*(u90+u2885*(u655*u2885+u552))+u2980*(u64+u298&
 &0*(u4125+u678))+u882)))
  u449=-3.4478602220868955d-2*u753
  u3903=3.4478602220868955d-2*u4167*(u772-u4003)
  u830=-1.49407276290432138d-1*u464
  u4129=-u4078
  value_1_=value_1_+(c(55)*(x*(u2886*(u3903+u846+u2885*(u437*u2885+u412&
 &9)+u2886*(u3969*u2886+u2885*(u437+u3792)+u434+u830))+u2885*(u760+u825&
 &*u2885)+u462+u449)))
  u4129=u894*(u900-u4003)
  u3903=-1.05646898223576648d-1*u60
  u830=-u493
  value_2_=value_2_+(c(55)*(z*(u2886*(u4129+u874+u2885*(u4104+u830)+u28&
 &86*(u469*u2886+u419+u3885+u3903))+u2885*(u609+u4105)+u4056+u92)))
  value_3_=value_3_+(c(55)*(x*(u2886*(u902+u870+u2885*(u3879*u2885+u411&
 &)+u2886*((u3858+u427+u3793)*u2886+u2885*(u88+u3794)+u712+u999))+u2885&
 &*(u957+u506*u2885)+u897+u891)))
  value_4_=value_4_+(c(55)*(z*(u2886*(u782+u3898+u2885*(u4106+u875)+u28&
 &86*((u955+u465)*u2886+u2885*(u972+u3874)+u3964+u953))+u2885*(u856+u41&
 &08)+u4072+u893)))
  value_5_=value_5_+(c(55)*(y*(u2886*(u917+u2980*(u4120+u645)+u2885*(u6&
 &30*u2885+u627)+u2886*(u2886*(u3882+u3991+u4126+u3931)+u2885*(u892+u38&
 &45)+u2980*(u62+u458)+u839))+u2885*(u958+u646*u2885)+u2980*(u648+u4121&
 &)+u497)))
  value_6_=value_6_+(c(55)*u4062)
  u4129=-9.47759034073127714d-3*u904
  u4064=3.38823854681143158d-1*u959
  u3903=-u829
  u449=4.34389557283516869d-2*u959
  u894=7.5294189929142924d-1*u968
  u807=-5.21267468740220243d-2*u753
  u457=-u4097
  u937=-7.90588994256000702d-1*u968
  value_7_=value_7_+(c(55)*(y*(u2886*(u4064+u2980*(u3842+u457)+u2885*(u&
 &3930+u943)+u2886*((u3879+u4096)*u2886+u2885*(u894+u903)+u2980*(u937+u&
 &4202)+u3903))+u2885*(u449+u4109)+u2980*(u807+u3932)+u4129)))
  u4129=u642*u2885
  value_8_=value_8_+(c(55)*(u934*(u2886*(u830+u3908+u767*u2885+(u426+u6&
 &28)*u2886)+u4129+u4116+u609)))
  value_9_=value_9_+(c(55)*(y*(u2886*(u3780+u2980*(u3925+u4157)+u2885*(&
 &u431*u2885+u475)+u2886*((u992+u4128)*u2886+u2885*(u606+u4127)+u446+u3&
 &784))+u2885*(u647+u3892*u2885)+u2980*(u641+u3926)+u588)))
  value_1_=value_1_+(c(56)*(u934*(u2886*(u3922+u539*u2885+(u827)*u2886)&
 &+u4188*u2885+u3923)))
  value_2_=value_2_+(c(56)*(y*(u2886*(u4110+u4129+u2886*((u845)*u2886+u&
 &637*u2885+u4111))+u505*u2885+u4112)))
  u894=-7.81901203110330365d-1*u753
  u4064=4.55505943633642423d6*u747
  u457=-2.25882569787428772d0*u968
  value_3_=value_3_+(c(56)*(u934*(u2886*(u4064+u4117+u61*u2885+u2886*(u&
 &408*u2886+u3887+u457))+u598*u2885+u4118+u894)))
  u894=2.23388946641306919d-2*u947
  u4064=-1.16720724620082865d0*u753
  u457=3.72947869237408871d6*u747
  u937=-1.3310258070711204d0*u968
  u807=(u2886*(u4064+u4113+u4205*u2885+u2886*(u2886*(u937+u818+u4107)+u&
 &3980*u2885+u4114+u457))+u554*u2885+u4115+u894)
  value_4_=value_4_+(c(56)*(y*u807))
  u3903=(2.08d2+2.1d1*pd2)
  u449=1.3245335165481581d-3*u3936*(8.9153215880999909d6*u658-u4003*u39&
 &03)
  u3994=-4.85662289400991303d-3*u4167*(3.92492180744517469d4*exppd2+6.2&
 &4072511166999363d7*u915)
  u4026=u743*(u429-u3862)
  u947=-4.92461561452605181d-1*u968
  u3898=3.64246717050743477d-1*u959
  u830=-2.09069199529897204d6*u747
  u870=1.02280478147848768d0*u968
  u4097=4.02056152942110007d4*u747
  u829=-6.31360976221288694d-2*u968
  u3966=1.45698686820297391d-1*u4167*(1.6132486683228555d6*u915+u453)
  u3779=-6.31360976221288694d-2*u4167*(4.71238426799570947d7*u915+u3872&
 &)
  u873=5.05088780977030955d-2*u4167*(u852+u3873*(5.4d1+pd2))
  u991=-u4097
  u895=6.31360976221288694d-2*u968
  value_5_=value_5_+(c(56)*(z*(u2886*(u3994+u2980*(u895*u2980+u3779)+u2&
 &885*(u829*u2885+u830)+u2886*(u2886*(u947+u810)+u2885*(u870+u4006)+u29&
 &80*(u873+u3888)+u4026))+u2885*(u3898+u4097*u2885)+u2980*(u3966+u991*u&
 &2980)+u449)))
  value_6_=value_6_+(c(56)*(x*u807))
  u830=3.90950601555165182d-1*u959
  u873=-2.27752971816821212d6*u747
  u870=+u873
  u829=-u715
  u449=5.99349925833740031d4*u747
  u3966=-3.90950601555165182d-1*u753
  u3779=-u873
  u873=-u449
  value_7_=value_7_+(c(56)*(z*(u2886*(u2980*(u4090+u3779)+u2885*(u3918*&
 &u2885+u870)+u2886*((u4096)*u2886+u2885*(u829+u903)+u2980*(u527+u4202)&
 &))+u2885*(u830+u449*u2885)+u2980*(u3966+u873*u2980))))
  value_8_=value_8_+(c(56)*(x*(u2886*(u3924+u553*u2885+u2886*((u776+u38&
 &75)*u2886+u671*u2885+u4036*u2980))+u3811*u2885+u616*u2980)))
  u830=8.61965055521723875d-2*u959
  u870=-8.61965055521723875d-2*u899
  u829=7.47036381452160691d-2*u60
  u3966=-8.61965055521723875d-2*u753
  u3779=-1.72393011104344775d-1*u523
  u873=3.73518190726080346d-1*u3827
  u449=-2.98814552580864276d-1*u654
  u895=-8.71542445027520807d-1*u968
  value_9_=value_9_+(c(56)*(z*(u2886*(u870+u2980*(u895*u2980+u873)+u288&
 &5*(u614*u2885+u666)+u2886*((u860+u45)*u2886+u3816+u2980*(u449+u3921)+&
 &u829))+u2885*(u3966+u82*u2885)+u2980*(u3779+u650*u2980)+u830)))
  if ( lmax .eq. 5 ) go to 100
  u830=A15_v(pd,exppd2,erfpd)
  u870=pd2*u830
  u829=6.36808684864285064d6*u870
  u3966=p**11
  u3779=u954*u3966
  u873=u3779*pd2
  u4097=u873*(u829-u4003*u3886)
  u895=1.60900143697388457d-1*u4097
  u449=6.36808684864285064d6*u830
  u4026=u3779*u4027
  u457=u4026*(u449-u3861)
  u991=3.48616978011008323d-1*u457
  u3994=u4026*(-u3861+u449)
  u947=-3.1375528020990749d0*u3994
  u3898=u4027*u830
  u4064=u3779*u3898
  u614=-3.33003478932827392d7*u4064
  u807=+u614
  u431=9.99010436798482176d7*u4064
  u4123=1.08257476426928461d8*u830
  u937=u4026*(u4123+u3872)
  u3906=1.04585093403302497d0*u937
  u430=-1.64348003919475352d0*u937
  u3886=+u430
  u417=2.05689205211164076d9*u830
  u3921=u4026*(u417-9.07496371663624206d2*u3857)
  u3871=-4.98024254301440461d-2*u3921
  u89=4.98024254301440461d-2*u3921
  u4036=u89*u2980
  u3824=u3871*u2980
  value_1_=value_1_+(c(57)*(u2934*((u991+u2980*(u2980*(u3906+u3824)+u80&
 &7))*u2885+u2980*(u947+u2980*(u2980*(u3886+u4036)+u431))+u895)))
  u895=4.87601068724199913d-2*u4097
  u947=3.5215632741192216d-2*u457
  u3886=-2.6411724555894162d0*u3994
  u953=-1.00915293476620229d7*u4064
  u874=+u953
  u475=1.31189881519606297d8*u4064
  u559=5.28234491117883239d-1*u937
  u875=-2.85246625203656949d0*u937
  u500=+u875
  u3883=-3.5215632741192216d-2*u3921
  u88=1.05646898223576648d-1*u3921
  u839=u88*u2980
  u3882=u3883*u2980
  value_2_=value_2_+(c(57)*(y*((u947+u2980*(u2980*(u559+u3882)+u874))*u&
 &2885+u2980*(u3886+u2980*(u2980*(u500+u839)+u475))+u895)*z))
  u3886=u873*(-u4003*(u3774+1.5d1)+u829)
  u500=-6.08145380196923616d-2*u3886
  u627=-7.90588994256000702d-1*u3994
  u955=1.31764832376000117d-1*u457
  u992=1.18588349138400105d0*u457
  u411=7.55180906550512439d7*u4064
  u856=-1.25863484425085406d7*u4064
  u4067=+u856
  u972=-3.77590453275256219d7*u4064
  u4078=+u972
  u493=-2.3717669827680021d0*u937
  u538=3.95294497128000351d-1*u937
  u628=6.21177066915429123d-1*u937
  u96=1.12941284893714386d-1*u3921
  u3895=-1.8823547482285731d-2*u3921
  u3879=u3895*u2980
  u4011=u96*u2980
  value_3_=value_3_+(c(57)*(u2934*((u627+u2980*(u2980*(u493+u4011)+u411&
 &))*u2886+(u955+u2980*(u2980*(u538+u3879)+u4067))*u2885+u2980*(u992+u2&
 &980*(u2980*(u628+u3879)+u4078))+u500)))
  u437=-1.84295880979078208d-2*u3886
  u551=-5.32410322828448158d-2*u3994
  u942=3.99307742121336118d-2*u457
  u956=9.98269355303340296d-1*u457
  u4107=1.52569582869849084d7*u4064
  u492=-1.14427187152386813d7*u4064
  u4014=+u492
  u464=-4.95851144327009522d7*u4064
  u4062=-7.98615484242672237d-1*u937
  u535=5.98961613182004178d-1*u937
  u539=1.07813090372760752d0*u937
  u92=5.32410322828448158d-2*u3921
  u3932=-3.99307742121336118d-2*u3921
  u3805=u3932*u2980
  u4159=u2980*(u539+u3805)
  u771=u92*u2980
  u598=(u2936*((u551+u2980*(u2980*(u4062+u771)+u4107))*u2886+(u942+u298&
 &0*(u2980*(u535+u3805)+u4014))*u2885+u2980*(u956+u2980*(u4159+u464))+u&
 &437))
  value_4_=value_4_+(c(57)*u598)
  u401=4.45766079404999545d7*u3898
  u3892=-7.06417875492350986d-4*u3779*(-u3975*(5.6d1*u4027+1.5d1*u41)+u&
 &401)
  u501=-1.94264915760396521d-2*u3886
  u629=-1.68362926992343652d-2*u3994
  u3783=1.9741069230792837d8*u870
  u776=1.94264915760396521d-3*u873
  u896=u776*(u3783-u4003*(7.5d1+6.2d1*pd2))
  u630=-5.68224878599159824d-2*u3994
  u469=7.57802334988499226d8*u870
  u810=2.38d2*pd2
  u897=u776*(u469-u4003*(6.15d2+u810))
  u993=1.1111953181494681d0*u457
  u4112=4.82467383530532009d6*u4064
  u554=1.40734719355006999d9*u830
  u846=-1.26272195244257739d-2*u4026
  u402=1.86036756191042962d4*exppd2
  u631=+u846*(u402+u554)
  u4129=1.62832741941554553d7*u4064
  u4156=1.93165301075499803d8*u830
  u616=7.71371915914080575d3*exppd2
  u632=-5.68224878599159824d-2*u4026*(u616+u4156)
  u3918=-5.78960860236638411d7*u4064
  u433=-2.52544390488515477d-1*u937
  u633=1.89408292866386608d-1*u4026*(6.17704424318356512d8*u830+u4023)
  u3915=-8.52337317898739736d-1*u937
  u462=6.31360976221288694d-3*u4026
  u634=u462*(4.59139061787149531d9*u830+u3902)
  u635=1.31323083054028048d0*u937
  u98=1.68362926992343652d-2*u3921
  u59=+u846*(1.8078998563297053d10*u830+u3872*(1.67d2+u4016))
  u99=5.68224878599159824d-2*u3921
  u434=-2.96739658824005686d-1*u937
  u3864=-5.05088780977030955d-2*u3921
  u3816=6.31360976221288694d-2*u3921
  u94=6.31360976221288694d-3*u3921
  u3881=u99*u2980
  u4015=u94*u2980
  u4002=u3816*u2980
  u3795=u3864*u2980
  u3933=u98*u2980
  value_5_=value_5_+(c(57)*(u2886*(u501+u2980*(u2980*(u3918+u2980*(u379&
 &5+u635))+u993)+(u2980*(u4112+u2980*(u3933+u433))+u629)*u2886)+u2885*(&
 &u896+u2980*(u2980*(u633+u2980*(u4002+u59))+u631)+(u2980*(u4129+u2980*&
 &(u3881+u3915))+u630)*u2885)+u2980*(u897+u2980*(u2980*(u634+u2980*(u40&
 &15+u434))+u632))+u3892))
  u502=-1.29007116685354746d-1*u3886
  u636=-3.7268722597991371d-1*u3994
  u948=2.79515419484935283d-1*u457
  u994=2.51563877536441755d0*u457
  u419=3.55995693362981196d7*u4064
  u3885=-2.66996770022235897d7*u4064
  u4056=-8.0099031006670769d7*u4064
  u4072=-1.11806167793974113d0*u937
  u552=8.38546258454805848d-1*u937
  u637=1.31771554900040919d0*u937
  value_6_=value_6_+(c(57)*(u2937*((u636+u2980*(u2980*(u4072+u771)+u419&
 &))*u2886+(u948+u2980*(u2980*(u552+u3805)+u3885))*u2885+u2980*(u994+u2&
 &980*(u2980*(u637+u3805)+u4056))+u502)))
  u4019=1.4d1*u4027
  u4006=1.5d1*u3996
  u902=u3779*(-u4003*(u4019+u4006)+u401)
  u426=-3.94899597530469881d-4*u902
  u957=3.04072690098461808d-2*u4097
  u498=4.45766079404999545d7*u870
  u4096=u873*(u498-u3988*(1.5d1+u407))
  u782=8.68779114567033738d-3*u4096
  u4115=-6.58824161880000585d-2*u3994
  u3780=-u3988*(u407-3.75d2)
  u3811=8.68779114567033738d-3*u873*(u498+u3780)
  u4117=-1.5811779885120014d0*u3994
  u4124=3.12036255583499681d8*u830
  u892=u4026*(u3819+u4124)
  u3997=-5.6470642446857193d-2*u892
  u4094=1.8879522663762811d7*u4064
  u4087=4.45766079404999545d7*u830
  u876=2.82353212234285965d-2*u4026
  u806=1.67886828757770478d4*exppd2
  u4082=u876*(u4087+u806)
  u952=u4026*(1.33729823821499863d8*u830+u3872)
  u465=8.47059636702857895d-1*u952
  u4120=-9.88236242820000877d-1*u937
  u832=5.79495903226499408d8*u830
  u3784=-2.82353212234285965d-2*u4026
  u90=+u3784*(-u4008+u832)
  u827=-1.5811779885120014d0*u937
  u3827=+u827
  u3851=3.78901167494249613d9*u830
  u413=u4026*(u3851-9.07496371663624206d2*u754)
  u60=-5.6470642446857193d-2*u413
  u845=6.58824161880000585d-2*u3921
  u3875=3.10588533457714561d-1*u937
  u95=5.6470642446857193d-2*u3921
  u888=-9.41177374114286549d-3*u3921
  u3887=u888*u2980
  u858=u95*u2980
  u772=u845*u2980
  value_7_=value_7_+(c(57)*((u957+u2980*(u2980*(u411+u2980*(u858+u3827)&
 &)+u4117))*u2886+u2885*(u782+u2980*(u2980*(u465+u2980*(u858+u60))+u399&
 &7)+(u2980*(u4094+u2980*(u772+u4120))+u4115)*u2885)+u2980*(u3811+u2980&
 &*(u2980*(u90+u2980*(u3887+u3875))+u4082))+u426))
  u426=1.13773582702313313d-1*u4097
  u957=7.39528287565036535d-1*u457
  u782=-2.21858486269510961d0*u3994
  u3997=-7.06407054336341601d7*u4064
  u3811=+u3997
  u4117=-u3997
  u3997=2.21858486269510961d0*u937
  u4082=-1.16211588045934313d0*u937
  u465=+u4082
  u4120=-1.05646898223576648d-1*u3921
  u90=3.5215632741192216d-2*u3921
  u3827=u4120*u2980
  u3875=u90*u2980
  value_8_=value_8_+(c(57)*(x*((u957+u2980*(u2980*(u3997+u3827)+u3811))&
 &*u2885+u2980*(u782+u2980*(u2980*(u465+u3875)+u4117))+u426)*z))
  u426=6.36808684864285064d6*u3898
  u782=4.0d0*u4027
  u465=u3779*(-u3988*(u782+1.5d1*u966)+u426)
  u4115=-2.0896122558102397d-3*u465
  u503=-3.4478602220868955d-2*u3886
  u604=-1.24506063575360115d-2*u3994
  u898=1.49407276290432138d-1*u4097
  u995=1.86759095363040173d0*u457
  u4083=3.56789441713743635d6*u4064
  u638=-1.60612822012214549d0*u3994
  u3874=-9.2765254845573345d7*u4064
  u818=-1.86759095363040173d-1*u937
  u868=3.68682423104201756d7*u4064
  u639=2.01699822992083387d0*u937
  u87=1.24506063575360115d-2*u3921
  u882=-4.85573647943904449d-1*u937
  u814=-7.47036381452160691d-2*u3921
  u4037=u87*u2980
  u4025=u814*u2980
  value_9_=value_9_+(c(57)*(u2885*(u503+u2980*(u2980*(u3874+u2980*(u402&
 &5+u639))+u995)+(u2980*(u4083+u2980*(u4037+u818))+u604)*u2885)+u2980*(&
 &u898+u2980*(u2980*(u868+u2980*(u4037+u882))+u638))+u4115))
  u4108=-1.39307483720682646d-3*u465
  u465=3.83095580231877278d-3*u4097
  u944=4.98024254301440461d-2*u457
  u486=7.27881602440566827d-2*u4097
  u504=-4.98024254301440461d-1*u3994
  u4181=-1.42715776685497454d7*u4064
  u4116=+u4181
  u566=-5.47826679731584507d-1*u3994
  u4175=3.80575404494659877d7*u4064
  u560=7.47036381452160692d-1*u937
  u3920=7.92865426030541411d6*u4064
  u4080=-1.09565335946316901d0*u937
  u3909=+u4080
  u3928=-4.98024254301440461d-2*u937
  u4130=u3928*u2980
  value_1_=value_1_+(c(58)*(u2885*(u465+u2980*(u2980*(u4175+u2980*(u403&
 &6+u3909))+u504)+(u2980*(u4116+u2980*(u3824+u560))+u944)*u2885)+u2980*&
 &(u486+u2980*(u2980*(u3920+u4130)+u566))+u4108))
  u4108=-5.28234491117883239d-1*u3994
  u465=-3.36384311588734096d6*u4064
  u486=+u465
  u504=-6.69097022082652103d-1*u3994
  u3909=5.71853329700847963d7*u4064
  u545=3.5215632741192216d-1*u937
  u3808=1.34553724635493638d7*u4064
  u4051=-2.00729106624795631d0*u937
  u4176=+u4051
  u893=-1.05646898223576648d-1*u937
  u795=(u2980*(u545+u3882)+u486)
  value_2_=value_2_+(c(58)*(x*(u2885*(u4108+u2980*(u2980*(u4176+u839)+u&
 &3909)+u795*u2885)+u2980*(u504+u2980*(u893*u2980+u3808))+u895)*z))
  u4108=u3779*(u401-u4003*(3.0d0*u4186+u4019))
  u4176=1.31633199176823294d-3*u4108
  u4099=1.4d1*pd2
  u504=-1.8823547482285731d-2*u873*(-u4003*(u4099+1.5d1)+u498)
  u885=-1.44796519094505623d-3*u873
  u505=+u885*(-u4003*(5.18d2*pd2+3.75d2)+1.64933449379849832d9*u870)
  u996=4.89412234539429006d-1*u4026*(7.88663063562691502d7*u830+u4047)
  u917=7.57802334988499226d8*u830
  u759=1.8823547482285731d-2*u4026
  u958=u759*(u917+1.3158697389122551d4*exppd2)
  u588=7.42943465674999241d7*u830
  u971=u4026*(u588+u3862)
  u848=-3.38823854681143158d0*u971
  u640=1.97647248564000175d0*u937
  u3969=-5.6470642446857193d-2*u4026
  u3868=+u3969*(4.6062494871849953d8*u830+u4008)
  u3843=1.28826396948044868d10*u830
  u4101=3.7647094964571462d-2*u4026*(u3843+u3872*(1.19d2+u3960))
  u4001=-1.31764832376000117d-1*u3921
  u641=1.31764832376000117d-1*u937
  u3961=u4001*u2980
  u4177=u2980*(u3961+u4101)
  value_3_=value_3_+(c(58)*(u2885*(u504+u2980*(u2980*(u848+u4177)+u996)&
 &+(u2980*(u4078+u2980*(u3961+u640))+u955)*u2885)+u2980*(u505+u2980*(u2&
 &980*(u3868+u641*u2980)+u958))+u4176))
  u80=1.84674518610642669d8*u870
  u3919=5.8d1*pd2
  u438=-6.14319603263594028d-3*u873*(-u4003*(u3919+4.5d1)+u80)
  u4109=u4026*(u4123+u4047)
  u959=5.32410322828448158d-2*u4109
  u541=-5.08565276232830279d6*u4064
  u3924=+u541
  u960=5.19100064757736954d-1*u457
  u753=-3.81423957174622709d6*u4064
  u3845=+u753
  u436=9.61581114145070446d8*u830
  u77=1.3310258070711204d-2*u4026
  u961=u77*(u436+1.40661937607861752d4*exppd2)
  u491=2.61091560794356876d8*u830
  u938=u4026*(u491+u4008)
  u4088=-2.66205161414224079d-1*u938
  u567=5.32410322828448158d-1*u937
  u398=-4.1956635289208498d7*u4064
  u540=3.99307742121336118d-1*u937
  u962=u4026*(u449+u3863)
  u424=-5.11113909915310231d0*u962
  u394=u4026*(3.35598176923478229d9*u830+u3872*(3.1d1+u3774))
  u3863=5.32410322828448158d-2*u394
  u482=-5.32410322828448158d-2*u3921
  u568=1.99653871060668059d-1*u937
  u4188=u482*u2980
  u879=(u2937*(u2886*(u959+u2980*(u2980*(u3863+u4188)+u4088)+(u2980*(u5&
 &67+u4188)+u3924)*u2886)+u2885*(u960+u2980*(u4159+u398)+(u2980*(u540+u&
 &3805)+u3845)*u2885)+u2980*(u961+u2980*(u568*u2980+u424))+u438))
  value_4_=value_4_+(c(58)*u879)
  u383=u873*(u829+4.53748185831812103d2*u3880)
  u3880=1.74838424184356869d-1*u383
  u997=2.52544390488515477d-1*u457
  u3922=1.60822461176844003d6*u4064
  u895=u4026*(-u3988+u449)
  u642=-1.01017756195406191d0*u895
  u3923=5.4277580647184851d6*u4064
  u643=-5.89270244473202781d-1*u4026*(-u3861+2.27431673165816094d7*u830&
 &)
  u3930=-2.73398184000634805d7*u4064
  u3839=-1.68362926992343652d-1*u937
  u825=1.24177693548535587d9*u830
  u69=6.31360976221288694d-2*u4026
  u644=u69*(u825+u3902)
  u3914=-5.68224878599159824d-1*u937
  u880=1.89408292866386608d-2*u4026
  u645=u880*(2.08660979073864073d9*u830+u3902)
  u646=9.59668683856358814d-1*u937
  u414=-2.52544390488515477d-2*u4026
  u61=+u414*(8.11931073201963457d9*u830+u3872*(7.5d1+u3992))
  u4122=-3.66189366208347442d-1*u937
  u600=u3839+u3933
  value_5_=value_5_+(c(58)*(u2934*(u2886*(u997+u2980*(u2980*(u646+u3795&
 &)+u3930)+(u2980*(u600)+u3922)*u2886)+u2885*(u642+u2980*(u2980*(u61+u4&
 &002)+u644)+(u2980*(u3914+u3881)+u3923)*u2885)+u2980*(u643+u2980*(u298&
 &0*(u4122+u4015)+u645))+u3880)))
  value_6_=value_6_+(c(58)*u598)
  u598=-u3861*(u4005-6.0d0)
  u792=u873*(u498+u598)
  u506=2.17194778641758435d-2*u792
  u949=-3.95294497128000351d-1*u3994
  u883=4.05241890368181404d6*u830
  u847=u4026*(-u3975+u883)
  u553=-1.65647217844114433d0*u847
  u4079=6.29317422125427032d6*u4064
  u3840=3.77186682573461153d7*u830
  u3962=u4026*(u3840+u3861)
  u494=-2.44706117269714503d-1*u3962
  u4110=-u972
  u972=u4026*(2.82318516956499712d8*u830+u4008)
  u454=2.82353212234285965d-1*u972
  u4155=-6.58824161880000585d-1*u937
  u3772=1.63447562448499833d8*u830
  u998=u4026*(u3772+u4008)
  u3917=8.47059636702857895d-2*u998
  u4179=-1.18588349138400105d0*u937
  u451=+u4179
  u3971=5.30461634491949458d9*u830
  u425=u4026*(u3971+u3872*(4.9d1+u4024))
  u62=-3.7647094964571462d-2*u425
  u569=1.12941284893714386d-1*u937
  value_7_=value_7_+(c(58)*(u2934*((u949+u2980*(u2980*(u451+u858)+u4110&
 &))*u2886+u2885*(u553+u2980*(u2980*(u62+u858)+u454)+(u2980*(u4155+u772&
 &)+u4079)*u2885)+u2980*(u494+u2980*(u2980*(u569+u3887)+u3917))+u506)))
  u506=-5.41778965249111015d-3*u3886
  u949=1.05646898223576648d-1*u457
  u553=-3.5215632741192216d-2*u3994
  u494=-3.02745880429860686d7*u4064
  u454=+u494
  u4155=1.68192155794367048d7*u4064
  u3917=1.58470347335364972d0*u937
  u451=-6.69097022082652103d-1*u937
  value_8_=value_8_+(c(58)*(y*((u949+u2980*(u2980*(u3917+u3827)+u454))*&
 &u2885+u2980*(u553+u2980*(u2980*(u451+u3875)+u4155))+u506)*z))
  u451=-4.59714696278252733d-2*u3886
  u950=3.98419403441152369d-1*u457
  u3838=1.18929813904581212d6*u4064
  u4106=-4.20218675796186947d7*u4064
  u855=-1.24506063575360115d-1*u937
  u647=1.44427033747417734d0*u937
  u499=-2.24110914435648207d-1*u937
  u429=u2980*(u855+u4037)
  u3855=(u429+u3838)*u2885
  value_9_=value_9_+(c(58)*(u2934*(u2885*(u950+u2980*(u2980*(u647+u4025&
 &)+u4106)+u3855)+u2980*(u950+u2980*(u2980*(u499+u4037)+u3838))+u451)))
  u899=1.53238232092750911d-2*u4097
  u393=-9.46246083172736875d-1*u3994
  u3812=5.23291181180157331d7*u4064
  u815=-1.24506063575360115d0*u937
  u3912=+u815
  value_1_=value_1_+(c(59)*(y*((u944+u2980*(u2980*(u560+u3824)+u4116))*&
 &u2885+u2980*(u393+u2980*(u2980*(u3912+u4036)+u3812))+u899)*z))
  u393=1.20993650124214162d8*u870
  u3912=8.12668447873666522d-3*u873
  u900=u3912*(u393-u3861*(4.5d1+u400))
  u570=-7.39528287565036535d-1*u3994
  u4165=u4026*(u4047+u4123)
  u605=-3.5215632741192216d-2*u4165
  u4042=-u465
  u465=7.52592082112336894d6*u830
  u4074=-1.5494878406124575d0*u4026
  u648=+u4074*(-u4003+u465)
  u606=1.7607816370596108d-1*u938
  u450=-3.5215632741192216d-1*u937
  u999=u4026*(2.9930008188621398d8*u830+u3872)
  u963=1.05646898223576648d-1*u999
  u822=-u3997
  u51=-3.5215632741192216d-2*u394
  u3833=-2.11293796447153296d-1*u937
  u3934=u3833*u2980
  value_2_=value_2_+(c(59)*(u2934*((u570+u2980*(u2980*(u822+u839)+u4117&
 &))*u2886+u2885*(u605+u2980*(u2980*(u51+u3875)+u606)+(u2980*(u450+u387&
 &5)+u4042)*u2885)+u2980*(u648+u2980*(u3934+u963))+u900)))
  u900=1.15837215275604498d-2*u4097
  u648=-1.12941284893714386d-1*u3994
  u963=1.8823547482285731d-2*u457
  u649=-3.20000307198857427d-1*u3994
  u396=3.23648959950219617d7*u4064
  u659=-5.39414933250366028d6*u4064
  u399=+u659
  u849=1.79804977750122009d6*u4064
  u3788=-1.69411927340571579d0*u937
  u571=2.82353212234285965d-1*u937
  u650=2.44706117269714503d-1*u937
  value_3_=value_3_+(c(59)*(u2936*((u648+u2980*(u2980*(u3788+u4011)+u39&
 &6))*u2886+(u963+u2980*(u2980*(u571+u3879)+u399))*u2885+u2980*(u649+u2&
 &980*(u2980*(u650+u3879)+u849))+u900)))
  u4086=7.0048955335071357d7*u870
  u420=9.21479404895391042d-3*u873
  u4093=1.1d1*pd2
  u833=u420*(u4086-u3861*(u4093-1.5d1))
  u546=-3.99307742121336118d-2*u3994
  u4118=-u541
  u541=-3.99307742121336118d-2*u4165
  u472=-u753
  u753=u4026*(u449+u3975)
  u572=-1.06482064565689632d0*u753
  u410=-6.35706595291037849d6*u4064
  u485=-5.32410322828448158d-1*u937
  u542=1.99653871060668059d-1*u938
  u811=-3.99307742121336118d-1*u937
  u490=1.19792322636400836d-1*u4026
  u543=u490*(1.29484432589071296d8*u830+u3872)
  u573=5.19100064757736954d-1*u937
  u39=-3.99307742121336118d-2*u394
  u86=3.99307742121336118d-2*u3921
  u821=-7.98615484242672237d-2*u937
  u463=u86*u2980
  u797=u2980*(u39+u463)
  u495=(u2934*(u2886*(u546+u2980*(u2980*(u573+u3805)+u410)+(u2980*(u485&
 &+u771)+u4118)*u2886)+u2885*(u541+u2980*(u797+u542)+(u2980*(u811+u463)&
 &+u472)*u2885)+u2980*(u572+u2980*(u821*u2980+u543))+u833))
  value_4_=value_4_+(c(59)*u495)
  u884=-u3988*(u3992-9.0d0)
  u4166=u873*(u829+u884)
  u901=1.16558949456237913d-1*u4166
  u651=-6.73451707969374606d-1*u753
  u3931=6.43289844707376012d6*u4064
  u652=-1.26272195244257739d-1*u3994
  u477=6.03084229413165011d5*u4064
  u653=-2.10453658740429565d-1*u4026*(3.9482138461585674d7*u830+u3861)
  u654=4.04071024781624764d0*u4026*(1.06134780810714177d7*u830+u3982)
  u3876=-6.73451707969374606d-1*u937
  u392=1.08555161294369702d7*u4064
  u4119=-6.31360976221288694d-2*u937
  u655=u880*(1.23753154425292731d9*u830+u3902)
  u801=5.41287382134642304d8*u830
  u929=-1.34690341593874921d-1*u4026
  u63=+u929*(u801-2.26874092915906051d2*u4020)
  u4020=6.73451707969374607d-2*u3921
  u416=-3.03053268586218573d-1*u937
  u480=-2.39917170964089704d-1*u937
  u4081=1.26272195244257739d-2*u3921
  u3849=u4020*u2980
  u4131=u4081*u2980
  value_5_=value_5_+(c(59)*(u2937*(u2886*(u651+u2980*(u63*u2980+u654)+(&
 &u2980*(u3876+u3849)+u3931)*u2886)+u2885*(u652+u2980*(u2980*(u416+u413&
 &1)+u392)+(u2980*(u4119+u4015)+u477)*u2885)+u2980*(u653+u2980*(u2980*(&
 &u480+u4015)+u655))+u901)))
  u466=-2.79236183301633649d-4*u902
  u902=3.07159801631797014d-3*u4097
  u903=u420*(u498-u3861*(u4005-4.5d1))
  u420=1.59723096848534447d-1*u457
  u3985=1.46465997518785565d8*u830
  u408=u4026*(-u4003+u3985)
  u607=-7.98615484242672237d-2*u408
  u4111=-u492
  u492=-3.99307742121336118d-2*u4026
  u656=+u492*(u4087-u484)
  u867=-2.28854374304773626d7*u4064
  u4103=3.18404342432142532d7*u830
  u939=u4026*(u4103+u3873)
  u608=2.39584645272801671d0*u939
  u421=-5.98961613182004178d-1*u937
  u657=3.99307742121336118d-2*u4026*(u4087+u3872)
  u658=7.98615484242672237d-1*u937
  u4043=u4026*(1.84037709925778383d9*u830+u3872*(1.7d1+pd2))
  u473=-7.98615484242672237d-2*u4043
  u860=u463+u473
  u3773=u2980*(u860)
  value_6_=value_6_+(c(59)*(u2886*(u902+u2980*(u2980*(u867+u2980*(u3805&
 &+u658))+u420)+(u2980*(u4107+u2980*(u771+u4062))+u551)*u2886)+u2885*(u&
 &942+u2980*(u2980*(u608+u3773)+u607)+(u2980*(u4111+u2980*(u463+u421))+&
 &u546)*u2885)+u2980*(u903+u2980*(u657*u2980+u656))+u466))
  u3860=2.9930008188621398d8*u870
  u3907=4.7d1*pd2
  u3925=u873*(u3860-u3861*(6.0d1+u3907))
  u4040=4.34389557283516869d-3*u3925
  u3996=2.12269561621428355d6*u830
  u56=u4026*(-u3975+u3996)
  u496=-4.06588625617371789d0*u56
  u985=-u659
  u659=-5.6470642446857193d-2*u3994
  u866=8.99024888750610046d5*u4064
  u3791=2.92931995037571129d7*u830
  u928=u4026*(-u3861+u3791)
  u850=-4.70588687057143275d-1*u928
  u564=1.9741069230792837d8*u830
  u951=u4026*(u564+u3872)
  u547=5.6470642446857193d-1*u951
  u4186=-5.6470642446857193d-1*u937
  u4076=3.59609955500244019d6*u4064
  u483=-9.4117737411428655d-2*u937
  u79=u4026*(3.33263211745642517d8*u830+u4008)
  u609=8.47059636702857895d-2*u79
  u4121=u4026*(u554+u3862*(2.6d1+pd2))
  u554=-2.25882569787428772d-1*u4121
  u3884=-5.6470642446857193d-2*u937
  u93=9.41177374114286549d-3*u3921
  u452=u93*u2980
  u3834=u2980*(u483+u452)
  u3980=u3884*u2980
  value_7_=value_7_+(c(59)*(u2937*(u2886*(u496+u2980*(u2980*(u554+u4011&
 &)+u547)+(u2980*(u4186+u858)+u985)*u2886)+u2885*(u659+u2980*(u3980+u40&
 &76)+(u3834+u866)*u2885)+u2980*(u850+u2980*(u888*u4191+u609))+u4040)))
  u4040=2.35619213399785474d8*u3898
  u496=7.4d1*u4027
  u850=-2.46263166022323189d-4*u3779*(-u4003*(u496+1.5d1*u668)+u4040)
  u547=1.89622637837188855d-2*u4097
  u554=-1.05646898223576648d-1*u3994
  u4125=8.85164071961356239d8*u870
  u4075=1.39d2*pd2
  u397=2.70889482624555508d-3*u873*(u4125-u3861*(1.65d2+u4075))
  u526=-9.86037716753382047d-1*u3994
  u660=-2.11293796447153296d-1*u408
  u408=-u494
  u494=8.27851290323570583d7*u830
  u523=u4026*(u3819+u494)
  u476=-1.7607816370596108d-1*u523
  u4071=4.70938036224227734d7*u4064
  u4070=6.33881389341459888d0*u939
  u4200=-u3917
  u3828=u4026*(2.56846169561928309d8*u830+u3872)
  u973=1.05646898223576648d-1*u3828
  u487=-9.86037716753382047d-1*u937
  u819=-2.11293796447153296d-1*u4043
  u453=-1.40862530964768864d-1*u937
  value_8_=value_8_+(c(59)*((u547+u2980*(u2980*(u4071+u2980*(u3875+u487&
 &))+u526))*u2886+u2885*(u949+u2980*(u2980*(u4070+u2980*(u839+u819))+u6&
 &60)+(u2980*(u408+u2980*(u839+u4200))+u554)*u2885)+u2980*(u397+u2980*(&
 &u2980*(u973+u453*u2980)+u476))+u850))
  u850=2.29857348139126367d-2*u4097
  u819=4.48221828871296415d-1*u457
  u397=-4.51933292837408604d7*u4064
  u526=2.0218068363778806d7*u4064
  u660=1.49407276290432138d0*u937
  u476=-3.73518190726080346d-1*u937
  value_9_=value_9_+(c(59)*(x*(u2885*(u819+u2980*(u2980*(u660+u4025)+u3&
 &97)+u3855)+u2980*(u566+u2980*(u2980*(u476+u4037)+u526))+u850)*z))
  u973=1.49407276290432138d-1*u457
  u487=-4.75719255618324846d6*u4064
  u453=+u487
  u3855=1.58573085206108282d6*u4064
  u547=4.98024254301440461d-1*u937
  u550=1.52230161797863951d7*u4064
  u851=-5.47826679731584507d-1*u937
  u4201=-1.49407276290432138d-1*u937
  u3959=(u2980*(u547+u3824)+u453)
  u524=u3959*u2885
  u4132=u4201*u2980
  value_1_=value_1_+(c(60)*(u2934*(u2885*(u973+u2980*(u2980*(u851+u4036&
 &)+u3855)+u524)+u2980*(u566+u2980*(u4132+u550))+u850)))
  u851=1.62533689574733305d-2*u4097
  u835=-6.72768623177468191d5*u4064
  u481=+u835
  u528=2.11293796447153296d-1*u937
  u746=2.69107449270987277d7*u4064
  u82=-1.23254714594172756d0*u937
  u4205=+u82
  u4104=-3.16940694670729944d-1*u937
  u412=(u528+u3882)
  u3972=u2980*u412
  value_2_=value_2_+(c(60)*(y*(u2885*(u553+u2980*(u2980*(u4205+u839)+u4&
 &155)+(u3972+u481)*u2885)+u2980*(u570+u2980*(u4104*u2980+u746))+u851)*&
 &z))
  u4205=u873*(-u3988*(u407+1.5d1)+u498)
  u507=-5.21267468740220243d-2*u4205
  u4169=u759*(1.38187484615549859d9*u830+u806)
  u4157=-2.82353212234285965d-1*u4026*(6.09213641853499378d8*u830+u4023&
 &)
  u661=1.31764832376000117d0*u937
  u522=-3.38823854681143158d-1*u4026*(u4156+u3872)
  u4156=u759*(2.3491872384643476d10*u830+u3872*(2.17d2+u844))
  u4158=u2980*(u4156+u3961)
  u3796=u538*u2980
  value_3_=value_3_+(c(60)*(u2934*(u2885*(u958+u2980*(u4158+u4157)+(u29&
 &80*(u661+u3961)+u4067)*u2885)+u2980*(u4169+u2980*(u3796+u522))+u507))&
 &)
  u3927=u873*(-u4003*(u4099+2.7d1)+u498)
  u439=-6.14319603263594028d-3*u3927
  u964=u490*(1.48588693134999848d7*u830-u3861)
  u922=-1.77997846681490598d6*u4064
  u4162=+u922
  u965=4.3923851633346973d-1*u4026*(2.02620945184090702d7*u830-u3861)
  u3856=+u492*(1.56018127791749841d9*u830+u3902)
  u574=5.59030838969870566d-1*u937
  u3858=2.22883039702499772d8*u830
  u945=u4026*(u3858+u3872)
  u4160=-1.59723096848534447d-1*u945
  u915=1.89450583747124806d10*u830
  u3957=u77*(u915+u3872*(1.75d2+u4016))
  u40=-9.31718064949784276d-2*u3921
  u575=2.79515419484935283d-1*u937
  u3935=u40*u2980
  u3786=u2980*(u3957+u3935)
  u742=(y*(u2885*(u964+u2980*(u3786+u3856)+(u2980*(u574+u3935)+u4162)*u&
 &2885)+u2980*(u965+u2980*(u575*u2980+u4160))+u439)*z)
  value_4_=value_4_+(c(60)*u742)
  u4=-5.65134300393880789d-3*u3779*(-u3981*(6.4d1*u4027+3.0d0*u3952)+u4&
 &26)
  u3817=3.18404342432142532d7*u870
  u766=-u4003*(1.0d1*pd2-3.0d0)
  u3844=u873*(u766+u3817)
  u508=-3.88529831520793042d-3*u3844
  u3991=u4026*(u449-u3988)
  u774=6.73451707969374607d-2*u3991
  u672=-3.21644922353688006d5*u4064
  u781=+u672
  u904=5.82794747281189563d-3*u873*(1.33729823821499863d8*u870-u4003*(2&
 &.9d1+4.2d1*pd2))
  u662=-1.32585805006470626d-1*u4026*(-u3861+1.36459003899489657d7*u830&
 &)
  u3842=1.08555161294369702d6*u4064
  u905=u776*(6.81385292804785018d8*u870-u4003*(1.23d2+2.14d2*pd2))
  u776=3.85973906824425607d6*u4064
  u3846=2.86563908188928279d8*u830
  u940=u4026*(u3846+u4008)
  u743=-5.05088780977030955d-2*u940
  u562=1.01017756195406191d-1*u937
  u4102=5.22183121588713752d7*u830
  u663=-5.68224878599159824d-1*u4026*(-u3861+u4102)
  u81=7.57633171465546432d-2*u4026
  u664=u81*(8.59691724566784836d8*u830+u4023)
  u820=-3.40934927159495895d-1*u937
  u4167=-2.10453658740429565d-3*u4026
  u665=+u4167*(5.76260196006401371d4*exppd2+3.40692646402392509d9*u830)
  u890=1.62386214640392691d9*u830
  u489=u4026*(u890-9.07496371663624206d2*u4018)
  u4161=3.36725853984687303d-2*u489
  u778=-1.68362926992343652d-2*u3921
  u666=1.32585805006470626d-1*u4026*(1.18719333392556001d9*u830+u3902)
  u4192=-6.31360976221288694d-3*u4026
  u64=+u4192*(4.38442779529060266d10*u830+u3872*(4.05d2+u4031))
  u3789=1.26272195244257739d-2*u4026
  u667=u3789*(9.48844940447784745d8*u830+u4023)
  u3848=4.87158643921178074d9*u830
  u3995=u4026*(u3848-9.07496371663624206d2*u714)
  u65=-5.05088780977030955d-2*u3995
  u520=1.13644975719831965d-1*u3921
  u747=-5.68224878599159824d-2*u937
  u4168=u778*u2980
  u3936=u520*u2980
  u4133=u562*u2980
  value_5_=value_5_+(c(60)*(u2886*(u508+u2980*(u2980*(u743+u4133)+u776)&
 &+u2886*((u781+u2980*(u4168+u562))*u2886+u2980*(u743+u2980*(u4168+u416&
 &1))+u774))+u2885*(u904+u2980*(u2980*(u666+u2980*(u3881+u65))+u663)+u2&
 &885*((u3842+u2980*(u3881+u820))*u2885+u2980*(u664+u2980*(u3936+u64))+&
 &u662))+u2980*(u905+u2980*(u2980*(u667+u747*u2980)+u665))+u4))
  value_6_=value_6_+(c(60)*u879)
  u879=2.63266398353646587d-4*u4108
  u405=1.44796519094505623d-3*u792
  u765=u4026*(-u3861+u4087)
  u422=-2.82353212234285965d-2*u765
  u45=1.25863484425085406d6*u4064
  u3954=4.90342687345499499d8*u870
  u754=7.7d1*pd2
  u3857=+u885*(-u3861*(u754+5.1d1)+u3954)
  u758=u759*(u3858+u484)
  u941=u4026*(u3858+u4008)
  u714=1.69411927340571579d-1*u941
  u4089=-3.95294497128000351d-1*u937
  u886=9.41177374114286549d-3*u4026
  u931=8.62121553080442996d3*exppd2
  u3893=u886*(u832+u931)
  u3967=u4026*(2.72412604080833055d8*u830+u4008)
  u4172=-2.54117891010857368d-1*u3967
  u387=u4026*(u3851-9.07496371663624206d2*u4033)
  u954=-2.82353212234285965d-2*u387
  u808=1.04012085194499894d8*u830
  u974=u4026*(u808+u3862)
  u440=-1.12941284893714386d-1*u974
  u4033=3.7647094964571462d-2*u425
  u966=6.58824161880000585d-2*u937
  u41=-6.58824161880000585d-2*u3921
  u4000=u41*u2980
  value_7_=value_7_+(c(60)*(u2885*(u405+u2980*(u2980*(u4172+u2980*(u400&
 &0+u4033))+u758)+u2885*((u45+u2980*(u772+u4089))*u2885+u2980*(u714+u95&
 &4*u2980)+u422))+u2980*(u3857+u2980*(u2980*(u440+u966*u2980)+u3893))+u&
 &879))
  u4172=3.87371960153114376d-1*u457
  u405=-u4155
  u422=1.05646898223576648d0*u937
  u3857=2.69107449270987277d6*u4064
  u758=-3.5215632741192216d-2*u937
  u714=(u2980*(u422+u3827)+u874)
  u3893=(u893+u3875)
  value_8_=value_8_+(c(60)*(x*(u2885*(u4172+u2980*(u2980*u3893+u405)+u7&
 &14*u2885)+u2980*(u553+u2980*(u758*u2980+u3857))+u506)*z))
  u4172=u3779*(u426-u3988*(1.5d1*u4061+u782))
  u954=6.96537418603413232d-4*u4172
  u440=-1.91547790115938639d-2*u3886
  u966=3.73518190726080346d-2*u457
  u879=2.37859627809162423d5*u4064
  u782=7.71937594167232714d-1*u457
  u3877=-7.47036381452160691d-2*u937
  u785=-2.25966646418704302d7*u4064
  u668=9.33795476815200865d-1*u937
  u3952=7.9286542603054141d5*u4064
  u669=9.96048508602880921d-2*u937
  u784=-1.24506063575360115d-2*u937
  u4134=u784*u2980
  value_9_=value_9_+(c(60)*(u2885*(u440+u2980*(u2980*(u785+u2980*(u4037&
 &+u669))+u782)+u2885*((u879+u2980*(u4037+u3877))*u2885+u2980*(u4116+u2&
 &980*(u4025+u668))+u966))+u2980*(u440+u2980*(u2980*(u3952+u4134)+u966)&
 &)+u954))
  u555=-4.98024254301440461d-2*u3994
  u441=-2.49012127150720231d-1*u3994
  u4113=-u4181
  u4181=5.70863106741989815d6*u4064
  u4105=-7.47036381452160692d-1*u937
  u779=u555+u2980*(u2980*(u4105+u4036)+u4113)
  value_1_=value_1_+(c(61)*(x*(u2885*(u779+u524)+u2980*(u441+u2980*(u41&
 &30+u4181))+u899)*z))
  u441=8.0d0*u4027
  u524=3.0d0*u563
  u506=u3779*(u426-u3975*(u524+u441))
  u563=5.91031598453575652d-3*u506
  u3913=u873*(u3987*(u4031+3.0d0)+u829)
  u509=-6.50134758298933218d-2*u3913
  u899=u4026*(u449-u4047)
  u3901=-3.5215632741192216d-2*u899
  u4114=-u835
  u835=-u3981*(6.4d1*pd2-2.7d1)
  u934=u873*(u835+u829)
  u584=-2.60053903319573287d-1*u934
  u3916=2.11293796447153296d-1*u4026
  u4090=u3916*(u4123-u4003)
  u877=u4039*u3966*pd2**3
  u3966=-1.91748353630823987d2*u877
  u4098=+u3966
  u580=5.9860016377242796d7*u830
  u461=u4026*(u580-u3861)
  u831=1.7607816370596108d-1*u461
  u4039=u4026*(1.2491247280030207d8*u830+u3872)
  u841=-1.37340967690649642d0*u4039
  u4130=+u841
  u4030=3.5215632741192216d-2*u4026
  u471=9.07496371663624206d2*u4032
  u4032=u4030*(u890+u471)
  u4091=9.55213027296427596d7*u830
  u946=u4026*(u4091+u3862)
  u4045=-2.11293796447153296d-1*u946
  u497=u4026*(u4123+u717*(1.6d1+pd2))
  u3963=3.38070074315445273d0*u497
  u52=-7.04312654823844319d-2*u3921
  u530=1.05646898223576648d-1*u937
  u3850=u530*u2980
  u3937=u52*u2980
  value_2_=value_2_+(c(61)*(u2885*(u509+u2980*(u2980*(u4130+u2980*(u382&
 &7+u3963))+u4090)+u2885*((u4114+u2980*(u3875+u3833))*u2885+u2980*(u409&
 &8+u2980*(u3937+u4032))+u3901))+u2980*(u584+u2980*(u2980*(u4045+u3850)&
 &+u831))+u563))
  u563=-2.89593038189011246d-3*u873
  u509=+u563*(-u3861*(1.13d2*pd2-3.0d1)+7.19593813896642122d8*u870)
  u4032=1.12941284893714386d-1*u4109
  u584=-1.07882986650073206d7*u4064
  u4090=+u584
  u4098=2.44706117269714503d-1*u457
  u831=-u849
  u4130=u4026*(u494-u3861)
  u3901=2.82353212234285965d-1*u4130
  u4045=-5.6470642446857193d-1*u938
  u4073=1.12941284893714386d0*u937
  u894=-1.9778547552513421d7*u4064
  u576=1.8823547482285731d-1*u937
  u3867=-6.77647709362286316d-1*u4026*(8.70305202647856254d7*u830+u3862&
 &)
  u4061=1.12941284893714386d-1*u394
  u406=-1.12941284893714386d-1*u3921
  u670=5.08235782021714737d-1*u937
  u671=3.57647402163428889d-1*u937
  u3869=u2980*(u670+u3879)
  u3998=u406*u2980
  value_3_=value_3_+(c(61)*(u2937*(u2886*(u4032+u2980*(u2980*(u4061+u39&
 &98)+u4045)+(u2980*(u4073+u3998)+u4090)*u2886)+u2885*(u4098+u2980*(u38&
 &69+u894)+(u2980*(u576+u3879)+u831)*u2885)+u2980*(u3901+u2980*(u671*u2&
 &980+u3867))+u509)))
  u3778=u3779*(-u3975*(u441+u524)+u426)
  u524=-2.23388946641306919d-3*u3778
  u441=-1.22863920652718806d-2*u3844
  u967=2.12964129131379263d-1*u3991
  u50=-1.01713055246566056d6*u4064
  u794=+u50
  u446=u873*(u498-u3975*(u3993-9.0d0))
  u834=1.22863920652718806d-2*u446
  u767=u4026*(u4047+u4103)
  u577=-3.99307742121336118d-2*u767
  u4084=7.62847914349245418d5*u4064
  u3993=u873*(u829+u835)
  u835=9.82911365221750445d-2*u3993
  u4126=1.22055666295879267d7*u4064
  u3794=-1.59723096848534447d-1*u940
  u531=3.19446193697068895d-1*u937
  u715=u4026*(-u4003+u491)
  u578=-7.98615484242672237d-2*u715
  u579=4.79169290545603342d-1*u946
  u3900=-2.39584645272801671d-1*u937
  u3965=u4026*(-u3861+u580)
  u580=-6.65512903535560198d-2*u3965
  u448=1.06482064565689632d-1*u489
  u72=u4026*(1.31910370436173335d8*u830+u3872)
  u581=8.38546258454805848d-1*u72
  u46=-3.99307742121336118d-2*u3995
  u582=7.98615484242672237d-2*u946
  u3908=u4026*(u801-2.26874092915906051d2*u809)
  u3911=-3.19446193697068895d-1*u3908
  u623=7.98615484242672237d-2*u3921
  u789=-3.99307742121336118d-2*u937
  u3938=u623*u2980
  u4135=u531*u2980
  u712=(u2886*(u441+u2980*(u2980*(u3794+u4135)+u4126)+u2886*((u794+u298&
 &0*(u4188+u531))*u2886+u2980*(u3794+u2980*(u4188+u448))+u967))+u2885*(&
 &u834+u2980*(u2980*(u581+u2980*(u463+u3911))+u578)+u2885*((u4084+u2980&
 &*(u463+u3900))*u2885+u2980*(u579+u2980*(u3938+u46))+u577))+u2980*(u83&
 &5+u2980*(u2980*(u582+u789*u2980)+u580))+u524)
  value_4_=value_4_+(c(61)*u712)
  u906=2.33117898912475825d-2*u4166
  u4183=1.68362926992343652d-2*u457
  u4184=-u672
  u672=-2.52544390488515477d-2*u765
  u427=-1.76781073341960834d-1*u4026
  u673=+u427*(2.63820740872346669d7*u830+u3861)
  u881=-8.04112305884220015d6*u4064
  u467=-1.01017756195406191d-1*u937
  u3910=3.78816585732773216d-2*u4026
  u674=u3910*(1.05073433002607036d9*u830+u3902)
  u675=u880*(1.00828041770178468d9*u830+u3902)
  u676=5.89270244473202781d-1*u937
  u3889=-1.01017756195406191d-1*u489
  u468=-2.14662731915238156d-1*u937
  value_5_=value_5_+(c(61)*(u2936*(u2886*(u4183+u2980*(u2980*(u676+u379&
 &5)+u881)+(u2980*(u467+u3933)+u4184)*u2886)+u2885*(u672+u2980*(u2980*(&
 &u3889+u4002)+u674)+(u2980*(u820+u3881)+u3842)*u2885)+u2980*(u673+u298&
 &0*(u2980*(u468+u4015)+u675))+u906)))
  value_6_=value_6_+(c(61)*u495)
  u495=1.59202171216071266d8*u870
  u793=1.44796519094505623d-3*u873
  u773=u793*(u495-u3861*(1.2d1+u4007))
  u3933=u4026*(-u4003+u4103)
  u711=-3.7647094964571462d-2*u3933
  u3964=5.5493328252459127d7*u830
  u787=-1.31764832376000117d-1*u4026*(-u3861+u3964)
  u4046=1.61824479975109808d7*u4064
  u975=u4026*(2.44109995864642608d8*u830+u4008)
  u442=1.69411927340571579d-1*u975
  u4195=u4026*(3.07790864351071114d8*u830+u4008)
  u677=8.47059636702857895d-2*u4195
  u458=-8.47059636702857895d-1*u937
  u388=u4026*(u801-1.13437046457953026d2*u4038)
  u611=-3.01176759716571696d-1*u388
  u853=-7.52941899291429239d-2*u937
  value_7_=value_7_+(c(61)*(u2936*((u659+u2980*(u2980*(u458+u858)+u4046&
 &))*u2886+u2885*(u711+u2980*(u2980*(u611+u858)+u442)+(u2980*(u4089+u77&
 &2)+u45)*u2885)+u2980*(u787+u2980*(u2980*(u853+u3887)+u677))+u773)))
  u773=1.46465997518785565d8*u870
  u711=-u3861*(u3814-3.0d0)
  u787=u873*(u773+u711)
  u442=1.35444741312277754d-2*u787
  u677=-2.46509429188345512d-1*u3994
  u611=-1.05646898223576648d-1*u4165
  u853=-u953
  u953=-7.0431265482384432d-1*u3933
  u479=2.35469018112113867d7*u4064
  u612=5.28234491117883239d-1*u938
  u4202=-u422
  u780=u4026*(1.80429127378214101d8*u830+u3872)
  u976=3.16940694670729944d-1*u780
  u3870=-7.39528287565036535d-1*u937
  u53=-1.05646898223576648d-1*u394
  u723=(u2980*(u4202+u839)+u853)
  u3939=u450*u2980
  value_8_=value_8_+(c(61)*(u2934*((u677+u2980*(u2980*(u3870+u3875)+u47&
 &9))*u2886+u2885*(u611+u2980*(u2980*(u53+u839)+u612)+u723*u2885)+u2980&
 &*(u953+u2980*(u3939+u976))+u442)))
  u442=-7.66191160463754555d-3*u3886
  u953=2.49012127150720231d-1*u457
  u976=-1.66501739466413696d7*u4064
  u840=+u976
  u3888=-u4083
  u678=9.96048508602880922d-1*u937
  u4058=u2980*(u3877+u4037)
  value_9_=value_9_+(c(61)*(y*(u2885*(u944+u2980*(u2980*(u678+u4025)+u8&
 &40)+(u4058+u879)*u2885)+u2980*(u953+u2980*(u429+u3888))+u442)*z))
  u429=2.73827734491642577d8*u870
  u404=u873*(u429-u3861*(3.0d1+u4010))
  u679=3.83095580231877278d-3*u404
  u583=-3.48616978011008323d-1*u3994
  u613=-4.98024254301440461d-2*u4165
  u887=-u487
  u487=1.20993650124214162d8*u830
  u889=u4026*(u3819+u487)
  u907=-9.96048508602880921d-2*u889
  u4199=-u614
  u614=2.49012127150720231d-1*u938
  u4196=-4.98024254301440461d-1*u937
  u435=2.10146866005214071d8*u830
  u968=u4026*(u435+u3872)
  u852=1.49407276290432138d-1*u968
  u4198=-u3906
  u54=-4.98024254301440461d-2*u394
  u908=-1.99209701720576184d-1*u937
  value_1_=value_1_+(c(62)*(u2934*((u583+u2980*(u2980*(u4198+u4036)+u41&
 &99))*u2886+u2885*(u613+u2980*(u2980*(u54+u4036)+u614)+(u2980*(u4196+u&
 &4036)+u887)*u2885)+u2980*(u907+u2980*(u908*u2980+u852))+u679)))
  u852=u873*(u3817-u3861*(1.8d1+u4021))
  u907=8.12668447873666522d-3*u852
  u770=u4026*(u3819+u4087)
  u908=-2.11293796447153296d-1*u770
  u679=1.05646898223576648d-1*u941
  u760=4.13925645161785292d8*u830
  u488=u4026*(u760+u3872)
  u680=1.05646898223576648d-1*u488
  u659=u4026*(2.70643691067321152d9*u830-9.07496371663624206d2*u872)
  u66=-3.5215632741192216d-2*u659
  u558=-4.22587592894306592d-1*u937
  u445=u2980*(u66+u3875)
  u4203=u2980*(u3833+u3875)
  u3940=u558*u2980
  value_2_=value_2_+(c(62)*(u2936*((u554+u2980*(u2980*(u4200+u839)+u408&
 &))*u2886+u2885*(u481+u2980*(u445+u679)+(u4203+u4114)*u2885)+u2980*(u9&
 &08+u2980*(u3940+u680))+u907)))
  u908=1.44796519094505623d-3*u3925
  u680=-5.45882876986286199d-1*u3994
  u3925=-u584
  u584=-1.8823547482285731d-2*u4165
  u3979=1.91042605459285519d7*u830
  u478=u4026*(-u4003+u3979)
  u681=-2.25882569787428772d-1*u478
  u3926=3.05668462175207416d7*u4064
  u3896=-u4073
  u585=9.4117737411428655d-2*u938
  u460=-1.8823547482285731d-1*u937
  u386=1.59202171216071266d8*u830
  u943=u4026*(u386+u3872)
  u615=5.6470642446857193d-2*u943
  u799=-2.82353212234285965d-1*u937
  u47=-1.8823547482285731d-2*u394
  u91=1.8823547482285731d-2*u3921
  u3781=-3.7647094964571462d-2*u937
  u3777=u91*u2980
  u804=u2980*(u47+u3777)
  value_3_=value_3_+(c(62)*(u2934*(u2886*(u680+u2980*(u2980*(u799+u3879&
 &)+u3926)+(u2980*(u3896+u4011)+u3925)*u2886)+u2885*(u584+u2980*(u804+u&
 &585)+(u2980*(u460+u3777)+u849)*u2885)+u2980*(u681+u2980*(u3781*u2980+&
 &u615))+u908)))
  u816=3.07159801631797014d-3*u792
  u586=-6.65512903535560198d-2*u3994
  u395=-u50
  u809=-u4084
  u50=u4026*(u4087-u4047)
  u587=-7.98615484242672237d-2*u50
  u3976=8.89989233407452989d6*u4064
  u3847=-3.19446193697068895d-1*u937
  u527=1.19792322636400836d-1*u941
  u556=u490*(u588+u3872)
  u588=6.65512903535560198d-2*u937
  u3821=-3.99307742121336118d-2*u659
  u862=u3821+u463
  u798=u2980*(u862)
  u4185=u2980*(u798+u527)
  u755=u2980*(u3900+u463)
  u4136=u556*u2980
  u599=(u2936*(u2886*(u586+u2980*(u2980*(u588+u3805)+u3976)+(u2980*(u38&
 &47+u771)+u395)*u2886)+u2885*(u809+u4185+(u755+u4084)*u2885)+u2980*(u5&
 &87+u4136)+u816))
  value_4_=value_4_+(c(62)*u599)
  u5=-5.73964523837535176d-4*u3779*(-u4003*(2.0d0*u4027+3.0d0*u4069)+u4&
 &26)
  u426=4.85662289400991303d-4*u873
  u909=u426*(6.55912945410213616d8*u870-u4003*(2.06d2*pd2-3.27d2))
  u682=-2.02035512390812382d-1*u895
  u4127=1.28657968941475202d6*u4064
  u3792=2.61091560794356876d8*u870
  u510=-4.85662289400991303d-4*u873*(-u4003*(8.2d1*pd2-9.0d0)+u3792)
  u3983=u4026*(u449-u4003)
  u891=3.78816585732773216d-2*u3983
  u3809=-1.20616845882633002d5*u4064
  u4085=+u3809
  u910=1.45698686820297391d-3*u873*(u829-u4003*(u3774-1.41d2))
  u3801=3.50244776675356785d8*u830
  u683=+u414*(u3801-u616)
  u616=1.21221307434487429d0*u939
  u842=-4.04071024781624764d-1*u937
  u474=8.08142049563249528d-1*u4026*(u449-u3981)
  u4018=4.77606513648213798d8*u830
  u828=-1.89408292866386608d-2*u4026
  u4210=+u828*(u4018+u4008)
  u532=3.78816585732773216d-2*u937
  u3981=2.18871805090046747d-1*u4026*(2.44926417255494255d6*u830-u3975)
  u933=7.32329987593927823d8*u830
  u684=u3910*(u933+u3902)
  u4069=-2.02035512390812382d-1*u4026*(u801-4.53748185831812103d2*u4017&
 &)
  u926=u4026*(5.20060425972499469d8*u830+u4023)
  u3831=-5.68224878599159824d-2*u926
  u3905=2.52544390488515477d-2*u4026
  u990=u3905*(u890+u3862*(3.0d1+pd2))
  u4100=-6.31360976221288694d-3*u3921
  u70=2.35619213399785474d8*u830
  u3899=+u4192*(u70+u4008)
  u805=-1.76781073341960834d-1*u937
  u857=u3789*(4.00552662779635305d9*u830+u3872*(3.7d1+u3774))
  u42=-1.89408292866386608d-2*u3921
  u589=6.31360976221288694d-3*u937
  u824=-1.26272195244257739d-2*u3921
  u3776=u4100*u2980
  u3818=u824*u2980
  u3977=u42*u2980
  u4137=u4069*u2980
  u4138=u589*u2980
  value_5_=value_5_+(c(62)*(u2886*(u909+u2980*(u2980*(u684+u2980*(u4015&
 &+u805))+u683)+u2886*((u4127+u2980*(u3849+u842))*u2886+u2980*(u616+u41&
 &37)+u682))+u2885*(u510+u2980*(u2980*(u3831+u2980*(u3818+u857))+u474)+&
 &u2885*((u4085+u2980*(u3776+u532))*u2885+u2980*(u4210+u2980*(u3977+u99&
 &0))+u891))+u2980*(u910+u2980*(u2980*(u3899+u4138)+u3981))+u5))
  u911=1.53579900815898507d-2*u792
  u601=8.9153215880999909d6*u830
  u685=-6.65512903535560198d-1*u4026*(-u4003+u601)
  u686=-1.3310258070711204d-1*u50
  u687=7.98615484242672237d-1*u971
  u812=-9.31718064949784276d-1*u937
  u447=-2.66205161414224079d-2*u4026
  u67=+u447*(u3851-9.07496371663624206d2*u3878)
  u859=9.31718064949784276d-2*u3921
  u3851=u859*u2980
  value_6_=value_6_+(c(62)*(u2937*(u2886*(u685+u2980*(u67*u2980+u687)+(&
 &u2980*(u812+u3851)+u3976)*u2886)+u2980*(u686+u4136)+u911)))
  u4066=2.61091560794356876d8*u3898
  u4065=8.2d1*u4027
  u6=-6.58165995884116468d-5*u3779*(-u4003*(u4065+1.5d1*u4077)+u4066)
  u3806=7.23982595472528115d-4*u873
  u423=u3806*(u469-u4003*(3.81d2+u810))
  u912=-5.6470642446857193d-2*u767
  u4028=1.07882986650073206d6*u4064
  u511=-2.17194778641758435d-3*u3844
  u4077=3.7647094964571462d-2*u3991
  u625=-1.79804977750122009d5*u4064
  u3823=+u625
  u4060=6.43176771712927914d8*u870
  u935=2.02d2*pd2
  u4171=u3806*(u4060-u4003*(3.45d2+u935))
  u3986=4.71238426799570947d7*u830
  u4164=-3.7647094964571462d-1*u4026
  u3978=+u4164*(u4047+u3986)
  u58=6.77647709362286316d-1*u946
  u762=-3.38823854681143158d-1*u937
  u4212=2.15765973300146411d6*u4064
  u4173=-2.82353212234285965d-2*u940
  u536=5.6470642446857193d-2*u937
  u913=-6.77647709362286316d-1*u56
  u56=2.01656083540356937d8*u830
  u4136=u4026*(u56+u3872)
  u3878=3.38823854681143158d-1*u4136
  u3973=-5.6470642446857193d-2*u3995
  u872=1.8823547482285731d-2*u489
  u4048=2.99674962916870015d5*u4064
  u688=-4.14118044610286082d-1*u937
  u617=9.41177374114286549d-3*u937
  u3852=u536*u2980
  value_7_=value_7_+(c(62)*(u2886*(u423+u2980*(u2980*(u3878+u2980*(u388&
 &7+u688))+u3978)+u2886*((u4028+u2980*(u858+u762))*u2886+u2980*(u58+u29&
 &80*(u4011+u3973))+u912))+u2885*(u511+u2980*(u2980*(u4173+u3852)+u4212&
 &)+u2885*((u3823+u2980*(u3887+u536))*u2885+u2980*(u4173+u2980*(u3887+u&
 &872))+u4077))+u2980*(u4171+u2980*(u2980*(u4048+u617*u2980)+u913))+u6)&
 &)
  u6=8.12668447873666522d-3*u404
  u423=1.42057322008186668d7*u830
  u912=-9.15606451270997615d-1*u4026*(-u4003+u423)
  u4171=-2.11293796447153296d-1*u889
  u3878=u4026*(u3979+u3982)
  u58=8.45175185788613183d0*u3878
  u913=3.16940694670729944d-1*u968
  u3973=u4026*(u417+u3872*(1.9d1+pd2))
  u417=-2.11293796447153296d-1*u3973
  u889=1.40862530964768864d-1*u3921
  u3978=u889*u2980
  value_8_=value_8_+(c(62)*(u2937*(u2886*(u912+u2980*(u2980*(u417+u3978&
 &)+u58)+u723*u2886)+u2980*(u4171+u2980*(u3940+u913))+u6)))
  u6=-8.7067177325426654d-5*u3779*(-u4003*(1.22d2*u4027+3.0d0*u3832)+3.&
 &88453297767213889d8*u3898)
  u912=6.70417265405785236d-3*u4097
  u4171=8.61965055521723874d-3*u873
  u913=u4171*(u498-u4003*(1.0d0+u4099))
  u417=1.27361736972857013d6*u830
  u688=-1.24506063575360115d-1*u4026*(-u4003+u417)
  u404=-u879
  u723=9.57738950579693194d-4*u873
  u914=u723*(1.4965004094310699d9*u870-u4003*(2.49d2+4.7d2*pd2))
  u689=-5.97629105161728553d-1*u4026*(-u3975+u4103)
  u418=3.73518190726080346d-2*u4026
  u690=u418*(u3846-u4008)
  u544=7.47036381452160691d-2*u937
  u691=-4.98024254301440461d0*u4026*(u3987+1.78306431761999818d6*u830)
  u4128=-u976
  u976=u4026*(5.04898314428111729d8*u830+u4023)
  u692=2.61462733508256242d-1*u976
  u3894=-u3862*(pd2-3.0d1)
  u4038=-4.98024254301440461d-2*u4026*(u3894+u890)
  u4059=-1.24506063575360115d-2*u3921
  u693=u418*(4.47888775021213828d8*u830+u4008)
  u418=-3.48616978011008323d-1*u937
  u3832=u4026*(3.57249672208863921d9*u830-9.07496371663624206d2*u3782)
  u68=-7.47036381452160691d-2*u3832
  u4017=6.22530317876800576d-2*u3921
  u3890=-8.71542445027520806d-2*u937
  u97=7.47036381452160691d-2*u3921
  u455=u4059*u2980
  u444=u2980*(u455+u544)
  u456=u97*u2980
  u3770=u4017*u2980
  value_9_=value_9_+(c(62)*((u912+u2980*(u2980*(u4128+u2980*(u4037+u418&
 &))+u583))*u2886+u2885*(u913+u2980*(u2980*(u692+u2980*(u456+u68))+u689&
 &)+u2885*((u404+u444)*u2885+u2980*(u690+u2980*(u3770+u4038))+u688))+u2&
 &980*(u914+u2980*(u2980*(u693+u3890*u2980)+u691))+u6))
  u988=-1.14928674069563183d-2*u3886
  u3782=9.96048508602880921d-2*u457
  u3793=-9.51438511236649692d5*u4064
  u407=+u3793
  u3859=1.14928674069563183d-2*u4097
  u529=2.98814552580864276d-1*u937
  u4012=-9.96048508602880921d-2*u3994
  u4044=-u3793
  u3793=-2.98814552580864276d-1*u937
  u4204=u2980*(u3824+u529)
  u3853=u3793*u2980
  u4139=u4044*u2980
  value_1_=value_1_+(c(63)*(u2885*(u988+u4191*(u3853+u4113)+u2885*((u40&
 &7+u4204)*u2885+u2980*(u4116+u89*u4191)+u3782))+u2980*(u3859+u2980*(u4&
 &139+u4012))))
  u988=-3.16940694670729944d-1*u3994
  u3859=-1.7607816370596108d-1*u3994
  u4012=-5.28234491117883239d-1*u937
  u3837=2.01830586953240457d6*u4064
  u4189=-6.33881389341459887d-1*u937
  u4140=u4189*u2980
  value_2_=value_2_+(c(63)*(x*(u2885*(u988+u2980*((u4012+u839)*u2885+u4&
 &140+u408)+(u3882+u530)*u800)+u2980*(u3859+u3837*u2980)+u851)*z))
  u988=7.89799195060939762d-4*u4108
  u3859=-4.34389557283516869d-3*u873
  u512=+u3859*(-u4003*(7.0d1*pd2+5.7d1)+2.22883039702499772d8*u870)
  u3984=3.7647094964571462d-1*u4026
  u977=u3984*(u601-u4003)
  u790=-2.51726968850170813d6*u4064
  u878=u4026*(u4087-u3861)
  u443=6.77647709362286316d-1*u878
  u525=6.68649119107499317d8*u830
  u57=-1.69411927340571579d-1*u4026
  u618=+u57*(u525+u4023)
  u590=7.90588994256000702d-1*u937
  u986=1.12941284893714386d-1*u413
  u4141=u590*u2980
  u4142=u790*u2980
  value_3_=value_3_+(c(63)*(u2885*(u512+u2980*(u2980*(u618+u4141)+u443)&
 &+u2885*((u790+u2980*(u3961+u590))*u2885+u2980*(u618+u2980*(u3961+u986&
 &))+u977))+u2980*(u512+u2980*(u4142+u977))+u988))
  u988=(x*(u2885*(u965+u2980*(u574*u2980+u3856)+u2885*((u575+u3935)*u28&
 &85+u3786+u4160))+u2980*(u964+u4162*u2980)+u439)*z)
  value_4_=value_4_+(c(63)*u988)
  u443=3.88529831520793042d-3*u873*(5.03078861042785201d8*u870-u3861*(2&
 &.7d1+7.9d1*pd2))
  u3786=u4026*(u4103+u3861)
  u977=5.05088780977030955d-2*u3786
  u986=-2.02035512390812382d-1*u939
  u618=5.05088780977030955d-2*u937
  u4187=-9.09159805758655719d-1*u478
  u694=u880*(1.77669623077135533d9*u830+u3902)
  u745=-1.70467463579747947d-1*u937
  u619=-4.54579902879327859d-1*u4026*(-u3861+3.75009558864523426d7*u830&
 &)
  u4206=-5.05088780977030955d-2*u941
  u4144=1.68362926992343652d-2*u659
  u4055=-1.45199419466179873d4*u3989
  u533=u69*(2.34982404714921189d9*u830+u4055)
  u69=+u4192*(3.99470088015366021d10*u830+u3872*(3.69d2+u4031))
  u4057=u880*(1.75971466584164106d9*u830+u3902)
  u871=(u618+u4168)*u2886
  u4174=(u745+u3881)*u2885
  u4143=u745*u2980
  value_5_=value_5_+(c(63)*(u2934*(u2886*(u977+u2980*(u4133+u4206)+u288&
 &6*(u871+u2980*(u4144+u4168)+u986))+u2885*(u4187+u2980*(u2980*(u69+u38&
 &81)+u533)+u2885*(u4174+u2980*(u69+u3936)+u694))+u2980*(u619+u2980*(u4&
 &143+u4057))+u443)))
  value_6_=value_6_+(c(63)*u742)
  u443=-1.8823547482285731d-2*u892
  u977=u4026*(u4124+u4008)
  u986=8.47059636702857895d-2*u977
  u4187=-1.97647248564000175d-1*u937
  u4206=u4026*(u4124+u3819)
  u619=1.8823547482285731d-2*u4206
  u694=u4026*(u3971-9.07496371663624206d2*u3825)
  u3825=-2.82353212234285965d-2*u694
  u533=-8.47059636702857895d-2*u977
  u742=2.82353212234285965d-2*u694
  u694=1.97647248564000175d-1*u937
  u4057=(u4187+u772)*u2885
  u4144=u694*u2980
  value_7_=value_7_+(c(63)*(u2934*(u2885*(u443+u4191*(u4000+u742)+u2885&
 &*(u4057+u3825*u2980+u986))+u2980*(u619+u2980*(u4144+u533)))))
  u443=-1.62533689574733305d-2*u3886
  u986=1.7607816370596108d-1*u457
  u3825=-u3837
  u742=3.16940694670729944d-1*u457
  u533=6.33881389341459887d-1*u937
  u750=u2980*(u533+u3827)
  value_8_=value_8_+(c(63)*(y*(u2885*(u986+u2980*(u2980*(u559+u3875)+u4&
 &54)+(u750+u3825)*u2885)+u2980*(u742+u893*u4191)+u443)*z))
  u986=-2.14073665028246181d6*u4064
  u742=-3.73518190726080346d-2*u937
  u619=-3.0921751615191115d7*u4064
  u695=4.85573647943904449d-1*u937
  u982=(u742+u4037)*u2885
  u4145=u742*u2980
  value_9_=value_9_+(c(63)*(u2934*(u2885*(u950+u2980*(u2980*(u695+u4037&
 &)+u619)+u2885*(u982+u2980*(u695+u4025)+u986))+u2980*(u950+u2980*(u414&
 &5+u986))+u451)))
  u619=-1.49407276290432138d-1*u3994
  u986=9.51438511236649692d6*u4064
  u783=-2.49012127150720231d-1*u937
  u428=(u529+u3824)
  u3830=u2980*u428
  value_1_=value_1_+(c(64)*(y*(u2885*(u944+u2980*(u2980*(u783+u4036)+u4&
 &53)+(u3830+u407)*u2885)+u2980*(u619+u2980*(u4132+u986)))*z))
  u783=-u3975*(u4016-5.0d0)
  u4132=-1.95040427489679965d-1*u873*(u783+u829)
  u513=7.0431265482384432d-1*u3991
  u695=7.0048955335071357d7*u830
  u978=u4026*(u695+u3872)
  u537=1.05646898223576648d-1*u978
  u4133=1.40862530964768864d-1*u4026*(u487-u3988)
  u384=-3.5215632741192216d-1*u938
  u768=u4026*(u4123+u3872*(1.0d0+u3774))
  u518=-3.5215632741192216d-2*u768
  u930=u4026*(1.54956779983642699d8*u830+u3872)
  u4190=-3.16940694670729944d-1*u930
  u3958=u4026*(8.98537054343506225d9*u830+u3872*(8.3d1+u3960))
  u3955=3.5215632741192216d-2*u3958
  u549=3.16940694670729944d-1*u937
  u3813=u3893*u2885
  u3941=u549*u2980
  value_2_=value_2_+(c(64)*(u2934*(u2885*(u513+u2980*(u2980*(u3955+u382&
 &7)+u384)+u2885*(u3813+u2980*(u518+u3937)+u537))+u2980*(u4133+u2980*(u&
 &3941+u4190))+u4132)))
  u4132=u873*(u598+u498)
  u513=-8.68779114567033738d-3*u4132
  u537=5.6470642446857193d-2*u878
  u4133=5.6470642446857193d-2*u4026
  u384=u4133*(u3858+u3861)
  u518=-1.69411927340571579d-1*u926
  u4190=-6.77647709362286316d-1*u971
  u3955=u759*(u915-9.07496371663624206d2*u409)
  u926=u2980*(u3955+u3961)
  value_3_=value_3_+(c(64)*(y*(u2885*(u537+u2980*(u926+u518)+(u2980*(u5&
 &90+u3961)+u790)*u2885)+u2980*(u384+u2980*(u3796+u4190))+u513)*z))
  u915=1.96582273044350089d-1*u3993
  u598=1.59723096848534447d-1*u3786
  u4194=-6.38892387394137789d-1*u939
  u557=1.59723096848534447d-1*u937
  u744=u4026*(u3987+u449)
  u696=-1.91667716218241337d0*u744
  u697=1.19792322636400836d-1*u951
  u4182=-1.19792322636400836d-1*u937
  u3948=2.3349651778357119d7*u830
  u698=-4.79169290545603342d-1*u4026*(-u3988+u3948)
  u769=-1.59723096848534447d-1*u941
  u748=5.32410322828448158d-2*u659
  u699=3.99307742121336118d-1*u938
  u3799=4.4385565335040669d9*u830
  u971=u4026*(u3799+u3872*(4.1d1+u3774))
  u48=-3.99307742121336118d-2*u971
  u700=u490*(1.88919909843071236d8*u830+u3872)
  u516=(u557+u4188)*u2886
  u718=(u4182+u463)*u2885
  u3942=u4182*u2980
  u4163=(u2934*(u2886*(u598+u2980*(u4135+u769)+u2886*(u516+u2980*(u748+&
 &u4188)+u4194))+u2885*(u696+u2980*(u2980*(u48+u463)+u699)+u2885*(u718+&
 &u2980*(u48+u3938)+u697))+u2980*(u698+u2980*(u3942+u700))+u915))
  value_4_=value_4_+(c(64)*u4163)
  u916=1.24329546086653774d-1*u3993
  u701=-3.36725853984687303d-1*u895
  u927=-5.05088780977030955d-2*u978
  u702=+u414*(-u3861+u70)
  u70=u4026*(1.17809606699892737d9*u830+u3902)
  u703=1.89408292866386608d-2*u70
  u704=-9.25996098457890084d-2*u4026*(-u3861+1.03626140537006388d8*u830&
 &)
  u705=1.68362926992343652d-1*u938
  u987=1.68362926992343652d-2*u768
  u706=3.78816585732773216d-2*u70
  u70=+u4192*(2.54405069603281883d10*u830+u3872*(2.35d2+u4031))
  u707=u880*(1.34791171629607005d9*u830+u3902)
  u71=-1.68362926992343652d-2*u3958
  u8=3.36725853984687303d-2*u3921
  u756=-5.11402390739243842d-1*u937
  u757=-1.57840244055322173d-1*u937
  u4135=5.05088780977030955d-2*u3921
  u3798=u4135*u2980
  u4146=u8*u2980
  value_5_=value_5_+(c(64)*(u2937*(u2886*(u701+u2980*(u2980*(u71+u3798)&
 &+u705)+u2886*(u871+u2980*(u987+u4146)+u927))+u2885*(u702+u2980*(u2980&
 &*(u756+u4015)+u706)+u2885*(u4174+u2980*(u70+u4002)+u703))+u2980*(u704&
 &+u2980*(u757*u2980+u707))+u916)))
  value_6_=value_6_+(c(64)*u712)
  u871=u873*(u711+u773)
  u711=-5.79186076378022492d-3*u871
  u712=5.6470642446857193d-2*u4109
  u4174=u4026*(u3985-u3819)
  u69=-3.7647094964571462d-2*u4174
  u708=8.47059636702857895d-2*u940
  u602=1.84674518610642669d8*u830
  u709=u4133*(u602+u4047)
  u55=-2.82353212234285965d-1*u938
  u620=5.6470642446857193d-1*u937
  u9=u4026*(u56+u4008)
  u56=1.69411927340571579d-1*u9
  u791=5.95416120348106535d9*u830
  u83=u4026*(u791-9.07496371663624206d2*u4197)
  u4197=-2.82353212234285965d-2*u83
  u3974=-8.47059636702857895d-2*u79
  u79=5.6470642446857193d-2*u394
  u4068=-5.6470642446857193d-2*u3921
  u432=-1.41176606117142983d-1*u937
  u936=1.78823701081714444d-1*u937
  u854=u4068*u2980
  value_7_=value_7_+(c(64)*(u2937*(u2886*(u712+u2980*(u2980*(u79+u854)+&
 &u55)+(u2980*(u620+u854)+u399)*u2886)+u2885*(u69+u2980*(u2980*(u432+u3&
 &887)+u56)+u2885*(u4057+u2980*(u4197+u858)+u708))+u2980*(u709+u2980*(u&
 &936*u2980+u3974))+u711)))
  u711=-1.97010532817858551d-3*u3778
  u712=8.66846344398577624d-2*u3993
  u69=-3.5215632741192216d-2*u4026
  u708=+u69*(u4047+u695)
  u4057=u873*(u829+u3987*(3.0d0+u4031))
  u936=2.16711586099644406d-2*u4057
  u56=u4026*(u564+u4003)
  u4197=-7.04312654823844319d-2*u56
  u79=2.53552555736583955d0*u939
  u4187=3.5215632741192216d-2*u3786
  u709=1.05646898223576648d-1*u945
  u55=-1.05646898223576648d-1*u659
  u979=u4026*(u695+u3862)
  u3968=-7.04312654823844319d-2*u979
  u925=u4026*(u4123+u3982*(8.0d0+pd2))
  u843=5.63450123859075456d-1*u925
  u4031=7.04312654823844319d-2*u3921
  u710=3.5215632741192216d-2*u937
  u4153=(u3837+u2980*(u839+u4189))
  u3943=u4031*u2980
  value_8_=value_8_+(c(64)*(u2885*(u712+u2980*(u2980*(u709+u2980*(u3882&
 &+u843))+u4197)+u2885*(u4153*u2885+u2980*(u79+u2980*(u3943+u55))+u708)&
 &)+u2980*(u936+u2980*(u2980*(u3968+u710*u2980)+u4187))+u711))
  u711=5.60277286089120519d-1*u937
  u712=1.86759095363040173d-1*u937
  value_9_=value_9_+(c(64)*(x*(u2885*(u953+u2980*(u2980*(u712+u4037)+u8&
 &40)+u2885*(u982+u2980*(u711+u4025)+u3888))+u2980*(u944+u2980*(u4134+u&
 &879))+u442)*z))
  u936=1.39307483720682646d-3*u506
  u4197=7.66191160463754555d-3*u4166
  u4187=-7.66191160463754555d-3*u873
  u3968=-u3988*(4.4d1*pd2-2.1d1)
  u843=+u4187*(u3968+u4086)
  u4134=u4026*(u4103-u4003)
  u982=9.96048508602880921d-2*u4134
  u980=u4026*(u4091+u3872)
  u514=2.98814552580864276d-1*u980
  u713=4.98024254301440461d-2*u4130
  u4207=u4026*(1.16748258891785595d8*u830+u3872)
  u515=-4.48221828871296415d-1*u4207
  u390=u4026*(u890-9.07496371663624206d2*u4022)
  u3970=-4.98024254301440461d-2*u390
  u981=u4026*(u4087+u3873)
  u4208=-1.99209701720576184d-1*u981
  u4022=u4026*(u917-4.53748185831812103d2*u802)
  u917=1.99209701720576184d-1*u4022
  u621=4.98024254301440461d-2*u937
  value_1_=value_1_+(c(65)*(u2885*(u4197+u2980*(u2980*(u515+u2980*(u382&
 &4+u917))+u982)+u2885*((u4044+u2980*(u4036+u3793))*u2885+u2980*(u514+u&
 &3970*u2980)+u407))+u2980*(u843+u2980*(u2980*(u4208+u621*u2980)+u713))&
 &+u936))
  u936=1.08257476426928461d8*u870
  u4197=-u3861*(1.7d1*pd2-9.0d0)
  u843=-1.62533689574733305d-2*u873*(u4197+u936)
  u982=1.05646898223576648d-1*u4109
  u514=-1.05646898223576648d-1*u3786
  u713=4.22587592894306592d-1*u939
  u515=u3916*(u4091-u4003)
  u3970=-5.28234491117883239d-1*u938
  u4208=u4026*(u3772+u3872)
  u917=-3.16940694670729944d-1*u4208
  u4209=1.05646898223576648d-1*u394
  u85=u2980*(u4209+u3827)
  value_2_=value_2_+(c(65)*(u2937*(u2886*(u982+u2980*(u85+u3970)+u714*u&
 &2886)+u2885*(u514+u2980*(u3934+u679)+u2885*(u3813+u445+u713))+u2980*(&
 &u515+u2980*(u3941+u917))+u843)))
  u3813=2.63266398353646588d-3*u506
  u514=-2.60633734370110121d-2*u3844
  u843=4.51765139574857544d-1*u3991
  u66=-u4212
  u917=5.79186076378022492d-3*u446
  u713=-1.8823547482285731d-2*u767
  u446=3.59609955500244018d5*u4064
  u445=7.6d1*pd2
  u515=+u563*(-u3988*(u445-1.5d1)+u393)
  u3865=2.58919167960175693d7*u4064
  u3772=-3.38823854681143158d-1*u940
  u563=6.77647709362286316d-1*u937
  u714=-3.7647094964571462d-2*u715
  u715=2.25882569787428772d-1*u946
  u4034=-1.12941284893714386d-1*u937
  u969=5.6470642446857193d-2*u457
  u679=2.25882569787428772d-1*u489
  u716=3.95294497128000351d-1*u72
  u72=-1.8823547482285731d-2*u3995
  u4053=u4026*(u449+u717)
  u717=3.01176759716571696d-1*u4053
  u415=-1.50588379858285848d-1*u3908
  u624=3.7647094964571462d-2*u3921
  u4035=-1.8823547482285731d-2*u937
  u3944=u624*u2980
  u4147=u563*u2980
  value_3_=value_3_+(c(65)*(u2886*(u514+u2980*(u2980*(u3772+u4147)+u386&
 &5)+u2886*((u66+u2980*(u3998+u563))*u2886+u2980*(u3772+u2980*(u3998+u6&
 &79))+u843))+u2885*(u917+u2980*(u2980*(u716+u2980*(u3777+u415))+u714)+&
 &u2885*((u446+u2980*(u3777+u4034))*u2885+u2980*(u715+u2980*(u3944+u72)&
 &)+u713))+u2980*(u515+u2980*(u2980*(u717+u4035*u2980)+u969))+u3813))
  u512=u873*(4.53748185831812103d2*u4004+u829)
  u836=-1.84295880979078208d-2*u512
  u4004=4.64870339950928097d8*u830
  u970=u77*(u4004-u4047)
  u77=u4026*(5.66759729529213707d8*u830+u4023)
  u3807=-3.99307742121336118d-2*u77
  u591=-1.19792322636400836d-1*u3786
  u592=4.79169290545603342d-1*u939
  u593=-1.3310258070711204d-1*u4026*(2.16514952853856922d7*u830+u4003)
  u869=-6.65512903535560198d-2*u938
  u752=2.66205161414224079d-2*u4026
  u3935=u752*(5.08810139206563766d9*u830+u3872*(4.7d1+u3992))
  u49=-1.3310258070711204d-2*u4026*(5.73764625062720842d9*u830+u3872*(5&
 &.3d1+u3960))
  u43=-1.3310258070711204d-2*u3921
  u3945=u3900*u2980
  u389=u2980*(u3945+u527)
  u4148=u43*u2980
  u749=(u2937*(u2886*(u970+u2980*(u2980*(u49+u463)+u869)+u2886*(u516+u2&
 &980*(u3935+u4148)+u3807))+u2885*(u591+u389+u2885*(u718+u798+u592))+u2&
 &980*(u593+u2980*(u3942+u543))+u836))
  value_4_=value_4_+(c(65)*u749)
  u516=-1.45698686820297391d-3*u873*(-u3861*(u3814+3.6d1)+u773)
  u718=+u4192*(6.55912945410213616d8*u830-u402)
  u719=6.06106537172437146d-1*u939
  u391=-2.02035512390812382d-1*u937
  u918=3.15680488110644347d-2*u4130
  u802=-7.57633171465546432d-2*u4026
  u517=+u802*(5.73127816377856558d7*u830+u3873)
  u534=1.89408292866386608d-2*u937
  u989=1.68362926992343652d-2*u4026*(u3985+1.7582742200982719d3*exppd2)
  u786=u880*(u825+u4055)
  u825=-9.4704146433193304d-2*u938
  u4055=5.5211312977733515d9*u830
  u4041=u462*(u4055+u3872*(5.1d1+u3774))
  u722=-3.78816585732773216d-2*u943
  u721=-1.32585805006470626d-1*u937
  u3814=u4026*(9.74317287842356148d8*u830+u3862*(1.8d1+pd2))
  u3810=5.05088780977030955d-2*u3814
  u720=3.15680488110644347d-2*u937
  u3953=(u391+u3849)*u2886
  u7=(u534+u3776)*u2885
  value_5_=value_5_+(c(65)*(u2934*(u2886*(u718+u2980*(u2980*(u721+u4015&
 &)+u786)+u2886*(u3953+u4137+u719))+u2885*(u918+u2980*(u2980*(u3810+u38&
 &18)+u825)+u2885*(u7+u2980*(u4041+u3977)+u517))+u2980*(u989+u2980*(u72&
 &0*u2980+u722))+u516)))
  value_6_=value_6_+(c(65)*u599)
  u918=u3806*(u4086-u3861*(2.4d1+u4093))
  u3806=6.88911213625908387d7*u830
  u989=-1.0352951115257152d-1*u4026*(-u3861+u3806)
  u786=1.69411927340571579d-1*u943
  u4137=-1.69411927340571579d-1*u937
  u721=u886*(u564-u3861)
  u722=+u3969*(u494+u3862)
  u564=2.82353212234285965d-2*u937
  u517=u4026*(-u4003+u449)
  u599=-5.6470642446857193d-2*u517
  u382=8.47059636702857895d-2*u4026*(6.26195206783213646d8*u830+u4023)
  u919=-4.70588687057143275d-2*u938
  u622=9.41177374114286549d-3*u394
  u761=-4.8000046079828614d-1*u937
  u594=9.4117737411428655d-2*u937
  u409=(u4137+u858)*u2886
  u3836=(u564+u3887)*u2885
  u3854=u594*u2980
  value_7_=value_7_+(c(65)*(u2934*(u2886*(u989+u2980*(u2980*(u761+u3887&
 &)+u382)+u2886*(u409+u2980*(u60+u4011)+u786))+u2885*(u721+u2980*(u3854&
 &+u919)+u2885*(u3836+u2980*(u622+u3887)+u722))+u2980*(u599+u617*u4191)&
 &+u918)))
  u918=8.12668447873666522d-3*u792
  u989=-1.21098352171944274d7*u4064
  u919=+u989
  u721=3.16940694670729944d-1*u941
  u622=3.16940694670729944d-1*u943
  u599=u2980*(u55+u839)
  u382=u2980*(u4189+u839)
  value_8_=value_8_+(c(65)*(u2936*((u553+u2980*(u2980*(u4012+u3875)+u85&
 &3))*u2886+u2885*(u3825+u2980*(u599+u721)+(u382+u3837)*u2885)+u2980*(u&
 &919+u2980*(u3940+u622))+u918)))
  u919=u723*(1.19083224069621307d9*u870-u3861*(1.87d2*pd2-6.0d1))
  u622=-8.71542445027520806d-2*u3994
  u60=4.11031060230583996d7*u830
  u786=u4026*(-u3861+u60)
  u722=-1.36956669932896127d-1*u786
  u761=u4026*(u449-u3873)
  u723=1.49407276290432138d-1*u761
  u558=3.73518190726080346d-2*u937
  u724=-2.19130671892633803d0*u744
  u3841=8.32508697332068481d6*u4064
  u725=3.11265158938400288d-1*u938
  u84=3.13946681638092536d9*u830
  u3956=-u3872*(u3774-2.9d1)
  u4178=u4026*(u3956+u84)
  u73=-1.24506063575360115d-2*u4178
  u983=u4026*(1.676929536809284d8*u830+u3872)
  u726=2.24110914435648207d-1*u983
  u3940=-2.61462733508256242d-1*u937
  u385=u4026*(1.19083224069621307d9*u830+u3873*(4.4d1+u4024))
  u74=-1.99209701720576184d-1*u385
  u3822=-2.36561520793184219d-1*u937
  u924=(u558+u455)*u2885
  value_9_=value_9_+(c(65)*(u2934*((u622+u2980*(u2980*(u3940+u4037)+u38&
 &41))*u2886+u2885*(u722+u2980*(u2980*(u74+u456)+u725)+u2885*(u924+u298&
 &0*(u73+u3770)+u723))+u2980*(u724+u2980*(u3822*u2980+u726))+u919)))
  u838=3.4478602220868955d-2*u383
  u837=-1.19525821032345711d0*u895
  u727=1.49407276290432138d-1*u941
  u4170=1.49407276290432138d-1*u945
  u75=-4.98024254301440461d-2*u659
  u751=u2980*(u75+u4036)
  u4211=u2980*(u3793+u4036)
  value_1_=value_1_+(c(66)*(u2936*((u779)*u2886+u2885*(u407+u2980*(u751&
 &+u727)+(u4211+u4044)*u2885)+u2980*(u837+u2980*(u3853+u4170))+u838)))
  u3820=1.40097910670142714d7*u830
  u4170=-5.28234491117883239d-1*u4026*(-u3861+u3820)
  u837=-1.05646898223576648d-1*u4026*(u4047+u3979)
  u779=6.04968250621070811d8*u830
  u728=1.05646898223576648d-1*u4026*(u779+u4008)
  u729=+u69*(u791-9.07496371663624206d2*u4049)
  u764=-8.45175185788613183d-1*u937
  u3946=u764*u2980
  u3947=u3857*u2980
  value_2_=value_2_+(c(66)*(u2934*(u2886*(u4170+u2980*(u3946+u728)+u288&
 &6*(u3893*u2886+u2980*(u729+u3978)+u709))+u2980*(u837+u3947)+u907)))
  u4170=4.34389557283516869d-3*u873
  u837=u4170*(u498-u3861*(3.3d1+u4005))
  u728=-2.07059022305143041d-1*u3994
  u4052=-u446
  u729=-1.12941284893714386d-1*u770
  u3893=-6.77647709362286316d-1*u937
  u565=5.6470642446857193d-2*u941
  u730=5.6470642446857193d-2*u945
  u44=-1.8823547482285731d-2*u659
  u861=u44+u3777
  u3866=u2980*(u861)
  u796=u2980*(u3866+u565)
  u3785=u2980*(u4034+u3777)
  value_3_=value_3_+(c(66)*(u2936*(u2886*(u728+u2980*(u2980*(u458+u3879&
 &)+u4110)+(u2980*(u3893+u4011)+u4212)*u2886)+u2885*(u4052+u796+(u3785+&
 &u446)*u2885)+u2980*(u729+u730*u2980)+u837)))
  u595=-1.19792322636400836d-1*u765
  u596=1.19792322636400836d-1*u945
  u709=-2.79515419484935283d-1*u937
  u597=-1.3310258070711204d-2*u50
  u50=-3.99307742121336118d-2*u413
  u775=(u709+u3851)*u2886
  u4149=u527*u2980
  u4150=u50*u2980
  u4151=u597*u2980
  u4193=(u2934*(u2886*(u595+u4149+u2886*(u775+u4150+u596))+u4151+u816))
  value_4_=value_4_+(c(66)*u4193)
  u920=u426*(1.06347050372335606d9*u870-u3861*(1.67d2*pd2-4.44d2))
  u731=-8.20769269087675302d-2*u4026*(9.65010083986647366d7*u830+u3861)
  u732=u3910*(6.81385292804785018d8*u830+u4023)
  u4005=1.89408292866386608d-2*u3786
  u4009=-7.57633171465546432d-2*u939
  u788=-8.41814634961718258d-3*u4026
  u733=+u788*(u3985-8.96152667017828903d3*exppd2)
  u734=9.4704146433193304d-2*u4026*(3.45999385442928218d8*u830+u4023)
  u76=+u414*(u4055-9.07496371663624206d2*u4013)
  u4055=-1.89408292866386608d-2*u941
  u907=6.31360976221288694d-3*u659
  u865=-6.06106537172437146d-1*u4053
  u4013=u462*(u3799+u3872*(4.1d1+u3992))
  u4053=u2980*(u907+u3776)
  u3799=u532*u2980
  u3800=u534*u2980
  value_5_=value_5_+(c(66)*(u2937*(u2886*(u731+u2980*(u2980*(u4013+u377&
 &6)+u734)+u2886*(u3953+u2980*(u76+u3818)+u732))+u2885*(u4005+u2980*(u3&
 &799+u4055)+u2885*(u7+u4053+u4009))+u2980*(u733+u2980*(u3800+u865))+u9&
 &20)))
  u3953=u3779*(-u4003*(u4019+3.0d0*u3835)+u401)
  u7=-2.79236183301633649d-4*u3953
  u521=u873*(u498-u4003*(2.7d1+u4099))
  u921=9.21479404895391042d-3*u521
  u735=-3.99307742121336118d-2*u770
  u3829=-u922
  u922=6.14319603263594028d-3*u792
  u736=-2.39584645272801671d-1*u765
  u737=2.39584645272801671d-1*u945
  u864=-5.59030838969870566d-1*u937
  value_6_=value_6_+(c(66)*(u2886*(u921+u2980*(u4149+u736)+u2886*((u382&
 &9+u2980*(u3851+u864))*u2886+u2980*(u737+u4150)+u735))+u2980*(u922+u41&
 &51)+u7))
  u763=6.51584335925275303d-3*u873*(u3783-u3861*(4.4d1+u3775))
  u4150=1.5219727568256413d9*u830
  u78=-9.41177374114286549d-3*u4026
  u3990=2.67711429640769141d4*exppd2
  u4049=+u78*(u3990+u4150)
  u923=1.69411927340571579d-1*u951
  u738=2.82353212234285965d-2*u3786
  u826=-1.12941284893714386d-1*u939
  u741=+u3969*(u803+u3985)
  u4149=u4026*(u933+u4023)
  u4151=1.41176606117142983d-1*u4149
  u739=-5.6470642446857193d-2*u971
  u984=-2.82353212234285965d-2*u941
  u740=9.41177374114286549d-3*u659
  u3835=6.83258915450463635d6*u4064
  u3891=-9.31765600373143684d-1*u937
  u4152=u564*u2980
  value_7_=value_7_+(c(66)*(u2937*(u2886*(u4049+u2980*(u2980*(u3891+u38&
 &87)+u4151)+u2886*(u409+u2980*(u739+u4011)+u923))+u2885*(u738+u2980*(u&
 &3852+u984)+u2885*(u3836+u2980*(u740+u3887)+u826))+u2980*(u741+u2980*(&
 &u4152+u3835))+u763)))
  u763=-7.38789498066969565d-4*u3953
  u4049=3.4d1*pd2
  u923=u3912*(u936-u4003*(2.1d1+u4049))
  u738=-3.5215632741192216d-2*u523
  u3836=1.46280320617259974d-1*u383
  u741=-6.33881389341459887d-1*u765
  u4151=6.33881389341459887d-1*u943
  u739=-1.05646898223576648d-1*u767
  u984=u4026*(u3801+u4008)
  u740=3.16940694670729944d-1*u984
  u409=-1.05646898223576648d-1*u413
  value_8_=value_8_+(c(66)*(u2886*(u923+u2980*(u2980*(u740+u3946)+u741)&
 &+u2886*(u4153*u2886+u2980*(u4151+u2980*(u3978+u409))+u738))+u2980*(u3&
 &836+u2980*(u3947+u739))+u763))
  u763=2.87321685173907958d-3*u873
  u923=u763*(u3954-u3861*(1.2d1+u754))
  u738=-1.24506063575360115d-2*u4026*(u931+6.9412146650207072d8*u830)
  u3836=7.13578883427487269d6*u4064
  u409=3.73518190726080346d-2*u3786
  u4151=-1.49407276290432138d-1*u939
  u739=-1.34466548661388924d0*u4026*(-u3988+1.20286084918809401d7*u830)
  u740=1.86759095363040173d-1*u77
  u3891=-3.73518190726080346d-2*u941
  u826=1.24506063575360115d-2*u659
  u741=1.79288731548518566d0*u4026*(u3948+u3982)
  u523=-3.73518190726080346d-2*u4026
  u77=+u523*(7.4697658734580638d9*u830+u3872*(6.9d1+u3992))
  u3982=8.71542445027520806d-2*u3921
  u4054=u2980*(u826+u455)
  u3801=u3982*u2980
  u3948=u544*u2980
  u4153=u3940*u2980
  value_9_=value_9_+(c(66)*(u2937*(u2886*(u738+u2980*(u2980*(u77+u3801)&
 &+u740)+(u2980*(u4105+u456)+u3836)*u2886)+u2885*(u409+u2980*(u3948+u38&
 &91)+u2885*(u924+u4054+u4151))+u2980*(u739+u2980*(u4153+u741))+u923)))
  u924=-2.29857348139126367d-2*u3886
  u4180=5.47826679731584507d-1*u457
  u4154=-u550
  u550=1.49407276290432138d-1*u937
  u3999=-u3855
  u3855=5.47826679731584507d-1*u937
  u863=(u550+u3824)
  u519=u863*u2885
  u3949=u4196*u2980
  u3950=u887*u2980
  value_1_=value_1_+(c(67)*(u2934*(u2885*(u4180+u2980*(u3949+u3999)+u28&
 &85*(u519+u2980*(u3855+u4036)+u4154))+u2980*(u619+u3950)+u924)))
  u924=5.41778965249111015d-3*u4097
  u3999=-u3857
  u932=-3.87371960153114376d-1*u3994
  u3855=u853*u2980
  u4154=u4202*u2980
  value_2_=value_2_+(c(67)*(y*(u2885*(u947+u2980*(u4154+u4155)+u2885*((&
 &u710+u3882)*u2885+u2980*(u530+u839)+u3999))+u2980*(u932+u3855)+u924)*&
 &z))
  u932=(u538+u3961)*u2885
  u4155=u4067*u2980
  value_3_=value_3_+(c(67)*(u2934*(u2885*(u4169+u2980*(u661*u2980+u4157&
 &)+u2885*(u932+u4158+u522))+u2980*(u958+u4155)+u507)))
  u4169=u482*u2885
  u4156=u540*u2980
  u4157=u3845*u2980
  u522=(u2936*(u2886*(u959+u2885*(u2885*(u3863+u4169)+u4088)+(u2885*(u5&
 &67+u4169)+u3924)*u2886)+u2885*(u961+u2980*(u4156+u398)+u2885*((u568+u&
 &3805)*u2885+u4159+u424))+u2980*(u960+u4157)+u438))
  value_4_=value_4_+(c(67)*u522)
  u507=u778*u2885
  u4158=u820*u2980
  u4159=u3842*u2980
  value_5_=value_5_+(c(67)*(u2886*(u508+u2885*(u2885*(u743+u562*u2885)+&
 &u776)+u2886*((u781+u2885*(u507+u562))*u2886+u2885*(u743+u2885*(u507+u&
 &4161))+u774))+u2885*(u905+u2980*(u2980*(u664+u4158)+u663)+u2885*(u288&
 &5*(u667+u2980*(u3936+u65)+(u3881+u747)*u2885)+u2980*(u666+u2980*(u388&
 &1+u64))+u665))+u2980*(u904+u2980*(u4159+u662))+u4))
  value_6_=value_6_+(c(67)*u988)
  u3957=-2.63266398353646587d-4*u3953
  u747=u793*(u3954-u3861*(5.1d1+u754))
  u3953=+u78*(u931+u832)
  u4161=1.12941284893714386d-1*u974
  u4162=-6.58824161880000585d-2*u937
  u781=-1.44796519094505623d-3*u4132
  u40=-1.8823547482285731d-2*u4026
  u743=+u40*(u484+u3858)
  u3954=2.54117891010857368d-1*u3967
  u3967=2.82353212234285965d-2*u878
  u439=-1.69411927340571579d-1*u941
  u3856=2.82353212234285965d-2*u387
  u387=-u45
  u4160=u387*u2980
  value_7_=value_7_+(c(67)*(u2885*(u747+u2980*(u2980*(u439+u3796)+u743)&
 &+u2885*(u2885*(u4161+u62*u2980+(u772+u4162)*u2885)+u2980*(u3954+u2980&
 &*(u4000+u3856))+u3953))+u2980*(u781+u2980*(u4160+u3967))+u3957))
  u4161=-u746
  u4162=-u82
  u62=(u549+u3827)
  value_8_=value_8_+(c(67)*(x*(u2885*(u957+u2980*(u3934+u405)+u2885*(u6&
 &2*u2885+u2980*(u4162+u3875)+u4161))+u2980*(u947+u4114*u2980)+u443)*z)&
 &)
  u4161=u3877*u2980
  u4162=u879*u2980
  value_9_=value_9_+(c(67)*(u2885*(u440+u2980*(u2980*(u4116+u4161)+u782&
 &)+u2885*(u2885*(u3952+u2980*(u4025+u669)+(u4037+u784)*u2885)+u2980*(u&
 &785+u2980*(u4037+u668))+u966))+u2980*(u440+u2980*(u4162+u966))+u954))
  u3954=-u986
  u3967=2.49012127150720231d-1*u937
  value_1_=value_1_+(c(68)*(x*(u2885*(u973+u2980*(u3853+u887)+u2885*(u5&
 &19+u2980*(u3967+u4036)+u3954))+u2980*(u555+u4139))*z))
  u3954=1.97010532817858551d-3*u506
  u3967=-2.16711586099644406d-2*u3913
  u3856=-3.5215632741192216d-2*u3786
  u78=7.04312654823844319d-2*u979
  u743=-8.66846344398577624d-2*u934
  u781=7.04312654823844319d-2*u56
  u440=-1.05646898223576648d-1*u945
  u782=-5.63450123859075456d-1*u925
  u954=u4030*(u695+u4047)
  u785=-u79
  u4139=1.05646898223576648d-1*u659
  u3952=u533*u2980
  u3953=u3825*u2980
  value_2_=value_2_+(c(68)*(u2885*(u3967+u2980*(u2980*(u785+u3952)+u781&
 &)+u2885*(u2885*(u78+u2980*(u3937+u782)+(u3875+u758)*u2885)+u2980*(u44&
 &0+u2980*(u3827+u4139))+u3856))+u2980*(u743+u2980*(u3953+u954))+u3954)&
 &)
  value_3_=value_3_+(c(68)*(x*(u2885*(u384+u2980*(u4141+u518)+u2885*(u9&
 &32+u926+u4190))+u2980*(u537+u4142)+u513)*z))
  u3954=u4084*u2980
  u3967=(u2886*(u441+u2885*(u2885*(u3794+u531*u2885)+u4126)+u2886*((u79&
 &4+u2885*(u4169+u531))*u2886+u2885*(u3794+u2885*(u4169+u448))+u967))+u&
 &2885*(u835+u2980*(u2980*(u579+u3945)+u578)+u2885*(u2885*(u582+u2980*(&
 &u3938+u3911)+(u463+u789)*u2885)+u2980*(u581+u2980*(u463+u46))+u580))+&
 &u2980*(u834+u2980*(u3954+u577))+u524)
  value_4_=value_4_+(c(68)*u3967)
  u3856=u4135*u2885
  value_5_=value_5_+(c(68)*(u2936*(u2886*(u701+u2885*(u2885*(u71+u3856)&
 &+u705)+u2886*((u618+u507)*u2886+u2885*(u987+u8*u2885)+u927))+u2885*(u&
 &704+u2980*(u2980*(u70+u3881)+u706)+u2885*((u757+u4015)*u2885+u2980*(u&
 &756+u4002)+u707))+u2980*(u702+u2980*(u4143+u703))+u916)))
  value_6_=value_6_+(c(68)*u4163)
  u748=5.79186076378022492d-3*u787
  u745=-5.6470642446857193d-2*u4165
  u781=+u3969*(u4047+u602)
  u782=2.82353212234285965d-1*u938
  u3955=-5.6470642446857193d-2*u394
  u785=-1.78823701081714444d-1*u937
  u537=3.7647094964571462d-2*u4174
  u384=-1.69411927340571579d-1*u9
  u598=1.41176606117142983d-1*u937
  u9=-8.47059636702857895d-2*u940
  u932=2.82353212234285965d-2*u83
  u988=u95*u2885
  value_7_=value_7_+(c(68)*(u2936*(u2886*(u745+u2885*(u2885*(u3955+u988&
 &)+u782)+(u2885*(u4186+u988)+u985)*u2886)+u2885*(u781+u2980*(u2980*(u9&
 &32+u4000)+u384)+u2885*((u785+u452)*u2885+u2980*(u598+u854)+u609))+u29&
 &80*(u537+u2980*(u4144+u9))+u748)))
  u3955=1.95040427489679965d-1*u873*(u829+u783)
  u745=-1.40862530964768864d-1*u4026*(-u3988+u487)
  u781=3.16940694670729944d-1*u930
  u782=-7.0431265482384432d-1*u895
  u930=3.5215632741192216d-1*u938
  u785=-3.5215632741192216d-2*u3958
  u537=-1.05646898223576648d-1*u978
  u384=3.5215632741192216d-2*u768
  u9=(u4104+u839)
  u932=u9*u2885
  value_8_=value_8_+(c(68)*(u2934*(u2885*(u745+u2980*(u2980*(u384+u3882&
 &)+u930)+u2885*(u932+u2980*(u785+u3943)+u781))+u2980*(u782+u2980*(u385&
 &0+u537))+u3955)))
  value_9_=value_9_+(c(68)*(y*(u2885*(u944+u2980*(u2980*(u711+u4037)+u8&
 &40)+u2885*((u784+u4037)*u2885+u2980*(u712+u4025)+u879))+u2980*(u953+u&
 &2980*(u4145+u3888))+u442)*z))
  u930=u4026*(u4087+u4003)
  u745=-9.96048508602880921d-2*u930
  u781=1.49407276290432138d-1*u952
  u782=9.96048508602880921d-2*u930
  u930=u4026*(2.27340700496549768d9*u830-9.07496371663624206d2*u459)
  u785=-4.98024254301440461d-2*u930
  u537=-1.49407276290432138d-1*u952
  u384=4.98024254301440461d-2*u930
  u930=(u4201+u4036)
  u787=u930*u2885
  u3955=u550*u2980
  value_1_=value_1_+(c(69)*(u2934*(u2885*(u745+u4191*(u3824+u384)+u2885&
 &*(u787+u785*u2980+u781))+u2980*(u782+u2980*(u3955+u537)))))
  u745=4.87601068724199913d-2*u383
  u781=-7.04312654823844319d-2*u4026*(-u4003+u695)
  u782=u4026*(u3985+u3872)
  u785=1.05646898223576648d-1*u782
  u537=3.16940694670729944d-1*u3786
  u384=-3.16940694670729944d-1*u941
  u748=-1.26776277868291978d0*u939
  u4194=+u748
  u459=u2980*(u4139+u3827)
  u772=u90*u2885
  u83=u2885*(u51+u772)
  value_2_=value_2_+(c(69)*(u2936*(u2886*(u605+u2885*(u83+u606)+(u2885*&
 &(u450+u772)+u4042)*u2886)+u2885*(u781+u2980*(u459+u384)+u2885*(u893*u&
 &2885+u750+u785))+u2980*(u537+u2980*(u3941+u4194))+u745)))
  u745=5.21267468740220243d-2*u383
  u781=3.38823854681143158d-1*u3786
  u785=-1.35529541872457263d0*u939
  u537=3.38823854681143158d-1*u937
  u4194=-9.03530279149715087d-1*u744
  u744=5.6470642446857193d-2*u951
  u790=-1.12941284893714386d-1*u3933
  u954=-3.38823854681143158d-1*u941
  u985=1.12941284893714386d-1*u659
  u743=1.8823547482285731d-1*u938
  u78=-1.8823547482285731d-2*u971
  u518=(u537+u3998)*u2886
  u750=(u3884+u3777)*u2885
  value_3_=value_3_+(c(69)*(u2934*(u2886*(u781+u2980*(u4147+u954)+u2886&
 &*(u518+u2980*(u985+u3998)+u785))+u2885*(u4194+u2980*(u2980*(u78+u3777&
 &)+u743)+u2885*(u750+u2980*(u78+u3944)+u744))+u2980*(u790+u2980*(u3980&
 &+u615))+u745)))
  u4147=u86*u2885
  u781=(u2936*(u2886*(u970+u2885*(u2885*(u49+u4147)+u869)+u2886*((u557+&
 &u4169)*u2886+u2885*(u3935+u43*u2885)+u3807))+u2885*(u593+u4185+u2885*&
 &(u4182*u2885+u755+u543))+u2980*(u591+u2980*(u3942+u592))+u836))
  value_4_=value_4_+(c(69)*u781)
  u785=-1.64830837614881897d-3*u3778
  u3942=3.37508602978071084d8*u870
  u744=1.29509943840264347d-3*u873*(u3942-u3861*(u4029-1.2d1))
  u790=-1.68362926992343652d-2*u770
  u954=u4026*(u494+u3872)
  u985=-1.68362926992343652d-2*u954
  u743=1.68362926992343652d-2*u937
  u78=u873*(u3817+u766)
  u4194=3.88529831520793042d-3*u78
  u745=-1.15792172047327682d7*u4064
  u742=1.51526634293109286d-1*u940
  u784=-3.78816585732773216d-2*u517
  u4185=-6.31360976221288694d-3*u954
  u926=2.59019887680528695d-3*u4166
  u4165=-3.36725853984687303d-2*u56
  u4190=5.05088780977030955d-2*u945
  u746=2.69380683187749843d-1*u925
  u925=u3905*(u491-u4003)
  u4186=-1.70467463579747947d-1*u4026*(1.45050867107976042d8*u830+u3872&
 &)
  u519=u3789*(2.4899219578193546d9*u830+u3872*(2.3d1+pd2))
  u986=-4.20907317480859129d-3*u3933
  u79=-5.05088780977030955d-2*u659
  u747=1.26272195244257739d-2*u3995
  u4163=u416*u2980
  value_5_=value_5_+(c(69)*(u2886*(u744+u2980*(u2980*(u616+u4163)+u4165&
 &)+u2885*(u2885*(u742+u416*u2885)+u745)+u2886*(u2886*(u985+u2980*(u414&
 &6+u746)+u2885*(u3856+u416)+(u4168+u743)*u2886)+u2885*(u742+u2885*(u38&
 &56+u3889))+u2980*(u4190+u2980*(u3798+u79))+u790))+u2885*(u4194+u2980*&
 &(u2980*(u4186+u2980*(u3776+u519))+u925)+u2885*(u2885*(u4185+u2980*(u3&
 &977+u519)+(u3776+u589)*u2885)+u2980*(u4186+u2980*(u3977+u747))+u784))&
 &+u2980*(u926+u2980*(u2980*(u4185+u4138)+u986))+u785))
  value_6_=value_6_+(c(69)*u749)
  u790=-u3861*(u4093-2.0d0)
  u744=u4170*(u4086+u790)
  u746=-1.29459583980087847d7*u4064
  u985=+u746
  u743=1.69411927340571579d-1*u940
  u4138=+u4164*(-u3975+u417)
  u745=-1.12941284893714386d-1*u489
  u4190=-9.41177374114286549d-3*u954
  u926=+u3859*(u790+u4086)
  u4170=-u746
  u790=-1.69411927340571579d-1*u940
  u746=1.8823547482285731d-2*u4043
  u4146=u4026*(u417-u3975)
  u986=3.7647094964571462d-1*u4146
  u785=1.12941284893714386d-1*u489
  u8=9.41177374114286549d-3*u954
  u925=-1.8823547482285731d-2*u4043
  u742=-9.41177374114286549d-3*u937
  u4164=u537*u2980
  value_7_=value_7_+(c(69)*(u2886*(u2980*(u2980*(u790+u4164)+u4170)+u28&
 &85*(u2885*(u743+u762*u2885)+u985)+u2886*((u2980*(u854+u537)+u2885*(u9&
 &88+u762))*u2886+u2885*(u743+u2885*(u988+u745))+u2980*(u790+u2980*(u85&
 &4+u785))))+u2885*(u744+u4191*(u2980*(u925+u452)+u598)+u2885*(u2885*(u&
 &4190+u2980*(u3887+u746)+(u3887+u617)*u2885)+u2980*(u432+u93*u4191)+u4&
 &138))+u2980*(u926+u2980*(u2980*(u8+u742*u2980)+u986))))
  u8=-4.87601068724199913d-2*u512
  u985=3.5215632741192216d-2*u4109
  u925=-3.16940694670729944d-1*u3786
  u744=-u748
  u4190=u4026*(u695-u4003)
  u926=7.04312654823844319d-2*u4190
  u4138=-1.7607816370596108d-1*u938
  u745=-1.05646898223576648d-1*u782
  u617=3.5215632741192216d-2*u394
  u743=u2980*(u617+u3882)
  value_8_=value_8_+(c(69)*(u2937*(u2886*(u985+u2980*(u743+u4138)+u795*&
 &u2886)+u2885*(u925+u2980*(u4140+u721)+u2885*(u932+u599+u744))+u2980*(&
 &u926+u2980*(u3850+u745))+u8)))
  u8=-1.39307483720682646d-3*u3778
  u925=7.66191160463754555d-3*u78
  u744=-7.47036381452160691d-2*u517
  u517=-1.24506063575360115d-2*u954
  u745=1.24506063575360115d-2*u937
  u926=-4.98024254301440461d-2*u56
  u746=2.61462733508256242d-1*u4026*(1.41007637362805978d8*u830+u3872)
  u986=-2.4901212715072023d-2*u4026*(-u3872*(pd2-1.0d0)+u4123)
  u56=-7.47036381452160691d-2*u659
  value_9_=value_9_+(c(69)*(u2885*(u925+u2980*(u2980*(u746+u2980*(u455+&
 &u986))+u926)+u2885*(u2885*(u517+u2980*(u3770+u986)+(u455+u745)*u2885)&
 &+u2980*(u746+u2980*(u3770+u56))+u744))+u2980*(u925+u2980*(u2980*(u517&
 &+u745*u2980)+u744))+u8))
  u926=-2.29857348139126367d-2*u3844
  u986=4.98024254301440461d-2*u4109
  u517=-1.49407276290432138d-1*u3786
  u782=5.97629105161728553d-1*u939
  u4170=1.99209701720576184d-1*u4026
  u785=u4170*(u4087-u3988)
  u745=-2.49012127150720231d-1*u938
  u746=-1.49407276290432138d-1*u943
  u4109=4.98024254301440461d-2*u394
  u744=u2980*(u4109+u3824)
  value_1_=value_1_+(c(70)*(u2937*(u2886*(u986+u2980*(u744+u745)+u3959*&
 &u2886)+u2885*(u517+u2980*(u3853+u727)+u2885*(u787+u751+u782))+u2980*(&
 &u785+u2980*(u3955+u746))+u926)))
  u926=2.95515799226787826d-3*u506
  u517=-1.62533689574733305d-2*u3844
  u782=2.81725061929537728d-1*u3991
  u785=-1.34553724635493638d6*u4064
  u75=+u785
  u748=8.12668447873666522d-3*u78
  u790=-8.0732234781296183d6*u4064
  u749=+u790
  u747=1.05646898223576648d-1*u940
  u4165=-1.40862530964768864d-1*u895
  u727=-7.04312654823844319d-2*u489
  u4185=-2.43800534362099957d-2*u3844
  u795=2.42196704343888549d7*u4064
  u3853=-3.16940694670729944d-1*u940
  u925=4.22587592894306592d-1*u3991
  u4186=2.11293796447153296d-1*u489
  value_2_=value_2_+(c(70)*(u2886*(u517+u2980*(u2980*(u3853+u3952)+u795&
 &)+u2885*(u2885*(u747+u3833*u2885)+u749)+u2886*((u75+u2980*(u3827+u533&
 &)+u2885*(u772+u3833))*u2886+u2885*(u747+u2885*(u772+u727))+u2980*(u38&
 &53+u2980*(u3827+u4186))+u782))+u2885*(u748+u2885*(u4114*u2885+u4165))&
 &+u2980*(u4185+u2980*(u3953+u925))+u926))
  u926=u873*(-u4003*(u3919-3.3d1)+u80)
  u517=-8.68779114567033738d-3*u926
  u782=u759*(1.07620667742064176d9*u830+u3819)
  u75=+u57*(3.20527038048356815d8*u830+u4008)
  u748=-5.6470642446857193d-2*u3786
  u749=2.25882569787428772d-1*u939
  u727=8.28236089220572163d-1*u4026*(9.84158876608440553d6*u830-u3988)
  u925=-4.70588687057143275d-1*u938
  u4185=3.01176759716571696d-1*u385
  u954=+u3969*(-u3872+u3979)
  u4186=u4026*(u84+u3956)
  u84=1.8823547482285731d-2*u4186
  u57=-9.4117737411428655d-2*u3921
  u3956=u4034*u2980
  u747=u2980*(u3956+u565)
  u4165=u57*u2980
  value_3_=value_3_+(c(70)*(u2937*(u2886*(u782+u2980*(u2980*(u84+u3777)&
 &+u925)+u2886*(u518+u2980*(u4185+u4165)+u75))+u2885*(u748+u747+u2885*(&
 &u750+u3866+u749))+u2980*(u727+u2980*(u3980+u954))+u517)))
  u3919=3.72314911068844866d-4*u506
  u518=-2.04773201087864676d-3*u926
  u3853=5.32410322828448158d-2*u4190
  u4140=+u447*(3.37508602978071084d8*u830+u4008)
  u750=5.32410322828448158d-2*u937
  u926=9.21479404895391042d-3*u78
  u4190=-9.15417497219094502d6*u4064
  u599=1.19792322636400836d-1*u940
  u751=-1.59723096848534447d-1*u895
  u4194=-7.98615484242672237d-2*u489
  u519=-1.02386600543932338d-3*u3927
  u932=u4026*(u487-u3861)
  u787=5.32410322828448158d-2*u932
  u3927=-2.79515419484935283d-1*u4026*(2.4107757355576506d8*u830+u4008)
  u784=u752*(6.38719110918877919d9*u830+u3872*(5.9d1+u3992))
  u752=-1.06482064565689632d-1*u753
  u753=u490*(u386+u4008)
  u386=-7.98615484242672237d-2*u4026*(u801+u3872*(5.0d0+pd2))
  u793=(u2886*(u518+u2980*(u2980*(u753+u3945)+u787)+u2885*(u2885*(u599+&
 &u3900*u2885)+u4190)+u2886*(u2886*(u4140+u2980*(u4148+u784)+u2885*(u41&
 &47+u3900)+(u4188+u750)*u2886)+u2885*(u599+u2885*(u4147+u4194))+u2980*&
 &(u3927+u2980*(u463+u386))+u3853))+u2885*(u926+u2885*(u4084*u2885+u751&
 &))+u2980*(u519+u2980*(u3954+u752))+u3919)
  value_4_=value_4_+(c(70)*u793)
  u769=1.61887429800330434d-4*u873
  u927=u769*(4.54044592308235251d9*u870-u3861*(7.13d2*pd2-4.44d2))
  u754=+u929*(2.05604646705039859d2*exppd2+u695)
  u755=u880*(1.36064788999335575d9*u830+u3902)
  u756=-9.4704146433193304d-2*u4026*(-u3861+4.20293732010428142d7*u830)
  u757=2.20976341677451043d-1*u938
  u80=+u846*(u3843+u3872*(1.19d2+u4016))
  u758=u81*(u3979-u3873)
  u81=+u4192*(-u3872*(u3774-4.9d1)+u3971)
  u3843=4.41952683354902085d-2*u3921
  u4192=2.10453658740429565d-3*u4026*(u436+1.95111719907679204d4*exppd2&
 &)
  u3969=-6.23187037060270512d6*u4064
  u4148=-5.05088780977030955d-2*u937
  u55=-3.78816585732773216d-2*u941
  u759=7.57633171465546432d-2*u937
  u721=+u828*(u760+u4008)
  u760=4.10384634543837651d-1*u937
  u783=1.26272195244257739d-2*u659
  u436=u2980*(u783+u3818)
  u3945=u4100*u2885
  u4166=u534*u2885
  value_5_=value_5_+(c(70)*(u2936*(u2886*(u754+u2980*(u2980*(u760+u4015&
 &)+u3969)+u2885*(u2885*(u81+u3945)+u757)+u2886*((u600+u3856)*u2886+u28&
 &85*(u80+u3843*u2885)+u2980*(u4148+u3795)+u755))+u2885*(u756+u2980*(u4&
 &36+u55)+u2885*(u4166+u2980*(u759+u3818)+u758))+u2980*(u4192+u2980*(u3&
 &800+u721))+u927)))
  value_6_=value_6_+(c(70)*u4193)
  u4193=1.41176606117142983d-1*u457
  u600=-5.6470642446857193d-1*u4026*(-u3988+u3979)
  u987=-1.0352951115257152d-1*u786
  u786=2.35294343528571638d-1*u938
  u45=-1.50588379858285848d-1*u385
  u768=1.12941284893714386d-1*u761
  u761=-9.41177374114286549d-3*u4178
  u4178=4.70588687057143275d-2*u3921
  u4191=-2.54117891010857368d-1*u3994
  u3856=2.69707466625183014d6*u4064
  u4188=-4.23529818351428947d-1*u937
  u3957=u888*u2885
  value_7_=value_7_+(c(70)*(u2936*(u2886*(u600+u2980*(u2980*(u4188+u388&
 &7)+u4094)+u2885*(u2885*(u761+u3957)+u786)+u2886*((u4137+u988)*u2886+u&
 &2885*(u45+u4178*u2885)+u2980*(u762+u858)+u609))+u2885*(u987+u2885*(u5&
 &64*u2885+u768))+u2980*(u4191+u2980*(u4152+u3856))+u4193)))
  u4152=u4026*(u4123+u3861)
  u4191=-1.05646898223576648d-1*u4152
  u987=u4026*(1.3797521505392843d8*u830+u3872)
  u4188=3.16940694670729944d-1*u987
  u45=2.65336952026785443d8*u830
  u564=3.16940694670729944d-1*u4026*(u45+u4008)
  u786=+u69*(9.20188549628891917d9*u830-9.07496371663624206d2*u4092)
  u4178=u9*u2886
  value_8_=value_8_+(c(70)*(u2934*(u2886*(u4191+u2980*(u3946+u564)+u288&
 &6*(u4178+u2980*(u786+u3978)+u4188))+u2980*(u708+u3947)+u918)))
  u4191=u763*(u829-u3861*(pd2-1.2d1))
  u4188=u4170*(u449+u3987)
  u3947=-u3838
  u786=1.24506063575360115d-2*u4130
  u564=-6.22530317876800576d-2*u938
  u761=1.24506063575360115d-1*u937
  u4193=1.24506063575360115d-2*u394
  u4170=-1.12055457217824104d-1*u4026
  u762=+u4170*(u695+u3861)
  u763=2.24110914435648207d-1*u941
  u4130=-4.48221828871296415d-1*u937
  u764=1.12055457217824104d-1*u940
  u3946=u2980*(u56+u456)
  u3857=u4059*u2885
  value_9_=value_9_+(c(70)*(u2936*(u2886*(u4188+u2980*(u2980*(u818+u403&
 &7)+u4083)+u2885*(u2885*(u4193+u3857)+u564)+(u2885*(u761+u3857)+u3947)&
 &*u2886)+u2885*(u786+u2980*(u3946+u763)+u2885*(u558*u2885+u2980*(u4130&
 &+u456)+u4151))+u2980*(u762+u2980*(u4153+u764))+u4191)))
  u600=-1.49407276290432138d-1*u765
  u765=1.49407276290432138d-1*u943
  u9=-4.98024254301440461d-2*u767
  u767=1.49407276290432138d-1*u984
  u988=-4.98024254301440461d-2*u413
  u768=1.90287702247329938d6*u4064
  u766=-5.97629105161728553d-1*u937
  u4153=9.96048508602880921d-2*u3921
  value_1_=value_1_+(c(71)*(u2934*(u2886*(u600+u2980*(u766*u2980+u767)+&
 &u2886*(u930*u2886+u2980*(u988+u4153*u2980)+u765))+u2980*(u9+u768*u298&
 &0)+u838)))
  u708=1.21900267181049978d-1*u383
  u918=3.3659887628540782d7*u830
  u767=-2.46509429188345512d-1*u4026*(-u3861+u918)
  u988=u4026*(u602+u3872)
  u768=1.05646898223576648d-1*u988
  u766=-6.33881389341459887d-1*u3994
  u4153=4.03661173906480915d7*u4064
  value_2_=value_2_+(c(71)*(u2936*(u2886*(u767+u2980*(u4154+u4153)+u288&
 &5*(u450*u2885+u606)+u2886*((u893+u772)*u2886+u83+u382+u768))+u2885*(u&
 &605+u4042*u2885)+u2980*(u766+u3855)+u708)))
  u767=-8.47059636702857895d-1*u4026*(-u3861+u601)
  u768=u4133*(u525+u3872)
  u838=-1.8823547482285731d-2*u770
  u600=u4133*(u525+u4008)
  u601=+u40*(1.13670350248274884d10*u830+u3872*(1.05d2+u3774))
  u4092=1.31764832376000117d-1*u3921
  value_3_=value_3_+(c(71)*(u2934*(u2886*(u767+u600*u2980+u2886*((u4089&
 &+u4092*u2980)*u2886+u601*u2980+u768))+u838*u2980+u837)))
  u838=1.02386600543932338d-3*u873*(8.08747029777642031d8*u870-u3861*(7&
 &.5d1+1.27d2*pd2))
  u600=-3.32756451767780099d-1*u928
  u601=u490*(2.18637648470071205d8*u830+u3872)
  u4092=-1.73033354919245651d-1*u937
  u602=-1.3310258070711204d-1*u3994
  u51=2.03426110493132112d7*u4064
  u4154=-8.78477032666939461d-1*u937
  u525=(u2936*(u2886*(u600+u2980*(u4156+u51)+u2885*(u811*u2885+u542)+u2&
 &886*((u4092+u771+u4147)*u2886+u2885*(u39+u4147)+u2980*(u4154+u3805)+u&
 &601))+u2885*(u541+u472*u2885)+u2980*(u602+u4157)+u838))
  value_4_=value_4_+(c(71)*u525)
  u9=-1.47170390727573122d-5*u3779*(-u4003*(3.58d2*u4027+3.0d0*u4063)+1&
 &.13988754590707026d9*u3898)
  u928=u769*(4.96073965509278065d9*u870-u4003*(1.558d3*pd2-2.49d2))
  u769=+u4167*(3.72073512382085924d4*exppd2+2.64275604218678301d9*u830)
  u770=u462*(1.8446224904902124d9*u830+u3902)
  u382=-6.73451707969374607d-2*u937
  u83=-1.45698686820297391d-3*u3844
  u930=1.44740215059159603d6*u4064
  u383=-1.89408292866386608d-2*u940
  u929=2.52544390488515477d-2*u3991
  u771=1.26272195244257739d-2*u489
  u772=1.13321200860231304d-3*u873*(u773-u4003*(4.6d1*pd2-1.23d2))
  u773=-2.73589756362558434d-1*u4026*(3.59551980531065567d7*u830+u3861)
  u774=1.70467463579747947d-1*u4026*(4.75247962963531261d8*u830+u4023)
  u82=+u788*(2.32753574317896191d10*u830+u3872*(2.15d2+u844))
  u3902=6.31360976221288694d-3*u3786
  u788=-7.57633171465546432d-2*u980
  u520=u462*(u3848-9.07496371663624206d2*u403)
  u3858=u4085*u2980
  u4167=u532*u2885
  u4168=u4085*u2885
  value_5_=value_5_+(c(71)*(u2886*(u928+u2980*(u2980*(u788+u3799)+u773)&
 &+u2885*(u2885*(u383+u4167)+u930)+u2886*(u2886*(u770+u2980*(u3818+u82)&
 &+u2885*(u3945+u532)+(u3849+u382)*u2886)+u2885*(u383+u2885*(u3945+u771&
 &))+u2980*(u774+u2980*(u3776+u520))+u769))+u2885*(u83+u2885*(u4168+u92&
 &9))+u2980*(u772+u2980*(u3858+u3902))+u9))
  u3849=3.12036255583499681d8*u870
  u383=u873*(u3849-u3861*(7.5d1+u777))
  u929=3.07159801631797014d-3*u383
  u771=+u492*(u931+u4124)
  u772=u490*(u4124+u3872)
  u773=-1.3310258070711204d-2*u892
  u774=1.99653871060668059d-1*u977
  u82=+u492*(u3971-9.07496371663624206d2*u3826)
  value_6_=value_6_+(c(71)*(u2937*(u2886*(u771+u774*u2980+u2886*(u775+u&
 &82*u2980+u772))+u773*u2980+u929)))
  u4063=1.46465997518785565d8*u3898
  u403=4.6d1*u4027
  u520=-1.9744979876523494d-4*u3779*(-u4003*(u403+3.0d0*u3904)+u4063)
  u492=5.41287382134642304d8*u870
  u3902=1.7d2*pd2
  u83=2.17194778641758435d-3*u873*(u492-u4003*(1.83d2+u3902))
  u3826=1.04618569656275403d8*u830
  u930=-6.58824161880000585d-2*u4026*(u3819+u3826)
  u788=4.26661818859070993d8*u830
  u3904=u4026*(u788+u4008)
  u931=2.82353212234285965d-2*u3904
  u775=1.08597389320879217d-2*u521
  u977=3.88453297767213889d8*u830
  u606=u4026*(u4050+u977)
  u521=-5.6470642446857193d-2*u606
  u777=4.23529818351428947d-1*u4026*(3.24772429280785383d8*u830+u4008)
  u778=-4.51765139574857544d-1*u3908
  u776=-8.47059636702857895d-2*u3994
  u490=2.69707466625183014d7*u4064
  u605=-1.55294266728857281d0*u937
  u4169=u3823*u2980
  value_7_=value_7_+(c(71)*(u2886*(u83+u2980*(u2980*(u490+u3852)+u521)+&
 &u2885*(u2885*(u4173+u536*u2885)+u4212)+u2886*(u2886*(u931+u2980*(u401&
 &1+u778)+u2885*(u3957+u536)+(u858+u3884)*u2886)+u2885*(u4173+u2885*(u3&
 &957+u872))+u2980*(u777+u2980*(u3887+u605))+u930))+u2885*(u511+u2885*(&
 &u3823*u2885+u4077))+u2980*(u775+u2980*(u4169+u776))+u520))
  u520=u3912*(u3792-u3861*(1.5d1+u823))
  u511=5.04548419546318166d7*u830
  u930=-4.57803225635498807d-1*u4026*(-u3861+u511)
  u872=u4026*(1.76183736145785534d8*u830+u3872)
  u775=3.16940694670729944d-1*u872
  u521=-3.5215632741192216d-2*u606
  u490=u4026*(3.03545473118642547d8*u830+u4008)
  u3957=5.28234491117883239d-1*u490
  u4077=u4026*(1.11505200719736315d10*u830+u3872*(1.03d2+u3960))
  u605=-3.5215632741192216d-2*u4077
  u606=-1.40862530964768864d0*u937
  u4157=+u606
  value_8_=value_8_+(c(71)*(u2937*(u2886*(u930+u2980*(u4157*u2980+u3957&
 &)+u2886*(u4178+u2980*(u605+u3978)+u775))+u2980*(u521+u3808*u2980)+u52&
 &0)))
  u520=-2.87321685173907958d-3*u3886
  u930=u4171*(5.73127816377856558d7*u870-u4003*(5.0d0+1.8d1*pd2))
  u775=+u523*(u4047+u4087)
  u521=-2.87321685173907958d-3*u3844
  u3957=2.85431553370994907d6*u4064
  u605=-3.73518190726080346d-2*u940
  u4157=4.98024254301440461d-2*u3991
  u609=2.4901212715072023d-2*u489
  u931=u4171*(u4086-u4003*(9.0d0+2.2d1*pd2))
  u4156=u4026*(-u3861+u494)
  u776=-2.24110914435648207d-1*u4156
  u777=1.12055457217824104d-1*u4026*(u779+u4023)
  u778=+u4170*(-u3861+u3979)
  u779=4.48221828871296415d-1*u943
  u83=+u523*(7.03673596775034996d9*u830+u3872*(6.5d1+u3992))
  u3908=1.66501739466413696d6*u4064
  u844=-5.22925467016512484d-1*u937
  u4170=u844*u2980
  u4171=u3908*u2980
  value_9_=value_9_+(c(71)*(u2886*(u930+u2980*(u2980*(u779+u4170)+u776)&
 &+u2885*(u2885*(u605+u544*u2885)+u3957)+u2886*((u3838+u2980*(u456+u413&
 &0)+u2885*(u3857+u544))*u2886+u2885*(u605+u2885*(u3857+u609))+u2980*(u&
 &777+u2980*(u3801+u83))+u775))+u2885*(u521+u2885*(u404*u2885+u4157))+u&
 &2980*(u931+u2980*(u4171+u778))+u520))
  u3857=1.39307483720682646d-3*u4172
  u523=-7.27881602440566827d-2*u3886
  u4178=-u3920
  u4172=-3.83095580231877278d-3*u3886
  u4173=4.98024254301440461d-1*u457
  u4174=-u4175
  u4175=-u4080
  u3958=u4113*u2980
  u3959=u4105*u2980
  u3960=u555*u2980
  value_1_=value_1_+(c(72)*(u2885*(u523+u2980*(u3958+u4173)+u2885*(u288&
 &5*(u4178+u2980*(u4036+u4175)+(u3824+u621)*u2885)+u2980*(u4174+u3959)+&
 &u4180))+u2980*(u4172+u3960)+u3857))
  u4175=6.69097022082652103d-1*u937
  u3920=u3883*u2885
  u4080=u839+u3920
  u4172=u408*u2980
  u4173=u4200*u2980
  u4174=u554*u2980
  value_2_=value_2_+(c(72)*(x*(u2885*(u947+u4172+u2885*(u2885*(u4175+u4&
 &080)+u4173+u405))+u4174+u924)*z))
  u4175=u955*u2980
  value_3_=value_3_+(c(72)*(u2885*(u505+u2980*(u4078*u2980+u996)+u2885*&
 &(u2885*(u3868+u4177+(u3961+u641)*u2885)+u2980*(u848+u640*u2980)+u958)&
 &)+u2980*(u504+u4175)+u4176))
  u641=u3932*u2885
  u4101=u3805+u641
  u848=u92*u2885
  u4176=u4014*u2980
  u4177=u535*u2980
  u4178=u942*u2980
  u504=(u2937*((u551+u2885*(u2885*(u4062+u848)+u4107))*u2886+u2885*(u95&
 &6+u4176+u2885*(u2885*(u539+u4101)+u4177+u464))+u4178+u437))
  value_4_=value_4_+(c(72)*u504)
  u640=u94*u2885
  u505=u4002+u640
  u4001=u3864*u2885
  u3805=u98*u2885
  value_5_=value_5_+(c(72)*(u2934*(u2886*(u997+u2885*(u2885*(u646+u4001&
 &)+u3930)+(u2885*(u3839+u3805)+u3922)*u2886)+u2885*(u643+u2980*(u3914*&
 &u2980+u644)+u2885*(u2885*(u4122+u505)+u2980*(u61+u3881)+u645))+u2980*&
 &(u642+u3923*u2980)+u3880)))
  value_6_=value_6_+(c(72)*u522)
  u522=-2.17194778641758435d-2*u4132
  u3961=3.95294497128000351d-1*u457
  u4180=2.44706117269714503d-1*u3962
  u3962=-8.47059636702857895d-2*u998
  u523=-u4179
  u398=1.65647217844114433d0*u4026*(u883-u3975)
  u4179=-2.82353212234285965d-1*u972
  u3924=-u4079
  u4079=6.58824161880000585d-1*u937
  u998=u93*u2885
  u883=u854+u998
  u3863=u4068*u2885
  value_7_=value_7_+(c(72)*(u2934*((u3961+u2885*(u2885*(u523+u3863)+u40&
 &78))*u2886+u2885*(u4180+u2980*(u4079*u2980+u4179)+u2885*(u2885*(u4034&
 &+u883)+u2980*(u4033+u4000)+u3962))+u2980*(u398+u3924*u2980)+u522)))
  u522=-4.87601068724199913d-2*u3886
  u4180=6.69097022082652103d-1*u457
  u4033=-u3808
  u523=5.28234491117883239d-1*u457
  u3961=-u3909
  u3962=-u4051
  u4179=u4042*u2980
  value_8_=value_8_+(c(72)*(y*(u2885*(u4180+u2980*(u3939+u3961)+u2885*(&
 &(u530+u3827)*u2885+u2980*(u3962+u3875)+u4033))+u2980*(u523+u4179)+u52&
 &2)*z))
  u3909=u87*u2885
  u4051=u4025+u3909
  u3961=u855*u2980
  u3962=u3838*u2980
  value_9_=value_9_+(c(72)*(u2934*(u2885*(u950+u2980*(u3961+u4106)+u288&
 &5*(u2885*(u499+u4051)+u2980*(u647+u4037)+u3838))+u2980*(u950+u3962)+u&
 &451)))
  u523=-1.53238232092750911d-2*u3886
  u647=-u4181
  value_1_=value_1_+(c(73)*(y*(u2885*(u953+u2980*(u3949+u4116)+u2885*((&
 &u621+u3824)*u2885+u2980*(u560+u4036)+u647))+u2980*(u944+u3950)+u523)*&
 &z))
  u647=-1.35444741312277754d-2*u871
  u4106=2.46509429188345512d-1*u457
  u4180=7.0431265482384432d-1*u4134
  u621=-u479
  u4181=-3.16940694670729944d-1*u780
  u780=7.39528287565036535d-1*u937
  u451=u2980*(u422*u2980+u3970)
  u871=u2980*(u982+u874*u2980)
  value_2_=value_2_+(c(73)*(u2934*((u4106+u2885*(u2885*(u780+u3920)+u62&
 &1))*u2886+u2885*(u4180+u451+u2885*((u545+u3827)*u2885+u85+u4181))+u87&
 &1+u647)))
  u647=u406*u2885
  u4180=u576*u2980
  u4181=u831*u2980
  value_3_=value_3_+(c(73)*(u2936*(u2886*(u4032+u2885*(u2885*(u4061+u64&
 &7)+u4045)+(u2885*(u4073+u647)+u4090)*u2886)+u2885*(u3901+u2980*(u4180&
 &+u894)+u2885*((u671+u3879)*u2885+u3869+u3867))+u2980*(u4098+u4181)+u5&
 &09)))
  u3867=u2980*(u811*u2980+u542)
  u3901=u2980*(u541+u472*u2980)
  u4032=(u2934*(u2886*(u546+u2885*(u2885*(u573+u641)+u410)+(u2885*(u485&
 &+u848)+u4118)*u2886)+u2885*(u572+u3867+u2885*((u821+u463)*u2885+u797+&
 &u543))+u3901+u833))
  value_4_=value_4_+(c(73)*u4032)
  value_5_=value_5_+(c(73)*(u2937*(u2886*(u4183+u2885*(u2885*(u676+u400&
 &1)+u881)+(u2885*(u467+u3805)+u4184)*u2886)+u2885*(u673+u2980*(u4158+u&
 &674)+u2885*(u2885*(u468+u505)+u2980*(u3889+u3881)+u675))+u2980*(u672+&
 &u4159)+u906)))
  value_6_=value_6_+(c(73)*u3967)
  u524=+u885*(-u3861*(u4007+1.2d1)+u495)
  u789=1.31764832376000117d-1*u4026*(u3964-u3861)
  u4184=-u4046
  u623=-8.47059636702857895d-2*u4195
  u794=8.47059636702857895d-1*u937
  u3964=7.52941899291429239d-2*u937
  u4183=3.7647094964571462d-2*u4134
  u3967=-1.69411927340571579d-1*u975
  u4195=3.01176759716571696d-1*u388
  value_7_=value_7_+(c(73)*(u2937*((u969+u2885*(u2885*(u794+u3863)+u418&
 &4))*u2886+u2885*(u789+u2980*(u3796+u3967)+u2885*(u2885*(u3964+u883)+u&
 &2980*(u4195+u4000)+u623))+u2980*(u4183+u4160)+u524)))
  u524=-5.91031598453575652d-3*u3778
  u789=2.60053903319573287d-1*u3993
  u4184=-1.7607816370596108d-1*u3965
  u623=2.11293796447153296d-1*u946
  u794=6.50134758298933218d-2*u4057
  u3993=-2.11293796447153296d-1*u4026
  u4183=+u3993*(-u4003+u4123)
  u3967=-u841
  u4195=-u3963
  u387=3.5215632741192216d-2*u899
  u3965=-u3966
  u3966=+u69*(u471+u890)
  u471=(u839+u893)
  u3963=u528*u2980
  u3964=u481*u2980
  value_8_=value_8_+(c(73)*(u2885*(u789+u2980*(u2980*(u3965+u3963)+u418&
 &3)+u2885*(u2885*(u623+u2980*(u3943+u4195)+u471*u2885)+u2980*(u3967+u2&
 &980*(u3882+u3966))+u4184))+u2980*(u794+u2980*(u3964+u387))+u524))
  value_9_=value_9_+(c(73)*(x*(u2885*(u953+u2980*(u4161+u840)+u2885*(u2&
 &885*(u855+u4051)+u2980*(u678+u4037)+u3888))+u2980*(u944+u4162)+u442)*&
 &z))
  u524=7.66191160463754555d-3*u873*(u4086+u3968)
  u789=-4.98024254301440461d-2*u4156
  u4184=1.99209701720576184d-1*u981
  u623=+u4187*(u884+u829)
  u794=-9.96048508602880921d-2*u3933
  u3968=4.48221828871296415d-1*u4207
  u3967=-1.99209701720576184d-1*u4022
  u4195=-2.98814552580864276d-1*u980
  u4187=4.98024254301440461d-2*u390
  u3965=u529*u2980
  u3966=u407*u2980
  value_1_=value_1_+(c(74)*(u2885*(u524+u2980*(u2980*(u4195+u3965)+u794&
 &)+u2885*(u2885*(u4184+u3967*u2980+(u4036+u3928)*u2885)+u2980*(u3968+u&
 &2980*(u3824+u4187))+u789))+u2980*(u623+u2980*(u3966+u4044))+u8))
  u524=-8.12668447873666522d-3*u4132
  u3967=-u989
  u4195=-3.16940694670729944d-1*u943
  u623=4.22587592894306592d-1*u937
  u794=u2980*(u3837+u3953)+u524
  u390=u2980*(u3952+u384)
  value_2_=value_2_+(c(74)*(u2937*((u947+u2885*(u2885*(u559+u3920)+u874&
 &))*u2886+u2885*(u3967+u390+u2885*((u623+u3827)*u2885+u459+u4195))+u79&
 &4)))
  u3967=u446*u2980
  value_3_=value_3_+(c(74)*(u2886*(u514+u2885*(u2885*(u3772+u563*u2885)&
 &+u3865)+u2886*((u66+u2885*(u647+u563))*u2886+u2885*(u3772+u2885*(u647&
 &+u679))+u843))+u2885*(u515+u2980*(u2980*(u715+u3956)+u714)+u2885*(u28&
 &85*(u717+u2980*(u3944+u415)+(u3777+u4035)*u2885)+u2980*(u716+u2980*(u&
 &3777+u72))+u969))+u2980*(u917+u2980*(u3967+u713))+u3813))
  u4195=u2980*(u809+u3954)+u816
  u3954=(u2937*(u2886*(u586+u2885*(u2885*(u588+u641)+u3976)+(u2885*(u38&
 &47+u848)+u395)*u2886)+u2885*(u587+u389+u2885*(u2980*(u862+u4147)+u556&
 &))+u4195))
  value_4_=value_4_+(c(74)*u3954)
  u789=5.38761366375499685d-1*u4026*(u449+5.49460693780709969d1*exppd2)
  u4184=-1.02524319000238052d7*u4064
  u4187=2.52544390488515477d-1*u937
  u3772=-3.78816585732773216d-2*u988
  u3968=1.70467463579747947d-1*u937
  u624=3.15680488110644347d-2*u4026*(u4102-u3861)
  u989=1.51526634293109286d-1*u941
  u387=-7.57633171465546432d-2*u981
  u981=u3798+u3805
  value_5_=value_5_+(c(74)*(u2934*(u2886*(u718+u2980*(u4163+u989)+u2885&
 &*(u2885*(u3968+u640)+u4184)+u2886*((u391+u981)*u2886+u2885*(u4187+u40&
 &01)+u2980*(u79+u3798)+u719))+u2885*(u789+u2980*(u2980*(u4041+u3776)+u&
 &825)+u2885*((u720+u3818)*u2885+u2980*(u3810+u3977)+u3772))+u2980*(u62&
 &4+u2980*(u3800+u387))+u516)))
  value_6_=value_6_+(c(74)*u781)
  u789=u873*(-u3861*(u4093+2.4d1)+u4086)
  u4184=-7.23982595472528115d-4*u789
  u4187=1.0352951115257152d-1*u4026
  u3772=u4187*(u3806-u3861)
  u3968=-1.69411927340571579d-1*u943
  u624=1.69411927340571579d-1*u937
  u989=-1.12941284893714386d-1*u4026
  u387=+u989*(u449+u3988)
  u3807=-1.52834231087603708d7*u4064
  u3935=-9.41177374114286549d-3*u4156
  u4207=5.6470642446857193d-2*u659
  u3810=4.70588687057143275d-2*u938
  u781=1.12941284893714386d-1*u939
  u3806=-9.41177374114286549d-3*u394
  u394=-2.82353212234285965d-2*u937
  u4156=(u624+u854)*u2886
  u4086=u2885*(u620+u3863)
  u4182=u394*u2980
  u4093=u2980*(u4182+u781)
  u4183=u742*u2885
  value_7_=value_7_+(c(74)*(u2934*(u2886*(u3772+u2980*(u4164+u439)+u288&
 &5*(u2885*(u598+u998)+u3807)+u2886*(u4156+u4086+u2980*(u4207+u854)+u39&
 &68))+u2885*(u387+u2980*(u2980*(u3806+u452)+u3810)+u2885*(u4183+u3834+&
 &u4028))+u2980*(u3935+u4093)+u4184)))
  u4184=1.62533689574733305d-2*u873*(u936+u4197)
  u3772=+u3993*(-u4003+u4091)
  u3968=3.16940694670729944d-1*u4208
  u387=1.05646898223576648d-1*u3786
  u4091=-1.05646898223576648d-1*u941
  u3807=-4.22587592894306592d-1*u939
  u4164=3.5215632741192216d-2*u659
  u3993=u2980*(u4164+u3882)
  u3806=u88*u2885
  u3810=u2885*(u53+u3806)
  value_8_=value_8_+(c(74)*(u2936*(u2886*(u611+u2885*(u3810+u612)+(u288&
 &5*(u4202+u3806)+u853)*u2886)+u2885*(u3772+u2980*(u3993+u4091)+u2885*(&
 &u4104*u2885+u3972+u3968))+u2980*(u387+u2980*(u3850+u3807))+u4184)))
  u4184=u558*u2980
  value_9_=value_9_+(c(74)*(u2934*((u622+u2885*(u2885*(u3940+u3909)+u38&
 &41))*u2886+u2885*(u724+u2980*(u2980*(u73+u455)+u725)+u2885*((u3822+u4&
 &56)*u2885+u2980*(u74+u3770)+u726))+u2980*(u722+u2980*(u4184+u723))+u9&
 &19)))
  u3772=2.29857348139126367d-2*u78
  u3968=-1.99209701720576184d-1*u4026*(-u3988+u4087)
  u387=1.49407276290432138d-1*u3786
  u722=-1.49407276290432138d-1*u941
  u3972=-5.97629105161728553d-1*u939
  u3822=4.98024254301440461d-2*u659
  u723=u2980*(u3822+u3824)
  u3807=u89*u2885
  u4197=u2885*(u54+u3807)
  value_1_=value_1_+(c(75)*(u2936*(u2886*(u613+u2885*(u4197+u614)+(u288&
 &5*(u4196+u3807)+u887)*u2886)+u2885*(u3968+u2980*(u723+u722)+u2885*(u4&
 &201*u2885+u3830+u765))+u2980*(u387+u2980*(u3955+u3972))+u3772)))
  u3772=1.05646898223576648d-1*u4152
  u3968=-3.16940694670729944d-1*u987
  u387=7.04312654823844319d-2*u457
  u3972=u62*u2886
  value_2_=value_2_+(c(75)*(u2934*(u2886*(u3772+u390+u2885*(u528*u2885+&
 &u4033)+u2886*(u3972+u2885*(u545+u3920)+u459+u3968))+u2885*(u387+u481*&
 &u2885)+u794)))
  u3772=u91*u2885
  u3968=u3884*u2885
  value_3_=value_3_+(c(75)*(u2936*(u2886*(u782+u2885*(u2885*(u84+u3772)&
 &+u925)+u2886*((u537+u647)*u2886+u2885*(u4185+u57*u2885)+u75))+u2885*(&
 &u727+u796+u2885*(u3968+u3785+u954))+u2980*(u748+u2980*(u3980+u749))+u&
 &517)))
  u517=2.66205161414224079d-2*u457
  u782=2.39584645272801671d-1*u937
  u796=u463+u848
  u3785=(u2934*(u2886*(u595+u389+u782*u800+u2886*((u709+u796)*u2886+u28&
 &85*(u811+u641)+u798+u596))+u2885*(u517+u809*u2885)+u4195))
  value_4_=value_4_+(c(75)*u3785)
  u4185=u3843*u2980
  value_5_=value_5_+(c(75)*(u2937*(u2886*(u754+u2980*(u2980*(u81+u3776)&
 &+u757)+u2885*(u2885*(u760+u640)+u3969)+u2886*((u3839+u981)*u2886+u288&
 &5*(u4148+u4001)+u2980*(u80+u4185)+u755))+u2885*(u4192+u2980*(u759*u29&
 &80+u55)+u2885*((u534+u3818)*u2885+u436+u721))+u2980*(u756+u2980*(u380&
 &0+u758))+u927)))
  value_6_=value_6_+(c(75)*u793)
  u387=-1.41176606117142983d-1*u3994
  u787=5.6470642446857193d-1*u4026*(u3979-u3988)
  u4192=2.54117891010857368d-1*u457
  u3853=-u4094
  u4094=-u3856
  u4194=4.23529818351428947d-1*u937
  u784=u4187*(u60-u3861)
  u386=-2.35294343528571638d-1*u938
  u783=1.50588379858285848d-1*u385
  u4190=+u989*(-u3873+u449)
  u3969=9.41177374114286549d-3*u4186
  u84=-4.70588687057143275d-2*u3921
  u4186=u84*u2980
  u4187=u394*u2885
  value_7_=value_7_+(c(75)*(u2937*(u2886*(u787+u2980*(u2980*(u3969+u452&
 &)+u386)+u2885*(u2885*(u4194+u998)+u3853)+u2886*(u4156+u2885*(u537+u38&
 &63)+u2980*(u783+u4186)+u3974))+u2885*(u4192+u2885*(u4187+u4094))+u298&
 &0*(u784+u2980*(u4182+u4190))+u387)))
  u387=-2.95515799226787826d-3*u3778
  u787=1.62533689574733305d-2*u78
  u4192=-2.81725061929537728d-1*u895
  u4194=-u785
  u784=2.43800534362099957d-2*u78
  u386=-u795
  u783=3.16940694670729944d-1*u940
  u4190=-4.22587592894306592d-1*u895
  u3969=-2.11293796447153296d-1*u489
  u989=-8.12668447873666522d-3*u3844
  u385=-u790
  u384=-1.05646898223576648d-1*u940
  u785=1.40862530964768864d-1*u3991
  u793=7.04312654823844319d-2*u489
  value_8_=value_8_+(c(75)*(u2886*(u787+u2980*(u2980*(u384+u3963)+u385)&
 &+u2885*(u2885*(u783+u4189*u2885)+u386)+u2886*((u4194+u2980*(u3882+u52&
 &8)+u2885*(u3806+u4189))*u2886+u2885*(u783+u2885*(u3806+u3969))+u2980*&
 &(u384+u2980*(u3882+u793))+u4192))+u2885*(u784+u2885*(u3837*u2885+u419&
 &0))+u2980*(u989+u2980*(u3964+u785))+u387))
  u384=u2980*(u4184+u4151)
  value_9_=value_9_+(c(75)*(u2937*(u2886*(u4188+u2980*(u2980*(u4193+u45&
 &5)+u564)+u2885*(u2885*(u818+u3909)+u4083)+(u2980*(u761+u455)+u3947)*u&
 &2886)+u2885*(u762+u2980*(u4130*u2980+u763)+u2885*((u3940+u456)*u2885+&
 &u3946+u764))+u2980*(u786+u384)+u4191)))
  u787=1.14928674069563183d-2*u78
  u385=-1.14172621348397963d7*u4064
  u786=+u385
  u784=1.49407276290432138d-1*u940
  u386=-1.99209701720576184d-1*u895
  u4190=-9.96048508602880921d-2*u489
  u3969=-1.14928674069563183d-2*u3844
  u989=-u385
  u385=-1.49407276290432138d-1*u940
  u4188=1.99209701720576184d-1*u3991
  u4193=9.96048508602880921d-2*u489
  value_1_=value_1_+(c(76)*(u2886*(u2980*(u2980*(u385+u3965)+u989)+u288&
 &5*(u2885*(u784+u3793*u2885)+u786)+u2886*((u4204+u2885*(u3807+u3793))*&
 &u2886+u2885*(u784+u2885*(u3807+u4190))+u2980*(u385+u2980*(u3824+u4193&
 &))))+u2885*(u787+u2885*(u4044*u2885+u386))+u2980*(u3969+u2980*(u3966+&
 &u4188))))
  u787=2.35619213399785474d8*u870
  u385=-8.12668447873666522d-3*u873
  u386=-u3861*(u817-1.5d1)
  u4193=+u385*(u386+u787)
  u4190=1.7607816370596108d-1*u932
  u3969=-3.16940694670729944d-1*u983
  u817=2.11293796447153296d-1*u457
  value_2_=value_2_+(c(76)*(u2937*(u2886*(u4190+u451+u2885*(u545*u2885+&
 &u4033)+u2886*(u3972+u2885*(u528+u3920)+u85+u3969))+u2885*(u817+u486*u&
 &2885)+u871+u4193)))
  u4193=6.31839356048751809d-3*u506
  u4190=-2.78009316661450796d-1*u934
  u3969=1.8823547482285731d-1*u461
  u989=-2.25882569787428772d-1*u946
  u4188=4.34389557283516869d-3*u78
  u786=-4.31531946600292822d6*u4064
  u785=5.6470642446857193d-2*u940
  u4191=-7.52941899291429239d-2*u895
  u932=-3.7647094964571462d-2*u489
  u871=u873*(-u4003*(u4099+1.0d0)+u498)
  u793=-1.30316867185055061d-2*u871
  u4099=u4026*(u4103-u3975)
  u790=9.03530279149715087d-1*u4099
  u4189=-3.95294497128000351d-1*u976
  u4192=1.12941284893714386d-1*u3832
  u451=u4026*(u417-u4003)
  u3970=1.8823547482285731d-1*u451
  u417=u4026*(-u4008+u3846)
  u4194=-5.6470642446857193d-2*u417
  u506=u4026*(u890+u3894)
  u387=7.52941899291429239d-2*u506
  value_3_=value_3_+(c(76)*(u2886*(u4190+u2980*(u2980*(u4194+u3956)+u79&
 &0)+u2885*(u2885*(u785+u4034*u2885)+u786)+u2886*(u2886*(u989+u2980*(u4&
 &165+u4192)+u2885*(u3772+u4034)+(u3998+u569)*u2886)+u2885*(u785+u2885*&
 &(u3772+u932))+u2980*(u4189+u2980*(u3777+u387))+u3969))+u2885*(u4188+u&
 &2885*(u446*u2885+u4191))+u2980*(u793+u2980*(u3967+u3970))+u4193))
  u3956=(u2937*(u2886*(u600+u3867+u2885*(u540*u2885+u51)+u2886*((u4092+&
 &u796)*u2886+u2885*(u4154+u641)+u797+u601))+u2885*(u602+u3845*u2885)+u&
 &3901+u838))
  value_4_=value_4_+(c(76)*u3956)
  u4190=-3.23774859600660869d-4*u873*(-u3861*(u4029+2.22d2)+u3942)
  u3942=3.24098634460261529d-1*u4026
  u989=u3942*(1.08340178853534212d7*u830-u3861)
  u387=-3.78816585732773216d-2*u4026
  u786=+u387*(u45+u3872)
  u4191=-1.3048126841906633d-1*u937
  u932=3.78816585732773216d-2*u457
  u4193=-1.20616845882633002d7*u4064
  u790=6.94497073843417563d-1*u937
  u4189=-u3809
  u4192=-3.78816585732773216d-2*u937
  u3970=7.57633171465546432d-2*u3983
  u4194=+u387*(u801+u4008)
  u387=2.02035512390812382d-1*u4026*(u801-1.13437046457953026d2*u470)
  u45=-6.31360976221288694d-2*u3921
  u3969=u98*u2886
  u4188=u45*u2980
  u793=u4188+u4001+u3969
  value_5_=value_5_+(c(76)*(u2934*(u2886*(u989+u2980*(u3799+u4194)+u288&
 &5*(u4192*u2885+u4193)+u2886*(u2886*(u4191+u793)+u2885*(u790+u640)+u29&
 &80*(u387+u3776)+u786))+u2885*(u932+u4189*u2885)+u2980*(u3970+u3858)+u&
 &4190)))
  value_6_=value_6_+(c(76)*u525)
  u3970=3.7647094964571462d-2*u457
  u989=7.5294189929142924d-1*u937
  u470=-u625
  u4190=-3.7647094964571462d-2*u3994
  u4191=-u856
  u932=-7.5294189929142924d-1*u937
  u856=u858+u3863
  u4189=u470*u2885
  value_7_=value_7_+(c(76)*(u2934*(u2886*(u2980*(u3852+u4191)+u2885*(u3&
 &968+u4067)+u2886*((u856)*u2886+u2885*(u989+u998)+u2980*(u932+u3887)))&
 &+u2885*(u3970+u4189)+u2980*(u4190+u4169))))
  u4191=u3912*(u787+u386)
  u3970=-1.7607816370596108d-1*u4026*(-u3861+u487)
  u989=3.16940694670729944d-1*u983
  u4190=-2.11293796447153296d-1*u3994
  value_8_=value_8_+(c(76)*(u2936*(u2886*(u3970+u2980*(u3939+u3808)+u28&
 &85*(u4202*u2885+u612)+u2886*((u4104+u3806)*u2886+u3810+u4203+u989))+u&
 &2885*(u611+u853*u2885)+u2980*(u4190+u4179)+u4191)))
  u3808=5.74643370347815916d-3*u792
  u3970=-7.47036381452160691d-2*u4152
  u989=2.24110914435648207d-1*u987
  u4190=-2.4901212715072023d-2*u3994
  u4191=-4.98024254301440461d-2*u3933
  u932=2.24110914435648207d-1*u975
  u786=-3.98419403441152369d-1*u388
  u3810=(u499+u456)*u2886
  value_9_=value_9_+(c(76)*(u2934*(u2886*(u3970+u2980*(u4170+u932)+u288&
 &5*(u3877*u2885+u887)+u2886*(u3810+u2885*(u855+u3909)+u2980*(u786+u380&
 &1)+u989))+u2885*(u4190+u879*u2885)+u2980*(u4191+u4171)+u3808)))
  u3808=u4026*(-u3861+u695)
  u4191=-1.49407276290432138d-1*u3808
  u989=u4026*(1.71938344913356967d8*u830+u3872)
  u3970=1.49407276290432138d-1*u989
  u4190=-2.98814552580864276d-1*u3994
  u4170=1.90287702247329938d7*u4064
  value_1_=value_1_+(c(77)*(u2936*(u2886*(u4191+u2980*(u3949+u4170)+u28&
 &85*(u4196*u2885+u614)+u2886*((u4201+u3807)*u2886+u4197+u4211+u3970))+&
 &u2885*(u613+u887*u2885)+u2980*(u4190+u3950)+u973)))
  u3793=2.01830586953240458d7*u4064
  u3970=u874*u2885
  u4190=u559*u2885
  u4191=u947*u2885
  value_2_=value_2_+(c(77)*(u2934*(u2886*(u766+u4172+u3970+u2886*((u383&
 &3+u4080)*u2886+u4190+u4173+u3793))+u4191+u4174+u851)))
  u932=2.17194778641758435d-2*u852
  u786=+u40*(u402+4.52134166253642395d8*u830)
  u787=u4133*(u788+u3872)
  u788=-6.77647709362286316d-1*u3994
  u402=6.83258915450463635d7*u4064
  u386=-2.25882569787428772d0*u937
  value_3_=value_3_+(c(77)*(u2936*(u2886*(u786+u2980*(u4180+u402)+u2885&
 &*(u460*u2885+u585)+u2886*((u4137+u4011+u3772)*u2886+u2885*(u47+u3772)&
 &+u2980*(u386+u3879)+u787))+u2885*(u584+u849*u2885)+u2980*(u788+u4181)&
 &+u932)))
  u852=1.02386600543932338d-2*u4097
  u625=-7.71994968101249829d-1*u3994
  u4011=4.83137012421188765d7*u4064
  u4180=-1.25116425864685317d0*u937
  u40=u4101+u92*u2886
  u4192=u4014*u2885
  u4193=u535*u2885
  u4194=u942*u2885
  u4181=(u2934*(u2886*(u625+u4176+u4192+u2886*(u2886*(u4180+u40)+u4193+&
 &u4177+u4011))+u4194+u4178+u852))
  value_4_=value_4_+(c(77)*u4181)
  u525=-4.85662289400991303d-3*u789
  u4176=6.31360976221288694d-2*u3786
  u789=u3910*(u435-u3872)
  u435=3.15680488110644347d-1*u457
  u4177=-3.13603799294845806d7*u4064
  u790=1.02280478147848768d0*u937
  u3910=2.52544390488515477d0*u4146
  u387=-6.31360976221288694d-2*u4026*(u933+u4008)
  u4146=u3905*(u791+u3872*(5.5d1+pd2))
  u3905=-u477
  u791=6.31360976221288694d-2*u937
  value_5_=value_5_+(c(77)*(u2937*(u2886*(u4176+u2980*(u791*u2980+u387)&
 &+u2885*(u4119*u2885+u4177)+u2886*(u2886*(u4122+u793)+u2885*(u790+u640&
 &)+u2980*(u4146+u3776)+u789))+u2885*(u435+u477*u2885)+u2980*(u3910+u39&
 &05*u2980)+u525)))
  u4178=-9.30787277672112164d-5*u3779*(-u4003*(9.8d1*u4027+1.5d1*u3929)&
 &+3.12036255583499681d8*u3898)
  u933=4.09546402175729352d-3*u873*(u3849-u3975*(1.185d3+3.92d2*pd2))
  u792=+u447*(1.20243269245430207d4*exppd2+u4124)
  u793=1.59723096848534447d-1*u4026*(u808+u3873)
  u388=-9.31718064949784276d-2*u937
  u934=1.02386600543932338d-3*u383
  u796=-5.32410322828448158d-2*u4026
  u3809=+u796*(6.91965983393513457d3*exppd2+u4124)
  u795=1.19792322636400836d0*u974
  u794=+u796*(u3971+u3862*(9.8d1+u4024))
  value_6_=value_6_+(c(77)*(u2886*(u933+u3809*u2980+u2886*(u2886*(u793+&
 &u794*u2980+(u3851+u388)*u2886)+u795*u2980+u792))+u934*u2980+u4178))
  u934=5.21267468740220243d-2*u4097
  u3809=-7.15294804326857777d-1*u3994
  u383=1.43843982200097607d7*u4064
  u794=3.38823854681143158d-1*u457
  u4024=-3.41629457725231818d7*u4064
  u795=+u4024
  u4124=-3.01176759716571696d-1*u3994
  u3862=3.17655460691882216d7*u4064
  u974=-1.0917657539725724d0*u937
  u796=+u974
  u3851=-u866
  u3971=u3851*u2980
  u4195=u483*u2885
  u4196=u866*u2885
  value_7_=value_7_+(c(77)*(u2937*(u2886*(u3809+u2980*(u3854+u3862)+u28&
 &85*(u4195+u795)+u2886*((u4034+u856)*u2886+u2885*(u4073+u998)+u2980*(u&
 &796+u3887)+u383))+u2885*(u794+u4196)+u2980*(u4124+u3971)+u934)))
  u3809=-8.12668447873666522d-3*u3886
  u794=3.25067379149466609d-1*u4057
  u795=-6.33881389341459887d-1*u478
  u4124=3.38070074315445273d0*u962
  u4073=1.18264470046224369d7*u830
  u796=-2.95811315026014614d0*u4026*(-u3988+u4073)
  u797=3.16940694670729944d0*u979
  u798=-4.22587592894306592d-1*u3814
  u85=-1.40862530964768864d-1*u3994
  u389=-2.11293796447153296d0*u937
  value_8_=value_8_+(c(77)*(u2886*(u794+u2980*(u4153*u2980+u796)+u2886*&
 &(u2886*(u4124+u2980*(u3978+u798)+u471*u2886)+u2980*(u797+u389*u2980)+&
 &u795))+u2980*(u708+u85*u2980)+u3809))
  u3978=2.24110914435648207d-1*u457
  u794=-2.24110914435648207d-1*u3808
  u795=2.24110914435648207d-1*u989
  u796=-7.47036381452160691d-2*u3994
  u797=-2.19130671892633803d0*u847
  u798=3.73518190726080346d-1*u972
  u85=-4.98024254301440461d-2*u425
  u389=-8.71542445027520807d-1*u937
  value_9_=value_9_+(c(77)*(u2937*(u2886*(u794+u2980*(u389*u2980+u798)+&
 &u2885*(u855*u2885+u887)+u2886*(u3810+u2885*(u3877+u3909)+u2980*(u85+u&
 &3801)+u795))+u2885*(u796+u3838*u2885)+u2980*(u797+u3841*u2980)+u3978)&
 &))
  u3809=-1.60900143697388457d-1*u3886
  u4200=3.1375528020990749d0*u457
  u3810=-u431
  u4197=-u430
  u3808=u3871*u2885
  u430=u4036+u3808
  value_1_=value_1_+(c(78)*(u2934*(u2885*(u4200+u4199*u2980+u2885*(u288&
 &5*(u4197+u430)+u4198*u2980+u3810))+u583*u2980+u3809)))
  u3809=-1.13773582702313313d-1*u3886
  u4200=2.21858486269510961d0*u457
  u3810=-u4082
  u4197=u4117*u2980
  u4198=u822*u2980
  u4199=u570*u2980
  value_2_=value_2_+(c(78)*(y*(u2885*(u4200+u4197+u2885*(u2885*(u3810+u&
 &4080)+u4198+u3811))+u4199+u3809)*z))
  u3809=u3895*u2885
  u4082=u3879+u3809
  u3810=u96*u2885
  value_3_=value_3_+(c(78)*(u2934*((u627+u2885*(u2885*(u493+u3810)+u411&
 &))*u2886+u2885*(u992+u4155+u2885*(u2885*(u628+u4082)+u3796+u4078))+u4&
 &175+u500)))
  u4200=u3885*u2980
  u4201=u552*u2980
  u4202=u948*u2980
  value_4_=value_4_+(c(78)*(u2936*((u636+u2885*(u2885*(u4072+u848)+u419&
 &))*u2886+u2885*(u994+u4200+u2885*(u2885*(u637+u4101)+u4201+u4056))+u4&
 &202+u502)))
  value_5_=value_5_+(c(78)*(u2886*(u501+u2885*(u2885*(u3918+u2885*(u400&
 &1+u635))+u993)+(u2885*(u4112+u2885*(u3805+u433))+u629)*u2886)+u2885*(&
 &u897+u2980*(u4129*u2980+u631)+u2885*(u2885*(u634+u2980*(u3881+u59)+u2&
 &885*(u640+u4002+u434))+u2980*(u633+u3915*u2980)+u632))+u2980*(u896+u6&
 &30*u2980)+u3892))
  value_6_=value_6_+(c(78)*u504)
  u437=3.94899597530469881d-4*u3779*(u401-u4003*(u4006+u4019))
  u4006=-3.04072690098461808d-2*u3886
  u3972=-8.68779114567033738d-3*u873*(u3780+u498)
  u464=1.5811779885120014d0*u457
  u3780=+u3784*(u806+u4087)
  u4203=-u411
  u4204=u876*(u832-u4008)
  u4208=-u827
  u4209=-3.10588533457714561d-1*u937
  u390=-8.68779114567033738d-3*u4205
  u4205=5.6470642446857193d-2*u4206
  u4206=-8.47059636702857895d-1*u952
  u4207=5.6470642446857193d-2*u413
  u4211=6.58824161880000585d-2*u457
  u827=9.88236242820000877d-1*u937
  value_7_=value_7_+(c(78)*((u4006+u2885*(u2885*(u4203+u2885*(u3863+u42&
 &08))+u464))*u2886+u2885*(u3972+u2980*(u3853*u2980+u4205)+u2885*(u2885&
 &*(u4204+u2980*(u4000+u4207)+u2885*(u998+u854+u4209))+u2980*(u4206+u82&
 &7*u2980)+u3780))+u2980*(u390+u4211*u2980)+u437))
  u4205=2.6411724555894162d0*u457
  u4206=-u475
  u4207=-u875
  u475=u4120*u2885
  u875=u3875+u475
  u4203=u4012*u2980
  u4204=u553*u2980
  value_8_=value_8_+(c(78)*(x*(u2885*(u4205+u3855+u2885*(u2885*(u4207+u&
 &875)+u4203+u4206))+u4204+u522)*z))
  value_9_=value_9_+(c(78)*(u2885*(u898+u2980*(u4083*u2980+u995)+u2885*&
 &(u2885*(u868+u2980*(u4037+u639)+u2885*(u3909+u4025+u882))+u2980*(u387&
 &4+u818*u2980)+u638))+u2980*(u503+u604*u2980)+u4115))
  u4205=9.46246083172736875d-1*u457
  u4206=-u3812
  u4207=-u815
  value_1_=value_1_+(c(79)*(x*(u2885*(u4205+u3958+u2885*(u2885*(u4207+u&
 &430)+u3959+u4206))+u3960+u523)*z))
  u4205=2.46263166022323189d-4*u3779*(u4040-u4003*(1.5d1*u4095+u496))
  u4206=-1.89622637837188855d-2*u3886
  u4207=-2.70889482624555508d-3*u873*(-u3861*(u4075+1.65d2)+u4125)
  u4075=9.86037716753382047d-1*u457
  u4040=u4026*(u494+u3819)
  u4125=1.7607816370596108d-1*u4040
  u4211=-u4071
  u3972=-1.05646898223576648d-1*u3828
  u4208=9.86037716753382047d-1*u937
  u4209=1.40862530964768864d-1*u937
  u4071=u3916*(u3985-u4003)
  u390=-u4070
  u4070=2.11293796447153296d-1*u4043
  u4043=u2980*(u454*u2980+u4071)
  u939=u2980*(u3827+u4070)
  u494=u2980*(u390+u3917*u2980)
  u3828=u2980*(u554+u949*u2980)
  value_2_=value_2_+(c(79)*((u4206+u2885*(u2885*(u4211+u2885*(u3920+u42&
 &08))+u4075))*u2886+u2885*(u4207+u4043+u2885*(u2885*(u3972+u939+(u3827&
 &+u4209)*u2885)+u494+u4125))+u3828+u4205))
  u4205=u399*u2980
  u4206=u571*u2980
  u4207=u963*u2980
  value_3_=value_3_+(c(79)*(u2937*((u648+u2885*(u2885*(u3788+u3810)+u39&
 &6))*u2886+u2885*(u649+u4205+u2885*(u2885*(u650+u4082)+u4206+u849))+u4&
 &207+u900)))
  u3788=u2980*(u4111*u2980+u607)
  u396=u2980*(u608+u421*u2980)
  u648=u2980*(u942+u546*u2980)
  value_4_=value_4_+(c(79)*(u2886*(u902+u2885*(u2885*(u867+u2885*(u641+&
 &u658))+u420)+(u2885*(u4107+u2885*(u848+u4062))+u551)*u2886)+u2885*(u9&
 &03+u3788+u2885*(u2885*(u657+u2980*(u4147+u860))+u396+u656))+u648+u466&
 &))
  u3972=u4020*u2885
  u4208=u4119*u2980
  u4209=u477*u2980
  value_5_=value_5_+(c(79)*(u2936*(u2886*(u651+u2885*(u63*u2885+u654)+(&
 &u2885*(u3876+u3972)+u3931)*u2886)+u2885*(u653+u2980*(u4208+u392)+u288&
 &5*(u2885*(u480+u4131+u640)+u2980*(u416+u4015)+u655))+u2980*(u652+u420&
 &9)+u901)))
  value_6_=value_6_+(c(79)*u4032)
  u4032=+u3859*(-u3861*(u3907+6.0d1)+u3860)
  u3859=u4026*(u3996-u3975)
  u3975=4.06588625617371789d0*u3859
  u3860=4.70588687057143275d-1*u4026*(u3791-u3861)
  u3791=-5.6470642446857193d-1*u951
  u3996=2.25882569787428772d-1*u4121
  u4211=-u4076
  value_7_=value_7_+(c(79)*(u2936*(u2886*(u3975+u2885*(u2885*(u3996+u64&
 &7)+u3791)+(u4086+u399)*u2886)+u2885*(u3860+u2980*(u3854+u4211)+u2885*&
 &(u93*u800+u2980*(u536+u3887)+u3974))+u2980*(u969+u3971)+u4032)))
  u4211=+u385*(-u3861*(u400+4.5d1)+u393)
  u3974=-u4074*(u465-u4003)
  u3975=-1.05646898223576648d-1*u999
  u999=u2980*(u545*u2980+u4138)
  u4074=u2980*(u985+u486*u2980)
  value_8_=value_8_+(c(79)*(u2934*((u957+u2885*(u2885*(u3997+u475)+u381&
 &1))*u2886+u2885*(u3974+u999+u2885*(u412*u2885+u743+u3975))+u4074+u421&
 &1)))
  value_9_=value_9_+(c(79)*(y*(u2885*(u566+u2980*(u3961+u397)+u2885*(u2&
 &885*(u476+u4051)+u2980*(u660+u4037)+u526))+u2980*(u819+u3962)+u850)*z&
 &))
  u4211=u873*(-u3861*(u4010+3.0d1)+u429)
  u429=-3.83095580231877278d-3*u4211
  u3975=u4026*(u487+u3819)
  u397=9.96048508602880921d-2*u3975
  u3974=-1.49407276290432138d-1*u968
  u526=1.99209701720576184d-1*u937
  u487=u2980*(u547*u2980+u745)
  u4010=u2980*(u986+u453*u2980)
  value_1_=value_1_+(c(80)*(u2934*((u991+u2885*(u2885*(u3906+u3808)+u80&
 &7))*u2886+u2885*(u397+u487+u2885*((u526+u3824)*u2885+u744+u3974))+u40&
 &10+u429)))
  u3824=-8.12668447873666522d-3*u4211
  u3906=9.15606451270997615d-1*u4026*(u423-u4003)
  u3974=2.11293796447153296d-1*u3975
  u526=-u58
  u4211=-3.16940694670729944d-1*u968
  u3975=2.11293796447153296d-1*u3973
  u58=-1.40862530964768864d-1*u3921
  u3921=u2885*(u422+u475)
  u3973=u58*u2885
  value_2_=value_2_+(c(80)*(u2936*(u2886*(u3906+u2885*(u2885*(u3975+u39&
 &73)+u526)+(u3921+u874)*u2886)+u2885*(u3974+u2885*(u623*u2885+u4211))+&
 &u3824)))
  u968=u2980*(u460*u2980+u585)
  u3824=u2980*(u584+u849*u2980)
  value_3_=value_3_+(c(80)*(u2934*(u2886*(u680+u2885*(u2885*(u799+u3809&
 &)+u3926)+(u2885*(u3896+u3810)+u3925)*u2886)+u2885*(u681+u968+u2885*((&
 &u3781+u3777)*u2885+u804+u615))+u3824+u908)))
  u3974=u859*u2885
  value_4_=value_4_+(c(80)*(u2936*(u2886*(u685+u2885*(u67*u2885+u687)+(&
 &u2885*(u812+u3974)+u3976)*u2886)+u2885*(u686+u556*u2885)+u911)))
  value_5_=value_5_+(c(80)*(u2886*(u909+u2885*(u2885*(u684+u2885*(u640+&
 &u805))+u683)+u2886*((u4127+u2885*(u3972+u842))*u2886+u2885*(u616+u406&
 &9*u2885)+u682))+u2885*(u910+u2980*(u2980*(u4210+u3799)+u474)+u2885*(u&
 &2885*(u3899+u2980*(u3977+u857)+(u3818+u589)*u2885)+u2980*(u3831+u2980&
 &*(u3776+u990))+u3981))+u2980*(u510+u2980*(u3858+u891))+u5))
  value_6_=value_6_+(c(80)*u3954)
  u526=6.58165995884116468d-5*u3779*(u4066-u4003*(1.5d1*u813+u4065))
  u4211=-7.23982595472528115d-4*u873
  u3975=+u4211*(-u4003*(u810+3.81d2)+u469)
  u990=u4026*(u4103+u4047)
  u3976=5.6470642446857193d-2*u990
  u4103=-u4028
  u4028=+u4211*(-u4003*(u935+3.45d2)+u4060)
  u4211=u3984*(u3986+u4047)
  u935=-6.77647709362286316d-1*u946
  u946=6.77647709362286316d-1*u3859
  u3986=-3.38823854681143158d-1*u4136
  u4136=5.6470642446857193d-2*u3995
  u3995=-u4048
  u4048=4.14118044610286082d-1*u937
  u3859=2.17194778641758435d-3*u78
  u799=2.82353212234285965d-2*u940
  u940=-3.7647094964571462d-2*u895
  u895=-1.8823547482285731d-2*u489
  u489=u2885*(u3863+u537)
  u4210=u470*u2980
  value_7_=value_7_+(c(80)*(u2886*(u3975+u2885*(u2885*(u3986+u2885*(u99&
 &8+u4048))+u4211)+u2886*((u4103+u489)*u2886+u2885*(u935+u2885*(u647+u4&
 &136))+u3976))+u2885*(u4028+u2980*(u2980*(u799+u3980)+u66)+u2885*(u288&
 &5*(u3995+u2980*(u452+u3884)+u4183)+u2980*(u799+u2980*(u452+u895))+u94&
 &6))+u2980*(u3859+u2980*(u4210+u940))+u526))
  u526=+u385*(-u3861*(u4021+1.8d1)+u3817)
  u4211=u3916*(u4087+u3819)
  u935=-1.05646898223576648d-1*u488
  u488=u2980*(u4114+u3964)+u526
  u3964=u2980*(u3963+u4091)
  value_8_=value_8_+(c(80)*(u2937*((u949+u2885*(u2885*(u3917+u475)+u454&
 &))*u2886+u2885*(u4211+u3964+u2885*((u623+u3882)*u2885+u3993+u935))+u4&
 &88)))
  u4211=u404*u2980
  value_9_=value_9_+(c(80)*((u912+u2885*(u2885*(u4128+u2885*(u3909+u418&
 &))+u583))*u2886+u2885*(u914+u2980*(u2980*(u690+u3948)+u689)+u2885*(u2&
 &885*(u693+u2980*(u3770+u68)+(u456+u3890)*u2885)+u2980*(u692+u2980*(u4&
 &55+u4038))+u691))+u2980*(u913+u2980*(u4211+u688))+u6))
  u935=-3.4478602220868955d-2*u512
  u3975=1.19525821032345711d0*u3991
  u3976=-1.49407276290432138d-1*u945
  u945=u2980*(u4044+u3966)+u935
  u3966=u2980*(u3965+u722)
  value_1_=value_1_+(c(81)*(u2937*((u944+u2885*(u2885*(u560+u3808)+u411&
 &6))*u2886+u2885*(u3975+u3966+u2885*(u428*u2885+u723+u3976))+u945)))
  u428=7.38789498066969565d-4*u4108
  u3965=+u385*(-u4003*(u4049+2.1d1)+u936)
  u4049=3.5215632741192216d-2*u4040
  u4108=-1.46280320617259974d-1*u512
  u3975=6.33881389341459887d-1*u878
  u3976=-6.33881389341459887d-1*u943
  u936=1.05646898223576648d-1*u990
  u943=-3.16940694670729944d-1*u984
  u984=1.05646898223576648d-1*u413
  u990=8.45175185788613183d-1*u937
  value_2_=value_2_+(c(81)*(u2886*(u3965+u2885*(u2885*(u943+u990*u2885)&
 &+u3975)+u2886*((u3825+u2885*(u475+u533))*u2886+u2885*(u3976+u2885*(u3&
 &973+u984))+u4049))+u2885*(u4108+u2885*(u3999*u2885+u936))+u428))
  u3965=u2980*(u4052+u3967)+u837
  value_3_=value_3_+(c(81)*(u2937*(u2886*(u728+u2885*(u2885*(u458+u3809&
 &)+u4110)+(u2885*(u3893+u3810)+u4212)*u2886)+u2885*(u729+u747+u2885*(u&
 &2980*(u861+u3772)+u730))+u3965)))
  value_4_=value_4_+(c(81)*(u2886*(u921+u2885*(u527*u2885+u736)+u2886*(&
 &(u3829+u2885*(u3974+u864))*u2886+u2885*(u737+u50*u2885)+u735))+u2885*&
 &(u922+u597*u2885)+u7))
  value_5_=value_5_+(c(81)*(u2936*(u2886*(u731+u2885*(u2885*(u4013+u394&
 &5)+u734)+u2886*((u391+u3972)*u2886+u2885*(u76+u824*u2885)+u732))+u288&
 &5*(u733+u2980*(u4053+u4055)+u2885*(u4166+u2980*(u532+u3776)+u865))+u2&
 &980*(u4005+u2980*(u3800+u4009))+u920)))
  value_6_=value_6_+(c(81)*u3785)
  u3785=-6.51584335925275303d-3*u873*(-u3861*(u3775+4.4d1)+u3783)
  u3775=u886*(u4150+u3990)
  u3990=-1.69411927340571579d-1*u951
  u3975=u4133*(u3985+u803)
  u3976=-1.41176606117142983d-1*u4149
  u391=5.6470642446857193d-2*u971
  u4212=-u3835
  u3835=9.31765600373143684d-1*u937
  u990=-2.82353212234285965d-2*u3786
  u936=2.82353212234285965d-2*u941
  u800=-9.41177374114286549d-3*u659
  value_7_=value_7_+(c(81)*(u2936*(u2886*(u3775+u2885*(u2885*(u3835+u99&
 &8)+u3976)+u2886*((u624+u3863)*u2886+u2885*(u391+u647)+u3990))+u2885*(&
 &u3975+u2980*(u2980*(u800+u452)+u936)+u2885*(u4187+u2980*(u3884+u452)+&
 &u4212))+u2980*(u990+u4093)+u3785)))
  u3975=5.28234491117883239d-1*u4026*(u3820-u3861)
  u3820=-u4153
  u4153=(u530+u3882)*u2886
  value_8_=value_8_+(c(81)*(u2934*(u2886*(u3975+u3964+u2885*(u533*u2885&
 &+u3820)+u2886*(u4153+u3921+u3993+u440))+u2885*(u817+u3825*u2885)+u488&
 &)))
  u3975=u3982*u2885
  u3976=u97*u2885
  value_9_=value_9_+(c(81)*(u2936*(u2886*(u738+u2885*(u2885*(u77+u3975)&
 &+u740)+(u2885*(u4105+u3976)+u3836)*u2886)+u2885*(u739+u2980*(u4054+u3&
 &891)+u2885*(u3940*u2885+u2980*(u544+u455)+u741))+u2980*(u409+u384)+u9&
 &23)))
  u990=1.49407276290432138d-1*u878
  u878=-u4170
  u4170=u863*u2886
  value_1_=value_1_+(c(82)*(u2934*(u2886*(u990+u3966+u2885*(u529*u2885+&
 &u878)+u2886*(u4170+u2885*(u547+u3808)+u723+u746))+u2885*(u3782+u407*u&
 &2885)+u945)))
  u990=+u385*(-u3861*(u823+1.5d1)+u3792)
  u936=4.57803225635498807d-1*u4026*(u511-u3861)
  u800=-3.16940694670729944d-1*u872
  u872=u4030*(u977+u4050)
  u391=-5.28234491117883239d-1*u490
  u490=3.5215632741192216d-2*u4077
  u4212=-u606
  value_2_=value_2_+(c(82)*(u2936*(u2886*(u936+u2885*(u4212*u2885+u391)&
 &+u2886*((u549+u475)*u2886+u2885*(u490+u3973)+u800))+u2885*(u872+u4033&
 &*u2885)+u990)))
  u990=-7.52941899291429239d-2*u3994
  u936=2.51726968850170813d7*u4064
  u800=-1.50588379858285848d0*u937
  u3973=u3777+u3810
  value_3_=value_3_+(c(82)*(u2934*(u2886*(u767+u747+u2885*(u569*u2885+u&
 &936)+u2886*((u4089+u3973)*u2886+u2885*(u800+u3809)+u3866+u768))+u2885&
 &*(u990+u4052*u2885)+u3965)))
  value_4_=value_4_+(c(82)*(u2936*(u2886*(u771+u774*u2885+u2886*((u709+&
 &u3974)*u2886+u82*u2885+u772))+u773*u2885+u929)))
  u709=-1.61887429800330434d-4*u873*(-u4003*(3.02d2*pd2+8.61d2)+9.61581&
 &114145070446d8*u870)
  u3974=u3942*(1.52999489220639918d7*u830-u3861)
  u3942=+u828*(8.17237812242499165d8*u830+u4023)
  u990=-9.25996098457890084d-2*u937
  u4023=4.41952683354902085d-2*u4026*(u4073-u3861)
  u4073=+u802*(u3846+u3872)
  u3846=u462*(1.7862483610443196d10*u830+u3872*(1.65d2+u3992))
  u3992=u426*(5.66759729529213707d8*u870-u4003*(9.0d0+1.78d2*pd2))
  u936=+u414*(u4047+5.28551208437356603d8*u830)
  u800=u880*(4.80790557072535223d9*u830-3.44848621232177198d4*u3989)
  u3989=+u846*(1.43982443647814853d10*u830+u3872*(1.33d2+u4016))
  u4016=+u427*(-u3988+9.09726692663264377d5*u830)
  u391=u880*(u4018-u4008)
  u4008=+u846*(-u3872*(pd2-4.5d1)+u3848)
  u4212=u45*u2885
  value_5_=value_5_+(c(82)*(u2886*(u928+u2980*(u2980*(u391+u3799)+u936)&
 &+u2885*(u2885*(u4073+u4167)+u3974)+u2886*(u2886*(u770+u2980*(u4185+u3&
 &989)+u2885*(u4212+u990)+(u3805+u3798+u382)*u2886)+u2885*(u3942+u2885*&
 &(u3945+u3846))+u2980*(u800+u2980*(u3776+u4008))+u769))+u2885*(u709+u2&
 &885*(u4168+u4023))+u2980*(u3992+u2980*(u3858+u4016))+u9))
  value_6_=value_6_+(c(82)*u3956)
  u3845=1.9744979876523494d-4*u3779*(u4063-u4003*(3.0d0*u3897+u403))
  u3897=-2.17194778641758435d-3*u873*(-u4003*(u3902+1.83d2)+u492)
  u492=6.58824161880000585d-2*u4026*(u3826+u3819)
  u990=-2.82353212234285965d-2*u3904
  u3904=-1.95475300777582591d-2*u3886
  u3826=8.47059636702857895d-1*u457
  u3819=-2.96678213287701315d7*u4064
  u3902=2.82353212234285965d-2*u457
  u936=7.05883030585714913d-1*u937
  u800=-6.51584335925275303d-3*u871
  u871=4.51765139574857544d-1*u4099
  u4099=-1.97647248564000175d-1*u976
  u391=5.6470642446857193d-2*u3832
  u3832=9.4117737411428655d-2*u451
  u451=-2.82353212234285965d-2*u417
  u417=3.7647094964571462d-2*u506
  value_7_=value_7_+(c(82)*(u2886*(u3897+u2980*(u2980*(u451+u3980)+u871&
 &)+u2885*(u2885*(u4090+u3968)+u3826)+u2886*(u2886*(u990+u2980*(u4186+u&
 &391)+u489+(u854+u536)*u2886)+u2885*(u3819+u2885*(u998+u936))+u2980*(u&
 &4099+u2980*(u452+u417))+u492))+u2885*(u3904+u2885*(u4189+u3902))+u298&
 &0*(u800+u2980*(u4210+u3832))+u3845))
  u391=-1.21900267181049978d-1*u512
  u936=2.46509429188345512d-1*u4026*(u918-u3861)
  u800=-1.05646898223576648d-1*u988
  u990=6.33881389341459887d-1*u457
  value_8_=value_8_+(c(82)*(u2937*(u2886*(u936+u999+u2885*(u422*u2885+u&
 &3820)+u2886*(u4153+u2885*(u533+u475)+u743+u800))+u2885*(u990+u3970)+u&
 &4074+u391)))
  value_9_=value_9_+(c(82)*(u2886*(u930+u2980*(u2980*(u605+u3948)+u3957&
 &)+u2885*(u2885*(u779+u844*u2885)+u776)+u2886*((u3838+u444+u2885*(u397&
 &6+u4130))*u2886+u2885*(u777+u2885*(u3975+u83))+u2980*(u605+u2980*(u45&
 &5+u609))+u775))+u2885*(u931+u2885*(u3908*u2885+u778))+u2980*(u521+u29&
 &80*(u4211+u4157))+u520))
  u391=1.49407276290432138d-1*u4026*(u695-u3861)
  u936=-1.49407276290432138d-1*u989
  u800=2.98814552580864276d-1*u457
  value_1_=value_1_+(c(83)*(u2937*(u2886*(u391+u487+u2885*(u547*u2885+u&
 &878)+u2886*(u4170+u2885*(u529+u3808)+u744+u936))+u2885*(u800+u453*u28&
 &85)+u4010+u619)))
  u391=8.12668447873666522d-3*u4097
  u3808=-3.25067379149466609d-1*u3913
  u800=6.33881389341459887d-1*u4026*(u3979-u4003)
  u936=-u4124
  u4124=-u3793
  value_2_=value_2_+(c(83)*(u2886*(u3808+u4043+u2885*(u3970+u990)+u2886&
 &*(u2886*(u936+u939+u2885*(u3920+u528)+(u3827+u530)*u2886)+u2885*(u412&
 &4+u4190)+u494+u800))+u2885*(u443+u4191)+u3828+u391))
  value_3_=value_3_+(c(83)*(u2937*(u2886*(u786+u968+u2885*(u576*u2885+u&
 &402)+u2886*((u4137+u3973)*u2886+u2885*(u386+u3809)+u804+u787))+u2885*&
 &(u788+u831*u2885)+u3824+u932)))
  value_4_=value_4_+(c(83)*(u2886*(u933+u3788+u2885*(u4192+u625)+u2886*&
 &(u2886*(u793+u3773+u2885*(u641+u4180)+(u848+u463+u388)*u2886)+u2885*(&
 &u4011+u4193)+u396+u792))+u2885*(u852+u4194)+u648+u4178))
  value_5_=value_5_+(c(83)*(u2936*(u2886*(u4176+u2980*(u4208+u4177)+u28&
 &85*(u791*u2885+u387)+u2886*(u2886*(u4122+u3795+u4212+u3969)+u2885*(u4&
 &146+u3945)+u2980*(u790+u4015)+u789))+u2885*(u3910+u3905*u2885)+u2980*&
 &(u435+u4209)+u525)))
  value_6_=value_6_+(c(83)*u4181)
  u4011=-5.21267468740220243d-2*u3886
  u4014=7.15294804326857777d-1*u457
  u4180=-u383
  u936=3.01176759716571696d-1*u457
  u800=-u3862
  u3862=-u974
  u391=-3.38823854681143158d-1*u3994
  u974=-u4024
  value_7_=value_7_+(c(83)*(u2936*(u2886*(u4014+u2980*(u3854+u974)+u288&
 &5*(u4195+u800)+u2886*((u569+u856)*u2886+u2885*(u3862+u998)+u2980*(u38&
 &96+u3887)+u4180))+u2885*(u936+u4196)+u2980*(u391+u3971)+u4011)))
  value_8_=value_8_+(c(83)*(u2934*(u2886*(u990+u3855+u454*u2885+u2886*(&
 &(u528+u875)*u2886+u3917*u2885+u4203+u4124))+u949*u2885+u4204+u443)))
  value_9_=value_9_+(c(83)*(u2936*(u2886*(u794+u2980*(u3961+u887)+u2885&
 &*(u389*u2885+u798)+u2886*((u499+u3976)*u2886+u2885*(u85+u3975)+u4058+&
 &u795))+u2885*(u797+u3841*u2885)+u2980*(u796+u3962)+u3978)))
  value_1_=value_1_+(c(84)*(u2934*(u2886*(u3958+u4116*u2885+u2886*((u43&
 &0)*u2886+u560*u2885+u3959))+u944*u2885+u3960)))
  value_2_=value_2_+(c(84)*(u2936*(u2886*(u4197+u621*u2885+u2886*((u408&
 &0)*u2886+u780*u2885+u4198))+u4106*u2885+u4199)))
  u936=-2.82353212234285965d0*u3994
  u800=1.40247882645095167d8*u4064
  u391=-3.04941469213028842d0*u937
  value_3_=value_3_+(c(84)*(u2934*(u2886*(u936+u4205+u399*u2885+u2886*(&
 &u2886*(u391+u4082+u96*u2886)+u571*u2885+u4206+u800))+u963*u2885+u4207&
 &+u934)))
  u936=1.72009488913806328d-1*u4097
  u800=-3.35418503381922339d0*u3994
  u2934=1.06798708008894359d8*u4064
  u391=-1.75695406533387892d0*u937
  u3895=(u2886*(u800+u4200+u3885*u2885+u2886*(u2886*(u391+u40)+u552*u28&
 &85+u4201+u2934))+u948*u2885+u4202+u936)
  value_4_=value_4_+(c(84)*(u2936*u3895))
  u2936=-2.6490670330963162d-4*u3779*(-u4003*(u4019+5.0d0*u3903)+u401)
  u401=2.33117898912475825d-2*u873*(u498+u3987*(1.035d3+1.12d2*pd2))
  u4200=-3.28307707635070121d-1*u4026*(2.94936320790677867d3*exppd2+u38&
 &40)
  u3779=5.05088780977030955d-2*u4026*(9.36108766750499044d8*u830-u3873)
  u4201=-6.43988195745714467d-1*u937
  u4019=-2.42831144700495651d-2*u3886
  u3903=1.30060361101585471d0*u457
  u4003=-6.39269283177954912d7*u4064
  u3873=1.37636692816240935d0*u937
  u3840=-6.31360976221288694d-3*u3994
  u40=1.80925268823949503d6*u4064
  u498=-9.4704146433193304d-2*u937
  u3987=-1.94264915760396521d-3*u873*(-u3988*(u445+3.75d2)+u393)
  u393=1.6415385381753506d-1*u4026*(7.29880723421372881d7*u830+u484)
  u484=-1.89408292866386608d-1*u4026*(u4004+u3872)
  u4004=u3789*(1.54808191290507699d10*u830+u3872*(1.43d2+u3774))
  u3774=6.31360976221288694d-3*u457
  u3789=-u40
  u3872=9.4704146433193304d-2*u937
  value_5_=value_5_+(c(84)*(u2886*(u401+u2980*(u3789*u2980+u393)+u2885*&
 &(u40*u2885+u3903)+u2886*(u2886*(u3779+u2980*(u3776+u4004)+u2885*(u640&
 &+u3873)+u2886*(u3969+u4001+u4188+u4201))+u2885*(u4003+u498*u2885)+u29&
 &80*(u484+u3872*u2980)+u4200))+u2885*(u4019+u3840*u2885)+u2980*(u3987+&
 &u3774*u2980)+u2936))
  value_6_=value_6_+(c(84)*(u2937*u3895))
  u552=-2.60633734370110121d-2*u3886
  u800=1.41176606117142983d0*u457
  u391=-7.01239413225475836d7*u4064
  u92=+u391
  u948=1.52470734606514421d0*u937
  u3885=-9.41177374114286549d-3*u3994
  u2934=2.60633734370110121d-2*u4097
  u3932=-1.41176606117142983d0*u3994
  u3895=-u391
  u936=-u948
  u391=9.41177374114286549d-3*u457
  value_7_=value_7_+(c(84)*(u2886*(u2980*(u4094*u2980+u3932)+u2885*(u38&
 &56*u2885+u800)+u2886*(u2886*(u2980*(u3887+u936)+u2885*(u998+u948)+(u3&
 &863+u858)*u2886)+u2885*(u92+u432*u2885)+u2980*(u3895+u598*u2980)))+u2&
 &885*(u552+u3885*u2885)+u2980*(u2934+u391*u2980)))
  value_8_=value_8_+(c(84)*(u2937*(u2886*(u479*u2980+u3811*u2885+u2886*&
 &((u875)*u2886+u3997*u2885+u3870*u2980))+u957*u2885+u677*u2980)))
  u2937=-5.74643370347815916d-3*u3886
  u3886=2.29857348139126367d-1*u4057
  u4057=-4.48221828871296415d-1*u478
  u478=2.39051642064691421d0*u962
  u962=5.74643370347815916d-3*u4097
  u4097=-2.24110914435648207d-1*u3994
  u3994=1.14928674069563183d-2*u4096
  u4096=-7.47036381452160691d-2*u892
  u892=1.12055457217824104d0*u952
  u952=-7.47036381452160691d-2*u413
  u413=2.49752609199620544d7*u4064
  u4064=-1.30731366754128121d0*u937
  value_9_=value_9_+(c(84)*(u2886*(u3886+u2980*(u413*u2980+u4096)+u2885&
 &*(u4083*u2885+u4097)+u2886*(u2886*(u478+u2980*(u3801+u952)+u2885*(u39&
 &09+u3877)+(u456+u3877)*u2886)+u2885*(u3836+u818*u2885)+u2980*(u892+u4&
 &064*u2980)+u4057))+u2885*(u962+u604*u2885)+u2980*(u3994+u622*u2980)+u&
 &2937))
  if ( lmax .eq. 6 ) go to 100

  ! set computation flag
  D_X_Y_shell_optv1_4 = .false.

  ! setting results
100  continue
  value(1:1+9*(NVEC-1):9) = value_1_
  value(2:2+9*(NVEC-1):9) = value_2_
  value(3:3+9*(NVEC-1):9) = value_3_
  value(4:4+9*(NVEC-1):9) = value_4_
  value(5:5+9*(NVEC-1):9) = value_5_
  value(6:6+9*(NVEC-1):9) = value_6_
  value(7:7+9*(NVEC-1):9) = value_7_
  value(8:8+9*(NVEC-1):9) = value_8_
  value(9:9+9*(NVEC-1):9) = value_9_

  contains

    include "Av_functions.h"

end function

  
end module
