module mod_D_X_Y_shell_opt_4
  
  use mod_A_functions
  
  implicit none
  
  private
  public :: D_X_Y_shell_opt_4
  
contains
  
  

!> Compute D X Y/R coulomb integral for the right hand shell 
!!
recursive function D_X_Y_shell_opt_4(a12,a3,r3,c,lmax,value)
  
  use mod_A_functions
  
  implicit none
  
  ! input arguments
  real(kind=8), intent(in)  :: a12      !< combined exponent = a1*a2/(a1+a2) of the two first solid harmonics 
  real(kind=8), intent(in)  :: a3       !< exponent of the third solid harmonics 
  real(kind=8), intent(in)  :: r3(3)    !< third solid harmonic center, taken with combined center (a1*R1+a2*R2)/(a1+a2) as origin
  real(kind=8), intent(in)  :: c(*)     !< coeffcient of the hermit decomposition of the product of orbital 1 and 2
  integer     , intent(in)  :: lmax     !< l of maximum non 0 coeff in hermit decomposition~%
  real(kind=8), intent(out) :: value(*) !< result value
  
  ! return value
  
  logical :: D_X_Y_shell_opt_4
  ! local variables
  real(kind=8) :: d,p,pd,pd2,erfpd,exppd2,x,y,z
  real(kind=8) :: u2855
  real(kind=8) :: u2856
  real(kind=8) :: u2857
  real(kind=8) :: u2858
  real(kind=8) :: u2867
  real(kind=8) :: u2966
  real(kind=8) :: u2983
  real(kind=8) :: u2984
  real(kind=8) :: u2985
  real(kind=8) :: u3621
  real(kind=8) :: u3754
  real(kind=8) :: u3755
  real(kind=8) :: u3756
  real(kind=8) :: u3757
  real(kind=8) :: u3758
  real(kind=8) :: u3759
  real(kind=8) :: u376
  real(kind=8) :: u3760
  real(kind=8) :: u3761
  real(kind=8) :: u3762
  real(kind=8) :: u3763
  real(kind=8) :: u3764
  real(kind=8) :: u3765
  real(kind=8) :: u3766
  real(kind=8) :: u3767
  real(kind=8) :: u3768
  real(kind=8) :: u3769
  real(kind=8) :: u377
  real(kind=8) :: u3770
  real(kind=8) :: u3771
  real(kind=8) :: u3772
  real(kind=8) :: u3773
  real(kind=8) :: u3774
  real(kind=8) :: u3775
  real(kind=8) :: u3776
  real(kind=8) :: u3777
  real(kind=8) :: u3778
  real(kind=8) :: u3779
  real(kind=8) :: u378
  real(kind=8) :: u3780
  real(kind=8) :: u3781
  real(kind=8) :: u3782
  real(kind=8) :: u3783
  real(kind=8) :: u3784
  real(kind=8) :: u3785
  real(kind=8) :: u3786
  real(kind=8) :: u3787
  real(kind=8) :: u3788
  real(kind=8) :: u3789
  real(kind=8) :: u379
  real(kind=8) :: u3790
  real(kind=8) :: u3791
  real(kind=8) :: u3792
  real(kind=8) :: u3793
  real(kind=8) :: u3794
  real(kind=8) :: u3795
  real(kind=8) :: u3796
  real(kind=8) :: u3797
  real(kind=8) :: u3798
  real(kind=8) :: u3799
  real(kind=8) :: u38
  real(kind=8) :: u380
  real(kind=8) :: u3800
  real(kind=8) :: u3801
  real(kind=8) :: u3802
  real(kind=8) :: u3803
  real(kind=8) :: u3804
  real(kind=8) :: u3805
  real(kind=8) :: u3806
  real(kind=8) :: u3807
  real(kind=8) :: u3808
  real(kind=8) :: u3809
  real(kind=8) :: u381
  real(kind=8) :: u3810
  real(kind=8) :: u3811
  real(kind=8) :: u3812
  real(kind=8) :: u3813
  real(kind=8) :: u3814
  real(kind=8) :: u3815
  real(kind=8) :: u3816
  real(kind=8) :: u3817
  real(kind=8) :: u3818
  real(kind=8) :: u3819
  real(kind=8) :: u382
  real(kind=8) :: u3820
  real(kind=8) :: u3821
  real(kind=8) :: u3822
  real(kind=8) :: u3823
  real(kind=8) :: u3824
  real(kind=8) :: u3825
  real(kind=8) :: u3826
  real(kind=8) :: u3827
  real(kind=8) :: u3828
  real(kind=8) :: u3829
  real(kind=8) :: u383
  real(kind=8) :: u3830
  real(kind=8) :: u3831
  real(kind=8) :: u3832
  real(kind=8) :: u3833
  real(kind=8) :: u3834
  real(kind=8) :: u3835
  real(kind=8) :: u3836
  real(kind=8) :: u3837
  real(kind=8) :: u3838
  real(kind=8) :: u3839
  real(kind=8) :: u384
  real(kind=8) :: u3840
  real(kind=8) :: u3841
  real(kind=8) :: u3842
  real(kind=8) :: u3843
  real(kind=8) :: u3844
  real(kind=8) :: u3845
  real(kind=8) :: u3846
  real(kind=8) :: u3847
  real(kind=8) :: u3848
  real(kind=8) :: u3849
  real(kind=8) :: u385
  real(kind=8) :: u3850
  real(kind=8) :: u3851
  real(kind=8) :: u3852
  real(kind=8) :: u3853
  real(kind=8) :: u3854
  real(kind=8) :: u3855
  real(kind=8) :: u3856
  real(kind=8) :: u3857
  real(kind=8) :: u3858
  real(kind=8) :: u3859
  real(kind=8) :: u386
  real(kind=8) :: u3860
  real(kind=8) :: u3861
  real(kind=8) :: u3862
  real(kind=8) :: u3863
  real(kind=8) :: u3864
  real(kind=8) :: u3865
  real(kind=8) :: u3866
  real(kind=8) :: u3867
  real(kind=8) :: u3868
  real(kind=8) :: u3869
  real(kind=8) :: u387
  real(kind=8) :: u3870
  real(kind=8) :: u3871
  real(kind=8) :: u3872
  real(kind=8) :: u3873
  real(kind=8) :: u3874
  real(kind=8) :: u3875
  real(kind=8) :: u3876
  real(kind=8) :: u3877
  real(kind=8) :: u3878
  real(kind=8) :: u3879
  real(kind=8) :: u388
  real(kind=8) :: u3880
  real(kind=8) :: u3881
  real(kind=8) :: u3882
  real(kind=8) :: u3883
  real(kind=8) :: u3884
  real(kind=8) :: u3885
  real(kind=8) :: u3886
  real(kind=8) :: u3887
  real(kind=8) :: u3888
  real(kind=8) :: u3889
  real(kind=8) :: u389
  real(kind=8) :: u3890
  real(kind=8) :: u3891
  real(kind=8) :: u3892
  real(kind=8) :: u3893
  real(kind=8) :: u3894
  real(kind=8) :: u3895
  real(kind=8) :: u3896
  real(kind=8) :: u3897
  real(kind=8) :: u3898
  real(kind=8) :: u3899
  real(kind=8) :: u39
  real(kind=8) :: u390
  real(kind=8) :: u3900
  real(kind=8) :: u3901
  real(kind=8) :: u3902
  real(kind=8) :: u3903
  real(kind=8) :: u3904
  real(kind=8) :: u3905
  real(kind=8) :: u3906
  real(kind=8) :: u3907
  real(kind=8) :: u3908
  real(kind=8) :: u3909
  real(kind=8) :: u391
  real(kind=8) :: u3910
  real(kind=8) :: u3911
  real(kind=8) :: u3912
  real(kind=8) :: u3913
  real(kind=8) :: u3914
  real(kind=8) :: u3915
  real(kind=8) :: u3916
  real(kind=8) :: u3917
  real(kind=8) :: u3918
  real(kind=8) :: u3919
  real(kind=8) :: u392
  real(kind=8) :: u3920
  real(kind=8) :: u3921
  real(kind=8) :: u3922
  real(kind=8) :: u3923
  real(kind=8) :: u3924
  real(kind=8) :: u3925
  real(kind=8) :: u3926
  real(kind=8) :: u3927
  real(kind=8) :: u3928
  real(kind=8) :: u3929
  real(kind=8) :: u393
  real(kind=8) :: u3930
  real(kind=8) :: u3931
  real(kind=8) :: u3932
  real(kind=8) :: u3933
  real(kind=8) :: u3934
  real(kind=8) :: u3935
  real(kind=8) :: u3936
  real(kind=8) :: u3937
  real(kind=8) :: u3938
  real(kind=8) :: u3939
  real(kind=8) :: u394
  real(kind=8) :: u3940
  real(kind=8) :: u3941
  real(kind=8) :: u3942
  real(kind=8) :: u3943
  real(kind=8) :: u3944
  real(kind=8) :: u3945
  real(kind=8) :: u3946
  real(kind=8) :: u3947
  real(kind=8) :: u3948
  real(kind=8) :: u3949
  real(kind=8) :: u395
  real(kind=8) :: u3950
  real(kind=8) :: u3951
  real(kind=8) :: u3952
  real(kind=8) :: u3953
  real(kind=8) :: u3954
  real(kind=8) :: u3955
  real(kind=8) :: u3956
  real(kind=8) :: u3957
  real(kind=8) :: u3958
  real(kind=8) :: u3959
  real(kind=8) :: u396
  real(kind=8) :: u3960
  real(kind=8) :: u3961
  real(kind=8) :: u3962
  real(kind=8) :: u3963
  real(kind=8) :: u3964
  real(kind=8) :: u3965
  real(kind=8) :: u3966
  real(kind=8) :: u3967
  real(kind=8) :: u3968
  real(kind=8) :: u3969
  real(kind=8) :: u397
  real(kind=8) :: u3970
  real(kind=8) :: u3971
  real(kind=8) :: u3972
  real(kind=8) :: u3973
  real(kind=8) :: u3974
  real(kind=8) :: u3975
  real(kind=8) :: u3976
  real(kind=8) :: u3977
  real(kind=8) :: u3978
  real(kind=8) :: u3979
  real(kind=8) :: u398
  real(kind=8) :: u3980
  real(kind=8) :: u3981
  real(kind=8) :: u3982
  real(kind=8) :: u3983
  real(kind=8) :: u3984
  real(kind=8) :: u3985
  real(kind=8) :: u3986
  real(kind=8) :: u3987
  real(kind=8) :: u3988
  real(kind=8) :: u3989
  real(kind=8) :: u399
  real(kind=8) :: u3990
  real(kind=8) :: u3991
  real(kind=8) :: u3992
  real(kind=8) :: u3993
  real(kind=8) :: u3994
  real(kind=8) :: u3995
  real(kind=8) :: u3996
  real(kind=8) :: u3997
  real(kind=8) :: u3998
  real(kind=8) :: u3999
  real(kind=8) :: u4
  real(kind=8) :: u40
  real(kind=8) :: u400
  real(kind=8) :: u4000
  real(kind=8) :: u4001
  real(kind=8) :: u4002
  real(kind=8) :: u4003
  real(kind=8) :: u4004
  real(kind=8) :: u4005
  real(kind=8) :: u4006
  real(kind=8) :: u4007
  real(kind=8) :: u4008
  real(kind=8) :: u4009
  real(kind=8) :: u401
  real(kind=8) :: u4010
  real(kind=8) :: u4011
  real(kind=8) :: u4012
  real(kind=8) :: u4013
  real(kind=8) :: u4014
  real(kind=8) :: u4015
  real(kind=8) :: u4016
  real(kind=8) :: u4017
  real(kind=8) :: u4018
  real(kind=8) :: u4019
  real(kind=8) :: u402
  real(kind=8) :: u4020
  real(kind=8) :: u4021
  real(kind=8) :: u4022
  real(kind=8) :: u4023
  real(kind=8) :: u4024
  real(kind=8) :: u4025
  real(kind=8) :: u4026
  real(kind=8) :: u4027
  real(kind=8) :: u4028
  real(kind=8) :: u4029
  real(kind=8) :: u403
  real(kind=8) :: u4030
  real(kind=8) :: u4031
  real(kind=8) :: u4032
  real(kind=8) :: u4033
  real(kind=8) :: u4034
  real(kind=8) :: u4035
  real(kind=8) :: u4036
  real(kind=8) :: u4037
  real(kind=8) :: u4038
  real(kind=8) :: u4039
  real(kind=8) :: u404
  real(kind=8) :: u4040
  real(kind=8) :: u4041
  real(kind=8) :: u4042
  real(kind=8) :: u4043
  real(kind=8) :: u4044
  real(kind=8) :: u4045
  real(kind=8) :: u4046
  real(kind=8) :: u4047
  real(kind=8) :: u4048
  real(kind=8) :: u4049
  real(kind=8) :: u405
  real(kind=8) :: u4050
  real(kind=8) :: u4051
  real(kind=8) :: u4052
  real(kind=8) :: u4053
  real(kind=8) :: u4054
  real(kind=8) :: u4055
  real(kind=8) :: u4056
  real(kind=8) :: u4057
  real(kind=8) :: u4058
  real(kind=8) :: u4059
  real(kind=8) :: u406
  real(kind=8) :: u4060
  real(kind=8) :: u4061
  real(kind=8) :: u4062
  real(kind=8) :: u4063
  real(kind=8) :: u4064
  real(kind=8) :: u4065
  real(kind=8) :: u4066
  real(kind=8) :: u4067
  real(kind=8) :: u4068
  real(kind=8) :: u4069
  real(kind=8) :: u407
  real(kind=8) :: u4070
  real(kind=8) :: u4071
  real(kind=8) :: u4072
  real(kind=8) :: u4073
  real(kind=8) :: u4074
  real(kind=8) :: u4075
  real(kind=8) :: u4076
  real(kind=8) :: u4077
  real(kind=8) :: u4078
  real(kind=8) :: u4079
  real(kind=8) :: u408
  real(kind=8) :: u4080
  real(kind=8) :: u4081
  real(kind=8) :: u4082
  real(kind=8) :: u4083
  real(kind=8) :: u4084
  real(kind=8) :: u4085
  real(kind=8) :: u4086
  real(kind=8) :: u4087
  real(kind=8) :: u4088
  real(kind=8) :: u4089
  real(kind=8) :: u409
  real(kind=8) :: u4090
  real(kind=8) :: u4091
  real(kind=8) :: u4092
  real(kind=8) :: u4093
  real(kind=8) :: u4094
  real(kind=8) :: u4095
  real(kind=8) :: u4096
  real(kind=8) :: u4097
  real(kind=8) :: u4098
  real(kind=8) :: u4099
  real(kind=8) :: u41
  real(kind=8) :: u410
  real(kind=8) :: u4100
  real(kind=8) :: u4101
  real(kind=8) :: u4102
  real(kind=8) :: u4103
  real(kind=8) :: u4104
  real(kind=8) :: u4105
  real(kind=8) :: u4106
  real(kind=8) :: u4107
  real(kind=8) :: u4108
  real(kind=8) :: u4109
  real(kind=8) :: u411
  real(kind=8) :: u4110
  real(kind=8) :: u4111
  real(kind=8) :: u4112
  real(kind=8) :: u4113
  real(kind=8) :: u4114
  real(kind=8) :: u4115
  real(kind=8) :: u4116
  real(kind=8) :: u4117
  real(kind=8) :: u4118
  real(kind=8) :: u4119
  real(kind=8) :: u412
  real(kind=8) :: u4120
  real(kind=8) :: u4121
  real(kind=8) :: u4122
  real(kind=8) :: u4123
  real(kind=8) :: u4124
  real(kind=8) :: u4125
  real(kind=8) :: u4126
  real(kind=8) :: u4127
  real(kind=8) :: u4128
  real(kind=8) :: u4129
  real(kind=8) :: u413
  real(kind=8) :: u4130
  real(kind=8) :: u4131
  real(kind=8) :: u4132
  real(kind=8) :: u4133
  real(kind=8) :: u4134
  real(kind=8) :: u4135
  real(kind=8) :: u4136
  real(kind=8) :: u4137
  real(kind=8) :: u4138
  real(kind=8) :: u4139
  real(kind=8) :: u414
  real(kind=8) :: u4140
  real(kind=8) :: u4141
  real(kind=8) :: u4142
  real(kind=8) :: u4143
  real(kind=8) :: u4144
  real(kind=8) :: u4145
  real(kind=8) :: u4146
  real(kind=8) :: u4147
  real(kind=8) :: u4148
  real(kind=8) :: u4149
  real(kind=8) :: u415
  real(kind=8) :: u4150
  real(kind=8) :: u4151
  real(kind=8) :: u4152
  real(kind=8) :: u4153
  real(kind=8) :: u4154
  real(kind=8) :: u4155
  real(kind=8) :: u4156
  real(kind=8) :: u4157
  real(kind=8) :: u4158
  real(kind=8) :: u4159
  real(kind=8) :: u416
  real(kind=8) :: u4160
  real(kind=8) :: u4161
  real(kind=8) :: u4162
  real(kind=8) :: u4163
  real(kind=8) :: u4164
  real(kind=8) :: u4165
  real(kind=8) :: u4166
  real(kind=8) :: u4167
  real(kind=8) :: u4169
  real(kind=8) :: u417
  real(kind=8) :: u4170
  real(kind=8) :: u4171
  real(kind=8) :: u4172
  real(kind=8) :: u4173
  real(kind=8) :: u4174
  real(kind=8) :: u4175
  real(kind=8) :: u4176
  real(kind=8) :: u4177
  real(kind=8) :: u4178
  real(kind=8) :: u4179
  real(kind=8) :: u418
  real(kind=8) :: u4180
  real(kind=8) :: u4181
  real(kind=8) :: u4182
  real(kind=8) :: u4183
  real(kind=8) :: u4184
  real(kind=8) :: u4185
  real(kind=8) :: u4186
  real(kind=8) :: u4187
  real(kind=8) :: u4188
  real(kind=8) :: u4189
  real(kind=8) :: u419
  real(kind=8) :: u4190
  real(kind=8) :: u4191
  real(kind=8) :: u4192
  real(kind=8) :: u4193
  real(kind=8) :: u4194
  real(kind=8) :: u4195
  real(kind=8) :: u4196
  real(kind=8) :: u4197
  real(kind=8) :: u4198
  real(kind=8) :: u4199
  real(kind=8) :: u42
  real(kind=8) :: u420
  real(kind=8) :: u4200
  real(kind=8) :: u4201
  real(kind=8) :: u4202
  real(kind=8) :: u4203
  real(kind=8) :: u4204
  real(kind=8) :: u4205
  real(kind=8) :: u4206
  real(kind=8) :: u4207
  real(kind=8) :: u4208
  real(kind=8) :: u4209
  real(kind=8) :: u421
  real(kind=8) :: u4210
  real(kind=8) :: u4211
  real(kind=8) :: u4212
  real(kind=8) :: u4213
  real(kind=8) :: u4214
  real(kind=8) :: u4215
  real(kind=8) :: u4216
  real(kind=8) :: u4217
  real(kind=8) :: u4218
  real(kind=8) :: u4219
  real(kind=8) :: u422
  real(kind=8) :: u4220
  real(kind=8) :: u4221
  real(kind=8) :: u4222
  real(kind=8) :: u4223
  real(kind=8) :: u4224
  real(kind=8) :: u4225
  real(kind=8) :: u423
  real(kind=8) :: u424
  real(kind=8) :: u425
  real(kind=8) :: u426
  real(kind=8) :: u427
  real(kind=8) :: u428
  real(kind=8) :: u429
  real(kind=8) :: u43
  real(kind=8) :: u430
  real(kind=8) :: u431
  real(kind=8) :: u432
  real(kind=8) :: u433
  real(kind=8) :: u434
  real(kind=8) :: u435
  real(kind=8) :: u436
  real(kind=8) :: u437
  real(kind=8) :: u438
  real(kind=8) :: u439
  real(kind=8) :: u44
  real(kind=8) :: u440
  real(kind=8) :: u441
  real(kind=8) :: u442
  real(kind=8) :: u443
  real(kind=8) :: u444
  real(kind=8) :: u445
  real(kind=8) :: u446
  real(kind=8) :: u447
  real(kind=8) :: u448
  real(kind=8) :: u449
  real(kind=8) :: u45
  real(kind=8) :: u450
  real(kind=8) :: u451
  real(kind=8) :: u452
  real(kind=8) :: u453
  real(kind=8) :: u454
  real(kind=8) :: u455
  real(kind=8) :: u456
  real(kind=8) :: u457
  real(kind=8) :: u458
  real(kind=8) :: u459
  real(kind=8) :: u46
  real(kind=8) :: u460
  real(kind=8) :: u461
  real(kind=8) :: u462
  real(kind=8) :: u463
  real(kind=8) :: u464
  real(kind=8) :: u465
  real(kind=8) :: u466
  real(kind=8) :: u467
  real(kind=8) :: u468
  real(kind=8) :: u469
  real(kind=8) :: u47
  real(kind=8) :: u470
  real(kind=8) :: u471
  real(kind=8) :: u472
  real(kind=8) :: u473
  real(kind=8) :: u474
  real(kind=8) :: u475
  real(kind=8) :: u476
  real(kind=8) :: u477
  real(kind=8) :: u478
  real(kind=8) :: u479
  real(kind=8) :: u48
  real(kind=8) :: u480
  real(kind=8) :: u481
  real(kind=8) :: u482
  real(kind=8) :: u483
  real(kind=8) :: u484
  real(kind=8) :: u485
  real(kind=8) :: u486
  real(kind=8) :: u487
  real(kind=8) :: u488
  real(kind=8) :: u489
  real(kind=8) :: u49
  real(kind=8) :: u490
  real(kind=8) :: u491
  real(kind=8) :: u492
  real(kind=8) :: u493
  real(kind=8) :: u494
  real(kind=8) :: u495
  real(kind=8) :: u496
  real(kind=8) :: u497
  real(kind=8) :: u498
  real(kind=8) :: u499
  real(kind=8) :: u5
  real(kind=8) :: u50
  real(kind=8) :: u500
  real(kind=8) :: u501
  real(kind=8) :: u502
  real(kind=8) :: u503
  real(kind=8) :: u504
  real(kind=8) :: u505
  real(kind=8) :: u506
  real(kind=8) :: u507
  real(kind=8) :: u508
  real(kind=8) :: u509
  real(kind=8) :: u51
  real(kind=8) :: u510
  real(kind=8) :: u511
  real(kind=8) :: u512
  real(kind=8) :: u513
  real(kind=8) :: u514
  real(kind=8) :: u515
  real(kind=8) :: u516
  real(kind=8) :: u517
  real(kind=8) :: u518
  real(kind=8) :: u519
  real(kind=8) :: u52
  real(kind=8) :: u520
  real(kind=8) :: u521
  real(kind=8) :: u522
  real(kind=8) :: u523
  real(kind=8) :: u524
  real(kind=8) :: u525
  real(kind=8) :: u526
  real(kind=8) :: u527
  real(kind=8) :: u528
  real(kind=8) :: u529
  real(kind=8) :: u53
  real(kind=8) :: u530
  real(kind=8) :: u531
  real(kind=8) :: u532
  real(kind=8) :: u533
  real(kind=8) :: u534
  real(kind=8) :: u535
  real(kind=8) :: u536
  real(kind=8) :: u537
  real(kind=8) :: u538
  real(kind=8) :: u539
  real(kind=8) :: u54
  real(kind=8) :: u540
  real(kind=8) :: u541
  real(kind=8) :: u542
  real(kind=8) :: u543
  real(kind=8) :: u544
  real(kind=8) :: u545
  real(kind=8) :: u546
  real(kind=8) :: u547
  real(kind=8) :: u548
  real(kind=8) :: u549
  real(kind=8) :: u55
  real(kind=8) :: u550
  real(kind=8) :: u551
  real(kind=8) :: u552
  real(kind=8) :: u553
  real(kind=8) :: u554
  real(kind=8) :: u555
  real(kind=8) :: u556
  real(kind=8) :: u557
  real(kind=8) :: u558
  real(kind=8) :: u559
  real(kind=8) :: u56
  real(kind=8) :: u560
  real(kind=8) :: u561
  real(kind=8) :: u562
  real(kind=8) :: u563
  real(kind=8) :: u564
  real(kind=8) :: u565
  real(kind=8) :: u566
  real(kind=8) :: u567
  real(kind=8) :: u568
  real(kind=8) :: u569
  real(kind=8) :: u57
  real(kind=8) :: u570
  real(kind=8) :: u571
  real(kind=8) :: u572
  real(kind=8) :: u573
  real(kind=8) :: u574
  real(kind=8) :: u575
  real(kind=8) :: u576
  real(kind=8) :: u577
  real(kind=8) :: u578
  real(kind=8) :: u579
  real(kind=8) :: u58
  real(kind=8) :: u580
  real(kind=8) :: u581
  real(kind=8) :: u582
  real(kind=8) :: u583
  real(kind=8) :: u584
  real(kind=8) :: u585
  real(kind=8) :: u586
  real(kind=8) :: u587
  real(kind=8) :: u588
  real(kind=8) :: u589
  real(kind=8) :: u59
  real(kind=8) :: u590
  real(kind=8) :: u591
  real(kind=8) :: u592
  real(kind=8) :: u593
  real(kind=8) :: u594
  real(kind=8) :: u595
  real(kind=8) :: u596
  real(kind=8) :: u597
  real(kind=8) :: u598
  real(kind=8) :: u599
  real(kind=8) :: u6
  real(kind=8) :: u60
  real(kind=8) :: u600
  real(kind=8) :: u601
  real(kind=8) :: u602
  real(kind=8) :: u603
  real(kind=8) :: u604
  real(kind=8) :: u605
  real(kind=8) :: u606
  real(kind=8) :: u607
  real(kind=8) :: u608
  real(kind=8) :: u609
  real(kind=8) :: u61
  real(kind=8) :: u610
  real(kind=8) :: u611
  real(kind=8) :: u612
  real(kind=8) :: u613
  real(kind=8) :: u614
  real(kind=8) :: u615
  real(kind=8) :: u616
  real(kind=8) :: u617
  real(kind=8) :: u618
  real(kind=8) :: u619
  real(kind=8) :: u62
  real(kind=8) :: u620
  real(kind=8) :: u621
  real(kind=8) :: u622
  real(kind=8) :: u623
  real(kind=8) :: u624
  real(kind=8) :: u625
  real(kind=8) :: u626
  real(kind=8) :: u627
  real(kind=8) :: u628
  real(kind=8) :: u629
  real(kind=8) :: u63
  real(kind=8) :: u630
  real(kind=8) :: u631
  real(kind=8) :: u632
  real(kind=8) :: u633
  real(kind=8) :: u634
  real(kind=8) :: u635
  real(kind=8) :: u636
  real(kind=8) :: u637
  real(kind=8) :: u638
  real(kind=8) :: u639
  real(kind=8) :: u64
  real(kind=8) :: u640
  real(kind=8) :: u641
  real(kind=8) :: u642
  real(kind=8) :: u643
  real(kind=8) :: u644
  real(kind=8) :: u645
  real(kind=8) :: u646
  real(kind=8) :: u647
  real(kind=8) :: u648
  real(kind=8) :: u649
  real(kind=8) :: u65
  real(kind=8) :: u650
  real(kind=8) :: u651
  real(kind=8) :: u653
  real(kind=8) :: u654
  real(kind=8) :: u655
  real(kind=8) :: u656
  real(kind=8) :: u657
  real(kind=8) :: u658
  real(kind=8) :: u659
  real(kind=8) :: u66
  real(kind=8) :: u660
  real(kind=8) :: u661
  real(kind=8) :: u662
  real(kind=8) :: u663
  real(kind=8) :: u664
  real(kind=8) :: u666
  real(kind=8) :: u667
  real(kind=8) :: u668
  real(kind=8) :: u67
  real(kind=8) :: u670
  real(kind=8) :: u671
  real(kind=8) :: u672
  real(kind=8) :: u673
  real(kind=8) :: u674
  real(kind=8) :: u675
  real(kind=8) :: u676
  real(kind=8) :: u677
  real(kind=8) :: u678
  real(kind=8) :: u679
  real(kind=8) :: u68
  real(kind=8) :: u680
  real(kind=8) :: u681
  real(kind=8) :: u682
  real(kind=8) :: u683
  real(kind=8) :: u684
  real(kind=8) :: u685
  real(kind=8) :: u686
  real(kind=8) :: u687
  real(kind=8) :: u688
  real(kind=8) :: u689
  real(kind=8) :: u69
  real(kind=8) :: u690
  real(kind=8) :: u691
  real(kind=8) :: u692
  real(kind=8) :: u693
  real(kind=8) :: u694
  real(kind=8) :: u695
  real(kind=8) :: u696
  real(kind=8) :: u697
  real(kind=8) :: u698
  real(kind=8) :: u699
  real(kind=8) :: u7
  real(kind=8) :: u70
  real(kind=8) :: u700
  real(kind=8) :: u701
  real(kind=8) :: u702
  real(kind=8) :: u703
  real(kind=8) :: u704
  real(kind=8) :: u705
  real(kind=8) :: u706
  real(kind=8) :: u707
  real(kind=8) :: u708
  real(kind=8) :: u709
  real(kind=8) :: u71
  real(kind=8) :: u710
  real(kind=8) :: u711
  real(kind=8) :: u712
  real(kind=8) :: u713
  real(kind=8) :: u714
  real(kind=8) :: u715
  real(kind=8) :: u716
  real(kind=8) :: u717
  real(kind=8) :: u718
  real(kind=8) :: u719
  real(kind=8) :: u72
  real(kind=8) :: u720
  real(kind=8) :: u721
  real(kind=8) :: u722
  real(kind=8) :: u723
  real(kind=8) :: u724
  real(kind=8) :: u725
  real(kind=8) :: u726
  real(kind=8) :: u727
  real(kind=8) :: u728
  real(kind=8) :: u729
  real(kind=8) :: u73
  real(kind=8) :: u731
  real(kind=8) :: u732
  real(kind=8) :: u733
  real(kind=8) :: u734
  real(kind=8) :: u735
  real(kind=8) :: u736
  real(kind=8) :: u737
  real(kind=8) :: u738
  real(kind=8) :: u739
  real(kind=8) :: u74
  real(kind=8) :: u740
  real(kind=8) :: u741
  real(kind=8) :: u743
  real(kind=8) :: u744
  real(kind=8) :: u746
  real(kind=8) :: u747
  real(kind=8) :: u748
  real(kind=8) :: u749
  real(kind=8) :: u75
  real(kind=8) :: u750
  real(kind=8) :: u751
  real(kind=8) :: u752
  real(kind=8) :: u753
  real(kind=8) :: u754
  real(kind=8) :: u755
  real(kind=8) :: u756
  real(kind=8) :: u757
  real(kind=8) :: u758
  real(kind=8) :: u759
  real(kind=8) :: u76
  real(kind=8) :: u760
  real(kind=8) :: u761
  real(kind=8) :: u762
  real(kind=8) :: u763
  real(kind=8) :: u764
  real(kind=8) :: u765
  real(kind=8) :: u766
  real(kind=8) :: u767
  real(kind=8) :: u768
  real(kind=8) :: u769
  real(kind=8) :: u77
  real(kind=8) :: u770
  real(kind=8) :: u771
  real(kind=8) :: u772
  real(kind=8) :: u773
  real(kind=8) :: u774
  real(kind=8) :: u775
  real(kind=8) :: u776
  real(kind=8) :: u777
  real(kind=8) :: u778
  real(kind=8) :: u779
  real(kind=8) :: u78
  real(kind=8) :: u780
  real(kind=8) :: u781
  real(kind=8) :: u782
  real(kind=8) :: u783
  real(kind=8) :: u784
  real(kind=8) :: u785
  real(kind=8) :: u786
  real(kind=8) :: u787
  real(kind=8) :: u788
  real(kind=8) :: u789
  real(kind=8) :: u79
  real(kind=8) :: u790
  real(kind=8) :: u791
  real(kind=8) :: u792
  real(kind=8) :: u793
  real(kind=8) :: u794
  real(kind=8) :: u795
  real(kind=8) :: u796
  real(kind=8) :: u797
  real(kind=8) :: u798
  real(kind=8) :: u799
  real(kind=8) :: u8
  real(kind=8) :: u80
  real(kind=8) :: u800
  real(kind=8) :: u801
  real(kind=8) :: u802
  real(kind=8) :: u803
  real(kind=8) :: u804
  real(kind=8) :: u805
  real(kind=8) :: u806
  real(kind=8) :: u807
  real(kind=8) :: u808
  real(kind=8) :: u809
  real(kind=8) :: u81
  real(kind=8) :: u810
  real(kind=8) :: u811
  real(kind=8) :: u812
  real(kind=8) :: u813
  real(kind=8) :: u814
  real(kind=8) :: u815
  real(kind=8) :: u816
  real(kind=8) :: u817
  real(kind=8) :: u818
  real(kind=8) :: u819
  real(kind=8) :: u82
  real(kind=8) :: u820
  real(kind=8) :: u821
  real(kind=8) :: u822
  real(kind=8) :: u823
  real(kind=8) :: u824
  real(kind=8) :: u825
  real(kind=8) :: u826
  real(kind=8) :: u827
  real(kind=8) :: u828
  real(kind=8) :: u829
  real(kind=8) :: u83
  real(kind=8) :: u830
  real(kind=8) :: u831
  real(kind=8) :: u832
  real(kind=8) :: u833
  real(kind=8) :: u834
  real(kind=8) :: u835
  real(kind=8) :: u836
  real(kind=8) :: u837
  real(kind=8) :: u838
  real(kind=8) :: u839
  real(kind=8) :: u84
  real(kind=8) :: u840
  real(kind=8) :: u841
  real(kind=8) :: u842
  real(kind=8) :: u843
  real(kind=8) :: u844
  real(kind=8) :: u845
  real(kind=8) :: u847
  real(kind=8) :: u848
  real(kind=8) :: u849
  real(kind=8) :: u85
  real(kind=8) :: u850
  real(kind=8) :: u851
  real(kind=8) :: u852
  real(kind=8) :: u853
  real(kind=8) :: u854
  real(kind=8) :: u855
  real(kind=8) :: u856
  real(kind=8) :: u857
  real(kind=8) :: u858
  real(kind=8) :: u859
  real(kind=8) :: u86
  real(kind=8) :: u860
  real(kind=8) :: u861
  real(kind=8) :: u862
  real(kind=8) :: u863
  real(kind=8) :: u864
  real(kind=8) :: u865
  real(kind=8) :: u866
  real(kind=8) :: u867
  real(kind=8) :: u868
  real(kind=8) :: u869
  real(kind=8) :: u87
  real(kind=8) :: u870
  real(kind=8) :: u871
  real(kind=8) :: u872
  real(kind=8) :: u873
  real(kind=8) :: u874
  real(kind=8) :: u875
  real(kind=8) :: u876
  real(kind=8) :: u877
  real(kind=8) :: u878
  real(kind=8) :: u879
  real(kind=8) :: u88
  real(kind=8) :: u880
  real(kind=8) :: u881
  real(kind=8) :: u882
  real(kind=8) :: u883
  real(kind=8) :: u884
  real(kind=8) :: u885
  real(kind=8) :: u886
  real(kind=8) :: u887
  real(kind=8) :: u888
  real(kind=8) :: u889
  real(kind=8) :: u89
  real(kind=8) :: u890
  real(kind=8) :: u891
  real(kind=8) :: u892
  real(kind=8) :: u893
  real(kind=8) :: u894
  real(kind=8) :: u895
  real(kind=8) :: u896
  real(kind=8) :: u897
  real(kind=8) :: u898
  real(kind=8) :: u899
  real(kind=8) :: u9
  real(kind=8) :: u90
  real(kind=8) :: u900
  real(kind=8) :: u901
  real(kind=8) :: u902
  real(kind=8) :: u903
  real(kind=8) :: u904
  real(kind=8) :: u905
  real(kind=8) :: u906
  real(kind=8) :: u907
  real(kind=8) :: u908
  real(kind=8) :: u909
  real(kind=8) :: u91
  real(kind=8) :: u910
  real(kind=8) :: u911
  real(kind=8) :: u912
  real(kind=8) :: u913
  real(kind=8) :: u914
  real(kind=8) :: u915
  real(kind=8) :: u916
  real(kind=8) :: u917
  real(kind=8) :: u918
  real(kind=8) :: u919
  real(kind=8) :: u92
  real(kind=8) :: u920
  real(kind=8) :: u921
  real(kind=8) :: u922
  real(kind=8) :: u923
  real(kind=8) :: u924
  real(kind=8) :: u925
  real(kind=8) :: u926
  real(kind=8) :: u927
  real(kind=8) :: u928
  real(kind=8) :: u929
  real(kind=8) :: u93
  real(kind=8) :: u930
  real(kind=8) :: u931
  real(kind=8) :: u932
  real(kind=8) :: u933
  real(kind=8) :: u934
  real(kind=8) :: u935
  real(kind=8) :: u936
  real(kind=8) :: u937
  real(kind=8) :: u938
  real(kind=8) :: u939
  real(kind=8) :: u94
  real(kind=8) :: u940
  real(kind=8) :: u941
  real(kind=8) :: u942
  real(kind=8) :: u943
  real(kind=8) :: u944
  real(kind=8) :: u945
  real(kind=8) :: u946
  real(kind=8) :: u947
  real(kind=8) :: u948
  real(kind=8) :: u949
  real(kind=8) :: u95
  real(kind=8) :: u950
  real(kind=8) :: u951
  real(kind=8) :: u952
  real(kind=8) :: u953
  real(kind=8) :: u954
  real(kind=8) :: u955
  real(kind=8) :: u956
  real(kind=8) :: u957
  real(kind=8) :: u958
  real(kind=8) :: u959
  real(kind=8) :: u96
  real(kind=8) :: u960
  real(kind=8) :: u961
  real(kind=8) :: u962
  real(kind=8) :: u963
  real(kind=8) :: u964
  real(kind=8) :: u965
  real(kind=8) :: u966
  real(kind=8) :: u967
  real(kind=8) :: u968
  real(kind=8) :: u969
  real(kind=8) :: u97
  real(kind=8) :: u970
  real(kind=8) :: u971
  real(kind=8) :: u972
  real(kind=8) :: u973
  real(kind=8) :: u974
  real(kind=8) :: u975
  real(kind=8) :: u976
  real(kind=8) :: u977
  real(kind=8) :: u978
  real(kind=8) :: u979
  real(kind=8) :: u98
  real(kind=8) :: u980
  real(kind=8) :: u981
  real(kind=8) :: u982
  real(kind=8) :: u983
  real(kind=8) :: u984
  real(kind=8) :: u985
  real(kind=8) :: u986
  real(kind=8) :: u987
  real(kind=8) :: u988
  real(kind=8) :: u989
  real(kind=8) :: u99
  real(kind=8) :: u990
  real(kind=8) :: u991
  real(kind=8) :: u992
  real(kind=8) :: u993
  real(kind=8) :: u994
  real(kind=8) :: u995
  real(kind=8) :: u996
  real(kind=8) :: u997
  real(kind=8) :: u998
  real(kind=8) :: u999

  
  ! init computation flag
  D_X_Y_shell_opt_4=.true.
  
  ! compute local quantities
  p=sqrt((a12*a3)/(a12+a3))
  d=sqrt(sum(r3**2))
  
  ! renormalize xyz
  x=0.0d0
  y=0.0d0
  z=0.0d0
  if (d>0.0d0) then
    x=r3(1)/d
    y=r3(2)/d
    z=r3(3)/d
  end if
  
  ! compute local quantities
  pd=p*d
  pd2=pd*pd
  erfpd =erf(pd)
  exppd2=exp(-pd2)
  
  ! computation section
  value(1:9)=0.0d0
  u2858=a3**(-4)
  u3621=exppd2*pd2
  u2966=-2.83592616144882564d1*u3621
  u3797=u2858*(3.2986722862692829d2*A5_(p,pd,erfpd,exppd2)+u2966)*p**5
  u382=-4.98024254301440461d-2*u3797
  u977=4.98024254301440461d-2*u3797
  u2855=x**2
  u2856=y**2
  u2983=x*y
  value(1)=value(1)+(c(1)*(u2983*(u382*u2856+u977*u2855)))
  u382=-3.5215632741192216d-2*u3797
  u977=1.05646898223576648d-1*u3797
  value(2)=value(2)+(c(1)*(y*(u382*u2856+u977*u2855)*z))
  u382=1.12941284893714386d-1*u3797
  u977=-1.8823547482285731d-2*u3797
  u2857=z**2
  value(3)=value(3)+(c(1)*(u2983*(u382*u2857+u977*u2856+u977*u2855)))
  u382=5.32410322828448158d-2*u3797
  u977=-3.99307742121336118d-2*u3797
  u2984=y*z
  u4011=(u382*u2857+u977*u2856+u977*u2855)
  value(4)=value(4)+(c(1)*(u2984*u4011))
  u2985=1.68362926992343652d-2*u3797
  u500=-5.05088780977030955d-2*u3797
  u3825=6.31360976221288694d-3*u3797
  u3774=1.26272195244257739d-2*u3797
  u2867=y**4
  value(5)=value(5)+(c(1)*(u2857*(u500*u2855+u500*u2856+u2985*u2857)+u3&
 &825*u2867+u2855*(u3825*u2855+u3774*u2856)))
  u2985=x*z
  value(6)=value(6)+(c(1)*(u2985*u4011))
  u3774=-5.6470642446857193d-2*u3797
  u977=9.41177374114286549d-3*u3797
  u500=5.6470642446857193d-2*u3797
  u382=-9.41177374114286549d-3*u3797
  u999=x**4
  value(7)=value(7)+(c(1)*((u500*u2855+u3774*u2856)*u2857+u977*u2867+u3&
 &82*u999))
  u977=-1.05646898223576648d-1*u3797
  u3774=3.5215632741192216d-2*u3797
  value(8)=value(8)+(c(1)*(x*(u977*u2856+u3774*u2855)*z))
  u977=1.24506063575360115d-2*u3797
  u3774=-7.47036381452160691d-2*u3797
  value(9)=value(9)+(c(1)*(u977*u2867+u2855*(u977*u2855+u3774*u2856)))
  if ( lmax .eq. 0 ) return
  u3774=A6_(p,pd,erfpd,exppd2)
  u977=p**6
  u382=exppd2*pd
  u500=-2.83592616144882564d1*u382
  u3797=u2858*(3.2986722862692829d2*u3774+u500)*u977
  u3825=4.98024254301440461d-2*u3797
  u4011=-1.49407276290432138d-1*u3797
  u926=2.0d0*pd2
  u4072=(9.0d0+u926)
  u4055=u2858*(2.96880505764235461d3*u3774+u500*u4072)*u977
  u3774=-4.98024254301440461d-2*u4055
  u977=4.98024254301440461d-2*u4055
  u3872=u977*u2855
  value(1)=value(1)+(c(2)*(y*((u3825+u3774*u2855)*u2856+u2855*(u4011+u3&
 &872))))
  u3825=-2.11293796447153296d-1*u3797
  u4011=-3.5215632741192216d-2*u4055
  u737=1.05646898223576648d-1*u4055
  u3990=u4011*u2856
  u3991=u737*u2855
  value(2)=value(2)+(c(2)*(u2983*(u3990+u3991+u3825)*z))
  u3825=-1.12941284893714386d-1*u3797
  u4037=1.8823547482285731d-2*u3797
  u3999=5.6470642446857193d-2*u3797
  u3876=1.12941284893714386d-1*u4055
  u4035=-1.8823547482285731d-2*u4055
  u3873=u4035*u2855
  u621=u3999+u3873
  value(3)=value(3)+(c(2)*(y*((u3825+u3876*u2855)*u2857+(u4037+u3873)*u&
 &2856+u2855*(u621))))
  u4000=7.98615484242672237d-2*u3797
  u950=5.32410322828448158d-2*u4055
  u780=-3.99307742121336118d-2*u4055
  u528=u2983*z
  u3824=u780*u2855
  u3992=u950*u2857
  u3993=u780*u2856
  u712=(u528*(u3992+u3993+u3824+u4000))
  value(4)=value(4)+(c(2)*u712)
  u990=1.01017756195406191d-1*u3797
  u3877=1.68362926992343652d-2*u4055
  u4002=-2.52544390488515477d-2*u3797
  u599=-5.05088780977030955d-2*u4055
  u513=6.31360976221288694d-3*u4055
  u504=1.26272195244257739d-2*u4055
  u3830=u599*u2855+u599*u2856+u3877*u2857
  u3878=u504*u2855+u513*u2856
  u3994=u513*u2855
  u4036=(u2857*(u990+u3830)+u2856*(u4002+u3878)+u2855*(u4002+u3994))
  value(5)=value(5)+(c(2)*(x*u4036))
  u668=-5.32410322828448158d-2*u3797
  u3768=3.99307742121336118d-2*u3797
  u4001=1.19792322636400836d-1*u3797
  value(6)=value(6)+(c(2)*(z*((u668+u950*u2855)*u2857+(u3768+u3824)*u28&
 &56+u2855*(u4001+u3824))))
  u3840=-5.6470642446857193d-2*u4055
  u638=9.41177374114286549d-3*u4055
  u741=3.7647094964571462d-2*u3797
  u969=5.6470642446857193d-2*u4055
  u4017=-9.41177374114286549d-3*u4055
  u3754=u969*u2855+u3840*u2856
  u3995=u4017*u2855
  value(7)=value(7)+(c(2)*(x*((u3825+u3754)*u2857+u638*u2867+u2855*(u74&
 &1+u3995))))
  u741=1.05646898223576648d-1*u3797
  u3869=-1.05646898223576648d-1*u3797
  u677=-1.05646898223576648d-1*u4055
  u3851=3.5215632741192216d-2*u4055
  u3874=u3851*u2855
  value(8)=value(8)+(c(2)*(((u741+u677*u2855)*u2856+u2855*(u3869+u3874)&
 &)*z))
  u3786=1.49407276290432138d-1*u3797
  u609=1.24506063575360115d-2*u4055
  u4013=-4.98024254301440461d-2*u3797
  u3778=-7.47036381452160691d-2*u4055
  u4055=u3778*u2855+u609*u2856
  u3875=u609*u2855
  value(9)=value(9)+(c(2)*(x*(u2856*(u3786+u4055)+u2855*(u4013+u3875)))&
 &)
  u3996=u3774*u2856
  value(1)=value(1)+(c(3)*(x*(u2856*(u3786+u3872+u3996)+u4013*u2855)))
  u470=u3991+u3990
  u3997=u3869*u2855
  value(2)=value(2)+(c(3)*((u2856*(u741+u470)+u3997)*z))
  u3998=u4035*u2856
  value(3)=value(3)+(c(3)*(x*((u3825+u3876*u2856)*u2857+u2856*(u621+u39&
 &98)+u4037*u2855)))
  u4037=u3824+u3993
  u3999=u3768*u2855
  value(4)=value(4)+(c(3)*(z*((u668+u950*u2856)*u2857+u2856*(u4001+u403&
 &7)+u3999)))
  value(5)=value(5)+(c(3)*(y*u4036))
  value(6)=value(6)+(c(3)*u712)
  u4002=1.12941284893714386d-1*u3797
  u668=-3.7647094964571462d-2*u3797
  u4000=u638*u2856
  value(7)=value(7)+(c(3)*(y*((u4002+u3754)*u2857+u2856*(u668+u4000)+u4&
 &017*u999)))
  u668=2.11293796447153296d-1*u3797
  u4001=u677*u2856
  value(8)=value(8)+(c(3)*(u2983*(u4001+u3874+u668)*z))
  value(9)=value(9)+(c(3)*(y*(u2856*(u4013+u4055)+u2855*(u3786+u3875)))&
 &)
  value(1)=value(1)+(c(4)*(u2983*(u3996+u3872)*z))
  u3774=3.5215632741192216d-2*u3797
  value(2)=value(2)+(c(4)*(y*((u470)*u2857+u3774*u2856+u3997)))
  u4011=-2.25882569787428772d-1*u3797
  value(3)=value(3)+(c(4)*(u528*(u3876*u2857+u3998+u3873+u4011)))
  u4011=-1.59723096848534447d-1*u3797
  u977=(u2857*(u4011+u4037+u3992)+u3768*u2856+u3999)
  value(4)=value(4)+(c(4)*(y*u977))
  u3774=-6.73451707969374607d-2*u3797
  value(5)=value(5)+(c(4)*(z*(u2857*(u3774+u3830)+u2856*(u990+u3878)+u2&
 &855*(u990+u3994))))
  value(6)=value(6)+(c(4)*(x*u977))
  value(7)=value(7)+(c(4)*(z*((u3754)*u2857+u2856*(u4002+u4000)+u2855*(&
 &u3825+u3995))))
  u3774=-3.5215632741192216d-2*u3797
  value(8)=value(8)+(c(4)*(x*((u3874+u4001)*u2857+u741*u2856+u3774*u285&
 &5)))
  value(9)=value(9)+(c(4)*((u609*u2867+u2855*(u3875+u3778*u2856))*z))
  if ( lmax .eq. 1 ) return
  u3774=A7_(p,pd,erfpd,exppd2)
  u780=p**7
  u599=u2858*u3774*u780
  u4011=9.85691283332452413d1*u599
  u3840=2.96880505764235461d3*u3774
  u3995=-5.67185232289765129d1*u3621
  u677=u2858*(u3840+u3995)*u780
  u3778=1.49407276290432138d-1*u677
  u950=-3.48616978011008323d-1*u677
  u513=(1.1d1+u926)
  u609=u2858*(3.26568556340659007d4*u3774+u3995*u513)*u780
  u977=-4.98024254301440461d-2*u609
  u737=4.98024254301440461d-2*u609
  u3876=u977*u2855
  u3877=u737*u2855
  value(1)=value(1)+(c(5)*(u2983*((u3778+u3876)*u2856+u2855*(u950+u3877&
 &)+u4011)))
  u950=6.96988990600847647d1*u599
  u504=3.5215632741192216d-2*u677
  u638=-5.28234491117883239d-1*u677
  u969=-3.5215632741192216d-2*u609
  u3851=1.05646898223576648d-1*u609
  u3825=u3851*u2855
  value(2)=value(2)+(c(5)*(y*((u504+u969*u2855)*u2856+u2855*(u638+u3825&
 &)+u950)*z))
  u638=-3.72556286454539258d1*u599
  u3869=+u638
  u4013=-3.38823854681143158d-1*u677
  u668=5.6470642446857193d-2*u677
  u3768=1.31764832376000117d-1*u677
  u990=1.12941284893714386d-1*u609
  u741=-1.8823547482285731d-2*u609
  u3786=u741*u2855
  u4002=u990*u2855
  value(3)=value(3)+(c(5)*(u2983*((u4013+u4002)*u2857+(u668+u3786)*u285&
 &6+u2855*(u3768+u3786)+u3869)))
  u3797=-2.6343707652568261d1*u599
  u4017=-5.32410322828448158d-2*u677
  u4036=3.99307742121336118d-2*u677
  u4035=1.99653871060668059d-1*u677
  u712=5.32410322828448158d-2*u609
  u3830=-3.99307742121336118d-2*u609
  u3754=u3830*u2855
  u3878=u712*u2855
  u4055=(u2984*((u4017+u3878)*u2857+(u4036+u3754)*u2856+u2855*(u4035+u3&
 &754)+u3797))
  value(4)=value(4)+(c(5)*u4055)
  u470=-1.66612236391446779d1*u599
  u4037=+u470
  u621=-1.68362926992343652d-2*u677
  u3992=-1.13437046457953026d2*u3621
  u4081=u2858*(6.92721180116549409d3*u3774+u3992)*u780
  u3831=2.52544390488515477d-2*u4081
  u605=-5.68224878599159824d-2*u677
  u3874=2.52544390488515477d-2*u2858
  u873=u3874*(7.58694625841935067d3*u3774+u3992)*u780
  u585=2.02035512390812382d-1*u677
  u3776=1.68362926992343652d-2*u609
  u3875=-1.26272195244257739d-2*u2858
  u3994=8.0d0*pd2
  u675=+u3875*(1.69221888285614213d5*u3774+u3995*(5.7d1+u3994))*u780
  u5=5.68224878599159824d-2*u609
  u3832=-1.07331365957619078d-1*u677
  u4112=-5.05088780977030955d-2*u609
  u755=6.31360976221288694d-2*u609
  u4079=6.31360976221288694d-3*u609
  u3826=u4079*u2855
  u3879=u4112*u2855
  u4003=u755*u2855
  u4004=u5*u2855
  value(5)=value(5)+(c(5)*(u2857*(u4037+u2855*(u3879+u585)+(u3776*u2855&
 &+u621)*u2857)+u2856*(u3831+u2855*(u4003+u675)+(u4004+u605)*u2856)+u28&
 &55*(u873+u2855*(u3826+u3832))+u4037))
  u3795=-7.9031122957704783d1*u599
  u4014=-1.59723096848534447d-1*u677
  u539=1.19792322636400836d-1*u677
  u3939=2.79515419484935283d-1*u677
  value(6)=value(6)+(c(5)*(u2985*((u4014+u3878)*u2857+(u539+u3754)*u285&
 &6+u2855*(u3939+u3754)+u3795)))
  u595=-1.86278143227269629d1*u599
  u4031=+u595
  u3899=5.58834429681808886d1*u599
  u919=u2858*(1.64933614313464145d3*u3774+u2966)*u780
  u575=1.12941284893714386d-1*u919
  u953=-6.58824161880000585d-2*u677
  u4000=-7.08981540362206411d0*u3621
  u3977=u2858*(3.2986722862692829d2*u3774+u4000)*u780
  u761=4.51765139574857544d-1*u3977
  u4001=(1.3d1+u926)
  u43=u3621*u4001
  u827=u2858*(3.85944657493506099d4*u3774-5.67185232289765129d1*u43)*u7&
 &80
  u3828=-5.6470642446857193d-2*u827
  u3829=6.58824161880000585d-2*u609
  u3763=2.82353212234285965d-2*u677
  u893=5.6470642446857193d-2*u609
  u4016=-9.41177374114286549d-3*u609
  u3787=u4016*u2855
  u3827=u893*u2855
  value(7)=value(7)+(c(5)*((u3899+u2855*(u3827+u4013))*u2857+u2856*(u57&
 &5+u2855*(u3827+u3828)+(u3829*u2855+u953)*u2856)+u2855*(u761+u2855*(u3&
 &787+u3763))+u4031))
  u4031=3.16940694670729944d-1*u677
  u3763=-2.46509429188345512d-1*u677
  u953=-1.05646898223576648d-1*u609
  u761=3.5215632741192216d-2*u609
  u3828=u953*u2855
  u3829=u761*u2855
  value(8)=value(8)+(c(5)*(x*((u4031+u3828)*u2856+u2855*(u3763+u3829)+u&
 &950)*z))
  u3763=-4.92845641666226206d1*u599
  u575=+u3763
  u608=-1.24506063575360115d-2*u677
  u444=-u3763
  u3763=3.73518190726080346d-1*u677
  u3850=1.24506063575360115d-2*u609
  u866=-1.12055457217824104d-1*u677
  u3887=-7.47036381452160691d-2*u609
  u8=u3850*u2855
  u3880=u3887*u2855
  value(9)=value(9)+(c(5)*(u2856*(u575+u2855*(u3880+u3763)+(u8+u608)*u2&
 &856)+u2855*(u444+u2855*(u8+u866))))
  u389=4.98024254301440461d-2*u677
  u663=-4.98024254301440461d-2*u677
  u3881=u663*u2855
  value(1)=value(1)+(c(6)*(u2856*(u575+u737*u999+(u3876+u389)*u2856)+u2&
 &855*(u444+u3881)))
  u3833=-1.05646898223576648d-1*u677
  u4005=u969*u2856
  u469=u3825+u4005
  u4006=u3833*u2855
  u628=u4006+u950
  value(2)=value(2)+(c(6)*(x*(u2856*(u3833+u469)+u628)*z))
  u674=-u638
  u638=-5.6470642446857193d-2*u4081
  u3873=2.25882569787428772d-1*u2858
  u4081=(7.0d0+pd2)
  u3999=u3621*u4081
  u4008=u3873*(2.07816354034964823d4*u3774-5.67185232289765129d1*u3999)&
 &*u780
  u654=-1.31764832376000117d-1*u609
  u4007=u654*u2855
  value(3)=value(3)+(c(6)*(u2856*(u638+u2855*(u4007+u4008)+(u4007+u3768&
 &)*u2856)+u2855*(u638+u3768*u2855)+u674))
  u4008=u3830*u2856
  u654=u3754+u4008
  u638=u4036*u2855+u3797
  u3882=u712*u2856
  u88=(u2985*((u4017+u3882)*u2857+u2856*(u4035+u654)+u638))
  value(4)=value(4)+(c(6)*u88)
  u844=-u470
  u470=-5.05088780977030955d-2*u677
  u520=1.26272195244257739d-2*u609
  u3883=u4112*u2856
  u4009=u3776*u2857
  u620=u3879+u3883+u4009
  u3884=u4079*u2856
  u683=u520*u2855+u3884
  value(5)=value(5)+(c(6)*(u2983*(u2857*(u585+u620)+u2856*(u470+u683)+u&
 &2855*(u470+u3826)+u844)))
  value(6)=value(6)+(c(6)*u4055)
  u844=-3.7647094964571462d-2*u677
  u470=-5.6470642446857193d-2*u609
  u4055=9.41177374114286549d-3*u609
  u3846=3.7647094964571462d-2*u677
  u3885=u470*u2856
  u997=u3827+u3885
  u3788=u4055*u2856
  value(7)=value(7)+(c(6)*(u2983*((u997)*u2857+u2856*(u844+u3788)+u2855&
 &*(u3846+u3787))))
  u844=-u950
  u3846=1.05646898223576648d-1*u677
  value(8)=value(8)+(c(6)*(y*((u3846+u3828)*u2856+u2855*(u3846+u3829)+u&
 &844)*z))
  u670=-u4011
  u4011=9.96048508602880921d-2*u677
  u3886=u3850*u2856
  u739=u3880+u3886
  value(9)=value(9)+(c(6)*(u2983*(u2856*(u4011+u739)+u2855*(u4011+u8)+u&
 &670)))
  u4011=-1.49407276290432138d-1*u677
  value(1)=value(1)+(c(7)*(y*((u389+u3876)*u2856+u2855*(u4011+u3877))*z&
 &))
  u4012=-2.11293796447153296d-1*u677
  value(2)=value(2)+(c(7)*(u2983*((u4012+u469)*u2857+u504*u2856+u628)))
  u504=7.45112572909078515d1*u599
  u950=-1.12941284893714386d-1*u677
  u628=1.8823547482285731d-2*u677
  u3755=-1.69411927340571579d-1*u677
  value(3)=value(3)+(c(7)*(u2984*((u950+u4002)*u2857+(u628+u3786)*u2856&
 &+u2855*(u3755+u3786)+u504)))
  u3918=-7.98615484242672237d-2*u677
  u414=u654+u712*u2857
  u959=(u2983*(u2857*(u3918+u414)+u4036*u2856+u638))
  value(4)=value(4)+(c(7)*u959)
  u638=-6.66448945565787114d1*u599
  u598=3.36725853984687303d-2*u677
  u786=7.57633171465546432d-2*u677
  u3969=(u2857*(u598+u620)+u2856*(u786+u683)+u2855*(u786+u3826)+u638)
  value(5)=value(5)+(c(7)*(u2985*u3969))
  u620=-1.31718538262841305d1*u599
  u683=6.58592691314206525d1*u599
  u988=-3.99307742121336118d-2*u677
  u3876=3.99307742121336118d-2*u2858
  u3901=u3876*(2.30907060038849803d3*u3774+u3995)*u780
  u610=(6.0d0+pd2)
  u4002=u3621*u610
  u771=u2858*(8.90641517292706383d3*u3774-2.83592616144882564d1*u4002)*&
 &u780
  u4015=-1.59723096848534447d-1*u771
  u3928=3.99307742121336118d-2*u609
  u4010=u3928*u2855
  u3898=u4010+u4015
  value(6)=value(6)+(c(7)*(u2857*(u683+u2855*(u3754+u3918)+(u3878+u4017&
 &)*u2857)+u2856*(u4036+u2855*(u3898)+(u4010+u988)*u2856)+u3901*u2855+u&
 &620))
  u966=1.12941284893714386d-1*u677
  u4224=-7.52941899291429239d-2*u677
  value(7)=value(7)+(c(7)*(u2985*((u950+u997)*u2857+u2856*(u966+u3788)+&
 &u2855*(u4224+u3787)+u504)))
  u4224=-3.48494495300423824d1*u599
  u3834=+u4224
  u3862=-u4224
  u4224=u2858*(3.62853951489621119d3*u3774+u3995)*u780
  u3774=1.05646898223576648d-1*u4224
  u975=-4.22587592894306592d-1*u771
  u4007=-1.40862530964768864d-1*u677
  value(8)=value(8)+(c(7)*((u3862+u2855*(u3829+u4012))*u2857+u2856*(u38&
 &46+u2855*(u3825+u975)+(u3825+u3833)*u2856)+u2855*(u3774+u4007*u2855)+&
 &u3834))
  value(9)=value(9)+(c(7)*(x*(u2856*(u3778+u739)+u2855*(u663+u8))*z))
  u3774=3.48616978011008323d-1*u677
  u975=u3877+u977*u2856
  value(1)=value(1)+(c(8)*(u2983*(u2856*(u3774+u975)+u4011*u2855+u670))&
 &)
  u670=2.46509429188345512d-1*u677
  u3774=-3.16940694670729944d-1*u677
  u4011=u3774*u2855
  value(2)=value(2)+(c(8)*(y*(u2856*(u670+u469)+u4011+u844)*z))
  u670=u3786+u741*u2856
  u4012=u990*u2856
  value(3)=value(3)+(c(8)*(u2983*((u4013+u4012)*u2857+u2856*(u3768+u670&
 &)+u668*u2855+u3869)))
  u4013=u539*u2855
  value(4)=value(4)+(c(8)*(u2984*((u4014+u3882)*u2857+u2856*(u3939+u654&
 &)+u4013+u3795)))
  value(5)=value(5)+(c(8)*(u2857*(u4037+u2856*(u3883+u585)+(u3776*u2856&
 &+u621)*u2857)+u2856*(u873+u2855*(u4004+u675)+u2856*(u3884+u4003+u3832&
 &))+u2855*(u3831+u605*u2855)+u4037))
  value(6)=value(6)+(c(8)*u88)
  u585=-u595
  u5=-u3899
  u3831=-4.51765139574857544d-1*u3977
  u873=3.38823854681143158d-1*u677
  u605=-2.82353212234285965d-2*u677
  u4035=-1.12941284893714386d-1*u919
  u4014=5.6470642446857193d-2*u827
  u668=6.58824161880000585d-2*u677
  u3939=-6.58824161880000585d-2*u609
  value(7)=value(7)+(c(8)*((u5+u2856*(u3885+u873))*u2857+u2856*(u3831+u&
 &2855*(u3939*u2855+u4014)+u2856*(u3788+u470*u2855+u605))+u2855*(u4035+&
 &u668*u2855)+u585))
  u585=5.28234491117883239d-1*u677
  u5=-3.5215632741192216d-2*u677
  u605=u3829+u953*u2856
  u4035=u5*u2855+u844
  value(8)=value(8)+(c(8)*(x*(u2856*(u585+u605)+u4035)*z))
  value(9)=value(9)+(c(8)*(u2856*(u444+u2855*(u8+u3763)+u2856*(u3886+u3&
 &880+u866))+u2855*(u575+u608*u2855)))
  value(1)=value(1)+(c(9)*(x*(u2856*(u3778+u975)+u3881)*z))
  u585=-1.05646898223576648d-1*u4224
  u608=2.11293796447153296d-1*u677
  u4014=1.40862530964768864d-1*u677
  u866=4.22587592894306592d-1*u771
  value(2)=value(2)+(c(9)*((u3834+u2856*(u4005+u608))*u2857+u2856*(u585&
 &+u2855*(u3828+u866)+(u3828+u4014)*u2856)+u2855*(u3833+u3846*u2855)+u3&
 &862))
  u585=u628*u2855+u504
  value(3)=value(3)+(c(9)*(u2985*((u950+u4012)*u2857+u2856*(u3755+u670)&
 &+u585)))
  value(4)=value(4)+(c(9)*(u2857*(u683+u2856*(u4008+u3918)+(u3882+u4017&
 &)*u2857)+u2856*(u3901+u2855*(u3928*u2856+u3898))+u2855*(u4036+u988*u2&
 &855)+u620))
  value(5)=value(5)+(c(9)*(u2984*u3969))
  value(6)=value(6)+(c(9)*u959)
  u4014=-u504
  u866=7.52941899291429239d-2*u677
  value(7)=value(7)+(c(9)*(u2984*((u966+u997)*u2857+u2856*(u866+u3788)+&
 &u2855*(u950+u3787)+u4014)))
  u4014=u3846*u2856
  value(8)=value(8)+(c(9)*(u2983*((u608+u605)*u2857+u4014+u4035)))
  value(9)=value(9)+(c(9)*(y*(u2856*(u663+u739)+u2855*(u3778+u8))*z))
  value(1)=value(1)+(c(10)*(u2983*((u975)*u2857+u389*u2856+u3881)))
  value(2)=value(2)+(c(10)*(u2984*((u469)*u2857+u4014+u4011)))
  u866=-5.6470642446857193d-1*u677
  value(3)=value(3)+(c(10)*(u2983*(u2857*(u866+u670+u990*u2857)+u628*u2&
 &856+u585)))
  u866=1.05374830610273044d2*u599
  u3755=-3.7268722597991371d-1*u677
  u3774=(u2857*(u3755+u414)+u539*u2856+u4013+u866)
  value(4)=value(4)+(c(10)*(u2984*u3774))
  u4011=4.16530590978616946d0*u599
  u977=6.24795886467925419d1*u599
  u969=-1.51526634293109286d-1*u677
  u4013=-2.52544390488515477d-2*u2858
  u663=+u4013*(u3840+u2966)*u780
  u950=2.65171610012941251d-1*u677
  u3918=6.31360976221288694d-3*u677
  u608=5.05088780977030955d-2*u771
  u520=-1.26272195244257739d-2*u609
  u4015=u520*u2855
  value(5)=value(5)+(c(10)*(u2857*(u977+u2855*(u3826+u950)+u2856*(u3884&
 &+u950)+u2857*(u4009+u3883+u3879+u969))+u2856*(u663+u2855*(u4015+u608)&
 &+(u4015+u3918)*u2856)+u2855*(u663+u3918*u2855)+u4011))
  value(6)=value(6)+(c(10)*(u2985*u3774))
  u712=2.82353212234285965d-1*u677
  u977=-9.41177374114286549d-3*u677
  u969=-2.82353212234285965d-1*u677
  u3776=9.41177374114286549d-3*u677
  value(7)=value(7)+(c(10)*(u2857*(u2855*(u3787+u969)+u2856*(u3788+u712&
 &)+(u3885+u3827)*u2857)+u2856*(u3869+u977*u2856)+u2855*(u674+u3776*u28&
 &55)))
  value(8)=value(8)+(c(10)*(u2985*((u605)*u2857+u4031*u2856+u4006)))
  u3776=-2.46422820833113103d1*u599
  u977=+u3776
  u969=-u3776
  u3776=1.49407276290432138d-1*u919
  u712=-7.47036381452160691d-2*u677
  u3827=-8.71542445027520806d-2*u677
  u953=-2.98814552580864276d-1*u771
  u520=7.47036381452160691d-2*u609
  u4016=u520*u2855
  value(9)=value(9)+(c(10)*((u969+u2855*(u8+u712)+u2856*(u3886+u712))*u&
 &2857+u2856*(u3776+u2855*(u4016+u953)+(u4016+u3827)*u2856)+u2855*(u377&
 &6+u3827*u2855)+u977))
  if ( lmax .eq. 2 ) return
  u977=A8_(p,pd,erfpd,exppd2)
  u969=p**8
  u953=u2858*u977*u969
  u520=-9.85691283332452413d1*u953
  u712=+u520
  u3850=2.96880505764235461d3*u977
  u3827=-5.67185232289765129d1*u382
  u3776=u2858*(u3850+u3827)*u969
  u893=-1.49407276290432138d-1*u3776
  u761=1.34466548661388924d0*u3776
  u4055=3.26568556340659007d4*u977
  u470=u2858*(u4055+u3827*u513)*u969
  u609=2.98814552580864276d-1*u470
  u771=-5.97629105161728553d-1*u470
  u4011=pd2**2
  u4006=4.0d0*u4011
  u663=(u926+1.1d1)
  u950=u2858*(4.24539123242856709d5*u977+u3827*(1.3d1*u663+u4006))*u969
  u3918=-4.98024254301440461d-2*u950
  u608=4.98024254301440461d-2*u950
  u866=u608*u2855
  u3755=u3918*u2855
  value(1)=value(1)+(c(11)*(y*((u893+u2855*(u3755+u609))*u2856+u2855*(u&
 &761+u2855*(u866+u771))+u712)))
  u771=1.26776277868291978d0*u3776
  u988=1.05646898223576648d-1*u470
  u3774=-9.50822084012189831d-1*u470
  u5=-3.5215632741192216d-2*u950
  u605=1.05646898223576648d-1*u950
  u3851=u605*u2855
  u3830=u5*u2855
  value(2)=value(2)+(c(11)*(u2983*((u988+u3830)*u2856+u2855*(u3774+u385&
 &1)+u771)*z))
  u3774=3.72556286454539258d1*u953
  u4112=3.38823854681143158d-1*u3776
  u4036=-5.6470642446857193d-2*u3776
  u539=-5.08235782021714737d-1*u3776
  u3846=-6.77647709362286316d-1*u470
  u3778=1.12941284893714386d-1*u470
  u4035=2.25882569787428772d-1*u470
  u628=1.12941284893714386d-1*u950
  u786=-1.8823547482285731d-2*u950
  u585=u786*u2855
  u389=u585+u4035
  u668=u2855*(u389)
  u3831=u628*u2855
  value(3)=value(3)+(c(11)*(y*((u4112+u2855*(u3831+u3846))*u2857+(u4036&
 &+u2855*(u585+u3778))*u2856+u2855*(u539+u668)+u3774)))
  u873=-4.79169290545603342d-1*u3776
  u3939=-1.59723096848534447d-1*u470
  u4007=1.19792322636400836d-1*u470
  u3763=3.59376967909202506d-1*u470
  u598=5.32410322828448158d-2*u950
  u3901=-3.99307742121336118d-2*u950
  u966=u3901*u2855
  u3832=u598*u2855
  u737=(u528*((u3939+u3832)*u2857+(u4007+u966)*u2856+u2855*(u3763+u966)&
 &+u873))
  value(4)=value(4)+(c(11)*u737)
  u990=-3.03053268586218573d-1*u2858
  u3928=+u990*(1.64933614313464145d3*u977+u500)*u969
  u677=-4.54579902879327859d-1*u3776
  u919=-5.05088780977030955d-2*u470
  u3977=u382*u4001
  u4224=u2858*(3.85944657493506099d4*u977-5.67185232289765129d1*u3977)*&
 &u969
  u3797=1.51526634293109286d-1*u4224
  u4037=-1.70467463579747947d-1*u470
  u844=1.48440252882117731d4*u977
  u3869=(1.0d1+pd2)
  u575=u382*u3869
  u3795=u2858*(u844-2.83592616144882564d1*u575)*u969
  u670=2.02035512390812382d-1*u3795
  u638=4.04071024781624764d-1*u470
  u620=1.68362926992343652d-2*u950
  u3786=1.6d1*u4011
  u444=+u3875*(2.51457788382307436d6*u977+u3827*(7.7d1*u663+u3786))*u96&
 &9
  u674=5.68224878599159824d-2*u950
  u683=-1.452130245308964d-1*u470
  u3862=-5.05088780977030955d-2*u950
  u675=6.31360976221288694d-2*u950
  u599=6.31360976221288694d-3*u950
  u3840=u599*u2855
  u3833=u3862*u2855
  u3834=u674*u2855
  u3887=u675*u2855
  u4017=u620*u2855
  value(5)=value(5)+(c(11)*(x*(u2857*(u677+u2855*(u3833+u638)+(u4017+u9&
 &19)*u2857)+u2856*(u3797+u2855*(u3887+u444)+(u3834+u4037)*u2856)+u2855&
 &*(u670+u2855*(u3840+u683))+u3928)))
  u780=7.9031122957704783d1*u953
  u513=1.59723096848534447d-1*u3776
  u4079=-1.19792322636400836d-1*u3776
  u595=-1.07813090372760752d0*u3776
  u3899=-3.19446193697068895d-1*u470
  u504=2.39584645272801671d-1*u470
  u741=4.79169290545603342d-1*u470
  u3969=u2855*(u966+u741)
  value(6)=value(6)+(c(11)*(z*((u513+u2855*(u3832+u3899))*u2857+(u4079+&
 &u2855*(u966+u504))*u2856+u2855*(u595+u3969)+u780)))
  u88=u2858*(2.52898208613978356d3*u977+u3827)*u969
  u959=-1.69411927340571579d-1*u88
  u827=8.47059636702857895d-1*u3776
  u414=8.90641517292706383d3*u977
  u469=u2858*(u414+u500*u610)*u969
  u610=6.77647709362286316d-1*u469
  u654=-1.97647248564000175d-1*u470
  u997=u2858*(u414+u3827*(3.0d0+pd2))*u969
  u414=1.12941284893714386d-1*u997
  u739=-5.6470642446857193d-1*u470
  u975=u2858*(5.55166545779120312d5*u977+u3827*(1.7d1*u663+u4006))*u969
  u3998=-5.6470642446857193d-2*u975
  u755=6.58824161880000585d-2*u950
  u3898=8.47059636702857895d-2*u470
  u621=5.6470642446857193d-2*u950
  u3768=-9.41177374114286549d-3*u950
  u8=u3768*u2855
  u4031=u621*u2855
  u4018=u755*u2855
  value(7)=value(7)+(c(11)*(x*((u827+u2855*(u4031+u739))*u2857+u2856*(u&
 &610+u2855*(u4031+u3998)+(u4018+u654)*u2856)+u2855*(u414+u2855*(u8+u38&
 &98))+u959)))
  u959=-6.96988990600847647d1*u953
  u827=+u959
  u610=-3.16940694670729944d-1*u3776
  u414=9.50822084012189831d-1*u3776
  u3998=6.33881389341459887d-1*u470
  u3916=-4.22587592894306592d-1*u470
  u973=-1.05646898223576648d-1*u950
  u612=3.5215632741192216d-2*u950
  u398=u612*u2855
  u3756=u973*u2855
  value(8)=value(8)+(c(11)*(((u610+u2855*(u3756+u3998))*u2856+u2855*(u4&
 &14+u2855*(u398+u3916))+u827)*z))
  u414=-8.9644365774259283d-1*u3776
  u456=-3.73518190726080346d-2*u470
  u94=5.97629105161728553d-1*u3776
  u462=6.72332743306944622d-1*u470
  u594=1.24506063575360115d-2*u950
  u4122=-1.86759095363040173d-1*u470
  u573=-7.47036381452160691d-2*u950
  u782=u594*u2855
  u3789=u573*u2855
  value(9)=value(9)+(c(11)*(x*(u2856*(u414+u2855*(u3789+u462)+(u782+u45&
 &6)*u2856)+u2855*(u94+u2855*(u782+u4122))+u712)))
  u4077=1.49407276290432138d-1*u470
  u70=3.48616978011008323d-1*u3776
  u664=-1.99209701720576184d-1*u470
  u877=-4.98024254301440461d-2*u470
  u4019=u877*u2855
  value(1)=value(1)+(c(12)*(x*(u2856*(u893+u2855*(u866+u664)+(u3755+u40&
 &77)*u2856)+u2855*(u70+u4019)+u712)))
  u70=1.05646898223576648d-1*u3776
  u4078=3.5215632741192216d-2*u470
  u3791=5.28234491117883239d-1*u3776
  u80=-1.05646898223576648d-1*u470
  u3888=u80*u2855
  value(2)=value(2)+(c(12)*((u2856*(u70+u2855*(u3851+u3916)+(u3830+u407&
 &8)*u2856)+u2855*(u3791+u3888)+u827)*z))
  u827=-1.41796308072441282d1*u382
  u4025=u2858*(7.69690200129499343d2*u977+u827)*u969
  u3757=1.35529541872457263d0*u4025
  u3793=2.28597989438461305d5*u977
  u3886=-5.6470642446857193d-2*u2858
  u4015=1.2d1*pd2
  u3787=(7.7d1+u4015)
  u3883=u382*u3787
  u657=+u3886*(u3793-5.67185232289765129d1*u3883)*u969
  u4023=3.95294497128000351d-1*u470
  u3826=-1.8823547482285731d-2*u2858
  u998=u382*(9.1d1+u4015)
  u3792=+u3826*(2.7016126024545427d5*u977-5.67185232289765129d1*u998)*u&
 &969
  u4016=3.0d0*u4011
  u756=3.01176759716571696d-1*u2858*(u3793+u500*(1.4d1*u663+u4016))*u96&
 &9
  u3793=-1.31764832376000117d-1*u950
  u847=1.31764832376000117d-1*u470
  u3889=u3793*u2855
  u684=u2855*(u3889+u756)
  value(3)=value(3)+(c(12)*(x*(u2856*(u657+u684+(u3889+u4023)*u2856)+u2&
 &855*(u3792+u847*u2855)+u3757)))
  u454=7.98615484242672237d-2*u2858
  u3788=-1.13437046457953026d2*u382
  u524=u454*(6.26747734391163751d3*u977+u3788)*u969
  u858=-5.32410322828448158d-2*u4224
  u3836=5.32410322828448158d-2*u470
  u55=-3.59376967909202506d-1*u3776
  u863=3.99307742121336118d-2*u470
  u3885=-3.99307742121336118d-2*u2858
  u476=+u3885*(1.57346668055044794d5*u977+u3827*(5.3d1+u3994))*u969
  u642=u2858*(u4055+u827*(4.0d0*u663+u4011))*u969
  u949=8.51856516525517053d-1*u642
  u976=-5.32410322828448158d-2*u950
  u3758=1.99653871060668059d-1*u470
  u4020=u976*u2855
  u909=(z*(u2857*(u858+u2855*(u4020+u949)+(u4020+u3836)*u2857)+u2856*(u&
 &55+u3969+(u966+u863)*u2856)+u2855*(u476+u3758*u2855)+u524))
  value(4)=value(4)+(c(12)*u909)
  u3969=-6.06106537172437146d-1*u4025
  u4025=-5.05088780977030955d-2*u3776
  u429=-1.68362926992343652d-2*u470
  u4014=(1.5d1+u926)
  u4197=u382*u4014
  u758=u2858*(4.45320758646353191d4*u977-5.67185232289765129d1*u4197)*u&
 &969
  u3893=5.05088780977030955d-2*u758
  u3839=-5.68224878599159824d-2*u470
  u4020=3.0d0*pd2
  u3765=(2.0d1+u4020)
  u3986=u382*u3765
  u3974=u2858*(u844-1.41796308072441282d1*u3986)*u969
  u844=4.04071024781624764d-1*u3974
  u4024=2.02035512390812382d-1*u470
  u445=+u3875*(2.44926417255494255d6*u977+u3827*(7.5d1*u663+u3786))*u96&
 &9
  u446=-2.33603561201876817d-1*u470
  value(5)=value(5)+(c(12)*(y*(u2857*(u4025+u2855*(u3833+u4024)+(u4017+&
 &u429)*u2857)+u2856*(u3893+u2855*(u3887+u445)+(u3834+u3839)*u2856)+u28&
 &55*(u844+u2855*(u3840+u446))+u3969)))
  value(6)=value(6)+(c(12)*u737)
  u737=-1.69411927340571579d-1*u3776
  u4033=1.69411927340571579d-1*u3776
  u3978=u2858*(u4055+u500*(2.2d1+u4020))*u969
  u4055=7.52941899291429239d-2*u3978
  u3794=-6.58824161880000585d-2*u470
  u760=5.04696859799200284d4*u977
  u606=u382*(1.7d1+u4020)
  u759=u2858*(u760-5.67185232289765129d1*u606)*u969
  u3837=1.12941284893714386d-1*u759
  u517=-3.38823854681143158d-1*u470
  u511=1.73081334860549274d6*u977
  u4009=1.2d1*u4011
  u3910=5.3d1*u663
  u448=u2858*(u511+u3827*(u3910+u4009))*u969
  u4034=-1.8823547482285731d-2*u448
  u4=-8.47059636702857895d-2*u470
  value(7)=value(7)+(c(12)*(y*((u4033+u2855*(u4031+u517))*u2857+u2856*(&
 &u4055+u2855*(u4031+u4034)+(u4018+u3794)*u2856)+u2855*(u3837+u2855*(u8&
 &+u4))+u737)))
  u4034=-4.22587592894306592d-1*u3776
  u3794=3.16940694670729944d-1*u470
  u3837=-3.5215632741192216d-2*u470
  value(8)=value(8)+(c(12)*(u2983*((u3794+u3756)*u2856+u2855*(u3837+u39&
 &8)+u4034)*z))
  u4034=-u520
  u4055=-9.96048508602880921d-2*u3776
  u520=-1.24506063575360115d-2*u470
  u4219=-5.97629105161728553d-1*u3776
  u685=3.237157652959363d-1*u470
  u650=3.73518190726080346d-2*u470
  u4133=(u782+u520)*u2856
  value(9)=value(9)+(c(12)*(y*(u2856*(u4055+u2855*(u3789+u685)+u4133)+u&
 &2855*(u4219+u2855*(u782+u650))+u4034)))
  u781=2.98814552580864276d-1*u3776
  u592=-3.48616978011008323d-1*u470
  value(1)=value(1)+(c(13)*(u2983*((u4077+u3755)*u2856+u2855*(u592+u866&
 &)+u781)*z))
  u781=u2858*(3.62853951489621119d3*u977+u3827)*u969
  u592=-1.05646898223576648d-1*u781
  u4127=3.16940694670729944d-1*u3776
  u597=3.5215632741192216d-2*u4224
  u4070=(1.7d1+u926)
  u3942=u382*u4070
  u971=u2858*(u760-5.67185232289765129d1*u3942)*u969
  u760=1.05646898223576648d-1*u971
  u3915=-6.33881389341459887d-1*u470
  u4083=-5.63450123859075456d-1*u642
  u698=-2.11293796447153296d-1*u470
  u4021=u698*u2855
  value(2)=value(2)+(c(13)*(y*((u4127+u2855*(u3851+u3915))*u2857+u2856*&
 &(u597+u2855*(u398+u4083)+(u398+u3837)*u2856)+u2855*(u760+u4021)+u592)&
 &))
  u592=5.6470642446857193d-1*u3776
  u597=5.6470642446857193d-2*u470
  u760=-9.4117737411428655d-2*u470
  value(3)=value(3)+(c(13)*(u528*((u517+u3831)*u2857+(u597+u585)*u2856+&
 &u2855*(u760+u585)+u592)))
  u4083=u2858*(2.74889357189106908d3*u977+u3827)*u969
  u722=-1.19792322636400836d-1*u4083
  u3950=1.99653871060668059d-1*u3776
  u3848=-5.32410322828448158d-2*u470
  u565=3.99307742121336118d-2*u4224
  u4215=-3.99307742121336118d-2*u470
  u4017=6.0d0*pd2
  u84=u382*(3.1d1+u4017)
  u4145=u3876*(9.20329567869129929d4*u977-5.67185232289765129d1*u84)*u9&
 &69
  u3809=-7.98615484242672237d-2*u470
  u644=-6.38892387394137789d-1*u642
  u601=3.99307742121336118d-2*u950
  u3790=u601*u2855
  u788=u3790+u644
  u4064=u2855*(u788)
  u3943=(y*(u2857*(u3950+u2855*(u966+u3809)+(u3832+u3848)*u2857)+u2856*&
 &(u565+u4064+(u3790+u4215)*u2856)+u2855*(u4145+u3809*u2855)+u722))
  value(4)=value(4)+(c(13)*u3943)
  u3760=-1.51526634293109286d-1*u88
  u843=1.34690341593874921d-1*u3974
  u3974=-6.73451707969374607d-2*u470
  u697=7.57633171465546432d-2*u3776
  u3796=-6.31360976221288694d-3*u470
  u600=1.63284278170329504d5*u977
  u682=u3874*(u600+u3827*(5.5d1+u4015))*u969
  u3879=-1.34690341593874921d-1*u2858
  u603=+u3879*(u600+u500*(1.0d1*u663+u4016))*u969
  u600=6.73451707969374607d-2*u950
  u607=-1.13644975719831965d-1*u470
  u924=-1.07331365957619078d-1*u470
  u721=1.26272195244257739d-2*u950
  u3890=u721*u2855
  value(5)=value(5)+(c(13)*(z*(u2857*(u843+u603*u2855+(u600*u2855+u3974&
 &)*u2857)+u2856*(u697+u2855*(u3890+u607)+(u3840+u3796)*u2856)+u2855*(u&
 &682+u2855*(u3840+u924))+u3760)))
  u697=-1.19792322636400836d-1*u2858
  u3796=+u697*(2.30907060038849803d3*u977+u3827)*u969
  u682=3.59376967909202506d-1*u3776
  u603=-1.19792322636400836d-1*u470
  u600=2.07816354034964823d4*u977
  u924=(7.0d0+u926)
  u962=u382*u924
  u4038=u3876*(u600-5.67185232289765129d1*u962)*u969
  u633=7.98615484242672237d-2*u470
  value(6)=value(6)+(c(13)*(x*(u2857*(u682+u2855*(u966+u633)+(u3832+u39&
 &39)*u2857)+u2856*(u4007+u4064+(u3790+u603)*u2856)+u4038*u2855+u3796))&
 &)
  u512=u2858*(3.40862802914492566d3*u977+u3827)*u969
  u3976=-1.69411927340571579d-1*u512
  u461=(8.0d0+pd2)
  u854=u382*u461
  u72=u2858*(u3850-7.08981540362206411d0*u854)*u969
  u3850=9.03530279149715087d-1*u72
  u3759=-5.6470642446857193d-2*u470
  u876=5.6470642446857193d-2*u3776
  u4075=-9.41177374114286549d-3*u470
  u755=2.25882569787428772d-1*u3978
  u3832=2.0d0*u4011
  u4018=-1.12941284893714386d-1*u2858
  u3884=+u4018*(3.59225411974724908d5*u977+u3827*(1.1d1*u663+u3832))*u9&
 &69
  u629=9.41177374114286549d-3*u950
  value(7)=value(7)+(c(13)*(z*(u2857*(u3850+u2855*(u3831+u3884)+(u4031+&
 &u3759)*u2857)+u2856*(u876+u3759*u2855+(u629*u2855+u4075)*u2856)+u2855&
 &*(u755+u2855*(u8+u654))+u3976)))
  u3976=u2858*(3.18871654339364014d3*u977+u3827)*u969
  u3850=-3.16940694670729944d-1*u3976
  u876=-3.16940694670729944d-1*u470
  u755=(4.9d1+u4017)
  u3884=u382*u755
  u654=u2858*(1.45471447824475376d5*u977-5.67185232289765129d1*u3884)*u&
 &969
  u3759=3.5215632741192216d-2*u654
  u4028=-3.5215632741192216d-1*u470
  u4039=-1.69035037157722637d0*u642
  u757=+u4039
  u3779=-1.40862530964768864d-1*u470
  u676=u2855*(u3851+u757)
  value(8)=value(8)+(c(13)*(x*((u3791+u2855*(u398+u4028))*u2857+u2856*(&
 &u3794+u676+(u3851+u876)*u2856)+u2855*(u3759+u3779*u2855)+u3850)))
  u3850=1.49407276290432138d-1*u3776
  u3759=3.73518190726080346d-1*u470
  u3779=-1.12055457217824104d-1*u470
  value(9)=value(9)+(c(13)*((u2856*(u893+u2855*(u3789+u3759)+u4133)+u28&
 &55*(u3850+u2855*(u782+u3779)))*z))
  u3791=-3.48616978011008323d-1*u3776
  u4133=4.98024254301440461d-2*u470
  u3946=1.99209701720576184d-1*u470
  u535=-1.49407276290432138d-1*u470
  u4080=(u3755+u4133)*u2856
  u3891=u535*u2855
  value(1)=value(1)+(c(14)*(y*(u2856*(u3791+u2855*(u866+u3946)+u4080)+u&
 &2855*(u3850+u3891)+u4034)))
  u3791=4.22587592894306592d-1*u3776
  u3835=u5*u2856
  u596=u3851+u3835
  u4022=u876*u2855
  value(2)=value(2)+(c(14)*(u2983*(u2856*(u4078+u596)+u4022+u3791)*z))
  value(3)=value(3)+(c(14)*(y*(u2856*(u3792+u684+(u3889+u847)*u2856)+u2&
 &855*(u657+u4023*u2855)+u3757)))
  u3791=u3901*u2856
  u657=u966+u3791
  u3792=u598*u2856
  u4023=u4007*u2855
  u847=(u528*((u3939+u3792)*u2857+u2856*(u3763+u657)+u4023+u873))
  value(4)=value(4)+(c(14)*u847)
  u3757=u599*u2856
  u684=u3757+u3887
  u3793=u3862*u2856
  u3892=u620*u2856
  value(5)=value(5)+(c(14)*(x*(u2857*(u4025+u2856*(u3793+u4024)+(u3892+&
 &u429)*u2857)+u2856*(u844+u2855*(u3834+u445)+u2856*(u684+u446))+u2855*&
 &(u3893+u3839*u2855)+u3969)))
  value(6)=value(6)+(c(14)*u909)
  u949=-1.12941284893714386d-1*u759
  u858=3.38823854681143158d-1*u470
  u476=-5.6470642446857193d-2*u950
  u429=-7.52941899291429239d-2*u3978
  u3839=1.8823547482285731d-2*u448
  u448=6.58824161880000585d-2*u470
  u446=-6.58824161880000585d-2*u950
  u3836=u476*u2856
  u863=u2856*(u3836+u858)
  u3758=u629*u2856
  u3893=u476*u2855
  u844=u3758+u3893
  u4024=u446*u2855
  value(7)=value(7)+(c(14)*(x*((u737+u863)*u2857+u2856*(u949+u2855*(u40&
 &24+u3839)+u2856*(u844+u3898))+u2855*(u429+u448*u2855)+u4033)))
  u429=-u959
  u3839=-5.28234491117883239d-1*u3776
  u448=-1.05646898223576648d-1*u3776
  u3898=4.22587592894306592d-1*u470
  u55=(u3756+u988)
  value(8)=value(8)+(c(14)*((u2856*(u3839+u2855*(u398+u3898)+u55*u2856)&
 &+u2855*(u448+u3837*u2855)+u429)*z))
  u3837=u594*u2856
  u3969=u3837+u3789
  u4025=u520*u2855
  value(9)=value(9)+(c(14)*(x*(u2856*(u4219+u2855*(u782+u685)+u2856*(u3&
 &969+u650))+u2855*(u4055+u4025)+u4034)))
  value(1)=value(1)+(c(15)*((u2856*(u893+u608*u999+u4080)+u2855*(u3850+&
 &u4019))*z))
  u685=3.16940694670729944d-1*u4083
  u650=1.03908177017482411d5*u977
  u4055=(3.5d1+u4017)
  u4219=u382*u4055
  u524=u2858*(u650-5.67185232289765129d1*u4219)*u969
  u4033=-1.05646898223576648d-1*u524
  u959=2.11293796447153296d-1*u470
  u909=3.5215632741192216d-1*u470
  u4080=-1.05646898223576648d-1*u4224
  u976=-u4039
  u4039=u2855*(u3756+u976)
  u458=u2856*(u3835+u959)
  value(2)=value(2)+(c(15)*(x*((u448+u458)*u2857+u2856*(u4033+u4039+(u3&
 &756+u909)*u2856)+u2855*(u4080+u988*u2855)+u685)))
  u685=3.38823854681143158d-1*u4083
  u4033=-1.12941284893714386d-1*u4224
  u4080=1.8823547482285731d-2*u470
  u4019=4.0d0*pd2
  u947=u2858*(6.8282516325774156d4*u977+u3827*(2.3d1+u4019))*u969
  u3983=-1.69411927340571579d-1*u947
  u4027=1.80706055829943018d0*u642
  u447=-1.12941284893714386d-1*u950
  u3900=3.57647402163428889d-1*u470
  u4026=u447*u2855
  value(3)=value(3)+(c(15)*(z*(u2857*(u4033+u2855*(u4026+u4027)+(u4026+&
 &u3778)*u2857)+u2856*(u737+u668+(u585+u4080)*u2856)+u2855*(u3983+u3900&
 &*u2855)+u685)))
  u685=(x*(u2857*(u3950+u2856*(u3791+u3809)+(u3792+u3848)*u2857)+u2856*&
 &(u4145+u4064+(u3790+u3809)*u2856)+u2855*(u565+u4215*u2855)+u722))
  value(4)=value(4)+(c(15)*u685)
  u3900=-3.53562146683921668d-1*u3776
  u4080=1.34690341593874921d-1*u470
  u4033=5.05088780977030955d-2*u470
  u4027=u620*u2857
  value(5)=value(5)+(c(15)*(u528*(u2857*(u4080+u3833+u3793+u4027)+u2856&
 &*(u4033+u3890+u3757)+u2855*(u4033+u3840)+u3900)))
  value(6)=value(6)+(c(15)*u3943)
  u3900=7.52941899291429239d-2*u470
  u4033=-7.52941899291429239d-2*u470
  value(7)=value(7)+(c(15)*(u528*((u4031+u3836)*u2857+u2856*(u3900+u375&
 &8)+u2855*(u4033+u8))))
  u3900=-3.16940694670729944d-1*u4083
  u4033=1.05646898223576648d-1*u4224
  u3983=1.05646898223576648d-1*u524
  u737=(u3851+u80)
  value(8)=value(8)+(c(15)*(y*((u70+u2855*(u398+u698))*u2857+u2856*(u40&
 &33+u676+u737*u2856)+u2855*(u3983+u4028*u2855)+u3900)))
  u3900=-2.98814552580864276d-1*u3776
  u3983=9.96048508602880921d-2*u470
  value(9)=value(9)+(c(15)*(u2983*(u2856*(u3983+u3789+u3837)+u2855*(u39&
 &83+u782)+u3900)*z))
  u3983=4.98024254301440461d-2*u4224
  u4028=1.49407276290432138d-1*u4224
  u4033=-2.98814552580864276d-1*u470
  u4083=-7.96838806882304737d-1*u642
  value(1)=value(1)+(c(16)*(y*((u3850+u2855*(u866+u4033))*u2857+u2856*(&
 &u3983+u2855*(u866+u4083)+(u866+u877)*u2856)+u2855*(u4028+u664*u2855)+&
 &u893)))
  u4083=6.33881389341459887d-1*u3776
  u4028=u988*u2856
  value(2)=value(2)+(c(16)*(u528*((u698+u596)*u2857+u4028+u4022+u4083))&
 &)
  u4083=+u3886*(4.28827397215006777d3*u977+u3827)*u969
  u3983=6.21177066915429123d-1*u3776
  u877=-1.12941284893714386d-1*u470
  u664=1.8823547482285731d-2*u4224
  u3943=-1.8823547482285731d-2*u470
  u668=5.6470642446857193d-2*u758
  u758=-3.01176759716571696d-1*u642
  u676=1.8823547482285731d-2*u950
  u596=-3.7647094964571462d-2*u470
  u3894=u676*u2855
  u533=u2855*(u3894+u758)
  value(3)=value(3)+(c(16)*(y*(u2857*(u3983+u2855*(u585+u739)+(u3831+u8&
 &77)*u2857)+u2856*(u664+u533+(u3894+u3943)*u2856)+u2855*(u668+u596*u28&
 &55)+u4083)))
  u89=7.98615484242672237d-2*u3776
  u779=-2.92825677555646487d-1*u470
  u4029=u598*u2857
  u4030=u4007*u2856
  u4113=(u528*(u2857*(u779+u657+u4029)+u4030+u4023+u89))
  value(4)=value(4)+(c(16)*u4113)
  u3831=3.78816585732773216d-2*u2858
  u657=u3831*(4.72809694365263882d3*u977+u3827)*u969
  u386=-3.40934927159495895d-1*u3776
  u968=u382*u4081
  u4081=u2858*(u600-5.67185232289765129d1*u968)*u969
  u790=-7.57633171465546432d-2*u4081
  u593=2.65171610012941251d-1*u470
  u631=3.15680488110644347d-2*u470
  u4134=-5.05088780977030955d-2*u3795
  u3795=2.39917170964089704d-1*u470
  u3975=2.02035512390812382d-1*u642
  u3799=-1.26272195244257739d-2*u950
  u770=6.31360976221288694d-3*u470
  u3949=u4027+u3793
  u723=u2857*(u3949+u3833+u919)
  u3895=u3799*u2855
  u3902=u2855*(u3895+u3975)
  value(5)=value(5)+(c(16)*(x*(u2857*(u386+u2855*(u3840+u3795)+u2856*(u&
 &3757+u593)+u723)+u2856*(u790+u3902+(u3895+u631)*u2856)+u2855*(u4134+u&
 &770*u2855)+u657)))
  u548=+u697*(3.84845100064749672d3*u977+u3827)*u969
  u4027=2.66205161414224079d-2*u2858
  u778=(3.5d1+u4020)
  u4140=u382*u778
  u757=u4027*(u650-5.67185232289765129d1*u4140)*u969
  u650=-9.31718064949784276d-2*u470
  u703=3.99307742121336118d-2*u524
  u3833=6.0d0*u4011
  u4023=-2.66205161414224079d-2*u2858
  u524=+u4023*(1.14298994719230652d6*u977+u3827*(3.5d1*u663+u3833))*u96&
 &9
  u882=9.31718064949784276d-2*u950
  value(6)=value(6)+(c(16)*(z*(u2857*(u757+u524*u2855+(u882*u2855+u650)&
 &*u2857)+u703*u2855+u548)))
  u703=-7.45112572909078515d1*u953
  u524=+u703
  u882=-1.12941284893714386d-1*u3776
  u4063=2.82353212234285965d-1*u470
  u787=7.52941899291429239d-2*u3776
  u655=-2.44706117269714503d-1*u470
  u408=9.41177374114286549d-3*u470
  u527=u3836+u4031
  u4031=u4075*u2856
  u4032=u408*u2855
  value(7)=value(7)+(c(16)*(x*(u2857*(u592+u2855*(u8+u655)+u2856*(u3758&
 &+u4063)+(u527+u877)*u2857)+u2856*(u882+u4031)+u2855*(u787+u4032)+u524&
 &)))
  u655=2.11293796447153296d-1*u4081
  u524=3.16940694670729944d-1*u4224
  u882=-2.11293796447153296d-1*u2858*(2.93911700706593106d5*u977+u3827*&
 &(9.0d0*u663+u3832))*u969
  u4063=1.40862530964768864d-1*u950
  value(8)=value(8)+(c(16)*(z*(u2857*(u655+u2855*(u4063*u2855+u882)+u73&
 &7*u2857)+u2855*(u524+u3916*u2855)+u610)))
  u655=-2.24110914435648207d-1*u3776
  u524=2.24110914435648207d-1*u3776
  u882=1.49407276290432138d-1*u759
  u4063=-7.47036381452160691d-2*u470
  u787=-2.36561520793184219d-1*u470
  u3916=9.96048508602880921d-2*u3978
  u3978=-1.24506063575360115d-1*u470
  u759=-1.19525821032345711d0*u642
  u663=7.47036381452160691d-2*u950
  u737=-8.71542445027520806d-2*u470
  u3889=u2856*(u3837+u4063)
  u3838=u663*u2855
  u4062=u2855*(u3838+u759)
  value(9)=value(9)+(c(16)*(x*((u524+u2855*(u782+u3978)+u3889)*u2857+u2&
 &856*(u882+u4062+(u3838+u787)*u2856)+u2855*(u3916+u737*u2855)+u655)))
  u3897=-u761
  u761=5.97629105161728553d-1*u470
  u3896=u3918*u2856
  value(1)=value(1)+(c(17)*(x*(u2856*(u3897+u4033*u2855+u2856*(u3896+u8&
 &66+u761))+u3850*u2855+u4034)))
  u3897=-9.50822084012189831d-1*u3776
  u761=u3835+u3851
  u4033=u3915*u2855
  u4034=u4127*u2855
  value(2)=value(2)+(c(17)*((u2856*(u3897+u4033+u2856*(u761+u3898))+u40&
 &34+u429)*z))
  u3897=u628*u2856
  u3898=u786*u2856
  value(3)=value(3)+(c(17)*(x*((u4112+u2856*(u3897+u3846))*u2857+u2856*&
 &(u539+u3778*u2855+u2856*(u3898+u389))+u4036*u2855+u3774)))
  u539=u3791+u966
  u4035=u504*u2855
  u4036=u4079*u2855
  value(4)=value(4)+(c(17)*(z*((u513+u2856*(u3792+u3899))*u2857+u2856*(&
 &u595+u4035+u2856*(u539+u741))+u4036+u780)))
  value(5)=value(5)+(c(17)*(y*(u2857*(u677+u2856*(u3793+u638)+(u3892+u9&
 &19)*u2857)+u2856*(u670+u2855*(u3834+u444)+u2856*(u684+u683))+u2855*(u&
 &3797+u4037*u2855)+u3928)))
  value(6)=value(6)+(c(17)*u847)
  u4037=1.69411927340571579d-1*u88
  u873=-8.47059636702857895d-1*u3776
  u3899=-1.12941284893714386d-1*u997
  u997=5.6470642446857193d-1*u470
  u670=-6.77647709362286316d-1*u469
  u3928=5.6470642446857193d-2*u975
  u677=1.97647248564000175d-1*u470
  value(7)=value(7)+(c(17)*(y*((u873+u2856*(u3836+u997))*u2857+u2856*(u&
 &3899+u2855*(u4024+u3928)+u2856*(u844+u4))+u2855*(u670+u677*u2855)+u40&
 &37)))
  u677=-u771
  u3928=9.50822084012189831d-1*u470
  u3899=u973*u2856
  u873=u398+u3899
  value(8)=value(8)+(c(17)*(u2983*(u2856*(u3928+u873)+u3888+u677)*z))
  value(9)=value(9)+(c(17)*(y*(u2856*(u94+u2855*(u782+u462)+u2856*(u396&
 &9+u4122))+u2855*(u414+u456*u2855)+u712)))
  u677=3.48616978011008323d-1*u470
  u3928=u866+u3896
  value(1)=value(1)+(c(18)*(u2983*(u2856*(u677+u3928)+u3891+u3900)*z))
  u677=3.16940694670729944d-1*u3976
  u4037=-3.5215632741192216d-2*u654
  u997=1.40862530964768864d-1*u470
  u595=u2855*(u876+u3794*u2855)
  value(2)=value(2)+(c(18)*(y*((u3839+u2856*(u3835+u909))*u2857+u2856*(&
 &u4037+u4039+(u3756+u997)*u2856)+u595+u677)))
  u677=u585+u3898
  u4037=u597*u2855
  value(3)=value(3)+(c(18)*(u528*((u517+u3897)*u2857+u2856*(u760+u677)+&
 &u4037+u592)))
  u997=u2855*(u4007+u603*u2855)
  value(4)=value(4)+(c(18)*(y*(u2857*(u682+u2856*(u3791+u633)+(u3792+u3&
 &939)*u2857)+u2856*(u4038+u2855*(u601*u2856+u788))+u997+u3796)))
  u3939=-2.77798829537367025d-1*u3776
  u670=4.41952683354902085d-2*u470
  u4=7.57633171465546432d-2*u947
  u3834=-8.08142049563249528d-1*u642
  u760=5.05088780977030955d-2*u950
  u909=-1.57840244055322173d-1*u470
  u4038=u760*u2855
  value(5)=value(5)+(c(18)*(z*(u2857*(u843+u2855*(u4038+u3834)+u2856*(u&
 &3793+u4080)+(u3892+u4038+u3974)*u2857)+u2856*(u3939+u2855*(u3840+u607&
 &)+u2856*(u3757+u3890+u670))+u2855*(u4+u909*u2855)+u3760)))
  value(6)=value(6)+(c(18)*u685)
  u3939=1.69411927340571579d-1*u512
  u4215=-9.03530279149715087d-1*u72
  u3848=-6.21177066915429123d-1*u3776
  u607=2.82353212234285965d-2*u470
  u4145=9.03530279149715087d-1*u642
  u565=1.78823701081714444d-1*u470
  value(7)=value(7)+(c(18)*(z*(u2857*(u4215+u2855*(u3893+u4145)+u863+(u&
 &3893+u597)*u2857)+u2856*(u3848+u2855*(u8+u597)+u2856*(u3758+u607))+u2&
 &855*(u949+u565*u2855)+u3939)))
  u4215=1.05646898223576648d-1*u781
  u3848=-1.05646898223576648d-1*u971
  u3939=-3.5215632741192216d-2*u4224
  u949=5.63450123859075456d-1*u642
  value(8)=value(8)+(c(18)*(x*((u610+u2856*(u3899+u3998))*u2857+u2856*(&
 &u3848+u2855*(u3830+u949)+(u3830+u959)*u2856)+u2855*(u3939+u4078*u2855&
 &)+u4215)))
  value(9)=value(9)+(c(18)*((u2856*(u3850+u2855*(u782+u3759)+u2856*(u39&
 &69+u3779))+u2855*(u893+u4025))*z))
  u3939=-1.49407276290432138d-1*u4224
  u520=-4.98024254301440461d-2*u4224
  u949=7.96838806882304737d-1*u642
  value(1)=value(1)+(c(19)*(x*((u893+u2856*(u3896+u609))*u2857+u2856*(u&
 &3939+u2855*(u3755+u949)+(u3755+u3946)*u2856)+u2855*(u520+u4133*u2855)&
 &+u3850)))
  u3939=-2.11293796447153296d-1*u4081
  u520=-6.33881389341459887d-1*u3776
  value(2)=value(2)+(c(19)*(z*(u2857*(u3939+u4039+u458+u55*u2857)+u2856&
 &*(u520+u4028)+u595+u4127)))
  value(3)=value(3)+(c(19)*(x*(u2857*(u3983+u2856*(u3898+u739)+(u3897+u&
 &877)*u2857)+u2856*(u668+u533+(u3894+u596)*u2856)+u2855*(u664+u3943*u2&
 &855)+u4083)))
  value(4)=value(4)+(c(19)*(z*(u2857*(u757+u4064+u2856*(u3791+u779)+(u3&
 &792+u3790+u650)*u2857)+u2856*(u89+u4030)+u997+u548)))
  value(5)=value(5)+(c(19)*(y*(u2857*(u386+u2855*(u3840+u593)+u2856*(u3&
 &757+u3795)+u723)+u2856*(u4134+u3902+(u3895+u770)*u2856)+u2855*(u790+u&
 &631*u2855)+u657)))
  value(6)=value(6)+(c(19)*u4113)
  u876=-u703
  u3939=-5.6470642446857193d-1*u3776
  u949=-7.52941899291429239d-2*u3776
  u3943=2.44706117269714503d-1*u470
  u596=1.12941284893714386d-1*u3776
  u919=-2.82353212234285965d-1*u470
  value(7)=value(7)+(c(19)*(y*(u2857*(u3939+u2855*(u8+u919)+u2856*(u375&
 &8+u3943)+(u527+u3778)*u2857)+u2856*(u949+u4031)+u2855*(u596+u4032)+u8&
 &76)))
  value(8)=value(8)+(c(19)*(u528*((u959+u873)*u2857+u3794*u2856+u3888+u&
 &520)))
  value(9)=value(9)+(c(19)*(y*((u524+u2855*(u782+u4063)+u2856*(u3837+u3&
 &978))*u2857+u2856*(u3916+u4062+(u3838+u737)*u2856)+u2855*(u882+u787*u&
 &2855)+u655)))
  value(1)=value(1)+(c(20)*(u528*((u3928)*u2857+u4077*u2856+u3891)))
  value(2)=value(2)+(c(20)*(y*(u2857*(u4033+u959*u2856+(u761)*u2857)+u4&
 &48*u2856+u4034)))
  u949=1.35529541872457263d0*u3776
  u3943=-1.01647156404342947d0*u470
  value(3)=value(3)+(c(20)*(u528*(u2857*(u3943+u677+u628*u2857)+u597*u2&
 &856+u4037+u949)))
  u949=-1.05374830610273044d2*u953
  u3943=1.43750787163681003d0*u3776
  u596=-6.38892387394137789d-1*u470
  u4031=(u2857*(u3943+u4035+u504*u2856+u2857*(u4029+u539+u596))+u4079*u&
 &2856+u4036+u949)
  value(4)=value(4)+(c(20)*(y*u4031))
  u4029=-3.78816585732773216d-2*u2858
  u3939=+u4029*(-u3827+5.49778714378213817d2*u977)*u969
  u80=u382*(pd2-2.5d1)
  u876=u3874*(7.42201264410588653d4*u977+5.67185232289765129d1*u80)*u96&
 &9
  u919=-2.39917170964089704d-1*u470
  u520=-6.43988195745714467d-1*u3776
  u3848=4.67207122403753633d-1*u470
  u4215=-1.89408292866386608d-2*u470
  u739=-1.51526634293109286d-1*u2858
  u3915=+u739*(u600+u500*(1.4d1+pd2))*u969
  u4032=1.26272195244257739d-2*u2858
  u4075=u4032*(u511+u3827*(u3910+u4006))*u969
  u535=-6.31360976221288694d-2*u950
  u877=1.89408292866386608d-2*u470
  u779=-6.31360976221288694d-3*u950
  value(5)=value(5)+(c(20)*(z*(u2857*(u876+u2855*(u779*u2855+u4075)+u28&
 &56*(u3757+u3848)+u2857*(u3949+u535*u2855+u919))+u2856*(u520+u4215*u28&
 &56)+u2855*(u3915+u877*u2855)+u3939)))
  value(6)=value(6)+(c(20)*(x*u4031))
  u779=-6.77647709362286316d-1*u3776
  u3939=5.08235782021714737d-1*u470
  u876=-2.82353212234285965d-2*u470
  u919=6.77647709362286316d-1*u3776
  u520=-5.08235782021714737d-1*u470
  value(7)=value(7)+(c(20)*(z*(u2857*(u2855*(u8+u520)+u2856*(u3758+u393&
 &9)+(u527)*u2857)+u2856*(u779+u876*u2856)+u2855*(u919+u607*u2855))))
  value(8)=value(8)+(c(20)*(x*(u2857*(u4021+u3998*u2856+(u3899+u398)*u2&
 &857)+u610*u2856+u70*u2855)))
  u779=1.49407276290432138d-1*u4081
  u3939=8.9644365774259283d-1*u469
  u876=-7.47036381452160691d-2*u975
  u919=-2.61462733508256242d-1*u470
  u520=8.71542445027520806d-2*u950
  value(9)=value(9)+(c(20)*(z*(u2857*(u779+u2855*(u520*u2855+u876)+u388&
 &9+(u3838+u4063)*u2857)+u2856*(u524+u456*u2856)+u2855*(u3939+u919*u285&
 &5)+u655)))
  if ( lmax .eq. 3 ) return
  u779=A9_(p,pd,erfpd,exppd2)
  u3939=p**9
  u876=u2858*u779*u3939
  u919=-8.87122154999207172d3*u876
  u520=+u919
  u3848=3.26568556340659007d4*u779
  u4215=u2858*(u3848+u3992)*u3939
  u3915=-7.47036381452160692d-1*u4215
  u4075=3.73518190726080346d0*u4215
  u535=4.24539123242856709d5*u779
  u877=u2858*(u535-1.13437046457953026d2*u43)*u3939
  u698=4.98024254301440461d-1*u877
  u4063=-8.9644365774259283d-1*u877
  u456=6.36808684864285064d6*u779
  u4031=(u926+1.3d1)
  u607=1.5d1*u4031
  u3894=(u607+u4006)
  u949=u2858*(u456+u3992*u3894)*u3939
  u3943=-4.98024254301440461d-2*u949
  u596=4.98024254301440461d-2*u949
  u650=u596*u2855
  u787=u3943*u2855
  value(1)=value(1)+(c(21)*(u2983*((u3915+u2855*(u787+u698))*u2856+u285&
 &5*(u4075+u2855*(u650+u4063))+u520)))
  u520=-3.76374054924457729d3*u876
  u4063=+u520
  u3978=-1.05646898223576648d-1*u4215
  u737=4.12022903071948927d0*u4215
  u3975=2.11293796447153296d-1*u877
  u4007=-1.47905657513007307d0*u877
  u988=+u4007
  u504=-3.5215632741192216d-2*u949
  u597=1.05646898223576648d-1*u949
  u3794=u504*u2855
  u959=(u3794+u3975)
  u3778=u597*u2855
  value(2)=value(2)+(c(21)*(y*((u3978+u2855*u959)*u2856+u2855*(u737+u28&
 &55*(u3778+u988))+u4063)*z))
  u988=3.35300657809085332d3*u876
  u565=1.69411927340571579d0*u4215
  u4145=-2.82353212234285965d-1*u4215
  u3998=-1.41176606117142983d0*u4215
  u4077=-1.12941284893714386d0*u877
  u4078=1.8823547482285731d-1*u877
  u4133=3.38823854681143158d-1*u877
  u609=1.12941284893714386d-1*u949
  u843=-1.8823547482285731d-2*u949
  u3834=u843*u2855
  u3759=u609*u2855
  value(3)=value(3)+(c(21)*(u2983*((u565+u2855*(u3759+u4077))*u2857+(u4&
 &145+u2855*(u3834+u4078))*u2856+u2855*(u3998+u2855*(u3834+u4133))+u988&
 &)))
  u3946=1.42256021323868609d3*u876
  u858=1.59723096848534447d-1*u4215
  u909=-1.19792322636400836d-1*u4215
  u4080=-1.55730019427321086d0*u4215
  u664=-3.19446193697068895d-1*u877
  u668=2.39584645272801671d-1*u877
  u757=5.59030838969870566d-1*u877
  u593=5.32410322828448158d-2*u949
  u631=-3.99307742121336118d-2*u949
  u4134=u631*u2855
  u3795=u593*u2855
  u770=(u2984*((u858+u2855*(u3795+u664))*u2857+(u909+u2855*(u4134+u668)&
 &)*u2856+u2855*(u4080+u2855*(u4134+u757))+u3946))
  value(4)=value(4)+(c(21)*u770)
  u741=1.48440252882117731d4*u779
  u408=2.83592616144882564d1*exppd2
  u882=3.36725853984687303d-2*u2858*(u741+u408)*u3939
  u3916=1.34955911477071891d3*u876
  u4=5.05088780977030955d-2*u4215
  u470=u2858*(3.85944657493506099d4*u779+u3992)*u3939
  u4224=-1.51526634293109286d-1*u470
  u670=1.70467463579747947d-1*u4215
  u971=u2858*(4.45320758646353191d4*u779+u3992)*u3939
  u601=-3.03053268586218573d-1*u971
  u947=-1.66679297722420215d0*u4215
  u4081=-1.01017756195406191d-1*u877
  u3888=1.6d1*pd2
  u997=u3831*(4.21273437679450119d6*u779+u3992*(1.29d2+u3888))*u3939
  u72=-3.40934927159495895d-1*u877
  u3838=6.31360976221288694d-3*u2858
  u4079=u3838*(8.98063529936812269d6*u779+u3992*(2.75d2+u3888))*u3939
  u873=6.56615415270140241d-1*u877
  u610=1.68362926992343652d-2*u949
  u3891=8.0d0*u4011
  u448=+u4013*(2.16514952853856922d7*u779+u3992*(5.1d1*u4031+u3891))*u3&
 &939
  u655=5.68224878599159824d-2*u949
  u539=-1.89408292866386608d-1*u877
  u3928=-5.05088780977030955d-2*u949
  u677=6.31360976221288694d-2*u949
  u595=6.31360976221288694d-3*u949
  u414=u595*u2855
  u3760=u655*u2855
  u3796=u3928*u2855
  u3839=u677*u2855
  u3900=u610*u2855
  value(5)=value(5)+(c(21)*(u2857*(u3916+u2855*(u2855*(u873+u3796)+u947&
 &)+(u2855*(u4081+u3900)+u4)*u2857)+u2856*(u4224+u2855*(u2855*(u448+u38&
 &39)+u997)+(u2855*(u72+u3760)+u670)*u2856)+u2855*(u601+u2855*(u2855*(u&
 &539+u414)+u4079))+u882))
  u4083=7.11280106619343048d3*u876
  u386=7.98615484242672237d-1*u4215
  u548=-5.98961613182004178d-1*u4215
  u462=-2.99480806591002089d0*u4215
  u3850=-5.32410322828448158d-1*u877
  u4127=3.99307742121336118d-1*u877
  u676=7.18753935818405013d-1*u877
  value(6)=value(6)+(c(21)*(u2985*((u386+u2855*(u3795+u3850))*u2857+(u5&
 &48+u2855*(u4134+u4127))*u2856+u2855*(u462+u2855*(u4134+u676))+u4083))&
 &)
  u70=6.8282516325774156d4*u779
  u620=-2.26874092915906051d2*exppd2
  u89=u2858*(u70+u620)*u3939
  u524=6.274515827428577d-3*u89
  u4112=-2.51475493356813999d3*u876
  u513=+u4112
  u94=u2858*(8.90641517292706383d3*u779+u2966)*u3939
  u3898=-6.77647709362286316d-1*u94
  u444=1.97647248564000175d-1*u4215
  u621=-2.26874092915906051d2*u3621
  u4122=u2858*(4.05736691211121797d4*u779+u621)*u3939
  u3776=-1.69411927340571579d-1*u4122
  u88=2.54117891010857368d0*u4215
  u781=9.47048813387911121d5*u779
  u512=(2.9d1+u4019)
  u3976=u2858*(u781+u3992*u512)*u3939
  u712=1.69411927340571579d-1*u3976
  u517=-3.95294497128000351d-1*u877
  u3774=9.79705669021977021d4*u779
  u780=-u3992*(u4019-3.0d0)
  u3758=-2.82353212234285965d-2*u2858
  u953=+u3758*(u780+u3774)*u3939
  u600=-8.47059636702857895d-1*u877
  u511=1.1d1*u4031
  u977=u2858*(4.6699303556714238d6*u779+u3992*(u511+u3832))*u3939
  u3835=-1.12941284893714386d-1*u977
  u722=6.58824161880000585d-2*u949
  u893=1.50588379858285848d-1*u877
  u603=5.6470642446857193d-2*u949
  u3809=-9.41177374114286549d-3*u949
  u429=u3809*u2855
  u43=u603*u2855
  u3901=u722*u2855
  value(7)=value(7)+(c(21)*((u513+u2855*(u2855*(u600+u43)+u88))*u2857+u&
 &2856*(u3898+u2855*(u2855*(u3835+u43)+u712)+(u2855*(u517+u3901)+u444)*&
 &u2856)+u2855*(u3776+u2855*(u2855*(u893+u429)+u953))+u524))
  u524=-6.27290091540762882d3*u876
  u513=+u524
  u3898=-1.58470347335364972d0*u4215
  u444=+u3898
  u3776=2.6411724555894162d0*u4215
  u712=1.05646898223576648d0*u877
  u517=-6.33881389341459887d-1*u877
  u953=-1.05646898223576648d-1*u949
  u600=3.5215632741192216d-2*u949
  u3835=u953*u2855
  u893=u600*u2855
  value(8)=value(8)+(c(21)*(x*((u444+u2855*(u3835+u712))*u2856+u2855*(u&
 &3776+u2855*(u893+u517))+u513)*z))
  u513=2.96880505764235461d3*u779
  u599=5.67185232289765129d1*exppd2
  u675=3.32016169534293641d-2*u2858*(u513+u599)*u3939
  u3974=2.66136646499762151d3*u876
  u786=3.73518190726080346d-2*u4215
  u973=-6.2098550849944502d3*u876
  u657=-2.9134418876634267d0*u4215
  u3862=-7.47036381452160691d-2*u877
  u469=1.53142458197692942d0*u4215
  u969=1.04585093403302497d0*u877
  u592=1.24506063575360115d-2*u949
  u771=-2.73913339865792253d-1*u877
  u703=-7.47036381452160691d-2*u949
  u3918=u592*u2855
  u3889=u703*u2855
  value(9)=value(9)+(c(21)*(u2856*(u3974+u2855*(u2855*(u969+u3889)+u657&
 &)+(u2855*(u3862+u3918)+u786)*u2856)+u2855*(u973+u2855*(u2855*(u771+u3&
 &918)+u469))+u675))
  u476=4.43561077499603586d2*u876
  u3768=-1.49407276290432138d-1*u4215
  u573=-3.99204969749643227d3*u876
  u685=+u573
  u642=4.48221828871296415d-1*u4215
  u847=2.98814552580864276d-1*u877
  u4113=5.97629105161728553d-1*u4215
  u4024=-4.48221828871296415d-1*u877
  u4064=-4.98024254301440461d-2*u877
  u4039=u4064*u2855
  value(1)=value(1)+(c(22)*(u2856*(u476+u2855*(u2855*(u4024+u650)+u642)&
 &+(u2855*(u847+u787)+u3768)*u2856)+u2855*(u685+u2855*(u4039+u4113))+u6&
 &75))
  u685=9.50822084012189831d-1*u4215
  u4024=1.05646898223576648d-1*u877
  u863=-8.45175185788613183d-1*u877
  u458=-1.05646898223576648d-1*u877
  u533=(u3794+u4024)
  u3902=u458*u2855
  value(2)=value(2)+(c(22)*(x*(u2856*(u685+u2855*(u3778+u863)+u533*u285&
 &6)+u2855*(u685+u3902)+u4063)*z))
  u685=2.07816354034964823d4*u779
  u4063=-5.0196126619428616d-2*u2858*(1.41796308072441282d1*exppd2+u685&
 &)*u3939
  u594=1.69411927340571579d-1*u2858
  u3755=u594*(7.6199329812820435d4*u779+u621)*u3939
  u4062=-3.95294497128000351d-1*u4215
  u629=-4.53748185831812103d2*u3621
  u56=u594*(1.45471447824475376d5*u779+u629)*u3939
  u663=-1.69411927340571579d-1*u2858
  u683=(6.3d1+u3994)
  u3799=+u663*(2.05738190494615174d6*u779+u3992*u683)*u3939
  u3969=7.90588994256000702d-1*u877
  u3899=u3621*(2.1d1+u926)
  u527=u2858*(6.85793968315383915d5*u779-1.13437046457953026d2*u3899)*u&
 &3939
  u684=-1.12941284893714386d-1*u527
  u844=4.9d1*u4031
  u3949=(u844+u3891)
  u950=5.6470642446857193d-2*u2858
  u761=u950*(2.08024170388999788d7*u779+u3992*u3949)*u3939
  u3846=-1.31764832376000117d-1*u949
  u3763=1.31764832376000117d-1*u877
  u3797=u3846*u2855
  u638=u2855*(u761+u3797)
  value(3)=value(3)+(c(22)*(u2856*(u3755+u2855*(u638+u3799)+(u2855*(u39&
 &69+u3797)+u4062)*u2856)+u2855*(u56+u2855*(u3763*u2855+u684))+u4063))
  u55=1.91667716218241337d0*u94
  u389=4.89852834510988511d5*u779
  u4021=u3621*u4014
  u788=u2858*(u389-1.13437046457953026d2*u4021)*u3939
  u654=-1.59723096848534447d-1*u788
  u682=1.59723096848534447d-1*u877
  u3950=-1.31771554900040919d0*u4215
  u3997=1.19792322636400836d-1*u877
  u782=u2858*(2.38395046128681075d6*u779+u3992*(7.3d1+u3994))*u3939
  u3779=-3.99307742121336118d-2*u782
  u723=2.12269561621428355d6*u779
  u3983=5.0d0*u4031
  u790=(u3983+u4011)
  u866=u2858*(u723+u3992*u790)*u3939
  u966=2.12964129131379263d-1*u866
  u398=-5.32410322828448158d-2*u949
  u633=6.38892387394137789d-1*u877
  u3910=1.99653871060668059d-1*u877
  u3851=u2855*(u4134+u633)
  u3840=u398*u2855
  u8=(u2985*(u2857*(u654+u2855*(u3840+u966)+(u3840+u682)*u2857)+u2856*(&
 &u3950+u3851+(u4134+u3997)*u2856)+u2855*(u3779+u3910*u2855)+u55))
  value(4)=value(4)+(c(22)*u8)
  u585=-1.81831961151731144d0*u94
  u5=-4.54579902879327859d-1*u4215
  u3911=-5.05088780977030955d-2*u877
  u4137=1.14298994719230652d6*u779
  u3804=(3.5d1+u4019)
  u612=u3621*u3804
  u913=u2858*(u4137-1.13437046457953026d2*u612)*u3939
  u3767=7.57633171465546432d-2*u913
  u3766=-1.70467463579747947d-1*u877
  u453=1.27361736972857013d6*u779
  u975=7.57633171465546432d-2*u2858
  u445=u975*(u453+u3992*(3.9d1+u4019))*u3939
  u538=4.04071024781624764d-1*u877
  u449=+u3875*(4.03312167080713874d7*u779+u3992*(9.5d1*u4031+u3786))*u3&
 &939
  u729=-2.71485219775154138d-1*u877
  value(5)=value(5)+(c(22)*(u2983*(u2857*(u5+u2855*(u3796+u538)+(u3900+&
 &u3911)*u2857)+u2856*(u3767+u2855*(u3839+u449)+(u3760+u3766)*u2856)+u2&
 &855*(u445+u2855*(u414+u729))+u585)))
  value(6)=value(6)+(c(22)*u770)
  u770=u2858*(3.13373867195581876d4*u779+u3992)*u3939
  u411=-5.08235782021714737d-1*u770
  u95=8.47059636702857895d-1*u4215
  u4038=8.16421390851647518d5*u779
  u615=(2.5d1+u4020)
  u605=u3621*u615
  u925=u2858*(u4038-1.13437046457953026d2*u605)*u3939
  u3896=1.12941284893714386d-1*u925
  u3964=-1.97647248564000175d-1*u877
  u4173=1.63284278170329504d5*u779
  u628=u3621*u3765
  u3765=u2858*(u4173-2.83592616144882564d1*u628)*u3939
  u423=4.51765139574857544d-1*u3765
  u49=-5.6470642446857193d-1*u877
  u881=6.5d1*u4031
  u801=(u881+u4009)
  u978=u2858*(2.75950430107856861d7*u779+u3992*u801)*u3939
  u836=-1.8823547482285731d-2*u978
  u3938=-2.82353212234285965d-2*u877
  u3945=u2855*(u429+u3938)
  value(7)=value(7)+(c(22)*(u2983*((u95+u2855*(u43+u49))*u2857+u2856*(u&
 &3896+u2855*(u43+u836)+(u3901+u3964)*u2856)+u2855*(u423+u3945)+u411)))
  u411=1.25458018308152576d3*u876
  u3896=-3.16940694670729944d-1*u4215
  u3964=6.33881389341459887d-1*u877
  u423=-2.11293796447153296d-1*u877
  u49=u2855*(u893+u423)
  value(8)=value(8)+(c(22)*(y*((u3896+u2855*(u3835+u3964))*u2856+u2855*&
 &(u3896+u49)+u411)*z))
  u836=4.43561077499603586d3*u876
  u702=-3.73518190726080346d-2*u877
  u980=6.22530317876800577d-1*u877
  u60=(u3918+u702)
  u3905=u60*u2856
  value(9)=value(9)+(c(22)*(u2983*(u2856*(u3915+u2855*(u3889+u980)+u390&
 &5)+u2855*(u3915+u2855*u60)+u836)))
  u60=-8.87122154999207171d2*u876
  u4082=+u60
  u465=1.34466548661388924d0*u4215
  u4111=-5.97629105161728553d-1*u877
  value(1)=value(1)+(c(23)*(y*((u3768+u2855*(u787+u847))*u2856+u2855*(u&
 &465+u2855*(u650+u4111))+u4082)*z))
  u4111=-3.16940694670729944d-1*u971
  u864=-u3898
  u3898=1.05646898223576648d-1*u788
  u4199=(2.5d1+u926)
  u3791=u3621*u4199
  u773=u2858*(u4038-1.13437046457953026d2*u3791)*u3939
  u4038=1.05646898223576648d-1*u773
  u483=-u712
  u762=-1.40862530964768864d-1*u866
  u3903=u423*u2855
  value(2)=value(2)+(c(23)*(u2983*((u864+u2855*(u3778+u483))*u2857+u285&
 &6*(u3898+u2855*(u893+u762)+(u893+u458)*u2856)+u2855*(u4038+u3903)+u41&
 &11)))
  u4111=-1.67650328904542666d3*u876
  u4038=3.38823854681143158d-1*u4215
  u3800=-5.6470642446857193d-2*u4215
  u3780=-6.77647709362286316d-1*u877
  u4128=1.12941284893714386d-1*u877
  value(3)=value(3)+(c(23)*(u2984*((u4038+u2855*(u3759+u3780))*u2857+(u&
 &3800+u2855*(u3834+u4128))*u2856+u2855*(u95+u843*u999)+u4111)))
  u91=-3.59376967909202506d-1*u2858
  u859=+u91*(2.86984488905427612d4*u779+u3992)*u3939
  u380=3.59376967909202506d-1*u4215
  u534=-1.59723096848534447d-1*u877
  u61=1.19792322636400836d-1*u788
  u4135=-1.19792322636400836d-1*u877
  u4114=7.98615484242672237d-2*u877
  u630=-1.59723096848534447d-1*u866
  u591=3.99307742121336118d-2*u949
  u3810=-7.98615484242672237d-2*u877
  u447=u591*u2855
  u981=u2855*(u447+u630)
  u3971=(u2983*(u2857*(u380+u2855*(u4134+u4114)+(u3795+u534)*u2857)+u28&
 &56*(u61+u981+(u447+u4135)*u2856)+u2855*(u3997+u3810*u2855)+u859))
  value(4)=value(4)+(c(23)*u3971)
  u776=-9.09159805758655719d-1*u2858*(1.28648219164502033d4*u779+u3995)&
 &*u3939
  u4004=4.04071024781624764d-1*u3765
  u4147=-2.02035512390812382d-1*u877
  u3801=3.03053268586218573d-1*u4215
  u3935=-1.89408292866386608d-2*u877
  u860=u2858*(u3774-5.67185232289765129d1*u4002)*u3939
  u460=6.06106537172437146d-1*u860
  u785=(1.0d1*u4031+u4016)
  u450=+u3879*(u723+u3995*u785)*u3939
  u723=6.73451707969374607d-2*u949
  u499=-1.6415385381753506d-1*u877
  u522=-1.452130245308964d-1*u877
  u678=1.26272195244257739d-2*u949
  u3904=u678*u2855
  u4040=u723*u2855
  value(5)=value(5)+(c(23)*(u2985*(u2857*(u4004+u450*u2855+(u4040+u4147&
 &)*u2857)+u2856*(u3801+u2855*(u3904+u499)+(u414+u3935)*u2856)+u2855*(u&
 &460+u2855*(u414+u522))+u776)))
  u4036=1.3310258070711204d-2*u2858
  u3757=-1.13437046457953026d2*exppd2
  u4126=u4036*(u685+u3757)*u3939
  u3806=-1.06692015992901457d3*u876
  u961=1.19792322636400836d-1*u4215
  u7=-2.39584645272801671d-1*u2858*(u685+u3992)*u3939
  u381=(2.7d1+u4019)
  u4193=u2858*(8.81735102119779319d5*u779+u3992*u381)*u3939
  u69=1.19792322636400836d-1*u4193
  u948=-2.39584645272801671d-1*u877
  u71=2.28597989438461305d5*u779
  u3897=u3621*u924
  u924=u3876*(u71-1.13437046457953026d2*u3897)*u3939
  u4131=2.79515419484935283d-1*u877
  u477=(2.1d1*u4031+u4006)
  u957=u2858*(8.9153215880999909d6*u779+u3992*u477)*u3939
  u406=-3.99307742121336118d-2*u957
  u468=u406+u447
  u667=u2855*(u468)
  value(6)=value(6)+(c(23)*(u2857*(u3806+u2855*(u2855*(u4131+u4134)+u96&
 &1)+(u2855*(u664+u3795)+u858)*u2857)+u2856*(u909+u2855*(u667+u69)+(u28&
 &55*(u948+u447)+u961)*u2856)+u2855*(u7+u924*u2855)+u4126))
  u681=u2858*(3.92542002066044665d4*u779+u3992)*u3939
  u831=-5.08235782021714737d-1*u681
  u3893=u3621*u3869
  u3869=u2858*(u4173-5.67185232289765129d1*u3893)*u3939
  u4173=6.77647709362286316d-1*u3869
  u932=-1.69411927340571579d-1*u877
  u98=1.69411927340571579d-1*u4215
  u598=u3621*(2.9d1+u4020)
  u774=u2858*(u781-1.13437046457953026d2*u598)*u3939
  u781=1.12941284893714386d-1*u774
  u636=(u607+u3832)
  u607=u2858*(u456+u3992*u636)*u3939
  u536=-1.12941284893714386d-1*u607
  u514=-5.6470642446857193d-2*u877
  u604=9.41177374114286549d-3*u949
  u952=-1.41176606117142983d-1*u877
  u3841=u604*u2855
  u4041=u514*u2855
  value(7)=value(7)+(c(23)*(u2985*(u2857*(u4173+u2855*(u3759+u536)+(u43&
 &+u932)*u2857)+u2856*(u98+u4041+(u3841+u3938)*u2856)+u2855*(u781+u2855&
 &*(u429+u952))+u831)))
  u831=8.60953466716282837d4*u779
  u4173=1.1738544247064072d-2*u2858*(u831-u3757)*u3939
  u781=-1.56822522885190721d3*u876
  u536=+u781
  u952=3.16940694670729944d-1*u4215
  u819=-6.33881389341459887d-1*u470
  u459=3.16940694670729944d-1*u4193
  u3980=7.51107679583515716d5*u779
  u870=u2858*(u3980+u3992*(2.3d1+u926))*u3939
  u3849=1.05646898223576648d-1*u870
  u4125=-5.28234491117883239d-1*u877
  u3944=-1.05646898223576648d-1*u957
  u4174=-1.40862530964768864d-1*u877
  value(8)=value(8)+(c(23)*((u536+u2855*(u2855*(u4125+u893)+u864))*u285&
 &7+u2856*(u3896+u2855*(u2855*(u3944+u3778)+u459)+(u2855*(u517+u3778)+u&
 &952)*u2856)+u2855*(u819+u2855*(u4174*u2855+u3849))+u4173))
  u4173=-8.9644365774259283d-1*u4215
  u536=6.72332743306944622d-1*u877
  u819=-1.86759095363040173d-1*u877
  value(9)=value(9)+(c(23)*(x*(u2856*(u4173+u2855*(u3889+u536)+u3905)+u&
 &2855*(u4113+u2855*(u3918+u819))+u4082)*z))
  u459=1.49407276290432138d-1*u877
  u3849=7.47036381452160692d-1*u4215
  u4125=-1.49407276290432138d-1*u877
  u3944=(u787+u459)
  u4174=u3944*u2856
  u3905=u4125*u2855
  value(1)=value(1)+(c(24)*(u2983*(u2856*(u3915+u596*u999+u4174)+u2855*&
 &(u3849+u3905))))
  u3764=-u411
  u463=-3.5215632741192216d-2*u4215
  u457=3.5215632741192216d-2*u877
  u4074=1.37340967690649642d0*u4215
  u4044=-2.81725061929537728d-1*u877
  u4065=-3.16940694670729944d-1*u877
  u4042=u4065*u2855
  value(2)=value(2)+(c(24)*(y*(u2856*(u463+u2855*(u3778+u4044)+(u3794+u&
 &457)*u2856)+u2855*(u4074+u4042)+u3764)*z))
  u463=u2858*(3.46360590058274705d4*u779+u3992)*u3939
  u4044=1.01647156404342947d0*u463
  u4203=-1.69411927340571579d-1*u913
  u784=3.95294497128000351d-1*u877
  u467=1.48588693134999848d7*u779
  u47=3.5d1*u4031
  u473=(u47+u3833)
  u967=u2858*(u467+u3992*u473)*u3939
  u4045=7.52941899291429239d-2*u967
  u4043=u784*u2855
  value(3)=value(3)+(c(24)*(u2983*(u2856*(u4203+u2855*(u3797+u4045)+(u3&
 &797+u784)*u2856)+u2855*(u4203+u4043)+u4044)))
  u4044=u398*u2856
  u4045=u3997*u2855
  u4203=(u2984*(u2857*(u654+u2856*(u4044+u966)+(u4044+u682)*u2857)+u285&
 &6*(u3779+u3851+(u4134+u3910)*u2856)+u2855*(u3950+u4045)+u55))
  value(4)=value(4)+(c(24)*u4203)
  u3851=1.12241951328229101d-2*u2858*(u70+u408)*u3939
  u408=u2858*(6.23449062104894468d4*u779+u621)*u3939
  u772=5.05088780977030955d-2*u408
  u378=2.93911700706593106d5*u779
  u83=(9.0d0+pd2)
  u3837=u3621*u83
  u4129=u2858*(u378-1.13437046457953026d2*u3837)*u3939
  u4073=-3.36725853984687303d-2*u4129
  u4037=1.68362926992343652d-2*u877
  u4003=u3838*(6.17214571483845523d6*u779+u3992*(1.89d2+u3888))*u3939
  u4047=-5.68224878599159824d-2*u877
  u4044=-5.05088780977030955d-2*u2858
  u3895=-9.07496371663624206d2*u3621
  u466=+u4044*(2.87974090591308397d5*u779+u3895)*u3939
  u852=-5.05088780977030955d-2*u4193
  u3951=1.68362926992343652d-2*u957
  u451=-1.68362926992343652d-2*u949
  u821=3.2d1*pd2
  u4061=u3831*(7.93561591907801387d6*u779+u3992*(2.43d2+u821))*u3939
  u4035=3.2d1*u4011
  u608=-6.31360976221288694d-3*u2858
  u452=+u608*(1.0316300694801418d8*u779+u3992*(2.43d2*u4031+u4035))*u39&
 &39
  u931=u3838*(5.91089086976592803d6*u779+u3992*(1.81d2+u3888))*u3939
  u407=1.01017756195406191d-1*u877
  u3906=1.13644975719831965d-1*u949
  u4046=u451*u2855
  value(5)=value(5)+(c(24)*(u2857*(u772+u2855*(u407*u2855+u852)+u2857*(&
 &(u4037+u4046)*u2857+u2855*(u3951+u4046)+u4073))+u2856*(u585+u2855*(u2&
 &855*(u452+u3760)+u4061)+u2856*((u4047+u3760)*u2856+u2855*(u452+u3906*&
 &u2855)+u4003))+u2855*(u466+u2855*(u4047*u2855+u931))+u3851))
  value(6)=value(6)+(c(24)*u8)
  u3851=-1.69411927340571579d-1*u463
  u772=2.82353212234285965d-2*u913
  u4073=-6.58824161880000585d-2*u877
  u4037=1.69411927340571579d-1*u463
  u4003=(u47+u3891)
  u931=u2858*(u467+u3992*u4003)*u3939
  u852=-2.82353212234285965d-2*u931
  u3951=-2.82353212234285965d-2*u913
  u4061=2.82353212234285965d-2*u931
  u931=6.58824161880000585d-2*u877
  u407=-6.58824161880000585d-2*u949
  u3906=u407*u2855
  value(7)=value(7)+(c(24)*(u2856*(u3851+u999*(u3906+u4061)+u2856*((u40&
 &73+u3901)*u2856+u852*u2855+u772))+u2855*(u4037+u2855*(u931*u2855+u395&
 &1))))
  u4061=-u4074
  u772=3.16940694670729944d-1*u877
  u931=3.5215632741192216d-2*u4215
  u3851=2.81725061929537728d-1*u877
  u852=-3.5215632741192216d-2*u877
  u3951=(u3835+u772)
  value(8)=value(8)+(c(24)*(x*(u2856*(u4061+u2855*(u893+u3851)+u3951*u2&
 &856)+u2855*(u931+u852*u2855)+u411)*z))
  u4061=-3.32016169534293641d-2*u2858*(u599+u513)*u3939
  u931=1.77424430999841434d3*u876
  u3851=-3.73518190726080346d-2*u4215
  u4073=-1.24506063575360115d-2*u877
  u4037=-1.56877640104953745d0*u4215
  u466=2.61462733508256242d-1*u877
  u4047=u4073*u2855
  value(9)=value(9)+(c(24)*(u2856*(u931+u2855*(u2855*(u466+u3918)+u4037&
 &)+u2856*((u4073+u3918)*u2856+u2855*(u466+u3889)+u3851))+u2855*(u931+u&
 &2855*(u4047+u3851))+u4061))
  u4037=3.48616978011008323d-1*u4215
  u466=-1.99209701720576184d-1*u877
  value(1)=value(1)+(c(25)*(x*(u2856*(u3768+u2855*(u650+u466)+u4174)+u2&
 &855*(u4037+u4039)+u4082)*z))
  u4037=-7.08981540362206411d0*exppd2
  u3851=u2858*(u513+u4037)*u3939
  u931=-1.87816707953025152d-1*u3851
  u513=4.22587592894306592d-1*u94
  u411=u2858*(u378+u3992*u4072)*u3939
  u378=3.5215632741192216d-2*u411
  u4072=u2858*(7.71889314987012199d3*u779+u2966)*u3939
  u4074=2.11293796447153296d0*u4072
  u8=-2.11293796447153296d-1*u4193
  u4174=3.0d0*u4031
  u42=u2858*(u453+u3992*(u4174+u4006))*u3939
  u569=-3.5215632741192216d-2*u42
  u4194=5.55166545779120312d5*u779
  u759=u3621*u4070
  u4070=u2858*(u4194-1.13437046457953026d2*u759)*u3939
  u3901=-1.05646898223576648d-1*u4070
  u3913=u2858*(8.06624334161427748d6*u779+u3992*(1.9d1*u4031+u4006))*u3&
 &939
  u3847=1.05646898223576648d-1*u3913
  u516=-7.04312654823844319d-2*u949
  u4048=u4024*u2855
  value(2)=value(2)+(c(25)*(u2856*(u513+u2855*(u2855*(u3847+u3835)+u8)+&
 &u2856*((u852+u893)*u2856+u2855*(u569+u516*u2855)+u378))+u2855*(u4074+&
 &u2855*(u4048+u3901))+u931))
  u931=2.03294312808685895d0*u2858
  u513=u931*(1.55037597454656296d4*u779+u3995)*u3939
  u378=-3.38823854681143158d-1*u788
  u8=-6.21177066915429123d-1*u4215
  u569=5.6470642446857193d-2*u877
  u3901=3.6d1*pd2
  u3847=+u3826*(9.2418901444406499d6*u779+u3992*(2.83d2+u3901))*u3939
  u516=4.51765139574857544d-1*u866
  u852=-1.12941284893714386d-1*u949
  u503=3.01176759716571696d-1*u877
  u4171=3.57647402163428889d-1*u877
  u849=u2855*(u3834+u503)
  u3842=u852*u2855
  value(3)=value(3)+(c(25)*(u2985*(u2857*(u378+u2855*(u3842+u516)+(u384&
 &2+u4133)*u2857)+u2856*(u8+u849+(u3834+u569)*u2856)+u2855*(u3847+u4171&
 &*u2855)+u513)))
  u731=7.09880430437930877d-2*u3851
  u489=1.59723096848534447d-1*u408
  u783=-1.06482064565689632d-1*u4129
  u3864=5.32410322828448158d-2*u877
  u410=u2858*(9.8960168588078487d2*u779-3.54490770181103205d0*u3621)*u3&
 &939
  u3782=-1.15000629730944802d1*u410
  u3805=3.99307742121336118d-2*u527
  u4130=-3.99307742121336118d-2*u877
  u875=-7.98615484242672237d-1*u4072
  u815=-1.59723096848534447d-1*u4193
  u763=5.32410322828448158d-2*u957
  u3912=2.39584645272801671d-1*u4193
  u4218=(2.7d1*u4031+u4006)
  u979=u2858*(1.14625563275571312d7*u779+u3992*u4218)*u3939
  u3979=-3.99307742121336118d-2*u979
  u519=3.99307742121336118d-2*u4070
  u639=3.19446193697068895d-1*u877
  u764=7.98615484242672237d-2*u949
  u3914=(u2857*(u489+u2855*(u639*u2855+u815)+u2857*((u3864+u3840)*u2857&
 &+u2855*(u763+u3840)+u783))+u2856*(u3782+u2855*(u2855*(u3979+u447)+u39&
 &12)+u2856*((u4130+u447)*u2856+u2855*(u3979+u764*u2855)+u3805))+u2855*&
 &(u875+u2855*(u4130*u2855+u519))+u731)
  value(4)=value(4)+(c(25)*u3914)
  u552=7.58694625841935067d3*u779
  u4143=-1.81831961151731144d0*u2858*(u552+u2966)*u3939
  u4033=9.0d0*pd2
  u4066=1.34690341593874921d-1*u2858*(u4194+u2966*(6.8d1+u4033))*u3939
  u412=-1.68362926992343652d-1*u877
  u3840=5.05088780977030955d-2*u2858
  u4067=u3840*(1.53487221480109733d6*u779+u3992*(4.7d1+u4017))*u3939
  u724=-2.02035512390812382d-1*u866
  u656=5.05088780977030955d-2*u949
  u426=-1.57840244055322173d-1*u877
  u4030=-5.05088780977030955d-2*u4215
  u867=1.34690341593874921d-1*u877
  u943=-1.38899414768683513d-1*u877
  u671=1.89408292866386608d-2*u877
  u4049=u656*u2856
  value(5)=value(5)+(c(25)*(u2984*(u2857*(u4066+u2855*(u3796+u867)+u285&
 &6*(u4049+u724)+(u4049+u3900+u412)*u2857)+u2856*(u4067+u2855*(u3904+u9&
 &43)+(u414+u426)*u2856)+u2855*(u4030+u2855*(u414+u671))+u4143)))
  value(6)=value(6)+(c(25)*u3971)
  u3971=-5.08235782021714737d-1*u4215
  u4034=u3621*u461
  u461=u2858*(u3848-1.41796308072441282d1*u4034)*u3939
  u3848=2.71059083744914526d0*u461
  u451=u3621*(7.0d1+u4033)
  u4071=u2858*(u4137-5.67185232289765129d1*u451)*u3939
  u3917=7.52941899291429239d-2*u4071
  u700=-2.25882569787428772d-1*u866
  u464=-1.78823701081714444d-1*u877
  u3790=7.34118351809143509d-1*u4215
  u3808=-3.38823854681143158d-1*u877
  u704=-9.4117737411428655d-2*u877
  u4050=u603*u2856
  value(7)=value(7)+(c(25)*(u2984*(u2857*(u3848+u2855*(u43+u3808)+u2856&
 &*(u4050+u700)+(u4050+u932)*u2857)+u2856*(u3917+u704*u2855+(u3841+u464&
 &)*u2856)+u2855*(u3790+u3945)+u3971)))
  u3971=-9.50822084012189831d-1*u770
  u3917=5.28234491117883239d-1*u4215
  u700=3.16940694670729944d-1*u788
  u464=4.7352440669395556d6*u779
  u3790=1.8d1*pd2
  u976=u3621*(1.45d2+u3790)
  u704=u2858*(u464-1.13437046457953026d2*u976)*u3939
  u3945=3.5215632741192216d-2*u704
  u728=-3.5215632741192216d-1*u877
  u765=-4.22587592894306592d-1*u866
  u3792=(u3778+u4065)
  u4051=u728*u2855
  value(8)=value(8)+(c(25)*(u2983*((u3917+u2855*(u893+u728))*u2857+u285&
 &6*(u700+u2855*(u3778+u765)+u3792*u2856)+u2855*(u3945+u4051)+u3971)))
  u3971=-u60
  u3945=-9.96048508602880921d-2*u4215
  u60=-5.97629105161728553d-1*u4215
  u4068=3.237157652959363d-1*u877
  u4069=3.73518190726080346d-2*u877
  value(9)=value(9)+(c(25)*(y*(u2856*(u3945+u2855*(u3889+u4068)+(u3918+&
 &u4073)*u2856)+u2855*(u60+u2855*(u3918+u4069))+u3971)*z))
  u727=-4.48221828871296415d-1*u463
  u912=1.49407276290432138d-1*u788
  u833=(5.5d1+u4017)
  u452=u3621*u833
  u3777=u2858*(1.79612705987362454d6*u779-1.13437046457953026d2*u452)*u&
 &3939
  u4124=4.98024254301440461d-2*u3777
  u4076=-4.98024254301440461d-1*u877
  u766=-1.99209701720576184d-1*u866
  value(1)=value(1)+(c(26)*(u2983*((u3849+u2855*(u650+u4076))*u2857+u28&
 &56*(u912+u2855*(u650+u766)+(u650+u4125)*u2856)+u2855*(u4124+u466*u285&
 &5)+u727)))
  u727=-3.16940694670729944d-1*u470
  u4124=2.11293796447153296d-1*u4129
  u466=1.90164416802437966d0*u4215
  u4052=u600*u2856
  value(2)=value(2)+(c(26)*(u2984*(u2857*(u4124+u2855*(u3778+u517)+u285&
 &6*(u4052+u762)+(u4052+u458)*u2857)+u2856*(u3898+u458*u2856)+u2855*(u4&
 &66+u4042)+u727)))
  u4124=+u663*(4.25528724928737494d4*u779+u3992)*u3939
  u3898=1.75058991585257298d0*u4215
  u4025=5.6470642446857193d-2*u788
  u4042=1.8823547482285731d-2*u2858
  u622=u4042*(1.73081334860549274d6*u779+u3992*(5.3d1+u4017))*u3939
  u4118=-4.89412234539429006d-1*u877
  u679=-7.52941899291429239d-2*u866
  u613=1.8823547482285731d-2*u949
  u651=-3.7647094964571462d-2*u877
  u3761=u613*u2855
  u3967=u2855*(u3761+u679)
  value(3)=value(3)+(c(26)*(u2983*(u2857*(u3898+u2855*(u3834+u4118)+(u3&
 &759+u3808)*u2857)+u2856*(u4025+u3967+(u3761+u514)*u2856)+u2855*(u622+&
 &u651*u2855)+u4124)))
  u532=+u91*(3.33165900913197573d4*u779+u3992)*u3939
  u91=u4027*(2.71051901762746976d6*u779+u3992*(8.3d1+u4033))*u3939
  u376=-1.73033354919245651d-1*u877
  u529=7.98615484242672237d-2*u4215
  u521=-2.92825677555646487d-1*u877
  u3907=u591*u2856
  u561=(u2984*(u2857*(u91+u2855*(u4134+u521)+u2856*(u3907+u630)+(u3907+&
 &u3795+u376)*u2857)+u2856*(u61+u4135*u2856)+u2855*(u529+u4045)+u532))
  value(4)=value(4)+(c(26)*u561)
  u4046=-1.81499274332724841d3*exppd2
  u640=7.01512195801431882d-4*u2858*(3.05786920937162525d5*u779+u4046)*&
 &u3939
  u4045=-1.89408292866386608d-2*u2858
  u4039=-1.81499274332724841d3*u3621
  u871=+u4045*(u389+u4039)*u3939
  u775=2.02035512390812382d-1*u3869
  u526=-6.73451707969374607d-2*u877
  u3795=1.89408292866386608d-2*u2858
  u383=u3795*(1.27658617478621248d5*u779+u629)*u3939
  u4050=u3621*(1.2d1+pd2)
  u4132=u2858*(u3774-2.83592616144882564d1*u4050)*u3939
  u725=-5.05088780977030955d-2*u4132
  u3828=6.31360976221288694d-3*u877
  u4052=-1.13644975719831965d-1*u2858
  u680=+u4052*(u552+u3995)*u3939
  u552=u3795*(5.71494973596153262d6*u779+u3992*(1.75d2+u821))*u3939
  u632=-5.68224878599159824d-2*u4193
  u732=3.3d1*u4031
  u767=u3838*(1.40097910670142714d7*u779+u3992*(u732+u4006))*u3939
  u3970=-6.31360976221288694d-3*u949
  u762=-5.05088780977030955d-2*u860
  u860=-9.4704146433193304d-2*u877
  u824=(6.0d0*u4031+u4011)
  u494=1.01017756195406191d-1*u2858
  u777=u494*(u453+u3995*u824)*u3939
  u453=-1.89408292866386608d-2*u949
  u941=-1.26272195244257739d-2*u949
  u3843=u941*u2855
  u3844=u3970*u2855
  u4053=u453*u2855
  u4054=u3828*u2855
  value(5)=value(5)+(c(26)*(u2857*(u871+u2855*(u2855*(u860+u414)+u552)+&
 &u2857*((u526+u4040)*u2857+u724*u2855+u775))+u2856*(u383+u2855*(u2855*&
 &(u777+u3843)+u632)+u2856*((u3828+u3844)*u2856+u2855*(u767+u4053)+u725&
 &))+u2855*(u680+u2855*(u4054+u762))+u640))
  u383=-3.59376967909202506d-1*u463
  u4040=u3621*u778
  u680=u454*(u4137-1.13437046457953026d2*u4040)*u3939
  u552=-2.79515419484935283d-1*u877
  u762=u3621*u4055
  u860=u3876*(u4137-1.13437046457953026d2*u762)*u3939
  u454=-2.66205161414224079d-2*u967
  u725=9.31718064949784276d-2*u949
  u4055=u725*u2855
  value(6)=value(6)+(c(26)*(u2985*(u2857*(u680+u454*u2855+(u4055+u552)*&
 &u2857)+u860*u2855+u383)))
  u778=1.57346668055044794d5*u779
  u722=4.53748185831812103d2*exppd2
  u576=3.1372579137142885d-3*u2858*(u778+u722)*u3939
  u409=u2858*(1.55367464683283225d5*u779+u629)*u3939
  u4110=-8.47059636702857895d-2*u409
  u3775=5.6470642446857193d-2*u527
  u747=2.82353212234285965d-2*u408
  u3900=-1.8823547482285731d-2*u4129
  u424=9.41177374114286549d-3*u877
  u4049=-5.6470642446857193d-1*u2858
  u768=+u4049*(1.12814592190409475d4*u779+u2966)*u3939
  u4057=8.47059636702857895d-2*u782
  u782=-5.6470642446857193d-2*u979
  u623=-2.82353212234285965d-2*u4193
  u634=9.41177374114286549d-3*u957
  u92=-5.36471103245143333d-1*u877
  u3845=u569*u2855
  u4056=u424*u2855
  value(7)=value(7)+(c(26)*(u2857*(u4110+u2855*(u2855*(u92+u429)+u4057)&
 &+u2857*((u514+u43)*u2857+u2855*(u782+u3759)+u3775))+u2856*(u747+u2855&
 &*(u3845+u623)+u2856*((u424+u429)*u2856+u2855*(u634+u429)+u3900))+u285&
 &5*(u768+u2855*(u4056+u98))+u576))
  u782=-9.50822084012189831d-1*u463
  u623=2.11293796447153296d-1*u925
  u634=1.05646898223576648d-1*u3777
  u92=-7.04312654823844319d-2*u967
  u4110=-4.22587592894306592d-1*u877
  u768=1.40862530964768864d-1*u949
  u4057=u768*u2855
  value(8)=value(8)+(c(26)*(u2985*(u2857*(u623+u2855*(u4057+u92)+u3792*&
 &u2857)+u2855*(u634+u4110*u2855)+u782)))
  u782=5.54451346874504482d2*u876
  u623=-u782
  u634=u2858*(1.33596227593905957d5*u779+u629)*u3939
  u92=-3.73518190726080346d-2*u634
  u3792=u2858*exppd2
  u4110=u3792*u4011*u3939
  u3775=2.82472002361899568d0*u4110
  u747=1.24506063575360115d-2*u877
  u3900=-7.47036381452160691d-2*u2858
  u98=+u3900*(1.9297232874675305d5*u779-6.23903755518741642d2*u3621)*u3&
 &939
  u576=5.60277286089120519d-1*u4215
  u430=1.86759095363040173d-1*u4193
  u497=-1.24506063575360115d-2*u2858
  u701=-1.5d1*u4001
  u455=+u497*(-u3992*(u4006+u701)+u456)*u3939
  u456=-1.24506063575360115d-2*u949
  u4198=2.98814552580864276d-1*u3869
  u726=-2.98814552580864276d-1*u866
  u769=6.22530317876800576d-2*u949
  u3860=-8.71542445027520806d-2*u877
  u645=7.47036381452160691d-2*u949
  u3762=u645*u2855
  u4058=u456*u2855
  u4059=u769*u2855
  value(9)=value(9)+(c(26)*((u623+u2855*(u2855*(u819+u3918)+u576))*u285&
 &7+u2856*(u92+u2855*(u2855*(u726+u3762)+u430)+u2856*((u747+u4058)*u285&
 &6+u2855*(u455+u4059)+u3775))+u2855*(u98+u2855*(u3860*u2855+u4198))+u7&
 &82))
  u3798=-u573
  u573=4.98024254301440461d-2*u877
  u413=-u476
  u476=-4.48221828871296415d-1*u4215
  u4060=4.48221828871296415d-1*u877
  u707=1.49407276290432138d-1*u4215
  u3807=-2.98814552580864276d-1*u877
  u3908=u3807*u2855
  u3909=u707*u2855
  value(1)=value(1)+(c(27)*(u2856*(u3798+u2855*(u3908+u476)+u2856*((u57&
 &3+u787)*u2856+u2855*(u4060+u650)+u60))+u2855*(u413+u3909)+u4061))
  u3798=u504*u2856
  u413=u3798+u3778
  u4060=u517*u2855
  u4061=u952*u2855
  value(2)=value(2)+(c(27)*(x*(u2856*(u952+u4060+u2856*(u413+u3975))+u4&
 &061+u3764)*z))
  value(3)=value(3)+(c(27)*(u2856*(u56+u2855*(u3969*u2855+u3799)+u2856*&
 &((u3763+u3797)*u2856+u638+u684))+u2855*(u3755+u4062*u2855)+u4063))
  u3799=u631*u2856
  u684=u3799+u4134
  u3763=u593*u2856
  u4062=u668*u2855
  u4063=u909*u2855
  u3755=(u2985*((u858+u2856*(u3763+u664))*u2857+u2856*(u4080+u4062+u285&
 &6*(u684+u757))+u4063+u3946))
  value(4)=value(4)+(c(27)*u3755)
  u4134=u595*u2856
  u3764=u3928*u2856
  u3846=u610*u2856
  value(5)=value(5)+(c(27)*(u2983*(u2857*(u5+u2856*(u3764+u538)+(u3846+&
 &u3911)*u2857)+u2856*(u445+u2855*(u3760+u449)+u2856*(u4134+u3839+u729)&
 &)+u2855*(u3767+u3766*u2855)+u585)))
  value(6)=value(6)+(c(27)*u4203)
  u654=5.08235782021714737d-1*u770
  u3950=-8.47059636702857895d-1*u4215
  u3911=-4.51765139574857544d-1*u3765
  u3766=5.6470642446857193d-1*u877
  u729=2.82353212234285965d-2*u877
  u682=-5.6470642446857193d-2*u949
  u585=-1.12941284893714386d-1*u925
  u3910=1.8823547482285731d-2*u978
  u3767=1.97647248564000175d-1*u877
  u445=u604*u2856
  u4203=u682*u2856
  u3765=u682*u2855
  value(7)=value(7)+(c(27)*(u2983*((u3950+u2856*(u4203+u3766))*u2857+u2&
 &856*(u3911+u2855*(u3906+u3910)+u2856*(u445+u3765+u729))+u2855*(u585+u&
 &3767*u2855)+u654)))
  u654=-u520
  u585=-9.50822084012189831d-1*u4215
  u3766=8.45175185788613183d-1*u877
  value(8)=value(8)+(c(27)*(y*(u2856*(u585+u2855*(u893+u3766)+(u3835+u4&
 &024)*u2856)+u2855*(u585+u3902)+u654)*z))
  u3766=u592*u2856
  u585=u3766+u3889
  u3910=u702*u2855
  value(9)=value(9)+(c(27)*(u2983*(u2856*(u3915+u2855*(u3918+u980)+u285&
 &6*(u585+u702))+u2855*(u3915+u3910)+u836)))
  u3911=-3.48616978011008323d-1*u4215
  u980=1.99209701720576184d-1*u877
  value(1)=value(1)+(c(28)*(y*(u2856*(u3911+u2855*(u650+u980)+(u787+u57&
 &3)*u2856)+u2855*(u707+u3905)+u3971)*z))
  u3911=9.50822084012189831d-1*u770
  u3950=-5.28234491117883239d-1*u4215
  u638=-3.5215632741192216d-2*u704
  u704=3.5215632741192216d-1*u877
  u5=-3.16940694670729944d-1*u788
  u520=4.22587592894306592d-1*u866
  u836=u2855*(u3835+u520)
  u3779=u2855*(u5+u772*u2855)
  value(2)=value(2)+(c(28)*(u2983*((u3950+u2856*(u3798+u704))*u2857+u28&
 &56*(u638+u836+(u3835+u704)*u2856)+u3779+u3911)))
  u3911=u852*u2856
  value(3)=value(3)+(c(28)*(u2984*(u2857*(u378+u2856*(u3911+u516)+(u391&
 &1+u4133)*u2857)+u2856*(u3847+u849+(u3834+u4171)*u2856)+u2855*(u8+u384&
 &5)+u513)))
  u8=u2855*(u61+u4135*u2855)
  u55=(u2983*(u2857*(u380+u2856*(u3799+u4114)+(u3763+u534)*u2857)+u2856&
 &*(u3997+u981+(u447+u3810)*u2856)+u8+u859))
  value(4)=value(4)+(c(28)*u55)
  u633=u4134+u3904
  u3847=u656*u2855
  value(5)=value(5)+(c(28)*(u2985*(u2857*(u4066+u2855*(u3847+u724)+u285&
 &6*(u3764+u867)+(u3846+u3847+u412)*u2857)+u2856*(u4030+u2855*(u414+u94&
 &3)+u2856*(u633+u671))+u2855*(u4067+u426*u2855)+u4143)))
  value(6)=value(6)+(c(28)*u3914)
  u4067=5.08235782021714737d-1*u4215
  u815=-u3848
  u638=1.69411927340571579d-1*u877
  u426=-7.34118351809143509d-1*u4215
  u943=-7.52941899291429239d-2*u4071
  u4066=2.25882569787428772d-1*u866
  u849=9.4117737411428655d-2*u877
  u783=1.78823701081714444d-1*u877
  value(7)=value(7)+(c(28)*(u2985*(u2857*(u815+u2855*(u3765+u4066)+u285&
 &6*(u4203+u4133)+(u3765+u638)*u2857)+u2856*(u426+u2855*(u429+u849)+u28&
 &56*(u445+u729))+u2855*(u943+u783*u2855)+u4067)))
  u4067=1.87816707953025152d-1*u3851
  u815=-u4074
  u426=1.05646898223576648d-1*u4070
  u943=-4.22587592894306592d-1*u94
  u4066=2.11293796447153296d-1*u4193
  u849=-1.05646898223576648d-1*u3913
  u783=-3.5215632741192216d-2*u411
  u3805=3.5215632741192216d-2*u42
  u867=7.04312654823844319d-2*u949
  u4074=(u458+u3778)
  value(8)=value(8)+(c(28)*(u2856*(u815+u2855*(u2855*(u3805+u3794)+u406&
 &6)+u2856*(u4074*u2856+u2855*(u849+u867*u2855)+u426))+u2855*(u943+u285&
 &5*(u457*u2855+u783))+u4067))
  value(9)=value(9)+(c(28)*(x*(u2856*(u60+u2855*(u3918+u4068)+u2856*(u5&
 &85+u4069))+u2855*(u3945+u4047)+u3971)*z))
  u867=u2858*(u741+u3995)*u3939
  u4068=-2.98814552580864276d-1*u867
  u943=4.98024254301440461d-2*u788
  u4047=2.98814552580864276d-1*u867
  u867=-4.98024254301440461d-2*u788
  value(1)=value(1)+(c(29)*(u2856*(u4068+u999*(u787+u596)+u2856*((u4064&
 &+u650)*u2856+u787+u943))+u2855*(u4047+u2855*(u573*u2855+u867))))
  u4068=u2858*(3.0677652262304331d4*u779+u3992)*u3939
  u943=9.50822084012189831d-1*u4068
  u4047=u3621*(2.3d1+u4020)
  u867=u2858*(u3980-1.13437046457953026d2*u4047)*u3939
  u60=-2.11293796447153296d-1*u867
  u783=-6.33881389341459887d-1*u4215
  value(2)=value(2)+(c(29)*(u2985*(u2857*(u60+u836+u2856*(u3798+u3975)+&
 &u3951*u2857)+u2856*(u783+u4024*u2856)+u3779+u943)))
  u943=-2.00784506477714464d-1*u3851
  u60=3.38823854681143158d-1*u408
  u783=-2.25882569787428772d-1*u4129
  u815=-5.42118167489829052d0*u410
  u426=1.8823547482285731d-2*u527
  u4067=-1.8823547482285731d-2*u877
  u4069=3.38823854681143158d-1*u2858*(1.68232286599733428d4*u779+u3995)&
 &*u3939
  u4064=-3.38823854681143158d-1*u4193
  u849=1.12941284893714386d-1*u957
  u4066=1.12941284893714386d-1*u4193
  u457=-1.8823547482285731d-2*u979
  u3805=u4042*(u3774+u3992*(3.0d0+u926))*u3939
  u3912=6.77647709362286316d-1*u877
  u4171=3.7647094964571462d-2*u949
  value(3)=value(3)+(c(29)*(u2857*(u60+u2855*(u3912*u2855+u4064)+u2857*&
 &((u4128+u3842)*u2857+u2855*(u849+u3842)+u783))+u2856*(u815+u2855*(u28&
 &55*(u457+u3761)+u4066)+u2856*((u4067+u3761)*u2856+u2855*(u457+u4171*u&
 &2855)+u426))+u2855*(u4069+u2855*(u4067*u2855+u3805))+u943))
  u943=(u2985*(u2857*(u91+u981+u2856*(u3799+u521)+(u3763+u447+u376)*u28&
 &57)+u2856*(u529+u3997*u2856)+u8+u532))
  value(4)=value(4)+(c(29)*u943)
  u60=1.13644975719831965d-1*u2858
  u783=u60*(4.18931380356198928d4*u779+u3992)*u3939
  u815=-8.71278147185378397d-1*u4215
  u426=5.05088780977030955d-2*u877
  u3842=u3621*(2.6d1+u4020)
  u4066=+u4044*(u535-5.67185232289765129d1*u3842)*u3939
  u4064=2.39917170964089704d-1*u877
  u849=3.15680488110644347d-2*u877
  u3805=5.05088780977030955d-2*u866
  u3912=u610*u2857
  u4171=u3912+u3764
  value(5)=value(5)+(c(29)*(u2983*(u2857*(u815+u2855*(u414+u4064)+u2856&
 &*(u4134+u4064)+u2857*(u4171+u3796+u426))+u2856*(u4066+u2855*(u3843+u3&
 &805)+(u3843+u849)*u2856)+u2855*(u4066+u849*u2855)+u783)))
  value(6)=value(6)+(c(29)*u561)
  u783=-7.52941899291429239d-2*u4215
  u815=2.44706117269714503d-1*u877
  u3805=-9.41177374114286549d-3*u877
  u4066=7.52941899291429239d-2*u4215
  u4064=-2.44706117269714503d-1*u877
  u849=u4203+u43
  value(7)=value(7)+(c(29)*(u2983*(u2857*(u2855*(u429+u4064)+u2856*(u44&
 &5+u815)+(u849)*u2857)+u2856*(u783+u3805*u2856)+u2855*(u4066+u4056))))
  u783=-9.50822084012189831d-1*u4068
  u815=2.11293796447153296d-1*u867
  u4066=6.33881389341459887d-1*u4215
  u4064=u597*u2856
  value(8)=value(8)+(c(29)*(u2984*(u2857*(u815+u49+u2856*(u4064+u765)+(&
 &u4064+u4065)*u2857)+u2856*(u700+u4065*u2856)+u2855*(u4066+u3902)+u783&
 &)))
  u783=-6.72332743306944622d-1*u770
  u815=3.73518190726080346d-1*u4215
  u4066=9.96048508602880921d-2*u4071
  u4065=-1.24506063575360115d-1*u877
  u424=-2.36561520793184219d-1*u877
  value(9)=value(9)+(c(29)*(u2983*((u815+u2855*(u3918+u4065)+u2856*(u37&
 &66+u4065))*u2857+u2856*(u4066+u2855*(u3762+u726)+(u3762+u424)*u2856)+&
 &u2855*(u4066+u424*u2855)+u783)))
  u783=2.39051642064691421d0*u461
  u815=8.9644365774259283d-1*u4215
  u4065=u596*u2856
  value(1)=value(1)+(c(30)*(u2984*(u2857*(u783+u2855*(u650+u3807)+u2856&
 &*(u4065+u766)+(u4065+u4125)*u2857)+u2856*(u912+u4125*u2856)+u2855*(u8&
 &15+u3905)+u476)))
  u815=-1.88187027462228865d3*u876
  u4125=+u815
  u424=1.26776277868291978d0*u4215
  u4066=u3975*u2856
  u4067=u3978*u2856
  value(2)=value(2)+(c(30)*(u2983*(u2857*(u424+u4060+u4066+(u413+u423)*&
 &u2857)+u4067+u4061+u4125)))
  u4125=+u663*(5.64072960952047376d4*u779+u3992)*u3939
  u517=1.12941284893714386d-1*u2858*(u4194+u3992*(1.7d1+pd2))*u3939
  u867=1.35529541872457263d0*u4215
  u700=-1.01647156404342947d0*u877
  u4068=u613*u2856
  u4069=u514*u2856
  value(3)=value(3)+(c(30)*(u2984*(u2857*(u517+u2855*(u3834+u700)+u2856&
 &*(u4068+u679)+(u4068+u3759+u932)*u2857)+u2856*(u4025+u4069)+u2855*(u8&
 &67+u3845)+u4125)))
  u912=-2.37093368873114349d2*u876
  u5=9.58338581091206684d-1*u4215
  u573=-5.59030838969870566d-1*u877
  u411=u593*u2857+u684
  u4070=u668*u2856
  u4071=u909*u2856
  u561=(u2983*(u2857*(u5+u4062+u4070+u2857*(u411+u573))+u4071+u4063+u91&
 &2))
  value(4)=value(4)+(c(30)*u561)
  u42=u60*(5.37683582661893113d4*u779+u3992)*u3939
  u60=u3621*(2.8d1+u4020)
  u8=-1.01017756195406191d-1*u2858*(u71-2.83592616144882564d1*u60)*u393&
 &9
  u3782=-1.13644975719831965d-1*u877
  u3980=-6.43988195745714467d-1*u4215
  u4143=4.67207122403753633d-1*u877
  u4062=-7.57633171465546432d-2*u2858
  u3945=(1.5d1+pd2)
  u4068=u3621*u3945
  u513=+u4062*(u389-1.13437046457953026d2*u4068)*u3939
  u489=5.5d1*u4031
  u770=u4032*(2.3349651778357119d7*u779+u3992*(u489+u4006))*u3939
  u408=-6.31360976221288694d-2*u949
  u4072=u408*u2855
  value(5)=value(5)+(c(30)*(u2985*(u2857*(u8+u2855*(u3844+u770)+u2856*(&
 &u4134+u4143)+u2857*(u4171+u4072+u3782))+u2856*(u3980+u3935*u2856)+u28&
 &55*(u513+u671*u2855)+u42)))
  u4171=1.03908177017482411d5*u779
  u410=4.43675269023706798d-3*u2858*(u4171+u722)*u3939
  u4130=+u697*(u4171+u621)*u3939
  u4171=(3.5d1+u926)
  u4060=u3621*u4171
  u981=u3876*(u4137-1.13437046457953026d2*u4060)*u3939
  u4061=-9.31718064949784276d-2*u877
  u4194=-1.19792322636400836d-1*u463
  u535=1.19792322636400836d-1*u913
  u913=(u47+u4006)
  u71=+u3885*(u467+u3992*u913)*u3939
  value(6)=value(6)+(c(30)*(u2857*(u4130+u535*u2855+u2857*((u4061+u4055&
 &)*u2857+u71*u2855+u981))+u4194*u2855+u410))
  u4194=-4.02360789370902398d3*u876
  u535=+u4194
  u71=1.01647156404342947d0*u4215
  u467=-1.12941284893714386d-1*u877
  u3779=-6.77647709362286316d-1*u4215
  u731=5.08235782021714737d-1*u877
  u4063=5.6470642446857193d-1*u4215
  u503=-4.70588687057143275d-1*u877
  u3913=u729*u2855
  u4073=u3938*u2856
  value(7)=value(7)+(c(30)*(u2985*(u2857*(u71+u2855*(u429+u503)+u2856*(&
 &u445+u731)+(u849+u467)*u2857)+u2856*(u3779+u4073)+u2855*(u4063+u3913)&
 &+u535)))
  u467=9.40935137311144324d2*u876
  u3779=u2858*(u70+u621)*u3939
  u731=-3.16940694670729944d-1*u3779
  u4063=u3621*(1.9d1+u926)
  u503=u2858*(6.20480257047252114d5*u779-1.13437046457953026d2*u4063)*u&
 &3939
  u70=1.05646898223576648d-1*u503
  u3851=(3.1d1+u4019)
  u3759=u3621*u3851
  u4030=3.16940694670729944d-1*u2858*(1.01236252465604292d6*u779-1.1343&
 &7046457953026d2*u3759)*u3939
  u476=(2.3d1*u4031+u4006)
  u4055=-1.05646898223576648d-1*u2858
  u398=+u4055*(9.76439983458570431d6*u779+u3992*u476)*u3939
  u520=4.22587592894306592d-1*u4215
  value(8)=value(8)+(c(30)*(u2857*(u731+u2855*(u863*u2855+u4030)+u2857*&
 &(u4074*u2857+u2855*(u398+u4057)+u70))+u2855*(u727+u520*u2855)+u467))
  u731=-6.72332743306944622d-1*u4215
  u70=3.58577463097037132d0*u461
  u4030=-2.24110914435648207d-1*u877
  u398=2.24110914435648207d-1*u4215
  u520=1.49407276290432138d-1*u925
  u458=-2.4901212715072023d-2*u978
  u863=-2.61462733508256242d-1*u877
  u727=8.71542445027520806d-2*u949
  u4074=u727*u2855
  value(9)=value(9)+(c(30)*(u2985*(u2857*(u70+u2855*(u4074+u458)+u2856*&
 &(u3766+u3862)+(u3762+u4030)*u2857)+u2856*(u398+u702*u2856)+u2855*(u52&
 &0+u863*u2855)+u731)))
  u3848=-u919
  u836=-u4075
  u4075=8.9644365774259283d-1*u877
  u3914=u3943*u2856
  u49=u3914+u650
  value(1)=value(1)+(c(31)*(u2983*(u2856*(u836+u4076*u2855+u2856*(u49+u&
 &4075))+u3849*u2855+u3848)))
  u3848=-u524
  u524=-u3776
  u4075=u483*u2855
  u4076=u864*u2855
  value(2)=value(2)+(c(31)*(y*(u2856*(u524+u4075+u2856*(u413+u3964))+u4&
 &076+u3848)*z))
  u3848=u843*u2856
  u524=u3848+u3834
  u3849=u609*u2856
  value(3)=value(3)+(c(31)*(u2983*((u565+u2856*(u3849+u4077))*u2857+u28&
 &56*(u3998+u4078*u2855+u2856*(u524+u4133))+u4145*u2855+u988)))
  u4077=u4127*u2855
  u4078=u548*u2855
  value(4)=value(4)+(c(31)*(u2984*((u386+u2856*(u3763+u3850))*u2857+u28&
 &56*(u462+u4077+u2856*(u684+u676))+u4078+u4083)))
  value(5)=value(5)+(c(31)*(u2857*(u3916+u2856*(u2856*(u873+u3764)+u947&
 &)+(u2856*(u4081+u3846)+u4)*u2857)+u2856*(u601+u2855*(u72*u2855+u997)+&
 &u2856*(u2856*(u539+u3839+u4134)+u2855*(u448+u3760)+u4079))+u2855*(u42&
 &24+u670*u2855)+u882))
  value(6)=value(6)+(c(31)*u3755)
  u4081=-6.274515827428577d-3*u89
  u670=-u4112
  u4224=1.69411927340571579d-1*u4122
  u3850=-u88
  u3760=2.82353212234285965d-2*u2858
  u386=u3760*(u3774+u780)*u3939
  u4079=8.47059636702857895d-1*u877
  u764=-1.50588379858285848d-1*u877
  u539=6.77647709362286316d-1*u94
  u4080=-1.69411927340571579d-1*u3976
  u565=1.12941284893714386d-1*u977
  u4=-1.97647248564000175d-1*u4215
  value(7)=value(7)+(c(31)*((u670+u2856*(u2856*(u4079+u4203)+u3850))*u2&
 &857+u2856*(u4224+u2855*(u4043+u4080)+u2856*(u2856*(u764+u3765+u445)+u&
 &2855*(u565+u3906)+u386))+u2855*(u539+u4*u2855)+u4081))
  u4080=-u737
  u565=-u4007
  u4224=1.05646898223576648d-1*u4215
  u3850=u953*u2856
  u386=u3850+u893
  u4079=u4224*u2855
  value(8)=value(8)+(c(31)*(x*(u2856*(u4080+u3903+u2856*(u386+u565))+u4&
 &079+u654)*z))
  value(9)=value(9)+(c(31)*(u2856*(u973+u2855*(u3862*u2855+u657)+u2856*&
 &(u2856*(u771+u3889+u3766)+u2855*(u969+u3918)+u469))+u2855*(u3974+u786&
 &*u2855)+u675))
  u4080=-u465
  u565=5.97629105161728553d-1*u877
  value(1)=value(1)+(c(32)*(x*(u2856*(u4080+u3908+u2856*(u49+u565))+u39&
 &09+u3971)*z))
  u4080=-1.1738544247064072d-2*u2858*(-u3757+u831)*u3939
  u565=-u781
  u4=6.33881389341459887d-1*u470
  u4081=-1.05646898223576648d-1*u870
  u670=5.28234491117883239d-1*u877
  u764=1.40862530964768864d-1*u877
  u539=-3.16940694670729944d-1*u4193
  u771=1.05646898223576648d-1*u957
  u4064=u2855*(u3964*u2855+u539)
  u462=u2855*(u771+u3835)
  u657=u2855*(u952+u3896*u2855)
  value(2)=value(2)+(c(32)*((u565+u2856*(u2856*(u670+u3798)+u444))*u285&
 &7+u2856*(u4+u4064+u2856*((u764+u3835)*u2856+u462+u4081))+u657+u4080))
  u4080=u3800*u2855
  value(3)=value(3)+(c(32)*(u2985*((u4038+u2856*(u3849+u3780))*u2857+u2&
 &856*(u95+u2855*(u3848+u4128)+u843*u2867)+u4080+u4111)))
  u764=u2855*(u948*u2855+u69)
  u95=u2855*(u909+u961*u2855)
  value(4)=value(4)+(c(32)*(u2857*(u3806+u2856*(u2856*(u4131+u3799)+u96&
 &1)+(u2856*(u664+u3763)+u858)*u2857)+u2856*(u7+u764+u2856*(u2855*(u468&
 &+u3907)+u924))+u95+u4126))
  u4081=u3935*u2855
  value(5)=value(5)+(c(32)*(u2984*(u2857*(u4004+u450*u2856+(u723*u2856+&
 &u4147)*u2857)+u2856*(u460+u2855*(u414+u499)+u2856*(u633+u522))+u2855*&
 &(u3801+u4081)+u776)))
  value(6)=value(6)+(c(32)*u55)
  u4114=5.08235782021714737d-1*u681
  u3780=-6.77647709362286316d-1*u3869
  u664=-1.12941284893714386d-1*u774
  u534=1.12941284893714386d-1*u607
  u4147=1.41176606117142983d-1*u877
  u460=-1.69411927340571579d-1*u4215
  value(7)=value(7)+(c(32)*(u2984*(u2857*(u3780+u2856*(u3911+u534)+(u42&
 &03+u638)*u2857)+u2856*(u664+u2855*(u429+u569)+u2856*(u445+u4147))+u28&
 &55*(u460+u3913)+u4114)))
  u3869=3.16940694670729944d-1*u971
  u3780=-1.05646898223576648d-1*u773
  u664=-1.05646898223576648d-1*u788
  u534=1.40862530964768864d-1*u866
  u773=u2855*(u3794+u534)
  u460=u2855*(u664+u4048)
  value(8)=value(8)+(c(32)*(u2983*((u444+u2856*(u3850+u712))*u2857+u285&
 &6*(u3780+u773+u959*u2856)+u460+u3869)))
  value(9)=value(9)+(c(32)*(y*(u2856*(u4113+u2855*(u3918+u536)+u2856*(u&
 &585+u819))+u2855*(u4173+u3910)+u4082)*z))
  u3869=4.48221828871296415d-1*u463
  u536=-4.98024254301440461d-2*u3777
  u3780=-1.49407276290432138d-1*u788
  u4147=1.99209701720576184d-1*u866
  u788=u2855*(u787+u4147)
  u774=u2855*(u3780+u459*u2855)
  value(1)=value(1)+(c(33)*(u2983*((u3915+u2856*(u3914+u698))*u2857+u28&
 &56*(u536+u788+(u787+u980)*u2856)+u774+u3869)))
  u3869=9.50822084012189831d-1*u463
  u536=-2.11293796447153296d-1*u925
  u4114=-1.05646898223576648d-1*u3777
  u4131=7.04312654823844319d-2*u967
  u681=4.22587592894306592d-1*u877
  u925=-1.40862530964768864d-1*u949
  value(2)=value(2)+(c(33)*(u2984*(u2857*(u536+u2856*(u925*u2856+u4131)&
 &+(u3850+u772)*u2857)+u2856*(u4114+u681*u2856)+u3869)))
  u925=u2855*(u4025+u4041)
  value(3)=value(3)+(c(33)*(u2983*(u2857*(u3898+u2856*(u3848+u4118)+(u3&
 &849+u3808)*u2857)+u2856*(u622+u3967+(u3761+u651)*u2856)+u925+u4124)))
  value(4)=value(4)+(c(33)*(u2984*(u2857*(u680+u454*u2856+(u725*u2856+u&
 &552)*u2857)+u860*u2856+u383)))
  u454=1.89408292866386608d-1*u2858
  u4114=u454*(2.15733167522011102d4*u779+u3995)*u3939
  u680=-7.76574000752185093d-1*u4215
  u860=-1.01017756195406191d-1*u4132
  u4041=2.08349122153025269d-1*u877
  u3808=-1.89408292866386608d-2*u634
  u681=1.51526634293109286d-1*u4193
  u536=-5.05088780977030955d-2*u957
  u622=1.43239448782705802d0*u4110
  u4131=-3.03053268586218573d-1*u877
  value(5)=value(5)+(c(33)*(u2857*(u871+u2855*(u4131*u2855+u681)+u2856*&
 &(u2856*(u4041+u4134)+u680)+u2857*((u526+u3847+u3846)*u2857+u2856*(u42&
 &6+u3764)+u2855*(u536+u3847)+u775))+u2856*(u4114+u2855*(u2855*(u767+u3&
 &844)+u632)+u2856*((u3828+u3843)*u2856+u2855*(u777+u4053)+u860))+u2855&
 &*(u3808+u2855*(u4054+u622))+u640))
  value(6)=value(6)+(c(33)*u943)
  u3828=-3.1372579137142885d-3*u2858*(u722+u778)*u3939
  u61=8.47059636702857895d-2*u409
  u4041=-5.6470642446857193d-2*u527
  u4004=u950*(u741-u3995)*u3939
  u3808=-1.60941330973543d0*u4215
  u681=2.82353212234285965d-2*u634
  u3997=-1.69411927340571579d-1*u4193
  u4048=5.6470642446857193d-2*u957
  u521=2.82353212234285965d-2*u4193
  u526=-2.13528763025153118d0*u4110
  u376=-9.41177374114286549d-3*u957
  value(7)=value(7)+(c(33)*(u2857*(u61+u2855*(u4133*u2855+u3997)+u2856*&
 &(u2856*(u3767+u445)+u3808)+u2857*((u569+u3765)*u2857+u2856*(u4133+u42&
 &03)+u2855*(u4048+u3765)+u4041))+u2856*(u4004+u2855*((u514+u3841)*u285&
 &6+u2855*(u376+u3841)+u521)+u3805*u2867)+u2855*(u681+u2855*(u3805*u285&
 &5+u526))+u3828))
  u3805=3.16940694670729944d-1*u470
  u4004=-2.11293796447153296d-1*u4129
  u526=-u466
  value(8)=value(8)+(c(33)*(u2985*(u2857*(u4004+u773+u2856*(u3850+u3964&
 &)+u533*u2857)+u2856*(u526+u772*u2856)+u460+u3805)))
  value(9)=value(9)+(c(33)*((u623+u2856*(u2856*(u819+u3766)+u576))*u285&
 &7+u2856*(u98+u2855*(u2855*(u455+u4058)+u430)+u2856*((u3860+u3762)*u28&
 &56+u2855*(u726+u4059)+u4198))+u2855*(u92+u2855*(u747*u2855+u3775))+u7&
 &82))
  u4004=-u783
  value(1)=value(1)+(c(34)*(u2985*(u2857*(u4004+u788+u2856*(u3914+u847)&
 &+u3944*u2857)+u2856*(u4173+u459*u2856)+u774+u642)))
  u4004=-u467
  u526=3.16940694670729944d-1*u3779
  u3805=-1.05646898223576648d-1*u503
  u681=-u815
  u3808=-u424
  value(2)=value(2)+(c(34)*(u2857*(u526+u4064+u2856*(u4066+u3808)+u2857&
 &*((u4024+u3835)*u2857+u2856*(u3975+u3798)+u462+u3805))+u2856*(u681+u4&
 &067)+u657+u4004))
  value(3)=value(3)+(c(34)*(u2985*(u2857*(u517+u3967+u2856*(u3848+u700)&
 &+(u3849+u3761+u932)*u2857)+u2856*(u867+u569*u2856)+u925+u4125)))
  value(4)=value(4)+(c(34)*(u2857*(u4130+u764+u2856*(u4070+u5)+u2857*((&
 &u4061+u447+u3763)*u2857+u2856*(u573+u3799)+u667+u981))+u2856*(u912+u4&
 &071)+u95+u410))
  value(5)=value(5)+(c(34)*(u2984*(u2857*(u8+u2855*(u414+u4143)+u2856*(&
 &u3970*u2856+u770)+u2857*(u3912+u408*u2856+u3796+u3782))+u2856*(u513+u&
 &671*u2856)+u2855*(u3980+u4081)+u42)))
  value(6)=value(6)+(c(34)*u561)
  u573=-u4194
  u526=-u71
  u3805=-5.6470642446857193d-1*u4215
  u948=4.70588687057143275d-1*u877
  u4004=6.77647709362286316d-1*u4215
  u8=-5.08235782021714737d-1*u877
  value(7)=value(7)+(c(34)*(u2984*(u2857*(u526+u2855*(u429+u8)+u2856*(u&
 &445+u948)+(u849+u4128)*u2857)+u2856*(u3805+u4073)+u2855*(u4004+u3913)&
 &+u573)))
  value(8)=value(8)+(c(34)*(u2983*(u2857*(u3808+u3903+u3964*u2856+(u386&
 &+u3975)*u2857)+u3896*u2856+u4079+u681)))
  value(9)=value(9)+(c(34)*(u2984*(u2857*(u70+u2855*(u3918+u3862)+u2856&
 &*(u727*u2856+u458)+(u645*u2856+u4030)*u2857)+u2856*(u520+u863*u2856)+&
 &u2855*(u398+u3910)+u731)))
  value(1)=value(1)+(c(35)*(u2983*(u2857*(u3908+u847*u2856+(u49)*u2857)&
 &+u3768*u2856+u3909)))
  value(2)=value(2)+(c(35)*(u2984*(u2857*(u4075+u704*u2856+(u413)*u2857&
 &)+u3950*u2856+u4076)))
  u3805=4.40471011085486105d0*u4215
  u948=-1.5811779885120014d0*u877
  value(3)=value(3)+(c(35)*(u2983*(u2857*(u3805+u4128*u2855+u4128*u2856&
 &+u2857*(u609*u2857+u524+u948))+u3800*u2856+u4080+u535)))
  u3805=-9.48373475492457397d3*u876
  u948=3.99307742121336119d0*u4215
  u483=-9.58338581091206684d-1*u877
  u8=(u2857*(u948+u4077+u4127*u2856+u2857*(u411+u483))+u548*u2856+u4078&
 &+u3805)
  value(4)=value(4)+(c(35)*(u2984*u8))
  u573=1.40302439160286376d-3*u2858*(u741-u4046)*u3939
  u526=+u4029*(-u621+1.53388261311521655d5*u779)*u3939
  u4004=u4032*(u464-u3992*(u926-1.45d2))*u3939
  u3808=-3.5777121985873026d-1*u877
  u681=1.91187541259185178d3*u876
  u534=-2.04560956295697537d0*u4215
  u4147=7.19751512892269111d-1*u877
  u702=1.89408292866386608d-2*u4215
  u700=-3.78816585732773216d-2*u877
  u4077=1.51526634293109286d-1*u2858
  u3935=u4077*(u685+u2966)*u3939
  u932=+u4029*(2.64520530635933796d6*u779+u3992*(8.1d1+u4019))*u3939
  u819=3.9d1*u4031
  u376=u3874*(1.65570258064714117d7*u779+u3992*(u819+u3832))*u3939
  u521=-1.89408292866386608d-2*u4215
  u3807=3.78816585732773216d-2*u877
  value(5)=value(5)+(c(35)*(u2857*(u526+u2855*(u3807*u2855+u932)+u2856*&
 &(u700*u2856+u534)+u2857*(u2857*(u3808+u4072+u3764+u3912)+u2856*(u4147&
 &+u4134)+u2855*(u376+u3844)+u4004))+u2856*(u681+u702*u2856)+u2855*(u39&
 &35+u521*u2855)+u573))
  value(6)=value(6)+(c(35)*(u2985*u8))
  u376=2.01180394685451199d3*u876
  u3935=-2.20235505542743053d0*u4215
  u4004=+u3935
  u534=2.82353212234285965d-2*u4215
  u4147=-u376
  u702=-u3935
  u3935=-7.90588994256000702d-1*u877
  u932=-2.82353212234285965d-2*u4215
  value(7)=value(7)+(c(35)*(u2857*(u2855*(u3845+u702)+u2856*(u4069+u400&
 &4)+u2857*((u43+u4203)*u2857+u2856*(u3969+u445)+u2855*(u3935+u429)))+u&
 &2856*(u376+u534*u2856)+u2855*(u4147+u932*u2855)))
  value(8)=value(8)+(c(35)*(u2985*(u2857*(u4051+u712*u2856+(u386)*u2857&
 &)+u444*u2856+u3917*u2855)))
  u4004=6.65341616249405378d2*u876
  u534=-2.24110914435648207d-1*u3779
  u4147=7.47036381452160691d-2*u503
  u702=-u4004
  u3935=-8.9644365774259283d-1*u94
  u932=2.24110914435648207d-1*u3976
  u376=-1.49407276290432138d-1*u977
  u521=2.61462733508256242d-1*u4215
  u3807=-5.22925467016512484d-1*u877
  value(9)=value(9)+(c(35)*(u2857*(u534+u2855*(u3807*u2855+u932)+u2856*&
 &(u3862*u2856+u642)+u2857*((u3862+u3762)*u2857+u2856*(u3862+u3766)+u28&
 &55*(u376+u4074)+u4147))+u2856*(u702+u786*u2856)+u2855*(u3935+u521*u28&
 &55)+u4004))
  if ( lmax .eq. 4 ) return
  u534=A10_(p,pd,erfpd,exppd2)
  u4147=p**10
  u702=u2858*u534*u4147
  u3935=8.87122154999207172d3*u702
  u932=3.26568556340659007d4*u534
  u376=u2858*(u932+u3788)*u4147
  u521=7.47036381452160692d-1*u376
  u3807=-1.41936912475910531d1*u376
  u4004=+u3807
  u3808=u2858*(4.24539123242856709d5*u534-1.13437046457953026d2*u3977)*&
 &u4147
  u681=-2.24110914435648207d0*u3808
  u408=+u681
  u4041=8.21740019597376761d0*u3808
  u483=6.36808684864285064d6*u534
  u948=u2858*(u483+u3788*u3894)*u4147
  u728=7.47036381452160692d-1*u948
  u3805=-1.24506063575360115d0*u948
  u573=+u3805
  u526=(u4006-u701)
  u4074=pd2**3
  u645=8.0d0*u4074
  u700=u2858*(1.08257476426928461d8*u534+u3788*(1.7d1*u526+u645))*u4147
  u8=-4.98024254301440461d-2*u700
  u3782=4.98024254301440461d-2*u700
  u4061=u3782*u2855
  u4030=u8*u2855
  value(1)=value(1)+(c(36)*(y*((u521+u2855*(u2855*(u728+u4030)+u408))*u&
 &2856+u2855*(u4004+u2855*(u2855*(u573+u4061)+u4041))+u3935)))
  u4004=-9.50822084012189831d0*u376
  u573=+u4004
  u863=-5.28234491117883239d-1*u3808
  u4048=1.00364553312397816d1*u3808
  u664=3.5215632741192216d-1*u948
  u3780=-2.11293796447153296d0*u948
  u764=+u3780
  u3997=-3.5215632741192216d-2*u700
  u668=1.05646898223576648d-1*u700
  u3975=u668*u2855
  u4024=u3997*u2855
  value(2)=value(2)+(c(36)*(u2983*((u863+u2855*(u4024+u664))*u2856+u285&
 &5*(u4048+u2855*(u3975+u764))+u573)*z))
  u573=-3.35300657809085332d3*u702
  u764=-1.69411927340571579d0*u376
  u61=2.82353212234285965d-1*u376
  u569=5.36471103245143333d0*u376
  u4133=5.08235782021714737d0*u3808
  u4127=-8.47059636702857895d-1*u3808
  u3964=-3.10588533457714561d0*u3808
  u4128=+u3964
  u459=-1.69411927340571579d0*u948
  u772=2.82353212234285965d-1*u948
  u729=4.70588687057143275d-1*u948
  u847=1.12941284893714386d-1*u700
  u671=-1.8823547482285731d-2*u700
  u4025=u671*u2855
  u4058=u729+u4025
  u3828=u847*u2855
  value(3)=value(3)+(c(36)*(y*((u764+u2855*(u2855*(u459+u3828)+u4133))*&
 &u2857+(u61+u2855*(u2855*(u772+u4025)+u4127))*u2856+u2855*(u569+u2855*&
 &(u2855*(u4058)+u4128))+u573)))
  u3969=3.59376967909202507d0*u376
  u69=7.98615484242672237d-1*u3808
  u704=-5.98961613182004178d-1*u3808
  u4131=-3.79342355015269313d0*u3808
  u536=-5.32410322828448158d-1*u948
  u622=3.99307742121336118d-1*u948
  u680=7.98615484242672237d-1*u948
  u860=5.32410322828448158d-2*u700
  u747=-3.99307742121336118d-2*u700
  u430=u747*u2855
  u3767=u860*u2855
  u4114=(u528*((u69+u2855*(u3767+u536))*u2857+(u704+u2855*(u430+u622))*&
 &u2856+u2855*(u4131+u2855*(u430+u680))+u3969))
  value(4)=value(4)+(c(36)*u4114)
  u4124=2.07816354034964823d4*u534
  u426=1.51526634293109286d0*u2858*(u4124+u3827)*u4147
  u4143=3.78816585732773216d0*u376
  u513=2.52544390488515477d-1*u3808
  u610=1.01236252465604292d6*u534
  u520=u2858*(u610+u3788*u3851)*u4147
  u877=-3.78816585732773216d-1*u520
  u949=8.52337317898739736d-1*u3808
  u788=2.97177386269999696d6*u534
  u527=-1.26272195244257739d-1*u2858
  u925=+u527*(u788+u3788*(9.1d1+u3994))*u4147
  u3777=-4.29325463830476312d0*u3808
  u3976=-1.68362926992343652d-1*u948
  u4075=6.31360976221288694d-2*u2858
  u773=u4075*(6.75017205956142168d7*u534+u3788*(1.59d2*u4031+u4035))*u4&
 &147
  u460=-5.68224878599159824d-1*u948
  u774=u3838*(1.93165301075499803d8*u534+u3788*(4.55d2*u4031+u4035))*u4&
 &147
  u775=9.59668683856358814d-1*u948
  u4198=1.68362926992343652d-2*u700
  u517=u2858*(2.10146866005214071d8*u534+u3788*(3.3d1*u526+u645))*u4147
  u70=-5.05088780977030955d-2*u517
  u4129=5.68224878599159824d-2*u700
  u461=-2.39917170964089704d-1*u948
  u3869=-5.05088780977030955d-2*u700
  u4132=6.31360976221288694d-2*u700
  u909=6.31360976221288694d-3*u700
  u3915=u4198*u2855
  u548=u3976+u3915
  u3896=u909*u2855
  u444=u4129*u2855
  u3768=u3869*u2855
  u3800=u4132*u2855
  value(5)=value(5)+(c(36)*(x*(u2857*(u4143+u2855*(u2855*(u775+u3768)+u&
 &3777)+(u2855*(u548)+u513)*u2857)+u2856*(u877+u2855*(u2855*(u70+u3800)&
 &+u773)+(u2855*(u460+u444)+u949)*u2856)+u2855*(u925+u2855*(u2855*(u461&
 &+u3896)+u774))+u426)))
  u859=-7.11280106619343048d3*u702
  u532=-7.98615484242672237d-1*u376
  u4064=5.98961613182004178d-1*u376
  u4173=1.13802706504580794d1*u376
  u4224=2.39584645272801671d0*u3808
  u601=-1.79688483954601253d0*u3808
  u591=-6.58857774500204595d0*u3808
  u462=-7.98615484242672237d-1*u948
  u657=5.98961613182004178d-1*u948
  u776=9.98269355303340296d-1*u948
  value(6)=value(6)+(c(36)*(z*((u532+u2855*(u2855*(u462+u3767)+u4224))*&
 &u2857+(u4064+u2855*(u2855*(u657+u430)+u601))*u2856+u2855*(u4173+u2855&
 &*(u2855*(u776+u430)+u591))+u859)))
  u7=u2858*(u4124+u3788)*u4147
  u4124=8.47059636702857895d-1*u7
  u871=-5.92941745692000526d0*u376
  u383=+u871
  u92=2.28597989438461305d5*u534
  u98=u2858*(u92-1.13437046457953026d2*u968)*u4147
  u3950=-1.69411927340571579d0*u98
  u4125=+u3950
  u656=9.88236242820000877d-1*u3808
  u4130=u2858*(u92+u3788*(7.0d0+u4019))*u4147
  u731=-2.82353212234285965d-1*u4130
  u924=5.92941745692000526d0*u3808
  u858=1.48588693134999848d7*u534
  u952=u2858*(u858+u3788*u4003)*u4147
  u961=2.82353212234285965d-1*u952
  u864=-6.58824161880000585d-1*u948
  u380=-u3788*(u3891-3.5d1*u4001)
  u529=+u3758*(u380+u858)*u4147
  u707=-1.18588349138400105d0*u948
  u3905=+u707
  u642=4.45766079404999545d7*u534
  u592=2.0d0*u4074
  u95=u2858*(u642+u3788*(7.0d0*u526+u592))*u4147
  u5=-2.25882569787428772d-1*u95
  u565=6.58824161880000585d-2*u700
  u4=2.25882569787428772d-1*u948
  u670=5.6470642446857193d-2*u700
  u386=-9.41177374114286549d-3*u700
  u4065=u386*u2855
  u4038=u670*u2855
  u3801=u565*u2855
  value(7)=value(7)+(c(36)*(x*((u383+u2855*(u2855*(u3905+u4038)+u924))*&
 &u2857+u2856*(u4125+u2855*(u2855*(u5+u4038)+u961)+(u2855*(u864+u3801)+&
 &u656)*u2856)+u2855*(u731+u2855*(u2855*(u4+u4065)+u529))+u4124)))
  u4124=6.27290091540762882d3*u702
  u383=1.58470347335364972d0*u376
  u4125=-1.00364553312397816d1*u376
  u656=+u4125
  u731=-4.75411042006094916d0*u3808
  u961=+u731
  u864=5.81057940229671564d0*u3808
  u529=1.58470347335364972d0*u948
  u3905=-8.803908185298054d-1*u948
  u5=-1.05646898223576648d-1*u700
  u4=3.5215632741192216d-2*u700
  u3904=u4*u2855
  u3898=u5*u2855
  value(8)=value(8)+(c(36)*(((u383+u2855*(u2855*(u529+u3898)+u961))*u28&
 &56+u2855*(u656+u2855*(u2855*(u3905+u3904)+u864))+u4124)*z))
  u656=1.33068323249881076d4*u702
  u3905=6.72332743306944622d0*u376
  u576=1.86759095363040173d-1*u3808
  u867=-8.21740019597376761d0*u376
  u42=-7.09684562379552657d0*u3808
  u398=-1.24506063575360115d-1*u948
  u539=3.17490462117168294d0*u3808
  u777=1.49407276290432138d0*u948
  u4215=1.24506063575360115d-2*u700
  u463=-3.73518190726080346d-1*u948
  u94=-7.47036381452160691d-2*u700
  u763=u4215*u2855
  u4145=u2855*(u398+u763)
  u675=u94*u2855
  value(9)=value(9)+(c(36)*(x*(u2856*(u3905+u2855*(u2855*(u777+u675)+u4&
 &2)+(u4145+u576)*u2856)+u2855*(u867+u2855*(u2855*(u463+u763)+u539))+u6&
 &56)))
  u767=-7.47036381452160692d-1*u376
  u595=-7.47036381452160692d-1*u3808
  u882=-3.73518190726080346d0*u376
  u4122=+u882
  u4126=-u681
  u681=4.98024254301440461d-1*u948
  u640=8.9644365774259283d-1*u3808
  u409=-7.47036381452160692d-1*u948
  u410=-4.98024254301440461d-2*u948
  u89=u767+u2855*(u2855*(u409+u4061)+u4126)
  u4082=u410*u2855
  value(1)=value(1)+(c(37)*(x*(u2856*(u89+(u2855*(u681+u4030)+u595)*u28&
 &56)+u2855*(u4122+u2855*(u4082+u640))+u3935)))
  u4122=3.76374054924457729d3*u702
  u640=-9.50822084012189831d-1*u376
  u912=-1.05646898223576648d-1*u3808
  u973=-4.12022903071948927d0*u376
  u4111=+u973
  u3806=3.48634764137802938d0*u3808
  u623=2.11293796447153296d-1*u948
  u535=1.47905657513007307d0*u3808
  u651=-1.37340967690649642d0*u948
  u447=+u651
  u3971=-1.05646898223576648d-1*u948
  u723=(u623+u4024)
  u3849=u2855*u723
  u4083=u3971*u2855
  value(2)=value(2)+(c(37)*((u2856*(u640+u2855*(u2855*(u447+u3975)+u380&
 &6)+(u3849+u912)*u2856)+u2855*(u4111+u2855*(u4083+u535))+u4122)*z))
  u640=u2858*(3.46360590058274705d4*u534+u3788)*u4147
  u4111=-1.69411927340571579d0*u640
  u447=2.82353212234285965d-1*u2858
  u3974=u447*(u788-1.13437046457953026d2*u998)*u4147
  u654=-1.97647248564000175d0*u3808
  u552=+u654
  u876=u2858*(2.05738190494615174d6*u534+u3788*u683)*u4147
  u3774=2.82353212234285965d-1*u876
  u741=3.26895124896999666d7*u534
  u3755=-2.82353212234285965d-1*u2858
  u831=7.7d1*u4031
  u464=+u3755*(u741+u3788*(u831+u3786))*u4147
  u778=1.31764832376000117d0*u948
  u770=7.0d0*u4031
  u779=-4.51765139574857544d-1*u2858*(u788+u3788*(u770+u4011))*u4147
  u727=1.6d1*u4074
  u4110=u950*(4.0118947146449959d8*u534+u3788*(6.3d1*u526+u727))*u4147
  u4003=-1.31764832376000117d-1*u700
  u780=1.31764832376000117d-1*u948
  u3860=u4003*u2855
  u3775=u2855*(u4110+u3860)
  value(3)=value(3)+(c(37)*(x*(u2856*(u3974+u2855*(u3775+u464)+(u2855*(&
 &u778+u3860)+u552)*u2856)+u2855*(u3774+u2855*(u780*u2855+u779))+u4111)&
 &))
  u4113=u2858*(8.90641517292706383d3*u534+u500)*u4147
  u757=-1.91667716218241337d0*u4113
  u873=4.89852834510988511d5*u534
  u3946=u2858*(u873-1.13437046457953026d2*u4197)*u4147
  u503=1.59723096848534447d-1*u3946
  u3938=-1.59723096848534447d-1*u3808
  u3847=1.31771554900040919d0*u376
  u3977=-1.19792322636400836d-1*u3808
  u3762=1.19792322636400836d-1*u2858
  u3834=u3762*(3.95147953172197399d6*u534+u3788*(1.21d2+u3888))*u4147
  u4197=-1.59723096848534447d-1*u952
  u631=3.19446193697068895d-1*u948
  u3939=-3.23439271118282256d0*u3808
  u632=2.39584645272801671d-1*u948
  u953=2.08024170388999788d7*u534
  u968=u2858*(u953+u3788*u3949)*u4147
  u726=-7.98615484242672237d-2*u968
  u3943=1.59202171216071266d8*u534
  u701=u2858*(u3943+u3788*(2.5d1*u526+u645))*u4147
  u843=5.32410322828448158d-2*u701
  u504=-5.32410322828448158d-2*u700
  u682=8.38546258454805848d-1*u948
  u683=1.99653871060668059d-1*u948
  u3851=u504*u2855
  u3809=u2855*(u843+u3851)
  u3928=u2855*(u682+u430)
  u737=(z*(u2857*(u503+u2855*(u3809+u4197)+(u2855*(u631+u3851)+u3938)*u&
 &2857)+u2856*(u3847+u2855*(u3928+u3939)+(u2855*(u632+u430)+u3977)*u285&
 &6)+u2855*(u3834+u2855*(u683*u2855+u726))+u757))
  value(4)=value(4)+(c(37)*u737)
  u4007=1.81831961151731144d0*u4113
  u4113=4.54579902879327859d-1*u376
  u4112=5.05088780977030955d-2*u3808
  u88=1.14298994719230652d6*u534
  u524=u2858*(u88+u3788*u3804)*u4147
  u3804=-7.57633171465546432d-2*u524
  u3776=1.70467463579747947d-1*u3808
  u712=-2.2728995143966393d-1*u876
  u876=-1.66679297722420215d0*u3808
  u465=-1.01017756195406191d-1*u948
  u781=u3831*(7.0048955335071357d7*u534+u3788*(1.65d2*u4031+u4035))*u41&
 &47
  u466=-3.40934927159495895d-1*u948
  u977=9.6d1*u4011
  u782=u3838*(2.89960221174871132d8*u534+u3788*(6.83d2*u4031+u977))*u41&
 &47
  u783=6.56615415270140241d-1*u948
  u815=9.55213027296427596d7*u534
  u424=1.5d1*u526
  u3766=4.0d0*u4074
  u4194=u2858*(u815+u3788*(u424+u3766))*u4147
  u71=-1.01017756195406191d-1*u4194
  u467=-3.15680488110644347d-1*u948
  value(5)=value(5)+(c(37)*(y*(u2857*(u4113+u2855*(u2855*(u783+u3768)+u&
 &876)+(u2855*(u465+u3915)+u4112)*u2857)+u2856*(u3804+u2855*(u2855*(u71&
 &+u3800)+u781)+(u2855*(u466+u444)+u3776)*u2856)+u2855*(u712+u2855*(u28&
 &55*(u467+u3896)+u782))+u4007)))
  value(6)=value(6)+(c(37)*u4114)
  u4114=u2858*(3.13373867195581876d4*u534+u3788)*u4147
  u698=5.08235782021714737d-1*u4114
  u423=-8.47059636702857895d-1*u376
  u4118=8.16421390851647518d5*u534
  u703=u2858*(u4118+u3788*u615)*u4147
  u615=-1.12941284893714386d-1*u703
  u852=1.97647248564000175d-1*u3808
  u3799=u2858*(7.02122396132416865d6*u534+u3788*(2.15d2+u3901))*u4147
  u561=-5.6470642446857193d-2*u3799
  u943=2.54117891010857368d0*u3808
  u600=2.4d1*u4011
  u980=u2858*(4.88219991729285216d7*u534+u3788*(1.15d2*u4031+u600))*u41&
 &47
  u630=5.6470642446857193d-2*u980
  u411=-3.95294497128000351d-1*u948
  u522=3.60858254756428203d7*u534
  u3779=8.5d1*u4031
  u981=u2858*(u522+u3788*(u3779+u600))*u4147
  u3967=2.82353212234285965d-2*u981
  u638=-8.47059636702857895d-1*u948
  u49=3.18404342432142532d7*u534
  u667=1.0d1*u526
  u604=3.0d0*u4074
  u849=u2858*(u49+u3827*(u667+u604))*u4147
  u445=-3.01176759716571696d-1*u849
  u836=3.7647094964571462d-2*u948
  value(7)=value(7)+(c(37)*(y*((u423+u2855*(u2855*(u638+u4038)+u943))*u&
 &2857+u2856*(u615+u2855*(u2855*(u445+u4038)+u630)+(u2855*(u411+u3801)+&
 &u852)*u2856)+u2855*(u561+u2855*(u2855*(u836+u4065)+u3967))+u698)))
  u698=1.05646898223576648d0*u376
  u423=-1.58470347335364972d0*u3808
  u615=+u423
  u852=5.28234491117883239d-1*u3808
  u561=1.05646898223576648d0*u948
  u630=-4.22587592894306592d-1*u948
  value(8)=value(8)+(c(37)*(u2983*((u615+u2855*(u3898+u561))*u2856+u285&
 &5*(u852+u2855*(u3904+u630))+u698)*z))
  u3967=-4.43561077499603586d3*u702
  u638=3.73518190726080346d-2*u3808
  u445=-u882
  u836=-2.61462733508256242d0*u3808
  u882=-7.47036381452160691d-2*u948
  u585=-5.60277286089120519d-1*u3808
  u784=9.96048508602880922d-1*u948
  u609=u2855*(u882+u763)
  u870=(u609+u638)*u2856
  value(9)=value(9)+(c(37)*(y*(u2856*(u521+u2855*(u2855*(u784+u675)+u83&
 &6)+u870)+u2855*(u445+u2855*(u4145+u585))+u3967)))
  u4145=-2.98814552580864277d0*u376
  u3998=+u4145
  u607=3.73518190726080346d0*u3808
  u3810=-8.9644365774259283d-1*u948
  value(1)=value(1)+(c(38)*(u2983*((u595+u2855*(u4030+u681))*u2856+u285&
 &5*(u607+u2855*(u4061+u3810))+u3998)*z))
  u3998=u2858*(4.45320758646353191d4*u534+u3788)*u4147
  u3810=3.16940694670729944d-1*u3998
  u866=-u383
  u499=-1.05646898223576648d-1*u3946
  u765=1.05646898223576648d-1*u3808
  u971=1.63284278170329504d5*u534
  u3970=u2858*(u971-5.67185232289765129d1*u575)*u4147
  u575=-2.53552555736583955d0*u3970
  u514=+u575
  u3980=-u731
  u731=1.05646898223576648d-1*u952
  u389=-2.11293796447153296d-1*u948
  u969=u2858*(u858+u3788*u913)*u4147
  u913=1.05646898223576648d-1*u969
  u412=-u529
  u55=-3.5215632741192216d-2*u701
  u3798=u2855*(u389+u3904)
  u3916=u389*u2855
  value(2)=value(2)+(c(38)*(y*((u866+u2855*(u2855*(u412+u3975)+u3980))*&
 &u2857+u2856*(u499+u2855*(u2855*(u55+u3904)+u731)+(u3798+u765)*u2856)+&
 &u2855*(u514+u2855*(u3916+u913))+u3810)))
  u3810=-2.25882569787428772d0*u376
  u514=1.69411927340571579d0*u3808
  u913=-2.82353212234285965d-1*u3808
  u3979=8.47059636702857895d-1*u3808
  u941=-1.12941284893714386d0*u948
  u413=+u941
  u684=1.8823547482285731d-1*u948
  u633=1.12941284893714386d-1*u948
  value(3)=value(3)+(c(38)*(u528*((u514+u2855*(u3828+u413))*u2857+(u913&
 &+u2855*(u4025+u684))*u2856+u2855*(u3979+u2855*(u4025+u633))+u3810)))
  u639=3.59376967909202506d-1*u2858
  u3949=u639*(2.86984488905427612d4*u534+u3788)*u4147
  u4203=-3.59376967909202506d-1*u376
  u4067=1.59723096848534447d-1*u3808
  u997=-1.19792322636400836d-1*u3946
  u516=1.19792322636400836d-1*u3808
  u378=5.55166545779120312d5*u534
  u519=-4.79169290545603342d-1*u2858*(u378-1.13437046457953026d2*u606)*&
 &u4147
  u3843=-3.19446193697068895d-1*u948
  u606=1.19792322636400836d-1*u952
  u919=-2.39584645272801671d-1*u948
  u685=u3876*(u953+u3788*(u844+u4009))*u4147
  u634=2.79515419484935283d-1*u948
  u3862=-3.99307742121336118d-2*u701
  u947=3.99307742121336118d-2*u700
  u4193=-7.98615484242672237d-2*u948
  u676=u947*u2855
  u875=u2855*(u3862+u676)
  u959=(y*(u2857*(u4203+u2855*(u2855*(u634+u430)+u516)+(u2855*(u3843+u3&
 &767)+u4067)*u2857)+u2856*(u997+u2855*(u875+u606)+(u2855*(u919+u676)+u&
 &516)*u2856)+u2855*(u519+u2855*(u4193*u2855+u685))+u3949))
  value(4)=value(4)+(c(38)*u959)
  u533=9.09159805758655719d-1*u2858*(1.28648219164502033d4*u534+u3827)*&
 &u4147
  u3944=u2858*(u971-2.83592616144882564d1*u3986)*u4147
  u3986=-4.04071024781624764d-1*u3944
  u3951=2.02035512390812382d-1*u3808
  u3864=-3.03053268586218573d-1*u376
  u988=1.89408292866386608d-2*u3808
  u998=+u990*(u610-1.13437046457953026d2*u84)*u4147
  u610=2.12269561621428355d6*u534
  u84=u2858*(u610+u3827*u785)*u4147
  u785=8.08142049563249528d-1*u84
  u468=-4.04071024781624764d-1*u948
  u91=7.95514830038823754d-1*u3808
  u3835=-3.78816585732773216d-2*u948
  u786=u3838*(1.71089266666871254d8*u534+u3788*(4.03d2*u4031+u977))*u41&
 &47
  u3978=5.0d0*u526
  u72=+u3879*(u49+u3788*(u3978+u604))*u4147
  u4137=6.73451707969374607d-2*u700
  u469=-2.2728995143966393d-1*u948
  u470=-1.89408292866386608d-1*u948
  u4135=1.26272195244257739d-2*u700
  u3917=u4137*u2855
  u3918=u4135*u2855
  value(5)=value(5)+(c(38)*(z*(u2857*(u3986+u2855*(u72*u2855+u785)+(u28&
 &55*(u468+u3917)+u3951)*u2857)+u2856*(u3864+u2855*(u2855*(u469+u3918)+&
 &u91)+(u2855*(u3835+u3896)+u988)*u2856)+u2855*(u998+u2855*(u2855*(u470&
 &+u3896)+u786))+u533)))
  u650=5.98961613182004178d-1*u7
  u429=-5.98961613182004178d-1*u376
  u3778=5.98961613182004178d-1*u3808
  u3889=-3.99307742121336118d-1*u2858
  u4134=+u3889*(u92-1.13437046457953026d2*u962)*u4147
  u92=-9.98269355303340296d-1*u3808
  u962=u2858*(1.40097910670142714d7*u534+u3788*(u732+u3891))*u4147
  u732=1.99653871060668059d-1*u962
  u414=-3.99307742121336118d-1*u948
  u787=u3876*(u788+u3788*(u770+u4006))*u4147
  u788=5.19100064757736954d-1*u948
  u770=1.71938344913356967d8*u534
  u893=2.7d1*u526
  u43=u2858*(u770+u3788*(u893+u645))*u4147
  u56=-3.99307742121336118d-2*u43
  u4057=u56+u676
  u538=u2855*(u4057)
  value(6)=value(6)+(c(38)*(x*(u2857*(u429+u2855*(u2855*(u788+u430)+u92&
 &)+(u2855*(u536+u3767)+u69)*u2857)+u2856*(u704+u2855*(u538+u732)+(u285&
 &5*(u414+u676)+u3778)*u2856)+u2855*(u4134+u787*u2855)+u650)))
  u706=u2858*(3.92542002066044665d4*u534+u3788)*u4147
  u377=5.08235782021714737d-1*u706
  u4170=-6.77647709362286316d-1*u3970
  u4183=1.69411927340571579d-1*u3808
  u3815=-1.69411927340571579d-1*u376
  u390=2.82353212234285965d-2*u3808
  u388=u2858*(9.56845870078130891d6*u534+u3788*(2.93d2+u3901))*u4147
  u760=-5.6470642446857193d-2*u388
  u96=1.06134780810714177d7*u534
  u3859=2.5d1*u4031
  u955=u2858*(u96+u3788*(u3859+u4006))*u4147
  u64=3.38823854681143158d-1*u955
  u3871=-3.38823854681143158d-1*u948
  u4072=3.38823854681143158d-1*u3808
  u4200=-5.6470642446857193d-2*u948
  u457=u2858*(1.99533387924142653d7*u534+u3788*(4.7d1*u4031+u3891))*u41&
 &47
  u3839=8.47059636702857895d-2*u457
  u4144=u2858*(u49+u3788*(u3978+u4074))*u4147
  u4167=-4.51765139574857544d-1*u4144
  u708=9.41177374114286549d-3*u700
  u946=-7.52941899291429239d-2*u948
  u3881=u2855*(u3871+u4038)
  u3769=u708*u2855
  u63=u2855*(u4200+u3769)
  u3852=u4200*u2855
  value(7)=value(7)+(c(38)*(z*(u2857*(u4170+u2855*(u2855*(u4167+u3828)+&
 &u64)+(u3881+u4183)*u2857)+u2856*(u3815+u2855*(u3852+u4072)+(u63+u390)&
 &*u2856)+u2855*(u760+u2855*(u2855*(u946+u4065)+u3839))+u377)))
  u377=u2858*(3.66152623775890402d4*u534+u3788)*u4147
  u4170=1.58470347335364972d0*u377
  u4183=-3.69764143782518268d0*u376
  u3815=+u4183
  u760=-u423
  u64=u2858*(u378-1.13437046457953026d2*u3942)*u4147
  u3839=-1.05646898223576648d0*u64
  u4167=+u3839
  u946=3.69764143782518268d0*u3808
  u378=5.28234491117883239d-1*u962
  u771=-u561
  u423=1.3160712820528558d7*u534
  u3942=3.1d1*u4031
  u743=u2858*(u423+u3788*(u3942+u4006))*u4147
  u3982=1.05646898223576648d-1*u743
  u4217=-7.39528287565036535d-1*u948
  u789=-1.05646898223576648d-1*u43
  u471=-1.40862530964768864d-1*u948
  value(8)=value(8)+(c(38)*(x*((u3815+u2855*(u2855*(u4217+u3904)+u946))&
 &*u2857+u2856*(u615+u2855*(u2855*(u789+u3975)+u378)+(u2855*(u771+u3975&
 &)+u760)*u2856)+u2855*(u4167+u2855*(u471*u2855+u3982))+u4170)))
  u3815=8.87122154999207171d2*u702
  u4167=8.9644365774259283d-1*u376
  u378=-2.09170186806604993d0*u376
  u3982=-2.9134418876634267d0*u3808
  u4217=1.53142458197692942d0*u3808
  u789=1.04585093403302497d0*u948
  u471=-2.73913339865792253d-1*u948
  value(9)=value(9)+(c(38)*((u2856*(u4167+u2855*(u2855*(u789+u675)+u398&
 &2)+u870)+u2855*(u378+u2855*(u2855*(u471+u763)+u4217))+u3815)*z))
  u870=-1.49407276290432138d-1*u3808
  u492=-2.24110914435648207d0*u376
  u4116=+u492
  u635=2.98814552580864276d-1*u948
  u637=1.49407276290432138d0*u3808
  u862=-2.49012127150720231d-1*u948
  u872=-1.49407276290432138d-1*u948
  u588=u2855*(u635+u4030)
  u4172=(u588+u870)*u2856
  u4084=u872*u2855
  value(1)=value(1)+(c(39)*(y*(u2856*(u521+u2855*(u2855*(u862+u4061)+u5&
 &95)+u4172)+u2855*(u4116+u2855*(u4084+u637)))))
  u4116=-3.16940694670729944d0*u376
  u862=+u4116
  u614=1.05646898223576648d-1*u948
  u438=2.6411724555894162d0*u3808
  u537=-7.0431265482384432d-1*u948
  u750=-3.16940694670729944d-1*u948
  value(2)=value(2)+(c(39)*(u2983*(u2856*(u852+u2855*(u3975+u537)+(u402&
 &4+u614)*u2856)+u2855*(u438+u750*u2855)+u862)*z))
  u537=-1.01647156404342947d0*u640
  u4097=1.69411927340571579d-1*u524
  u812=-3.95294497128000351d-1*u3808
  u4051=2.4d1*pd2
  u3764=(1.75d2+u4051)
  u4098=u594*(5.71494973596153262d6*u534+u3788*u3764)*u4147
  u593=4.8d1*u4011
  u472=+u3886*(1.04012085194499894d8*u534+u3788*(2.45d2*u4031+u593))*u4&
 &147
  u733=7.90588994256000702d-1*u948
  u3761=u2858*(u858+u3788*u473)*u4147
  u473=-2.25882569787428772d-1*u3761
  u3845=4.8d1*u4074
  u746=u4042*(1.11441519851249886d9*u534+u3788*(1.75d2*u526+u3845))*u41&
 &47
  u646=3.95294497128000351d-1*u948
  u567=u2855*(u746+u3860)
  u3919=u646*u2855
  value(3)=value(3)+(c(39)*(y*(u2856*(u4097+u2855*(u567+u472)+(u2855*(u&
 &733+u3860)+u812)*u2856)+u2855*(u4098+u2855*(u3919+u473))+u537)))
  u648=2.39584645272801671d-1*u524
  u963=+u3885*(u642+u3788*(1.05d2*u4031+u3786))*u4147
  u642=2.22883039702499772d8*u534
  u85=3.5d1*u526
  u915=u2858*(u642+u3788*(u85+u645))*u4147
  u4100=5.32410322828448158d-2*u915
  u44=-9.31718064949784276d-2*u700
  u4085=u44*u2855
  u3929=(u2983*(u2856*(u963+u2855*(u4085+u4100)+(u4085+u634)*u2856)+u28&
 &55*(u963+u634*u2855)+u648)*z)
  value(4)=value(4)+(c(39)*u3929)
  u3818=-7.08981540362206411d0*u382
  u738=1.69709830408282401d1*u2858*(2.12057504117311044d3*u534+u3818)*u&
 &4147
  u4076=(2.5d1+u4019)
  u897=u2858*(u4118+u3788*u4076)*u4147
  u840=1.51526634293109286d-1*u897
  u515=u2858*(u610+u3788*u790)*u4147
  u790=-2.02035512390812382d-1*u515
  u686=5.05088780977030955d-2*u948
  u823=+u4062*(1.02869095247307587d7*u534+u3788*(3.15d2+4.4d1*pd2))*u41&
 &47
  u791=u3795*(u815+u3788*(2.25d2*u4031+u4035))*u4147
  u3863=-1.70467463579747947d-1*u948
  u830=4.08210695425823759d6*u534
  u986=+u4062*(u830+u3788*(1.25d2+u3888))*u4147
  u474=-5.05088780977030955d-2*u952
  u616=1.68362926992343652d-2*u701
  u39=-1.68362926992343652d-2*u700
  u427=5.15815034740070902d8*u534
  u4012=1.215d3*u4031
  u792=u4032*(u427+u3788*(u4012+2.56d2*u4011))*u4147
  u4085=6.4d1*u4074
  u73=+u608*(2.00594735732249795d9*u534+u3788*(3.15d2*u526+u4085))*u414&
 &7
  u793=u3838*(9.97666939620713267d7*u534+u3788*(2.35d2*u4031+u4035))*u4&
 &147
  u794=1.01017756195406191d-1*u948
  u74=+u608*(1.89132179404692664d9*u534+u3788*(2.97d2*u526+u4085))*u414&
 &7
  u4101=1.13644975719831965d-1*u700
  u475=-5.68224878599159824d-2*u948
  u3920=u39*u2855
  u4086=u4101*u2855
  value(5)=value(5)+(c(39)*(x*(u2857*(u840+u2855*(u794*u2855+u474)+u285&
 &7*((u686+u3920)*u2857+u2855*(u616+u3920)+u790))+u2856*(u823+u2855*(u2&
 &855*(u74+u444)+u792)+u2856*((u3863+u444)*u2856+u2855*(u73+u4086)+u791&
 &))+u2855*(u986+u2855*(u475*u2855+u793))+u738)))
  value(6)=value(6)+(c(39)*u737)
  u737=-3.38823854681143158d-1*u640
  u848=u2858*(u88-1.13437046457953026d2*u4219)*u4147
  u4219=-5.6470642446857193d-2*u848
  u4096=8.47059636702857895d-2*u952
  u415=-1.97647248564000175d-1*u948
  u951=3.95294497128000351d-1*u3946
  u795=-1.12941284893714386d-1*u952
  u3858=u2858*(u642+u3788*(u85+u727))*u4147
  u450=-2.82353212234285965d-2*u3858
  u796=-2.82353212234285965d-2*u968
  u936=3.12036255583499681d8*u534
  u853=4.9d1*u526
  u4099=u2858*(u936+u3788*(u853+u727))*u4147
  u647=2.82353212234285965d-2*u4099
  u3784=6.58824161880000585d-2*u948
  u40=-6.58824161880000585d-2*u700
  u3802=u40*u2855
  value(7)=value(7)+(c(39)*(x*(u2856*(u4219+u2855*(u2855*(u647+u3802)+u&
 &795)+u2856*((u415+u3801)*u2856+u450*u2855+u4096))+u2855*(u951+u2855*(&
 &u3784*u2855+u796))+u737)))
  u737=-1.25458018308152576d3*u702
  u4219=+u737
  u4096=1.37340967690649642d0*u376
  u951=-3.16940694670729944d-1*u3808
  u795=3.16940694670729944d-1*u376
  u450=-2.21858486269510961d0*u3808
  u796=+u450
  u647=6.33881389341459887d-1*u948
  u3784=2.11293796447153296d-1*u3808
  u416=-3.5215632741192216d-2*u948
  u3861=u2855*(u647+u3898)
  value(8)=value(8)+(c(39)*((u2856*(u4096+u2855*(u2855*(u614+u3904)+u79&
 &6)+(u3861+u951)*u2856)+u2855*(u795+u2855*(u416*u2855+u3784))+u4219)*z&
 &))
  u4219=-3.73518190726080346d-2*u948
  u795=5.60277286089120519d-1*u948
  u796=1.86759095363040173d-1*u948
  u816=-1.24506063575360115d-2*u948
  u4087=u816*u2855
  value(9)=value(9)+(c(39)*(x*(u2856*(u445+u2855*(u2855*(u796+u763)+u83&
 &6)+u2856*((u4219+u763)*u2856+u2855*(u795+u675)+u585))+u2855*(u521+u28&
 &55*(u4087+u638))+u3967)))
  u508=1.49407276290432138d-1*u376
  u531=-1.34466548661388924d0*u376
  u4169=+u531
  u3947=4.48221828871296415d-1*u3808
  u987=5.97629105161728553d-1*u3808
  u417=-4.48221828871296415d-1*u948
  value(1)=value(1)+(c(40)*((u2856*(u508+u2855*(u2855*(u417+u4061)+u394&
 &7)+u4172)+u2855*(u4169+u2855*(u4082+u987))+u3815)*z))
  u508=u2858*(1.48440252882117731d4*u534+u3827)*u4147
  u4169=-2.53552555736583955d0*u508
  u987=+u4169
  u4172=2.53552555736583955d0*u3944
  u3756=u2858*(u610+u3788*(u3983+u4006))*u4147
  u3983=1.05646898223576648d-1*u3756
  u3820=2.75950430107856861d7*u534
  u837=u2858*(u3820+u3788*(u881+u3786))*u4147
  u4136=-2.11293796447153296d-1*u837
  u509=-u3788*(u645-5.0d0*u3894)
  u960=3.5215632741192216d-2*u2858
  u581=u960*(u49+u509)*u4147
  u970=u2858*(9.76439983458570431d6*u534+u3788*u476)*u4147
  u476=-1.05646898223576648d-1*u970
  u839=1.05646898223576648d-1*u701
  u75=-7.04312654823844319d-2*u700
  u3993=u2855*(u839+u3898)
  u3853=u614*u2855
  u4088=u75*u2855
  value(2)=value(2)+(c(40)*(x*(u2856*(u4172+u2855*(u3993+u4136)+u2856*(&
 &(u3971+u3904)*u2856+u2855*(u581+u4088)+u3983))+u2855*(u4172+u2855*(u3&
 &853+u476))+u987)))
  u987=-u931*(1.55037597454656296d4*u534+u3827)*u4147
  u3983=3.38823854681143158d-1*u3946
  u4136=-3.38823854681143158d-1*u3808
  u581=6.21177066915429123d-1*u376
  u931=-5.6470642446857193d-2*u3808
  u709=u594*(5.12712633454834641d6*u534+u3788*(1.57d2+u4051))*u4147
  u4192=-3.38823854681143158d-1*u952
  u687=6.77647709362286316d-1*u948
  u4078=-1.52470734606514421d0*u3808
  u958=u2858*(8.9153215880999909d6*u534+u3788*u477)*u4147
  u477=-3.38823854681143158d-1*u958
  u758=1.12941284893714386d-1*u701
  u641=-1.12941284893714386d-1*u700
  u797=3.57647402163428889d-1*u948
  u3854=u641*u2855
  u4119=u2855*(u758+u3854)
  u3952=u2855*(u646+u4025)
  value(3)=value(3)+(c(40)*(z*(u2857*(u3983+u2855*(u4119+u4192)+(u2855*&
 &(u687+u3854)+u4136)*u2857)+u2856*(u581+u2855*(u3952+u4078)+(u2855*(u6&
 &33+u4025)+u931)*u2856)+u2855*(u709+u2855*(u797*u2855+u477))+u987)))
  u4121=9.58338581091206684d-1*u508
  u4148=4.79169290545603342d-1*u897
  u688=-6.38892387394137789d-1*u515
  u689=1.59723096848534447d-1*u948
  u4082=1.1d1*pd2
  u818=u2858*(u88+u3827*(7.0d1+u4082))*u4147
  u453=-4.79169290545603342d-1*u818
  u690=1.19792322636400836d-1*u955
  u744=-1.19792322636400836d-1*u948
  u4053=-9.58338581091206684d-1*u3944
  u562=5.73127816377856558d7*u534
  u82=u2858*(u562+u3788*(1.35d2*u4031+u4035))*u4147
  u691=7.98615484242672237d-2*u82
  u3841=-3.99307742121336118d-2*u915
  u692=3.99307742121336118d-2*u970
  u45=-3.99307742121336118d-2*u517
  u984=7.98615484242672237d-2*u700
  u4026=-3.99307742121336118d-2*u948
  u4089=u984*u2855
  u814=(x*(u2857*(u4148+u2855*(u631*u2855+u4197)+u2857*((u689+u3851)*u2&
 &857+u3809+u688))+u2856*(u453+u2855*(u2855*(u45+u676)+u691)+u2856*((u7&
 &44+u676)*u2856+u2855*(u3841+u4089)+u690))+u2855*(u4053+u2855*(u4026*u&
 &2855+u692))+u4121))
  value(4)=value(4)+(c(40)*u814)
  u3809=u2858*(u971+u3788*(5.0d0+pd2))*u4147
  u3940=-9.09159805758655719d-1*u3809
  u506=-2.52544390488515477d-1*u3808
  u418=-5.05088780977030955d-2*u948
  u798=7.57633171465546432d-2*u952
  u799=2.52544390488515477d-2*u981
  u800=3.36725853984687303d-1*u948
  u3844=3.2d1*u4074
  u76=+u3875*(6.04968250621070811d8*u534+u3788*(9.5d1*u526+u3844))*u414&
 &7
  value(5)=value(5)+(c(40)*(u528*(u2857*(u506+u2855*(u3768+u800)+(u3915&
 &+u418)*u2857)+u2856*(u798+u2855*(u3800+u76)+(u444+u3863)*u2856)+u2855&
 &*(u799+u2855*(u3896+u3863))+u3940)))
  value(6)=value(6)+(c(40)*u959)
  u959=u2858*(u830+u3788*(1.25d2+u3790))*u4147
  u644=-5.6470642446857193d-2*u959
  u982=u2858*(u96+u3788*(u3859+u3833))*u4147
  u93=1.12941284893714386d-1*u982
  u385=6.77647709362286316d-1*u515
  u510=-5.6470642446857193d-1*u948
  u603=2.4d1*u4074
  u842=u2858*(4.13925645161785292d8*u534+u3788*(6.5d1*u526+u603))*u4147
  u3966=-1.8823547482285731d-2*u842
  u478=-1.41176606117142983d-1*u948
  value(7)=value(7)+(c(40)*(u528*((u3979+u2855*(u4038+u510))*u2857+u285&
 &6*(u93+u2855*(u4038+u3966)+(u3801+u415)*u2856)+u2855*(u385+u2855*(u40&
 &65+u478))+u644)))
  u644=9.50822084012189831d-1*u4114
  u93=-5.28234491117883239d-1*u376
  u385=-3.16940694670729944d-1*u3946
  u510=3.16940694670729944d-1*u3808
  u3966=-5.0710511147316791d0*u3944
  u802=+u3966
  u734=3.16940694670729944d-1*u952
  u613=-6.33881389341459887d-1*u948
  u851=u2858*(u3820+u3788*u801)*u4147
  u801=1.05646898223576648d-1*u851
  u479=-5.28234491117883239d-1*u948
  u46=-1.05646898223576648d-1*u701
  u99=-3.5215632741192216d-1*u948
  u50=u2855*(u613+u3975)
  u907=(u50+u510)
  u3921=u99*u2855
  value(8)=value(8)+(c(40)*(y*((u93+u2855*(u2855*(u479+u3904)+u760))*u2&
 &857+u2856*(u385+u2855*(u2855*(u46+u3975)+u734)+u907*u2856)+u2855*(u80&
 &2+u2855*(u3921+u801))+u644)))
  u644=1.49407276290432138d0*u376
  u802=6.22530317876800577d-1*u948
  u97=(u763+u4219)
  value(9)=value(9)+(c(40)*(u2983*(u2856*(u595+u2855*(u675+u802)+u97*u2&
 &856)+u2855*(u595+u2855*u97)+u644)*z))
  u97=4.48221828871296415d-1*u640
  u985=-1.49407276290432138d-1*u3946
  u841=1.49407276290432138d-1*u3808
  u954=-8.9644365774259283d-1*u3946
  u735=1.49407276290432138d-1*u952
  u677=-2.98814552580864276d-1*u948
  u803=1.49407276290432138d-1*u955
  u57=-4.98024254301440461d-2*u701
  u480=-1.99209701720576184d-1*u948
  u3763=u2855*(u677+u4061)
  value(1)=value(1)+(c(41)*(y*((u89)*u2857+u2856*(u985+u2855*(u2855*(u5&
 &7+u4061)+u735)+(u3763+u841)*u2856)+u2855*(u954+u2855*(u480*u2855+u803&
 &))+u97)))
  u97=u2858*(u4118+u3788*u4199)*u4147
  u954=-3.16940694670729944d-1*u97
  u4199=1.91042605459285519d7*u534
  u89=4.5d1*u4031
  u3948=u2858*(u4199+u3788*(u89+u4006))*u4147
  u989=1.05646898223576648d-1*u3948
  u857=u2858*(u49+u3788*(u3978+u592))*u4147
  u3978=-1.40862530964768864d-1*u857
  u3922=u630*u2855
  value(2)=value(2)+(c(41)*(u528*((u760+u2855*(u3975+u771))*u2857+u2856&
 &*(u614+u2855*(u3904+u3978)+(u3904+u3971)*u2856)+u2855*(u989+u3922)+u9&
 &54)))
  u954=u594*(4.25528724928737494d4*u534+u3788)*u4147
  u989=-1.75058991585257298d0*u376
  u3978=-5.6470642446857193d-2*u3946
  u625=5.6470642446857193d-2*u3808
  u964=u2858*(u932-1.41796308072441282d1*u854)*u4147
  u932=-5.42118167489829052d0*u964
  u854=3.21882661947086d0*u3808
  u481=-6.77647709362286316d-1*u948
  u666=5.6470642446857193d-2*u952
  u874=-1.12941284893714386d-1*u948
  u804=5.6470642446857193d-2*u958
  u41=-1.8823547482285731d-2*u701
  u523=1.8823547482285731d-2*u700
  u419=-3.7647094964571462d-2*u948
  u3973=u523*u2855
  u62=u2855*(u41+u3973)
  value(3)=value(3)+(c(41)*(y*(u2857*(u989+u2855*(u2855*(u411+u4025)+u8&
 &54)+(u2855*(u481+u3828)+u4072)*u2857)+u2856*(u3978+u2855*(u62+u666)+(&
 &u2855*(u874+u3973)+u625)*u2856)+u2855*(u932+u2855*(u419*u2855+u804))+&
 &u954)))
  u4138=-1.19792322636400836d-1*u848
  u4141=9.98269355303340296d-1*u3808
  u3989=-1.59723096848534447d-1*u948
  u693=1.19792322636400836d-1*u948
  u658=u3876*(u858+u3788*(u47+u4009))*u4147
  u3914=-1.3310258070711204d-1*u948
  u47=-1.59723096848534447d-1*u857
  u978=u676+u47
  u4090=u658*u2855
  u90=(u528*(u2857*(u4141+u2855*(u430+u3914)+(u3767+u3989)*u2857)+u2856&
 &*(u693+u2855*(u978)+(u676+u744)*u2856)+u4090+u4138))
  value(4)=value(4)+(c(41)*u90)
  u898=u3795*(1.9297232874675305d5*u534-1.24780751103748328d3*u382)*u41&
 &47
  u4069=(5.15d2+9.6d1*pd2)
  u4189=+u4045*(1.68182806515439389d7*u534+u3788*u4069)*u4147
  u551=6.06106537172437146d-1*u515
  u482=-2.02035512390812382d-1*u948
  u3991=u3795*(6.69465540498350965d6*u534+u3788*(2.05d2+u821))*u4147
  u4202=u2858*(u483+u3788*u636)*u4147
  u483=-3.78816585732773216d-2*u4202
  u636=1.89408292866386608d-2*u948
  u387=u4032*(u4118+1.13437046457953026d2*u80)*u4147
  u80=8.27851290323570583d7*u534
  u458=6.4d1*u4011
  u805=u3795*(u80+u3788*(1.95d2*u4031+u458))*u4147
  u77=-2.02035512390812382d-1*u857
  u564=1.84674518610642669d8*u534
  u86=1.04d2*u4011
  u484=+u608*(u564+u3788*(4.35d2*u4031+u86))*u4147
  u820=2.86563908188928279d8*u534
  u3803=4.5d1*u526
  u861=u3838*(u820+u3788*(u3803+u645))*u4147
  u905=-6.31360976221288694d-3*u700
  u485=-1.26272195244257739d-2*u948
  u486=-1.32585805006470626d-1*u948
  u428=2.52544390488515477d-2*u4194
  u48=-1.89408292866386608d-2*u700
  u694=6.31360976221288694d-3*u948
  u753=-1.26272195244257739d-2*u700
  u3770=u753*u2855
  u724=u2855*(u428+u3770)
  u3829=u905*u2855
  u3923=u48*u2855
  u4091=u694*u2855
  value(5)=value(5)+(c(41)*(x*(u2857*(u4189+u2855*(u2855*(u486+u3896)+u&
 &805)+u2857*((u482+u3917)*u2857+u77*u2855+u551))+u2856*(u3991+u2855*(u&
 &724+u484)+u2856*((u636+u3829)*u2856+u2855*(u861+u3923)+u483))+u2855*(&
 &u387+u2855*(u4091+u485))+u898)))
  u4146=3.59376967909202506d-1*u640
  u3915=-7.98615484242672237d-2*u2858
  u449=+u3915*(u88-1.13437046457953026d2*u4140)*u4147
  u4028=2.79515419484935283d-1*u3808
  u4140=-2.39584645272801671d-1*u848
  u806=1.59723096848534447d-1*u3761
  u3761=-5.59030838969870566d-1*u948
  u440=1.2d1*u4074
  u498=u2858*(u642+u3788*(u85+u440))*u4147
  u78=-2.66205161414224079d-2*u498
  u441=9.31718064949784276d-2*u700
  u3924=u441*u2855
  value(6)=value(6)+(c(41)*(z*(u2857*(u449+u2855*(u78*u2855+u806)+(u285&
 &5*(u3761+u3924)+u4028)*u2857)+u2855*(u4140+u4090)+u4146)))
  u559=u2858*(4.05736691211121797d4*u534+u3788)*u4147
  u835=4.23529818351428947d-1*u559
  u3912=7.2d1*pd2
  u4090=(5.95d2+u3912)
  u3785=u2858*(1.94308291022692109d7*u534+u3788*u4090)*u4147
  u856=-2.82353212234285965d-2*u3785
  u3825=1.69411927340571579d-1*u955
  u674=-1.69411927340571579d-1*u948
  u502=8.47059636702857895d-2*u897
  u807=-1.12941284893714386d-1*u515
  u695=2.82353212234285965d-2*u948
  u6=-5.6470642446857193d-1*u3970
  u809=1.25239041356642729d8*u534
  u885=2.95d2*u4031
  u4204=u2858*(u809+u3788*(u885+u593))*u4147
  u965=2.82353212234285965d-2*u4204
  u659=-5.6470642446857193d-2*u915
  u487=-2.82353212234285965d-2*u952
  u832=9.41177374114286549d-3*u701
  u3941=1.12941284893714386d-1*u3808
  u4142=-4.8000046079828614d-1*u948
  u624=5.6470642446857193d-2*u948
  u736=9.41177374114286549d-3*u948
  u3855=u624*u2855
  u4092=u736*u2855
  value(7)=value(7)+(c(41)*(x*(u2857*(u856+u2855*(u2855*(u4142+u4065)+u&
 &965)+u2857*((u674+u4038)*u2857+u2855*(u659+u3828)+u3825))+u2856*(u502&
 &+u2855*(u3855+u487)+u2856*((u695+u4065)*u2856+u2855*(u832+u4065)+u807&
 &))+u2855*(u6+u2855*(u4092+u3941))+u835)))
  u835=9.50822084012189831d-1*u640
  u856=-2.11293796447153296d-1*u703
  u3825=-1.90164416802437966d0*u3946
  u6=+u3825
  u965=2.53552555736583955d0*u515
  u659=3.16940694670729944d-1*u955
  u4142=-2.11293796447153296d-1*u4194
  u845=1.40862530964768864d-1*u700
  u3925=u845*u2855
  value(8)=value(8)+(c(41)*(z*(u2857*(u856+u2855*(u2855*(u4142+u3925)+u&
 &965)+u907*u2857)+u2855*(u6+u2855*(u3922+u659))+u835)))
  u856=1.00849911496041693d0*u2858*(3.35365015770710428d4*u534+u3788)*u&
 &4147
  u6=-1.30731366754128121d0*u376
  u659=u2858*(3.42896984157691957d6*u534+u3788*(1.05d2+u3888))*u4147
  u4142=-1.12055457217824104d-1*u659
  u907=-5.0d0*u4001
  u957=-u3788*(u3832+u907)
  u869=7.47036381452160691d-2*u2858
  u808=u869*(u610+u957)*u4147
  u696=3.73518190726080346d-2*u948
  u4213=1.9d1*pd2
  u4139=+u3900*(4.7352440669395556d6*u534+u3788*(1.45d2+u4213))*u4147
  u525=1.30731366754128121d0*u3808
  u934=7.2d1*u4011
  u983=u2858*(u809+u3788*(u885+u934))*u4147
  u809=3.73518190726080346d-2*u983
  u885=-u3788*(u645-3.5d1*u3894)
  u79=+u497*(u885+u642)*u4147
  u530=-1.24506063575360115d-2*u700
  u810=u869*(1.14625563275571312d7*u534+u3788*u4218)*u4147
  u4218=-2.61462733508256242d-1*u948
  u869=u2858*(u80+u3788*(1.3d1*u526+u3766))*u4147
  u80=-1.49407276290432138d-1*u869
  u850=6.22530317876800576d-2*u700
  u488=-8.71542445027520806d-2*u948
  u4179=7.47036381452160691d-2*u700
  u756=u4179*u2855
  u3771=u530*u2855
  u3856=u850*u2855
  value(9)=value(9)+(c(41)*(x*((u6+u2855*(u2855*(u4218+u763)+u525))*u28&
 &57+u2856*(u4142+u2855*(u2855*(u80+u756)+u809)+u2856*((u696+u3771)*u28&
 &56+u2855*(u79+u3856)+u808))+u2855*(u4139+u2855*(u488*u2855+u810))+u85&
 &6)))
  u3927=-u492
  u492=-u637
  u637=1.49407276290432138d-1*u948
  u4008=7.47036381452160692d-1*u3808
  u817=2.49012127150720231d-1*u948
  u3870=(u637+u4030)
  u3857=u677*u2855
  u3926=u841*u2855
  value(1)=value(1)+(c(42)*(x*(u2856*(u3927+u2855*(u3857+u4008)+u2856*(&
 &u3870*u2856+u2855*(u817+u4061)+u492))+u2855*(u767+u3926))))
  u3927=-u737
  u492=-3.16940694670729944d-1*u376
  u817=-2.11293796447153296d-1*u3808
  u737=3.5215632741192216d-2*u948
  u4095=-u4096
  u4096=-u450
  u450=(u3971+u3975)
  u4093=u613*u2855
  u4094=u510*u2855
  value(2)=value(2)+(c(42)*((u2856*(u492+u2855*(u4093+u4096)+u2856*((u7&
 &37+u4024)*u2856+u2855*u450+u817))+u2855*(u4095+u4094)+u3927)*z))
  value(3)=value(3)+(c(42)*(x*(u2856*(u4098+u2855*(u733*u2855+u472)+u28&
 &56*((u646+u3860)*u2856+u567+u473))+u2855*(u4097+u812*u2855)+u537)))
  u3927=u504*u2856
  u567=u2856*(u843+u3927)
  u4095=u632*u2855
  u4096=u3977*u2855
  u472=(z*(u2857*(u503+u2856*(u567+u4197)+(u2856*(u631+u3927)+u3938)*u2&
 &857)+u2856*(u3834+u2855*(u4095+u3939)+u2856*((u683+u430)*u2856+u3928+&
 &u726))+u2855*(u3847+u4096)+u757))
  value(4)=value(4)+(c(42)*u472)
  u4097=u39*u2856
  u4098=u3863*u2855
  value(5)=value(5)+(c(42)*(y*(u2857*(u840+u2856*(u794*u2856+u474)+u285&
 &7*((u686+u4097)*u2857+u2856*(u616+u4097)+u790))+u2856*(u986+u2855*(u2&
 &855*(u73+u444)+u792)+u2856*((u475+u444)*u2856+u2855*(u74+u4086)+u793)&
 &)+u2855*(u823+u2855*(u4098+u791))+u738)))
  value(6)=value(6)+(c(42)*u3929)
  u407=3.38823854681143158d-1*u640
  u986=-3.95294497128000351d-1*u3946
  u648=2.82353212234285965d-2*u968
  u840=-6.58824161880000585d-2*u948
  u746=5.6470642446857193d-2*u848
  u4100=1.12941284893714386d-1*u952
  u616=-2.82353212234285965d-2*u4099
  u4101=-8.47059636702857895d-2*u952
  u537=2.82353212234285965d-2*u3858
  u738=1.97647248564000175d-1*u948
  u4099=u738*u2855
  value(7)=value(7)+(c(42)*(y*(u2856*(u986+u2855*(u2855*(u537+u3802)+u4&
 &100)+u2856*((u840+u3801)*u2856+u616*u2855+u648))+u2855*(u746+u2855*(u&
 &4099+u4101))+u407)))
  u407=-u4116
  u746=-u438
  u648=3.16940694670729944d-1*u948
  u4100=7.0431265482384432d-1*u948
  value(8)=value(8)+(c(42)*(u2983*(u2856*(u746+u2855*(u3904+u4100)+(u38&
 &98+u648)*u2856)+u2855*(u863+u4083)+u407)*z))
  u746=(u816+u763)*u2856
  u4100=u4219*u2855
  value(9)=value(9)+(c(42)*(y*(u2856*(u521+u2855*(u2855*(u795+u763)+u83&
 &6)+u2856*(u746+u2855*(u796+u675)+u638))+u2855*(u445+u2855*(u4100+u585&
 &))+u3967)))
  value(1)=value(1)+(c(43)*(u2983*(u2856*(u595+u3782*u999+(u4030+u637)*&
 &u2856)+u2855*(u4008+u4084))*z))
  u616=-8.45175185788613183d-1*u508
  u4101=u3792*pd*pd2*u4147
  u537=7.98951473461766613d0*u4101
  u4084=1.3d1*pd2
  u840=u2858*(u971+u3818*(8.0d1+u4084))*u4147
  u565=3.38070074315445273d0*u840
  u971=u2858*(u49+u3788*(7.5d1*u4031+u3786))*u4147
  u3818=-7.04312654823844319d-2*u971
  u986=u2858*(u562+u3788*(9.0d0*u526+u645))*u4147
  u794=-3.5215632741192216d-2*u986
  u438=u2858*(2.3349651778357119d7*u534+u3788*(u489+u4009))*u4147
  u489=-1.05646898223576648d-1*u438
  u3929=3.5215632741192216d-2*u842
  u3928=u648*u2855
  value(2)=value(2)+(c(43)*(y*(u2856*(u537+u2855*(u2855*(u3929+u3898)+u&
 &3818)+u2856*((u416+u3904)*u2856+u2855*(u794+u4088)+u737))+u2855*(u565&
 &+u2855*(u3928+u489))+u616)))
  u3929=3.38823854681143158d-1*u848
  u3818=-1.69411927340571579d-1*u952
  u616=7.52941899291429239d-2*u498
  value(3)=value(3)+(c(43)*(u2983*(u2856*(u3818+u2855*(u3860+u616)+(u38&
 &60+u646)*u2856)+u2855*(u3818+u3919)+u3929)*z))
  u3929=u744*u2855
  u616=(y*(u2857*(u4148+u2856*(u631*u2856+u4197)+u2857*((u689+u3927)*u2&
 &857+u567+u688))+u2856*(u4053+u2855*(u2855*(u3841+u676)+u691)+u2856*((&
 &u4026+u676)*u2856+u2855*(u45+u4089)+u692))+u2855*(u453+u2855*(u3929+u&
 &690))+u4121))
  value(4)=value(4)+(c(43)*u616)
  u794=4.04071024781624764d-1*u2858
  u567=u794*(5.04696859799200284d4*u534-1.84335200494173667d2*u382)*u41&
 &47
  u492=4.4086755105988966d6*u534
  u4116=u2858*(u492+u3788*(1.35d2+u3888))*u4147
  u812=-1.68362926992343652d-2*u4116
  u823=1.27361736972857013d6*u534
  u963=u2858*(u823+u3788*(u4174+u3832))*u4147
  u4174=-3.36725853984687303d-2*u963
  u434=1.68362926992343652d-2*u948
  u811=+u739*(2.57989159509120616d6*u534+u3788*(7.9d1+u4015))*u4147
  u739=1.51526634293109286d-1*u952
  u420=-3.03053268586218573d-1*u948
  u972=u3795*(7.08980335815570705d7*u534+u3788*(1.67d2*u4031+u4035))*u4&
 &147
  u58=-5.05088780977030955d-2*u701
  u933=5.05088780977030955d-2*u700
  u490=-1.57840244055322173d-1*u948
  u579=-2.02035512390812382d-1*u2858
  u505=+u579*(1.40424479226483373d6*u534+u3827*(8.6d1+u4084))*u4147
  u81=3.36725853984687303d-2*u971
  u813=1.68362926992343652d-2*u986
  u59=1.02280478147848768d0*u3808
  u649=u3838*(1.95712535814956943d8*u534+u3788*(4.61d2*u4031+u977))*u41&
 &47
  u425=-1.68362926992343652d-2*u842
  u792=3.36725853984687303d-2*u700
  u3930=u933*u2856
  u44=u2856*(u58+u3930)
  u3931=u933*u2855
  value(5)=value(5)+(c(43)*(z*(u2857*(u812+u2855*(u2855*(u425+u3931)+u8&
 &1)+u2856*(u44+u739)+u2857*((u434+u3920)*u2857+u2856*(u420+u3930)+u285&
 &5*(u813+u792*u2855)+u4174))+u2856*(u811+u2855*(u2855*(u3863+u3896)+u5&
 &9)+u2856*((u490+u3896)*u2856+u2855*(u3863+u3918)+u972))+u2855*(u505+u&
 &2855*(u490*u2855+u649))+u567)))
  value(6)=value(6)+(c(43)*u814)
  u567=3.75553839791757858d6*u534
  u39=(1.15d2+u3790)
  u649=u2858*(u567+u3788*u39)*u4147
  u434=-1.12941284893714386d-1*u649
  u811=1.69411927340571579d-1*u952
  u972=u2858*(5.30673904053570887d7*u534+u3788*(1.25d2*u4031+u600))*u41&
 &47
  u505=2.82353212234285965d-2*u972
  u81=-5.6470642446857193d-2*u701
  u813=-1.78823701081714444d-1*u948
  u59=1.12941284893714386d-1*u649
  u649=3.38823854681143158d-1*u948
  u425=-2.82353212234285965d-2*u972
  u792=5.6470642446857193d-2*u701
  u4174=-5.6470642446857193d-2*u700
  u812=1.41176606117142983d-1*u948
  u814=1.78823701081714444d-1*u948
  u3772=u4174*u2855
  u4088=u2855*(u792+u3772)
  u3932=u670*u2856
  u699=u2856*(u81+u3932)
  value(7)=value(7)+(c(43)*(z*(u2857*(u2855*(u4088+u3818)+u2856*(u699+u&
 &811)+(u2856*(u3871+u3932)+u2855*(u649+u3772))*u2857)+u2856*(u434+u999&
 &*(u4065+u812)+u2856*((u813+u3769)*u2856+u478*u2855+u505))+u2855*(u59+&
 &u2855*(u814*u2855+u425)))))
  u434=8.45175185788613183d-1*u508
  u505=-u565
  u813=1.05646898223576648d-1*u438
  u565=-u537
  u425=7.04312654823844319d-2*u971
  u814=-3.5215632741192216d-2*u842
  u438=3.5215632741192216d-2*u986
  u840=7.04312654823844319d-2*u700
  u986=(u750+u3975)
  u4101=u840*u2855
  value(8)=value(8)+(c(43)*(x*(u2856*(u505+u2855*(u2855*(u438+u4024)+u4&
 &25)+u2856*(u986*u2856+u2855*(u814+u4101)+u813))+u2855*(u565+u2855*(u7&
 &37*u2855+u416))+u434)))
  u434=-u3815
  u505=5.97629105161728553d-1*u376
  u425=-3.73518190726080346d-2*u3808
  u438=-1.56877640104953745d0*u3808
  u814=2.61462733508256242d-1*u948
  value(9)=value(9)+(c(43)*((u2856*(u505+u2855*(u2855*(u814+u763)+u438)&
 &+u2856*(u746+u2855*(u814+u675)+u425))+u2855*(u505+u2855*(u4087+u425))&
 &+u434)*z))
  u438=-5.97629105161728553d-1*u508
  u425=-2.98814552580864276d-1*u3809
  u4089=7.0d0*pd2
  u505=u2858*(1.46955850353296553d6*u534+u3788*(4.5d1+u4089))*u4147
  u537=9.96048508602880921d-2*u505
  u565=u2858*(u815+u3788*(u424+u645))*u4147
  u815=-4.98024254301440461d-2*u565
  u424=-4.98024254301440461d-2*u958
  u59=1.33729823821499863d8*u534
  u746=u2858*(u59+u3788*(2.1d1*u526+u645))*u4147
  u816=4.98024254301440461d-2*u746
  u740=4.98024254301440461d-2*u948
  value(1)=value(1)+(c(44)*(x*(u2856*(u425+u2855*(u2855*(u816+u4030)+u4&
 &80)+u2856*((u872+u4061)*u2856+u815*u2855+u637))+u2855*(u537+u2855*(u7&
 &40*u2855+u424))+u438)))
  u438=-1.26776277868291978d0*u508
  u425=+u438
  u537=2.11293796447153296d-1*u3946
  u815=-3.16940694670729944d-1*u897
  u424=4.22587592894306592d-1*u515
  u816=9.50822084012189831d-1*u897
  u421=-3.16940694670729944d-1*u952
  u566=-1.26776277868291978d0*u515
  u737=+u566
  u3933=u4*u2856
  u865=u2856*(u55+u3933)
  value(2)=value(2)+(c(44)*(z*(u2857*(u537+u2855*(u3993+u421)+u2856*(u8&
 &65+u731)+(u2856*(u389+u3933)+u3861+u817)*u2857)+u2856*(u815+u2856*(u3&
 &971*u2856+u424))+u2855*(u816+u2855*(u3928+u737))+u425)))
  u425=-1.12941284893714386d0*u508
  u537=1.01647156404342947d0*u897
  u815=-1.35529541872457263d0*u515
  u424=-2.25882569787428772d-1*u818
  u816=5.6470642446857193d-2*u955
  u737=1.8823547482285731d-1*u2858*(2.93911700706593106d5*u534+u3788*u8&
 &3)*u4147
  u817=3.7647094964571462d-2*u82
  u82=-1.8823547482285731d-2*u915
  u818=u4042*(3.82085210918571038d6*u534+u3788*(9.0d0*u4031+u4006))*u41&
 &47
  u83=-1.8823547482285731d-2*u517
  u3861=3.7647094964571462d-2*u700
  u491=-1.8823547482285731d-2*u948
  u4102=u3861*u2855
  value(3)=value(3)+(c(44)*(x*(u2857*(u537+u2855*(u687*u2855+u4192)+u28&
 &57*((u649+u3854)*u2857+u4119+u815))+u2856*(u424+u2855*(u2855*(u83+u39&
 &73)+u817)+u2856*((u4200+u3973)*u2856+u2855*(u82+u4102)+u816))+u2855*(&
 &u737+u2855*(u491*u2855+u818))+u425)))
  u4119=-1.59723096848534447d-1*u508
  u3920=2.2d1*pd2
  u4083=u4027*(u492+u3788*(1.35d2+u3920))*u4147
  u492=+u4023*(1.65570258064714117d7*u534+u3788*(u819+u3891))*u4147
  u819=5.32410322828448158d-2*u948
  u710=-3.59376967909202506d-1*u897
  u838=4.79169290545603342d-1*u515
  u4115=3.99307742121336118d-2*u524
  u490=-1.3310258070711204d-2*u2858
  u822=8.8d1*u4011
  u493=+u490*(1.46465997518785565d8*u534+u3788*(3.45d2*u4031+u822))*u41&
 &47
  u474=u4027*(u820+u3788*(u3803+u727))*u4147
  u820=3.19446193697068895d-1*u84
  u84=+u490*(u642+u3788*(u85+u603))*u4147
  u85=-1.3310258070711204d-2*u700
  u3803=u947*u2856
  u4056=u2856*(u3862+u3803)
  u568=(z*(u2857*(u4083+u2855*(u2855*(u84+u676)+u493)+u2856*(u4056+u606&
 &)+u2857*((u819+u3851)*u2857+u2856*(u919+u3803)+u2855*(u474+u85*u2855)&
 &+u492))+u2856*(u710+u2856*(u744*u2856+u838))+u2855*(u4115+u2855*(u392&
 &9+u820))+u4119))
  value(4)=value(4)+(c(44)*u568)
  u87=-1.47468160395338933d3*u382
  u518=u3795*(3.69121428833532757d5*u534+u87)*u4147
  u3851=-5.68224878599159824d-2*u2858
  u422=+u3851*(6.49871427117911424d6*u534+u3788*(1.99d2+u821))*u4147
  u821=u494*(u423+u3788*(u3942+u3833))*u4147
  u423=-6.31360976221288694d-3*u4116
  u494=-1.26272195244257739d-2*u963
  u3942=3.59225411974724908d5*u534
  u480=5.1d1*pd2
  u478=u794*(u3942-3.54490770181103205d0*u382*(3.52d2+u480))*u4147
  u794=-7.76574000752185093d-1*u3808
  u495=+u608*(u3943+u3788*(3.75d2*u4031+u822))*u4147
  u3943=u3838*(2.48355387097071175d8*u534+u3788*(3.9d1*u526+u645))*u414&
 &7
  u822=+u990*(u823+u3827*u824)*u4147
  u823=2.08349122153025269d-1*u948
  u824=3.15680488110644347d-2*u948
  value(5)=value(5)+(c(44)*(y*(u2857*(u422+u2855*(u2855*(u823+u3896)+u7&
 &94)+u2856*(u420*u2856+u739)+u2857*((u548+u3930)*u2857+u44+u2855*(u686&
 &+u3768)+u821))+u2856*(u423+u2855*(u724+u495)+u2856*((u694+u3829)*u285&
 &6+u2855*(u3943+u3923)+u494))+u2855*(u478+u2855*(u824*u2855+u822))+u51&
 &8)))
  value(6)=value(6)+(c(44)*u90)
  u990=3.67059175904571754d-1*u376
  u90=1.57079475599856982d7*u534
  u724=u2858*(u90+u3788*(4.81d2+u3912))*u4147
  u44=-2.82353212234285965d-2*u724
  u548=1.69411927340571579d-1*u958
  u507=-9.41177374114286549d-3*u4116
  u4120=-1.8823547482285731d-2*u963
  u3927=5.0d0*pd2
  u4117=u2858*(u3942+u3827*(2.2d1+u3927))*u4147
  u3942=1.12941284893714386d-1*u4117
  u544=1.60941330973543d0*u3808
  u496=-4.70588687057143275d-2*u962
  u826=9.41177374114286549d-3*u43
  u930=-1.12941284893714386d-1*u3808
  u825=9.4117737411428655d-2*u948
  u4103=u825*u2855
  value(7)=value(7)+(c(44)*(y*(u2857*(u44+u2855*(u2855*(u415+u4065)+u54&
 &4)+u2856*(u3871*u2856+u811)+u2857*((u674+u3932)*u2857+u699+u3881+u548&
 &))+u2856*(u507+u2855*(u4103+u496)+u2856*((u736+u4065)*u2856+u2855*(u8&
 &26+u4065)+u4120))+u2855*(u3942+u2855*(u4092+u930))+u990)))
  u4120=-3.16940694670729944d-1*u848
  u507=-4.22587592894306592d-1*u857
  value(8)=value(8)+(c(44)*(u528*((u852+u2855*(u3904+u99))*u2857+u2856*&
 &(u648+u2855*(u3975+u507)+(u3975+u750)*u2856)+u2855*(u813+u3922)+u4120&
 &)))
  u630=3.73518190726080346d-2*u2858
  u507=u630*(4.00788682781717872d5*u534+u87)*u4147
  u4120=-1.86759095363040173d-1*u376
  u3942=-1.24506063575360115d-2*u4116
  u496=-2.4901212715072023d-2*u963
  u826=1.24506063575360115d-2*u948
  u990=-5.97629105161728553d-1*u2858
  u44=+u990*(u4118+u827*(2.0d2+3.1d1*pd2))*u4147
  u4118=5.60277286089120519d-1*u3808
  u548=1.24506063575360115d-2*u2858
  u827=u548*(1.9741069230792837d8*u534+u3788*(4.65d2*u4031+u86))*u4147
  u86=+u497*(-u3788*(u645-9.0d0*u3894)+u562)*u4147
  u811=8.9644365774259283d-1*u515
  u497=-1.86759095363040173d-1*u948
  u87=-4.98024254301440461d-2*u498
  u498=-2.36561520793184219d-1*u948
  value(9)=value(9)+(c(44)*(y*((u4120+u2855*(u2855*(u497+u763)+u4118))*&
 &u2857+u2856*(u3942+u2855*(u2855*(u87+u756)+u827)+u2856*((u826+u3771)*&
 &u2856+u2855*(u86+u3856)+u496))+u2855*(u44+u2855*(u498*u2855+u811))+u5&
 &07)))
  u562=-4.48221828871296415d-1*u3946
  u3881=-4.98024254301440461d-1*u948
  u699=-1.99209701720576184d-1*u857
  value(1)=value(1)+(c(45)*(u528*((u4008+u2855*(u4061+u3881))*u2857+u28&
 &56*(u637+u2855*(u4061+u699)+(u4061+u872)*u2856)+u2855*(u803+u3857)+u5&
 &62)))
  u562=u2858*(3.85944657493506099d4*u534+u3788)*u4147
  u699=3.16940694670729944d-1*u562
  u3830=-3.16940694670729944d-1*u520
  u828=1.05646898223576648d-1*u970
  u384=-1.90164416802437966d0*u376
  u829=+u384
  u705=3.80328833604875932d0*u3808
  value(2)=value(2)+(c(45)*(y*(u2857*(u3830+u2855*(u4093+u705)+u2856*(u&
 &389*u2856+u731)+u2857*((u3971+u3933)*u2857+u865+u50+u828))+u2856*(u49&
 &9+u765*u2856)+u2855*(u829+u4094)+u699)))
  u699=+u663*(u88+u3788*u4171)*u4147
  u3830=-u3964
  u829=5.6470642446857193d-2*u969
  u499=-9.4117737411428655d-1*u948
  u88=-7.52941899291429239d-2*u857
  u81=u3973+u88
  value(3)=value(3)+(c(45)*(u528*(u2857*(u3830+u2855*(u4025+u499)+(u382&
 &8+u3871)*u2857)+u2856*(u624+u2855*(u81)+(u3973+u4200)*u2856)+u829*u28&
 &55+u699)))
  u3964=u639*(3.33165900913197573d4*u534+u3788)*u4147
  u4171=+u697*(2.90646015143186516d6*u534+u3788*(8.9d1+u4015))*u4147
  u697=u3876*(u741+u3788*(u831+u4009))*u4147
  u4087=-1.73033354919245651d-1*u948
  u75=-7.98615484242672237d-2*u376
  u639=9.58338581091206684d-1*u3808
  u50=(y*(u2857*(u4171+u2855*(u4095+u639)+u2856*(u919*u2856+u606)+u2857&
 &*((u4087+u3767+u3803)*u2857+u4056+u2855*(u3761+u430)+u697))+u2856*(u9&
 &97+u516*u2856)+u2855*(u75+u4096)+u3964))
  value(4)=value(4)+(c(45)*u50)
  u865=u3795*(8.85693508863302459d5*u534-3.5165484401965438d3*u382)*u41&
 &47
  u382=+u4044*(u830+u500*(5.0d2+7.1d1*pd2))*u4147
  u830=u3838*(9.12759114972141925d7*u534+u3788*(2.15d2*u4031+u4035))*u4&
 &147
  u500=-6.73451707969374607d-2*u948
  u4056=5.68224878599159824d-2*u897
  u831=-1.89408292866386608d-2*u952
  u741=3.78816585732773216d-2*u948
  u501=-7.57633171465546432d-2*u515
  u966=6.31360976221288694d-3*u701
  u855=-4.41952683354902085d-2*u2858*(u567+u3788*(1.15d2+u4051))*u4147
  u567=u3838*(4.05434862696928157d8*u534+u3788*(9.55d2*u4031+2.72d2*u40&
 &11))*u4147
  u4175=+u4013*(3.37508602978071084d8*u534+u3788*(5.3d1*u526+u727))*u41&
 &47
  u868=+u4045*(u96+u3788*(u3859+u3891))*u4147
  u96=6.31360976221288694d-3*u3858
  u3858=u905*u2856
  u3859=u636*u2855
  u4104=u636*u2856
  value(5)=value(5)+(c(45)*(z*(u2857*(u382+u2855*(u2855*(u96+u3829)+u56&
 &7)+u2856*(u2856*(u966+u3858)+u831)+u2857*((u500+u3917)*u2857+u2856*(u&
 &741+u3858)+u2855*(u4175+u3770)+u830))+u2856*(u4056+u2856*(u4104+u501)&
 &)+u2855*(u855+u2855*(u3859+u868))+u865)))
  u4056=-3.59376967909202506d-1*u524
  u831=1.19792322636400836d-1*u969
  u501=-2.79515419484935283d-1*u948
  u966=-3.99307742121336118d-2*u848
  value(6)=value(6)+(c(45)*(x*(u2857*(u4056+u606*u2855+u2857*((u501+u39&
 &24)*u2857+u3841*u2855+u831))+u966*u2855+u4146)))
  u855=u2858*(3.79347312920967534d4*u534+u3788)*u4147
  u567=7.62353673032572105d-1*u855
  u4175=u2858*(5.06181262328021461d6*u534+u3788*(1.55d2+u4213))*u4147
  u868=-5.6470642446857193d-2*u4175
  u96=u2858*(2.25005735318714056d7*u534+u3788*(5.3d1*u4031+u3891))*u414&
 &7
  u3930=2.82353212234285965d-2*u96
  u834=-4.23529818351428947d-1*u524
  u3917=1.41176606117142983d-1*u2858
  u524=u3917*(3.69349037221285337d7*u534+u3788*(8.7d1*u4031+u3786))*u41&
 &47
  u586=-5.6470642446857193d-2*u517
  u517=1.27058945505428684d0*u3808
  u379=-9.88236242820000877d-1*u948
  u4105=u386*u2856
  u4106=u695*u2855
  value(7)=value(7)+(c(45)*(z*(u2857*(u868+u2855*(u2855*(u379+u4065)+u5&
 &24)+u2856*(u2856*(u832+u4105)+u487)+u2857*((u4200+u4038)*u2857+u2856*&
 &(u624+u4105)+u2855*(u586+u3828)+u3930))+u2856*(u502+u2856*(u695*u2856&
 &+u807))+u2855*(u834+u2855*(u4106+u517))+u567)))
  u4105=(8.5d1+u4015)
  u868=u2858*(2.77583272889560156d6*u534+u3788*u4105)*u4147
  u3930=-3.16940694670729944d-1*u868
  u834=1.79612705987362454d6*u534
  u524=u2858*(u834+u3788*u833)*u4147
  u586=-1.05646898223576648d-1*u524
  u517=1.05646898223576648d-1*u972
  u379=u2858*(5.41287382134642304d8*u534+u3788*(8.5d1*u526+u603))*u4147
  u502=-3.5215632741192216d-2*u379
  u832=4.22587592894306592d-1*u3808
  u833=-8.45175185788613183d-1*u948
  value(8)=value(8)+(c(45)*(x*(u2857*(u3930+u2855*(u833*u2855+u517)+u28&
 &57*(u986*u2857+u2855*(u502+u3925)+u801))+u2855*(u586+u832*u2855)+u835&
 &)))
  u3930=5.60277286089120519d-1*u376
  u586=-2.98814552580864277d0*u964
  u517=1.12055457217824104d-1*u897
  u502=-3.73518190726080346d-2*u952
  u832=7.47036381452160691d-2*u948
  u833=-1.49407276290432138d-1*u515
  u964=1.24506063575360115d-2*u701
  u3906=-3.36166371653472311d-1*u2858
  u487=(5.5d1+u3994)
  u567=+u3906*(u834+u3788*u487)*u4147
  u834=1.12055457217824104d-1*u971
  u971=1.12055457217824104d-1*u2858
  u835=u971*(u4199+u3788*(u89+u3891))*u4147
  u4199=3.50244776675356785d8*u534
  u801=-3.73518190726080346d-2*u2858
  u89=+u801*(u4199+u3788*(5.5d1*u526+u727))*u4147
  u986=8.71542445027520806d-2*u700
  u3934=u986*u2855
  u4107=u530*u2856
  u4108=u4218*u2855
  value(9)=value(9)+(c(45)*(z*(u2857*(u586+u2855*(u2855*(u89+u3934)+u83&
 &4)+u2856*(u2856*(u964+u4107)+u502)+(u2856*(u832+u4107)+u2855*(u417+u7&
 &56)+u576)*u2857)+u2856*(u517+u2856*(u696*u2856+u833))+u2855*(u567+u28&
 &55*(u4108+u835))+u3930)))
  u4123=-u3935
  u4095=-8.9644365774259283d-1*u3808
  u803=(u740+u4030)*u2856
  u3935=u3881*u2855
  u3936=u4008*u2855
  value(1)=value(1)+(c(46)*(y*(u2856*(u445+u2855*(u3935+u408)+u2856*(u8&
 &03+u2855*(u728+u4061)+u4095))+u2855*(u521+u3936)+u4123)))
  u4095=-u698
  u698=4.22587592894306592d-1*u948
  u3773=u3997*u2856
  u4094=u3773+u3975
  u3937=u760*u2855
  u4109=u771*u2855
  value(2)=value(2)+(c(46)*(u2983*(u2856*(u863+u4109+u2856*(u4094+u698)&
 &)+u3937+u4095)*z))
  value(3)=value(3)+(c(46)*(y*(u2856*(u3774+u2855*(u778*u2855+u464)+u28&
 &56*((u780+u3860)*u2856+u3775+u779))+u2855*(u3974+u552*u2855)+u4111)))
  u4095=u747*u2856
  u3974=u4095+u430
  u3774=u860*u2856
  u4110=u622*u2855
  u4111=u704*u2855
  u464=(u528*((u69+u2856*(u3774+u536))*u2857+u2856*(u4131+u4110+u2856*(&
 &u3974+u680))+u4111+u3969))
  value(4)=value(4)+(c(46)*u464)
  u4003=u909*u2856
  u814=u3800+u4003
  u3775=u3869*u2856
  u3860=u4198*u2856
  value(5)=value(5)+(c(46)*(x*(u2857*(u4113+u2856*(u2856*(u783+u3775)+u&
 &876)+(u2856*(u465+u3860)+u4112)*u2857)+u2856*(u712+u2855*(u466*u2855+&
 &u781)+u2856*(u2856*(u467+u814)+u2855*(u71+u444)+u782))+u2855*(u3804+u&
 &3776*u2855)+u4007)))
  value(6)=value(6)+(c(46)*u472)
  u3938=-5.08235782021714737d-1*u4114
  u3939=8.47059636702857895d-1*u376
  u3804=5.6470642446857193d-2*u3799
  u712=-u943
  u3922=-2.82353212234285965d-2*u981
  u503=8.47059636702857895d-1*u948
  u3834=1.12941284893714386d-1*u703
  u4112=-5.6470642446857193d-2*u980
  u3776=3.01176759716571696d-1*u849
  u467=-1.97647248564000175d-1*u3808
  u779=u708*u2856
  u726=u3772+u779
  u4007=u4174*u2856
  value(7)=value(7)+(c(46)*(x*((u3939+u2856*(u2856*(u503+u4007)+u712))*&
 &u2857+u2856*(u3804+u2855*(u3919+u4112)+u2856*(u2856*(u419+u726)+u2855&
 &*(u3776+u3802)+u3922))+u2855*(u3834+u467*u2855)+u3938)))
  u3776=-u4122
  u467=-u973
  u3938=-u535
  u3939=9.50822084012189831d-1*u376
  u3804=-u3806
  u712=-u651
  u3834=(u614+u3898)
  u4112=u765*u2855
  value(8)=value(8)+(c(46)*((u2856*(u467+u2855*(u3916+u3804)+u2856*(u38&
 &34*u2856+u2855*(u712+u3904)+u3938))+u2855*(u3939+u4112)+u3776)*z))
  u3776=u4215*u2856
  u467=u675+u3776
  u3938=u882*u2855
  u3939=u638*u2855
  value(9)=value(9)+(c(46)*(x*(u2856*(u445+u2855*(u3938+u836)+u2856*(u2&
 &856*(u398+u467)+u2855*(u784+u763)+u585))+u2855*(u521+u3939)+u3967)))
  u3804=-u531
  u712=-5.97629105161728553d-1*u3808
  u3922=-1.49407276290432138d-1*u376
  u503=-4.48221828871296415d-1*u3808
  u836=4.48221828871296415d-1*u948
  value(1)=value(1)+(c(47)*((u2856*(u3804+u2855*(u3857+u503)+u2856*(u80&
 &3+u2855*(u836+u4061)+u712))+u2855*(u3922+u3926)+u434)*z))
  u3804=-9.50822084012189831d-1*u4114
  u712=5.28234491117883239d-1*u376
  u3922=-u3966
  u503=-1.05646898223576648d-1*u851
  u836=5.28234491117883239d-1*u948
  u851=3.16940694670729944d-1*u3946
  u803=u2855*(u647*u2855+u421)
  u4122=u2855*(u851+u951*u2855)
  value(2)=value(2)+(c(47)*(x*((u712+u2856*(u2856*(u836+u3773)+u615))*u&
 &2857+u2856*(u3922+u803+u2856*((u664+u3898)*u2856+u3993+u503))+u4122+u&
 &3804)))
  u3804=u641*u2856
  u3922=u2856*(u758+u3804)
  u4113=u633*u2855
  u4114=u931*u2855
  value(3)=value(3)+(c(47)*(z*(u2857*(u3983+u2856*(u3922+u4192)+(u2856*&
 &(u687+u3804)+u4136)*u2857)+u2856*(u709+u2855*(u4113+u4078)+u2856*((u7&
 &97+u4025)*u2856+u3952+u477))+u2855*(u581+u4114)+u987)))
  u790=u2855*(u919*u2855+u606)
  u4097=u2855*(u997+u516*u2855)
  u683=(x*(u2857*(u4203+u2856*(u2856*(u634+u4095)+u516)+(u2856*(u3843+u&
 &3774)+u4067)*u2857)+u2856*(u519+u790+u2856*((u4193+u676)*u2856+u875+u&
 &685))+u4097+u3949))
  value(4)=value(4)+(c(47)*u683)
  value(5)=value(5)+(c(47)*(u528*(u2857*(u506+u2856*(u3775+u800)+(u3860&
 &+u418)*u2857)+u2856*(u799+u2855*(u444+u76)+u2856*(u4003+u3800+u3863))&
 &+u2855*(u798+u4098)+u3940)))
  value(6)=value(6)+(c(47)*u616)
  u3932=5.6470642446857193d-2*u959
  u4121=-6.77647709362286316d-1*u515
  u506=5.6470642446857193d-1*u948
  u4026=-1.12941284893714386d-1*u982
  u3940=1.8823547482285731d-2*u842
  value(7)=value(7)+(c(47)*(u528*((u4127+u2856*(u4007+u506))*u2857+u285&
 &6*(u4121+u2855*(u3802+u3940)+u2856*(u779+u3772+u812))+u2855*(u4026+u4&
 &099)+u3932)))
  u3940=-u4169
  u4121=-u4172
  u506=2.11293796447153296d-1*u837
  u4026=-1.05646898223576648d-1*u3756
  u4099=-3.5215632741192216d-2*u2858
  u3932=+u4099*(u509+u49)*u4147
  value(8)=value(8)+(c(47)*(y*(u2856*(u4121+u2855*(u2855*(u3932+u4024)+&
 &u506)+u2856*(u450*u2856+u2855*(u46+u4101)+u828))+u2855*(u4121+u2855*(&
 &u3853+u4026))+u3940)))
  value(9)=value(9)+(c(47)*(u2983*(u2856*(u595+u2855*(u763+u802)+u2856*&
 &(u3776+u675+u4219))+u2855*(u595+u4100)+u644)*z))
  u840=5.97629105161728553d-1*u508
  u506=-9.96048508602880921d-2*u505
  u4026=4.98024254301440461d-2*u958
  u3932=2.98814552580864276d-1*u3809
  u837=1.99209701720576184d-1*u948
  u4121=-4.98024254301440461d-2*u746
  u505=4.98024254301440461d-2*u565
  u3940=u637*u2855
  value(1)=value(1)+(c(48)*(y*(u2856*(u506+u2855*(u2855*(u505+u4030)+u8&
 &37)+u2856*((u410+u4061)*u2856+u4121*u2855+u4026))+u2855*(u3932+u2855*&
 &(u3940+u872))+u840)))
  u4026=3.16940694670729944d-1*u848
  u3932=4.22587592894306592d-1*u857
  value(2)=value(2)+(c(48)*(u528*((u863+u2856*(u3773+u664))*u2857+u2856&
 &*(u489+u2855*(u3898+u3932)+(u3898+u698)*u2856)+u2855*(u750+u3928)+u40&
 &26)))
  value(3)=value(3)+(c(48)*(y*(u2857*(u537+u2856*(u687*u2856+u4192)+u28&
 &57*((u649+u3804)*u2857+u3922+u815))+u2856*(u737+u2855*(u2855*(u82+u39&
 &73)+u817)+u2856*((u491+u3973)*u2856+u2855*(u83+u4102)+u818))+u2855*(u&
 &424+u2855*(u3852+u816))+u425)))
  u4026=(u528*(u2857*(u4141+u2856*(u4095+u3914)+(u3774+u3989)*u2857)+u2&
 &856*(u658+u2855*(u3803+u978))+u2855*(u693+u3929)+u4138))
  value(4)=value(4)+(c(48)*u4026)
  u3932=u2855*(u428+u3923)
  u3923=u3931+u3860
  value(5)=value(5)+(c(48)*(x*(u2857*(u422+u2855*(u420*u2855+u739)+u285&
 &6*(u2856*(u823+u4003)+u794)+u2857*((u3976+u3923)*u2857+u2856*(u686+u3&
 &775)+u2855*(u58+u3931)+u821))+u2856*(u478+u2855*(u2855*(u3943+u3829)+&
 &u495)+u2856*((u824+u3770)*u2856+u3932+u822))+u2855*(u423+u2855*(u4091&
 &+u494))+u518)))
  value(6)=value(6)+(c(48)*u568)
  u4119=-3.67059175904571754d-1*u376
  u425=2.82353212234285965d-2*u724
  u506=-1.69411927340571579d-1*u958
  u838=1.69411927340571579d-1*u948
  u424=-1.12941284893714386d-1*u4117
  u4121=-u544
  u422=-9.41177374114286549d-3*u948
  u423=9.41177374114286549d-3*u4116
  u478=4.70588687057143275d-2*u962
  u504=-9.4117737411428655d-2*u948
  u820=1.8823547482285731d-2*u963
  u4116=-9.41177374114286549d-3*u43
  u3922=u2856*(u649+u4007)
  u4115=u422*u2856
  value(7)=value(7)+(c(48)*(x*(u2857*(u425+u2855*(u649*u2855+u3818)+u28&
 &56*(u2856*(u738+u779)+u4121)+u2857*((u838+u3772)*u2857+u3922+u4088+u5&
 &06))+u2856*(u424+u2855*(u2855*(u4116+u3769)+u478)+u2856*(u4115+u2855*&
 &(u504+u3769)+u3941))+u2855*(u423+u2855*(u422*u2855+u820))+u4119)))
  u4119=-u438
  u445=-2.11293796447153296d-1*u3946
  u506=-9.50822084012189831d-1*u897
  u424=-u566
  u4116=3.16940694670729944d-1*u897
  u423=-1.05646898223576648d-1*u952
  u478=-4.22587592894306592d-1*u515
  u820=3.5215632741192216d-2*u701
  u963=u2855*(u820+u4024)
  u3941=u668*u2856
  u4121=u2856*(u46+u3941)
  value(8)=value(8)+(c(48)*(z*(u2857*(u445+u2855*(u963+u423)+u2856*(u41&
 &21+u734)+(u2856*(u613+u3941)+u3849+u3784)*u2857)+u2856*(u506+u2856*(u&
 &750*u2856+u424))+u2855*(u4116+u2855*(u3853+u478))+u4119)))
  value(9)=value(9)+(c(48)*(x*((u4120+u2856*(u2856*(u497+u3776)+u4118))&
 &*u2857+u2856*(u44+u2855*(u2855*(u86+u3771)+u827)+u2856*((u498+u756)*u&
 &2856+u2855*(u87+u3856)+u811))+u2855*(u3942+u2855*(u826*u2855+u496))+u&
 &507)))
  u4119=-4.48221828871296415d-1*u897
  u445=5.97629105161728553d-1*u515
  u506=4.48221828871296415d-1*u897
  u424=-1.49407276290432138d-1*u952
  u478=-5.97629105161728553d-1*u515
  u87=4.98024254301440461d-2*u701
  u811=u2855*(u87+u4030)
  u3942=u3782*u2856
  u4120=u2856*(u57+u3942)
  value(1)=value(1)+(c(49)*(z*(u2857*(u2855*(u811+u424)+u2856*(u4120+u7&
 &35)+(u2856*(u677+u3942)+u588)*u2857)+u2856*(u4119+u2856*(u872*u2856+u&
 &445))+u2855*(u506+u2855*(u3940+u478)))))
  u4119=u2858*(3.0677652262304331d4*u534+u3788)*u4147
  u445=-9.50822084012189831d-1*u4119
  u506=u2858*(2.51457788382307436d6*u534-1.13437046457953026d2*u3883)*u&
 &4147
  u478=3.16940694670729944d-1*u506
  u3883=u2858*(2.58968865178142593d7*u534+u3788*(6.1d1*u4031+u4009))*u4&
 &147
  u4118=-1.05646898223576648d-1*u3883
  u4116=6.33881389341459887d-1*u376
  u425=-1.26776277868291978d0*u3808
  u505=+u425
  u737=u2856*(u623+u3773)
  value(2)=value(2)+(c(49)*(x*(u2857*(u478+u803+u2856*(u623*u2856+u505)&
 &+u2857*((u648+u3898)*u2857+u737+u3993+u4118))+u2856*(u4116+u912*u2856&
 &)+u4122+u445)))
  u445=-2.71059083744914526d0*u508
  u478=2.71059083744914526d0*u3944
  u4118=-1.12941284893714386d-1*u970
  u4116=-1.69411927340571579d-1*u897
  u505=2.25882569787428772d-1*u515
  u839=1.69411927340571579d-1*u659
  u851=-5.6470642446857193d-2*u983
  u3944=2.25882569787428772d-1*u869
  u724=u2858*(u957+u610)*u4147
  u4117=-1.12941284893714386d-1*u724
  u44=u2858*(u642+u885)*u4147
  u84=1.8823547482285731d-2*u44
  u840=-9.4117737411428655d-2*u700
  u3861=u523*u2856
  u803=u2856*(u41+u3861)
  u3943=u4200*u2856
  value(3)=value(3)+(c(49)*(z*(u2857*(u478+u2855*(u2855*(u84+u3973)+u85&
 &1)+u2856*(u803+u666)+u2857*((u633+u3854)*u2857+u2856*(u874+u3861)+u28&
 &55*(u3944+u840*u2855)+u4118))+u2856*(u4116+u2856*(u3943+u505))+u2855*&
 &(u839+u2855*(u3852+u4117))+u445)))
  u445=u676+u3774
  u478=(x*(u2857*(u4171+u790+u2856*(u632*u2856+u639)+u2857*((u4087+u445&
 &)*u2857+u2856*(u3761+u4095)+u875+u697))+u2856*(u75+u3977*u2856)+u4097&
 &+u3964))
  value(4)=value(4)+(c(49)*u478)
  u851=u4032*(8.65406674302746369d6*u534+u3788*(2.65d2+u3790))*u4147
  u84=-1.07331365957619078d0*u3808
  u505=-7.57633171465546432d-2*u4202
  u839=4.41952683354902085d-1*u948
  u4117=5.05088780977030955d-2*u857
  u3944=u4198*u2857
  value(5)=value(5)+(c(49)*(u528*(u2857*(u84+u2855*(u3896+u839)+u2856*(&
 &u4003+u839)+u2857*(u3944+u3775+u3768+u418))+u2856*(u505+u2855*(u3770+&
 &u4117)+(u3770+u636)*u2856)+u2855*(u505+u3859)+u851)))
  value(6)=value(6)+(c(49)*u50)
  u851=-5.6470642446857193d-1*u3808
  u84=-2.82353212234285965d-2*u948
  u4117=5.6470642446857193d-1*u3808
  u840=-4.70588687057143275d-1*u948
  u4116=u84*u2856
  value(7)=value(7)+(c(49)*(u528*(u2857*(u2855*(u4065+u840)+u2856*(u779&
 &+u729)+(u4007+u4038)*u2857)+u2856*(u851+u4116)+u2855*(u4117+u4106))))
  u851=9.50822084012189831d-1*u4119
  u4117=-3.16940694670729944d-1*u506
  u840=1.05646898223576648d-1*u3883
  u4118=-6.33881389341459887d-1*u376
  u3883=-u425
  value(8)=value(8)+(c(49)*(y*(u2857*(u4117+u2855*(u3916+u3883)+u2856*(&
 &u613*u2856+u734)+u2857*((u750+u3941)*u2857+u4121+u3798+u840))+u2856*(&
 &u385+u510*u2856)+u2855*(u4118+u4112)+u851)))
  u3883=-2.24110914435648207d-1*u848
  u851=3.73518190726080346d-1*u3808
  u840=1.49407276290432138d-1*u982
  u4117=-2.98814552580864276d-1*u857
  value(9)=value(9)+(c(49)*(u528*((u851+u2855*(u763+u398)+u2856*(u3776+&
 &u398))*u2857+u2856*(u840+u2855*(u756+u4117)+(u756+u4218)*u2856)+u2855&
 &*(u840+u4108)+u3883)))
  u3883=4.48221828871296415d-1*u376
  u851=u2858*(8.81735102119779319d5*u534+u3788*u381)*u4147
  u4117=-4.48221828871296415d-1*u851
  u4118=1.49407276290432138d-1*u958
  u781=-8.9644365774259283d-1*u376
  u4108=1.79288731548518566d0*u3808
  value(1)=value(1)+(c(50)*(y*(u2857*(u4117+u2855*(u3857+u4108)+u2856*(&
 &u677*u2856+u735)+u2857*((u872+u3942)*u2857+u4120+u3763+u4118))+u2856*&
 &(u985+u841*u2856)+u2855*(u781+u3926)+u3883)))
  u3883=2.11293796447153296d0*u3808
  u4117=u664*u2856
  u4118=u863*u2856
  value(2)=value(2)+(c(50)*(u528*(u2857*(u3883+u4109+u4117+(u4094+u389)&
 &*u2857)+u4118+u3937+u862)))
  u781=u594*(5.64072960952047376d4*u534+u3788)*u4147
  u985=+u663*(1.73081334860549274d6*u534+u3788*(5.3d1+u4019))*u4147
  u510=u2858*(u953+u3788*(u844+u4006))*u4147
  u841=5.6470642446857193d-2*u510
  u819=-1.35529541872457263d0*u376
  u848=4.40471011085486105d0*u3808
  u506=-1.5811779885120014d0*u948
  value(3)=value(3)+(c(50)*(y*(u2857*(u985+u2855*(u4113+u848)+u2856*(u8&
 &74*u2856+u666)+u2857*((u674+u3828+u3861)*u2857+u803+u2855*(u506+u4025&
 &)+u841))+u2856*(u3978+u625*u2856)+u2855*(u819+u4114)+u781)))
  u862=-1.99653871060668059d0*u376
  u803=3.19446193697068895d0*u3808
  u425=-8.78477032666939461d-1*u948
  u4119=u622*u2856
  u4120=u860*u2857
  u4121=u704*u2856
  u496=(u528*(u2857*(u803+u4110+u4119+u2857*(u4120+u3974+u425))+u4121+u&
 &4111+u862))
  value(4)=value(4)+(c(50)*u496)
  u738=+u4052*(5.37683582661893113d4*u534+u3788)*u4147
  u518=u4032*(1.08094192148758131d7*u534+u3788*(3.31d2+u3901))*u4147
  u507=+u3875*(4.6699303556714238d6*u534+u3788*(u511+u4009))*u4147
  u508=-2.31499024614472521d-1*u948
  u984=6.43988195745714467d-1*u376
  u970=-2.04560956295697537d0*u3808
  u842=7.19751512892269111d-1*u948
  u758=u975*(u873+u3788*u3945)*u4147
  u509=+u4029*(u522+u3788*(u3779+u3891))*u4147
  u3847=2.02035512390812382d-1*u2858
  u642=u3847*(u49+u3827*(u667+u4074))*u4147
  u49=-6.31360976221288694d-2*u700
  u4097=-1.89408292866386608d-2*u3808
  u497=u49*u2855+u3775+u3944
  value(5)=value(5)+(c(50)*(x*(u2857*(u518+u2855*(u741*u2855+u509)+u285&
 &6*(u3835*u2856+u970)+u2857*(u2857*(u508+u497)+u2856*(u842+u4003)+u285&
 &5*(u642+u3829)+u507))+u2856*(u984+u988*u2856)+u2855*(u758+u4097*u2855&
 &)+u738)))
  u522=5.98961613182004178d-1*u2858*(4.84904826081584586d4*u534+u3788)*&
 &u4147
  u3942=1.60018592606922914d6*u534
  u4113=-1.99653871060668059d-1*u2858
  u3827=(4.9d1+u4019)
  u667=+u4113*(u3942+u3788*u3827)*u4147
  u843=3.99307742121336118d-2*u510
  u510=-9.31718064949784276d-2*u948
  u4114=+u4113*(u3942-1.13437046457953026d2*u3884)*u4147
  u3942=1.99653871060668059d-1*u968
  u3884=+u3885*(u936+u3788*(u853+u645))*u4147
  value(6)=value(6)+(c(50)*(z*(u2857*(u667+u3942*u2855+u2857*((u510+u39&
 &24)*u2857+u3884*u2855+u843))+u4114*u2855+u522)))
  u4114=4.02360789370902398d3*u702
  u3924=-4.40471011085486105d0*u376
  u3884=+u3924
  u3942=1.5811779885120014d0*u3808
  u853=6.77647709362286316d-1*u376
  u4113=-2.20235505542743053d0*u3808
  u474=+u4113
  u844=-5.6470642446857193d-1*u376
  u511=-u654
  u3779=-7.5294189929142924d-1*u948
  u381=-2.82353212234285965d-2*u3808
  u492=u4038+u4007
  u3945=u381*u2855
  u4122=u390*u2856
  value(7)=value(7)+(c(50)*(x*(u2857*(u3884+u2855*(u3855+u511)+u2856*(u&
 &3943+u474)+u2857*((u874+u492)*u2857+u2856*(u733+u779)+u2855*(u3779+u4&
 &065)+u3942))+u2856*(u853+u4122)+u2855*(u844+u3945)+u4114)))
  u3884=u2858*(9.47048813387911121d5*u534+u3788*u512)*u4147
  u853=-5.28234491117883239d-1*u3884
  u474=1.05646898223576648d-1*u955
  u844=-1.58470347335364972d0*u3946
  u511=5.28234491117883239d-1*u2858*(u90+u3788*(3.7d1*u4031+u3891))*u41&
 &47
  u3779=+u4055*(u564+u3788*(2.9d1*u526+u645))*u4147
  u564=-1.40862530964768864d0*u948
  value(8)=value(8)+(c(50)*(z*(u2857*(u853+u2855*(u564*u2855+u511)+u285&
 &7*(u450*u2857+u2855*(u3779+u3925)+u474))+u2855*(u844+u3883*u2855)+u38&
 &3)))
  u853=6.72332743306944622d-1*u376
  u474=-6.72332743306944622d-1*u851
  u844=2.24110914435648207d-1*u958
  u511=-2.24110914435648207d-1*u948
  u3779=-2.24110914435648207d-1*u376
  u564=-1.49407276290432138d-1*u703
  u845=7.47036381452160691d-2*u980
  u90=-3.98419403441152369d-1*u849
  u849=2.61462733508256242d-1*u3808
  u512=-5.22925467016512484d-1*u948
  u790=u2856*(u882+u3776)
  value(9)=value(9)+(c(50)*(x*(u2857*(u474+u2855*(u512*u2855+u845)+u285&
 &6*(u882*u2856+u3947)+u2857*((u511+u756)*u2857+u790+u2855*(u90+u3934)+&
 &u844))+u2856*(u3779+u638*u2856)+u2855*(u564+u849*u2855)+u853)))
  u3806=-u3807
  u3807=-u4041
  u4038=-u3805
  u3805=u8*u2856
  value(1)=value(1)+(c(51)*(x*(u2856*(u3806+u4126*u2855+u2856*(u2856*(u&
 &4038+u4061+u3805)+u409*u2855+u3807))+u767*u2855+u4123)))
  u3806=-u4124
  u3807=-u4125
  u4038=-u864
  u4126=8.803908185298054d-1*u948
  u767=u3975+u3773
  u4123=u3980*u2855
  u4124=u412*u2855
  u4125=u866*u2855
  value(2)=value(2)+(c(51)*((u2856*(u3807+u4123+u2856*(u2856*(u4126+u76&
 &7)+u4124+u4038))+u4125+u3806)*z))
  u3806=u847*u2856
  u3807=u671*u2856
  value(3)=value(3)+(c(51)*(x*((u764+u2856*(u2856*(u459+u3806)+u4133))*&
 &u2857+u2856*(u569+u4127*u2855+u2856*(u2856*(u4058+u3807)+u772*u2855+u&
 &4128))+u61*u2855+u573)))
  u4038=u430+u4095
  u4126=u601*u2855
  u4127=u657*u2855
  u4128=u4064*u2855
  value(4)=value(4)+(c(51)*(z*((u532+u2856*(u2856*(u462+u3774)+u4224))*&
 &u2857+u2856*(u4173+u4126+u2856*(u2856*(u776+u4038)+u4127+u591))+u4128&
 &+u859)))
  value(5)=value(5)+(c(51)*(y*(u2857*(u4143+u2856*(u2856*(u775+u3775)+u&
 &3777)+(u2856*(u3976+u3860)+u513)*u2857)+u2856*(u925+u2855*(u460*u2855&
 &+u773)+u2856*(u2856*(u461+u814)+u2855*(u70+u444)+u774))+u2855*(u877+u&
 &949*u2855)+u426)))
  value(6)=value(6)+(c(51)*u464)
  u4131=-8.47059636702857895d-1*u7
  u3925=-u871
  u421=2.82353212234285965d-1*u4130
  u3777=-u924
  u4041=u3760*(u858+u380)*u4147
  u4129=-u707
  u4132=-2.25882569787428772d-1*u948
  u4133=-u3950
  u513=-2.82353212234285965d-1*u952
  u464=2.25882569787428772d-1*u95
  u4083=-9.88236242820000877d-1*u3808
  u4130=6.58824161880000585d-1*u948
  value(7)=value(7)+(c(51)*(y*((u3925+u2856*(u2856*(u4129+u4007)+u3777)&
 &)*u2857+u2856*(u421+u2855*(u4130*u2855+u513)+u2856*(u2856*(u4132+u726&
 &)+u2855*(u464+u3802)+u4041))+u2855*(u4133+u4083*u2855)+u4131)))
  u4083=-u4004
  u4130=-u4048
  u4131=-u3780
  u3777=u5*u2856
  u4041=u3777+u3904
  u4129=u852*u2855
  value(8)=value(8)+(c(51)*(u2983*(u2856*(u4130+u3921+u2856*(u4041+u413&
 &1))+u4129+u4083)*z))
  value(9)=value(9)+(c(51)*(y*(u2856*(u867+u2855*(u398*u2855+u42)+u2856&
 &*(u2856*(u463+u467)+u2855*(u777+u763)+u539))+u2855*(u3905+u576*u2855)&
 &+u656)))
  u4083=-u4145
  u4130=-u607
  u4131=8.9644365774259283d-1*u948
  u464=u3805+u4061
  value(1)=value(1)+(c(52)*(u2983*(u2856*(u4130+u3935+u2856*(u464+u4131&
 &))+u3936+u4083)*z))
  u4083=-u4170
  u4130=-u4183
  u4131=-u3839
  u3925=-u946
  u777=-1.05646898223576648d-1*u743
  u4132=7.39528287565036535d-1*u948
  u4133=1.40862530964768864d-1*u948
  u513=-5.28234491117883239d-1*u962
  u421=1.05646898223576648d-1*u43
  u839=u2855*(u561*u2855+u513)
  u505=u2855*(u421+u3898)
  u733=u2855*(u760+u615*u2855)
  value(2)=value(2)+(c(52)*(y*((u4130+u2856*(u2856*(u4132+u3773)+u3925)&
 &)*u2857+u2856*(u4131+u839+u2856*((u4133+u3898)*u2856+u505+u777))+u733&
 &+u4083)))
  u4083=u3807+u4025
  u4130=u684*u2855
  u4131=u913*u2855
  value(3)=value(3)+(c(52)*(u528*((u514+u2856*(u3806+u413))*u2857+u2856&
 &*(u3979+u4130+u2856*(u4083+u633))+u4131+u3810)))
  u3925=u2855*(u414*u2855+u732)
  u777=u2855*(u704+u3778*u2855)
  value(4)=value(4)+(c(52)*(y*(u2857*(u429+u2856*(u2856*(u788+u4095)+u9&
 &2)+(u2856*(u536+u3774)+u69)*u2857)+u2856*(u4134+u3925+u2856*(u2855*(u&
 &4057+u3803)+u787))+u777+u650)))
  u4132=u4137*u2856
  u4133=u3835*u2855
  u4134=u988*u2855
  value(5)=value(5)+(c(52)*(z*(u2857*(u3986+u2856*(u72*u2856+u785)+(u28&
 &56*(u468+u4132)+u3951)*u2857)+u2856*(u998+u2855*(u4133+u91)+u2856*(u2&
 &856*(u470+u3918+u4003)+u2855*(u469+u3896)+u786))+u2855*(u3864+u4134)+&
 &u533)))
  value(6)=value(6)+(c(52)*u683)
  u519=-5.08235782021714737d-1*u706
  u683=6.77647709362286316d-1*u3970
  u772=-1.69411927340571579d-1*u3808
  u92=5.6470642446857193d-2*u388
  u3803=-3.38823854681143158d-1*u955
  u72=-8.47059636702857895d-2*u457
  u676=4.51765139574857544d-1*u4144
  u463=7.52941899291429239d-2*u948
  u514=1.69411927340571579d-1*u376
  value(7)=value(7)+(c(52)*(z*(u2857*(u683+u2856*(u2856*(u676+u3804)+u3&
 &803)+(u3922+u772)*u2857)+u2856*(u92+u2855*(u3855+u4136)+u2856*(u2856*&
 &(u463+u779)+u2855*(u624+u4065)+u72))+u2855*(u514+u3945)+u519)))
  u676=-3.16940694670729944d-1*u3998
  u463=-u575
  u514=-1.05646898223576648d-1*u969
  u92=1.05646898223576648d-1*u3946
  u3803=u2855*(u623*u2855+u423)
  u72=u2855*(u92+u912*u2855)
  value(8)=value(8)+(c(52)*(x*((u383+u2856*(u2856*(u529+u3777)+u961))*u&
 &2857+u2856*(u463+u3803+u2856*(u723*u2856+u963+u514))+u72+u676)))
  value(9)=value(9)+(c(52)*((u2856*(u378+u2855*(u3938+u3982)+u2856*(u28&
 &56*(u471+u467)+u2855*(u789+u763)+u4217))+u2855*(u4167+u3939)+u3815)*z&
 &))
  u676=-4.48221828871296415d-1*u640
  u463=8.9644365774259283d-1*u3946
  u514=-1.49407276290432138d-1*u955
  u519=1.49407276290432138d-1*u3946
  u683=u2855*(u635*u2855+u424)
  u772=u2855*(u519+u870*u2855)
  value(1)=value(1)+(c(53)*(x*((u521+u2856*(u2856*(u728+u3805)+u408))*u&
 &2857+u2856*(u463+u683+u2856*((u837+u4030)*u2856+u811+u514))+u772+u676&
 &)))
  u676=-9.50822084012189831d-1*u640
  u463=2.11293796447153296d-1*u703
  u86=-u3825
  u4136=-u965
  u40=-3.16940694670729944d-1*u955
  u459=2.11293796447153296d-1*u4194
  u91=-1.40862530964768864d-1*u700
  u794=u2856*(u647+u3777)
  u4135=u91*u2856
  value(2)=value(2)+(c(53)*(z*(u2857*(u463+u2856*(u2856*(u459+u4135)+u4&
 &136)+(u794+u951)*u2857)+u2856*(u86+u2856*(u698*u2856+u40))+u676)))
  u463=u2855*(u874*u2855+u666)
  u86=u2855*(u3978+u625*u2855)
  value(3)=value(3)+(c(53)*(x*(u2857*(u989+u2856*(u2856*(u411+u3807)+u8&
 &54)+(u2856*(u481+u3806)+u4072)*u2857)+u2856*(u932+u463+u2856*((u419+u&
 &3973)*u2856+u62+u804))+u86+u954)))
  u4136=u441*u2856
  value(4)=value(4)+(c(53)*(z*(u2857*(u449+u2856*(u78*u2856+u806)+(u285&
 &6*(u3761+u4136)+u4028)*u2857)+u2856*(u4140+u658*u2856)+u4146)))
  value(5)=value(5)+(c(53)*(y*(u2857*(u4189+u2856*(u2856*(u486+u4003)+u&
 &805)+u2857*((u482+u4132)*u2857+u77*u2856+u551))+u2856*(u387+u2855*(u2&
 &855*(u861+u3829)+u484)+u2856*((u694+u3770)*u2856+u3932+u485))+u2855*(&
 &u3991+u2855*(u3859+u483))+u898)))
  value(6)=value(6)+(c(53)*u4026)
  u47=-4.23529818351428947d-1*u559
  u459=2.82353212234285965d-2*u3785
  u40=-1.69411927340571579d-1*u955
  u4138=5.6470642446857193d-1*u3970
  u3898=-2.82353212234285965d-2*u4204
  u788=5.6470642446857193d-2*u915
  u4026=4.8000046079828614d-1*u948
  u4140=-8.47059636702857895d-2*u897
  u4141=2.82353212234285965d-2*u952
  u3905=1.12941284893714386d-1*u515
  u515=-9.41177374114286549d-3*u701
  u4137=u84*u2855
  value(7)=value(7)+(c(53)*(y*(u2857*(u459+u2856*(u2856*(u4026+u779)+u3&
 &898)+u2857*((u838+u4007)*u2857+u2856*(u788+u3804)+u40))+u2856*(u4138+&
 &u2855*(u2855*(u515+u3769)+u4141)+u2856*(u4115+u63+u930))+u2855*(u4140&
 &+u2855*(u4137+u3905))+u47)))
  u4138=3.16940694670729944d-1*u97
  u459=-1.05646898223576648d-1*u3948
  u40=1.40862530964768864d-1*u857
  value(8)=value(8)+(c(53)*(u528*((u615+u2856*(u3777+u561))*u2857+u2856&
 &*(u459+u2855*(u4024+u40)+(u4024+u698)*u2856)+u2855*(u3971+u3853)+u413&
 &8)))
  u4138=u696*u2855
  value(9)=value(9)+(c(53)*(y*((u6+u2856*(u2856*(u4218+u3776)+u525))*u2&
 &857+u2856*(u4139+u2855*(u2855*(u79+u3771)+u809)+u2856*((u488+u756)*u2&
 &856+u2855*(u80+u3856)+u810))+u2855*(u4142+u2855*(u4138+u808))+u856)))
  u459=4.48221828871296415d-1*u3946
  u40=1.99209701720576184d-1*u857
  value(1)=value(1)+(c(54)*(u528*((u595+u2856*(u3805+u681))*u2857+u2856&
 &*(u514+u2855*(u4030+u40)+(u4030+u635)*u2856)+u2855*(u872+u3940)+u459)&
 &))
  u459=3.16940694670729944d-1*u868
  u40=1.05646898223576648d-1*u524
  u809=-1.05646898223576648d-1*u972
  u810=3.5215632741192216d-2*u379
  u4141=-4.22587592894306592d-1*u3808
  u4140=8.45175185788613183d-1*u948
  value(2)=value(2)+(c(54)*(y*(u2857*(u459+u2856*(u4140*u2856+u809)+u28&
 &57*((u648+u3777)*u2857+u2856*(u810+u4135)+u503))+u2856*(u40+u4141*u28&
 &56)+u676)))
  value(3)=value(3)+(c(54)*(u528*(u2857*(u3830+u2856*(u3807+u499)+(u380&
 &6+u3871)*u2857)+u2856*(u829+u2855*(u3861+u81))+u2855*(u624+u3852)+u69&
 &9)))
  value(4)=value(4)+(c(54)*(y*(u2857*(u4056+u606*u2856+u2857*((u501+u41&
 &36)*u2857+u3841*u2856+u831))+u966*u2856+u4146)))
  u459=(7.55d2+u3912)
  u40=u3838*(2.4655926003719755d7*u534+u3788*u459)*u4147
  u809=+u608*(1.54956779983642699d8*u534+u3788*(3.65d2*u4031+u593))*u41&
 &47
  u3898=2.52544390488515477d-2*u948
  u47=+u4045*(u3820+u3788*(u881+u3891))*u4147
  u4140=u3838*(7.32329987593927823d8*u534+u3788*(1.15d2*u526+u727))*u41&
 &47
  u4141=+u4045*(1.45323007571593258d7*u534+u3788*(4.45d2+6.8d1*pd2))*u4&
 &147
  u810=u3838*(u427+u3788*(u4012+2.96d2*u4011))*u4147
  u515=+u3875*(6.68649119107499317d8*u534+u3788*(1.05d2*u526+u3844))*u4&
 &147
  u788=u975*(u610-u3788*(u4011+u907))*u4147
  u4026=+u608*(-u3788*(u645-5.5d1*u3894)+u4199)*u4147
  u4142=4.41952683354902085d-2*u700
  u4139=u49*u2856
  value(5)=value(5)+(c(54)*(z*(u2857*(u382+u2855*(u2855*(u4026+u3829)+u&
 &810)+u2856*(u2856*(u4140+u3858)+u809)+u2857*((u500+u3923)*u2857+u2856&
 &*(u3898+u4139)+u2855*(u515+u4142*u2855)+u830))+u2856*(u40+u2856*(u410&
 &4+u47))+u2855*(u4141+u2855*(u3859+u788))+u865)))
  value(6)=value(6)+(c(54)*u478)
  u40=-7.62353673032572105d-1*u855
  u809=5.6470642446857193d-2*u4175
  u4026=-2.82353212234285965d-2*u96
  u810=3.81176836516286053d0*u376
  u4140=-2.96470872846000263d0*u3808
  u4141=-4.23529818351428947d-1*u3808
  u788=4.23529818351428947d-1*u948
  u515=8.47059636702857895d-2*u659
  u47=-2.82353212234285965d-2*u983
  u3898=1.12941284893714386d-1*u869
  u4142=-5.6470642446857193d-2*u724
  u676=9.41177374114286549d-3*u44
  u3940=-4.70588687057143275d-2*u700
  value(7)=value(7)+(c(54)*(z*(u2857*(u809+u2855*(u2855*(u676+u3769)+u4&
 &7)+u2856*(u2856*(u788+u779)+u4140)+u2857*((u624+u3772)*u2857+u3922+u2&
 &855*(u3898+u3940*u2855)+u4026))+u2856*(u810+u2856*(u4116+u4141))+u285&
 &5*(u515+u2855*(u4137+u4142))+u40)))
  u4140=-3.16940694670729944d-1*u562
  u4141=3.16940694670729944d-1*u520
  u676=-u384
  u3898=-u705
  value(8)=value(8)+(c(54)*(x*(u2857*(u4141+u3803+u2856*(u647*u2856+u38&
 &98)+u2857*((u614+u4024)*u2857+u794+u963+u476))+u2856*(u676+u951*u2856&
 &)+u72+u4140)))
  u4140=u986*u2856
  u4141=u4179*u2856
  value(9)=value(9)+(c(54)*(z*(u2857*(u586+u2855*(u2855*(u964+u3771)+u5&
 &02)+u2856*(u2856*(u89+u4140)+u834)+(u2856*(u417+u4141)+u2855*(u832+u3&
 &771)+u576)*u2857)+u2856*(u567+u2856*(u4218*u2856+u835))+u2855*(u517+u&
 &2855*(u4138+u833))+u3930)))
  u676=-4.48221828871296415d-1*u376
  u3898=4.48221828871296415d-1*u851
  u809=-1.49407276290432138d-1*u958
  u4142=-u4108
  value(1)=value(1)+(c(55)*(x*(u2857*(u3898+u683+u2856*(u635*u2856+u414&
 &2)+u2857*(u3870*u2857+u2856*(u635+u3805)+u811+u809))+u2856*(u4167+u87&
 &0*u2856)+u772+u676)))
  u4142=5.28234491117883239d-1*u3884
  u3898=-1.05646898223576648d-1*u955
  u809=-u3883
  value(2)=value(2)+(c(55)*(z*(u2857*(u4142+u839+u2856*(u4117+u809)+u28&
 &57*(u3834*u2857+u737+u505+u3898))+u2856*(u407+u4118)+u733+u866)))
  value(3)=value(3)+(c(55)*(x*(u2857*(u985+u463+u2856*(u633*u2856+u848)&
 &+u2857*((u674+u3973+u3806)*u2857+u2856*(u506+u3807)+u62+u841))+u2856*&
 &(u819+u931*u2856)+u86+u781)))
  value(4)=value(4)+(c(55)*(z*(u2857*(u667+u3925+u2856*(u4119+u803)+u28&
 &57*((u510+u445)*u2857+u2856*(u425+u4095)+u538+u843))+u2856*(u862+u412&
 &1)+u777+u522)))
  value(5)=value(5)+(c(55)*(y*(u2857*(u518+u2855*(u4133+u970)+u2856*(u7&
 &41*u2856+u509)+u2857*(u2857*(u508+u3768+u4139+u3944)+u2856*(u642+u385&
 &8)+u2855*(u842+u3896)+u507))+u2856*(u758+u4097*u2856)+u2855*(u984+u41&
 &34)+u738)))
  value(6)=value(6)+(c(55)*u496)
  u4142=-u4114
  u947=-u3924
  u3898=-u3942
  u676=5.6470642446857193d-1*u376
  u4137=7.5294189929142924d-1*u948
  u810=-6.77647709362286316d-1*u376
  u970=-u4113
  u515=-7.90588994256000702d-1*u948
  value(7)=value(7)+(c(55)*(y*(u2857*(u947+u2855*(u3855+u970)+u2856*(u3&
 &943+u552)+u2857*((u633+u492)*u2857+u2856*(u4137+u779)+u2855*(u515+u40&
 &65)+u3898))+u2856*(u676+u4122)+u2855*(u810+u3945)+u4142)))
  u4142=u615*u2856
  value(8)=value(8)+(c(55)*(u528*(u2857*(u809+u3921+u561*u2856+(u4041+u&
 &623)*u2857)+u4142+u4129+u407)))
  value(9)=value(9)+(c(55)*(y*(u2857*(u474+u2855*(u3938+u3947)+u2856*(u&
 &512*u2856+u845)+u2857*((u511+u4141)*u2857+u2856*(u90+u4140)+u609+u844&
 &))+u2856*(u564+u849*u2856)+u2855*(u3779+u3939)+u853)))
  value(1)=value(1)+(c(56)*(u528*(u2857*(u3935+u681*u2856+(u464)*u2857)&
 &+u595*u2856+u3936)))
  value(2)=value(2)+(c(56)*(y*(u2857*(u4123+u4142+u2857*((u767)*u2857+u&
 &836*u2856+u4124))+u712*u2856+u4125)))
  u4137=-1.01647156404342947d1*u376
  u810=1.07294220649028667d1*u3808
  u970=-2.25882569787428772d0*u948
  value(3)=value(3)+(c(56)*(u528*(u2857*(u810+u4130+u684*u2856+u2857*(u&
 &847*u2857+u4083+u970))+u913*u2856+u4131+u4137)))
  u4137=9.48373475492457397d3*u702
  u810=-1.51736942006107725d1*u376
  u970=8.78477032666939461d0*u3808
  u515=-1.3310258070711204d0*u948
  u947=(u2857*(u810+u4126+u601*u2856+u2857*(u2857*(u515+u4038+u4120)+u6&
 &57*u2856+u4127+u970))+u4064*u2856+u4128+u4137)
  value(4)=value(4)+(c(56)*(y*u947))
  u3898=u454*(6.23449062104894468d4*u534-u3788)*u4147
  u676=-6.31360976221288694d-2*u2858
  u552=+u676*(-u3788*(u4019-1.47d2)+4.8005577782076874d6*u534)*u4147
  u809=u4032*(u59-u3788*(u4006-3.15d2*u4001))*u4147
  u3936=-4.92461561452605181d-1*u948
  u4026=4.7352073216596652d0*u376
  u40=-4.92461561452605181d0*u3808
  u3940=1.02280478147848768d0*u948
  u4097=9.4704146433193304d-2*u3808
  u788=-6.31360976221288694d-2*u948
  u4119=3.78816585732773216d-1*u2858
  u47=u4119*(6.20480257047252114d5*u534+u3788*(1.9d1+pd2))*u4147
  u4116=+u676*(4.71238426799570947d7*u534+u3788*(1.11d2*u4031+u3891))*u&
 &4147
  u474=u3840*(u770+u3788*(u893+u592))*u4147
  u758=-9.4704146433193304d-2*u3808
  u421=6.31360976221288694d-2*u948
  value(5)=value(5)+(c(56)*(z*(u2857*(u552+u2855*(u421*u2855+u4116)+u28&
 &56*(u788*u2856+u40)+u2857*(u2857*(u3936+u497)+u2856*(u3940+u4003)+u28&
 &55*(u474+u3829)+u809))+u2856*(u4026+u4097*u2856)+u2855*(u47+u758*u285&
 &5)+u3898)))
  value(6)=value(6)+(c(56)*(x*u947))
  u4116=5.08235782021714737d0*u376
  u47=-5.36471103245143333d0*u3808
  u758=+u47
  u4026=-u941
  u40=1.41176606117142983d-1*u3808
  u3940=-u4116
  u4097=-u47
  u47=-1.41176606117142983d-1*u3808
  value(7)=value(7)+(c(56)*(z*(u2857*(u2855*(u4103+u4097)+u2856*(u504*u&
 &2856+u758)+u2857*((u492)*u2857+u2856*(u4026+u779)+u2855*(u413+u4065))&
 &)+u2856*(u4116+u40*u2856)+u2855*(u3940+u47*u2855))))
  value(8)=value(8)+(c(56)*(x*(u2857*(u3937+u961*u2856+u2857*((u3904+u3&
 &777)*u2857+u529*u2856+u479*u2855))+u383*u2856+u93*u2855)))
  u758=1.12055457217824104d0*u376
  u4026=-3.73518190726080346d-1*u3884
  u40=7.47036381452160691d-2*u955
  u3940=-u758
  u4097=-2.24110914435648207d0*u98
  u47=3.73518190726080346d-1*u952
  u4116=-2.98814552580864276d-1*u95
  u474=-8.71542445027520807d-1*u948
  value(9)=value(9)+(c(56)*(z*(u2857*(u4026+u2855*(u474*u2855+u47)+u285&
 &6*(u398*u2856+u4008)+u2857*((u882+u756)*u2857+u790+u2855*(u4116+u3934&
 &)+u40))+u2856*(u3940+u576*u2856)+u2855*(u4097+u525*u2855)+u758)))
  if ( lmax .eq. 5 ) return
  u4026=A11_(p,pd,erfpd,exppd2)
  u40=p**11
  u3940=u2858*u4026*u40
  u4097=1.02462608902408428d6*u3940
  u47=4.24539123242856709d5*u4026
  u4116=u2858*(u47+u621)*u40
  u474=5.22925467016512484d0*u4116
  u758=-4.70632920314861236d1*u4116
  u421=+u758
  u809=6.36808684864285064d6*u4026
  u3936=u2858*(u809-2.26874092915906051d2*u4021)*u40
  u747=-5.22925467016512484d0*u3936
  u4137=+u747
  u970=1.56877640104953745d1*u3936
  u810=1.08257476426928461d8*u4026
  u3898=(u926+1.5d1)
  u552=1.7d1*u3898
  u788=(u552+u4006)
  u947=u2858*(u810+u621*u788)*u40
  u803=1.04585093403302497d0*u947
  u44=-1.64348003919475352d0*u947
  u790=+u44
  u3782=2.05689205211164076d9*u4026
  u3803=(u4006+1.7d1*u4014)
  u72=1.9d1*u3803
  u4108=u2858*(u3782+u621*(u72+u645))*u40
  u4114=-4.98024254301440461d-2*u4108
  u948=4.98024254301440461d-2*u4108
  u4056=u948*u2855
  u966=u4114*u2855
  value(1)=value(1)+(c(57)*(u2983*((u474+u2855*(u2855*(u803+u966)+u4137&
 &))*u2856+u2855*(u421+u2855*(u2855*(u790+u4056)+u970))+u4097)))
  u421=3.10508595312677627d5*u3940
  u790=5.28234491117883239d-1*u4116
  u3883=-3.9617586833841243d1*u4116
  u829=+u3883
  u3852=-1.58470347335364972d0*u3936
  u514=+u3852
  u478=2.06011451535974463d1*u3936
  u667=5.28234491117883239d-1*u947
  u496=-2.85246625203656949d0*u947
  u767=+u496
  u4=-3.5215632741192216d-2*u4108
  u3934=1.05646898223576648d-1*u4108
  u724=u3934*u2855
  u601=u4*u2855
  value(2)=value(2)+(c(57)*(y*((u790+u2855*(u2855*(u667+u601)+u514))*u2&
 &856+u2855*(u829+u2855*(u2855*(u767+u724)+u478))+u421)*z))
  u829=-3.87272259769493558d5*u3940
  u767=-1.18588349138400105d1*u4116
  u4118=1.97647248564000175d0*u4116
  u683=1.77882523707600158d1*u4116
  u3884=1.18588349138400105d1*u3936
  u497=-1.97647248564000175d0*u3936
  u839=+u497
  u985=-5.92941745692000526d0*u3936
  u4031=+u985
  u516=-2.3717669827680021d0*u947
  u638=3.95294497128000351d-1*u947
  u847=6.21177066915429123d-1*u947
  u505=1.12941284893714386d-1*u4108
  u3922=-1.8823547482285731d-2*u4108
  u760=u3922*u2855
  u4008=u505*u2855
  value(3)=value(3)+(c(57)*(u2983*((u767+u2855*(u2855*(u516+u4008)+u388&
 &4))*u2857+(u4118+u2855*(u2855*(u638+u760)+u839))*u2856+u2855*(u683+u2&
 &855*(u2855*(u847+u760)+u4031))+u829)))
  u3905=-1.17361217592191603d5*u3940
  u622=-7.98615484242672237d-1*u4116
  u41=5.98961613182004178d-1*u4116
  u513=1.49740403295501045d1*u4116
  u825=2.39584645272801671d0*u3936
  u733=-1.79688483954601253d0*u3936
  u4103=+u733
  u49=-7.78650097136605431d0*u3936
  u811=-7.98615484242672237d-1*u947
  u625=5.98961613182004178d-1*u947
  u639=1.07813090372760752d0*u947
  u4139=5.32410322828448158d-2*u4108
  u4021=-3.99307742121336118d-2*u4108
  u933=u4021*u2855
  u4171=u2855*(u639+u933)
  u3778=u4139*u2855
  u3925=(u2984*((u622+u2855*(u2855*(u811+u3778)+u825))*u2857+(u41+u2855&
 &*(u2855*(u625+u933)+u4103))*u2856+u2855*(u513+u2855*(u4171+u49))+u390&
 &5))
  value(4)=value(4)+(c(57)*u3925)
  u3930=7.6199329812820435d4*u4026
  u3932=-4.13254457163025327d-1*u2858*(u599+u3930)*u40
  u88=-1.23709585520649233d5*u3940
  u4007=-2.52544390488515477d-1*u4116
  u3771=u4119*(1.01236252465604292d6*u4026+u629)*u40
  u4119=-8.52337317898739736d-1*u4116
  u4138=1.13644975719831965d0*u2858*(1.2953886068179474d6*u4026+u629)*u&
 &40
  u417=1.66679297722420215d1*u4116
  u485=7.57633171465546432d-1*u3936
  u777=-1.89408292866386608d-1*u2858
  u955=+u777*(9.38231462366713327d7*u4026+u621*(2.21d2+u4051))*u40
  u4001=2.55701195369621921d0*u3936
  u772=-9.4704146433193304d-2*u2858
  u90=+u772*(1.15899180645299882d8*u4026+u621*(2.73d2+u3888))*u40
  u763=-9.09159805758655719d0*u3936
  u517=-2.52544390488515477d-1*u947
  u848=u454*(6.17704424318356512d8*u4026+u621*(9.7d1*u3898+u3786))*u40
  u518=-8.52337317898739736d-1*u947
  u849=u3838*(4.59139061787149531d9*u4026+u621*(7.21d2*u3898+u4035))*u4&
 &0
  u850=1.31323083054028048d0*u947
  u413=1.68362926992343652d-2*u4108
  u92=+u3875*(1.8078998563297053d10*u4026+u621*(1.67d2*u3803+u3844))*u4&
 &0
  u504=5.68224878599159824d-2*u4108
  u519=-2.96739658824005686d-1*u947
  u3806=-5.05088780977030955d-2*u4108
  u4003=6.31360976221288694d-2*u4108
  u3937=6.31360976221288694d-3*u4108
  u4028=u504*u2855
  u3830=u3937*u2855
  u3997=u4003*u2855
  u3808=u3806*u2855
  u3946=u413*u2855
  value(5)=value(5)+(c(57)*(u2857*(u88+u2855*(u2855*(u763+u2855*(u3808+&
 &u850))+u417)+(u2855*(u485+u2855*(u3946+u517))+u4007)*u2857)+u2856*(u3&
 &771+u2855*(u2855*(u848+u2855*(u3997+u92))+u955)+(u2855*(u4001+u2855*(&
 &u4028+u518))+u4119)*u2856)+u2855*(u4138+u2855*(u2855*(u849+u2855*(u38&
 &30+u519))+u90))+u3932))
  u4135=-8.2152852314534122d5*u3940
  u4121=-5.59030838969870566d0*u4116
  u794=4.19273129227402924d0*u4116
  u502=3.77345816304662632d1*u4116
  u4041=5.59030838969870566d0*u3936
  u463=-4.19273129227402924d0*u3936
  u842=-1.25781938768220877d1*u3936
  u520=-1.11806167793974113d0*u947
  u659=8.38546258454805848d-1*u947
  u851=1.31771554900040919d0*u947
  value(6)=value(6)+(c(57)*(u2985*((u4121+u2855*(u2855*(u520+u3778)+u40&
 &41))*u2857+(u794+u2855*(u2855*(u659+u933)+u463))*u2856+u2855*(u502+u2&
 &855*(u2855*(u851+u933)+u842))+u4135)))
  u737=2.28597989438461305d5*u4026
  u501=u2858*(u737-u722)*u40
  u963=-7.70054215184416268d-2*u501
  u86=1.93636129884746779d5*u3940
  u464=u2858*(u737+u3992)*u40
  u87=1.69411927340571579d0*u464
  u781=-9.88236242820000877d-1*u4116
  u486=u2858*(u3930+u3992)*u40
  u3930=5.08235782021714737d0*u486
  u75=-2.37176698276800211d1*u4116
  u492=+u75
  u3998=2.08024170388999788d7*u4026
  u3939=u2858*(u3998+u621*u755)*u40
  u77=-8.47059636702857895d-1*u3939
  u445=2.96470872846000263d0*u3936
  u3924=2.97177386269999696d6*u4026
  u4038=-u621*(u4019-7.0d0)
  u3973=4.23529818351428947d-1*u2858
  u508=u3973*(u3924+u4038)*u40
  u964=u2858*(1.33729823821499863d8*u4026+u621*(2.1d1*u3898+u4006))*u40
  u755=8.47059636702857895d-1*u964
  u4083=-9.88236242820000877d-1*u947
  u658=5.79495903226499408d8*u4026
  u634=-u621*(u3891-9.1d1*u4014)
  u819=+u3758*(u634+u658)*u40
  u738=-1.5811779885120014d0*u947
  u820=+u738
  u3864=3.78901167494249613d9*u4026
  u693=3.5d1*u3803
  u3896=u2858*(u3864+u621*(u693+u645))*u40
  u93=-5.6470642446857193d-2*u3896
  u3788=6.58824161880000585d-2*u4108
  u84=3.10588533457714561d-1*u947
  u91=5.6470642446857193d-2*u4108
  u6=-9.41177374114286549d-3*u4108
  u765=u6*u2855
  u644=u91*u2855
  u3779=u3788*u2855
  value(7)=value(7)+(c(57)*((u86+u2855*(u2855*(u3884+u2855*(u644+u820))&
 &+u492))*u2857+u2856*(u87+u2855*(u2855*(u755+u2855*(u644+u93))+u77)+(u&
 &2855*(u445+u2855*(u3779+u4083))+u781)*u2856)+u2855*(u3930+u2855*(u285&
 &5*(u819+u2855*(u765+u84))+u508))+u963))
  u963=7.2452005572958113d5*u3940
  u781=1.1092924313475548d1*u4116
  u492=-3.32787729404266441d1*u4116
  u77=+u492
  u819=-1.1092924313475548d1*u3936
  u755=+u819
  u4083=-u819
  u819=2.21858486269510961d0*u947
  u820=-1.16211588045934313d0*u947
  u84=+u820
  u508=-1.05646898223576648d-1*u4108
  u952=3.5215632741192216d-2*u4108
  u457=u508*u2855
  u449=u952*u2855
  value(8)=value(8)+(c(57)*(x*((u781+u2855*(u2855*(u819+u457)+u755))*u2&
 &856+u2855*(u77+u2855*(u2855*(u84+u449)+u4083))+u963)*z))
  u77=3.26568556340659007d4*u4026
  u84=u2858*(-u3757+u77)*u40
  u696=-4.07474389882996741d-1*u84
  u503=-2.19562733362303775d5*u3940
  u507=-1.86759095363040173d-1*u4116
  u419=9.51438511236649692d5*u3940
  u832=2.8013864304456026d1*u4116
  u3942=5.60277286089120519d-1*u3936
  u414=-2.40919233018321823d1*u4116
  u4113=-1.45672094383171335d1*u3936
  u4072=-1.86759095363040173d-1*u947
  u4120=5.78953195625424536d0*u3936
  u852=2.01699822992083387d0*u947
  u3881=1.24506063575360115d-2*u4108
  u521=-4.85573647943904449d-1*u947
  u3805=-7.47036381452160691d-2*u4108
  u4064=u3881*u2855
  u3835=u3805*u2855
  value(9)=value(9)+(c(57)*(u2856*(u503+u2855*(u2855*(u4113+u2855*(u383&
 &5+u852))+u832)+(u2855*(u3942+u2855*(u4064+u4072))+u507)*u2856)+u2855*&
 &(u419+u2855*(u2855*(u4120+u2855*(u4064+u521))+u414))+u696))
  u58=-2.71649593255331161d-1*u84
  u84=2.43958592624781972d4*u3940
  u506=7.47036381452160692d-1*u4116
  u85=4.63521325987085747d5*u3940
  u3931=-7.47036381452160692d0*u4116
  u909=+u3931
  u3763=-2.24110914435648207d0*u3936
  u831=+u3763
  u4167=-8.21740019597376761d0*u4116
  u4025=+u4167
  u420=5.97629105161728553d0*u3936
  u668=7.47036381452160692d-1*u947
  u494=1.24506063575360115d0*u3936
  u495=-1.09565335946316901d0*u947
  u3758=+u495
  u426=-4.98024254301440461d-2*u947
  u4143=u426*u2855
  value(1)=value(1)+(c(58)*(u2856*(u84+u2855*(u2855*(u420+u2855*(u4056+&
 &u3758))+u909)+(u2855*(u831+u2855*(u966+u668))+u506)*u2856)+u2855*(u85&
 &+u2855*(u2855*(u494+u4143)+u4025))+u58))
  u58=-7.9235173667682486d0*u4116
  u909=+u58
  u3758=-5.28234491117883239d-1*u3936
  u686=-1.00364553312397816d1*u4116
  u4136=+u686
  u739=8.97998634900401507d0*u3936
  u650=3.5215632741192216d-1*u947
  u821=2.11293796447153296d0*u3936
  u822=-2.00729106624795631d0*u947
  u3818=+u822
  u467=-1.05646898223576648d-1*u947
  u4094=(u2855*(u650+u601)+u3758)
  value(2)=value(2)+(c(58)*(x*(u2856*(u909+u2855*(u2855*(u3818+u724)+u7&
 &39)+u4094*u2856)+u2855*(u4136+u2855*(u467*u2855+u821))+u421)*z))
  u909=1.14298994719230652d6*u4026
  u4136=u2858*(u909-u620)*u40
  u3818=5.13369476789610845d-2*u4136
  u510=-8.47059636702857895d-1*u2858*(9.90591287566665655d5*u4026+u629)&
 &*u40
  u56=-2.54117891010857368d0*u2858*(9.39791734358118698d5*u4026+u629)*u&
 &40
  u4126=5.6470642446857193d-1*u2858*(6.83507988420999302d7*u4026+u621*(&
 &1.61d2+u3790))*u40
  u515=u447*(5.05201556658999484d7*u4026+u621*(1.19d2+u4015))*u40
  u712=2.22883039702499772d8*u4026
  u726=3.5d1*u3898
  u984=u2858*(u712+u621*(u726+u3833))*u40
  u522=-1.12941284893714386d0*u984
  u853=1.97647248564000175d0*u947
  u523=+u3826*(1.38187484615549859d9*u4026+u621*(2.17d2*u3898+u600))*u4&
 &0
  u805=1.28826396948044868d10*u4026
  u785=1.19d2*u3803
  u741=3.7647094964571462d-2*u2858*(u805+u621*(u785+u603))*u40
  u695=-1.31764832376000117d-1*u4108
  u854=1.31764832376000117d-1*u947
  u3991=u695*u2855
  u4189=u2855*(u3991+u741)
  value(3)=value(3)+(c(58)*(u2856*(u510+u2855*(u2855*(u522+u4189)+u4126&
 &)+(u2855*(u4031+u2855*(u3991+u853))+u4118)*u2856)+u2855*(u56+u2855*(u&
 &2855*(u523+u854*u2855)+u515))+u3818))
  u3768=-1.19792322636400836d0*u2858*(9.47048813387911121d5*u4026+u629)&
 &*u40
  u8=u2858*(7.21716509512856406d6*u4026-2.26874092915906051d2*u759)*u40
  u657=7.98615484242672237d-1*u8
  u694=-7.98615484242672237d-1*u3936
  u614=7.78650097136605431d0*u4116
  u493=-5.98961613182004178d-1*u3936
  u564=1.99653871060668059d-1*u2858
  u509=u564*(6.41054076096713631d7*u4026+u621*(1.51d2+u3888))*u40
  u3807=2.61091560794356876d8*u4026
  u845=4.1d1*u3898
  u951=u2858*(u3807+u621*(u845+u3891))*u40
  u376=-2.66205161414224079d-1*u951
  u699=5.32410322828448158d-1*u947
  u3974=-6.58857774500204595d0*u3936
  u640=3.99307742121336118d-1*u947
  u759=u2858*(u809+u2966*(8.0d0*u3898+u4011))*u40
  u700=-5.11113909915310231d0*u759
  u500=u2858*(3.35598176923478229d9*u4026+u621*(3.1d1*u3803+u645))*u40
  u4128=5.32410322828448158d-2*u500
  u7=-5.32410322828448158d-2*u4108
  u701=1.99653871060668059d-1*u947
  u932=u7*u2855
  u3964=(u2985*(u2857*(u657+u2855*(u2855*(u4128+u932)+u376)+(u2855*(u69&
 &9+u932)+u694)*u2857)+u2856*(u614+u2855*(u4171+u3974)+(u2855*(u640+u93&
 &3)+u493)*u2856)+u2855*(u509+u2855*(u701*u2855+u700))+u3768))
  value(4)=value(4)+(c(58)*u3964)
  u764=u2858*(4.89852834510988511d5*u4026+u621)*u40
  u833=2.2728995143966393d0*u764
  u834=3.78816585732773216d0*u4116
  u4127=2.52544390488515477d-1*u3936
  u488=2.12269561621428355d6*u4026
  u633=u2858*(u488-1.13437046457953026d2*u3893)*u40
  u3776=-3.03053268586218573d0*u633
  u729=8.52337317898739736d-1*u3936
  u4193=-2.52544390488515477d-1*u2858*(5.30673904053570887d7*u4026+u621&
 &*(1.25d2+u4015))*u40
  u3975=-4.29325463830476312d0*u3936
  u377=-1.68362926992343652d-1*u947
  u559=1.24177693548535587d9*u4026
  u562=1.95d2*u3898
  u855=u4075*(u559+u621*(u562+u4035))*u40
  u524=-5.68224878599159824d-1*u947
  u856=u3838*(6.25982937221592218d9*u4026+u621*(9.83d2*u3898+u977))*u40
  u857=9.59668683856358814d-1*u947
  u94=+u4013*(8.11931073201963457d9*u4026+u621*(7.5d1*u3803+u727))*u40
  u378=-3.66189366208347442d-1*u947
  u930=u377+u3946
  value(5)=value(5)+(c(58)*(u2983*(u2857*(u834+u2855*(u2855*(u857+u3808&
 &)+u3975)+(u2855*(u930)+u4127)*u2857)+u2856*(u3776+u2855*(u2855*(u94+u&
 &3997)+u855)+(u2855*(u524+u4028)+u729)*u2856)+u2855*(u4193+u2855*(u285&
 &5*(u378+u3830)+u856))+u833)))
  value(6)=value(6)+(c(58)*u3925)
  u3925=u2858*(3.80996649064102175d5*u4026+u621)*u40
  u823=2.54117891010857368d0*u3925
  u824=-5.92941745692000526d0*u4116
  u43=+u824
  u3834=u2858*(u3924-5.67185232289765129d1*u60)*u40
  u3923=-2.25882569787428772d0*u3834
  u915=+u3923
  u3929=9.88236242820000877d-1*u3936
  u806=3.26895124896999666d7*u4026
  u4129=u2858*(u806+u621*u3787)*u40
  u814=-2.82353212234285965d-1*u4129
  u483=-u985
  u985=u2858*(8.46955550869499135d8*u4026+u621*(1.33d2*u3898+u600))*u40
  u81=9.4117737411428655d-2*u985
  u779=-6.58824161880000585d-1*u947
  u897=4.90342687345499499d8*u4026
  u441=7.7d1*u3898
  u498=u2858*(u897+u621*(u441+u600))*u40
  u4194=2.82353212234285965d-2*u498
  u4192=-1.18588349138400105d0*u947
  u4061=+u4192
  u3853=5.30461634491949458d9*u4026
  u685=4.9d1*u3803
  u3756=u2858*(u3853+u621*(u685+u440))*u40
  u95=-3.7647094964571462d-2*u3756
  u702=1.12941284893714386d-1*u947
  value(7)=value(7)+(c(58)*(u2983*((u43+u2855*(u2855*(u4061+u644)+u483)&
 &)*u2857+u2856*(u915+u2855*(u2855*(u95+u644)+u81)+(u2855*(u779+u3779)+&
 &u3929)*u2856)+u2855*(u814+u2855*(u2855*(u702+u765)+u4194))+u823)))
  u43=-3.45009550347419585d4*u3940
  u915=+u43
  u3929=1.58470347335364972d0*u4116
  u814=-5.28234491117883239d-1*u4116
  u81=-4.75411042006094916d0*u3936
  u779=+u81
  u4194=2.6411724555894162d0*u3936
  u4061=1.58470347335364972d0*u947
  u3861=-6.69097022082652103d-1*u947
  value(8)=value(8)+(c(58)*(y*((u3929+u2855*(u2855*(u4061+u457)+u779))*&
 &u2856+u2855*(u814+u2855*(u2855*(u3861+u449)+u4194))+u915)*z))
  u3861=-2.92750311149738367d5*u3940
  u472=5.97629105161728553d0*u4116
  u957=1.86759095363040173d-1*u3936
  u3777=-6.59882136949408611d0*u3936
  u3829=-1.24506063575360115d-1*u947
  u858=1.44427033747417734d0*u947
  u610=-2.24110914435648207d-1*u947
  u606=u2855*(u3829+u4064)
  u651=(u606+u957)*u2856
  value(9)=value(9)+(c(58)*(u2983*(u2856*(u472+u2855*(u2855*(u858+u3835&
 &)+u3777)+u651)+u2855*(u472+u2855*(u2855*(u610+u4064)+u957))+u3861)))
  u4102=9.75834370499127889d4*u3940
  u491=-1.41936912475910531d1*u4116
  u411=+u491
  u687=8.21740019597376761d0*u3936
  u815=-1.24506063575360115d0*u947
  u808=+u815
  value(1)=value(1)+(c(59)*(y*((u506+u2855*(u2855*(u668+u966)+u831))*u2&
 &856+u2855*(u411+u2855*(u2855*(u808+u4056)+u687))+u4102)*z))
  u411=6.20480257047252114d5*u4026
  u808=u2858*(u411+u621)*u40
  u816=1.58470347335364972d0*u808
  u3858=-u781
  u682=-5.28234491117883239d-1*u8
  u69=5.28234491117883239d-1*u3936
  u481=5.51900860215713722d6*u4026
  u4095=u2858*(u481+u621*(1.3d1+pd2))*u40
  u817=-2.11293796447153296d0*u4095
  u972=+u817
  u743=1.7607816370596108d-1*u951
  u4088=-3.5215632741192216d-1*u947
  u3787=u2858*(2.9930008188621398d8*u4026+u621*(4.7d1*u3898+u4006))*u40
  u732=1.05646898223576648d-1*u3787
  u427=-u819
  u59=-3.5215632741192216d-2*u500
  u407=-2.11293796447153296d-1*u947
  u3947=u407*u2855
  value(2)=value(2)+(c(59)*(u2983*((u3858+u2855*(u2855*(u427+u724)+u408&
 &3))*u2857+u2856*(u682+u2855*(u2855*(u59+u449)+u743)+(u2855*(u4088+u44&
 &9)+u69)*u2856)+u2855*(u972+u2855*(u3947+u732))+u816)))
  u972=7.3766144717998773d4*u3940
  u732=-1.69411927340571579d0*u4116
  u664=2.82353212234285965d-1*u4116
  u843=-4.8000046079828614d0*u4116
  u4024=5.08235782021714737d0*u3936
  u3828=-8.47059636702857895d-1*u3936
  u681=2.82353212234285965d-1*u3936
  u525=-1.69411927340571579d0*u947
  u703=2.82353212234285965d-1*u947
  u859=2.44706117269714503d-1*u947
  value(3)=value(3)+(c(59)*(u2984*((u732+u2855*(u2855*(u525+u4008)+u402&
 &4))*u2857+(u664+u2855*(u2855*(u703+u760)+u3828))*u2856+u2855*(u843+u2&
 &855*(u2855*(u859+u760)+u681))+u972)))
  u3773=1.79688483954601253d0*u2858
  u734=u3773*(3.59225411974724908d5*u4026+u621)*u40
  u425=-5.98961613182004178d-1*u4116
  u3772=7.98615484242672237d-1*u3936
  u3802=-5.98961613182004178d-1*u8
  u757=5.98961613182004178d-1*u3936
  u585=u2858*(u488-5.67185232289765129d1*u628)*u40
  u628=-3.19446193697068895d0*u585
  u723=-9.98269355303340296d-1*u3936
  u379=-5.32410322828448158d-1*u947
  u641=1.99653871060668059d-1*u951
  u4087=-3.99307742121336118d-1*u947
  u642=u3876*(3.88453297767213889d8*u4026+u621*(6.1d1*u3898+u4009))*u40
  u704=5.19100064757736954d-1*u947
  u828=-3.99307742121336118d-2*u500
  u792=3.99307742121336118d-2*u4108
  u380=-7.98615484242672237d-2*u947
  u3904=u792*u2855
  u4092=u2855*(u828+u3904)
  u787=(u2983*(u2857*(u425+u2855*(u2855*(u704+u933)+u723)+(u2855*(u379+&
 &u3778)+u3772)*u2857)+u2856*(u3802+u2855*(u4092+u641)+(u2855*(u4087+u3&
 &904)+u757)*u2856)+u2855*(u628+u2855*(u380*u2855+u642))+u734))
  value(4)=value(4)+(c(59)*u787)
  u698=1.63284278170329504d5*u4026
  u4098=u2858*(u698+u3992)*u40
  u835=4.54579902879327859d0*u4098
  u479=-2.02035512390812382d0*u585
  u4130=1.01017756195406191d0*u3936
  u635=-1.89408292866386608d0*u4116
  u4093=9.4704146433193304d-2*u3936
  u893=6.58035641026427899d7*u4026
  u5=+u527*(u893+u621*(1.55d2+u4051))*u40
  u986=3.18404342432142532d7*u4026
  u860=1.34690341593874921d0*u2858*(u986+u3992*(1.0d1*u3898+u4016))*u40
  u526=-6.73451707969374606d-1*u947
  u4131=1.70467463579747947d0*u3936
  u381=-6.31360976221288694d-2*u947
  u861=u3838*(3.71259463275878192d9*u4026+u621*(5.83d2*u3898+u977))*u40
  u79=5.41287382134642304d8*u4026
  u3918=5.0d0*u3803
  u96=+u3879*(u79+u621*(u3918+u604))*u40
  u4004=6.73451707969374607d-2*u4108
  u4057=-3.03053268586218573d-1*u947
  u527=-2.39917170964089704d-1*u947
  u46=1.26272195244257739d-2*u4108
  u3862=u4004*u2855
  u4144=u46*u2855
  value(5)=value(5)+(c(59)*(u2985*(u2857*(u479+u2855*(u96*u2855+u860)+(&
 &u2855*(u526+u3862)+u4130)*u2857)+u2856*(u635+u2855*(u2855*(u4057+u414&
 &4)+u4131)+(u2855*(u381+u3830)+u4093)*u2856)+u2855*(u5+u2855*(u2855*(u&
 &527+u3830)+u861))+u835)))
  u89=-5.44510557438185616d-2*u501
  u812=1.95602029320319338d4*u3940
  u416=u3773*(u737+u621)*u40
  u4067=2.39584645272801671d0*u4116
  u949=u2858*(9.76439983458570431d6*u4026-2.26874092915906051d2*u4047)*&
 &u40
  u4047=-1.19792322636400836d0*u949
  u470=-u733
  u733=-5.98961613182004178d-1*u2858
  u826=+u733*(u3924-2.26874092915906051d2*u3897)*u40
  u3879=-3.59376967909202507d0*u3936
  u4012=u2858*(u986+u621*(5.0d0*u3898+u4011))*u40
  u986=2.39584645272801671d0*u4012
  u428=-5.98961613182004178d-1*u947
  u4115=4.45766079404999545d7*u4026
  u586=7.0d0*u3898
  u862=u3876*(u4115+u621*(u586+u4006))*u40
  u863=7.98615484242672237d-1*u947
  u499=u2858*(1.84037709925778383d9*u4026+u621*(1.7d1*u3803+u3766))*u40
  u60=-7.98615484242672237d-2*u499
  u4065=u3904+u60
  u595=u2855*(u4065)
  value(6)=value(6)+(c(59)*(u2857*(u812+u2855*(u2855*(u3879+u2855*(u933&
 &+u863))+u4067)+(u2855*(u825+u2855*(u3778+u811))+u622)*u2857)+u2856*(u&
 &41+u2855*(u2855*(u986+u595)+u4047)+(u2855*(u470+u2855*(u3904+u428))+u&
 &425)*u2856)+u2855*(u416+u2855*(u862*u2855+u826))+u89))
  u697=u2858*(5.11624071600365778d5*u4026+u621)*u40
  u82=2.54117891010857368d0*u697
  u3916=1.27361736972857013d6*u4026
  u4106=u2858*(u3916-5.67185232289765129d1*u4050)*u40
  u83=-6.77647709362286316d0*u4106
  u73=+u83
  u3784=8.47059636702857895d-1*u3936
  u76=-8.47059636702857895d-1*u4116
  u4110=1.41176606117142983d-1*u3936
  u907=4.88219991729285216d7*u4026
  u4030=u2858*(u907+u621*(1.15d2+u4015))*u40
  u802=-2.82353212234285965d-1*u4030
  u961=u2858*(1.9741069230792837d8*u4026+u621*(3.1d1*u3898+u4006))*u40
  u3928=5.6470642446857193d-1*u961
  u528=-5.6470642446857193d-1*u947
  u637=5.6470642446857193d-1*u3936
  u382=-9.4117737411428655d-2*u947
  u4200=u2858*(9.9978963523692755d8*u4026+u621*(1.57d2*u3898+u600))*u40
  u744=2.82353212234285965d-2*u4200
  u74=u2858*(1.40734719355006999d9*u4026+u621*(1.3d1*u3803+u592))*u40
  u674=-2.25882569787428772d-1*u74
  u978=-5.6470642446857193d-2*u947
  u804=9.41177374114286549d-3*u4108
  u3825=u804*u2855
  u4048=u2855*(u382+u3825)
  u3780=u978*u2855
  value(7)=value(7)+(c(59)*(u2985*(u2857*(u73+u2855*(u2855*(u674+u4008)&
 &+u3928)+(u2855*(u528+u644)+u3784)*u2857)+u2856*(u76+u2855*(u3780+u637&
 &)+(u4048+u4110)*u2856)+u2855*(u802+u2855*(u6*u999+u744))+u82)))
  u73=1.20830365846043833d6*u4026
  u802=-4.80213173743530218d-2*u2858*(u722+u73)*u40
  u3928=1.20753342621596855d5*u3940
  u637=-u3929
  u674=u2858*(5.04366992570573355d5*u4026+u621)*u40
  u489=4.75411042006094916d0*u674
  u410=-1.47905657513007307d1*u4116
  u3986=+u410
  u3933=-3.16940694670729944d0*u949
  u962=+u3933
  u793=-u81
  u81=2.75950430107856861d7*u4026
  u4141=u2858*(u81+u621*(6.5d1+u4017))*u40
  u684=-5.28234491117883239d-1*u4141
  u756=7.39528287565036535d0*u3936
  u4100=6.33881389341459888d0*u4012
  u383=-u4061
  u4086=u2858*(7.70538508685784927d8*u4026+u621*(1.21d2*u3898+u4009))*u&
 &40
  u864=3.5215632741192216d-2*u4086
  u529=-9.86037716753382047d-1*u947
  u4101=-2.11293796447153296d-1*u499
  u71=-1.40862530964768864d-1*u947
  value(8)=value(8)+(c(59)*((u3928+u2855*(u2855*(u756+u2855*(u449+u529)&
 &)+u3986))*u2857+u2856*(u3929+u2855*(u2855*(u4100+u2855*(u724+u4101))+&
 &u962)+(u2855*(u793+u2855*(u724+u383))+u637)*u2856)+u2855*(u489+u2855*&
 &(u2855*(u864+u71*u2855)+u684))+u802))
  u802=1.46375155574869183d5*u3940
  u3986=6.72332743306944622d0*u4116
  u962=-7.09684562379552657d0*u3936
  u684=3.17490462117168294d0*u3936
  u864=1.49407276290432138d0*u947
  u529=-3.73518190726080346d-1*u947
  value(9)=value(9)+(c(59)*(x*(u2856*(u3986+u2855*(u2855*(u864+u3835)+u&
 &962)+u651)+u2855*(u4025+u2855*(u2855*(u529+u4064)+u684))+u802)*z))
  u4101=2.24110914435648207d0*u4116
  u796=-7.47036381452160692d-1*u3936
  u71=2.49012127150720231d-1*u3936
  u651=4.98024254301440461d-1*u947
  u654=2.39051642064691421d0*u3936
  u4050=-5.47826679731584507d-1*u947
  u782=-1.49407276290432138d-1*u947
  u919=(u2855*(u651+u966)+u796)
  u830=u919*u2856
  u4145=u782*u2855
  value(1)=value(1)+(c(60)*(u2983*(u2856*(u4101+u2855*(u2855*(u4050+u40&
 &56)+u71)+u830)+u2855*(u4025+u2855*(u4145+u654))+u802)))
  u71=1.03502865104225876d5*u3940
  u4050=-1.05646898223576648d-1*u3936
  u607=2.11293796447153296d-1*u947
  u575=4.22587592894306592d0*u3936
  u941=-1.23254714594172756d0*u947
  u4170=+u941
  u778=-3.16940694670729944d-1*u947
  u3839=(u607+u601)
  u946=u2855*u3839
  value(2)=value(2)+(c(60)*(y*(u2856*(u814+u2855*(u2855*(u4170+u724)+u4&
 &194)+(u946+u4050)*u2856)+u2855*(u3858+u2855*(u778*u2855+u575))+u71)*z&
 &))
  u4170=-1.01647156404342947d1*u464
  u3815=u447*(9.21249897436999059d7*u4026+u621*(2.17d2+u4051))*u40
  u530=-9.4117737411428655d-2*u2858*(1.82764092556049813d9*u4026+u621*(&
 &2.87d2*u3898+u593))*u40
  u865=1.31764832376000117d0*u947
  u531=+u4018*(u658+u621*(9.1d1*u3898+u4009))*u40
  u4169=u4042*(2.3491872384643476d10*u4026+u621*(2.17d2*u3803+u3845))*u&
 &40
  u4172=u2855*(u4169+u3991)
  u3809=u638*u2855
  value(3)=value(3)+(c(60)*(u2983*(u2856*(u515+u2855*(u4172+u530)+(u285&
 &5*(u865+u3991)+u839)*u2856)+u2855*(u3815+u2855*(u3809+u531))+u4170)))
  u3993=u2858*(u909+u629)*u40
  u3966=-2.39584645272801671d-1*u3993
  u4180=u3876*(u4115+u621*(1.05d2+u3994))*u40
  u965=-2.79515419484935283d-1*u3936
  u438=7.42943465674999241d7*u4026
  u566=u3762*(u438+u621*(1.75d2+u3888))*u40
  u544=1.56018127791749841d9*u4026
  u841=2.45d2*u3898
  u384=+u3885*(u544+u621*(u841+u4035))*u40
  u705=5.59030838969870566d-1*u947
  u959=u2858*(u712+u621*(u726+u4006))*u40
  u385=-1.59723096848534447d-1*u959
  u4202=1.89450583747124806d10*u4026
  u885=1.75d2*u3803
  u565=u4036*(u4202+u621*(u885+u3844))*u40
  u936=-9.31718064949784276d-2*u4108
  u706=2.79515419484935283d-1*u947
  u3948=u936*u2855
  u551=u2855*(u565+u3948)
  u997=(y*(u2856*(u4180+u2855*(u551+u384)+(u2855*(u705+u3948)+u965)*u28&
 &56)+u2855*(u566+u2855*(u706*u2855+u385))+u3966)*z)
  value(4)=value(4)+(c(60)*u997)
  u631=5.44280927234431679d4*u4026
  u912=-6.61207131460840523d-1*u2858*(-u4037+u631)*u40
  u476=8.16421390851647518d5*u4026
  u649=u2858*(u476+u629)*u40
  u913=-1.51526634293109286d-1*u649
  u3978=2.02035512390812382d-1*u633
  u4215=-5.05088780977030955d-2*u3936
  u4199=8.33396488612101076d-1*u2858*(9.35173593157341702d5*u4026+u629)&
 &*u40
  u671=9.55213027296427596d7*u4026
  u931=+u4045*(u671+u621*(2.25d2+u3888))*u40
  u3888=1.70467463579747947d-1*u3936
  u567=u975*(1.74714177642252569d7*u4026-8.62121553080442996d3*u3621)*u&
 &40
  u3967=6.06106537172437146d-1*u3936
  u422=2.86563908188928279d8*u4026
  u70=4.5d1*u3898
  u953=u2858*(u422+u621*(u70+u3891))*u40
  u386=-5.05088780977030955d-2*u953
  u670=1.01017756195406191d-1*u947
  u708=8.70305202647856254d7*u4026
  u623=u2858*(u708+u621*(2.05d2+u4051))*u40
  u4174=-3.40934927159495895d-1*u623
  u866=u975*(8.59691724566784836d8*u4026+u621*(1.35d2*u3898+u3786))*u40
  u387=-3.40934927159495895d-1*u947
  u4175=+u608*(1.1356421546746417d9*u4026+u621*(2.675d3+2.72d2*pd2))*u4&
 &0
  u80=1.62386214640392691d9*u4026
  u434=1.5d1*u3803
  u4104=u2858*(u80+u621*(u434+u3766))*u40
  u3970=3.36725853984687303d-2*u4104
  u4203=-1.68362926992343652d-2*u4108
  u867=u3795*(8.31035333747892008d9*u4026+u621*(1.305d3*u3898+2.24d2*u4&
 &011))*u40
  u97=+u608*(4.38442779529060266d10*u4026+u621*(4.05d2*u3803+u4085))*u4&
 &0
  u868=u4032*(9.48844940447784745d8*u4026+u621*(1.49d2*u3898+u3786))*u4&
 &0
  u475=4.87158643921178074d9*u4026
  u3871=4.5d1*u3803
  u3843=u2858*(u475+u621*(u3871+u645))*u40
  u98=-5.05088780977030955d-2*u3843
  u3869=1.13644975719831965d-1*u4108
  u532=-5.68224878599159824d-2*u947
  u3810=u4203*u2855
  u3949=u3869*u2855
  u4146=u670*u2855
  value(5)=value(5)+(c(60)*(u2857*(u913+u2855*(u2855*(u386+u4146)+u3967&
 &)+u2857*((u4215+u2855*(u3810+u670))*u2857+u2855*(u386+u2855*(u3810+u3&
 &970))+u3978))+u2856*(u4199+u2855*(u2855*(u867+u2855*(u4028+u98))+u417&
 &4)+u2856*((u3888+u2855*(u4028+u387))*u2856+u2855*(u866+u2855*(u3949+u&
 &97))+u931))+u2855*(u567+u2855*(u2855*(u868+u532*u2855)+u4175))+u912))
  value(6)=value(6)+(c(60)*u3964)
  u3964=1.02673895357922169d-2*u4136
  u4179=1.69411927340571579d-1*u3925
  u998=1.48588693134999848d7*u4026
  u473=u2858*(u998-2.26874092915906051d2*u612)*u40
  u45=-8.47059636702857895d-2*u473
  u4219=1.97647248564000175d-1*u3936
  u569=u2858*(4.65662571078347103d5*u4026+u621)*u40
  u3969=-1.52470734606514421d0*u569
  u612=+u3969
  u4173=u2858*(u438+u621*(1.75d2+u3790))*u40
  u533=5.6470642446857193d-2*u4173
  u954=u2858*(u712+u621*(u726+u3891))*u40
  u898=1.69411927340571579d-1*u954
  u388=-3.95294497128000351d-1*u947
  u989=4.8d1*pd2
  u4198=u2858*(1.93165301075499803d8*u4026+u621*(4.55d2+u989))*u40
  u818=2.82353212234285965d-2*u4198
  u471=2.4517134367274975d9*u4026
  u773=3.85d2*u3898
  u4197=u2858*(u471+u621*(u773+u934))*u40
  u3971=-2.82353212234285965d-2*u4197
  u934=u2858*(u3864+u621*(u693+u727))*u40
  u573=-2.82353212234285965d-2*u934
  u838=3.12036255583499681d8*u4026
  u968=4.9d1*u3898
  u987=u2858*(u838+u621*(u968+u3833))*u40
  u4204=-3.7647094964571462d-2*u987
  u675=3.7647094964571462d-2*u3756
  u735=6.58824161880000585d-2*u947
  u609=-6.58824161880000585d-2*u4108
  u3872=u609*u2855
  value(7)=value(7)+(c(60)*(u2856*(u4179+u2855*(u2855*(u3971+u2855*(u38&
 &72+u675))+u533)+u2856*((u4219+u2855*(u3779+u388))*u2856+u2855*(u898+u&
 &573*u2855)+u45))+u2855*(u612+u2855*(u2855*(u4204+u735*u2855)+u818))+u&
 &3964))
  u3964=5.81057940229671564d0*u4116
  u735=-u4194
  u45=1.05646898223576648d0*u947
  u612=4.22587592894306592d-1*u3936
  u533=-3.5215632741192216d-2*u947
  u898=(u2855*(u45+u457)+u514)
  u818=(u467+u449)
  value(8)=value(8)+(c(60)*(x*(u2856*(u3964+u2855*(u2855*u818+u735)+u89&
 &8*u2856)+u2855*(u814+u2855*(u533*u2855+u612))+u915)*z))
  u3971=u2858*(u77-u3757)*u40
  u573=1.3582479662766558d-1*u3971
  u4204=-1.21979296312390986d5*u3940
  u4179=5.60277286089120519d-1*u4116
  u77=3.73518190726080346d-2*u3936
  u915=1.15790639125084907d1*u4116
  u827=-7.47036381452160691d-2*u947
  u576=-3.54842281189776328d0*u3936
  u869=9.33795476815200865d-1*u947
  u3820=1.24506063575360115d-1*u3936
  u870=9.96048508602880921d-2*u947
  u389=-1.24506063575360115d-2*u947
  u4147=u389*u2855
  value(9)=value(9)+(c(60)*(u2856*(u4204+u2855*(u2855*(u576+u2855*(u406&
 &4+u870))+u915)+u2856*((u77+u2855*(u4064+u827))*u2856+u2855*(u831+u285&
 &5*(u3835+u869))+u4179))+u2855*(u4204+u2855*(u2855*(u3820+u4147)+u4179&
 &))+u573))
  u3926=-7.47036381452160692d-1*u4116
  u881=-3.73518190726080346d0*u4116
  u680=+u881
  u3757=-u3763
  u3763=8.9644365774259283d-1*u3936
  u780=-7.47036381452160692d-1*u947
  u561=u3926+u2855*(u2855*(u780+u4056)+u3757)
  value(1)=value(1)+(c(61)*(x*(u2856*(u561+u830)+u2855*(u680+u2855*(u41&
 &43+u3763))+u4102)*z))
  u680=u2858*(u698-u599)*u40
  u3763=2.30502323396894504d-1*u680
  u830=u2858*(u631+u2966)*u40
  u631=-7.60657667209751865d0*u830
  u4036=+u631
  u3880=u2858*(u488+u621*(5.0d0+u926))*u40
  u418=-1.05646898223576648d-1*u3880
  u648=1.05646898223576648d-1*u3936
  u983=u2858*(u698-9.21676002470868334d1*u3621)*u40
  u688=-1.01421022294633582d1*u983
  u837=+u688
  u924=3.60858254756428203d7*u4026
  u3855=u2858*(u924+u621*(8.5d1+u4082))*u40
  u3979=6.33881389341459887d-1*u3855
  u4053=u3792*u4074*u40
  u689=-1.91748353630823987d2*u4053
  u840=+u689
  u905=2.8d1*pd2
  u511=u2858*(9.97666939620713267d7*u4026+u621*(2.35d2+u905))*u40
  u3798=1.05646898223576648d-1*u511
  u3792=u2858*(u80+u621*(2.55d2*u3898+5.2d1*u4011))*u40
  u450=-1.05646898223576648d-1*u3792
  u882=-1.5d1*u788
  u3856=-u621*(u645+u882)
  u534=u960*(u80+u3856)*u40
  u3775=1.5d1*u3898
  u960=u2858*(u671+u621*(u3775+u3832))*u40
  u3800=-2.11293796447153296d-1*u960
  u453=u2858*(u810+u3995*(4.0d0*u3803+u4074))*u40
  u691=3.38070074315445273d0*u453
  u61=-7.04312654823844319d-2*u4108
  u615=1.05646898223576648d-1*u947
  u3863=u615*u2855
  u3950=u61*u2855
  value(2)=value(2)+(c(61)*(u2856*(u4036+u2855*(u2855*(u450+u2855*(u457&
 &+u691))+u3979)+u2856*((u648+u2855*(u449+u407))*u2856+u2855*(u840+u285&
 &5*(u3950+u534))+u418))+u2855*(u837+u2855*(u2855*(u3800+u3863)+u3798))&
 &+u3763))
  u3763=-5.08235782021714737d0*u2858*(4.10024965183271864d5*u4026+u621)&
 &*u40
  u4036=1.69411927340571579d0*u8
  u418=-1.69411927340571579d0*u3936
  u837=+u418
  u3979=3.67059175904571754d0*u4116
  u840=-2.82353212234285965d-1*u3936
  u3798=u2858*(u81+u621*(6.5d1+u3994))*u40
  u450=8.47059636702857895d-1*u3798
  u534=-5.6470642446857193d-1*u951
  u3800=1.12941284893714386d0*u947
  u599=-3.10588533457714561d0*u3936
  u707=1.8823547482285731d-1*u947
  u535=-2.25882569787428772d-1*u2858*(u3807+u621*(u845+u3833))*u40
  u3807=1.12941284893714386d-1*u500
  u845=-1.12941284893714386d-1*u4108
  u871=5.08235782021714737d-1*u947
  u872=3.57647402163428889d-1*u947
  u677=u2855*(u871+u760)
  u3996=u845*u2855
  value(3)=value(3)+(c(61)*(u2985*(u2857*(u4036+u2855*(u2855*(u3807+u39&
 &96)+u534)+(u2855*(u3800+u3996)+u837)*u2857)+u2856*(u3979+u2855*(u677+&
 &u599)+(u2855*(u707+u760)+u840)*u2856)+u2855*(u450+u2855*(u872*u2855+u&
 &535))+u3763)))
  u4217=-8.71216891901096985d-2*u680
  u3944=-4.79169290545603342d-1*u649
  u512=6.38892387394137789d-1*u633
  u78=-1.59723096848534447d-1*u3936
  u398=u2858*(1.03908177017482411d5*u4026+u3995)*u40
  u3980=5.27086219600163676d0*u398
  u750=1.06134780810714177d7*u4026
  u3938=u2858*(u750-2.26874092915906051d2*u3791)*u40
  u3791=-1.19792322636400836d-1*u3938
  u4091=1.19792322636400836d-1*u3936
  u3761=3.83335432436482674d0*u983
  u3785=1.91667716218241337d0*u3936
  u613=-1.59723096848534447d-1*u953
  u616=3.19446193697068895d-1*u947
  u568=u2858*(u708+u621*(2.05d2+2.7d1*pd2))*u40
  u3841=-2.39584645272801671d-1*u568
  u708=4.79169290545603342d-1*u960
  u791=-2.39584645272801671d-1*u947
  u3833=-3.99307742121336118d-2*u511
  u4124=1.06482064565689632d-1*u4104
  u980=9.23372593053213343d8*u4026
  u943=2.8d1*u4011
  u429=1.45d2*u3898
  u4218=u2858*(u980+u621*(u429+u943))*u40
  u709=1.19792322636400836d-1*u4218
  u50=-3.99307742121336118d-2*u3843
  u710=7.98615484242672237d-2*u960
  u4183=u2858*(u79+u621*(u3918+u4074))*u40
  u3935=-3.19446193697068895d-1*u4183
  u4037=7.98615484242672237d-2*u4108
  u390=-3.99307742121336118d-2*u947
  u3951=u4037*u2855
  u4148=u616*u2855
  u3982=(u2857*(u3944+u2855*(u2855*(u613+u4148)+u3785)+u2857*((u78+u285&
 &5*(u932+u616))*u2857+u2855*(u613+u2855*(u932+u4124))+u512))+u2856*(u3&
 &980+u2855*(u2855*(u709+u2855*(u3904+u3935))+u3841)+u2856*((u4091+u285&
 &5*(u3904+u791))*u2856+u2855*(u708+u2855*(u3951+u50))+u3791))+u2855*(u&
 &3761+u2855*(u2855*(u710+u390*u2855)+u3833))+u4217)
  value(4)=value(4)+(c(61)*u3982)
  u444=9.09159805758655719d-1*u4098
  u423=2.52544390488515477d-1*u4116
  u424=5.05088780977030955d-2*u3936
  u3860=-7.57633171465546432d-2*u473
  u581=6.15581728702142228d7*u4026
  u3849=+u4062*(u581+u621*(1.45d2+u4051))*u40
  u4143=-1.26272195244257739d0*u3936
  u536=-1.01017756195406191d-1*u947
  u873=u3831*(1.05073433002607036d9*u4026+u621*(1.65d2*u3898+u4035))*u4&
 &0
  u874=u3838*(3.02484125310535405d9*u4026+u621*(4.75d2*u3898+u977))*u40
  u875=5.89270244473202781d-1*u947
  u62=-1.01017756195406191d-1*u4104
  u537=-2.14662731915238156d-1*u947
  value(5)=value(5)+(c(61)*(u2984*(u2857*(u423+u2855*(u2855*(u875+u3808&
 &)+u4143)+(u2855*(u536+u3946)+u424)*u2857)+u2856*(u3860+u2855*(u2855*(&
 &u62+u3997)+u873)+(u2855*(u387+u4028)+u3888)*u2856)+u2855*(u3849+u2855&
 &*(u2855*(u537+u3830)+u874))+u444)))
  value(6)=value(6)+(c(61)*u787)
  u787=u2858*(4.53567439362026399d5*u4026+u621)*u40
  u3878=5.08235782021714737d-1*u787
  u958=u2858*(u750-2.26874092915906051d2*u605)*u40
  u3977=-1.12941284893714386d-1*u958
  u3943=u2858*(1.29484432589071296d8*u4026+u621*(3.05d2+u3901))*u40
  u770=-5.6470642446857193d-2*u3943
  u753=2.54117891010857368d0*u3936
  u4107=7.32329987593927823d8*u4026
  u3767=1.15d2*u3898
  u988=u2858*(u4107+u621*(u3767+u600))*u40
  u408=5.6470642446857193d-2*u988
  u605=u2858*(u980+u621*(u429+u600))*u40
  u980=2.82353212234285965d-2*u605
  u429=-8.47059636702857895d-1*u947
  u3774=1.0d1*u3803
  u4132=u2858*(u79+u3992*(u3774+u604))*u40
  u3887=-3.01176759716571696d-1*u4132
  u774=-7.52941899291429239d-2*u947
  value(7)=value(7)+(c(61)*(u2984*((u76+u2855*(u2855*(u429+u644)+u753))&
 &*u2857+u2856*(u3977+u2855*(u2855*(u3887+u644)+u408)+(u2855*(u388+u377&
 &9)+u4219)*u2856)+u2855*(u770+u2855*(u2855*(u774+u765)+u980))+u3878)))
  u3878=u2858*(4.17282044213064287d5*u4026+u621)*u40
  u3977=4.75411042006094916d0*u3878
  u770=-3.69764143782518268d0*u4116
  u408=+u770
  u980=-1.58470347335364972d0*u8
  u4219=+u980
  u774=-u3852
  u3887=-2.11293796447153296d0*u958
  u3852=+u3887
  u76=3.69764143782518268d0*u3936
  u746=5.28234491117883239d-1*u951
  u446=-u45
  u775=8.5d1*u3898
  u3976=u2858*(u79+u621*(u775+u4009))*u40
  u4117=1.05646898223576648d-1*u3976
  u538=-7.39528287565036535d-1*u947
  u63=-1.05646898223576648d-1*u500
  u588=(u2855*(u446+u724)+u774)
  u3952=u4088*u2855
  value(8)=value(8)+(c(61)*(u2983*((u408+u2855*(u2855*(u538+u449)+u76))&
 &*u2857+u2856*(u4219+u2855*(u2855*(u63+u724)+u746)+u588*u2856)+u2855*(&
 &u3852+u2855*(u3952+u4117))+u3977)))
  u3852=-4.87917185249563944d4*u3940
  u4117=-u881
  u881=-2.61462733508256242d0*u3936
  u3799=+u881
  u591=-5.60277286089120519d-1*u3936
  u876=9.96048508602880922d-1*u947
  u4078=u2855*(u827+u4064)
  value(9)=value(9)+(c(61)*(y*(u2856*(u506+u2855*(u2855*(u876+u3835)+u3&
 &799)+(u4078+u77)*u2856)+u2855*(u4117+u2855*(u606+u591))+u3852)*z))
  u606=u2858*(4.68081597421611244d5*u4026+u621)*u40
  u799=2.24110914435648207d0*u606
  u4133=-u474
  u4224=-7.47036381452160692d-1*u8
  u789=7.47036381452160692d-1*u3936
  u3983=8.06624334161427748d6*u4026
  u844=u2858*(u3983-2.26874092915906051d2*u4063)*u40
  u800=-1.49407276290432138d0*u844
  u656=+u800
  u3946=-u747
  u747=2.49012127150720231d-1*u951
  u771=-4.98024254301440461d-1*u947
  u3770=2.10146866005214071d8*u4026
  u973=u2858*(u3770+u621*(3.3d1*u3898+u4006))*u40
  u42=1.49407276290432138d-1*u973
  u539=-u803
  u64=-4.98024254301440461d-2*u500
  u4140=-1.99209701720576184d-1*u947
  value(1)=value(1)+(c(62)*(u2983*((u4133+u2855*(u2855*(u539+u4056)+u39&
 &46))*u2857+u2856*(u4224+u2855*(u2855*(u64+u4056)+u747)+(u2855*(u771+u&
 &4056)+u789)*u2856)+u2855*(u656+u2855*(u4140*u2855+u42))+u799)))
  u656=u2858*(u476+u621)*u40
  u42=3.16940694670729944d-1*u656
  u4140=u2858*(u998-2.26874092915906051d2*u4060)*u40
  u476=-6.33881389341459887d-1*u4140
  u877=1.05646898223576648d-1*u954
  u925=4.13925645161785292d8*u4026
  u3870=6.5d1*u3898
  u4060=u2858*(u925+u621*(u3870+u4006))*u40
  u4058=1.05646898223576648d-1*u4060
  u632=u2858*(2.70643691067321152d9*u4026+u621*(2.5d1*u3803+u645))*u40
  u99=-3.5215632741192216d-2*u632
  u430=-4.22587592894306592d-1*u947
  u3989=u2855*(u99+u449)
  u3914=u2855*(u407+u449)
  u3953=u430*u2855
  value(2)=value(2)+(c(62)*(u2984*((u637+u2855*(u2855*(u383+u724)+u793)&
 &)*u2857+u2856*(u4050+u2855*(u3989+u877)+(u3914+u648)*u2856)+u2855*(u4&
 &76+u2855*(u3953+u4058))+u42)))
  u476=8.47059636702857895d-1*u697
  u4058=-8.18824315479429298d0*u4116
  u725=-u418
  u692=-2.82353212234285965d-1*u8
  u418=u2858*(3.82085210918571038d6*u4026-2.26874092915906051d2*u3837)*&
 &u40
  u412=-1.12941284893714386d0*u418
  u4123=4.8000046079828614d0*u3936
  u431=-u3800
  u711=9.4117737411428655d-2*u951
  u391=-1.8823547482285731d-1*u947
  u3984=1.59202171216071266d8*u4026
  u4221=2.5d1*u3898
  u956=u2858*(u3984+u621*(u4221+u4006))*u40
  u748=5.6470642446857193d-2*u956
  u540=-2.82353212234285965d-1*u947
  u51=-1.8823547482285731d-2*u500
  u982=1.8823547482285731d-2*u4108
  u541=-3.7647094964571462d-2*u947
  u3894=u982*u2855
  u415=u2855*(u51+u3894)
  value(3)=value(3)+(c(62)*(u2983*(u2857*(u4058+u2855*(u2855*(u540+u760&
 &)+u4123)+(u2855*(u431+u4008)+u725)*u2857)+u2856*(u692+u2855*(u415+u71&
 &1)+(u2855*(u391+u3894)+u681)*u2856)+u2855*(u412+u2855*(u541*u2855+u74&
 &8))+u476)))
  u666=3.59376967909202506d-1*u3925
  u4134=-9.98269355303340296d-1*u4116
  u4063=1.59723096848534447d-1*u3936
  u3837=-1.19792322636400836d-1*u3936
  u54=u2858*(u998-2.26874092915906051d2*u762)*u40
  u762=-2.39584645272801671d-1*u54
  u3854=1.39757709742467641d0*u3936
  u392=-3.19446193697068895d-1*u947
  u602=1.19792322636400836d-1*u954
  u660=u3876*(u712+u621*(u726+u4009))*u40
  u712=6.65512903535560198d-2*u947
  u726=-3.99307742121336118d-2*u632
  u766=u726+u3904
  u3848=u2855*(u766)
  u4201=u2855*(u3848+u602)
  u4205=u2855*(u791+u3904)
  u4149=u660*u2855
  u587=(u2984*(u2857*(u4134+u2855*(u2855*(u712+u933)+u3854)+(u2855*(u39&
 &2+u3778)+u4063)*u2857)+u2856*(u3837+u4201+(u4205+u4091)*u2856)+u2855*&
 &(u762+u4149)+u666))
  value(4)=value(4)+(c(62)*u587)
  u409=-9.07496371663624206d2*exppd2
  u647=-1.20532550005882387d-2*u2858*(3.03242230887754792d5*u4026+u409)&
 &*u40
  u928=5.68224878599159824d-2*u2858
  u589=-3.62998548665449682d3*u3621
  u4112=u928*(5.60609355051464629d6*u4026+u589)*u40
  u3804=-6.06106537172437146d-1*u633
  u4111=2.02035512390812382d-1*u3936
  u721=+u4045*(6.69465540498350965d6*u4026+u589)*u40
  u589=u2858*(u809-2.26874092915906051d2*u4068)*u40
  u3945=3.78816585732773216d-2*u589
  u4068=-1.89408292866386608d-2*u3936
  u3941=u928*(u698+u3895)*u40
  u3797=+u4062*(1.16748258891785595d8*u4026+u621*(2.75d2+u989))*u40
  u989=1.21221307434487429d0*u4012
  u542=-4.04071024781624764d-1*u947
  u66=2.1d1*pd2
  u4125=2.42442614868974858d0*u2858*(u488+u4000*(1.6d2+u66))*u40
  u3907=4.77606513648213798d8*u4026
  u543=+u4045*(u3907+u621*(7.5d1*u3898+u3891))*u40
  u617=3.78816585732773216d-2*u947
  u4142=u3840*(u750+u3995*(1.0d2+u4033))*u40
  u878=u3831*(u4107+u621*(u3767+u4035))*u40
  u65=+u579*(u79+u621*(u3918+u592))*u40
  u579=u2858*(u544+u621*(u841+u593))*u40
  u544=-1.89408292866386608d-2*u579
  u841=u3874*(u80+u621*(u434+u592))*u40
  u728=-6.31360976221288694d-3*u4108
  u545=+u608*(2.35619213399785474d8*u4026+u621*(3.7d1*u3898+u3891))*u40
  u546=-1.76781073341960834d-1*u947
  u4043=u4032*(4.00552662779635305d9*u4026+u621*(3.7d1*u3803+u645))*u40
  u3877=-1.89408292866386608d-2*u4108
  u713=6.31360976221288694d-3*u947
  u482=-1.26272195244257739d-2*u4108
  u761=u728*u2855
  u3882=u482*u2855
  u3781=u3877*u2855
  u4150=u65*u2855
  u4151=u713*u2855
  value(5)=value(5)+(c(62)*(u2857*(u4112+u2855*(u2855*(u878+u2855*(u383&
 &0+u546))+u3797)+u2857*((u4111+u2855*(u3862+u542))*u2857+u2855*(u989+u&
 &4150)+u3804))+u2856*(u721+u2855*(u2855*(u544+u2855*(u3882+u4043))+u41&
 &25)+u2856*((u4068+u2855*(u761+u617))*u2856+u2855*(u543+u2855*(u3781+u&
 &841))+u3945))+u2855*(u3941+u2855*(u2855*(u545+u4151)+u4142))+u647))
  u636=1.79688483954601253d0*u3925
  u929=u2858*(u998-2.26874092915906051d2*u4040)*u40
  u461=-3.99307742121336118d-1*u929
  u57=-3.99307742121336118d-1*u54
  u879=2.66205161414224079d-1*u984
  u547=-9.31718064949784276d-1*u947
  u4040=+u4023*(u3864+u621*(u693+u440))*u40
  u693=9.31718064949784276d-2*u4108
  u3864=u693*u2855
  value(6)=value(6)+(c(62)*(u2985*(u2857*(u461+u2855*(u4040*u2855+u879)&
 &+(u2855*(u547+u3864)+u3854)*u2857)+u2855*(u57+u4149)+u636)))
  u969=1.33893108099670193d6*u4026
  u477=-1.28342369197402711d-2*u2858*(-u409+u969)*u40
  u484=u2858*(2.15898101136324566d6*u4026+u3895)*u40
  u624=2.54117891010857368d-1*u484
  u646=-1.69411927340571579d-1*u3938
  u836=1.69411927340571579d-1*u3936
  u440=-8.47059636702857895d-2*u649
  u4190=1.12941284893714386d-1*u633
  u4149=-2.82353212234285965d-2*u3936
  u4054=u2858*(1.09944747301355199d6*u4026+u629)*u40
  u813=4.23529818351428947d-1*u4054
  u67=7.85397377999284912d7*u4026
  u4122=u2858*(u67+u621*(1.85d2+u3790))*u40
  u797=-2.25882569787428772d-1*u4122
  u3897=6.77647709362286316d-1*u960
  u393=-3.38823854681143158d-1*u947
  u468=3.38823854681143158d-1*u3936
  u394=-2.82353212234285965d-2*u953
  u626=5.6470642446857193d-2*u947
  u4059=-1.12941284893714386d0*u4106
  u4185=+u4059
  u4186=6.04968250621070811d8*u4026
  u4188=9.5d1*u3898
  u3790=u2858*(u4186+u621*(u4188+u4009))*u40
  u3824=1.12941284893714386d-1*u3790
  u3988=-5.6470642446857193d-2*u3843
  u4184=1.8823547482285731d-2*u4104
  u3859=4.70588687057143275d-2*u3936
  u469=-4.14118044610286082d-1*u947
  u749=9.41177374114286549d-3*u947
  u3865=u626*u2855
  value(7)=value(7)+(c(62)*(u2857*(u624+u2855*(u2855*(u3824+u2855*(u765&
 &+u469))+u797)+u2857*((u836+u2855*(u644+u393))*u2857+u2855*(u3897+u285&
 &5*(u4008+u3988))+u646))+u2856*(u440+u2855*(u2855*(u394+u3865)+u468)+u&
 &2856*((u4149+u2855*(u765+u626))*u2856+u2855*(u394+u2855*(u765+u4184))&
 &+u4190))+u2855*(u813+u2855*(u2855*(u3859+u749*u2855)+u4185))+u477))
  u477=4.75411042006094916d0*u606
  u624=u2858*(1.23116345740428446d7*u4026-2.26874092915906051d2*u598)*u&
 &40
  u646=-1.05646898223576648d0*u624
  u813=+u646
  u797=-3.16940694670729944d0*u844
  u3897=+u797
  u4185=1.91042605459285519d7*u4026
  u598=u2858*(u4185+u3992*(6.0d0*u3898+u4011))*u40
  u3988=8.45175185788613183d0*u598
  u3859=3.16940694670729944d-1*u973
  u469=u2858*(u3782+u621*(u72+u3766))*u40
  u3824=-2.11293796447153296d-1*u469
  u72=1.40862530964768864d-1*u4108
  u3782=u72*u2855
  value(8)=value(8)+(c(62)*(u2985*(u2857*(u813+u2855*(u2855*(u3824+u378&
 &2)+u3988)+u588*u2857)+u2855*(u3897+u2855*(u3953+u3859))+u477)))
  u813=-3.39561991569163951d-3*u2858*(-u409+9.96034096839009972d6*u4026&
 &)*u40
  u3897=4.26927537093368451d4*u3940
  u3859=u2858*(3.42896984157691957d6*u4026+u4039)*u40
  u3824=1.12055457217824104d-1*u3859
  u588=-u621*(pd2-5.0d0)
  u967=+u3900*(u588+u488)*u40
  u3900=-3.73518190726080346d-2*u3936
  u3769=7.84388200524768726d-1*u2858*(1.82722882714416349d6*u4026+u3895&
 &)*u40
  u3921=u2858*(u750+u3995*(1.0d2+u4084))*u40
  u4039=-1.79288731548518566d0*u3921
  u940=-u621*(u3891-4.5d1*u4014)
  u880=u630*(u422+u940)*u40
  u643=7.47036381452160691d-2*u947
  u678=2.9d1*pd2
  u55=+u990*(u998+u2966*(2.8d2+u678))*u40
  u4005=-u881
  u990=u2858*(3.5342882009967821d9*u4026+u621*(5.55d2*u3898+1.12d2*u401&
 &1))*u40
  u881=3.73518190726080346d-2*u990
  u4206=-u621*(u592+u882)
  u768=-4.98024254301440461d-2*u2858*(u4206+u80)*u40
  u3793=-1.24506063575360115d-2*u4108
  u882=u548*(1.34366632506364148d9*u4026+u621*(2.11d2*u3898+u600))*u40
  u548=-3.48616978011008323d-1*u947
  u48=u2858*(3.57249672208863921d9*u4026+u621*(3.3d1*u3803+u645))*u40
  u3857=-7.47036381452160691d-2*u48
  u4096=6.22530317876800576d-2*u4108
  u549=-8.71542445027520806d-2*u947
  u979=7.47036381452160691d-2*u4108
  u3789=u3793*u2855
  u4080=u2855*(u3789+u643)
  u3836=u979*u2855
  u3783=u4096*u2855
  value(9)=value(9)+(c(62)*((u3897+u2855*(u2855*(u4005+u2855*(u4064+u54&
 &8))+u4133))*u2857+u2856*(u3824+u2855*(u2855*(u881+u2855*(u3836+u3857)&
 &)+u4039)+u2856*((u3900+u4080)*u2856+u2855*(u880+u2855*(u3783+u768))+u&
 &967))+u2855*(u3769+u2855*(u2855*(u882+u549*u2855)+u55))+u813))
  u611=-7.31875777874345917d4*u3940
  u578=+u611
  u740=1.49407276290432138d0*u4116
  u798=-1.49407276290432138d-1*u3936
  u918=-u611
  u611=2.98814552580864276d-1*u947
  u4158=-u740
  u3801=1.49407276290432138d-1*u3936
  u4210=-2.98814552580864276d-1*u947
  u4223=u2855*(u966+u611)
  u3866=u4210*u2855
  u4152=u3801*u2855
  value(1)=value(1)+(c(63)*(u2856*(u578+u999*(u3866+u3757)+u2856*((u798&
 &+u4223)*u2856+u2855*(u831+u948*u999)+u740))+u2855*(u918+u2855*(u4152+&
 &u4158))))
  u578=-4.75411042006094916d0*u4116
  u918=+u578
  u4158=-2.6411724555894162d0*u4116
  u786=+u4158
  u395=-5.28234491117883239d-1*u947
  u4207=3.16940694670729944d-1*u3936
  u3985=-6.33881389341459887d-1*u947
  u4153=u3985*u2855
  value(2)=value(2)+(c(63)*(x*(u2856*(u918+u2855*((u395+u724)*u2856+u41&
 &53+u793)+(u601+u615)*u2867)+u2855*(u786+u4207*u2855)+u71)*z))
  u918=3.08021686073766507d-2*u4136
  u786=-5.08235782021714737d-1*u2858*(1.90498324532051088d6*u4026+u3895&
 &)*u40
  u991=2.25882569787428772d-1*u929
  u929=-3.95294497128000351d-1*u3936
  u4166=2.03294312808685895d0*u473
  u403=6.68649119107499317d8*u4026
  u582=1.05d2*u3898
  u550=+u663*(u403+u621*(u582+u3786))*u40
  u714=7.90588994256000702d-1*u947
  u884=1.12941284893714386d-1*u3896
  u4154=u714*u2855
  u4155=u929*u2855
  value(3)=value(3)+(c(63)*(u2856*(u786+u2855*(u2855*(u550+u4154)+u4166&
 &)+u2856*((u929+u2855*(u3991+u714))*u2856+u2855*(u550+u2855*(u3991+u88&
 &4))+u991))+u2855*(u786+u2855*(u4155+u991))+u918))
  u918=(x*(u2856*(u566+u2855*(u705*u2855+u384)+u2856*((u706+u3948)*u285&
 &6+u551+u385))+u2855*(u4180+u965*u2855)+u3966)*z)
  value(4)=value(4)+(c(63)*u918)
  u4166=-6.57934869456127549d3*u3621
  u884=u4077*(1.28994579754560308d7*u4026+u4166)*u40
  u4077=u2858*(u750+u621*u4076)*u40
  u991=1.51526634293109286d-1*u4077
  u786=-2.02035512390812382d-1*u4012
  u750=5.05088780977030955d-2*u947
  u551=-4.54579902879327859d0*u418
  u437=u3795*(1.77669623077135533d9*u4026+u621*(2.79d2*u3898+u4035))*u4&
 &0
  u396=-1.70467463579747947d-1*u947
  u917=3.37508602978071084d8*u4026
  u4157=+u4044*(u917+u621*(7.95d2+8.8d1*pd2))*u40
  u397=-5.05088780977030955d-2*u954
  u618=1.68362926992343652d-2*u632
  u3919=u4075*(2.34982404714921189d9*u4026+u621*(3.69d2*u3898+u458))*u4&
 &0
  u993=+u608*(3.99470088015366021d10*u4026+u621*(3.69d2*u3803+u4085))*u&
 &40
  u883=u3838*(5.27914399752492318d9*u4026+u621*(8.29d2*u3898+u977))*u40
  u895=(u750+u3810)*u2857
  u4208=(u396+u4028)*u2856
  u4156=u396*u2855
  value(5)=value(5)+(c(63)*(u2983*(u2857*(u991+u2855*(u4146+u397)+u2857&
 &*(u895+u2855*(u618+u3810)+u786))+u2856*(u551+u2855*(u2855*(u993+u4028&
 &)+u3919)+u2856*(u4208+u2855*(u993+u3949)+u437))+u2855*(u4157+u2855*(u&
 &4156+u883))+u884)))
  value(6)=value(6)+(c(63)*u997)
  u884=-2.82353212234285965d-1*u3939
  u991=u2858*(u838+u621*(u968+u3891))*u40
  u786=8.47059636702857895d-2*u991
  u551=-1.97647248564000175d-1*u947
  u437=2.82353212234285965d-1*u3939
  u883=u2858*(u3853+u621*(u685+u727))*u40
  u397=-2.82353212234285965d-2*u883
  u618=-8.47059636702857895d-2*u991
  u997=2.82353212234285965d-2*u883
  u883=1.97647248564000175d-1*u947
  u993=(u551+u3779)*u2856
  u4157=u883*u2855
  value(7)=value(7)+(c(63)*(u2983*(u2856*(u884+u999*(u3872+u997)+u2856*&
 &(u993+u397*u2855+u786))+u2855*(u437+u2855*(u4157+u618)))))
  u786=-u71
  u997=-u4158
  u3919=-3.16940694670729944d-1*u3936
  u397=-u578
  u618=6.33881389341459887d-1*u947
  u578=u2855*(u618+u457)
  value(8)=value(8)+(c(63)*(y*(u2856*(u997+u2855*(u2855*(u667+u449)+u77&
 &9)+(u578+u3919)*u2856)+u2855*(u397+u467*u999)+u786)*z))
  u997=-3.36166371653472311d-1*u3936
  u397=-3.73518190726080346d-2*u947
  u437=-4.8557364794390445d0*u3936
  u884=4.85573647943904449d-1*u947
  u4195=(u397+u4064)*u2856
  u4158=u397*u2855
  value(9)=value(9)+(c(63)*(u2983*(u2856*(u472+u2855*(u2855*(u884+u4064&
 &)+u437)+u2856*(u4195+u2855*(u884+u3835)+u997))+u2855*(u472+u2855*(u41&
 &58+u997))+u3861)))
  u884=-u4101
  u997=1.49407276290432138d0*u3936
  u437=-2.49012127150720231d-1*u947
  u448=(u611+u966)
  u3908=u2855*u448
  value(1)=value(1)+(c(64)*(y*(u2856*(u506+u2855*(u2855*(u437+u4056)+u7&
 &96)+(u3908+u798)*u2856)+u2855*(u884+u2855*(u4145+u997)))*z))
  u437=u2858*(9.79705669021977021d4*u4026+u3995)*u40
  u4178=-1.26776277868291978d1*u437
  u4222=+u4178
  u443=2.11293796447153296d0*u633
  u3961=7.0048955335071357d7*u4026
  u679=1.1d1*u3898
  u992=u2858*(u3961+u621*(u679+u4006))*u40
  u3968=1.05646898223576648d-1*u992
  u3817=u2858*(u3983+u3992*(3.8d1+u3927))*u40
  u3983=2.11293796447153296d0*u3817
  u436=-3.5215632741192216d-1*u951
  u572=u2858*(u810+u621*(u4006*(u926+1.0d0)+u552))*u40
  u552=-3.5215632741192216d-2*u572
  u784=4.64870339950928097d8*u4026
  u462=7.3d1*u3898
  u627=u2858*(u784+u621*(u462+u4009))*u40
  u553=-1.05646898223576648d-1*u627
  u3981=u2858*(8.98537054343506225d9*u4026+u621*(8.3d1*u3803+u603))*u40
  u914=3.5215632741192216d-2*u3981
  u653=3.16940694670729944d-1*u947
  u597=u818*u2856
  u3954=u653*u2855
  value(2)=value(2)+(c(64)*(u2983*(u2856*(u443+u2855*(u2855*(u914+u457)&
 &+u436)+u2856*(u597+u2855*(u552+u3950)+u3968))+u2855*(u3983+u2855*(u39&
 &54+u553))+u4222)))
  u4222=-1.01647156404342947d0*u3925
  u3968=1.69411927340571579d-1*u473
  u436=u594*(u438+u621*u3764)*u40
  u552=-5.6470642446857193d-2*u579
  u553=-2.25882569787428772d-1*u984
  u914=u4042*(u4202+u621*(u885+u3845))*u40
  u579=u2855*(u914+u3991)
  value(3)=value(3)+(c(64)*(y*(u2856*(u3968+u2855*(u579+u552)+(u2855*(u&
 &714+u3991)+u929)*u2856)+u2855*(u436+u2855*(u3809+u553))+u4222)*z))
  u438=7.66670864872965347d0*u983
  u4202=4.79169290545603342d-1*u4077
  u885=-6.38892387394137789d-1*u4012
  u661=1.59723096848534447d-1*u947
  u942=u2858*(u47-2.83592616144882564d1*u4034)*u40
  u570=-2.87501574327362005d1*u942
  u886=1.19792322636400836d-1*u961
  u3910=-1.19792322636400836d-1*u947
  u923=-1.59723096848534447d-1*u2858*(u3961+u3992*(3.3d2+4.1d1*pd2))*u4&
 &0
  u554=-1.59723096848534447d-1*u954
  u584=5.32410322828448158d-2*u632
  u887=3.99307742121336118d-1*u951
  u3812=4.4385565335040669d9*u4026
  u3813=4.1d1*u3803
  u4034=u2858*(u3812+u621*(u3813+u645))*u40
  u52=-3.99307742121336118d-2*u4034
  u399=5.66759729529213707d8*u4026
  u4146=8.9d1*u3898
  u888=u3876*(u399+u621*(u4146+u4009))*u40
  u401=(u661+u932)*u2857
  u596=(u3910+u3904)*u2856
  u3955=u3910*u2855
  u577=(u2983*(u2857*(u4202+u2855*(u4148+u554)+u2857*(u401+u2855*(u584+&
 &u932)+u885))+u2856*(u570+u2855*(u2855*(u52+u3904)+u887)+u2856*(u596+u&
 &2855*(u52+u3951)+u886))+u2855*(u923+u2855*(u3955+u888))+u438))
  value(4)=value(4)+(c(64)*u577)
  u439=4.84885229737949717d0*u983
  u935=-1.01017756195406191d0*u633
  u555=-5.05088780977030955d-2*u992
  u3816=+u4062*(u67+u621*(1.85d2+u4051))*u40
  u67=u2858*(1.17809606699892737d9*u4026+u621*(1.85d2*u3898+u4035))*u40
  u889=1.89408292866386608d-2*u67
  u939=+u4013*(3.79962515302356755d8*u4026+u621*(8.95d2+1.12d2*pd2))*u4&
 &0
  u890=1.68362926992343652d-1*u951
  u4214=1.68362926992343652d-2*u572
  u891=3.78816585732773216d-2*u67
  u67=+u608*(2.54405069603281883d10*u4026+u621*(2.35d2*u3803+u4085))*u4&
 &0
  u892=u3838*(4.04373514888821016d9*u4026+u621*(6.35d2*u3898+u977))*u40
  u4177=-1.68362926992343652d-2*u3981
  u571=3.36725853984687303d-2*u4108
  u556=-5.11402390739243842d-1*u947
  u557=-1.57840244055322173d-1*u947
  u4085=5.05088780977030955d-2*u4108
  u3811=u4085*u2855
  u4159=u571*u2855
  value(5)=value(5)+(c(64)*(u2985*(u2857*(u935+u2855*(u2855*(u4177+u381&
 &1)+u890)+u2857*(u895+u2855*(u4214+u4159)+u555))+u2856*(u3816+u2855*(u&
 &2855*(u556+u3830)+u891)+u2856*(u4208+u2855*(u67+u3997)+u889))+u2855*(&
 &u939+u2855*(u557*u2855+u892))+u439)))
  value(6)=value(6)+(c(64)*u3982)
  u3982=-2.03294312808685895d0*u3878
  u895=+u3982
  u4208=8.47059636702857895d-1*u8
  u405=u2858*(u907+u621*u39)*u40
  u402=-1.12941284893714386d-1*u405
  u3965=8.47059636702857895d-2*u953
  u4209=u2858*(u581-2.26874092915906051d2*u976)*u40
  u581=1.69411927340571579d-1*u4209
  u719=-2.82353212234285965d-1*u951
  u751=5.6470642446857193d-1*u947
  u718=u2858*(u4186+u621*(u4188+u600))*u40
  u4220=5.6470642446857193d-2*u718
  u550=5.95416120348106535d9*u4026
  u3948=5.5d1*u3803
  u916=u2858*(u550+u621*(u3948+u727))*u40
  u896=-2.82353212234285965d-2*u916
  u432=-2.82353212234285965d-2*u4200
  u4200=5.6470642446857193d-2*u500
  u9=-5.6470642446857193d-2*u4108
  u433=-1.41176606117142983d-1*u947
  u4076=1.78823701081714444d-1*u947
  u4022=u9*u2855
  value(7)=value(7)+(c(64)*(u2985*(u2857*(u4208+u2855*(u2855*(u4200+u40&
 &22)+u719)+(u2855*(u751+u4022)+u3828)*u2857)+u2856*(u402+u2855*(u2855*&
 &(u433+u765)+u4220)+u2856*(u993+u2855*(u896+u644)+u3965))+u2855*(u581+&
 &u2855*(u4076*u2855+u432))+u895)))
  u895=-7.68341077989648348d-2*u680
  u4208=3.38070074315445273d0*u983
  u402=2.3349651778357119d7*u4026
  u3965=u2858*(u402-2.26874092915906051d2*u452)*u40
  u581=-1.05646898223576648d-1*u3965
  u719=2.53552555736583955d0*u830
  u4220=u2858*(u893+u621*(1.55d2+u66))*u40
  u896=-2.11293796447153296d-1*u4220
  u4200=2.53552555736583955d0*u4012
  u4076=1.05646898223576648d-1*u4077
  u893=1.05646898223576648d-1*u959
  u66=-1.05646898223576648d-1*u632
  u993=u2858*(u3961+u621*(u679+u3832))*u40
  u679=-7.04312654823844319d-2*u993
  u574=u2858*(u810+u3992*(2.0d0*u3803+u4074))*u40
  u4075=5.63450123859075456d-1*u574
  u4071=7.04312654823844319d-2*u4108
  u894=3.5215632741192216d-2*u947
  u769=(u4207+u2855*(u724+u3985))
  u3956=u4071*u2855
  value(8)=value(8)+(c(64)*(u2856*(u4208+u2855*(u2855*(u893+u2855*(u601&
 &+u4075))+u896)+u2856*(u769*u2856+u2855*(u4200+u2855*(u3956+u66))+u581&
 &))+u2855*(u719+u2855*(u2855*(u679+u894*u2855)+u4076))+u895))
  u895=5.60277286089120519d-1*u947
  u896=1.86759095363040173d-1*u947
  value(9)=value(9)+(c(64)*(x*(u2856*(u4117+u2855*(u2855*(u896+u4064)+u&
 &3799)+u2856*(u4195+u2855*(u895+u3835)+u591))+u2855*(u506+u2855*(u4147&
 &+u77))+u3852)*z))
  u4076=5.43299186510662321d-2*u680
  u679=2.98814552580864276d-1*u4098
  u4075=u2858*(1.99569673319291615d5*u4026+u3992)*u40
  u4195=-2.68933097322777849d0*u4075
  u4035=+u4195
  u4148=2.98814552580864276d-1*u958
  u994=u2858*(u671+u621*(u3775+u4006))*u40
  u452=2.98814552580864276d-1*u994
  u3764=1.49407276290432138d-1*u3798
  u3775=3.50244776675356785d8*u4026
  u776=5.5d1*u3898
  u442=u2858*(u3775+u621*(u776+u4009))*u40
  u4147=-1.49407276290432138d-1*u442
  u754=u2858*(u80+u621*(u434+u645))*u40
  u3846=-4.98024254301440461d-2*u754
  u455=u2858*(u4115+u621*(u586+u4011))*u40
  u434=-1.99209701720576184d-1*u455
  u586=u2858*(7.57802334988499226d8*u4026+u621*(7.0d0*u3803+u592))*u40
  u3823=1.99209701720576184d-1*u586
  u752=4.98024254301440461d-2*u947
  value(1)=value(1)+(c(65)*(u2856*(u679+u2855*(u2855*(u4147+u2855*(u966&
 &+u3823))+u4148)+u2856*((u3801+u2855*(u4056+u4210))*u2856+u2855*(u452+&
 &u3846*u2855)+u798))+u2855*(u4035+u2855*(u2855*(u434+u752*u2855)+u3764&
 &))+u4076))
  u4076=u2858*(3.96547532699371651d5*u4026+u621)*u40
  u679=-4.43716972539021921d0*u4076
  u4035=+u679
  u3823=-u980
  u452=-3.16940694670729944d-1*u4077
  u3764=4.22587592894306592d-1*u4012
  u4147=u2858*(u671+u621*(2.25d2+u678))*u40
  u3846=2.11293796447153296d-1*u4147
  u434=-5.28234491117883239d-1*u951
  u4148=u2858*(u897+u621*(u441+u4009))*u40
  u39=-1.05646898223576648d-1*u4148
  u3821=1.05646898223576648d-1*u500
  u976=u2855*(u3821+u457)
  value(2)=value(2)+(c(65)*(u2985*(u2857*(u3823+u2855*(u976+u434)+u898*&
 &u2857)+u2856*(u452+u2855*(u3947+u877)+u2856*(u597+u3989+u3764))+u2855&
 &*(u3846+u2855*(u3954+u39))+u4035)))
  u4035=1.02673895357922169d-1*u680
  u452=-1.01647156404342947d0*u649
  u3764=1.35529541872457263d0*u633
  u3846=-3.38823854681143158d-1*u3936
  u39=2.48470826766171649d0*u398
  u3909=-5.6470642446857193d-2*u3938
  u99=5.6470642446857193d-2*u3936
  u597=-1.69411927340571579d0*u2858*(2.06826752349084038d5*u4026+u3992)&
 &*u40
  u441=4.06588625617371789d0*u3936
  u398=-3.38823854681143158d-1*u953
  u671=6.77647709362286316d-1*u947
  u877=-1.12941284893714386d-1*u568
  u897=2.25882569787428772d-1*u960
  u980=-1.12941284893714386d-1*u947
  u551=8.47059636702857895d-1*u4116
  u3989=2.25882569787428772d-1*u4104
  u898=5.6470642446857193d-2*u4218
  u4218=-1.8823547482285731d-2*u3843
  u465=u2858*(u809+u3995*(4.0d0*u3898+u4011))*u40
  u4196=3.01176759716571696d-1*u465
  u3893=-1.50588379858285848d-1*u4183
  u3850=3.7647094964571462d-2*u4108
  u558=-1.8823547482285731d-2*u947
  u3957=u3850*u2855
  u4160=u671*u2855
  value(3)=value(3)+(c(65)*(u2857*(u452+u2855*(u2855*(u398+u4160)+u441)&
 &+u2857*((u3846+u2855*(u3996+u671))*u2857+u2855*(u398+u2855*(u3996+u39&
 &89))+u3764))+u2856*(u39+u2855*(u2855*(u898+u2855*(u3894+u3893))+u877)&
 &+u2856*((u99+u2855*(u3894+u980))*u2856+u2855*(u897+u2855*(u3957+u4218&
 &))+u3909))+u2855*(u597+u2855*(u2855*(u4196+u558*u2855)+u551))+u4035))
  u406=-2.39584645272801671d-1*u764
  u984=u564*(3.09913559967285398d7*u4026+u621*(7.3d1+1.0d1*pd2))*u40
  u564=u2858*(u399+u621*(u4146+u3786))*u40
  u399=-3.99307742121336118d-2*u564
  u4146=-3.59376967909202506d-1*u4077
  u974=4.79169290545603342d-1*u4012
  u3819=+u3915*(u924+u621*(8.5d1+u4084))*u40
  u400=-6.65512903535560198d-2*u951
  u4084=u4027*(5.08810139206563766d9*u4026+u621*(4.7d1*u3803+u727))*u40
  u53=+u490*(5.73764625062720842d9*u4026+u621*(5.3d1*u3803+u603))*u40
  u38=-1.3310258070711204d-2*u4108
  u3958=u791*u2855
  u3987=u2855*(u3958+u602)
  u4161=u38*u2855
  u938=(u2985*(u2857*(u984+u2855*(u2855*(u53+u3904)+u400)+u2857*(u401+u&
 &2855*(u4084+u4161)+u399))+u2856*(u4146+u3987+u2856*(u596+u3848+u974))&
 &+u2855*(u3819+u2855*(u3955+u642))+u406))
  value(4)=value(4)+(c(65)*u938)
  u3890=-3.97757415019411877d-1*u2858*(5.3650548541679694d5*u4026+u621)&
 &*u40
  u596=+u4045*(2.18637648470071205d8*u4026+u621*u4069)*u40
  u4069=6.06106537172437146d-1*u4012
  u401=-2.02035512390812382d-1*u947
  u4176=9.4704146433193304d-2*u3798
  u937=+u4062*(5.73127816377856558d7*u4026+u621*(9.0d0*u3898+u4011))*u4&
 &0
  u619=1.89408292866386608d-2*u947
  u900=u3840*(u907+u3995*(4.6d2+u480))*u40
  u4070=u3795*(u559+u621*(u562+u458))*u40
  u559=-9.4704146433193304d-2*u951
  u562=5.5211312977733515d9*u4026
  u4044=5.1d1*u3803
  u4042=u3838*(u562+u621*(u4044+u645))*u40
  u480=-3.78816585732773216d-2*u956
  u4216=-1.32585805006470626d-1*u947
  u944=u2858*(9.74317287842356148d8*u4026+u621*(9.0d0*u3803+u592))*u40
  u592=5.05088780977030955d-2*u944
  u899=3.15680488110644347d-2*u947
  u920=(u401+u3862)*u2857
  u807=(u619+u761)*u2856
  value(5)=value(5)+(c(65)*(u2983*(u2857*(u596+u2855*(u2855*(u4216+u383&
 &0)+u4070)+u2857*(u920+u4150+u4069))+u2856*(u4176+u2855*(u2855*(u592+u&
 &3882)+u559)+u2856*(u807+u2855*(u4042+u3781)+u937))+u2855*(u900+u2855*&
 &(u899*u2855+u480))+u3890)))
  value(6)=value(6)+(c(65)*u587)
  u480=u2858*(5.98709019957874846d5*u4026+u621)*u40
  u937=8.47059636702857895d-2*u480
  u4216=u2858*(2.52600778329499742d8*u4026+u621*u4090)*u40
  u4070=-2.82353212234285965d-2*u4216
  u587=1.69411927340571579d-1*u956
  u456=-1.69411927340571579d-1*u947
  u900=u3917*(1.3160712820528558d7*u4026-2.26874092915906051d2*u3759)*u&
 &40
  u4176=8.27851290323570583d7*u4026
  u4211=1.3d1*u3898
  u3917=+u3886*(u4176+u621*(u4211+u3832))*u40
  u672=2.82353212234285965d-2*u947
  u921=-5.6470642446857193d-2*u589
  u4181=u3760*(1.87858562034964094d9*u4026+u621*(2.95d2*u3898+u593))*u4&
 &0
  u3754=-4.70588687057143275d-2*u951
  u3759=9.41177374114286549d-3*u500
  u4010=-4.8000046079828614d-1*u947
  u715=9.4117737411428655d-2*u947
  u731=(u456+u644)*u2857
  u4079=(u672+u765)*u2856
  u3867=u715*u2855
  value(7)=value(7)+(c(65)*(u2983*(u2857*(u4070+u2855*(u2855*(u4010+u76&
 &5)+u4181)+u2857*(u731+u2855*(u93+u4008)+u587))+u2856*(u900+u2855*(u38&
 &67+u3754)+u2856*(u4079+u2855*(u3759+u765)+u3917))+u2855*(u921+u749*u9&
 &99)+u937)))
  u937=9.50822084012189831d-1*u3925
  u4070=-1.90164416802437966d0*u3936
  u3754=+u4070
  u900=3.16940694670729944d-1*u954
  u3759=3.16940694670729944d-1*u956
  u921=u2855*(u66+u724)
  u4181=u2855*(u3985+u724)
  value(8)=value(8)+(c(65)*(u2984*((u814+u2855*(u2855*(u395+u449)+u774)&
 &)*u2857+u2856*(u3919+u2855*(u921+u900)+(u4181+u4207)*u2856)+u2855*(u3&
 &754+u2855*(u3953+u3759))+u937)))
  u3754=2.80138643044560259d0*u2858*(4.07122133571354896d5*u4026+u621)*&
 &u40
  u3759=-1.30731366754128121d0*u4116
  u4010=3.01422777502428264d7*u4026
  u3917=u2858*(u4010+u621*(7.1d1+u3994))*u40
  u430=-1.86759095363040173d-1*u3917
  u587=-u621*(u4011-u4014)
  u404=u2858*(u809+u587)*u40
  u3994=1.49407276290432138d-1*u404
  u662=3.73518190726080346d-2*u947
  u3832=-3.28696007838950704d1*u942
  u93=1.30731366754128121d0*u3936
  u901=3.11265158938400288d-1*u951
  u927=-u621*(u645-2.9d1*u788)
  u68=3.13946681638092536d9*u4026
  u3822=u2858*(u927+u68)*u40
  u3953=-1.24506063575360115d-2*u3822
  u995=u2858*(5.03078861042785201d8*u4026+u621*(7.9d1*u3898+u4009))*u40
  u902=7.47036381452160691d-2*u995
  u490=-2.61462733508256242d-1*u947
  u4066=u2858*(1.19083224069621307d9*u4026+u621*(1.1d1*u3803+u604))*u40
  u604=-1.99209701720576184d-1*u4066
  u560=-2.36561520793184219d-1*u947
  u3990=(u662+u3789)*u2856
  value(9)=value(9)+(c(65)*(u2983*((u3759+u2855*(u2855*(u490+u4064)+u93&
 &))*u2857+u2856*(u430+u2855*(u2855*(u604+u3836)+u901)+u2856*(u3990+u28&
 &55*(u3953+u3783)+u3994))+u2855*(u3832+u2855*(u560*u2855+u902))+u3754)&
 &))
  u720=4.48221828871296415d-1*u764
  u590=-3.58577463097037132d0*u633
  u904=+u590
  u903=1.49407276290432138d-1*u954
  u4090=1.49407276290432138d-1*u959
  u583=-4.98024254301440461d-2*u632
  u580=u2855*(u583+u4056)
  u4073=u2855*(u4210+u4056)
  value(1)=value(1)+(c(66)*(u2984*((u561)*u2857+u2856*(u798+u2855*(u580&
 &+u903)+(u4073+u3801)*u2856)+u2855*(u904+u2855*(u3866+u4090))+u720)))
  u3927=u2858*(u402+u621*(5.5d1+u4019))*u40
  u4090=-3.16940694670729944d-1*u3927
  u904=+u4055*(u4185+u621*(4.5d1+u926))*u40
  u926=1.05646898223576648d-1*u2858
  u561=u926*(u4186+u621*(u4188+u3891))*u40
  u673=+u4099*(u550+u621*(u3948+u645))*u40
  u435=-8.45175185788613183d-1*u947
  u3959=u435*u2855
  u3960=u612*u2855
  value(2)=value(2)+(c(66)*(u2983*(u2857*(u4090+u2855*(u3959+u561)+u285&
 &7*(u818*u2857+u2855*(u673+u3782)+u893))+u2855*(u904+u3960)+u42)))
  u4090=u594*(u909+u621)*u40
  u818=-3.10588533457714561d0*u4116
  u594=-5.6470642446857193d-2*u3936
  u893=-3.38823854681143158d-1*u4140
  u561=-6.77647709362286316d-1*u947
  u673=5.6470642446857193d-2*u954
  u904=5.6470642446857193d-2*u959
  u42=-1.8823547482285731d-2*u632
  u3902=u42+u3894
  u736=u2855*(u3902)
  u945=u2855*(u736+u673)
  u3892=u2855*(u980+u3894)
  value(3)=value(3)+(c(66)*(u2984*(u2857*(u818+u2855*(u2855*(u429+u760)&
 &+u483)+(u2855*(u561+u4008)+u468)*u2857)+u2856*(u594+u945+(u3892+u99)*&
 &u2856)+u2855*(u893+u904*u2855)+u4090)))
  u4150=-3.59376967909202506d-1*u473
  u716=1.19792322636400836d-1*u959
  u795=-2.79515419484935283d-1*u947
  u3911=-3.99307742121336118d-2*u54
  u54=-3.99307742121336118d-2*u3896
  u4182=(u795+u3864)*u2857
  u4162=u602*u2855
  u4163=u54*u2855
  u4164=u3911*u2855
  u4212=(u2983*(u2857*(u4150+u4162+u2857*(u4182+u4163+u716))+u4164+u666&
 &))
  value(4)=value(4)+(c(66)*u4212)
  u460=u928*(9.08949148481500903d6*u4026+u4166)*u40
  u4055=+u772*(8.36342072788427717d7*u4026+u621*(1.97d2+u905))*u40
  u905=u3831*(6.81385292804785018d8*u4026+u621*(1.07d2*u3898+u3786))*u4&
 &0
  u690=5.68224878599159824d-2*u4077
  u906=-7.57633171465546432d-2*u4012
  u4109=+u4013*(u907+u3992*(2.3d2+5.7d1*pd2))*u40
  u907=3.15680488110644347d-2*u2858*(1.03799815632878465d9*u4026+u621*(&
 &1.63d2*u3898+u593))*u40
  u4145=+u4013*(u562+u621*(u4044+u727))*u40
  u562=-1.89408292866386608d-2*u954
  u4044=6.31360976221288694d-3*u632
  u908=-6.06106537172437146d-1*u465
  u465=u3838*(u3812+u621*(u3813+u727))*u40
  u3796=u2855*(u4044+u761)
  u3812=u617*u2855
  u3813=u619*u2855
  value(5)=value(5)+(c(66)*(u2985*(u2857*(u4055+u2855*(u2855*(u465+u761&
 &)+u907)+u2857*(u920+u2855*(u4145+u3882)+u905))+u2856*(u690+u2855*(u38&
 &12+u562)+u2856*(u807+u3796+u906))+u2855*(u4109+u2855*(u3813+u908))+u4&
 &60)))
  u920=u2858*(-u620+u909)*u40
  u3903=-1.08902111487637123d-2*u920
  u807=3.59376967909202506d-1*u3993
  u4225=-1.19792322636400836d-1*u4140
  u458=2.79515419484935283d-1*u3936
  u466=7.18753935818405013d-1*u3925
  u620=-7.18753935818405013d-1*u473
  u909=2.39584645272801671d-1*u959
  u563=-5.59030838969870566d-1*u947
  value(6)=value(6)+(c(66)*(u2857*(u807+u2855*(u4162+u620)+u2857*((u458&
 &+u2855*(u3864+u563))*u2857+u2855*(u909+u4163)+u4225))+u2855*(u466+u41&
 &64)+u3903))
  u4164=u2858*(1.51854378698406438d7*u4026+u4166)*u40
  u4163=8.47059636702857895d-2*u4164
  u4162=u2858*(1.01464850455042754d8*u4026+u621*(2.39d2+u4051))*u40
  u4013=-1.41176606117142983d-1*u4162
  u4187=1.69411927340571579d-1*u961
  u3814=8.47059636702857895d-2*u4077
  u783=-1.12941284893714386d-1*u4012
  u3972=1.46465997518785565d8*u4026
  u4051=u2858*(u3972+u621*(3.45d2+u678))*u40
  u678=-5.6470642446857193d-2*u4051
  u3794=u2858*(u4107+u621*(u3767+u3786))*u40
  u772=1.41176606117142983d-1*u3794
  u910=-5.6470642446857193d-2*u4034
  u717=-2.82353212234285965d-2*u954
  u996=9.41177374114286549d-3*u632
  u981=1.07294220649028667d0*u3936
  u911=-9.31765600373143684d-1*u947
  u4165=u672*u2855
  value(7)=value(7)+(c(66)*(u2985*(u2857*(u4013+u2855*(u2855*(u911+u765&
 &)+u772)+u2857*(u731+u2855*(u910+u4008)+u4187))+u2856*(u3814+u2855*(u3&
 &865+u717)+u2856*(u4079+u2855*(u996+u765)+u783))+u2855*(u678+u2855*(u4&
 &165+u981))+u4163)))
  u4163=-2.88127904246118131d-2*u920
  u4013=u2858*(9.25277576298533854d5*u4026+u629)*u40
  u4187=9.50822084012189831d-1*u4013
  u4079=-1.05646898223576648d-1*u4141
  u783=1.90164416802437966d0*u764
  u678=-1.90164416802437966d0*u473
  u731=+u678
  u772=6.33881389341459887d-1*u956
  u717=-3.16940694670729944d-1*u3938
  u996=u2858*(u3775+u621*(u776+u3891))*u40
  u911=3.16940694670729944d-1*u996
  u910=-1.05646898223576648d-1*u3896
  value(8)=value(8)+(c(66)*(u2857*(u4187+u2855*(u2855*(u911+u3959)+u731&
 &)+u2857*(u769*u2857+u2855*(u772+u2855*(u3782+u910))+u4079))+u2855*(u7&
 &83+u2855*(u3960+u717))+u4163))
  u4163=u971*(1.25728894191153718d7*u4026+u4166)*u40
  u4187=-1.86759095363040173d-1*u2858
  u772=+u4187*(4.62747644334713813d7*u4026+u621*(1.09d2+u4015))*u40
  u776=1.12055457217824104d0*u3936
  u4079=1.12055457217824104d-1*u4077
  u717=-1.49407276290432138d-1*u4012
  u731=-1.49407276290432138d-1*u2858*(u810+u3992*(5.1d2+5.9d1*pd2))*u40
  u910=1.86759095363040173d-1*u564
  u564=-3.73518190726080346d-2*u954
  u769=1.24506063575360115d-2*u632
  u911=5.97629105161728553d-1*u2858*(u3961+u3992*(2.2d1*u3898+u4016))*u&
 &40
  u3775=+u801*(7.4697658734580638d9*u4026+u621*(6.9d1*u3803+u727))*u40
  u3765=8.71542445027520806d-2*u4108
  u4081=u2855*(u769+u3789)
  u3814=u3765*u2855
  u3961=u643*u2855
  u4166=u490*u2855
  value(9)=value(9)+(c(66)*(u2985*(u2857*(u772+u2855*(u2855*(u3775+u381&
 &4)+u910)+(u2855*(u780+u3836)+u776)*u2857)+u2856*(u4079+u2855*(u3961+u&
 &564)+u2856*(u3990+u4081+u717))+u2855*(u731+u2855*(u4166+u911))+u4163)&
 &))
  u3990=-u802
  u4191=-u4167
  u4167=-u654
  u654=1.49407276290432138d-1*u947
  u3913=-2.49012127150720231d-1*u3936
  u3868=5.47826679731584507d-1*u947
  u655=(u654+u966)
  u922=u655*u2856
  u3962=u771*u2855
  u3963=u789*u2855
  value(1)=value(1)+(c(67)*(u2983*(u2856*(u4191+u2855*(u3962+u3913)+u28&
 &56*(u922+u2855*(u3868+u4056)+u4167))+u2855*(u884+u3963)+u3990)))
  u3990=-u43
  u3913=-4.22587592894306592d-1*u3936
  u43=-u3964
  u3868=u774*u2855
  u4167=u446*u2855
  value(2)=value(2)+(c(67)*(y*(u2856*(u790+u2855*(u4167+u4194)+u2856*((&
 &u894+u601)*u2856+u2855*(u615+u724)+u3913))+u2855*(u43+u3868)+u3990)*z&
 &))
  u43=(u638+u3991)*u2856
  u894=u839*u2855
  value(3)=value(3)+(c(67)*(u2983*(u2856*(u3815+u2855*(u865*u2855+u530)&
 &+u2856*(u43+u4172+u531))+u2855*(u515+u894)+u4170)))
  u3815=u7*u2856
  u4169=u640*u2855
  u4170=u493*u2855
  u4194=(u2984*(u2857*(u657+u2856*(u2856*(u4128+u3815)+u376)+(u2856*(u6&
 &99+u3815)+u694)*u2857)+u2856*(u509+u2855*(u4169+u3974)+u2856*((u701+u&
 &933)*u2856+u4171+u700))+u2855*(u614+u4170)+u3768))
  value(4)=value(4)+(c(67)*u4194)
  u3964=u4203*u2856
  u4171=u387*u2855
  u4172=u3888*u2855
  value(5)=value(5)+(c(67)*(u2857*(u913+u2856*(u2856*(u386+u670*u2856)+&
 &u3967)+u2857*((u4215+u2856*(u3964+u670))*u2857+u2856*(u386+u2856*(u39&
 &64+u3970))+u3978))+u2856*(u567+u2855*(u2855*(u866+u4171)+u4174)+u2856&
 &*(u2856*(u868+u2855*(u3949+u98)+(u4028+u532)*u2856)+u2855*(u867+u2855&
 &*(u4028+u97))+u4175))+u2855*(u4199+u2855*(u4172+u931))+u912))
  value(6)=value(6)+(c(67)*u918)
  u965=-1.02673895357922169d-2*u920
  u913=-u3969
  u3966=-2.82353212234285965d-2*u4198
  u4174=3.7647094964571462d-2*u987
  u4175=-6.58824161880000585d-2*u947
  u912=-1.69411927340571579d-1*u3925
  u566=-5.6470642446857193d-2*u4173
  u3967=2.82353212234285965d-2*u4197
  u4197=8.47059636702857895d-2*u473
  u565=-1.69411927340571579d-1*u954
  u3869=2.82353212234285965d-2*u934
  u3978=-1.97647248564000175d-1*u3936
  u4173=u3978*u2855
  value(7)=value(7)+(c(67)*(u2856*(u913+u2855*(u2855*(u565+u3809)+u566)&
 &+u2856*(u2856*(u4174+u95*u2855+(u3779+u4175)*u2856)+u2855*(u3967+u285&
 &5*(u3872+u3869))+u3966))+u2855*(u912+u2855*(u4173+u4197))+u965))
  u4174=-u575
  u4175=-u941
  u3779=(u653+u457)
  value(8)=value(8)+(c(67)*(x*(u2856*(u781+u2855*(u3947+u735)+u2856*(u3&
 &779*u2856+u2855*(u4175+u449)+u4174))+u2855*(u790+u648*u2855)+u786)*z)&
 &)
  u4174=u827*u2855
  u4175=u77*u2855
  value(9)=value(9)+(c(67)*(u2856*(u4204+u2855*(u2855*(u831+u4174)+u915&
 &)+u2856*(u2856*(u3820+u2855*(u3835+u870)+(u4064+u389)*u2856)+u2855*(u&
 &576+u2855*(u4064+u869))+u4179))+u2855*(u4204+u2855*(u4175+u4179))+u57&
 &3))
  u3967=-u997
  u4197=2.49012127150720231d-1*u947
  value(1)=value(1)+(c(68)*(x*(u2856*(u4101+u2855*(u3866+u789)+u2856*(u&
 &922+u2855*(u4197+u4056)+u3967))+u2855*(u3926+u4152))*z))
  u3967=7.68341077989648348d-2*u680
  u4197=-u719
  u3869=-1.05646898223576648d-1*u4077
  u965=7.04312654823844319d-2*u993
  u913=-u4208
  u912=2.11293796447153296d-1*u4220
  u566=-1.05646898223576648d-1*u959
  u576=-5.63450123859075456d-1*u574
  u3820=1.05646898223576648d-1*u3965
  u934=-u4200
  u4152=1.05646898223576648d-1*u632
  u3965=u618*u2855
  u3966=u3919*u2855
  value(2)=value(2)+(c(68)*(u2856*(u4197+u2855*(u2855*(u934+u3965)+u912&
 &)+u2856*(u2856*(u965+u2855*(u3950+u576)+(u449+u533)*u2856)+u2855*(u56&
 &6+u2855*(u457+u4152))+u3869))+u2855*(u913+u2855*(u3966+u3820))+u3967)&
 &)
  value(3)=value(3)+(c(68)*(x*(u2856*(u436+u2855*(u4154+u552)+u2856*(u4&
 &3+u579+u553))+u2855*(u3968+u4155)+u4222)*z))
  u3967=u4091*u2855
  u4197=(u2857*(u3944+u2856*(u2856*(u613+u616*u2856)+u3785)+u2857*((u78&
 &+u2856*(u3815+u616))*u2857+u2856*(u613+u2856*(u3815+u4124))+u512))+u2&
 &856*(u3761+u2855*(u2855*(u708+u3958)+u3841)+u2856*(u2856*(u710+u2855*&
 &(u3951+u3935)+(u3904+u390)*u2856)+u2855*(u709+u2855*(u3904+u50))+u383&
 &3))+u2855*(u3980+u2855*(u3967+u3791))+u4217)
  value(4)=value(4)+(c(68)*u4197)
  u3869=u4085*u2856
  value(5)=value(5)+(c(68)*(u2984*(u2857*(u935+u2856*(u2856*(u4177+u386&
 &9)+u890)+u2857*((u750+u3964)*u2857+u2856*(u4214+u571*u2856)+u555))+u2&
 &856*(u939+u2855*(u2855*(u67+u4028)+u891)+u2856*((u557+u3830)*u2856+u2&
 &855*(u556+u3997)+u892))+u2855*(u3816+u2855*(u4156+u889))+u439)))
  value(6)=value(6)+(c(68)*u577)
  u929=-u3982
  u436=-8.47059636702857895d-1*u8
  u912=-1.69411927340571579d-1*u4209
  u576=2.82353212234285965d-1*u951
  u3820=-5.6470642446857193d-2*u500
  u3968=-1.78823701081714444d-1*u947
  u67=1.12941284893714386d-1*u405
  u4177=-5.6470642446857193d-2*u718
  u718=1.41176606117142983d-1*u947
  u4214=-8.47059636702857895d-2*u953
  u4198=2.82353212234285965d-2*u916
  u3816=u91*u2856
  value(7)=value(7)+(c(68)*(u2984*(u2857*(u436+u2856*(u2856*(u3820+u381&
 &6)+u576)+(u2856*(u528+u3816)+u3784)*u2857)+u2856*(u912+u2855*(u2855*(&
 &u4198+u3872)+u4177)+u2856*((u3968+u3825)*u2856+u2855*(u718+u4022)+u74&
 &4))+u2855*(u67+u2855*(u4157+u4214))+u929)))
  u929=-u4178
  u4214=-u3983
  u912=1.05646898223576648d-1*u627
  u576=-u443
  u3820=3.5215632741192216d-1*u951
  u627=-3.5215632741192216d-2*u3981
  u67=-1.05646898223576648d-1*u992
  u4177=3.5215632741192216d-2*u572
  u931=(u778+u724)
  u4198=u931*u2856
  value(8)=value(8)+(c(68)*(u2983*(u2856*(u4214+u2855*(u2855*(u4177+u60&
 &1)+u3820)+u2856*(u4198+u2855*(u627+u3956)+u912))+u2855*(u576+u2855*(u&
 &3863+u67))+u929)))
  value(9)=value(9)+(c(68)*(y*(u2856*(u506+u2855*(u2855*(u895+u4064)+u3&
 &799)+u2856*((u389+u4064)*u2856+u2855*(u896+u3835)+u77))+u2855*(u4117+&
 &u2855*(u4158+u591))+u3852)*z))
  u929=u2858*(u3924-2.26874092915906051d2*u3999)*u40
  u67=-1.49407276290432138d0*u929
  u912=+u67
  u576=1.49407276290432138d-1*u964
  u3820=-u67
  u929=u2858*(2.27340700496549768d9*u4026+u621*(2.1d1*u3803+u645))*u40
  u67=-4.98024254301440461d-2*u929
  u4177=-1.49407276290432138d-1*u964
  u627=4.98024254301440461d-2*u929
  u929=(u782+u4056)
  u4214=u929*u2856
  u3968=u654*u2855
  value(1)=value(1)+(c(69)*(u2983*(u2856*(u912+u999*(u966+u627)+u2856*(&
 &u4214+u67*u2855+u576))+u2855*(u3820+u2855*(u3968+u4177)))))
  u912=6.33881389341459887d-1*u764
  u576=u2858*(u402+u621*(5.5d1+u4089))*u40
  u3820=-2.11293796447153296d-1*u576
  u67=u2858*(u3972+u621*(2.3d1*u3898+u4006))*u40
  u4177=1.05646898223576648d-1*u67
  u627=9.50822084012189831d-1*u4077
  u436=-3.16940694670729944d-1*u954
  u572=-1.26776277868291978d0*u4012
  u570=+u572
  u4089=u2855*(u4152+u457)
  u3817=u952*u2856
  u405=u2856*(u59+u3817)
  value(2)=value(2)+(c(69)*(u2984*(u2857*(u682+u2856*(u405+u743)+(u2856&
 &*(u4088+u3817)+u69)*u2857)+u2856*(u3820+u2855*(u4089+u436)+u2856*(u46&
 &7*u2856+u578+u4177))+u2855*(u627+u2855*(u3954+u570))+u912)))
  u912=6.77647709362286316d-1*u764
  u3820=1.01647156404342947d0*u4077
  u4177=-1.35529541872457263d0*u4012
  u627=3.38823854681143158d-1*u947
  u570=-1.35529541872457263d1*u942
  u934=5.6470642446857193d-2*u961
  u965=-3.38823854681143158d-1*u958
  u575=-3.38823854681143158d-1*u954
  u913=1.12941284893714386d-1*u632
  u4222=1.8823547482285731d-1*u951
  u914=-1.8823547482285731d-2*u4034
  u443=(u627+u3996)*u2857
  u3983=(u978+u3894)*u2856
  value(3)=value(3)+(c(69)*(u2983*(u2857*(u3820+u2855*(u4160+u575)+u285&
 &7*(u443+u2855*(u913+u3996)+u4177))+u2856*(u570+u2855*(u2855*(u914+u38&
 &94)+u4222)+u2856*(u3983+u2855*(u914+u3957)+u934))+u2855*(u965+u2855*(&
 &u3780+u748))+u912)))
  u4160=u792*u2856
  u3820=(u2984*(u2857*(u984+u2856*(u2856*(u53+u4160)+u400)+u2857*((u661&
 &+u3815)*u2857+u2856*(u4084+u38*u2856)+u399))+u2856*(u3819+u4201+u2856&
 &*(u3910*u2856+u4205+u642))+u2855*(u4146+u2855*(u3955+u974))+u406))
  value(4)=value(4)+(c(69)*u3820)
  u4177=-6.42840266698039397d-2*u680
  u570=1.06068644005176501d0*u2858*(4.12098416334641128d5*u4026+u621)*u&
 &40
  u934=-5.05088780977030955d-2*u4140
  u965=u2858*(u4176+u621*(u4211+u4006))*u40
  u575=-1.68362926992343652d-2*u965
  u913=1.68362926992343652d-2*u947
  u4222=1.51526634293109286d-1*u649
  u914=-1.81831961151731144d0*u3936
  u912=1.51526634293109286d-1*u953
  u4199=-3.78816585732773216d-2*u589
  u567=-6.31360976221288694d-3*u965
  u915=1.01017756195406191d-1*u4098
  u438=-1.01017756195406191d-1*u4220
  u439=5.05088780977030955d-2*u959
  u569=2.69380683187749843d-1*u574
  u573=7.57633171465546432d-2*u568
  u568=+u4045*(1.30545780397178438d9*u4026+u621*(2.05d2*u3898+3.6d1*u40&
 &11))*u40
  u437=u4032*(2.4899219578193546d9*u4026+u621*(2.3d1*u3803+u3766))*u40
  u916=-1.26272195244257739d-2*u958
  u3815=-5.05088780977030955d-2*u632
  u3981=1.26272195244257739d-2*u3843
  u4176=u4057*u2855
  value(5)=value(5)+(c(69)*(u2857*(u570+u2855*(u2855*(u989+u4176)+u438)&
 &+u2856*(u2856*(u912+u4057*u2856)+u914)+u2857*(u2857*(u575+u2855*(u415&
 &9+u569)+u2856*(u3869+u4057)+(u3810+u913)*u2857)+u2856*(u912+u2856*(u3&
 &869+u62))+u2855*(u439+u2855*(u3811+u3815))+u934))+u2856*(u4222+u2855*&
 &(u2855*(u568+u2855*(u761+u437))+u573)+u2856*(u2856*(u567+u2855*(u3781&
 &+u437)+(u761+u713)*u2856)+u2855*(u568+u2855*(u3781+u3981))+u4199))+u2&
 &855*(u915+u2855*(u2855*(u567+u4151)+u916))+u4177))
  value(6)=value(6)+(c(69)*u938)
  u439=u2858*(5.38838117962087362d6*u4026-2.94936320790677867d3*u3621)*&
 &u40
  u570=5.6470642446857193d-2*u439
  u569=-2.03294312808685895d0*u3936
  u575=+u569
  u913=1.69411927340571579d-1*u953
  u3810=u2858*(u488+u3995*(2.0d1+pd2))*u40
  u914=-2.25882569787428772d-1*u3810
  u4199=-1.12941284893714386d-1*u4104
  u915=-9.41177374114286549d-3*u965
  u438=-5.6470642446857193d-2*u439
  u439=-u569
  u569=-1.69411927340571579d-1*u953
  u573=1.8823547482285731d-2*u499
  u934=2.25882569787428772d-1*u3810
  u3981=1.12941284893714386d-1*u4104
  u571=9.41177374114286549d-3*u965
  u4222=-1.8823547482285731d-2*u499
  u437=-9.41177374114286549d-3*u947
  u4177=u627*u2855
  value(7)=value(7)+(c(69)*(u2857*(u2855*(u2855*(u569+u4177)+u439)+u285&
 &6*(u2856*(u913+u393*u2856)+u575)+u2857*((u2855*(u4022+u627)+u2856*(u3&
 &816+u393))*u2857+u2856*(u913+u2856*(u3816+u4199))+u2855*(u569+u2855*(&
 &u4022+u3981))))+u2856*(u570+u999*(u2855*(u4222+u3825)+u718)+u2856*(u2&
 &856*(u915+u2855*(u765+u573)+(u765+u749)*u2856)+u2855*(u433+u804*u999)&
 &+u914))+u2855*(u438+u2855*(u2855*(u571+u437*u2855)+u934))))
  u3981=-6.33881389341459887d-1*u764
  u749=5.28234491117883239d-1*u8
  u573=-9.50822084012189831d-1*u4077
  u575=-u572
  u570=2.11293796447153296d-1*u576
  u438=-1.7607816370596108d-1*u951
  u914=-1.05646898223576648d-1*u67
  u569=3.5215632741192216d-2*u500
  u913=u2855*(u569+u601)
  value(8)=value(8)+(c(69)*(u2985*(u2857*(u749+u2855*(u913+u438)+u4094*&
 &u2857)+u2856*(u573+u2855*(u4153+u900)+u2856*(u4198+u921+u575))+u2855*&
 &(u570+u2855*(u3863+u914))+u3981)))
  u3981=-5.43299186510662321d-2*u680
  u573=2.98814552580864276d-1*u649
  u575=-7.47036381452160691d-2*u589
  u570=-1.24506063575360115d-2*u965
  u914=1.24506063575360115d-2*u947
  u4199=-1.49407276290432138d-1*u4220
  u915=u630*(9.87053461539641849d8*u4026+u621*(1.55d2*u3898+u943))*u40
  u4222=-2.4901212715072023d-2*u2858*(-u621*(u3766-u788)+u810)*u40
  u67=-7.47036381452160691d-2*u632
  value(9)=value(9)+(c(69)*(u2856*(u573+u2855*(u2855*(u915+u2855*(u3789&
 &+u4222))+u4199)+u2856*(u2856*(u570+u2855*(u3783+u4222)+(u3789+u914)*u&
 &2856)+u2855*(u915+u2855*(u3783+u67))+u575))+u2855*(u573+u2855*(u2855*&
 &(u570+u914*u2855)+u575))+u3981))
  u4199=-8.9644365774259283d-1*u649
  u570=7.47036381452160692d-1*u8
  u575=-4.48221828871296415d-1*u4077
  u4222=5.97629105161728553d-1*u4012
  u914=u2858*(u998-1.13437046457953026d2*u451)*u40
  u573=5.97629105161728553d-1*u914
  u439=-2.49012127150720231d-1*u951
  u571=-1.49407276290432138d-1*u956
  u451=4.98024254301440461d-2*u500
  u8=u2855*(u451+u966)
  value(1)=value(1)+(c(70)*(u2985*(u2857*(u570+u2855*(u8+u439)+u919*u28&
 &57)+u2856*(u575+u2855*(u3866+u903)+u2856*(u4214+u580+u4222))+u2855*(u&
 &573+u2855*(u3968+u571))+u4199)))
  u4199=1.15251161698447252d-1*u680
  u575=-6.33881389341459887d-1*u649
  u4222=8.45175185788613183d-1*u633
  u573=-2.11293796447153296d-1*u3936
  u941=3.16940694670729944d-1*u649
  u3866=-1.26776277868291978d0*u3936
  u934=+u3866
  u916=1.05646898223576648d-1*u953
  u4198=-4.22587592894306592d-1*u633
  u574=-7.04312654823844319d-2*u4104
  u943=-9.50822084012189831d-1*u649
  u4214=3.80328833604875932d0*u3936
  u572=-3.16940694670729944d-1*u953
  u4211=1.26776277868291978d0*u633
  u4203=2.11293796447153296d-1*u4104
  value(2)=value(2)+(c(70)*(u2857*(u575+u2855*(u2855*(u572+u3965)+u4214&
 &)+u2856*(u2856*(u916+u407*u2856)+u934)+u2857*((u573+u2855*(u457+u618)&
 &+u2856*(u3817+u407))*u2857+u2856*(u916+u2856*(u3817+u574))+u2855*(u57&
 &2+u2855*(u457+u4203))+u4222))+u2856*(u941+u2856*(u648*u2856+u4198))+u&
 &2855*(u943+u2855*(u3966+u4211))+u4199))
  u4199=u2858*(1.57841468897985187d6*u4026+u3895)*u40
  u572=-1.01647156404342947d0*u4199
  u4222=u447*(7.17471118280427839d7*u4026+u621*(1.69d2+u3920))*u40
  u573=+u3886*(9.61581114145070446d8*u4026+u621*(1.51d2*u3898+u600))*u4&
 &0
  u3920=-1.69411927340571579d-1*u4077
  u934=2.25882569787428772d-1*u4012
  u4198=u3873*(u924+u3992*(1.7d2+u4213))*u40
  u574=-4.70588687057143275d-1*u951
  u943=3.01176759716571696d-1*u4066
  u4203=-3.0d0*u4014
  u575=+u3886*(-u621*(u4006+u4203)+u4185)*u40
  u941=u2858*(u68+u927)*u40
  u3873=1.8823547482285731d-2*u941
  u68=-9.4117737411428655d-2*u4108
  u3969=u980*u2855
  u3886=u2855*(u3969+u673)
  u4178=u68*u2855
  value(3)=value(3)+(c(70)*(u2985*(u2857*(u4222+u2855*(u2855*(u3873+u38&
 &94)+u574)+u2857*(u443+u2855*(u943+u4178)+u573))+u2856*(u3920+u3886+u2&
 &856*(u3983+u736+u934))+u2855*(u4198+u2855*(u3780+u575))+u572)))
  u443=1.45202815316849498d-2*u680
  u3983=-2.39584645272801671d-1*u4199
  u4199=1.59723096848534447d-1*u576
  u576=+u4023*(u917+u621*(5.3d1*u3898+u3891))*u40
  u917=5.32410322828448158d-2*u947
  u3982=3.59376967909202506d-1*u649
  u4208=-1.43750787163681003d0*u3936
  u719=1.19792322636400836d-1*u953
  u4200=-4.79169290545603342d-1*u633
  u447=-7.98615484242672237d-2*u4104
  u4202=-3.99307742121336118d-2*u3993
  u4023=u2858*(4.03312167080713874d7*u4026+u621*(9.5d1+u4015))*u40
  u4209=1.59723096848534447d-1*u4023
  u577=+u3885*(1.68754301489035542d9*u4026+u621*(2.65d2*u3898+5.6d1*u40&
 &11))*u40
  u938=u4027*(6.38719110918877919d9*u4026+u621*(5.9d1*u3803+u727))*u40
  u936=-3.19446193697068895d-1*u585
  u918=u3762*(u3984+u621*(u4221+u3891))*u40
  u3984=+u3915*(u79+u621*(u3918+u3766))*u40
  u4221=(u2857*(u3983+u2855*(u2855*(u918+u3958)+u4209)+u2856*(u2856*(u7&
 &19+u791*u2856)+u4208)+u2857*(u2857*(u576+u2855*(u4161+u938)+u2856*(u4&
 &160+u791)+(u932+u917)*u2857)+u2856*(u719+u2856*(u4160+u447))+u2855*(u&
 &577+u2855*(u3904+u3984))+u4199))+u2856*(u3982+u2856*(u4091*u2856+u420&
 &0))+u2855*(u4202+u2855*(u3967+u936))+u443)
  value(4)=value(4)+(c(70)*u4221)
  u4201=1.87514209937722742d0*u2858*(3.91992223351666451d5*u4026+u621)*&
 &u40
  u999=-2.02035512390812382d0*u2858*(4.6699303556714238d6*u4026+u4000*(&
 &3.52d2+4.5d1*pd2))*u40
  u919=u3838*(4.08194366998006726d9*u4026+u621*(6.41d2*u3898+u977))*u40
  u935=+u4045*(u3770+u621*(4.95d2+5.6d1*pd2))*u40
  u920=2.20976341677451043d-1*u951
  u4000=+u3875*(u805+u621*(u785+u3844))*u40
  u785=u975*(u4185-u621*(u4011+u4203))*u40
  u4203=+u608*(-u621*(u645-4.9d1*u788)+u3853)*u40
  u975=4.41952683354902085d-2*u4108
  u585=u3838*(3.20527038048356815d8*u4026+u621*u459)*u40
  u4205=-9.78609513142997475d-1*u3936
  u578=-5.05088780977030955d-2*u947
  u579=-3.78816585732773216d-2*u954
  u921=7.57633171465546432d-2*u947
  u580=+u4045*(u925+u621*(u3870+u3891))*u40
  u922=4.10384634543837651d-1*u947
  u43=1.26272195244257739d-2*u632
  u4045=u2855*(u43+u3882)
  u3784=u728*u2856
  u4179=u619*u2856
  value(5)=value(5)+(c(70)*(u2984*(u2857*(u999+u2855*(u2855*(u922+u3830&
 &)+u4205)+u2856*(u2856*(u4203+u3784)+u920)+u2857*((u930+u3869)*u2857+u&
 &2856*(u4000+u975*u2856)+u2855*(u578+u3808)+u919))+u2856*(u935+u2855*(&
 &u4045+u579)+u2856*(u4179+u2855*(u921+u3882)+u785))+u2855*(u585+u2855*&
 &(u3813+u580))+u4201)))
  value(6)=value(6)+(c(70)*u4212)
  u4212=2.11764909175714474d0*u4116
  u584=u2858*(u4185+u3992*(9.0d1+u4082))*u40
  u923=-5.6470642446857193d-1*u584
  u4220=-1.41176606117142983d-1*u3917
  u997=2.35294343528571638d-1*u951
  u927=-1.50588379858285848d-1*u4066
  u998=1.12941284893714386d-1*u404
  u404=-9.41177374114286549d-3*u3822
  u3822=4.70588687057143275d-2*u4108
  u4215=-3.81176836516286053d0*u4116
  u4204=+u4215
  u4082=4.23529818351428947d-1*u3936
  u3972=-4.23529818351428947d-1*u947
  u3970=u6*u2856
  value(7)=value(7)+(c(70)*(u2984*(u2857*(u923+u2855*(u2855*(u3972+u765&
 &)+u445)+u2856*(u2856*(u404+u3970)+u997)+u2857*((u456+u3816)*u2857+u28&
 &56*(u927+u3822*u2856)+u2855*(u393+u644)+u744))+u2856*(u4220+u2856*(u6&
 &72*u2856+u998))+u2855*(u4204+u2855*(u4165+u4082))+u4212)))
  u3816=u2858*(u924+u621*u4105)*u40
  u4220=-3.16940694670729944d-1*u3816
  u997=u2858*(u925+u621*(u3870+u4009))*u40
  u4165=1.05646898223576648d-1*u997
  u393=7.9601085608035633d8*u4026
  u672=1.25d2*u3898
  u3822=u926*(u393+u621*(u672+u600))*u40
  u4204=+u4099*(9.20188549628891917d9*u4026+u621*(8.5d1*u3803+u603))*u4&
 &0
  u3972=u931*u2857
  value(8)=value(8)+(c(70)*(u2983*(u2857*(u4220+u2855*(u3959+u3822)+u28&
 &57*(u3972+u2855*(u4204+u3782)+u4165))+u2855*(u581+u3960)+u937)))
  u4220=u971*(u698+u621)*u40
  u4165=2.98814552580864277d0*u942
  u3822=-1.86759095363040173d-1*u3936
  u4204=3.73518190726080346d-2*u3798
  u581=-6.22530317876800576d-2*u951
  u923=1.24506063575360115d-1*u947
  u942=1.24506063575360115d-2*u500
  u937=+u3906*(u402+u621*u487)*u40
  u924=2.24110914435648207d-1*u954
  u402=-4.48221828871296415d-1*u947
  u925=1.12055457217824104d-1*u953
  u4105=u2855*(u67+u3836)
  u3870=u3793*u2856
  value(9)=value(9)+(c(70)*(u2984*(u2857*(u4165+u2855*(u2855*(u4072+u40&
 &64)+u3942)+u2856*(u2856*(u942+u3870)+u581)+(u2856*(u923+u3870)+u3822)&
 &*u2857)+u2856*(u4204+u2855*(u4105+u924)+u2856*(u662*u2856+u2855*(u402&
 &+u3836)+u717))+u2855*(u937+u2855*(u4166+u925))+u4220)))
  u3906=-4.48221828871296415d-1*u473
  u926=1.49407276290432138d-1*u956
  u487=-1.49407276290432138d-1*u3938
  u698=1.49407276290432138d-1*u996
  u998=-4.98024254301440461d-2*u3896
  u927=2.98814552580864276d-1*u3936
  u4166=-5.97629105161728553d-1*u947
  u404=9.96048508602880921d-2*u4108
  value(1)=value(1)+(c(71)*(u2983*(u2857*(u3906+u2855*(u4166*u2855+u698&
 &)+u2857*(u929*u2857+u2855*(u998+u404*u2855)+u926))+u2855*(u487+u927*u&
 &2855)+u720)))
  u3906=1.58470347335364972d0*u764
  u4166=u2858*(1.57079475599856982d7*u4026+u621*(3.7d1+u4019))*u40
  u404=-5.28234491117883239d-1*u4166
  u998=u2858*(1.84674518610642669d8*u4026+u621*(2.9d1*u3898+u4006))*u40
  u927=1.05646898223576648d-1*u998
  u487=-9.50822084012189831d0*u4116
  u3869=+u487
  u698=6.33881389341459888d0*u3936
  value(2)=value(2)+(c(71)*(u2984*(u2857*(u404+u2855*(u4167+u698)+u2856&
 &*(u4088*u2856+u743)+u2857*((u467+u3817)*u2857+u405+u4181+u927))+u2856&
 &*(u682+u69*u2856)+u2855*(u3869+u3868)+u3906)))
  u682=+u663*(u4115+u621*(1.05d2+u4019))*u40
  u927=u950*(u403+u621*(u582+u4006))*u40
  u612=-5.6470642446857193d-2*u4140
  u404=u950*(u403+u621*(u582+u3891))*u40
  u720=+u3826*(1.13670350248274884d10*u4026+u621*(1.05d2*u3803+u645))*u&
 &40
  u403=1.31764832376000117d-1*u4108
  value(3)=value(3)+(c(71)*(u2983*(u2857*(u682+u404*u2855+u2857*((u388+&
 &u403*u2855)*u2857+u720*u2855+u927))+u612*u2855+u4090)))
  u612=u3773*(4.60824518391818821d5*u4026+u621)*u40
  u3773=-1.99653871060668059d-1*u4030
  u720=u3876*(6.55912945410213616d8*u4026+u621*(1.03d2*u3898+u4009))*u4&
 &0
  u403=-1.73033354919245651d-1*u947
  u3826=-1.99653871060668059d0*u4116
  u3876=3.19446193697068895d0*u3936
  u404=-8.78477032666939461d-1*u947
  u3958=(u2984*(u2857*(u3773+u2855*(u4169+u3876)+u2856*(u4087*u2856+u64&
 &1)+u2857*((u403+u3778+u4160)*u2857+u2856*(u828+u4160)+u2855*(u404+u93&
 &3)+u720))+u2856*(u3802+u757*u2856)+u2855*(u3826+u4170)+u612))
  value(4)=value(4)+(c(71)*u3958)
  u3778=-5.73964523837535176d-4*u2858*(2.92278857924889811d7*u4026-1.45&
 &199419466179873d4*exppd2)*u40
  u3817=-7.71371915914080575d3*u3621
  u4169=u928*(1.41331614105207426d7*u4026+u3817)*u40
  u4170=+u608*(8.80918680728927672d8*u4026+u621*(2.075d3+2.22d2*pd2))*u&
 &40
  u928=2.10453658740429565d-3*u2858*(5.53386747147063721d9*u4026+u621*(&
 &8.69d2*u3898+u977))*u40
  u582=-6.73451707969374607d-2*u947
  u4115=-5.68224878599159824d-2*u649
  u405=2.2728995143966393d-1*u3936
  u583=-1.89408292866386608d-2*u953
  u931=7.57633171465546432d-2*u633
  u929=1.26272195244257739d-2*u4104
  u939=1.32585805006470626d-1*u2858*(1.25184613263919286d6*u4026+u3895)&
 &*u40
  u930=+u3875*(7.79029291150642061d8*u4026+u621*(1.835d3+2.88d2*pd2))*u&
 &40
  u4213=u3838*(1.2831695000015344d10*u4026+u621*(2.015d3*u3898+4.32d2*u&
 &4011))*u40
  u977=-8.41814634961718258d-3*u2858*(2.32753574317896191d10*u4026+u621&
 &*(2.15d2*u3803+u3845))*u40
  u932=1.89408292866386608d-2*u4077
  u743=-7.57633171465546432d-2*u994
  u4167=u3838*(u475+u621*(u3871+u727))*u40
  u3871=u4068*u2855
  u4180=u617*u2856
  u4181=u4068*u2856
  value(5)=value(5)+(c(71)*(u2857*(u4169+u2855*(u2855*(u743+u3812)+u930&
 &)+u2856*(u2856*(u583+u4180)+u405)+u2857*(u2857*(u928+u2855*(u3882+u97&
 &7)+u2856*(u3784+u617)+(u3862+u582)*u2857)+u2856*(u583+u2856*(u3784+u9&
 &29))+u2855*(u4213+u2855*(u761+u4167))+u4170))+u2856*(u4115+u2856*(u41&
 &81+u931))+u2855*(u939+u2855*(u3871+u932))+u3778))
  u743=u2858*(5.33395308689743045d5*u4026+u621)*u40
  u583=1.79688483954601253d0*u743
  u3862=+u733*(u3998+u621*u3827)*u40
  u929=u3762*(u838+u621*(u968+u4006))*u40
  u3762=-1.99653871060668059d-1*u3939
  u930=1.99653871060668059d-1*u991
  u991=+u3885*(u3853+u621*(u685+u645))*u40
  value(6)=value(6)+(c(71)*(u2985*(u2857*(u3862+u930*u2855+u2857*(u4182&
 &+u991*u2855+u929))+u3762*u2855+u583)))
  u977=3.75553839791757858d6*u4026
  u3827=-7.70054215184416268d-3*u2858*(-u4046+u977)*u40
  u733=u2858*(1.98273766349685826d6*u4026+u3895)*u40
  u4167=5.92941745692000526d-1*u733
  u3885=u2858*(2.44109995864642608d8*u4026+u621*(5.75d2+5.8d1*pd2))*u40
  u931=-2.82353212234285965d-2*u3885
  u939=4.26661818859070993d8*u4026
  u4213=6.7d1*u3898
  u435=u2858*(u939+u621*(u4213+u3891))*u40
  u3845=2.82353212234285965d-2*u435
  u4115=4.23529818351428947d-1*u3993
  u3959=u2858*(2.58968865178142593d7*u4026+u621*(6.1d1+u4017))*u40
  u744=-8.47059636702857895d-1*u3959
  u932=u3973*(3.24772429280785383d8*u4026+u621*(5.1d1*u3898+u3891))*u40
  u405=-4.51765139574857544d-1*u4183
  u4183=-1.27058945505428684d0*u4116
  u663=4.23529818351428948d0*u3936
  u59=-1.55294266728857281d0*u947
  u4182=u4149*u2855
  value(7)=value(7)+(c(71)*(u2857*(u4167+u2855*(u2855*(u663+u3865)+u744&
 &)+u2856*(u2856*(u394+u626*u2856)+u468)+u2857*(u2857*(u3845+u2855*(u40&
 &08+u405)+u2856*(u3970+u626)+(u644+u978)*u2857)+u2856*(u394+u2856*(u39&
 &70+u4184))+u2855*(u932+u2855*(u765+u59))+u931))+u2856*(u440+u2856*(u4&
 &149*u2856+u4190))+u2855*(u4115+u2855*(u4182+u4183))+u3827))
  u3827=u2858*(4.46310360332233976d5*u4026+u621)*u40
  u4167=4.75411042006094916d0*u3827
  u3970=u2858*(4.37275296940142411d7*u4026+u621*(1.03d2+u4015))*u40
  u3845=-5.28234491117883239d-1*u3970
  u4015=u2858*(5.28551208437356603d8*u4026+u621*(8.3d1*u3898+u4009))*u4&
 &0
  u663=1.05646898223576648d-1*u4015
  u744=-5.28234491117883239d-1*u3959
  u394=u2858*(9.10636419355927641d8*u4026+u621*(1.43d2*u3898+u600))*u40
  u931=1.7607816370596108d-1*u394
  u4115=u2858*(1.11505200719736315d10*u4026+u621*(1.03d2*u3803+u603))*u&
 &40
  u603=-3.5215632741192216d-2*u4115
  u59=-1.40862530964768864d0*u947
  u405=+u59
  value(8)=value(8)+(c(71)*(u2985*(u2857*(u3845+u2855*(u405*u2855+u931)&
 &+u2857*(u3972+u2855*(u603+u3782)+u663))+u2855*(u744+u821*u2855)+u4167&
 &)))
  u3845=-1.82968944468586479d4*u3940
  u663=5.60277286089120519d-1*u2858*(8.81735102119779319d5*u4026+u629)*&
 &u40
  u629=+u4187*(8.9153215880999909d6*u4026-2.26874092915906051d2*u3899)*&
 &u40
  u603=-1.12055457217824104d-1*u649
  u3899=4.48221828871296415d-1*u3936
  u405=-3.73518190726080346d-2*u953
  u744=1.49407276290432138d-1*u633
  u600=2.4901212715072023d-2*u4104
  u459=3.36166371653472311d-1*u2858*(1.79612705987362454d6*u4026+u3895)&
 &*u40
  u3960=-6.72332743306944622d-1*u3798
  u931=u971*(u4186+u621*(u4188+u3786))*u40
  u3786=-1.12055457217824104d-1*u2858*(u4185+u621*(4.5d1+u4019))*u40
  u932=4.48221828871296415d-1*u956
  u971=+u801*(7.03673596775034996d9*u4026+u621*(6.5d1*u3803+u727))*u40
  u801=2.61462733508256242d-1*u3936
  u440=-5.22925467016512484d-1*u947
  u4183=u440*u2855
  u4184=u801*u2855
  value(9)=value(9)+(c(71)*(u2857*(u663+u2855*(u2855*(u932+u4183)+u3960&
 &)+u2856*(u2856*(u405+u643*u2856)+u3899)+u2857*((u957+u2855*(u3836+u40&
 &2)+u2856*(u3870+u643))*u2857+u2856*(u405+u2856*(u3870+u600))+u2855*(u&
 &931+u2855*(u3814+u971))+u629))+u2856*(u603+u2856*(u3900*u2856+u744))+&
 &u2855*(u459+u2855*(u4184+u3786))+u3845))
  u805=2.71649593255331161d-1*u3971
  u3870=-u85
  u4190=-u494
  u4185=-u84
  u4186=-u3931
  u4187=-u420
  u4188=-u495
  u3971=u3757*u2855
  u3972=u780*u2855
  u3973=u3926*u2855
  value(1)=value(1)+(c(72)*(u2856*(u3870+u2855*(u3971+u4186)+u2856*(u28&
 &56*(u4190+u2855*(u4056+u4188)+(u966+u752)*u2856)+u2855*(u4187+u3972)+&
 &u4191))+u2855*(u4185+u3973)+u805))
  u4188=6.69097022082652103d-1*u947
  u420=u4*u2856
  u495=u724+u420
  u4185=u793*u2855
  u4186=u383*u2855
  u4187=u637*u2855
  value(2)=value(2)+(c(72)*(x*(u2856*(u790+u4185+u2856*(u2856*(u4188+u4&
 &95)+u4186+u735))+u4187+u3990)*z))
  u4188=u4118*u2855
  value(3)=value(3)+(c(72)*(u2856*(u56+u2855*(u4031*u2855+u4126)+u2856*&
 &(u2856*(u523+u4189+(u3991+u854)*u2856)+u2855*(u522+u853*u2855)+u515))&
 &+u2855*(u510+u4188)+u3818))
  u523=u4021*u2856
  u510=u933+u523
  u56=u4139*u2856
  u4189=u4103*u2855
  u4190=u625*u2855
  u4191=u41*u2855
  u4126=(u2985*((u622+u2856*(u2856*(u811+u56)+u825))*u2857+u2856*(u513+&
 &u4189+u2856*(u2856*(u639+u510)+u4190+u49))+u4191+u3905))
  value(4)=value(4)+(c(72)*u4126)
  u522=u3937*u2856
  u515=u3997+u522
  u695=u3806*u2856
  u3818=u413*u2856
  value(5)=value(5)+(c(72)*(u2983*(u2857*(u834+u2856*(u2856*(u857+u695)&
 &+u3975)+(u2856*(u377+u3818)+u4127)*u2857)+u2856*(u4193+u2855*(u524*u2&
 &855+u855)+u2856*(u2856*(u378+u515)+u2855*(u94+u4028)+u856))+u2855*(u3&
 &776+u729*u2855)+u833)))
  value(6)=value(6)+(c(72)*u4194)
  u4194=-u823
  u3974=-u824
  u4193=2.82353212234285965d-1*u4129
  u3975=-2.82353212234285965d-2*u498
  u933=-u4192
  u376=-u3923
  u4192=-9.4117737411428655d-2*u985
  u823=-9.88236242820000877d-1*u3936
  u3923=6.58824161880000585d-1*u947
  u4129=u804*u2856
  u498=u4022+u4129
  u824=u9*u2856
  value(7)=value(7)+(c(72)*(u2983*((u3974+u2856*(u2856*(u933+u824)+u403&
 &1))*u2857+u2856*(u4193+u2855*(u3923*u2855+u4192)+u2856*(u2856*(u980+u&
 &498)+u2855*(u675+u3872)+u3975))+u2855*(u376+u823*u2855)+u4194)))
  u3923=-u421
  u4193=-u686
  u686=-u821
  u4194=-u58
  u3974=-u739
  u3975=-u822
  u4192=u69*u2855
  value(8)=value(8)+(c(72)*(y*(u2856*(u4193+u2855*(u3952+u3974)+u2856*(&
 &(u615+u457)*u2856+u2855*(u3975+u449)+u686))+u2855*(u4194+u4192)+u3923&
 &)*z))
  u739=u3881*u2856
  u822=u3835+u739
  u3974=u3829*u2855
  u3975=u957*u2855
  value(9)=value(9)+(c(72)*(u2983*(u2856*(u472+u2855*(u3974+u3777)+u285&
 &6*(u2856*(u610+u822)+u2855*(u858+u4064)+u957))+u2855*(u472+u3975)+u38&
 &61)))
  u3777=-u4102
  u472=-8.9644365774259283d-1*u3936
  value(1)=value(1)+(c(73)*(y*(u2856*(u4117+u2855*(u3962+u831)+u2856*((&
 &u752+u966)*u2856+u2855*(u668+u4056)+u472))+u2855*(u506+u3963)+u3777)*&
 &z))
  u472=-u3977
  u3878=-u770
  u4193=-u3887
  u752=-u76
  u4194=-1.05646898223576648d-1*u3976
  u933=7.39528287565036535d-1*u947
  u4102=u2855*(u45*u2855+u434)
  u3861=u2855*(u3823+u514*u2855)
  value(2)=value(2)+(c(73)*(u2983*((u3878+u2856*(u2856*(u933+u420)+u752&
 &))*u2857+u2856*(u4193+u4102+u2856*((u650+u457)*u2856+u976+u4194))+u38&
 &61+u472)))
  u3887=u845*u2856
  u4193=u707*u2855
  u4194=u840*u2855
  value(3)=value(3)+(c(73)*(u2984*(u2857*(u4036+u2856*(u2856*(u3807+u38&
 &87)+u534)+(u2856*(u3800+u3887)+u837)*u2857)+u2856*(u450+u2855*(u4193+&
 &u599)+u2856*((u872+u760)*u2856+u677+u535))+u2855*(u3979+u4194)+u3763)&
 &))
  u3807=u2855*(u4087*u2855+u641)
  u4036=u2855*(u3802+u757*u2855)
  u3763=(u2983*(u2857*(u425+u2856*(u2856*(u704+u523)+u723)+(u2856*(u379&
 &+u56)+u3772)*u2857)+u2856*(u628+u3807+u2856*((u380+u3904)*u2856+u4092&
 &+u642))+u4036+u734))
  value(4)=value(4)+(c(73)*u3763)
  value(5)=value(5)+(c(73)*(u2985*(u2857*(u423+u2856*(u2856*(u875+u695)&
 &+u4143)+(u2856*(u536+u3818)+u424)*u2857)+u2856*(u3849+u2855*(u4171+u8&
 &73)+u2856*(u2856*(u537+u515)+u2855*(u62+u4028)+u874))+u2855*(u3860+u4&
 &172)+u444)))
  value(6)=value(6)+(c(73)*u4197)
  u3977=-5.08235782021714737d-1*u787
  u4197=5.6470642446857193d-2*u3943
  u3785=-u753
  u753=-2.82353212234285965d-2*u605
  u3976=8.47059636702857895d-1*u947
  u3979=7.52941899291429239d-2*u947
  u3980=1.12941284893714386d-1*u958
  u605=-5.6470642446857193d-2*u988
  u4217=3.01176759716571696d-1*u4132
  value(7)=value(7)+(c(73)*(u2985*((u551+u2856*(u2856*(u3976+u824)+u378&
 &5))*u2857+u2856*(u4197+u2855*(u3809+u605)+u2856*(u2856*(u3979+u498)+u&
 &2855*(u4217+u3872)+u753))+u2855*(u3980+u4173)+u3977)))
  u4173=-2.30502323396894504d-1*u680
  u4197=-u688
  u3785=-1.05646898223576648d-1*u511
  u753=2.11293796447153296d-1*u960
  u688=-u631
  u631=-6.33881389341459887d-1*u3855
  u3980=1.05646898223576648d-1*u3792
  u3792=-u691
  u4217=1.05646898223576648d-1*u3880
  u3978=-u689
  u3979=+u4099*(u3856+u80)*u40
  u3856=(u724+u467)
  u3976=u607*u2855
  u3977=u4050*u2855
  value(8)=value(8)+(c(73)*(u2856*(u4197+u2855*(u2855*(u3978+u3976)+u63&
 &1)+u2856*(u2856*(u753+u2855*(u3956+u3792)+u3856*u2856)+u2855*(u3980+u&
 &2855*(u601+u3979))+u3785))+u2855*(u688+u2855*(u3977+u4217))+u4173))
  value(9)=value(9)+(c(73)*(x*(u2856*(u4117+u2855*(u4174+u3799)+u2856*(&
 &u2856*(u3829+u822)+u2855*(u876+u4064)+u591))+u2855*(u506+u4175)+u3852&
 &)*z))
  u3799=-u4195
  u4197=-1.49407276290432138d-1*u3798
  u3785=1.99209701720576184d-1*u455
  u753=-2.98814552580864276d-1*u4098
  u4217=-2.98814552580864276d-1*u958
  u4195=1.49407276290432138d-1*u442
  u3980=-1.99209701720576184d-1*u586
  u442=-2.98814552580864276d-1*u994
  u586=4.98024254301440461d-2*u754
  u3978=u611*u2855
  u3979=u798*u2855
  value(1)=value(1)+(c(74)*(u2856*(u3799+u2855*(u2855*(u442+u3978)+u421&
 &7)+u2856*(u2856*(u3785+u3980*u2855+(u4056+u426)*u2856)+u2855*(u4195+u&
 &2855*(u966+u586))+u4197))+u2855*(u753+u2855*(u3979+u3801))+u3981))
  u426=-9.50822084012189831d-1*u3925
  u3980=-u4070
  u4197=-3.16940694670729944d-1*u956
  u753=4.22587592894306592d-1*u947
  u4217=u2855*(u4207+u3966)+u426
  u3966=u2855*(u3965+u436)
  value(2)=value(2)+(c(74)*(u2985*((u790+u2856*(u2856*(u667+u420)+u514)&
 &)*u2857+u2856*(u3980+u3966+u2856*((u753+u457)*u2856+u4089+u4197))+u42&
 &17)))
  u3980=u99*u2855
  value(3)=value(3)+(c(74)*(u2857*(u452+u2856*(u2856*(u398+u671*u2856)+&
 &u441)+u2857*((u3846+u2856*(u3887+u671))*u2857+u2856*(u398+u2856*(u388&
 &7+u3989))+u3764))+u2856*(u597+u2855*(u2855*(u897+u3969)+u877)+u2856*(&
 &u2856*(u4196+u2855*(u3957+u3893)+(u3894+u558)*u2856)+u2855*(u898+u285&
 &5*(u3894+u4218))+u551))+u2855*(u39+u2855*(u3980+u3909))+u4035))
  u4218=u2855*(u3837+u3967)+u666
  u3967=(u2985*(u2857*(u4134+u2856*(u2856*(u712+u523)+u3854)+(u2856*(u3&
 &92+u56)+u4063)*u2857)+u2856*(u762+u3987+u2856*(u2855*(u766+u4160)+u66&
 &0))+u4218))
  value(4)=value(4)+(c(74)*u3967)
  u4197=8.08142049563249528d0*u2858*(u47-1.77245385090551603d0*u3621*(1&
 &.28d2+1.5d1*pd2))*u40
  u3785=-1.60997048936428617d0*u3936
  u3981=2.52544390488515477d-1*u947
  u754=-3.78816585732773216d-2*u998
  u586=1.70467463579747947d-1*u947
  u441=1.89408292866386608d-2*u623
  u442=1.51526634293109286d-1*u954
  u3989=-7.57633171465546432d-2*u455
  u3621=u3811+u3818
  value(5)=value(5)+(c(74)*(u2983*(u2857*(u596+u2855*(u4176+u442)+u2856&
 &*(u2856*(u586+u522)+u3785)+u2857*((u401+u3621)*u2857+u2856*(u3981+u69&
 &5)+u2855*(u3815+u3811)+u4069))+u2856*(u4197+u2855*(u2855*(u4042+u761)&
 &+u559)+u2856*((u899+u3882)*u2856+u2855*(u592+u3781)+u754))+u2855*(u44&
 &1+u2855*(u3813+u3989))+u3890)))
  value(6)=value(6)+(c(74)*u3820)
  u4197=-8.47059636702857895d-2*u480
  u3785=2.82353212234285965d-2*u4216
  u3981=-1.69411927340571579d-1*u956
  u754=1.69411927340571579d-1*u947
  u586=+u4049*(u3916-1.13437046457953026d2*u4002)*u40
  u441=-2.4000023039914307d0*u3936
  u442=-2.82353212234285965d-2*u3798
  u3989=5.6470642446857193d-2*u632
  u3820=4.70588687057143275d-2*u951
  u4002=1.12941284893714386d-1*u4012
  u3819=-9.41177374114286549d-3*u500
  u500=-2.82353212234285965d-2*u947
  u4216=(u754+u4022)*u2857
  u3798=u2856*(u751+u824)
  u4195=u500*u2855
  u3916=u2855*(u4195+u4002)
  u4196=u437*u2856
  value(7)=value(7)+(c(74)*(u2983*(u2857*(u3785+u2855*(u4177+u565)+u285&
 &6*(u2856*(u718+u4129)+u441)+u2857*(u4216+u3798+u2855*(u3989+u4022)+u3&
 &981))+u2856*(u586+u2855*(u2855*(u3819+u3825)+u3820)+u2856*(u4196+u404&
 &8+u836))+u2855*(u442+u3916)+u4197)))
  u4197=-u679
  u3785=-2.11293796447153296d-1*u4147
  u3981=1.05646898223576648d-1*u4148
  u586=3.16940694670729944d-1*u4077
  u441=-1.05646898223576648d-1*u954
  u442=-4.22587592894306592d-1*u4012
  u4148=3.5215632741192216d-2*u632
  u4147=u2855*(u4148+u601)
  u3819=u3934*u2856
  u3989=u2856*(u63+u3819)
  value(8)=value(8)+(c(74)*(u2984*(u2857*(u4219+u2856*(u3989+u746)+(u28&
 &56*(u446+u3819)+u774)*u2857)+u2856*(u3785+u2855*(u4147+u441)+u2856*(u&
 &778*u2856+u946+u3981))+u2855*(u586+u2855*(u3863+u442))+u4197)))
  u4197=u662*u2855
  value(9)=value(9)+(c(74)*(u2983*((u3759+u2856*(u2856*(u490+u739)+u93)&
 &)*u2857+u2856*(u3832+u2855*(u2855*(u3953+u3789)+u901)+u2856*((u560+u3&
 &836)*u2856+u2855*(u604+u3783)+u902))+u2855*(u430+u2855*(u4197+u3994))&
 &+u3754)))
  u3785=8.9644365774259283d-1*u649
  u3981=-5.97629105161728553d-1*u914
  u586=4.48221828871296415d-1*u4077
  u442=-1.49407276290432138d-1*u954
  u914=-5.97629105161728553d-1*u4012
  u3754=4.98024254301440461d-2*u632
  u3832=u2855*(u3754+u966)
  u3820=u948*u2856
  u946=u2856*(u64+u3820)
  value(1)=value(1)+(c(75)*(u2984*(u2857*(u4224+u2856*(u946+u747)+(u285&
 &6*(u771+u3820)+u789)*u2857)+u2856*(u3981+u2855*(u3832+u442)+u2856*(u7&
 &82*u2856+u3908+u926))+u2855*(u586+u2855*(u3968+u914))+u3785)))
  u3785=3.16940694670729944d-1*u3816
  u3981=-1.05646898223576648d-1*u997
  u586=1.05646898223576648d0*u4116
  u3908=u3779*u2857
  value(2)=value(2)+(c(75)*(u2983*(u2857*(u3785+u3966+u2856*(u607*u2856&
 &+u686)+u2857*(u3908+u2856*(u650+u420)+u4089+u3981))+u2856*(u586+u4050&
 &*u2856)+u4217)))
  u3785=u982*u2856
  u3981=u978*u2856
  value(3)=value(3)+(c(75)*(u2984*(u2857*(u4222+u2856*(u2856*(u3873+u37&
 &85)+u574)+u2857*((u627+u3887)*u2857+u2856*(u943+u68*u2856)+u573))+u28&
 &56*(u4198+u945+u2856*(u3981+u3892+u575))+u2855*(u3920+u2855*(u3780+u9&
 &34))+u572)))
  u573=3.99307742121336118d-1*u4116
  u934=2.39584645272801671d-1*u947
  u943=u3904+u56
  u3873=(u2983*(u2857*(u4150+u3987+u934*u2867+u2857*((u795+u943)*u2857+&
 &u2856*(u4087+u523)+u3848+u716))+u2856*(u573+u3837*u2856)+u4218))
  value(4)=value(4)+(c(75)*u3873)
  u4198=u975*u2855
  value(5)=value(5)+(c(75)*(u2985*(u2857*(u999+u2855*(u2855*(u4203+u761&
 &)+u920)+u2856*(u2856*(u922+u522)+u4205)+u2857*((u377+u3621)*u2857+u28&
 &56*(u578+u695)+u2855*(u4000+u4198)+u919))+u2856*(u585+u2855*(u921*u28&
 &55+u579)+u2856*((u619+u3882)*u2856+u4045+u580))+u2855*(u935+u2855*(u3&
 &813+u785))+u4201)))
  value(6)=value(6)+(c(75)*u4221)
  u938=-u4212
  u443=5.6470642446857193d-1*u584
  u4202=-u4215
  u3987=-u445
  u445=-4.23529818351428947d-1*u3936
  u4205=4.23529818351428947d-1*u947
  u936=1.41176606117142983d-1*u3917
  u999=-2.35294343528571638d-1*u951
  u935=1.50588379858285848d-1*u4066
  u586=+u4018*(u587+u809)*u40
  u4203=9.41177374114286549d-3*u941
  u4018=-4.70588687057143275d-2*u4108
  u4199=u4018*u2855
  u4200=u500*u2856
  value(7)=value(7)+(c(75)*(u2985*(u2857*(u443+u2855*(u2855*(u4203+u382&
 &5)+u999)+u2856*(u2856*(u4205+u4129)+u3987)+u2857*(u4216+u2856*(u627+u&
 &824)+u2855*(u935+u4199)+u432))+u2856*(u4202+u2856*(u4200+u445))+u2855&
 &*(u936+u2855*(u4195+u586))+u938)))
  u938=-1.15251161698447252d-1*u680
  u443=6.33881389341459887d-1*u649
  u4202=-8.45175185788613183d-1*u633
  u587=2.11293796447153296d-1*u3936
  u936=9.50822084012189831d-1*u649
  u999=-u4214
  u935=3.16940694670729944d-1*u953
  u586=-u4211
  u4203=-2.11293796447153296d-1*u4104
  u4205=-3.16940694670729944d-1*u649
  u585=-u3866
  u584=-1.05646898223576648d-1*u953
  u4201=4.22587592894306592d-1*u633
  u3983=7.04312654823844319d-2*u4104
  value(8)=value(8)+(c(75)*(u2857*(u443+u2855*(u2855*(u584+u3976)+u585)&
 &+u2856*(u2856*(u935+u3985*u2856)+u999)+u2857*((u587+u2855*(u601+u607)&
 &+u2856*(u3819+u3985))*u2857+u2856*(u935+u2856*(u3819+u4203))+u2855*(u&
 &584+u2855*(u601+u3983))+u4202))+u2856*(u936+u2856*(u4207*u2856+u586))&
 &+u2855*(u4205+u2855*(u3977+u4201))+u938))
  u584=u2855*(u4197+u717)
  value(9)=value(9)+(c(75)*(u2985*(u2857*(u4165+u2855*(u2855*(u942+u378&
 &9)+u581)+u2856*(u2856*(u4072+u739)+u3942)+(u2855*(u923+u3789)+u3822)*&
 &u2857)+u2856*(u937+u2855*(u402*u2855+u924)+u2856*((u490+u3836)*u2856+&
 &u4105+u925))+u2855*(u4204+u584)+u4220)))
  u443=4.48221828871296415d-1*u649
  u585=-1.79288731548518566d0*u3936
  u587=+u585
  u936=1.49407276290432138d-1*u953
  u999=-5.97629105161728553d-1*u633
  u586=-9.96048508602880921d-2*u4104
  u4203=-4.48221828871296415d-1*u649
  u4205=-u585
  u585=-1.49407276290432138d-1*u953
  u942=5.97629105161728553d-1*u633
  u3983=9.96048508602880921d-2*u4104
  value(1)=value(1)+(c(76)*(u2857*(u2855*(u2855*(u585+u3978)+u4205)+u28&
 &56*(u2856*(u936+u4210*u2856)+u587)+u2857*((u4223+u2856*(u3820+u4210))&
 &*u2857+u2856*(u936+u2856*(u3820+u586))+u2855*(u585+u2855*(u966+u3983)&
 &)))+u2856*(u443+u2856*(u3801*u2856+u999))+u2855*(u4203+u2855*(u3979+u&
 &942))))
  u443=u2858*(4.02767886153479442d5*u4026+u621)*u40
  u587=-4.75411042006094916d0*u443
  u942=+u587
  u3983=5.28234491117883239d-1*u4023
  u999=-1.05646898223576648d-1*u995
  u585=3.16940694670729944d0*u4116
  value(2)=value(2)+(c(76)*(u2985*(u2857*(u3983+u4102+u2856*(u650*u2856&
 &+u686)+u2857*(u3908+u2856*(u607+u420)+u976+u999))+u2856*(u585+u3758*u&
 &2856)+u3861+u942)))
  u942=2.46417348859013206d-1*u680
  u3983=-1.08423633497965811d1*u983
  u999=1.12941284893714386d-1*u511
  u586=-2.25882569787428772d-1*u960
  u4203=1.69411927340571579d-1*u649
  u4205=-6.77647709362286316d-1*u3936
  u937=5.6470642446857193d-2*u953
  u4204=-2.25882569787428772d-1*u633
  u941=-3.7647094964571462d-2*u4104
  u4202=-1.69411927340571579d-1*u3859
  u4201=2.71059083744914526d0*u3921
  u4207=-5.6470642446857193d-2*u990
  u3821=1.12941284893714386d-1*u48
  u511=u2858*(u488+u588)*u40
  u938=1.12941284893714386d-1*u511
  u680=u2858*(u940+u422)*u40
  u43=-5.6470642446857193d-2*u680
  u983=u2858*(u80+u4206)*u40
  u3982=7.52941899291429239d-2*u983
  value(3)=value(3)+(c(76)*(u2857*(u3983+u2855*(u2855*(u43+u3969)+u4201&
 &)+u2856*(u2856*(u937+u980*u2856)+u4205)+u2857*(u2857*(u586+u2855*(u41&
 &78+u3821)+u2856*(u3785+u980)+(u3996+u702)*u2857)+u2856*(u937+u2856*(u&
 &3785+u941))+u2855*(u4207+u2855*(u3894+u3982))+u999))+u2856*(u4203+u28&
 &56*(u99*u2856+u4204))+u2855*(u4202+u2855*(u3980+u938))+u942))
  u3969=(u2985*(u2857*(u3773+u3807+u2856*(u640*u2856+u3876)+u2857*((u40&
 &3+u943)*u2857+u2856*(u404+u523)+u4092+u720))+u2856*(u3826+u493*u2856)&
 &+u4036+u612))
  value(4)=value(4)+(c(76)*u3969)
  u3983=+u4052*(9.61562971447495965d5*u4026+u621)*u40
  u999=u4032*(2.78073125724071145d8*u4026+u621*(6.55d2+u3901))*u40
  u586=+u3875*(u393+u621*(u672+u4009))*u40
  u4203=-1.3048126841906633d-1*u947
  u4205=5.68224878599159824d-1*u4116
  u4204=-1.89408292866386608d0*u3936
  u940=6.94497073843417563d-1*u947
  u4202=1.89408292866386608d-2*u3936
  u4206=-3.78816585732773216d-2*u947
  u4207=7.57633171465546432d-2*u589
  u3821=+u4029*(u79+u621*(u775+u3891))*u40
  u938=u3847*(u79+u3992*(u3774+u4074))*u40
  u43=-6.31360976221288694d-2*u4108
  u3982=u413*u2857
  u4201=u43*u2855
  u942=u4201+u695+u3982
  value(5)=value(5)+(c(76)*(u2983*(u2857*(u999+u2855*(u3812+u3821)+u285&
 &6*(u4206*u2856+u4204)+u2857*(u2857*(u4203+u942)+u2856*(u940+u522)+u28&
 &55*(u938+u761)+u586))+u2856*(u4205+u4202*u2856)+u2855*(u4207+u3871)+u&
 &3983)))
  value(6)=value(6)+(c(76)*u3958)
  u4203=5.6470642446857193d-1*u4116
  u999=7.5294189929142924d-1*u947
  u3774=2.82353212234285965d-2*u3936
  u3983=-5.6470642446857193d-1*u4116
  u938=-u497
  u4204=-7.5294189929142924d-1*u947
  u497=u644+u824
  u4202=u3774*u2856
  value(7)=value(7)+(c(76)*(u2983*(u2857*(u2855*(u3865+u938)+u2856*(u39&
 &81+u839)+u2857*((u497)*u2857+u2856*(u999+u4129)+u2855*(u4204+u765)))+&
 &u2856*(u4203+u4202)+u2855*(u3983+u4182))))
  u4204=-u587
  u443=-5.28234491117883239d-1*u4023
  u999=1.05646898223576648d-1*u995
  u3983=-u585
  value(8)=value(8)+(c(76)*(u2984*(u2857*(u443+u2855*(u3952+u821)+u2856&
 &*(u446*u2856+u746)+u2857*((u778+u3819)*u2857+u3989+u3914+u999))+u2856&
 &*(u4219+u774*u2856)+u2855*(u3983+u4192)+u4204)))
  u4204=6.72332743306944622d-1*u3925
  u443=-2.24110914435648207d-1*u3816
  u999=7.47036381452160691d-2*u997
  u3983=-3.73518190726080346d-1*u4116
  u938=-1.49407276290432138d-1*u958
  u4203=7.47036381452160691d-2*u988
  u586=-3.98419403441152369d-1*u4132
  u3821=(u610+u3836)*u2857
  value(9)=value(9)+(c(76)*(u2983*(u2857*(u443+u2855*(u4183+u4203)+u285&
 &6*(u827*u2856+u789)+u2857*(u3821+u2856*(u3829+u739)+u2855*(u586+u3814&
 &)+u999))+u2856*(u3983+u77*u2856)+u2855*(u938+u4184)+u4204)))
  u4183=u2858*(1.40097910670142714d7*u4026+u621*(3.3d1+u4019))*u40
  u4204=-7.47036381452160692d-1*u4183
  u999=u2858*(1.71938344913356967d8*u4026+u621*(2.7d1*u3898+u4006))*u40
  u3983=1.49407276290432138d-1*u999
  u4184=-4.48221828871296415d0*u4116
  u4203=+u4184
  u77=2.98814552580864277d0*u3936
  value(1)=value(1)+(c(77)*(u2984*(u2857*(u4204+u2855*(u3962+u77)+u2856&
 &*(u771*u2856+u747)+u2857*((u782+u3820)*u2857+u946+u4073+u3983))+u2856&
 &*(u4224+u789*u2856)+u2855*(u4203+u3963)+u4101)))
  u3820=3.16940694670729944d0*u3936
  u3983=u514*u2856
  u4203=u667*u2856
  u4204=u790*u2856
  value(2)=value(2)+(c(77)*(u2983*(u2857*(u3869+u4185+u3983+u2857*((u40&
 &7+u495)*u2857+u4203+u4186+u3820))+u4204+u4187+u71)))
  u383=8.47059636702857895d-1*u656
  u407=+u3755*(u4010+u621*(7.1d1+u4019))*u40
  u938=u950*(u939+u621*(u4213+u4006))*u40
  u950=-1.01647156404342947d1*u4116
  u3755=1.07294220649028667d1*u3936
  u586=-2.25882569787428772d0*u947
  value(3)=value(3)+(c(77)*(u2984*(u2857*(u407+u2855*(u4193+u3755)+u285&
 &6*(u391*u2856+u711)+u2857*((u456+u4008+u3785)*u2857+u2856*(u51+u3785)&
 &+u2855*(u586+u760)+u938))+u2856*(u692+u681*u2856)+u2855*(u950+u4194)+&
 &u383)))
  u4194=6.52006764401064461d4*u3940
  u4193=-1.15799245215187474d1*u4116
  u4008=7.58684710030538625d0*u3936
  u443=-1.25116425864685317d0*u947
  u3869=u510+u4139*u2857
  u4205=u4103*u2856
  u4206=u625*u2856
  u4207=u41*u2856
  u4010=(u2983*(u2857*(u4193+u4189+u4205+u2857*(u2857*(u443+u3869)+u420&
 &6+u4190+u4008))+u4207+u4191+u4194))
  value(4)=value(4)+(c(77)*u4010)
  u4189=-5.68224878599159824d-1*u480
  u480=1.89408292866386608d-1*u4077
  u939=u3831*(u3770-u621*(u4006-3.3d1*u4014))*u40
  u4190=4.7352073216596652d0*u4116
  u3770=-4.92461561452605181d0*u3936
  u940=1.02280478147848768d0*u947
  u3831=1.51526634293109286d0*u3810
  u587=+u676*(u4107+u621*(u3767+u3891))*u40
  u3767=u3874*(u550+u621*(u3948+u3766))*u40
  u3948=-9.4704146433193304d-2*u3936
  u941=6.31360976221288694d-2*u947
  value(5)=value(5)+(c(77)*(u2985*(u2857*(u480+u2855*(u941*u2855+u587)+&
 &u2856*(u381*u2856+u3770)+u2857*(u2857*(u378+u942)+u2856*(u940+u522)+u&
 &2855*(u3767+u761)+u939))+u2856*(u4190+u4093*u2856)+u2855*(u3831+u3948&
 &*u2855)+u4189)))
  u550=-1.81503519146061872d-2*u2858*(-u4046+1.60018592606922914d6*u402&
 &6)*u40
  u3874=7.18753935818405013d0*u2858*(1.77798436229914348d5*u4026+u3995)&
 &*u40
  u3995=+u3889*(u3998+u621*(4.9d1+u4020))*u40
  u942=5.32410322828448158d-2*u2858*(u838+u621*(u968+u4016))*u40
  u588=-9.31718064949784276d-2*u947
  u4016=5.98961613182004178d-1*u743
  u838=-7.98615484242672237d-1*u2858*(u3998+u3992*(9.8d1+u4033))*u40
  u3822=3.99307742121336118d-1*u987
  u987=-5.32410322828448158d-2*u2858*(u3853+u621*(u685+6.0d0*u4074))*u4&
 &0
  value(6)=value(6)+(c(77)*(u2857*(u3874+u838*u2855+u2857*(u2857*(u942+&
 &u987*u2855+(u3864+u588)*u2857)+u3822*u2855+u3995))+u4016*u2855+u550))
  u3853=3.31947651230994478d5*u3940
  u4074=-1.07294220649028667d1*u4116
  u3822=+u4074
  u685=2.25882569787428772d0*u3936
  u987=5.08235782021714737d0*u4116
  u4033=-5.36471103245143333d0*u3936
  u4210=+u4033
  u3998=-4.51765139574857544d0*u4116
  u838=+u3998
  u743=4.98824008280571871d0*u3936
  u4016=-1.0917657539725724d0*u947
  u968=+u4016
  u3864=-1.41176606117142983d-1*u3936
  u3984=u3864*u2855
  u4208=u382*u2856
  u4209=u4110*u2856
  value(7)=value(7)+(c(77)*(u2985*(u2857*(u3822+u2855*(u3867+u743)+u285&
 &6*(u4208+u4210)+u2857*((u980+u497)*u2857+u2856*(u3800+u4129)+u2855*(u&
 &968+u765)+u685))+u2856*(u987+u4209)+u2855*(u838+u3984)+u3853)))
  u3800=-5.17514325521129378d4*u3940
  u4210=+u3800
  u838=3.80328833604875932d1*u830
  u968=-3.16940694670729944d0*u418
  u3822=+u968
  u3889=3.38070074315445273d0*u759
  u943=-6.33881389341459888d0*u2858*(u481-1.13437046457953026d2*u3842)*&
 &u40
  u3823=3.16940694670729944d0*u993
  u4213=-4.22587592894306592d-1*u944
  u944=-2.11293796447153296d0*u4116
  u589=-2.11293796447153296d0*u947
  value(8)=value(8)+(c(77)*(u2857*(u838+u2855*(u698*u2855+u943)+u2857*(&
 &u2857*(u3889+u2855*(u3782+u4213)+u3856*u2857)+u2855*(u3823+u589*u2855&
 &)+u3822))+u2855*(u3906+u944*u2855)+u4210))
  u72=3.36166371653472311d0*u4116
  u3782=-1.12055457217824104d0*u4183
  u943=2.24110914435648207d-1*u999
  u3856=-1.12055457217824104d0*u4116
  u467=-2.98814552580864277d0*u3834
  u944=1.24506063575360115d-1*u985
  u985=-4.98024254301440461d-2*u3756
  u589=-8.71542445027520807d-1*u947
  value(9)=value(9)+(c(77)*(u2985*(u2857*(u3782+u2855*(u589*u2855+u944)&
 &+u2856*(u3829*u2856+u789)+u2857*(u3821+u2856*(u827+u739)+u2855*(u985+&
 &u3814)+u943))+u2856*(u3856+u957*u2856)+u2855*(u467+u93*u2855)+u72)))
  u3822=-u4097
  u4213=-u758
  u3823=-u970
  u4210=-u44
  u3821=u4114*u2856
  u44=u4056+u3821
  value(1)=value(1)+(c(78)*(u2983*(u2856*(u4213+u3946*u2855+u2856*(u285&
 &6*(u4210+u44)+u539*u2855+u3823))+u4133*u2855+u3822)))
  u3822=-u963
  u4213=-u492
  u3823=-u820
  u4210=u4083*u2855
  u4211=u427*u2855
  u4212=u3858*u2855
  value(2)=value(2)+(c(78)*(y*(u2856*(u4213+u4210+u2856*(u2856*(u3823+u&
 &495)+u4211+u755))+u4212+u3822)*z))
  u3822=u3922*u2856
  u820=u760+u3822
  u3823=u505*u2856
  value(3)=value(3)+(c(78)*(u2983*((u767+u2856*(u2856*(u516+u3823)+u388&
 &4))*u2857+u2856*(u683+u894+u2856*(u2856*(u847+u820)+u3809+u4031))+u41&
 &88+u829)))
  u4213=u463*u2855
  u4214=u659*u2855
  u4215=u794*u2855
  value(4)=value(4)+(c(78)*(u2984*((u4121+u2856*(u2856*(u520+u56)+u4041&
 &))*u2857+u2856*(u502+u4213+u2856*(u2856*(u851+u510)+u4214+u842))+u421&
 &5+u4135)))
  value(5)=value(5)+(c(78)*(u2857*(u88+u2856*(u2856*(u763+u2856*(u695+u&
 &850))+u417)+(u2856*(u485+u2856*(u3818+u517))+u4007)*u2857)+u2856*(u41&
 &38+u2855*(u4001*u2855+u955)+u2856*(u2856*(u849+u2855*(u4028+u92)+u285&
 &6*(u522+u3997+u519))+u2855*(u848+u518*u2855)+u90))+u2855*(u3771+u4119&
 &*u2855)+u3932))
  value(6)=value(6)+(c(78)*u4126)
  u3989=7.70054215184416268d-2*u501
  u3985=-u86
  u4221=-u3930
  u3930=-u75
  u4224=-4.23529818351428947d-1*u2858*(u4038+u3924)*u40
  u4216=-u3884
  u4217=u3760*(u658+u634)*u40
  u4222=-u738
  u4223=-3.10588533457714561d-1*u947
  u945=-u87
  u4218=8.47059636702857895d-1*u3939
  u4219=-8.47059636702857895d-1*u964
  u4220=5.6470642446857193d-2*u3896
  u946=9.88236242820000877d-1*u4116
  u69=9.88236242820000877d-1*u947
  value(7)=value(7)+(c(78)*((u3985+u2856*(u2856*(u4216+u2856*(u824+u422&
 &2))+u3930))*u2857+u2856*(u4221+u2855*(u3987*u2855+u4218)+u2856*(u2856&
 &*(u4217+u2855*(u3872+u4220)+u2856*(u4129+u4022+u4223))+u2855*(u4219+u&
 &69*u2855)+u4224))+u2855*(u945+u946*u2855)+u3989))
  u4218=-u3883
  u4219=-u478
  u4220=-u496
  u478=u508*u2856
  u496=u449+u478
  u4216=u395*u2855
  u4217=u814*u2855
  value(8)=value(8)+(c(78)*(x*(u2856*(u4218+u3868+u2856*(u2856*(u4220+u&
 &496)+u4216+u4219))+u4217+u3923)*z))
  value(9)=value(9)+(c(78)*(u2856*(u419+u2855*(u3942*u2855+u832)+u2856*&
 &(u2856*(u4120+u2855*(u4064+u852)+u2856*(u739+u3835+u521))+u2855*(u411&
 &3+u4072*u2855)+u414))+u2855*(u503+u507*u2855)+u696))
  u4218=-u491
  u4219=-u687
  u4220=-u815
  value(1)=value(1)+(c(79)*(x*(u2856*(u4218+u3971+u2856*(u2856*(u4220+u&
 &44)+u3972+u4219))+u3973+u3777)*z))
  u4218=4.80213173743530218d-2*u2858*(u73+u722)*u40
  u4219=-u3928
  u4220=-u489
  u4223=-u410
  u3987=5.28234491117883239d-1*u4141
  u946=-u756
  u3985=-3.5215632741192216d-2*u4086
  u4221=9.86037716753382047d-1*u947
  u4222=1.40862530964768864d-1*u947
  u4086=-u3933
  u945=-u4100
  u4012=2.11293796447153296d-1*u499
  u499=u2855*(u779*u2855+u4086)
  u4100=u2855*(u457+u4012)
  u949=u2855*(u945+u4061*u2855)
  u3933=u2855*(u637+u3929*u2855)
  value(2)=value(2)+(c(79)*((u4219+u2856*(u2856*(u946+u2856*(u420+u4221&
 &))+u4223))*u2857+u2856*(u4220+u499+u2856*(u2856*(u3985+u4100+(u457+u4&
 &222)*u2856)+u949+u3987))+u3933+u4218))
  u4218=u3828*u2855
  u4219=u703*u2855
  u4220=u664*u2855
  value(3)=value(3)+(c(79)*(u2985*((u732+u2856*(u2856*(u525+u3823)+u402&
 &4))*u2857+u2856*(u843+u4218+u2856*(u2856*(u859+u820)+u4219+u681))+u42&
 &20+u972)))
  u4024=u2855*(u470*u2855+u4047)
  u525=u2855*(u986+u428*u2855)
  u732=u2855*(u41+u425*u2855)
  value(4)=value(4)+(c(79)*(u2857*(u812+u2856*(u2856*(u3879+u2856*(u523&
 &+u863))+u4067)+(u2856*(u825+u2856*(u56+u811))+u622)*u2857)+u2856*(u41&
 &6+u4024+u2856*(u2856*(u862+u2855*(u4160+u4065))+u525+u826))+u732+u89)&
 &)
  u3985=u4004*u2856
  u4221=u381*u2855
  u4222=u4093*u2855
  value(5)=value(5)+(c(79)*(u2984*(u2857*(u479+u2856*(u96*u2856+u860)+(&
 &u2856*(u526+u3985)+u4130)*u2857)+u2856*(u5+u2855*(u4221+u4131)+u2856*&
 &(u2856*(u527+u4144+u522)+u2855*(u4057+u3830)+u861))+u2855*(u635+u4222&
 &)+u835)))
  value(6)=value(6)+(c(79)*u3763)
  u4223=-u82
  u3987=-u83
  u946=2.82353212234285965d-1*u4030
  u69=-5.6470642446857193d-1*u961
  u3989=2.25882569787428772d-1*u74
  u4224=-5.6470642446857193d-1*u3936
  value(7)=value(7)+(c(79)*(u2984*(u2857*(u3987+u2856*(u2856*(u3989+u38&
 &87)+u69)+(u3798+u3828)*u2857)+u2856*(u946+u2855*(u3867+u4224)+u2856*(&
 &u804*u2867+u2855*(u626+u765)+u432))+u2855*(u551+u3984)+u4223)))
  u4224=-u816
  u4223=-u817
  u3987=-1.05646898223576648d-1*u3787
  u3787=u2855*(u650*u2855+u438)
  u4095=u2855*(u749+u3758*u2855)
  value(8)=value(8)+(c(79)*(u2983*((u781+u2856*(u2856*(u819+u478)+u755)&
 &)*u2857+u2856*(u4223+u3787+u2856*(u3839*u2856+u913+u3987))+u4095+u422&
 &4)))
  value(9)=value(9)+(c(79)*(y*(u2856*(u4025+u2855*(u3974+u962)+u2856*(u&
 &2856*(u529+u822)+u2855*(u864+u4064)+u684))+u2855*(u3986+u3975)+u802)*&
 &z))
  u4224=-u799
  u4223=-u800
  u3987=-1.49407276290432138d-1*u973
  u946=1.99209701720576184d-1*u947
  u800=u2855*(u651*u2855+u439)
  u799=u2855*(u570+u796*u2855)
  value(1)=value(1)+(c(80)*(u2983*((u474+u2856*(u2856*(u803+u3821)+u413&
 &7))*u2857+u2856*(u4223+u800+u2856*((u946+u966)*u2856+u8+u3987))+u799+&
 &u4224)))
  u4224=-u477
  u4223=-u646
  u3987=-u797
  u946=-u3988
  u3988=-3.16940694670729944d-1*u973
  u3989=2.11293796447153296d-1*u469
  u69=-1.40862530964768864d-1*u4108
  u4108=u2856*(u45+u478)
  u3986=u69*u2856
  value(2)=value(2)+(c(80)*(u2984*(u2857*(u4223+u2856*(u2856*(u3989+u39&
 &86)+u946)+(u4108+u514)*u2857)+u2856*(u3987+u2856*(u753*u2856+u3988))+&
 &u4224)))
  u469=u2855*(u391*u2855+u711)
  u973=u2855*(u692+u681*u2855)
  value(3)=value(3)+(c(80)*(u2983*(u2857*(u4058+u2856*(u2856*(u540+u382&
 &2)+u4123)+(u2856*(u431+u3823)+u725)*u2857)+u2856*(u412+u469+u2856*((u&
 &541+u3894)*u2856+u415+u748))+u973+u476)))
  u3987=u693*u2856
  value(4)=value(4)+(c(80)*(u2984*(u2857*(u461+u2856*(u4040*u2856+u879)&
 &+(u2856*(u547+u3987)+u3854)*u2857)+u2856*(u57+u660*u2856)+u636)))
  value(5)=value(5)+(c(80)*(u2857*(u4112+u2856*(u2856*(u878+u2856*(u522&
 &+u546))+u3797)+u2857*((u4111+u2856*(u3985+u542))*u2857+u2856*(u989+u6&
 &5*u2856)+u3804))+u2856*(u3941+u2855*(u2855*(u543+u3812)+u4125)+u2856*&
 &(u2856*(u545+u2855*(u3781+u4043)+(u3882+u713)*u2856)+u2855*(u544+u285&
 &5*(u761+u841))+u4142))+u2855*(u721+u2855*(u3871+u3945))+u647))
  value(6)=value(6)+(c(80)*u3967)
  u3854=1.28342369197402711d-2*u2858*(u969-u409)*u40
  u3988=-2.54117891010857368d-1*u484
  u3989=1.69411927340571579d-1*u3938
  u484=-1.69411927340571579d-1*u3936
  u409=-4.23529818351428947d-1*u4054
  u4224=2.25882569787428772d-1*u4122
  u4122=-6.77647709362286316d-1*u960
  u960=-u4059
  u4106=-1.12941284893714386d-1*u3790
  u3790=5.6470642446857193d-2*u3843
  u4059=-4.70588687057143275d-2*u3936
  u3843=4.14118044610286082d-1*u947
  u4054=8.47059636702857895d-2*u649
  u946=2.82353212234285965d-2*u953
  u953=-1.12941284893714386d-1*u633
  u649=-1.8823547482285731d-2*u4104
  u4104=u2856*(u824+u627)
  u4223=u3774*u2855
  value(7)=value(7)+(c(80)*(u2857*(u3988+u2856*(u2856*(u4106+u2856*(u41&
 &29+u3843))+u4224)+u2857*((u484+u4104)*u2857+u2856*(u4122+u2856*(u3887&
 &+u3790))+u3989))+u2856*(u409+u2855*(u2855*(u946+u3780)+u3846)+u2856*(&
 &u2856*(u4059+u2855*(u3825+u978)+u4196)+u2855*(u946+u2855*(u3825+u649)&
 &)+u960))+u2855*(u4054+u2855*(u4223+u953))+u3854))
  u3790=-3.16940694670729944d-1*u656
  u4224=6.33881389341459887d-1*u4140
  u3988=-1.05646898223576648d-1*u4060
  u4060=u2855*(u648+u3977)+u3790
  u3977=u2855*(u3976+u441)
  value(8)=value(8)+(c(80)*(u2985*((u3929+u2856*(u2856*(u4061+u478)+u77&
 &9))*u2857+u2856*(u4224+u3977+u2856*((u753+u601)*u2856+u4147+u3988))+u&
 &4060)))
  u4224=u3900*u2855
  value(9)=value(9)+(c(80)*((u3897+u2856*(u2856*(u4005+u2856*(u739+u548&
 &))+u4133))*u2857+u2856*(u3769+u2855*(u2855*(u880+u3961)+u4039)+u2856*&
 &(u2856*(u882+u2855*(u3783+u3857)+(u3836+u549)*u2856)+u2855*(u881+u285&
 &5*(u3789+u768))+u55))+u2855*(u3824+u2855*(u4224+u967))+u813))
  u3769=-4.48221828871296415d-1*u764
  u3988=-u590
  u3989=-1.49407276290432138d-1*u959
  u959=u2855*(u3801+u3979)+u3769
  u3979=u2855*(u3978+u442)
  value(1)=value(1)+(c(81)*(u2985*((u506+u2856*(u2856*(u668+u3821)+u831&
 &))*u2857+u2856*(u3988+u3979+u2856*(u448*u2856+u3832+u3989))+u959)))
  u448=2.88127904246118131d-2*u4136
  u4136=-9.50822084012189831d-1*u4013
  u4013=1.05646898223576648d-1*u4141
  u4141=-u783
  u3988=-u678
  u3989=-6.33881389341459887d-1*u956
  u956=3.16940694670729944d-1*u3938
  u3938=-3.16940694670729944d-1*u996
  u996=1.05646898223576648d-1*u3896
  u590=8.45175185788613183d-1*u947
  value(2)=value(2)+(c(81)*(u2857*(u4136+u2856*(u2856*(u3938+u590*u2856&
 &)+u3988)+u2857*((u3919+u2856*(u478+u618))*u2857+u2856*(u3989+u2856*(u&
 &3986+u996))+u4013))+u2856*(u4141+u2856*(u3913*u2856+u956))+u448))
  u3913=u2855*(u594+u3980)+u4090
  value(3)=value(3)+(c(81)*(u2985*(u2857*(u818+u2856*(u2856*(u429+u3822&
 &)+u483)+(u2856*(u561+u3823)+u468)*u2857)+u2856*(u893+u3886+u2856*(u28&
 &55*(u3902+u3785)+u904))+u3913)))
  value(4)=value(4)+(c(81)*(u2857*(u807+u2856*(u602*u2856+u620)+u2857*(&
 &(u458+u2856*(u3987+u563))*u2857+u2856*(u909+u54*u2856)+u4225))+u2856*&
 &(u466+u3911*u2856)+u3903))
  value(5)=value(5)+(c(81)*(u2984*(u2857*(u4055+u2856*(u2856*(u465+u378&
 &4)+u907)+u2857*((u401+u3985)*u2857+u2856*(u4145+u482*u2856)+u905))+u2&
 &856*(u4109+u2855*(u3796+u562)+u2856*(u4179+u2855*(u617+u761)+u908))+u&
 &2855*(u690+u2855*(u3813+u906))+u460)))
  value(6)=value(6)+(c(81)*u3873)
  u2867=-8.47059636702857895d-2*u4164
  u4164=1.41176606117142983d-1*u4162
  u4162=-1.69411927340571579d-1*u961
  u3988=5.6470642446857193d-2*u4051
  u3989=-1.41176606117142983d-1*u3794
  u3794=5.6470642446857193d-2*u4034
  u4034=-u981
  u981=9.31765600373143684d-1*u947
  u590=-8.47059636702857895d-2*u4077
  u4077=2.82353212234285965d-2*u954
  u4225=-9.41177374114286549d-3*u632
  value(7)=value(7)+(c(81)*(u2984*(u2857*(u4164+u2856*(u2856*(u981+u412&
 &9)+u3989)+u2857*((u754+u824)*u2857+u2856*(u3794+u3887)+u4162))+u2856*&
 &(u3988+u2855*(u2855*(u4225+u3825)+u4077)+u2856*(u4200+u2855*(u978+u38&
 &25)+u4034))+u2855*(u590+u3916)+u2867)))
  u3988=3.16940694670729944d-1*u3927
  u3927=-u698
  u698=(u615+u601)*u2857
  value(8)=value(8)+(c(81)*(u2983*(u2857*(u3988+u3977+u2856*(u618*u2856&
 &+u3927)+u2857*(u698+u4108+u4147+u566))+u2856*(u585+u3919*u2856)+u4060&
 &)))
  u3988=u3765*u2856
  u3989=u979*u2856
  value(9)=value(9)+(c(81)*(u2984*(u2857*(u772+u2856*(u2856*(u3775+u398&
 &8)+u910)+(u2856*(u780+u3989)+u776)*u2857)+u2856*(u731+u2855*(u4081+u5&
 &64)+u2856*(u490*u2856+u2855*(u643+u3789)+u911))+u2855*(u4079+u584)+u4&
 &163)))
  u590=4.48221828871296415d-1*u473
  u473=-u77
  u77=u655*u2857
  value(1)=value(1)+(c(82)*(u2983*(u2857*(u590+u3979+u2856*(u611*u2856+&
 &u473)+u2857*(u77+u2856*(u651+u3821)+u3832+u571))+u2856*(u740+u798*u28&
 &56)+u959)))
  u590=-u4167
  u3827=5.28234491117883239d-1*u3970
  u3970=-1.05646898223576648d-1*u4015
  u4015=5.28234491117883239d-1*u3959
  u4225=-1.7607816370596108d-1*u394
  u394=3.5215632741192216d-2*u4115
  u4115=-u59
  value(2)=value(2)+(c(82)*(u2984*(u2857*(u3827+u2856*(u4115*u2856+u422&
 &5)+u2857*((u653+u478)*u2857+u2856*(u394+u3986)+u3970))+u2856*(u4015+u&
 &686*u2856)+u590)))
  u4225=-1.12941284893714386d0*u4116
  u590=3.95294497128000351d0*u3936
  u3827=-1.50588379858285848d0*u947
  u394=u3894+u3823
  value(3)=value(3)+(c(82)*(u2983*(u2857*(u682+u3886+u2856*(u702*u2856+&
 &u590)+u2857*((u388+u394)*u2857+u2856*(u3827+u3822)+u736+u927))+u2856*&
 &(u4225+u594*u2856)+u3913)))
  value(4)=value(4)+(c(82)*(u2984*(u2857*(u3862+u930*u2856+u2857*((u795&
 &+u3987)*u2857+u991*u2856+u929))+u3762*u2856+u583)))
  u3762=+u3851*(2.73954733374663945d6*u4026+u3895)*u40
  u3851=u4032*(3.92698688999642456d8*u4026+u621*(9.25d2+u3912))*u40
  u3912=+u608*(u471+u621*(u773+u593))*u40
  u471=-9.25996098457890084d-2*u947
  u593=u3795*(u81+u621*(6.5d1+u4019))*u40
  u4019=+u4062*(u422+u621*(u70+u4006))*u40
  u422=u3838*(1.7862483610443196d10*u4026+u621*(1.65d2*u3803+u727))*u40
  u3838=u3795*(1.45323007571593258d7*u4026+u3817)*u40
  u727=+u4062*(1.76183736145785534d8*u4026+u621*(4.15d2+5.4d1*pd2))*u40
  u3817=u3795*(4.80790557072535223d9*u4026+u621*(7.55d2*u3898+1.52d2*u4&
 &011))*u40
  u3898=+u3875*(1.43982443647814853d10*u4026+u621*(1.33d2*u3803+u3844))&
 &*u40
  u3844=+u4062*(-u3992*(pd2-1.0d1)+u488)*u40
  u590=u3795*(u3907-u621*(u3891-7.5d1*u4014))*u40
  u3795=+u3875*(-u621*(u3766-4.5d1*u788)+u475)*u40
  u4225=u43*u2856
  value(5)=value(5)+(c(82)*(u2857*(u4169+u2855*(u2855*(u590+u3812)+u727&
 &)+u2856*(u2856*(u4019+u4180)+u3851)+u2857*(u2857*(u928+u2855*(u4198+u&
 &3898)+u2856*(u4225+u471)+(u3818+u3811+u582)*u2857)+u2856*(u3912+u2856&
 &*(u3784+u422))+u2855*(u3817+u2855*(u761+u3795))+u4170))+u2856*(u3762+&
 &u2856*(u4181+u593))+u2855*(u3838+u2855*(u3871+u3844))+u3778))
  value(6)=value(6)+(c(82)*u3969)
  u3826=7.70054215184416268d-3*u2858*(u977-u4046)*u40
  u3773=-5.92941745692000526d-1*u733
  u3802=2.82353212234285965d-2*u3885
  u4046=-2.82353212234285965d-2*u435
  u3876=-1.24480369211622929d5*u3940
  u435=1.27058945505428684d1*u4116
  u3885=-4.65882800186571842d0*u3936
  u733=4.23529818351428947d-1*u4116
  u3969=7.05883030585714913d-1*u947
  u977=-8.47059636702857895d-2*u3859
  u3859=1.35529541872457263d0*u3921
  u3921=-2.82353212234285965d-2*u990
  u590=5.6470642446857193d-2*u48
  u48=5.6470642446857193d-2*u511
  u511=-2.82353212234285965d-2*u680
  u680=3.7647094964571462d-2*u983
  value(7)=value(7)+(c(82)*(u2857*(u3773+u2855*(u2855*(u511+u3780)+u385&
 &9)+u2856*(u2856*(u837+u3981)+u435)+u2857*(u2857*(u4046+u2855*(u4199+u&
 &590)+u4104+(u4022+u626)*u2857)+u2856*(u3885+u2856*(u4129+u3969))+u285&
 &5*(u3921+u2855*(u3825+u680))+u3802))+u2856*(u3876+u2856*(u4202+u733))&
 &+u2855*(u977+u2855*(u4223+u48))+u3826))
  u590=-u3906
  u3906=5.28234491117883239d-1*u4166
  u4166=-1.05646898223576648d-1*u998
  u998=-u487
  value(8)=value(8)+(c(82)*(u2985*(u2857*(u3906+u3787+u2856*(u45*u2856+&
 &u3927)+u2857*(u698+u2856*(u618+u478)+u913+u4166))+u2856*(u998+u3983)+&
 &u4095+u590)))
  value(9)=value(9)+(c(82)*(u2857*(u663+u2855*(u2855*(u405+u3961)+u3899&
 &)+u2856*(u2856*(u932+u440*u2856)+u3960)+u2857*((u957+u4080+u2856*(u39&
 &89+u402))*u2857+u2856*(u931+u2856*(u3988+u971))+u2855*(u405+u2855*(u3&
 &789+u600))+u629))+u2856*(u459+u2856*(u801*u2856+u3786))+u2855*(u603+u&
 &2855*(u4224+u744))+u3845))
  u590=7.47036381452160692d-1*u4183
  u4183=-1.49407276290432138d-1*u999
  u999=-u4184
  value(1)=value(1)+(c(83)*(u2985*(u2857*(u590+u800+u2856*(u651*u2856+u&
 &473)+u2857*(u77+u2856*(u611+u3821)+u8+u4183))+u2856*(u999+u796*u2856)&
 &+u799+u884)))
  u590=-u3800
  u3800=-u838
  u838=-u968
  u968=-u3889
  u3889=-u3820
  value(2)=value(2)+(c(83)*(u2857*(u3800+u499+u2856*(u3983+u998)+u2857*&
 &(u2857*(u968+u4100+u2856*(u420+u607)+(u457+u615)*u2857)+u2856*(u3889+&
 &u4203)+u949+u838))+u2856*(u786+u4204)+u3933+u590))
  value(3)=value(3)+(c(83)*(u2985*(u2857*(u407+u469+u2856*(u707*u2856+u&
 &3755)+u2857*((u456+u394)*u2857+u2856*(u586+u3822)+u415+u938))+u2856*(&
 &u950+u840*u2856)+u973+u383)))
  value(4)=value(4)+(c(83)*(u2857*(u3874+u4024+u2856*(u4205+u4193)+u285&
 &7*(u2857*(u942+u595+u2856*(u523+u443)+(u56+u3904+u588)*u2857)+u2856*(&
 &u4008+u4206)+u525+u3995))+u2856*(u4194+u4207)+u732+u550))
  value(5)=value(5)+(c(83)*(u2984*(u2857*(u480+u2855*(u4221+u3770)+u285&
 &6*(u941*u2856+u587)+u2857*(u2857*(u378+u3808+u4225+u3982)+u2856*(u376&
 &7+u3784)+u2855*(u940+u3830)+u939))+u2856*(u3831+u3948*u2856)+u2855*(u&
 &4190+u4222)+u4189)))
  value(6)=value(6)+(c(83)*u4010)
  u4008=-u3853
  u4010=-u4074
  u4074=-u685
  u685=-u3998
  u3998=-u743
  u743=-u4016
  u590=-u987
  u987=-u4033
  value(7)=value(7)+(c(83)*(u2984*(u2857*(u4010+u2855*(u3867+u987)+u285&
 &6*(u4208+u3998)+u2857*((u702+u497)*u2857+u2856*(u743+u4129)+u2855*(u4&
 &31+u765)+u4074))+u2856*(u685+u4209)+u2855*(u590+u3984)+u4008)))
  value(8)=value(8)+(c(83)*(u2983*(u2857*(u998+u3868+u779*u2856+u2857*(&
 &(u607+u496)*u2857+u4061*u2856+u4216+u3889))+u3929*u2856+u4217+u786)))
  value(9)=value(9)+(c(83)*(u2984*(u2857*(u3782+u2855*(u3974+u789)+u285&
 &6*(u589*u2856+u944)+u2857*((u610+u3989)*u2857+u2856*(u985+u3988)+u407&
 &8+u943))+u2856*(u467+u93*u2856)+u2855*(u3856+u3975)+u72)))
  value(1)=value(1)+(c(84)*(u2983*(u2857*(u3971+u831*u2856+u2857*((u44)&
 &*u2857+u668*u2856+u3972))+u506*u2856+u3973)))
  value(2)=value(2)+(c(84)*(u2984*(u2857*(u4210+u752*u2856+u2857*((u495&
 &)*u2857+u933*u2856+u4211))+u3878*u2856+u4212)))
  u3878=-4.23529818351428948d1*u4116
  u3858=2.20235505542743053d1*u3936
  u590=-3.04941469213028842d0*u947
  value(3)=value(3)+(c(84)*(u2983*(u2857*(u3878+u4218+u3828*u2856+u2857&
 &*(u2857*(u590+u820+u505*u2857)+u703*u2856+u4219+u3858))+u664*u2856+u4&
 &220+u3853)))
  u3922=1.09537136419378829d6*u3940
  u3858=-5.03127755072883509d1*u4116
  u2983=1.6770925169096117d1*u3936
  u590=-1.75695406533387892d0*u947
  u3878=(u2857*(u3858+u4213+u463*u2856+u2857*(u2857*(u590+u3869)+u659*u&
 &2856+u4214+u2983))+u794*u2856+u4215+u3922)
  value(4)=value(4)+(c(84)*(u2984*u3878))
  u4213=-1.72189357151260553d-2*u2858*(3.62998548665449682d3*exppd2+6.8&
 &5793968315383915d5*u4026)*u40
  u2984=4.54579902879327859d0*u2858*(u737-u2966)*u40
  u3869=-3.78816585732773216d-1*u2858
  u4214=+u3869*(-u621*(pd2-7.7d1)+u806)*u40
  u4215=u3840*(9.36108766750499044d8*u4026-u621*(u4011-1.47d2*u4014))*u&
 &40
  u4014=-6.43988195745714467d-1*u947
  u3828=-1.54636981900811541d5*u3940
  u806=1.95090541652378206d1*u4116
  u2858=-1.00386395219184902d1*u3936
  u737=1.37636692816240935d0*u947
  u2966=-9.4704146433193304d-2*u4116
  u3840=2.84112439299579912d-1*u3936
  u4011=-9.4704146433193304d-2*u947
  u3853=+u3869*(u411+u3992)*u40
  u3869=u454*(6.32563293631856497d7*u4026+u621*(1.49d2+u4017))*u40
  u4017=+u777*(u784+u621*(u462+u4006))*u40
  u4006=u4032*(1.54808191290507699d10*u4026+u621*(1.43d2*u3803+u645))*u&
 &40
  u3803=9.4704146433193304d-2*u4116
  u4026=-2.84112439299579912d-1*u3936
  u40=9.4704146433193304d-2*u947
  value(5)=value(5)+(c(84)*(u2857*(u2984+u2855*(u4026*u2855+u3869)+u285&
 &6*(u3840*u2856+u806)+u2857*(u2857*(u4215+u2855*(u761+u4006)+u2856*(u5&
 &22+u737)+u2857*(u3982+u695+u4201+u4014))+u2856*(u2858+u4011*u2856)+u2&
 &855*(u4017+u40*u2855)+u4214))+u2856*(u3828+u2966*u2856)+u2855*(u3853+&
 &u3803*u2855)+u4213))
  value(6)=value(6)+(c(84)*(u2985*u3878))
  u2983=-1.65973825615497239d5*u3940
  u659=+u2983
  u4139=2.11764909175714474d1*u4116
  u4021=-1.10117752771371526d1*u3936
  u3858=+u4021
  u590=1.52470734606514421d0*u947
  u3922=-1.41176606117142983d-1*u4116
  u3878=-u2983
  u2983=-u4139
  u794=-u4021
  u463=-u590
  u4021=1.41176606117142983d-1*u4116
  value(7)=value(7)+(c(84)*(u2857*(u2855*(u445*u2855+u2983)+u2856*(u408&
 &2*u2856+u4139)+u2857*(u2857*(u2855*(u765+u463)+u2856*(u4129+u590)+(u8&
 &24+u644)*u2857)+u2856*(u3858+u433*u2856)+u2855*(u794+u718*u2855)))+u2&
 &856*(u659+u3922*u2856)+u2855*(u3878+u4021*u2855)))
  value(8)=value(8)+(c(84)*(u2985*(u2857*(u76*u2855+u755*u2856+u2857*((&
 &u496)*u2857+u819*u2856+u538*u2855))+u781*u2856+u408*u2855)))
  u408=-3.65937888937172958d4*u3940
  u496=+u408
  u2985=2.68933097322777849d1*u830
  u830=-2.24110914435648207d0*u418
  u418=2.39051642064691421d0*u759
  u759=-u408
  u3940=-u72
  u4116=2.24110914435648207d0*u464
  u464=-1.12055457217824104d0*u3939
  u3939=1.12055457217824104d0*u964
  u964=-7.47036381452160691d-2*u3896
  u3896=3.92194100262384363d0*u3936
  u3936=-1.30731366754128121d0*u947
  value(9)=value(9)+(c(84)*(u2857*(u2985+u2855*(u3896*u2855+u464)+u2856&
 &*(u3942*u2856+u3940)+u2857*(u2857*(u418+u2855*(u3814+u964)+u2856*(u73&
 &9+u827)+(u3836+u827)*u2857)+u2856*(u776+u4072*u2856)+u2855*(u3939+u39&
 &36*u2855)+u830))+u2856*(u759+u507*u2856)+u2855*(u4116+u3759*u2855)+u4&
 &96))
  if ( lmax .eq. 6 ) return
  ! if we get here, something went wrong...
  D_X_Y_shell_opt_4=.false.
end function

  
end module
