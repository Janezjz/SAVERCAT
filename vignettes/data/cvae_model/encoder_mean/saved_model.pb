��
��
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring �
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape�"serve*2.0.02v2.0.0-rc2-26-g64c3d382ca8�

p
	e1/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape:
��*
shared_name	e1/kernel
i
e1/kernel/Read/ReadVariableOpReadVariableOp	e1/kernel*
dtype0* 
_output_shapes
:
��
g
e1/biasVarHandleOp*
shared_name	e1/bias*
dtype0*
_output_shapes
: *
shape:�
`
e1/bias/Read/ReadVariableOpReadVariableOpe1/bias*
dtype0*
_output_shapes	
:�
k
	be1/gammaVarHandleOp*
shared_name	be1/gamma*
dtype0*
_output_shapes
: *
shape:�
d
be1/gamma/Read/ReadVariableOpReadVariableOp	be1/gamma*
dtype0*
_output_shapes	
:�
w
be1/moving_meanVarHandleOp*
dtype0*
_output_shapes
: *
shape:�* 
shared_namebe1/moving_mean
p
#be1/moving_mean/Read/ReadVariableOpReadVariableOpbe1/moving_mean*
dtype0*
_output_shapes	
:�

be1/moving_varianceVarHandleOp*
shape:�*$
shared_namebe1/moving_variance*
dtype0*
_output_shapes
: 
x
'be1/moving_variance/Read/ReadVariableOpReadVariableOpbe1/moving_variance*
dtype0*
_output_shapes	
:�
p
	e2/kernelVarHandleOp*
shape:
��*
shared_name	e2/kernel*
dtype0*
_output_shapes
: 
i
e2/kernel/Read/ReadVariableOpReadVariableOp	e2/kernel*
dtype0* 
_output_shapes
:
��
g
e2/biasVarHandleOp*
shape:�*
shared_name	e2/bias*
dtype0*
_output_shapes
: 
`
e2/bias/Read/ReadVariableOpReadVariableOpe2/bias*
dtype0*
_output_shapes	
:�
k
	be2/gammaVarHandleOp*
dtype0*
_output_shapes
: *
shape:�*
shared_name	be2/gamma
d
be2/gamma/Read/ReadVariableOpReadVariableOp	be2/gamma*
dtype0*
_output_shapes	
:�
w
be2/moving_meanVarHandleOp*
shape:�* 
shared_namebe2/moving_mean*
dtype0*
_output_shapes
: 
p
#be2/moving_mean/Read/ReadVariableOpReadVariableOpbe2/moving_mean*
dtype0*
_output_shapes	
:�

be2/moving_varianceVarHandleOp*
shape:�*$
shared_namebe2/moving_variance*
dtype0*
_output_shapes
: 
x
'be2/moving_variance/Read/ReadVariableOpReadVariableOpbe2/moving_variance*
dtype0*
_output_shapes	
:�
w
z_mean/kernelVarHandleOp*
shape:	�*
shared_namez_mean/kernel*
dtype0*
_output_shapes
: 
p
!z_mean/kernel/Read/ReadVariableOpReadVariableOpz_mean/kernel*
dtype0*
_output_shapes
:	�
n
z_mean/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namez_mean/bias
g
z_mean/bias/Read/ReadVariableOpReadVariableOpz_mean/bias*
dtype0*
_output_shapes
:
h
bz/gammaVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_name
bz/gamma
a
bz/gamma/Read/ReadVariableOpReadVariableOpbz/gamma*
dtype0*
_output_shapes
:
t
bz/moving_meanVarHandleOp*
shared_namebz/moving_mean*
dtype0*
_output_shapes
: *
shape:
m
"bz/moving_mean/Read/ReadVariableOpReadVariableOpbz/moving_mean*
dtype0*
_output_shapes
:
|
bz/moving_varianceVarHandleOp*
dtype0*
_output_shapes
: *
shape:*#
shared_namebz/moving_variance
u
&bz/moving_variance/Read/ReadVariableOpReadVariableOpbz/moving_variance*
dtype0*
_output_shapes
:

NoOpNoOp
�,
ConstConst"/device:CPU:0*�+
value�+B�+ B�+
�
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer_with_weights-1
layer-4
layer-5
layer_with_weights-2
layer-6
layer_with_weights-3
layer-7
	layer_with_weights-4
	layer-8

layer_with_weights-5

layer-9
trainable_variables
	variables
regularization_losses
	keras_api

signatures
R
trainable_variables
	variables
regularization_losses
	keras_api
R
trainable_variables
	variables
regularization_losses
	keras_api
R
trainable_variables
	variables
regularization_losses
	keras_api
x

activation

kernel
bias
trainable_variables
 	variables
!regularization_losses
"	keras_api
�
#axis
	$gamma
%moving_mean
&moving_variance
'trainable_variables
(	variables
)regularization_losses
*	keras_api
R
+trainable_variables
,	variables
-regularization_losses
.	keras_api
x

activation

/kernel
0bias
1trainable_variables
2	variables
3regularization_losses
4	keras_api
�
5axis
	6gamma
7moving_mean
8moving_variance
9trainable_variables
:	variables
;regularization_losses
<	keras_api
x

activation

=kernel
>bias
?trainable_variables
@	variables
Aregularization_losses
B	keras_api
�
Caxis
	Dgamma
Emoving_mean
Fmoving_variance
Gtrainable_variables
H	variables
Iregularization_losses
J	keras_api
?
0
1
$2
/3
04
65
=6
>7
D8
n
0
1
$2
%3
&4
/5
06
67
78
89
=10
>11
D12
E13
F14
 
�
Klayer_regularization_losses
trainable_variables

Llayers
	variables
Mnon_trainable_variables
Nmetrics
regularization_losses
 
 
 
 
�
Olayer_regularization_losses
trainable_variables

Players
	variables
Qnon_trainable_variables
Rmetrics
regularization_losses
 
 
 
�
Slayer_regularization_losses
trainable_variables

Tlayers
	variables
Unon_trainable_variables
Vmetrics
regularization_losses
 
 
 
�
Wlayer_regularization_losses
trainable_variables

Xlayers
	variables
Ynon_trainable_variables
Zmetrics
regularization_losses
R
[trainable_variables
\	variables
]regularization_losses
^	keras_api
US
VARIABLE_VALUE	e1/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEe1/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
�
_layer_regularization_losses
trainable_variables

`layers
 	variables
anon_trainable_variables
bmetrics
!regularization_losses
 
TR
VARIABLE_VALUE	be1/gamma5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUE
`^
VARIABLE_VALUEbe1/moving_mean;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
hf
VARIABLE_VALUEbe1/moving_variance?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

$0

$0
%1
&2
 
�
clayer_regularization_losses
'trainable_variables

dlayers
(	variables
enon_trainable_variables
fmetrics
)regularization_losses
 
 
 
�
glayer_regularization_losses
+trainable_variables

hlayers
,	variables
inon_trainable_variables
jmetrics
-regularization_losses
US
VARIABLE_VALUE	e2/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEe2/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE

/0
01

/0
01
 
�
klayer_regularization_losses
1trainable_variables

llayers
2	variables
mnon_trainable_variables
nmetrics
3regularization_losses
 
TR
VARIABLE_VALUE	be2/gamma5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUE
`^
VARIABLE_VALUEbe2/moving_mean;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
hf
VARIABLE_VALUEbe2/moving_variance?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

60

60
71
82
 
�
olayer_regularization_losses
9trainable_variables

players
:	variables
qnon_trainable_variables
rmetrics
;regularization_losses
YW
VARIABLE_VALUEz_mean/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEz_mean/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE

=0
>1

=0
>1
 
�
slayer_regularization_losses
?trainable_variables

tlayers
@	variables
unon_trainable_variables
vmetrics
Aregularization_losses
 
SQ
VARIABLE_VALUEbz/gamma5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUE
_]
VARIABLE_VALUEbz/moving_mean;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
ge
VARIABLE_VALUEbz/moving_variance?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

D0

D0
E1
F2
 
�
wlayer_regularization_losses
Gtrainable_variables

xlayers
H	variables
ynon_trainable_variables
zmetrics
Iregularization_losses
 
F
0
1
2
3
4
5
6
7
	8

9
*
%0
&1
72
83
E4
F5
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
�
{layer_regularization_losses
[trainable_variables

|layers
\	variables
}non_trainable_variables
~metrics
]regularization_losses
 

0
 
 
 
 

%0
&1
 
 
 
 
 
 

0
 
 
 
 

70
81
 
 

0
 
 
 
 

E0
F1
 
 
 
 
 *
dtype0*
_output_shapes
: 
t
serving_default_BPlaceholder*
dtype0*'
_output_shapes
:���������*
shape:���������
|
serving_default_x_inputPlaceholder*
shape:����������*
dtype0*(
_output_shapes
:����������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_Bserving_default_x_input	e1/kernele1/biasbe1/moving_variance	be1/gammabe1/moving_mean	e2/kernele2/biasbe2/moving_variance	be2/gammabe2/moving_meanz_mean/kernelz_mean/biasbz/moving_variancebz/gammabz/moving_mean*+
_gradient_op_typePartitionedCall-4725*+
f&R$
"__inference_signature_wrapper_4110*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������
O
saver_filenamePlaceholder*
dtype0*
_output_shapes
: *
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenamee1/kernel/Read/ReadVariableOpe1/bias/Read/ReadVariableOpbe1/gamma/Read/ReadVariableOp#be1/moving_mean/Read/ReadVariableOp'be1/moving_variance/Read/ReadVariableOpe2/kernel/Read/ReadVariableOpe2/bias/Read/ReadVariableOpbe2/gamma/Read/ReadVariableOp#be2/moving_mean/Read/ReadVariableOp'be2/moving_variance/Read/ReadVariableOp!z_mean/kernel/Read/ReadVariableOpz_mean/bias/Read/ReadVariableOpbz/gamma/Read/ReadVariableOp"bz/moving_mean/Read/ReadVariableOp&bz/moving_variance/Read/ReadVariableOpConst*+
_gradient_op_typePartitionedCall-4762*&
f!R
__inference__traced_save_4761*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*
_output_shapes
: 
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filename	e1/kernele1/bias	be1/gammabe1/moving_meanbe1/moving_variance	e2/kernele2/bias	be2/gammabe2/moving_meanbe2/moving_variancez_mean/kernelz_mean/biasbz/gammabz/moving_meanbz/moving_variance*+
_gradient_op_typePartitionedCall-4820*)
f$R"
 __inference__traced_restore_4819*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*
_output_shapes
: ��	
�'
�
<__inference_bz_layer_call_and_return_conditional_losses_3696

inputs"
assignmovingavg_bz_moving_mean(
$assignmovingavg_1_bz_moving_variance)
%batchnorm_mul_readvariableop_bz_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: 
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
	keep_dims(*
T0*
_output_shapes

:d
moments/StopGradientStopGradientmoments/mean:output:0*
_output_shapes

:*
T0�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
	keep_dims(*
T0*
_output_shapes

:m
moments/SqueezeSqueezemoments/mean:output:0*
_output_shapes
:*
squeeze_dims
 *
T0s
moments/Squeeze_1Squeezemoments/variance:output:0*
_output_shapes
:*
squeeze_dims
 *
T0�
AssignMovingAvg/decayConst*
valueB
 *
�#<*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
: y
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_bz_moving_mean*
dtype0*
_output_shapes
:�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
_output_shapes
:�
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
_output_shapes
:*
T0*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean�
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_bz_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
dtype0*
_output_shapes
: *
valueB
 *
�#<*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp$assignmovingavg_1_bz_moving_variance*
dtype0*
_output_shapes
:�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
_output_shapes
:�
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
_output_shapes
:*
T0*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance�
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp$assignmovingavg_1_bz_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
dtype0*
_output_shapes
 T
batchnorm/add/yConst*
dtype0*
_output_shapes
: *
valueB
 *o�:q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
_output_shapes
:*
T0~
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:t
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:���������S
batchnorm/NegNegmoments/Squeeze:output:0*
_output_shapes
:*
T0a
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes
:t
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�
�
!__inference_bz_layer_call_fn_4682

inputs*
&statefulpartitionedcall_bz_moving_mean.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs&statefulpartitionedcall_bz_moving_mean*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3697*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3696*
Tout
2**
config_proto

CPU

GPU 2J 8�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�
q
E__inference_concatenate_layer_call_and_return_conditional_losses_4379
inputs_0
inputs_1
identityM
concat/axisConst*
value	B :*
dtype0*
_output_shapes
: x
concatConcatV2inputs_0inputs_1concat/axis:output:0*
T0*
N*(
_output_shapes
:����������X
IdentityIdentityconcat:output:0*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*:
_input_shapes)
':����������:���������:( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1
�<
�
 __inference__traced_restore_4819
file_prefix
assignvariableop_e1_kernel
assignvariableop_1_e1_bias 
assignvariableop_2_be1_gamma&
"assignvariableop_3_be1_moving_mean*
&assignvariableop_4_be1_moving_variance 
assignvariableop_5_e2_kernel
assignvariableop_6_e2_bias 
assignvariableop_7_be2_gamma&
"assignvariableop_8_be2_moving_mean*
&assignvariableop_9_be2_moving_variance%
!assignvariableop_10_z_mean_kernel#
assignvariableop_11_z_mean_bias 
assignvariableop_12_bz_gamma&
"assignvariableop_13_bz_moving_mean*
&assignvariableop_14_bz_moving_variance
identity_16��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_2�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�	RestoreV2�RestoreV2_1�
RestoreV2/tensor_namesConst"/device:CPU:0*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:�
RestoreV2/shape_and_slicesConst"/device:CPU:0*1
value(B&B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:�
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*P
_output_shapes>
<:::::::::::::::*
dtypes
2L
IdentityIdentityRestoreV2:tensors:0*
T0*
_output_shapes
:v
AssignVariableOpAssignVariableOpassignvariableop_e1_kernelIdentity:output:0*
dtype0*
_output_shapes
 N

Identity_1IdentityRestoreV2:tensors:1*
T0*
_output_shapes
:z
AssignVariableOp_1AssignVariableOpassignvariableop_1_e1_biasIdentity_1:output:0*
dtype0*
_output_shapes
 N

Identity_2IdentityRestoreV2:tensors:2*
T0*
_output_shapes
:|
AssignVariableOp_2AssignVariableOpassignvariableop_2_be1_gammaIdentity_2:output:0*
dtype0*
_output_shapes
 N

Identity_3IdentityRestoreV2:tensors:3*
_output_shapes
:*
T0�
AssignVariableOp_3AssignVariableOp"assignvariableop_3_be1_moving_meanIdentity_3:output:0*
dtype0*
_output_shapes
 N

Identity_4IdentityRestoreV2:tensors:4*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOp&assignvariableop_4_be1_moving_varianceIdentity_4:output:0*
dtype0*
_output_shapes
 N

Identity_5IdentityRestoreV2:tensors:5*
T0*
_output_shapes
:|
AssignVariableOp_5AssignVariableOpassignvariableop_5_e2_kernelIdentity_5:output:0*
dtype0*
_output_shapes
 N

Identity_6IdentityRestoreV2:tensors:6*
T0*
_output_shapes
:z
AssignVariableOp_6AssignVariableOpassignvariableop_6_e2_biasIdentity_6:output:0*
dtype0*
_output_shapes
 N

Identity_7IdentityRestoreV2:tensors:7*
_output_shapes
:*
T0|
AssignVariableOp_7AssignVariableOpassignvariableop_7_be2_gammaIdentity_7:output:0*
dtype0*
_output_shapes
 N

Identity_8IdentityRestoreV2:tensors:8*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOp"assignvariableop_8_be2_moving_meanIdentity_8:output:0*
dtype0*
_output_shapes
 N

Identity_9IdentityRestoreV2:tensors:9*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp&assignvariableop_9_be2_moving_varianceIdentity_9:output:0*
dtype0*
_output_shapes
 P
Identity_10IdentityRestoreV2:tensors:10*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp!assignvariableop_10_z_mean_kernelIdentity_10:output:0*
dtype0*
_output_shapes
 P
Identity_11IdentityRestoreV2:tensors:11*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpassignvariableop_11_z_mean_biasIdentity_11:output:0*
dtype0*
_output_shapes
 P
Identity_12IdentityRestoreV2:tensors:12*
T0*
_output_shapes
:~
AssignVariableOp_12AssignVariableOpassignvariableop_12_bz_gammaIdentity_12:output:0*
dtype0*
_output_shapes
 P
Identity_13IdentityRestoreV2:tensors:13*
_output_shapes
:*
T0�
AssignVariableOp_13AssignVariableOp"assignvariableop_13_bz_moving_meanIdentity_13:output:0*
dtype0*
_output_shapes
 P
Identity_14IdentityRestoreV2:tensors:14*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp&assignvariableop_14_bz_moving_varianceIdentity_14:output:0*
dtype0*
_output_shapes
 �
RestoreV2_1/tensor_namesConst"/device:CPU:0*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPH*
dtype0*
_output_shapes
:t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:�
RestoreV2_1	RestoreV2file_prefix!RestoreV2_1/tensor_names:output:0%RestoreV2_1/shape_and_slices:output:0
^RestoreV2"/device:CPU:0*
dtypes
2*
_output_shapes
:1
NoOpNoOp"/device:CPU:0*
_output_shapes
 �
Identity_15Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: �
Identity_16IdentityIdentity_15:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9
^RestoreV2^RestoreV2_1*
T0*
_output_shapes
: "#
identity_16Identity_16:output:0*Q
_input_shapes@
>: :::::::::::::::2*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122
RestoreV2_1RestoreV2_12*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92
	RestoreV2	RestoreV2:+ '
%
_user_specified_namefile_prefix: : : : : : : : :	 :
 : : : : : 
�(
�
?__inference_model_layer_call_and_return_conditional_losses_3983
x_input
b(
$e1_statefulpartitionedcall_e1_kernel&
"e1_statefulpartitionedcall_e1_bias3
/be1_statefulpartitionedcall_be1_moving_variance)
%be1_statefulpartitionedcall_be1_gamma/
+be1_statefulpartitionedcall_be1_moving_mean(
$e2_statefulpartitionedcall_e2_kernel&
"e2_statefulpartitionedcall_e2_bias3
/be2_statefulpartitionedcall_be2_moving_variance)
%be2_statefulpartitionedcall_be2_gamma/
+be2_statefulpartitionedcall_be2_moving_mean0
,z_mean_statefulpartitionedcall_z_mean_kernel.
*z_mean_statefulpartitionedcall_z_mean_bias1
-bz_statefulpartitionedcall_bz_moving_variance'
#bz_statefulpartitionedcall_bz_gamma-
)bz_statefulpartitionedcall_bz_moving_mean
identity��be1/StatefulPartitionedCall�be2/StatefulPartitionedCall�bz/StatefulPartitionedCall�e1/StatefulPartitionedCall�e2/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
concatenate/PartitionedCallPartitionedCallx_inputb**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3755*N
fIRG
E__inference_concatenate_layer_call_and_return_conditional_losses_3748*
Tout
2�
e1/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0$e1_statefulpartitionedcall_e1_kernel"e1_statefulpartitionedcall_e1_bias*+
_gradient_op_typePartitionedCall-3779*E
f@R>
<__inference_e1_layer_call_and_return_conditional_losses_3773*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
be1/StatefulPartitionedCallStatefulPartitionedCall#e1/StatefulPartitionedCall:output:0/be1_statefulpartitionedcall_be1_moving_variance%be1_statefulpartitionedcall_be1_gamma+be1_statefulpartitionedcall_be1_moving_mean*+
_gradient_op_typePartitionedCall-3454*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3453*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2�
dropout/PartitionedCallPartitionedCall$be1/StatefulPartitionedCall:output:0*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3850*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3838*
Tout
2**
config_proto

CPU

GPU 2J 8�
e2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0$e2_statefulpartitionedcall_e2_kernel"e2_statefulpartitionedcall_e2_bias*+
_gradient_op_typePartitionedCall-3872*E
f@R>
<__inference_e2_layer_call_and_return_conditional_losses_3866*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2�
be2/StatefulPartitionedCallStatefulPartitionedCall#e2/StatefulPartitionedCall:output:0/be2_statefulpartitionedcall_be2_moving_variance%be2_statefulpartitionedcall_be2_gamma+be2_statefulpartitionedcall_be2_moving_mean*+
_gradient_op_typePartitionedCall-3592*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3591*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
z_mean/StatefulPartitionedCallStatefulPartitionedCall$be2/StatefulPartitionedCall:output:0,z_mean_statefulpartitionedcall_z_mean_kernel*z_mean_statefulpartitionedcall_z_mean_bias**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3920*I
fDRB
@__inference_z_mean_layer_call_and_return_conditional_losses_3914*
Tout
2�
bz/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0-bz_statefulpartitionedcall_bz_moving_variance#bz_statefulpartitionedcall_bz_gamma)bz_statefulpartitionedcall_bz_moving_mean*+
_gradient_op_typePartitionedCall-3730*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3729*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:����������
IdentityIdentity#bz/StatefulPartitionedCall:output:0^be1/StatefulPartitionedCall^be2/StatefulPartitionedCall^bz/StatefulPartitionedCall^e1/StatefulPartitionedCall^e2/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2:
be1/StatefulPartitionedCallbe1/StatefulPartitionedCall2:
be2/StatefulPartitionedCallbe2/StatefulPartitionedCall28
bz/StatefulPartitionedCallbz/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall28
e1/StatefulPartitionedCalle1/StatefulPartitionedCall28
e2/StatefulPartitionedCalle2/StatefulPartitionedCall: : : : : :	 :
 : : : : : : :' #
!
_user_specified_name	x_input:!

_user_specified_nameB: : 
�
�
$__inference_model_layer_call_fn_4351
inputs_0
inputs_1%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias+
'statefulpartitionedcall_be1_moving_mean/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias+
'statefulpartitionedcall_be2_moving_mean/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias*
&statefulpartitionedcall_bz_moving_mean.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias'statefulpartitionedcall_be1_moving_mean+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias'statefulpartitionedcall_be2_moving_mean+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias&statefulpartitionedcall_bz_moving_mean*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma*+
_gradient_op_typePartitionedCall-4016*H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_4015*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall: : : :	 :
 : : : : : : :( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1: : : : 
�
�
$__inference_model_layer_call_fn_4372
inputs_0
inputs_1%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma+
'statefulpartitionedcall_be1_moving_mean%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma+
'statefulpartitionedcall_be2_moving_mean)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma*
&statefulpartitionedcall_bz_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma'statefulpartitionedcall_be1_moving_mean!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma'statefulpartitionedcall_be2_moving_mean%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma&statefulpartitionedcall_bz_moving_mean**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*+
_gradient_op_typePartitionedCall-4069*H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_4068*
Tout
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall: : : : : : : :	 :
 : : : : : : :( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1
�	
�
<__inference_e2_layer_call_and_return_conditional_losses_3866

inputs#
matmul_readvariableop_e2_kernel"
biasadd_readvariableop_e2_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpw
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_e2_kernel*
dtype0* 
_output_shapes
:
��j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_e2_bias*
dtype0*
_output_shapes	
:�w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������n
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:�����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�'
�
__inference__traced_save_4761
file_prefix(
$savev2_e1_kernel_read_readvariableop&
"savev2_e1_bias_read_readvariableop(
$savev2_be1_gamma_read_readvariableop.
*savev2_be1_moving_mean_read_readvariableop2
.savev2_be1_moving_variance_read_readvariableop(
$savev2_e2_kernel_read_readvariableop&
"savev2_e2_bias_read_readvariableop(
$savev2_be2_gamma_read_readvariableop.
*savev2_be2_moving_mean_read_readvariableop2
.savev2_be2_moving_variance_read_readvariableop,
(savev2_z_mean_kernel_read_readvariableop*
&savev2_z_mean_bias_read_readvariableop'
#savev2_bz_gamma_read_readvariableop-
)savev2_bz_moving_mean_read_readvariableop1
-savev2_bz_moving_variance_read_readvariableop
savev2_1_const

identity_1��MergeV2Checkpoints�SaveV2�SaveV2_1�
StringJoin/inputs_1Const"/device:CPU:0*<
value3B1 B+_temp_52d141997efe4f7db51688bd00b3253b/part*
dtype0*
_output_shapes
: s

StringJoin
StringJoinfile_prefixStringJoin/inputs_1:output:0"/device:CPU:0*
N*
_output_shapes
: L

num_shardsConst*
value	B :*
dtype0*
_output_shapes
: f
ShardedFilename/shardConst"/device:CPU:0*
value	B : *
dtype0*
_output_shapes
: �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*�
value�B�B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUE�
SaveV2/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*1
value(B&B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0$savev2_e1_kernel_read_readvariableop"savev2_e1_bias_read_readvariableop$savev2_be1_gamma_read_readvariableop*savev2_be1_moving_mean_read_readvariableop.savev2_be1_moving_variance_read_readvariableop$savev2_e2_kernel_read_readvariableop"savev2_e2_bias_read_readvariableop$savev2_be2_gamma_read_readvariableop*savev2_be2_moving_mean_read_readvariableop.savev2_be2_moving_variance_read_readvariableop(savev2_z_mean_kernel_read_readvariableop&savev2_z_mean_bias_read_readvariableop#savev2_bz_gamma_read_readvariableop)savev2_bz_moving_mean_read_readvariableop-savev2_bz_moving_variance_read_readvariableop"/device:CPU:0*
dtypes
2*
_output_shapes
 h
ShardedFilename_1/shardConst"/device:CPU:0*
value	B :*
dtype0*
_output_shapes
: �
ShardedFilename_1ShardedFilenameStringJoin:output:0 ShardedFilename_1/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2_1/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*1
value(B&B_CHECKPOINTABLE_OBJECT_GRAPHq
SaveV2_1/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B �
SaveV2_1SaveV2ShardedFilename_1:filename:0SaveV2_1/tensor_names:output:0"SaveV2_1/shape_and_slices:output:0savev2_1_const^SaveV2"/device:CPU:0*
dtypes
2*
_output_shapes
 �
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0ShardedFilename_1:filename:0^SaveV2	^SaveV2_1"/device:CPU:0*
N*
_output_shapes
:*
T0�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix	^SaveV2_1"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
_output_shapes
: *
T0s

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints^SaveV2	^SaveV2_1*
T0*
_output_shapes
: "!

identity_1Identity_1:output:0*�
_input_shapesy
w: :
��:�:�:�:�:
��:�:�:�:�:	�::::: 2
SaveV2_1SaveV2_12(
MergeV2CheckpointsMergeV2Checkpoints2
SaveV2SaveV2: : : : : : : :	 :
 : : : : : : :+ '
%
_user_specified_namefile_prefix: 
�
�
$__inference_model_layer_call_fn_4034
x_input
b%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias+
'statefulpartitionedcall_be1_moving_mean/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias+
'statefulpartitionedcall_be2_moving_mean/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias*
&statefulpartitionedcall_bz_moving_mean.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallx_inputb!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias'statefulpartitionedcall_be1_moving_mean+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias'statefulpartitionedcall_be2_moving_mean+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias&statefulpartitionedcall_bz_moving_mean*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*+
_gradient_op_typePartitionedCall-4016*H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_4015�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall: : : : : :' #
!
_user_specified_name	x_input:!

_user_specified_nameB: : : : : : : :	 :
 : 
�(
�
=__inference_be1_layer_call_and_return_conditional_losses_4438

inputs#
assignmovingavg_be1_moving_mean)
%assignmovingavg_1_be1_moving_variance*
&batchnorm_mul_readvariableop_be1_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 ZN
LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(e
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes
:	��
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*(
_output_shapes
:����������*
T0l
"moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
_output_shapes
:	�*
	keep_dims(*
T0n
moments/SqueezeSqueezemoments/mean:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:�t
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes	
:�*
squeeze_dims
 �
AssignMovingAvg/decayConst*
valueB
 *
�#<*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
: {
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_be1_moving_mean*
dtype0*
_output_shapes	
:��
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
_output_shapes	
:��
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
_output_shapes	
:��
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_be1_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
valueB
 *
�#<*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance*
dtype0*
_output_shapes
: �
 AssignMovingAvg_1/ReadVariableOpReadVariableOp%assignmovingavg_1_be1_moving_variance*
dtype0*
_output_shapes	
:��
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance*
_output_shapes	
:��
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
_output_shapes	
:�*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance�
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp%assignmovingavg_1_be1_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
dtype0*
_output_shapes
 *8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_varianceT
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: r
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
_output_shapes	
:�*
T0Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:�d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:����������T
batchnorm/NegNegmoments/Squeeze:output:0*
T0*
_output_shapes	
:�b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
_output_shapes	
:�*
T0u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*(
_output_shapes
:����������*
T0�
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp: : :& "
 
_user_specified_nameinputs: 
�(
�
=__inference_be2_layer_call_and_return_conditional_losses_4563

inputs#
assignmovingavg_be2_moving_mean)
%assignmovingavg_1_be2_moving_variance*
&batchnorm_mul_readvariableop_be2_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(e
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes
:	��
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*(
_output_shapes
:����������l
"moments/variance/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(n
moments/SqueezeSqueezemoments/mean:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:�t
moments/Squeeze_1Squeezemoments/variance:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:��
AssignMovingAvg/decayConst*
valueB
 *
�#<*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
: {
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_be2_moving_mean*
dtype0*
_output_shapes	
:��
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
_output_shapes	
:�*
T0*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean�
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
_output_shapes	
:��
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_be2_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
valueB
 *
�#<*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance*
dtype0*
_output_shapes
: �
 AssignMovingAvg_1/ReadVariableOpReadVariableOp%assignmovingavg_1_be2_moving_variance*
dtype0*
_output_shapes	
:��
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance*
_output_shapes	
:��
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance*
_output_shapes	
:��
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp%assignmovingavg_1_be2_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
dtype0*
_output_shapes
 *8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_varianceT
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: r
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
_output_shapes	
:�*
T0Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:�d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:����������T
batchnorm/NegNegmoments/Squeeze:output:0*
T0*
_output_shapes	
:�b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�
_
&__inference_dropout_layer_call_fn_4505

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3842*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3831�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*'
_input_shapes
:����������22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs
�
�
=__inference_be2_layer_call_and_return_conditional_losses_4584

inputs0
,batchnorm_readvariableop_be2_moving_variance*
&batchnorm_mul_readvariableop_be2_gamma.
*batchnorm_readvariableop_1_be2_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp,batchnorm_readvariableop_be2_moving_variance*
dtype0*
_output_shapes	
:�T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: x
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:�d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
batchnorm/ReadVariableOp_1ReadVariableOp*batchnorm_readvariableop_1_be2_moving_mean*
dtype0*
_output_shapes	
:�^
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
_output_shapes	
:�*
T0b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*3
_input_shapes"
 :����������:::24
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_1:& "
 
_user_specified_nameinputs: : : 
�'
�
<__inference_bz_layer_call_and_return_conditional_losses_4653

inputs"
assignmovingavg_bz_moving_mean(
$assignmovingavg_1_bz_moving_variance)
%batchnorm_mul_readvariableop_bz_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
	keep_dims(*
T0*
_output_shapes

:m
moments/SqueezeSqueezemoments/mean:output:0*
_output_shapes
:*
squeeze_dims
 *
T0s
moments/Squeeze_1Squeezemoments/variance:output:0*
squeeze_dims
 *
T0*
_output_shapes
:�
AssignMovingAvg/decayConst*
valueB
 *
�#<*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
: y
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_bz_moving_mean*
dtype0*
_output_shapes
:�
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
_output_shapes
:*
T0*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean�
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
_output_shapes
:�
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_bz_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*1
_class'
%#loc:@AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
valueB
 *
�#<*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
dtype0*
_output_shapes
: �
 AssignMovingAvg_1/ReadVariableOpReadVariableOp$assignmovingavg_1_bz_moving_variance*
dtype0*
_output_shapes
:�
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
_output_shapes
:�
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
_output_shapes
:�
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp$assignmovingavg_1_bz_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*7
_class-
+)loc:@AssignMovingAvg_1/bz/moving_variance*
dtype0*
_output_shapes
 T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: q
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
_output_shapes
:*
T0~
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:t
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
_output_shapes
:*
T0c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*'
_output_shapes
:���������*
T0S
batchnorm/NegNegmoments/Squeeze:output:0*
T0*
_output_shapes
:a
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes
:t
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp: : :& "
 
_user_specified_nameinputs: 
�
B
&__inference_dropout_layer_call_fn_4510

inputs
identity�
PartitionedCallPartitionedCallinputs**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3850*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3838*
Tout
2a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*'
_input_shapes
:����������:& "
 
_user_specified_nameinputs
�
�
"__inference_be2_layer_call_fn_4600

inputs/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma+
'statefulpartitionedcall_be2_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma'statefulpartitionedcall_be2_moving_mean*+
_gradient_op_typePartitionedCall-3592*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3591*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*3
_input_shapes"
 :����������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�
`
A__inference_dropout_layer_call_and_return_conditional_losses_4495

inputs
identity�Q
dropout/rateConst*
valueB
 *���=*
dtype0*
_output_shapes
: C
dropout/ShapeShapeinputs*
_output_shapes
:*
T0_
dropout/random_uniform/minConst*
valueB
 *    *
dtype0*
_output_shapes
: _
dropout/random_uniform/maxConst*
valueB
 *  �?*
dtype0*
_output_shapes
: �
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*
dtype0*(
_output_shapes
:�����������
dropout/random_uniform/subSub#dropout/random_uniform/max:output:0#dropout/random_uniform/min:output:0*
T0*
_output_shapes
: �
dropout/random_uniform/mulMul-dropout/random_uniform/RandomUniform:output:0dropout/random_uniform/sub:z:0*
T0*(
_output_shapes
:�����������
dropout/random_uniformAdddropout/random_uniform/mul:z:0#dropout/random_uniform/min:output:0*(
_output_shapes
:����������*
T0R
dropout/sub/xConst*
dtype0*
_output_shapes
: *
valueB
 *  �?b
dropout/subSubdropout/sub/x:output:0dropout/rate:output:0*
T0*
_output_shapes
: V
dropout/truediv/xConst*
valueB
 *  �?*
dtype0*
_output_shapes
: h
dropout/truedivRealDivdropout/truediv/x:output:0dropout/sub:z:0*
_output_shapes
: *
T0�
dropout/GreaterEqualGreaterEqualdropout/random_uniform:z:0dropout/rate:output:0*
T0*(
_output_shapes
:����������b
dropout/mulMulinputsdropout/truediv:z:0*
T0*(
_output_shapes
:����������p
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*(
_output_shapes
:����������*

SrcT0
j
dropout/mul_1Muldropout/mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:����������Z
IdentityIdentitydropout/mul_1:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*'
_input_shapes
:����������:& "
 
_user_specified_nameinputs
�
�
"__inference_be1_layer_call_fn_4475

inputs/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma+
'statefulpartitionedcall_be1_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma'statefulpartitionedcall_be1_moving_mean*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3454*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3453�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�
�
%__inference_z_mean_layer_call_fn_4618

inputs)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3920*I
fDRB
@__inference_z_mean_layer_call_and_return_conditional_losses_3914*
Tout
2**
config_proto

CPU

GPU 2J 8�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*/
_input_shapes
:����������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�
�
<__inference_bz_layer_call_and_return_conditional_losses_4674

inputs/
+batchnorm_readvariableop_bz_moving_variance)
%batchnorm_mul_readvariableop_bz_gamma-
)batchnorm_readvariableop_1_bz_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 Z N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp+batchnorm_readvariableop_bz_moving_variance*
dtype0*
_output_shapes
:T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: w
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:~
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:t
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:����������
batchnorm/ReadVariableOp_1ReadVariableOp)batchnorm_readvariableop_1_bz_moving_mean*
dtype0*
_output_shapes
:]
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes
:a
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes
:t
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_124
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�
�
"__inference_be2_layer_call_fn_4592

inputs+
'statefulpartitionedcall_be2_moving_mean/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs'statefulpartitionedcall_be2_moving_mean+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3559*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3558�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*3
_input_shapes"
 :����������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�
o
E__inference_concatenate_layer_call_and_return_conditional_losses_3748

inputs
inputs_1
identityM
concat/axisConst*
dtype0*
_output_shapes
: *
value	B :v
concatConcatV2inputsinputs_1concat/axis:output:0*
T0*
N*(
_output_shapes
:����������X
IdentityIdentityconcat:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*:
_input_shapes)
':����������:���������:& "
 
_user_specified_nameinputs:&"
 
_user_specified_nameinputs
�
`
A__inference_dropout_layer_call_and_return_conditional_losses_3831

inputs
identity�Q
dropout/rateConst*
dtype0*
_output_shapes
: *
valueB
 *���=C
dropout/ShapeShapeinputs*
_output_shapes
:*
T0_
dropout/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *    _
dropout/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *  �?�
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*
dtype0*(
_output_shapes
:�����������
dropout/random_uniform/subSub#dropout/random_uniform/max:output:0#dropout/random_uniform/min:output:0*
T0*
_output_shapes
: �
dropout/random_uniform/mulMul-dropout/random_uniform/RandomUniform:output:0dropout/random_uniform/sub:z:0*
T0*(
_output_shapes
:�����������
dropout/random_uniformAdddropout/random_uniform/mul:z:0#dropout/random_uniform/min:output:0*
T0*(
_output_shapes
:����������R
dropout/sub/xConst*
valueB
 *  �?*
dtype0*
_output_shapes
: b
dropout/subSubdropout/sub/x:output:0dropout/rate:output:0*
T0*
_output_shapes
: V
dropout/truediv/xConst*
valueB
 *  �?*
dtype0*
_output_shapes
: h
dropout/truedivRealDivdropout/truediv/x:output:0dropout/sub:z:0*
T0*
_output_shapes
: �
dropout/GreaterEqualGreaterEqualdropout/random_uniform:z:0dropout/rate:output:0*
T0*(
_output_shapes
:����������b
dropout/mulMulinputsdropout/truediv:z:0*
T0*(
_output_shapes
:����������p
dropout/CastCastdropout/GreaterEqual:z:0*

SrcT0
*

DstT0*(
_output_shapes
:����������j
dropout/mul_1Muldropout/mul:z:0dropout/Cast:y:0*(
_output_shapes
:����������*
T0Z
IdentityIdentitydropout/mul_1:z:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*'
_input_shapes
:����������:& "
 
_user_specified_nameinputs
�
�
!__inference_e2_layer_call_fn_4528

inputs%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias*+
_gradient_op_typePartitionedCall-3872*E
f@R>
<__inference_e2_layer_call_and_return_conditional_losses_3866*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*/
_input_shapes
:����������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
�	
�
<__inference_e2_layer_call_and_return_conditional_losses_4521

inputs#
matmul_readvariableop_e2_kernel"
biasadd_readvariableop_e2_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpw
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_e2_kernel*
dtype0* 
_output_shapes
:
��j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*(
_output_shapes
:����������*
T0r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_e2_bias*
dtype0*
_output_shapes	
:�w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������n
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:�����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�
�
"__inference_be1_layer_call_fn_4467

inputs+
'statefulpartitionedcall_be1_moving_mean/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs'statefulpartitionedcall_be1_moving_mean+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3421*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3420*
Tout
2**
config_proto

CPU

GPU 2J 8�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�
�
!__inference_bz_layer_call_fn_4690

inputs.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma*
&statefulpartitionedcall_bz_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma&statefulpartitionedcall_bz_moving_mean*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3730*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3729*
Tout
2**
config_proto

CPU

GPU 2J 8�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : : 
�Y
�
__inference__wrapped_model_3319
x_input
b,
(model_e1_matmul_readvariableop_e1_kernel+
'model_e1_biasadd_readvariableop_e1_bias:
6model_be1_batchnorm_readvariableop_be1_moving_variance4
0model_be1_batchnorm_mul_readvariableop_be1_gamma8
4model_be1_batchnorm_readvariableop_1_be1_moving_mean,
(model_e2_matmul_readvariableop_e2_kernel+
'model_e2_biasadd_readvariableop_e2_bias:
6model_be2_batchnorm_readvariableop_be2_moving_variance4
0model_be2_batchnorm_mul_readvariableop_be2_gamma8
4model_be2_batchnorm_readvariableop_1_be2_moving_mean4
0model_z_mean_matmul_readvariableop_z_mean_kernel3
/model_z_mean_biasadd_readvariableop_z_mean_bias8
4model_bz_batchnorm_readvariableop_bz_moving_variance2
.model_bz_batchnorm_mul_readvariableop_bz_gamma6
2model_bz_batchnorm_readvariableop_1_bz_moving_mean
identity��"model/be1/batchnorm/ReadVariableOp�$model/be1/batchnorm/ReadVariableOp_1�&model/be1/batchnorm/mul/ReadVariableOp�"model/be2/batchnorm/ReadVariableOp�$model/be2/batchnorm/ReadVariableOp_1�&model/be2/batchnorm/mul/ReadVariableOp�!model/bz/batchnorm/ReadVariableOp�#model/bz/batchnorm/ReadVariableOp_1�%model/bz/batchnorm/mul/ReadVariableOp�model/e1/BiasAdd/ReadVariableOp�model/e1/MatMul/ReadVariableOp�model/e2/BiasAdd/ReadVariableOp�model/e2/MatMul/ReadVariableOp�#model/z_mean/BiasAdd/ReadVariableOp�"model/z_mean/MatMul/ReadVariableOp_
model/concatenate/concat/axisConst*
value	B :*
dtype0*
_output_shapes
: �
model/concatenate/concatConcatV2x_inputb&model/concatenate/concat/axis:output:0*
T0*
N*(
_output_shapes
:�����������
model/e1/MatMul/ReadVariableOpReadVariableOp(model_e1_matmul_readvariableop_e1_kernel*
dtype0* 
_output_shapes
:
���
model/e1/MatMulMatMul!model/concatenate/concat:output:0&model/e1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
model/e1/BiasAdd/ReadVariableOpReadVariableOp'model_e1_biasadd_readvariableop_e1_bias*
dtype0*
_output_shapes	
:��
model/e1/BiasAddBiasAddmodel/e1/MatMul:product:0'model/e1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
model/e1/leaky_re_lu/LeakyRelu	LeakyRelumodel/e1/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������X
model/be1/LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: X
model/be1/LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z|
model/be1/LogicalAnd
LogicalAndmodel/be1/LogicalAnd/x:output:0model/be1/LogicalAnd/y:output:0*
_output_shapes
: �
"model/be1/batchnorm/ReadVariableOpReadVariableOp6model_be1_batchnorm_readvariableop_be1_moving_variance*
dtype0*
_output_shapes	
:�^
model/be1/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: �
model/be1/batchnorm/addAddV2*model/be1/batchnorm/ReadVariableOp:value:0"model/be1/batchnorm/add/y:output:0*
T0*
_output_shapes	
:�e
model/be1/batchnorm/RsqrtRsqrtmodel/be1/batchnorm/add:z:0*
T0*
_output_shapes	
:��
&model/be1/batchnorm/mul/ReadVariableOpReadVariableOp0model_be1_batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:��
model/be1/batchnorm/mulMulmodel/be1/batchnorm/Rsqrt:y:0.model/be1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:��
model/be1/batchnorm/mul_1Mul,model/e1/leaky_re_lu/LeakyRelu:activations:0model/be1/batchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
$model/be1/batchnorm/ReadVariableOp_1ReadVariableOp4model_be1_batchnorm_readvariableop_1_be1_moving_mean*
dtype0*
_output_shapes	
:�r
model/be1/batchnorm/NegNeg,model/be1/batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes	
:��
model/be1/batchnorm/mul_2Mulmodel/be1/batchnorm/Neg:y:0model/be1/batchnorm/mul:z:0*
T0*
_output_shapes	
:��
model/be1/batchnorm/add_1AddV2model/be1/batchnorm/mul_1:z:0model/be1/batchnorm/mul_2:z:0*
T0*(
_output_shapes
:����������t
model/dropout/IdentityIdentitymodel/be1/batchnorm/add_1:z:0*
T0*(
_output_shapes
:�����������
model/e2/MatMul/ReadVariableOpReadVariableOp(model_e2_matmul_readvariableop_e2_kernel*
dtype0* 
_output_shapes
:
���
model/e2/MatMulMatMulmodel/dropout/Identity:output:0&model/e2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
model/e2/BiasAdd/ReadVariableOpReadVariableOp'model_e2_biasadd_readvariableop_e2_bias*
dtype0*
_output_shapes	
:��
model/e2/BiasAddBiasAddmodel/e2/MatMul:product:0'model/e2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:�����������
model/e2/leaky_re_lu/LeakyRelu	LeakyRelumodel/e2/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������X
model/be2/LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 Z X
model/be2/LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: |
model/be2/LogicalAnd
LogicalAndmodel/be2/LogicalAnd/x:output:0model/be2/LogicalAnd/y:output:0*
_output_shapes
: �
"model/be2/batchnorm/ReadVariableOpReadVariableOp6model_be2_batchnorm_readvariableop_be2_moving_variance*
dtype0*
_output_shapes	
:�^
model/be2/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: �
model/be2/batchnorm/addAddV2*model/be2/batchnorm/ReadVariableOp:value:0"model/be2/batchnorm/add/y:output:0*
_output_shapes	
:�*
T0e
model/be2/batchnorm/RsqrtRsqrtmodel/be2/batchnorm/add:z:0*
_output_shapes	
:�*
T0�
&model/be2/batchnorm/mul/ReadVariableOpReadVariableOp0model_be2_batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:��
model/be2/batchnorm/mulMulmodel/be2/batchnorm/Rsqrt:y:0.model/be2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:��
model/be2/batchnorm/mul_1Mul,model/e2/leaky_re_lu/LeakyRelu:activations:0model/be2/batchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
$model/be2/batchnorm/ReadVariableOp_1ReadVariableOp4model_be2_batchnorm_readvariableop_1_be2_moving_mean*
dtype0*
_output_shapes	
:�r
model/be2/batchnorm/NegNeg,model/be2/batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes	
:��
model/be2/batchnorm/mul_2Mulmodel/be2/batchnorm/Neg:y:0model/be2/batchnorm/mul:z:0*
T0*
_output_shapes	
:��
model/be2/batchnorm/add_1AddV2model/be2/batchnorm/mul_1:z:0model/be2/batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
"model/z_mean/MatMul/ReadVariableOpReadVariableOp0model_z_mean_matmul_readvariableop_z_mean_kernel*
dtype0*
_output_shapes
:	��
model/z_mean/MatMulMatMulmodel/be2/batchnorm/add_1:z:0*model/z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
#model/z_mean/BiasAdd/ReadVariableOpReadVariableOp/model_z_mean_biasadd_readvariableop_z_mean_bias*
dtype0*
_output_shapes
:�
model/z_mean/BiasAddBiasAddmodel/z_mean/MatMul:product:0+model/z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
"model/z_mean/leaky_re_lu/LeakyRelu	LeakyRelumodel/z_mean/BiasAdd:output:0*
alpha%
�#<*'
_output_shapes
:���������W
model/bz/LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 Z W
model/bz/LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: y
model/bz/LogicalAnd
LogicalAndmodel/bz/LogicalAnd/x:output:0model/bz/LogicalAnd/y:output:0*
_output_shapes
: �
!model/bz/batchnorm/ReadVariableOpReadVariableOp4model_bz_batchnorm_readvariableop_bz_moving_variance*
dtype0*
_output_shapes
:]
model/bz/batchnorm/add/yConst*
dtype0*
_output_shapes
: *
valueB
 *o�:�
model/bz/batchnorm/addAddV2)model/bz/batchnorm/ReadVariableOp:value:0!model/bz/batchnorm/add/y:output:0*
T0*
_output_shapes
:b
model/bz/batchnorm/RsqrtRsqrtmodel/bz/batchnorm/add:z:0*
T0*
_output_shapes
:�
%model/bz/batchnorm/mul/ReadVariableOpReadVariableOp.model_bz_batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:�
model/bz/batchnorm/mulMulmodel/bz/batchnorm/Rsqrt:y:0-model/bz/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:�
model/bz/batchnorm/mul_1Mul0model/z_mean/leaky_re_lu/LeakyRelu:activations:0model/bz/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
#model/bz/batchnorm/ReadVariableOp_1ReadVariableOp2model_bz_batchnorm_readvariableop_1_bz_moving_mean*
dtype0*
_output_shapes
:o
model/bz/batchnorm/NegNeg+model/bz/batchnorm/ReadVariableOp_1:value:0*
_output_shapes
:*
T0|
model/bz/batchnorm/mul_2Mulmodel/bz/batchnorm/Neg:y:0model/bz/batchnorm/mul:z:0*
T0*
_output_shapes
:�
model/bz/batchnorm/add_1AddV2model/bz/batchnorm/mul_1:z:0model/bz/batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitymodel/bz/batchnorm/add_1:z:0#^model/be1/batchnorm/ReadVariableOp%^model/be1/batchnorm/ReadVariableOp_1'^model/be1/batchnorm/mul/ReadVariableOp#^model/be2/batchnorm/ReadVariableOp%^model/be2/batchnorm/ReadVariableOp_1'^model/be2/batchnorm/mul/ReadVariableOp"^model/bz/batchnorm/ReadVariableOp$^model/bz/batchnorm/ReadVariableOp_1&^model/bz/batchnorm/mul/ReadVariableOp ^model/e1/BiasAdd/ReadVariableOp^model/e1/MatMul/ReadVariableOp ^model/e2/BiasAdd/ReadVariableOp^model/e2/MatMul/ReadVariableOp$^model/z_mean/BiasAdd/ReadVariableOp#^model/z_mean/MatMul/ReadVariableOp*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2@
model/e2/MatMul/ReadVariableOpmodel/e2/MatMul/ReadVariableOp2B
model/e2/BiasAdd/ReadVariableOpmodel/e2/BiasAdd/ReadVariableOp2H
"model/be1/batchnorm/ReadVariableOp"model/be1/batchnorm/ReadVariableOp2J
#model/z_mean/BiasAdd/ReadVariableOp#model/z_mean/BiasAdd/ReadVariableOp2L
$model/be2/batchnorm/ReadVariableOp_1$model/be2/batchnorm/ReadVariableOp_12L
$model/be1/batchnorm/ReadVariableOp_1$model/be1/batchnorm/ReadVariableOp_12H
"model/z_mean/MatMul/ReadVariableOp"model/z_mean/MatMul/ReadVariableOp2H
"model/be2/batchnorm/ReadVariableOp"model/be2/batchnorm/ReadVariableOp2P
&model/be2/batchnorm/mul/ReadVariableOp&model/be2/batchnorm/mul/ReadVariableOp2B
model/e1/BiasAdd/ReadVariableOpmodel/e1/BiasAdd/ReadVariableOp2N
%model/bz/batchnorm/mul/ReadVariableOp%model/bz/batchnorm/mul/ReadVariableOp2@
model/e1/MatMul/ReadVariableOpmodel/e1/MatMul/ReadVariableOp2F
!model/bz/batchnorm/ReadVariableOp!model/bz/batchnorm/ReadVariableOp2P
&model/be1/batchnorm/mul/ReadVariableOp&model/be1/batchnorm/mul/ReadVariableOp2J
#model/bz/batchnorm/ReadVariableOp_1#model/bz/batchnorm/ReadVariableOp_1:' #
!
_user_specified_name	x_input:!

_user_specified_nameB: : : : : : : :	 :
 : : : : : : 
�
_
A__inference_dropout_layer_call_and_return_conditional_losses_4500

inputs

identity_1O
IdentityIdentityinputs*(
_output_shapes
:����������*
T0\

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:����������"!

identity_1Identity_1:output:0*'
_input_shapes
:����������:& "
 
_user_specified_nameinputs
�
�
<__inference_bz_layer_call_and_return_conditional_losses_3729

inputs/
+batchnorm_readvariableop_bz_moving_variance)
%batchnorm_mul_readvariableop_bz_gamma-
)batchnorm_readvariableop_1_bz_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp+batchnorm_readvariableop_bz_moving_variance*
dtype0*
_output_shapes
:T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: w
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:P
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:~
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:t
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:c
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*'
_output_shapes
:����������
batchnorm/ReadVariableOp_1ReadVariableOp)batchnorm_readvariableop_1_bz_moving_mean*
dtype0*
_output_shapes
:]
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
_output_shapes
:*
T0a
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes
:t
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*'
_output_shapes
:���������*
T0�
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*2
_input_shapes!
:���������:::28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_124
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp: : :& "
 
_user_specified_nameinputs: 
�	
�
@__inference_z_mean_layer_call_and_return_conditional_losses_3914

inputs'
#matmul_readvariableop_z_mean_kernel&
"biasadd_readvariableop_z_mean_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpz
MatMul/ReadVariableOpReadVariableOp#matmul_readvariableop_z_mean_kernel*
dtype0*
_output_shapes
:	�i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������u
BiasAdd/ReadVariableOpReadVariableOp"biasadd_readvariableop_z_mean_bias*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0m
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*'
_output_shapes
:����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�P
�

?__inference_model_layer_call_and_return_conditional_losses_4330
inputs_0
inputs_1&
"e1_matmul_readvariableop_e1_kernel%
!e1_biasadd_readvariableop_e1_bias4
0be1_batchnorm_readvariableop_be1_moving_variance.
*be1_batchnorm_mul_readvariableop_be1_gamma2
.be1_batchnorm_readvariableop_1_be1_moving_mean&
"e2_matmul_readvariableop_e2_kernel%
!e2_biasadd_readvariableop_e2_bias4
0be2_batchnorm_readvariableop_be2_moving_variance.
*be2_batchnorm_mul_readvariableop_be2_gamma2
.be2_batchnorm_readvariableop_1_be2_moving_mean.
*z_mean_matmul_readvariableop_z_mean_kernel-
)z_mean_biasadd_readvariableop_z_mean_bias2
.bz_batchnorm_readvariableop_bz_moving_variance,
(bz_batchnorm_mul_readvariableop_bz_gamma0
,bz_batchnorm_readvariableop_1_bz_moving_mean
identity��be1/batchnorm/ReadVariableOp�be1/batchnorm/ReadVariableOp_1� be1/batchnorm/mul/ReadVariableOp�be2/batchnorm/ReadVariableOp�be2/batchnorm/ReadVariableOp_1� be2/batchnorm/mul/ReadVariableOp�bz/batchnorm/ReadVariableOp�bz/batchnorm/ReadVariableOp_1�bz/batchnorm/mul/ReadVariableOp�e1/BiasAdd/ReadVariableOp�e1/MatMul/ReadVariableOp�e2/BiasAdd/ReadVariableOp�e2/MatMul/ReadVariableOp�z_mean/BiasAdd/ReadVariableOp�z_mean/MatMul/ReadVariableOpY
concatenate/concat/axisConst*
value	B :*
dtype0*
_output_shapes
: �
concatenate/concatConcatV2inputs_0inputs_1 concatenate/concat/axis:output:0*
N*(
_output_shapes
:����������*
T0}
e1/MatMul/ReadVariableOpReadVariableOp"e1_matmul_readvariableop_e1_kernel*
dtype0* 
_output_shapes
:
���
	e1/MatMulMatMulconcatenate/concat:output:0 e1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������x
e1/BiasAdd/ReadVariableOpReadVariableOp!e1_biasadd_readvariableop_e1_bias*
dtype0*
_output_shapes	
:��

e1/BiasAddBiasAdde1/MatMul:product:0!e1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������t
e1/leaky_re_lu/LeakyRelu	LeakyRelue1/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������R
be1/LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 Z R
be1/LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Zj
be1/LogicalAnd
LogicalAndbe1/LogicalAnd/x:output:0be1/LogicalAnd/y:output:0*
_output_shapes
: �
be1/batchnorm/ReadVariableOpReadVariableOp0be1_batchnorm_readvariableop_be1_moving_variance*
dtype0*
_output_shapes	
:�X
be1/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: �
be1/batchnorm/addAddV2$be1/batchnorm/ReadVariableOp:value:0be1/batchnorm/add/y:output:0*
_output_shapes	
:�*
T0Y
be1/batchnorm/RsqrtRsqrtbe1/batchnorm/add:z:0*
_output_shapes	
:�*
T0�
 be1/batchnorm/mul/ReadVariableOpReadVariableOp*be1_batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:��
be1/batchnorm/mulMulbe1/batchnorm/Rsqrt:y:0(be1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:��
be1/batchnorm/mul_1Mul&e1/leaky_re_lu/LeakyRelu:activations:0be1/batchnorm/mul:z:0*(
_output_shapes
:����������*
T0�
be1/batchnorm/ReadVariableOp_1ReadVariableOp.be1_batchnorm_readvariableop_1_be1_moving_mean*
dtype0*
_output_shapes	
:�f
be1/batchnorm/NegNeg&be1/batchnorm/ReadVariableOp_1:value:0*
_output_shapes	
:�*
T0n
be1/batchnorm/mul_2Mulbe1/batchnorm/Neg:y:0be1/batchnorm/mul:z:0*
T0*
_output_shapes	
:��
be1/batchnorm/add_1AddV2be1/batchnorm/mul_1:z:0be1/batchnorm/mul_2:z:0*(
_output_shapes
:����������*
T0h
dropout/IdentityIdentitybe1/batchnorm/add_1:z:0*(
_output_shapes
:����������*
T0}
e2/MatMul/ReadVariableOpReadVariableOp"e2_matmul_readvariableop_e2_kernel*
dtype0* 
_output_shapes
:
���
	e2/MatMulMatMuldropout/Identity:output:0 e2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������x
e2/BiasAdd/ReadVariableOpReadVariableOp!e2_biasadd_readvariableop_e2_bias*
dtype0*
_output_shapes	
:��

e2/BiasAddBiasAdde2/MatMul:product:0!e2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������t
e2/leaky_re_lu/LeakyRelu	LeakyRelue2/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������R
be2/LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: R
be2/LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Zj
be2/LogicalAnd
LogicalAndbe2/LogicalAnd/x:output:0be2/LogicalAnd/y:output:0*
_output_shapes
: �
be2/batchnorm/ReadVariableOpReadVariableOp0be2_batchnorm_readvariableop_be2_moving_variance*
dtype0*
_output_shapes	
:�X
be2/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: �
be2/batchnorm/addAddV2$be2/batchnorm/ReadVariableOp:value:0be2/batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Y
be2/batchnorm/RsqrtRsqrtbe2/batchnorm/add:z:0*
T0*
_output_shapes	
:��
 be2/batchnorm/mul/ReadVariableOpReadVariableOp*be2_batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:��
be2/batchnorm/mulMulbe2/batchnorm/Rsqrt:y:0(be2/batchnorm/mul/ReadVariableOp:value:0*
_output_shapes	
:�*
T0�
be2/batchnorm/mul_1Mul&e2/leaky_re_lu/LeakyRelu:activations:0be2/batchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
be2/batchnorm/ReadVariableOp_1ReadVariableOp.be2_batchnorm_readvariableop_1_be2_moving_mean*
dtype0*
_output_shapes	
:�f
be2/batchnorm/NegNeg&be2/batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes	
:�n
be2/batchnorm/mul_2Mulbe2/batchnorm/Neg:y:0be2/batchnorm/mul:z:0*
T0*
_output_shapes	
:��
be2/batchnorm/add_1AddV2be2/batchnorm/mul_1:z:0be2/batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
z_mean/MatMul/ReadVariableOpReadVariableOp*z_mean_matmul_readvariableop_z_mean_kernel*
dtype0*
_output_shapes
:	��
z_mean/MatMulMatMulbe2/batchnorm/add_1:z:0$z_mean/MatMul/ReadVariableOp:value:0*'
_output_shapes
:���������*
T0�
z_mean/BiasAdd/ReadVariableOpReadVariableOp)z_mean_biasadd_readvariableop_z_mean_bias*
dtype0*
_output_shapes
:�
z_mean/BiasAddBiasAddz_mean/MatMul:product:0%z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������{
z_mean/leaky_re_lu/LeakyRelu	LeakyReluz_mean/BiasAdd:output:0*
alpha%
�#<*'
_output_shapes
:���������Q
bz/LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: Q
bz/LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Zg
bz/LogicalAnd
LogicalAndbz/LogicalAnd/x:output:0bz/LogicalAnd/y:output:0*
_output_shapes
: �
bz/batchnorm/ReadVariableOpReadVariableOp.bz_batchnorm_readvariableop_bz_moving_variance*
dtype0*
_output_shapes
:W
bz/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: �
bz/batchnorm/addAddV2#bz/batchnorm/ReadVariableOp:value:0bz/batchnorm/add/y:output:0*
T0*
_output_shapes
:V
bz/batchnorm/RsqrtRsqrtbz/batchnorm/add:z:0*
_output_shapes
:*
T0�
bz/batchnorm/mul/ReadVariableOpReadVariableOp(bz_batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:}
bz/batchnorm/mulMulbz/batchnorm/Rsqrt:y:0'bz/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:�
bz/batchnorm/mul_1Mul*z_mean/leaky_re_lu/LeakyRelu:activations:0bz/batchnorm/mul:z:0*
T0*'
_output_shapes
:����������
bz/batchnorm/ReadVariableOp_1ReadVariableOp,bz_batchnorm_readvariableop_1_bz_moving_mean*
dtype0*
_output_shapes
:c
bz/batchnorm/NegNeg%bz/batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes
:j
bz/batchnorm/mul_2Mulbz/batchnorm/Neg:y:0bz/batchnorm/mul:z:0*
T0*
_output_shapes
:}
bz/batchnorm/add_1AddV2bz/batchnorm/mul_1:z:0bz/batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitybz/batchnorm/add_1:z:0^be1/batchnorm/ReadVariableOp^be1/batchnorm/ReadVariableOp_1!^be1/batchnorm/mul/ReadVariableOp^be2/batchnorm/ReadVariableOp^be2/batchnorm/ReadVariableOp_1!^be2/batchnorm/mul/ReadVariableOp^bz/batchnorm/ReadVariableOp^bz/batchnorm/ReadVariableOp_1 ^bz/batchnorm/mul/ReadVariableOp^e1/BiasAdd/ReadVariableOp^e1/MatMul/ReadVariableOp^e2/BiasAdd/ReadVariableOp^e2/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::24
e1/MatMul/ReadVariableOpe1/MatMul/ReadVariableOp2D
 be1/batchnorm/mul/ReadVariableOp be1/batchnorm/mul/ReadVariableOp24
e2/MatMul/ReadVariableOpe2/MatMul/ReadVariableOp2<
be1/batchnorm/ReadVariableOpbe1/batchnorm/ReadVariableOp26
e1/BiasAdd/ReadVariableOpe1/BiasAdd/ReadVariableOp2@
be2/batchnorm/ReadVariableOp_1be2/batchnorm/ReadVariableOp_12B
bz/batchnorm/mul/ReadVariableOpbz/batchnorm/mul/ReadVariableOp2@
be1/batchnorm/ReadVariableOp_1be1/batchnorm/ReadVariableOp_12<
z_mean/MatMul/ReadVariableOpz_mean/MatMul/ReadVariableOp2<
be2/batchnorm/ReadVariableOpbe2/batchnorm/ReadVariableOp2:
bz/batchnorm/ReadVariableOpbz/batchnorm/ReadVariableOp2>
bz/batchnorm/ReadVariableOp_1bz/batchnorm/ReadVariableOp_126
e2/BiasAdd/ReadVariableOpe2/BiasAdd/ReadVariableOp2D
 be2/batchnorm/mul/ReadVariableOp be2/batchnorm/mul/ReadVariableOp2>
z_mean/BiasAdd/ReadVariableOpz_mean/BiasAdd/ReadVariableOp:( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1: : : : : : : :	 :
 : : : : : : 
�)
�
?__inference_model_layer_call_and_return_conditional_losses_3952
x_input
b(
$e1_statefulpartitionedcall_e1_kernel&
"e1_statefulpartitionedcall_e1_bias/
+be1_statefulpartitionedcall_be1_moving_mean3
/be1_statefulpartitionedcall_be1_moving_variance)
%be1_statefulpartitionedcall_be1_gamma(
$e2_statefulpartitionedcall_e2_kernel&
"e2_statefulpartitionedcall_e2_bias/
+be2_statefulpartitionedcall_be2_moving_mean3
/be2_statefulpartitionedcall_be2_moving_variance)
%be2_statefulpartitionedcall_be2_gamma0
,z_mean_statefulpartitionedcall_z_mean_kernel.
*z_mean_statefulpartitionedcall_z_mean_bias-
)bz_statefulpartitionedcall_bz_moving_mean1
-bz_statefulpartitionedcall_bz_moving_variance'
#bz_statefulpartitionedcall_bz_gamma
identity��be1/StatefulPartitionedCall�be2/StatefulPartitionedCall�bz/StatefulPartitionedCall�dropout/StatefulPartitionedCall�e1/StatefulPartitionedCall�e2/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
concatenate/PartitionedCallPartitionedCallx_inputb*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3755*N
fIRG
E__inference_concatenate_layer_call_and_return_conditional_losses_3748�
e1/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0$e1_statefulpartitionedcall_e1_kernel"e1_statefulpartitionedcall_e1_bias*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3779*E
f@R>
<__inference_e1_layer_call_and_return_conditional_losses_3773*
Tout
2**
config_proto

CPU

GPU 2J 8�
be1/StatefulPartitionedCallStatefulPartitionedCall#e1/StatefulPartitionedCall:output:0+be1_statefulpartitionedcall_be1_moving_mean/be1_statefulpartitionedcall_be1_moving_variance%be1_statefulpartitionedcall_be1_gamma*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3421*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3420�
dropout/StatefulPartitionedCallStatefulPartitionedCall$be1/StatefulPartitionedCall:output:0**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3842*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3831*
Tout
2�
e2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0$e2_statefulpartitionedcall_e2_kernel"e2_statefulpartitionedcall_e2_bias*+
_gradient_op_typePartitionedCall-3872*E
f@R>
<__inference_e2_layer_call_and_return_conditional_losses_3866*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
be2/StatefulPartitionedCallStatefulPartitionedCall#e2/StatefulPartitionedCall:output:0+be2_statefulpartitionedcall_be2_moving_mean/be2_statefulpartitionedcall_be2_moving_variance%be2_statefulpartitionedcall_be2_gamma*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3558*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3559�
z_mean/StatefulPartitionedCallStatefulPartitionedCall$be2/StatefulPartitionedCall:output:0,z_mean_statefulpartitionedcall_z_mean_kernel*z_mean_statefulpartitionedcall_z_mean_bias**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3920*I
fDRB
@__inference_z_mean_layer_call_and_return_conditional_losses_3914*
Tout
2�
bz/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0)bz_statefulpartitionedcall_bz_moving_mean-bz_statefulpartitionedcall_bz_moving_variance#bz_statefulpartitionedcall_bz_gamma**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2*+
_gradient_op_typePartitionedCall-3697*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3696*
Tout
2�
IdentityIdentity#bz/StatefulPartitionedCall:output:0^be1/StatefulPartitionedCall^be2/StatefulPartitionedCall^bz/StatefulPartitionedCall ^dropout/StatefulPartitionedCall^e1/StatefulPartitionedCall^e2/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall28
bz/StatefulPartitionedCallbz/StatefulPartitionedCall28
e1/StatefulPartitionedCalle1/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall28
e2/StatefulPartitionedCalle2/StatefulPartitionedCall2:
be1/StatefulPartitionedCallbe1/StatefulPartitionedCall2:
be2/StatefulPartitionedCallbe2/StatefulPartitionedCall:' #
!
_user_specified_name	x_input:!

_user_specified_nameB: : : : : : : :	 :
 : : : : : : 
�
V
*__inference_concatenate_layer_call_fn_4385
inputs_0
inputs_1
identity�
PartitionedCallPartitionedCallinputs_0inputs_1**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3755*N
fIRG
E__inference_concatenate_layer_call_and_return_conditional_losses_3748*
Tout
2a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*:
_input_shapes)
':����������:���������:( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1
�
�
=__inference_be1_layer_call_and_return_conditional_losses_3453

inputs0
,batchnorm_readvariableop_be1_moving_variance*
&batchnorm_mul_readvariableop_be1_gamma.
*batchnorm_readvariableop_1_be1_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp,batchnorm_readvariableop_be1_moving_variance*
dtype0*
_output_shapes	
:�T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: x
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
_output_shapes	
:�*
T0d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
batchnorm/ReadVariableOp_1ReadVariableOp*batchnorm_readvariableop_1_be1_moving_mean*
dtype0*
_output_shapes	
:�^
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
_output_shapes	
:�*
T0b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*(
_output_shapes
:����������*
T0�
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_124
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�
_
A__inference_dropout_layer_call_and_return_conditional_losses_3838

inputs

identity_1O
IdentityIdentityinputs*
T0*(
_output_shapes
:����������\

Identity_1IdentityIdentity:output:0*(
_output_shapes
:����������*
T0"!

identity_1Identity_1:output:0*'
_input_shapes
:����������:& "
 
_user_specified_nameinputs
�(
�
=__inference_be1_layer_call_and_return_conditional_losses_3420

inputs#
assignmovingavg_be1_moving_mean)
%assignmovingavg_1_be1_moving_variance*
&batchnorm_mul_readvariableop_be1_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(e
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes
:	��
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*(
_output_shapes
:����������l
"moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(n
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes	
:�*
squeeze_dims
 t
moments/Squeeze_1Squeezemoments/variance:output:0*
_output_shapes	
:�*
squeeze_dims
 *
T0�
AssignMovingAvg/decayConst*
valueB
 *
�#<*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
: {
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_be1_moving_mean*
dtype0*
_output_shapes	
:��
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
_output_shapes	
:��
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
_output_shapes	
:�*
T0*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean�
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_be1_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*2
_class(
&$loc:@AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
dtype0*
_output_shapes
: *
valueB
 *
�#<*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp%assignmovingavg_1_be1_moving_variance*
dtype0*
_output_shapes	
:��
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance*
_output_shapes	
:��
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance*
_output_shapes	
:��
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp%assignmovingavg_1_be1_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*8
_class.
,*loc:@AssignMovingAvg_1/be1/moving_variance*
dtype0*
_output_shapes
 T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: r
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
_output_shapes	
:�*
T0d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:����������T
batchnorm/NegNegmoments/Squeeze:output:0*
_output_shapes	
:�*
T0b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp: : :& "
 
_user_specified_nameinputs: 
�
�
=__inference_be2_layer_call_and_return_conditional_losses_3591

inputs0
,batchnorm_readvariableop_be2_moving_variance*
&batchnorm_mul_readvariableop_be2_gamma.
*batchnorm_readvariableop_1_be2_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
dtype0
*
_output_shapes
: *
value	B
 Z^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp,batchnorm_readvariableop_be2_moving_variance*
dtype0*
_output_shapes	
:�T
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: x
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
_output_shapes	
:�*
T0�
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:�d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
batchnorm/ReadVariableOp_1ReadVariableOp*batchnorm_readvariableop_1_be2_moving_mean*
dtype0*
_output_shapes	
:�^
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes	
:�b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*3
_input_shapes"
 :����������:::28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_124
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp: : :& "
 
_user_specified_nameinputs: 
�	
�
@__inference_z_mean_layer_call_and_return_conditional_losses_4611

inputs'
#matmul_readvariableop_z_mean_kernel&
"biasadd_readvariableop_z_mean_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpz
MatMul/ReadVariableOpReadVariableOp#matmul_readvariableop_z_mean_kernel*
dtype0*
_output_shapes
:	�i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������u
BiasAdd/ReadVariableOpReadVariableOp"biasadd_readvariableop_z_mean_bias*
dtype0*
_output_shapes
:v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������m
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*'
_output_shapes
:����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�(
�
=__inference_be2_layer_call_and_return_conditional_losses_3558

inputs#
assignmovingavg_be2_moving_mean)
%assignmovingavg_1_be2_moving_variance*
&batchnorm_mul_readvariableop_be2_gamma
identity��#AssignMovingAvg/AssignSubVariableOp�AssignMovingAvg/ReadVariableOp�%AssignMovingAvg_1/AssignSubVariableOp� AssignMovingAvg_1/ReadVariableOp�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: h
moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(e
moments/StopGradientStopGradientmoments/mean:output:0*
_output_shapes
:	�*
T0�
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*(
_output_shapes
:����������l
"moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(n
moments/SqueezeSqueezemoments/mean:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:�t
moments/Squeeze_1Squeezemoments/variance:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:��
AssignMovingAvg/decayConst*
valueB
 *
�#<*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
: {
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_be2_moving_mean*
dtype0*
_output_shapes	
:��
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
_output_shapes	
:��
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
_output_shapes	
:�*
T0*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean�
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_be2_moving_meanAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*2
_class(
&$loc:@AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
 �
AssignMovingAvg_1/decayConst*
dtype0*
_output_shapes
: *
valueB
 *
�#<*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance�
 AssignMovingAvg_1/ReadVariableOpReadVariableOp%assignmovingavg_1_be2_moving_variance*
dtype0*
_output_shapes	
:��
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
_output_shapes	
:�*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance�
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
_output_shapes	
:�*
T0*8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_variance�
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp%assignmovingavg_1_be2_moving_varianceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
dtype0*
_output_shapes
 *8
_class.
,*loc:@AssignMovingAvg_1/be2/moving_varianceT
batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: r
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
_output_shapes	
:�*
T0d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*(
_output_shapes
:����������*
T0T
batchnorm/NegNegmoments/Squeeze:output:0*
_output_shapes	
:�*
T0b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp^AssignMovingAvg/ReadVariableOp&^AssignMovingAvg_1/AssignSubVariableOp!^AssignMovingAvg_1/ReadVariableOp^batchnorm/mul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*3
_input_shapes"
 :����������:::2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�)
�
?__inference_model_layer_call_and_return_conditional_losses_4015

inputs
inputs_1(
$e1_statefulpartitionedcall_e1_kernel&
"e1_statefulpartitionedcall_e1_bias/
+be1_statefulpartitionedcall_be1_moving_mean3
/be1_statefulpartitionedcall_be1_moving_variance)
%be1_statefulpartitionedcall_be1_gamma(
$e2_statefulpartitionedcall_e2_kernel&
"e2_statefulpartitionedcall_e2_bias/
+be2_statefulpartitionedcall_be2_moving_mean3
/be2_statefulpartitionedcall_be2_moving_variance)
%be2_statefulpartitionedcall_be2_gamma0
,z_mean_statefulpartitionedcall_z_mean_kernel.
*z_mean_statefulpartitionedcall_z_mean_bias-
)bz_statefulpartitionedcall_bz_moving_mean1
-bz_statefulpartitionedcall_bz_moving_variance'
#bz_statefulpartitionedcall_bz_gamma
identity��be1/StatefulPartitionedCall�be2/StatefulPartitionedCall�bz/StatefulPartitionedCall�dropout/StatefulPartitionedCall�e1/StatefulPartitionedCall�e2/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
concatenate/PartitionedCallPartitionedCallinputsinputs_1*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3755*N
fIRG
E__inference_concatenate_layer_call_and_return_conditional_losses_3748�
e1/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0$e1_statefulpartitionedcall_e1_kernel"e1_statefulpartitionedcall_e1_bias*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3779*E
f@R>
<__inference_e1_layer_call_and_return_conditional_losses_3773*
Tout
2**
config_proto

CPU

GPU 2J 8�
be1/StatefulPartitionedCallStatefulPartitionedCall#e1/StatefulPartitionedCall:output:0+be1_statefulpartitionedcall_be1_moving_mean/be1_statefulpartitionedcall_be1_moving_variance%be1_statefulpartitionedcall_be1_gamma**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3421*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3420*
Tout
2�
dropout/StatefulPartitionedCallStatefulPartitionedCall$be1/StatefulPartitionedCall:output:0*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3842*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3831�
e2/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0$e2_statefulpartitionedcall_e2_kernel"e2_statefulpartitionedcall_e2_bias*+
_gradient_op_typePartitionedCall-3872*E
f@R>
<__inference_e2_layer_call_and_return_conditional_losses_3866*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:�����������
be2/StatefulPartitionedCallStatefulPartitionedCall#e2/StatefulPartitionedCall:output:0+be2_statefulpartitionedcall_be2_moving_mean/be2_statefulpartitionedcall_be2_moving_variance%be2_statefulpartitionedcall_be2_gamma**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3559*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3558*
Tout
2�
z_mean/StatefulPartitionedCallStatefulPartitionedCall$be2/StatefulPartitionedCall:output:0,z_mean_statefulpartitionedcall_z_mean_kernel*z_mean_statefulpartitionedcall_z_mean_bias*+
_gradient_op_typePartitionedCall-3920*I
fDRB
@__inference_z_mean_layer_call_and_return_conditional_losses_3914*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:����������
bz/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0)bz_statefulpartitionedcall_bz_moving_mean-bz_statefulpartitionedcall_bz_moving_variance#bz_statefulpartitionedcall_bz_gamma**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3697*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3696*
Tout
2�
IdentityIdentity#bz/StatefulPartitionedCall:output:0^be1/StatefulPartitionedCall^be2/StatefulPartitionedCall^bz/StatefulPartitionedCall ^dropout/StatefulPartitionedCall^e1/StatefulPartitionedCall^e2/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall28
e1/StatefulPartitionedCalle1/StatefulPartitionedCall28
e2/StatefulPartitionedCalle2/StatefulPartitionedCall2:
be1/StatefulPartitionedCallbe1/StatefulPartitionedCall2:
be2/StatefulPartitionedCallbe2/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall28
bz/StatefulPartitionedCallbz/StatefulPartitionedCall:& "
 
_user_specified_nameinputs:&"
 
_user_specified_nameinputs: : : : : : : :	 :
 : : : : : : 
�
�
$__inference_model_layer_call_fn_4087
x_input
b%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma+
'statefulpartitionedcall_be1_moving_mean%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma+
'statefulpartitionedcall_be2_moving_mean)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma*
&statefulpartitionedcall_bz_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallx_inputb!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma'statefulpartitionedcall_be1_moving_mean!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma'statefulpartitionedcall_be2_moving_mean%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma&statefulpartitionedcall_bz_moving_mean*+
_gradient_op_typePartitionedCall-4069*H
fCRA
?__inference_model_layer_call_and_return_conditional_losses_4068*
Tout
2**
config_proto

CPU

GPU 2J 8*'
_output_shapes
:���������*
Tin
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:
 : : : : : : :' #
!
_user_specified_name	x_input:!

_user_specified_nameB: : : : : : : :	 
�
�
!__inference_e1_layer_call_fn_4403

inputs%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3779*E
f@R>
<__inference_e1_layer_call_and_return_conditional_losses_3773*
Tout
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*/
_input_shapes
:����������::22
StatefulPartitionedCallStatefulPartitionedCall:& "
 
_user_specified_nameinputs: : 
��
�
?__inference_model_layer_call_and_return_conditional_losses_4250
inputs_0
inputs_1&
"e1_matmul_readvariableop_e1_kernel%
!e1_biasadd_readvariableop_e1_bias'
#be1_assignmovingavg_be1_moving_mean-
)be1_assignmovingavg_1_be1_moving_variance.
*be1_batchnorm_mul_readvariableop_be1_gamma&
"e2_matmul_readvariableop_e2_kernel%
!e2_biasadd_readvariableop_e2_bias'
#be2_assignmovingavg_be2_moving_mean-
)be2_assignmovingavg_1_be2_moving_variance.
*be2_batchnorm_mul_readvariableop_be2_gamma.
*z_mean_matmul_readvariableop_z_mean_kernel-
)z_mean_biasadd_readvariableop_z_mean_bias%
!bz_assignmovingavg_bz_moving_mean+
'bz_assignmovingavg_1_bz_moving_variance,
(bz_batchnorm_mul_readvariableop_bz_gamma
identity��'be1/AssignMovingAvg/AssignSubVariableOp�"be1/AssignMovingAvg/ReadVariableOp�)be1/AssignMovingAvg_1/AssignSubVariableOp�$be1/AssignMovingAvg_1/ReadVariableOp� be1/batchnorm/mul/ReadVariableOp�'be2/AssignMovingAvg/AssignSubVariableOp�"be2/AssignMovingAvg/ReadVariableOp�)be2/AssignMovingAvg_1/AssignSubVariableOp�$be2/AssignMovingAvg_1/ReadVariableOp� be2/batchnorm/mul/ReadVariableOp�&bz/AssignMovingAvg/AssignSubVariableOp�!bz/AssignMovingAvg/ReadVariableOp�(bz/AssignMovingAvg_1/AssignSubVariableOp�#bz/AssignMovingAvg_1/ReadVariableOp�bz/batchnorm/mul/ReadVariableOp�e1/BiasAdd/ReadVariableOp�e1/MatMul/ReadVariableOp�e2/BiasAdd/ReadVariableOp�e2/MatMul/ReadVariableOp�z_mean/BiasAdd/ReadVariableOp�z_mean/MatMul/ReadVariableOpY
concatenate/concat/axisConst*
value	B :*
dtype0*
_output_shapes
: �
concatenate/concatConcatV2inputs_0inputs_1 concatenate/concat/axis:output:0*
T0*
N*(
_output_shapes
:����������}
e1/MatMul/ReadVariableOpReadVariableOp"e1_matmul_readvariableop_e1_kernel*
dtype0* 
_output_shapes
:
���
	e1/MatMulMatMulconcatenate/concat:output:0 e1/MatMul/ReadVariableOp:value:0*(
_output_shapes
:����������*
T0x
e1/BiasAdd/ReadVariableOpReadVariableOp!e1_biasadd_readvariableop_e1_bias*
dtype0*
_output_shapes	
:��

e1/BiasAddBiasAdde1/MatMul:product:0!e1/BiasAdd/ReadVariableOp:value:0*(
_output_shapes
:����������*
T0t
e1/leaky_re_lu/LeakyRelu	LeakyRelue1/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������R
be1/LogicalAnd/xConst*
value	B
 Z*
dtype0
*
_output_shapes
: R
be1/LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: j
be1/LogicalAnd
LogicalAndbe1/LogicalAnd/x:output:0be1/LogicalAnd/y:output:0*
_output_shapes
: l
"be1/moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
be1/moments/meanMean&e1/leaky_re_lu/LeakyRelu:activations:0+be1/moments/mean/reduction_indices:output:0*
_output_shapes
:	�*
	keep_dims(*
T0m
be1/moments/StopGradientStopGradientbe1/moments/mean:output:0*
_output_shapes
:	�*
T0�
be1/moments/SquaredDifferenceSquaredDifference&e1/leaky_re_lu/LeakyRelu:activations:0!be1/moments/StopGradient:output:0*(
_output_shapes
:����������*
T0p
&be1/moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
be1/moments/varianceMean!be1/moments/SquaredDifference:z:0/be1/moments/variance/reduction_indices:output:0*
	keep_dims(*
T0*
_output_shapes
:	�v
be1/moments/SqueezeSqueezebe1/moments/mean:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:�|
be1/moments/Squeeze_1Squeezebe1/moments/variance:output:0*
squeeze_dims
 *
T0*
_output_shapes	
:��
be1/AssignMovingAvg/decayConst*
valueB
 *
�#<*6
_class,
*(loc:@be1/AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
: �
"be1/AssignMovingAvg/ReadVariableOpReadVariableOp#be1_assignmovingavg_be1_moving_mean*
dtype0*
_output_shapes	
:��
be1/AssignMovingAvg/subSub*be1/AssignMovingAvg/ReadVariableOp:value:0be1/moments/Squeeze:output:0*
T0*6
_class,
*(loc:@be1/AssignMovingAvg/be1/moving_mean*
_output_shapes	
:��
be1/AssignMovingAvg/mulMulbe1/AssignMovingAvg/sub:z:0"be1/AssignMovingAvg/decay:output:0*
T0*6
_class,
*(loc:@be1/AssignMovingAvg/be1/moving_mean*
_output_shapes	
:��
'be1/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp#be1_assignmovingavg_be1_moving_meanbe1/AssignMovingAvg/mul:z:0#^be1/AssignMovingAvg/ReadVariableOp*6
_class,
*(loc:@be1/AssignMovingAvg/be1/moving_mean*
dtype0*
_output_shapes
 �
be1/AssignMovingAvg_1/decayConst*
dtype0*
_output_shapes
: *
valueB
 *
�#<*<
_class2
0.loc:@be1/AssignMovingAvg_1/be1/moving_variance�
$be1/AssignMovingAvg_1/ReadVariableOpReadVariableOp)be1_assignmovingavg_1_be1_moving_variance*
dtype0*
_output_shapes	
:��
be1/AssignMovingAvg_1/subSub,be1/AssignMovingAvg_1/ReadVariableOp:value:0be1/moments/Squeeze_1:output:0*
T0*<
_class2
0.loc:@be1/AssignMovingAvg_1/be1/moving_variance*
_output_shapes	
:��
be1/AssignMovingAvg_1/mulMulbe1/AssignMovingAvg_1/sub:z:0$be1/AssignMovingAvg_1/decay:output:0*
T0*<
_class2
0.loc:@be1/AssignMovingAvg_1/be1/moving_variance*
_output_shapes	
:��
)be1/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp)be1_assignmovingavg_1_be1_moving_variancebe1/AssignMovingAvg_1/mul:z:0%^be1/AssignMovingAvg_1/ReadVariableOp*<
_class2
0.loc:@be1/AssignMovingAvg_1/be1/moving_variance*
dtype0*
_output_shapes
 X
be1/batchnorm/add/yConst*
dtype0*
_output_shapes
: *
valueB
 *o�:~
be1/batchnorm/addAddV2be1/moments/Squeeze_1:output:0be1/batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Y
be1/batchnorm/RsqrtRsqrtbe1/batchnorm/add:z:0*
T0*
_output_shapes	
:��
 be1/batchnorm/mul/ReadVariableOpReadVariableOp*be1_batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:��
be1/batchnorm/mulMulbe1/batchnorm/Rsqrt:y:0(be1/batchnorm/mul/ReadVariableOp:value:0*
_output_shapes	
:�*
T0�
be1/batchnorm/mul_1Mul&e1/leaky_re_lu/LeakyRelu:activations:0be1/batchnorm/mul:z:0*
T0*(
_output_shapes
:����������\
be1/batchnorm/NegNegbe1/moments/Squeeze:output:0*
T0*
_output_shapes	
:�n
be1/batchnorm/mul_2Mulbe1/batchnorm/Neg:y:0be1/batchnorm/mul:z:0*
T0*
_output_shapes	
:��
be1/batchnorm/add_1AddV2be1/batchnorm/mul_1:z:0be1/batchnorm/mul_2:z:0*
T0*(
_output_shapes
:����������Y
dropout/dropout/rateConst*
valueB
 *���=*
dtype0*
_output_shapes
: \
dropout/dropout/ShapeShapebe1/batchnorm/add_1:z:0*
T0*
_output_shapes
:g
"dropout/dropout/random_uniform/minConst*
valueB
 *    *
dtype0*
_output_shapes
: g
"dropout/dropout/random_uniform/maxConst*
valueB
 *  �?*
dtype0*
_output_shapes
: �
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*
dtype0*(
_output_shapes
:�����������
"dropout/dropout/random_uniform/subSub+dropout/dropout/random_uniform/max:output:0+dropout/dropout/random_uniform/min:output:0*
T0*
_output_shapes
: �
"dropout/dropout/random_uniform/mulMul5dropout/dropout/random_uniform/RandomUniform:output:0&dropout/dropout/random_uniform/sub:z:0*
T0*(
_output_shapes
:�����������
dropout/dropout/random_uniformAdd&dropout/dropout/random_uniform/mul:z:0+dropout/dropout/random_uniform/min:output:0*
T0*(
_output_shapes
:����������Z
dropout/dropout/sub/xConst*
valueB
 *  �?*
dtype0*
_output_shapes
: z
dropout/dropout/subSubdropout/dropout/sub/x:output:0dropout/dropout/rate:output:0*
_output_shapes
: *
T0^
dropout/dropout/truediv/xConst*
dtype0*
_output_shapes
: *
valueB
 *  �?�
dropout/dropout/truedivRealDiv"dropout/dropout/truediv/x:output:0dropout/dropout/sub:z:0*
T0*
_output_shapes
: �
dropout/dropout/GreaterEqualGreaterEqual"dropout/dropout/random_uniform:z:0dropout/dropout/rate:output:0*
T0*(
_output_shapes
:�����������
dropout/dropout/mulMulbe1/batchnorm/add_1:z:0dropout/dropout/truediv:z:0*
T0*(
_output_shapes
:�����������
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

SrcT0
*

DstT0*(
_output_shapes
:�����������
dropout/dropout/mul_1Muldropout/dropout/mul:z:0dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:����������}
e2/MatMul/ReadVariableOpReadVariableOp"e2_matmul_readvariableop_e2_kernel*
dtype0* 
_output_shapes
:
���
	e2/MatMulMatMuldropout/dropout/mul_1:z:0 e2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������x
e2/BiasAdd/ReadVariableOpReadVariableOp!e2_biasadd_readvariableop_e2_bias*
dtype0*
_output_shapes	
:��

e2/BiasAddBiasAdde2/MatMul:product:0!e2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������t
e2/leaky_re_lu/LeakyRelu	LeakyRelue2/BiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:����������R
be2/LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 ZR
be2/LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: j
be2/LogicalAnd
LogicalAndbe2/LogicalAnd/x:output:0be2/LogicalAnd/y:output:0*
_output_shapes
: l
"be2/moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
be2/moments/meanMean&e2/leaky_re_lu/LeakyRelu:activations:0+be2/moments/mean/reduction_indices:output:0*
_output_shapes
:	�*
	keep_dims(*
T0m
be2/moments/StopGradientStopGradientbe2/moments/mean:output:0*
T0*
_output_shapes
:	��
be2/moments/SquaredDifferenceSquaredDifference&e2/leaky_re_lu/LeakyRelu:activations:0!be2/moments/StopGradient:output:0*
T0*(
_output_shapes
:����������p
&be2/moments/variance/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
be2/moments/varianceMean!be2/moments/SquaredDifference:z:0/be2/moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	�*
	keep_dims(v
be2/moments/SqueezeSqueezebe2/moments/mean:output:0*
T0*
_output_shapes	
:�*
squeeze_dims
 |
be2/moments/Squeeze_1Squeezebe2/moments/variance:output:0*
T0*
_output_shapes	
:�*
squeeze_dims
 �
be2/AssignMovingAvg/decayConst*
valueB
 *
�#<*6
_class,
*(loc:@be2/AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
: �
"be2/AssignMovingAvg/ReadVariableOpReadVariableOp#be2_assignmovingavg_be2_moving_mean*
dtype0*
_output_shapes	
:��
be2/AssignMovingAvg/subSub*be2/AssignMovingAvg/ReadVariableOp:value:0be2/moments/Squeeze:output:0*
T0*6
_class,
*(loc:@be2/AssignMovingAvg/be2/moving_mean*
_output_shapes	
:��
be2/AssignMovingAvg/mulMulbe2/AssignMovingAvg/sub:z:0"be2/AssignMovingAvg/decay:output:0*
T0*6
_class,
*(loc:@be2/AssignMovingAvg/be2/moving_mean*
_output_shapes	
:��
'be2/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp#be2_assignmovingavg_be2_moving_meanbe2/AssignMovingAvg/mul:z:0#^be2/AssignMovingAvg/ReadVariableOp*6
_class,
*(loc:@be2/AssignMovingAvg/be2/moving_mean*
dtype0*
_output_shapes
 �
be2/AssignMovingAvg_1/decayConst*
valueB
 *
�#<*<
_class2
0.loc:@be2/AssignMovingAvg_1/be2/moving_variance*
dtype0*
_output_shapes
: �
$be2/AssignMovingAvg_1/ReadVariableOpReadVariableOp)be2_assignmovingavg_1_be2_moving_variance*
dtype0*
_output_shapes	
:��
be2/AssignMovingAvg_1/subSub,be2/AssignMovingAvg_1/ReadVariableOp:value:0be2/moments/Squeeze_1:output:0*
T0*<
_class2
0.loc:@be2/AssignMovingAvg_1/be2/moving_variance*
_output_shapes	
:��
be2/AssignMovingAvg_1/mulMulbe2/AssignMovingAvg_1/sub:z:0$be2/AssignMovingAvg_1/decay:output:0*
T0*<
_class2
0.loc:@be2/AssignMovingAvg_1/be2/moving_variance*
_output_shapes	
:��
)be2/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp)be2_assignmovingavg_1_be2_moving_variancebe2/AssignMovingAvg_1/mul:z:0%^be2/AssignMovingAvg_1/ReadVariableOp*<
_class2
0.loc:@be2/AssignMovingAvg_1/be2/moving_variance*
dtype0*
_output_shapes
 X
be2/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: ~
be2/batchnorm/addAddV2be2/moments/Squeeze_1:output:0be2/batchnorm/add/y:output:0*
T0*
_output_shapes	
:�Y
be2/batchnorm/RsqrtRsqrtbe2/batchnorm/add:z:0*
T0*
_output_shapes	
:��
 be2/batchnorm/mul/ReadVariableOpReadVariableOp*be2_batchnorm_mul_readvariableop_be2_gamma*
dtype0*
_output_shapes	
:��
be2/batchnorm/mulMulbe2/batchnorm/Rsqrt:y:0(be2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:��
be2/batchnorm/mul_1Mul&e2/leaky_re_lu/LeakyRelu:activations:0be2/batchnorm/mul:z:0*
T0*(
_output_shapes
:����������\
be2/batchnorm/NegNegbe2/moments/Squeeze:output:0*
T0*
_output_shapes	
:�n
be2/batchnorm/mul_2Mulbe2/batchnorm/Neg:y:0be2/batchnorm/mul:z:0*
_output_shapes	
:�*
T0�
be2/batchnorm/add_1AddV2be2/batchnorm/mul_1:z:0be2/batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
z_mean/MatMul/ReadVariableOpReadVariableOp*z_mean_matmul_readvariableop_z_mean_kernel*
dtype0*
_output_shapes
:	��
z_mean/MatMulMatMulbe2/batchnorm/add_1:z:0$z_mean/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
z_mean/BiasAdd/ReadVariableOpReadVariableOp)z_mean_biasadd_readvariableop_z_mean_bias*
dtype0*
_output_shapes
:�
z_mean/BiasAddBiasAddz_mean/MatMul:product:0%z_mean/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������{
z_mean/leaky_re_lu/LeakyRelu	LeakyReluz_mean/BiasAdd:output:0*
alpha%
�#<*'
_output_shapes
:���������Q
bz/LogicalAnd/xConst*
dtype0
*
_output_shapes
: *
value	B
 ZQ
bz/LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: g
bz/LogicalAnd
LogicalAndbz/LogicalAnd/x:output:0bz/LogicalAnd/y:output:0*
_output_shapes
: k
!bz/moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:�
bz/moments/meanMean*z_mean/leaky_re_lu/LeakyRelu:activations:0*bz/moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(j
bz/moments/StopGradientStopGradientbz/moments/mean:output:0*
_output_shapes

:*
T0�
bz/moments/SquaredDifferenceSquaredDifference*z_mean/leaky_re_lu/LeakyRelu:activations:0 bz/moments/StopGradient:output:0*'
_output_shapes
:���������*
T0o
%bz/moments/variance/reduction_indicesConst*
dtype0*
_output_shapes
:*
valueB: �
bz/moments/varianceMean bz/moments/SquaredDifference:z:0.bz/moments/variance/reduction_indices:output:0*
_output_shapes

:*
	keep_dims(*
T0s
bz/moments/SqueezeSqueezebz/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 y
bz/moments/Squeeze_1Squeezebz/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 �
bz/AssignMovingAvg/decayConst*
valueB
 *
�#<*4
_class*
(&loc:@bz/AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
: 
!bz/AssignMovingAvg/ReadVariableOpReadVariableOp!bz_assignmovingavg_bz_moving_mean*
dtype0*
_output_shapes
:�
bz/AssignMovingAvg/subSub)bz/AssignMovingAvg/ReadVariableOp:value:0bz/moments/Squeeze:output:0*
_output_shapes
:*
T0*4
_class*
(&loc:@bz/AssignMovingAvg/bz/moving_mean�
bz/AssignMovingAvg/mulMulbz/AssignMovingAvg/sub:z:0!bz/AssignMovingAvg/decay:output:0*
T0*4
_class*
(&loc:@bz/AssignMovingAvg/bz/moving_mean*
_output_shapes
:�
&bz/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp!bz_assignmovingavg_bz_moving_meanbz/AssignMovingAvg/mul:z:0"^bz/AssignMovingAvg/ReadVariableOp*4
_class*
(&loc:@bz/AssignMovingAvg/bz/moving_mean*
dtype0*
_output_shapes
 �
bz/AssignMovingAvg_1/decayConst*
valueB
 *
�#<*:
_class0
.,loc:@bz/AssignMovingAvg_1/bz/moving_variance*
dtype0*
_output_shapes
: �
#bz/AssignMovingAvg_1/ReadVariableOpReadVariableOp'bz_assignmovingavg_1_bz_moving_variance*
dtype0*
_output_shapes
:�
bz/AssignMovingAvg_1/subSub+bz/AssignMovingAvg_1/ReadVariableOp:value:0bz/moments/Squeeze_1:output:0*
_output_shapes
:*
T0*:
_class0
.,loc:@bz/AssignMovingAvg_1/bz/moving_variance�
bz/AssignMovingAvg_1/mulMulbz/AssignMovingAvg_1/sub:z:0#bz/AssignMovingAvg_1/decay:output:0*
T0*:
_class0
.,loc:@bz/AssignMovingAvg_1/bz/moving_variance*
_output_shapes
:�
(bz/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp'bz_assignmovingavg_1_bz_moving_variancebz/AssignMovingAvg_1/mul:z:0$^bz/AssignMovingAvg_1/ReadVariableOp*:
_class0
.,loc:@bz/AssignMovingAvg_1/bz/moving_variance*
dtype0*
_output_shapes
 W
bz/batchnorm/add/yConst*
valueB
 *o�:*
dtype0*
_output_shapes
: z
bz/batchnorm/addAddV2bz/moments/Squeeze_1:output:0bz/batchnorm/add/y:output:0*
T0*
_output_shapes
:V
bz/batchnorm/RsqrtRsqrtbz/batchnorm/add:z:0*
T0*
_output_shapes
:�
bz/batchnorm/mul/ReadVariableOpReadVariableOp(bz_batchnorm_mul_readvariableop_bz_gamma*
dtype0*
_output_shapes
:}
bz/batchnorm/mulMulbz/batchnorm/Rsqrt:y:0'bz/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:�
bz/batchnorm/mul_1Mul*z_mean/leaky_re_lu/LeakyRelu:activations:0bz/batchnorm/mul:z:0*
T0*'
_output_shapes
:���������Y
bz/batchnorm/NegNegbz/moments/Squeeze:output:0*
T0*
_output_shapes
:j
bz/batchnorm/mul_2Mulbz/batchnorm/Neg:y:0bz/batchnorm/mul:z:0*
_output_shapes
:*
T0}
bz/batchnorm/add_1AddV2bz/batchnorm/mul_1:z:0bz/batchnorm/mul_2:z:0*
T0*'
_output_shapes
:����������
IdentityIdentitybz/batchnorm/add_1:z:0(^be1/AssignMovingAvg/AssignSubVariableOp#^be1/AssignMovingAvg/ReadVariableOp*^be1/AssignMovingAvg_1/AssignSubVariableOp%^be1/AssignMovingAvg_1/ReadVariableOp!^be1/batchnorm/mul/ReadVariableOp(^be2/AssignMovingAvg/AssignSubVariableOp#^be2/AssignMovingAvg/ReadVariableOp*^be2/AssignMovingAvg_1/AssignSubVariableOp%^be2/AssignMovingAvg_1/ReadVariableOp!^be2/batchnorm/mul/ReadVariableOp'^bz/AssignMovingAvg/AssignSubVariableOp"^bz/AssignMovingAvg/ReadVariableOp)^bz/AssignMovingAvg_1/AssignSubVariableOp$^bz/AssignMovingAvg_1/ReadVariableOp ^bz/batchnorm/mul/ReadVariableOp^e1/BiasAdd/ReadVariableOp^e1/MatMul/ReadVariableOp^e2/BiasAdd/ReadVariableOp^e2/MatMul/ReadVariableOp^z_mean/BiasAdd/ReadVariableOp^z_mean/MatMul/ReadVariableOp*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2H
"be1/AssignMovingAvg/ReadVariableOp"be1/AssignMovingAvg/ReadVariableOp2R
'be2/AssignMovingAvg/AssignSubVariableOp'be2/AssignMovingAvg/AssignSubVariableOp2D
 be1/batchnorm/mul/ReadVariableOp be1/batchnorm/mul/ReadVariableOp24
e2/MatMul/ReadVariableOpe2/MatMul/ReadVariableOp2F
!bz/AssignMovingAvg/ReadVariableOp!bz/AssignMovingAvg/ReadVariableOp26
e1/BiasAdd/ReadVariableOpe1/BiasAdd/ReadVariableOp2V
)be1/AssignMovingAvg_1/AssignSubVariableOp)be1/AssignMovingAvg_1/AssignSubVariableOp2B
bz/batchnorm/mul/ReadVariableOpbz/batchnorm/mul/ReadVariableOp2L
$be1/AssignMovingAvg_1/ReadVariableOp$be1/AssignMovingAvg_1/ReadVariableOp2H
"be2/AssignMovingAvg/ReadVariableOp"be2/AssignMovingAvg/ReadVariableOp2R
'be1/AssignMovingAvg/AssignSubVariableOp'be1/AssignMovingAvg/AssignSubVariableOp2<
z_mean/MatMul/ReadVariableOpz_mean/MatMul/ReadVariableOp2L
$be2/AssignMovingAvg_1/ReadVariableOp$be2/AssignMovingAvg_1/ReadVariableOp2T
(bz/AssignMovingAvg_1/AssignSubVariableOp(bz/AssignMovingAvg_1/AssignSubVariableOp2J
#bz/AssignMovingAvg_1/ReadVariableOp#bz/AssignMovingAvg_1/ReadVariableOp2V
)be2/AssignMovingAvg_1/AssignSubVariableOp)be2/AssignMovingAvg_1/AssignSubVariableOp2P
&bz/AssignMovingAvg/AssignSubVariableOp&bz/AssignMovingAvg/AssignSubVariableOp26
e2/BiasAdd/ReadVariableOpe2/BiasAdd/ReadVariableOp2D
 be2/batchnorm/mul/ReadVariableOp be2/batchnorm/mul/ReadVariableOp2>
z_mean/BiasAdd/ReadVariableOpz_mean/BiasAdd/ReadVariableOp24
e1/MatMul/ReadVariableOpe1/MatMul/ReadVariableOp: : : : : :	 :
 : : : : : : :( $
"
_user_specified_name
inputs/0:($
"
_user_specified_name
inputs/1: : 
�	
�
<__inference_e1_layer_call_and_return_conditional_losses_4396

inputs#
matmul_readvariableop_e1_kernel"
biasadd_readvariableop_e1_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpw
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_e1_kernel*
dtype0* 
_output_shapes
:
��j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_e1_bias*
dtype0*
_output_shapes	
:�w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*(
_output_shapes
:����������*
T0n
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:�����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�	
�
<__inference_e1_layer_call_and_return_conditional_losses_3773

inputs#
matmul_readvariableop_e1_kernel"
biasadd_readvariableop_e1_bias
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpw
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_e1_kernel*
dtype0* 
_output_shapes
:
��j
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:����������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_e1_bias*
dtype0*
_output_shapes	
:�w
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*(
_output_shapes
:����������*
T0n
leaky_re_lu/LeakyRelu	LeakyReluBiasAdd:output:0*
alpha%
�#<*(
_output_shapes
:�����������
IdentityIdentity#leaky_re_lu/LeakyRelu:activations:0^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*
T0*(
_output_shapes
:����������"
identityIdentity:output:0*/
_input_shapes
:����������::20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:& "
 
_user_specified_nameinputs: : 
�(
�
?__inference_model_layer_call_and_return_conditional_losses_4068

inputs
inputs_1(
$e1_statefulpartitionedcall_e1_kernel&
"e1_statefulpartitionedcall_e1_bias3
/be1_statefulpartitionedcall_be1_moving_variance)
%be1_statefulpartitionedcall_be1_gamma/
+be1_statefulpartitionedcall_be1_moving_mean(
$e2_statefulpartitionedcall_e2_kernel&
"e2_statefulpartitionedcall_e2_bias3
/be2_statefulpartitionedcall_be2_moving_variance)
%be2_statefulpartitionedcall_be2_gamma/
+be2_statefulpartitionedcall_be2_moving_mean0
,z_mean_statefulpartitionedcall_z_mean_kernel.
*z_mean_statefulpartitionedcall_z_mean_bias1
-bz_statefulpartitionedcall_bz_moving_variance'
#bz_statefulpartitionedcall_bz_gamma-
)bz_statefulpartitionedcall_bz_moving_mean
identity��be1/StatefulPartitionedCall�be2/StatefulPartitionedCall�bz/StatefulPartitionedCall�e1/StatefulPartitionedCall�e2/StatefulPartitionedCall�z_mean/StatefulPartitionedCall�
concatenate/PartitionedCallPartitionedCallinputsinputs_1**
config_proto

CPU

GPU 2J 8*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3755*N
fIRG
E__inference_concatenate_layer_call_and_return_conditional_losses_3748*
Tout
2�
e1/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0$e1_statefulpartitionedcall_e1_kernel"e1_statefulpartitionedcall_e1_bias*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3779*E
f@R>
<__inference_e1_layer_call_and_return_conditional_losses_3773*
Tout
2**
config_proto

CPU

GPU 2J 8�
be1/StatefulPartitionedCallStatefulPartitionedCall#e1/StatefulPartitionedCall:output:0/be1_statefulpartitionedcall_be1_moving_variance%be1_statefulpartitionedcall_be1_gamma+be1_statefulpartitionedcall_be1_moving_mean*
Tout
2**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3454*F
fAR?
=__inference_be1_layer_call_and_return_conditional_losses_3453�
dropout/PartitionedCallPartitionedCall$be1/StatefulPartitionedCall:output:0*
Tin
2*(
_output_shapes
:����������*+
_gradient_op_typePartitionedCall-3850*J
fERC
A__inference_dropout_layer_call_and_return_conditional_losses_3838*
Tout
2**
config_proto

CPU

GPU 2J 8�
e2/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0$e2_statefulpartitionedcall_e2_kernel"e2_statefulpartitionedcall_e2_bias**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3872*E
f@R>
<__inference_e2_layer_call_and_return_conditional_losses_3866*
Tout
2�
be2/StatefulPartitionedCallStatefulPartitionedCall#e2/StatefulPartitionedCall:output:0/be2_statefulpartitionedcall_be2_moving_variance%be2_statefulpartitionedcall_be2_gamma+be2_statefulpartitionedcall_be2_moving_mean**
config_proto

CPU

GPU 2J 8*(
_output_shapes
:����������*
Tin
2*+
_gradient_op_typePartitionedCall-3592*F
fAR?
=__inference_be2_layer_call_and_return_conditional_losses_3591*
Tout
2�
z_mean/StatefulPartitionedCallStatefulPartitionedCall$be2/StatefulPartitionedCall:output:0,z_mean_statefulpartitionedcall_z_mean_kernel*z_mean_statefulpartitionedcall_z_mean_bias*
Tout
2**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3920*I
fDRB
@__inference_z_mean_layer_call_and_return_conditional_losses_3914�
bz/StatefulPartitionedCallStatefulPartitionedCall'z_mean/StatefulPartitionedCall:output:0-bz_statefulpartitionedcall_bz_moving_variance#bz_statefulpartitionedcall_bz_gamma)bz_statefulpartitionedcall_bz_moving_mean**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-3730*E
f@R>
<__inference_bz_layer_call_and_return_conditional_losses_3729*
Tout
2�
IdentityIdentity#bz/StatefulPartitionedCall:output:0^be1/StatefulPartitionedCall^be2/StatefulPartitionedCall^bz/StatefulPartitionedCall^e1/StatefulPartitionedCall^e2/StatefulPartitionedCall^z_mean/StatefulPartitionedCall*'
_output_shapes
:���������*
T0"
identityIdentity:output:0*v
_input_shapese
c:����������:���������:::::::::::::::2:
be1/StatefulPartitionedCallbe1/StatefulPartitionedCall2:
be2/StatefulPartitionedCallbe2/StatefulPartitionedCall28
bz/StatefulPartitionedCallbz/StatefulPartitionedCall28
e1/StatefulPartitionedCalle1/StatefulPartitionedCall2@
z_mean/StatefulPartitionedCallz_mean/StatefulPartitionedCall28
e2/StatefulPartitionedCalle2/StatefulPartitionedCall:& "
 
_user_specified_nameinputs:&"
 
_user_specified_nameinputs: : : : : : : :	 :
 : : : : : : 
�
�
=__inference_be1_layer_call_and_return_conditional_losses_4459

inputs0
,batchnorm_readvariableop_be1_moving_variance*
&batchnorm_mul_readvariableop_be1_gamma.
*batchnorm_readvariableop_1_be1_moving_mean
identity��batchnorm/ReadVariableOp�batchnorm/ReadVariableOp_1�batchnorm/mul/ReadVariableOpN
LogicalAnd/xConst*
value	B
 Z *
dtype0
*
_output_shapes
: N
LogicalAnd/yConst*
value	B
 Z*
dtype0
*
_output_shapes
: ^

LogicalAnd
LogicalAndLogicalAnd/x:output:0LogicalAnd/y:output:0*
_output_shapes
: �
batchnorm/ReadVariableOpReadVariableOp,batchnorm_readvariableop_be1_moving_variance*
dtype0*
_output_shapes	
:�T
batchnorm/add/yConst*
dtype0*
_output_shapes
: *
valueB
 *o�:x
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
_output_shapes	
:�*
T0Q
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:��
batchnorm/mul/ReadVariableOpReadVariableOp&batchnorm_mul_readvariableop_be1_gamma*
dtype0*
_output_shapes	
:�u
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:�d
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:�����������
batchnorm/ReadVariableOp_1ReadVariableOp*batchnorm_readvariableop_1_be1_moving_mean*
dtype0*
_output_shapes	
:�^
batchnorm/NegNeg"batchnorm/ReadVariableOp_1:value:0*
T0*
_output_shapes	
:�b
batchnorm/mul_2Mulbatchnorm/Neg:y:0batchnorm/mul:z:0*
T0*
_output_shapes	
:�u
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/mul_2:z:0*
T0*(
_output_shapes
:�����������
IdentityIdentitybatchnorm/add_1:z:0^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/mul/ReadVariableOp*(
_output_shapes
:����������*
T0"
identityIdentity:output:0*3
_input_shapes"
 :����������:::28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_124
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:& "
 
_user_specified_nameinputs: : : 
�
�
"__inference_signature_wrapper_4110
b
x_input%
!statefulpartitionedcall_e1_kernel#
statefulpartitionedcall_e1_bias/
+statefulpartitionedcall_be1_moving_variance%
!statefulpartitionedcall_be1_gamma+
'statefulpartitionedcall_be1_moving_mean%
!statefulpartitionedcall_e2_kernel#
statefulpartitionedcall_e2_bias/
+statefulpartitionedcall_be2_moving_variance%
!statefulpartitionedcall_be2_gamma+
'statefulpartitionedcall_be2_moving_mean)
%statefulpartitionedcall_z_mean_kernel'
#statefulpartitionedcall_z_mean_bias.
*statefulpartitionedcall_bz_moving_variance$
 statefulpartitionedcall_bz_gamma*
&statefulpartitionedcall_bz_moving_mean
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallx_inputb!statefulpartitionedcall_e1_kernelstatefulpartitionedcall_e1_bias+statefulpartitionedcall_be1_moving_variance!statefulpartitionedcall_be1_gamma'statefulpartitionedcall_be1_moving_mean!statefulpartitionedcall_e2_kernelstatefulpartitionedcall_e2_bias+statefulpartitionedcall_be2_moving_variance!statefulpartitionedcall_be2_gamma'statefulpartitionedcall_be2_moving_mean%statefulpartitionedcall_z_mean_kernel#statefulpartitionedcall_z_mean_bias*statefulpartitionedcall_bz_moving_variance statefulpartitionedcall_bz_gamma&statefulpartitionedcall_bz_moving_mean**
config_proto

CPU

GPU 2J 8*
Tin
2*'
_output_shapes
:���������*+
_gradient_op_typePartitionedCall-4092*(
f#R!
__inference__wrapped_model_3319*
Tout
2�
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:���������"
identityIdentity:output:0*v
_input_shapese
c:���������:����������:::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall: : : : : : : :	 :
 : : : : : : :! 

_user_specified_nameB:'#
!
_user_specified_name	x_input"wL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*�
serving_default�
<
x_input1
serving_default_x_input:0����������
/
B*
serving_default_B:0���������6
bz0
StatefulPartitionedCall:0���������tensorflow/serving/predict*>
__saved_model_init_op%#
__saved_model_init_op

NoOp:��
�L
layer-0
layer-1
layer-2
layer_with_weights-0
layer-3
layer_with_weights-1
layer-4
layer-5
layer_with_weights-2
layer-6
layer_with_weights-3
layer-7
	layer_with_weights-4
	layer-8

layer_with_weights-5

layer-9
trainable_variables
	variables
regularization_losses
	keras_api

signatures
_default_save_signature
+�&call_and_return_all_conditional_losses
�__call__"�H
_tf_keras_model�H{"class_name": "Model", "name": "model", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3000], "dtype": "float32", "sparse": false, "name": "x_input"}, "name": "x_input", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": [null, 4], "dtype": "float32", "sparse": false, "name": "B"}, "name": "B", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": null, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["x_input", 0, 0, {}], ["B", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "e1", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "e1", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "be1", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "be1", "inbound_nodes": [[["e1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["be1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "e2", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "e2", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "be2", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "be2", "inbound_nodes": [[["e2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 30, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_mean", "inbound_nodes": [[["be2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "bz", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "bz", "inbound_nodes": [[["z_mean", 0, 0, {}]]]}], "input_layers": [["x_input", 0, 0], ["B", 0, 0]], "output_layers": [["bz", 0, 0]]}, "input_spec": [null, null], "keras_version": "2.2.4-tf", "backend": "tensorflow", "model_config": {"class_name": "Model", "config": {"name": "model", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 3000], "dtype": "float32", "sparse": false, "name": "x_input"}, "name": "x_input", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": [null, 4], "dtype": "float32", "sparse": false, "name": "B"}, "name": "B", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": null, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["x_input", 0, 0, {}], ["B", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "e1", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "e1", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "be1", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "be1", "inbound_nodes": [[["e1", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["be1", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "e2", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "e2", "inbound_nodes": [[["dropout", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "be2", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "be2", "inbound_nodes": [[["e2", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 30, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "z_mean", "inbound_nodes": [[["be2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "bz", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "bz", "inbound_nodes": [[["z_mean", 0, 0, {}]]]}], "input_layers": [["x_input", 0, 0], ["B", 0, 0]], "output_layers": [["bz", 0, 0]]}}}
�
trainable_variables
	variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "InputLayer", "name": "x_input", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": [null, 3000], "config": {"batch_input_shape": [null, 3000], "dtype": "float32", "sparse": false, "name": "x_input"}}
�
trainable_variables
	variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "InputLayer", "name": "B", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": [null, 4], "config": {"batch_input_shape": [null, 4], "dtype": "float32", "sparse": false, "name": "B"}}
�
trainable_variables
	variables
regularization_losses
	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Concatenate", "name": "concatenate", "trainable": null, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "concatenate", "trainable": null, "dtype": "float32", "axis": -1}}
�

activation

kernel
bias
trainable_variables
 	variables
!regularization_losses
"	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "e1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "e1", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 3004}}}}
�
#axis
	$gamma
%moving_mean
&moving_variance
'trainable_variables
(	variables
)regularization_losses
*	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "BatchNormalization", "name": "be1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "be1", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {"1": 128}}}}
�
+trainable_variables
,	variables
-regularization_losses
.	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dropout", "name": "dropout", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.1, "noise_shape": null, "seed": null}}
�

activation

/kernel
0bias
1trainable_variables
2	variables
3regularization_losses
4	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "e2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "e2", "trainable": true, "dtype": "float32", "units": 128, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 128}}}}
�
5axis
	6gamma
7moving_mean
8moving_variance
9trainable_variables
:	variables
;regularization_losses
<	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "BatchNormalization", "name": "be2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "be2", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {"1": 128}}}}
�

activation

=kernel
>bias
?trainable_variables
@	variables
Aregularization_losses
B	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "Dense", "name": "z_mean", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "z_mean", "trainable": true, "dtype": "float32", "units": 30, "activation": {"class_name": "LeakyReLU", "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}, "use_bias": true, "kernel_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 128}}}}
�
Caxis
	Dgamma
Emoving_mean
Fmoving_variance
Gtrainable_variables
H	variables
Iregularization_losses
J	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "BatchNormalization", "name": "bz", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "config": {"name": "bz", "trainable": true, "dtype": "float32", "axis": [1], "momentum": 0.99, "epsilon": 0.001, "center": false, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 2, "max_ndim": null, "min_ndim": null, "axes": {"1": 30}}}}
_
0
1
$2
/3
04
65
=6
>7
D8"
trackable_list_wrapper
�
0
1
$2
%3
&4
/5
06
67
78
89
=10
>11
D12
E13
F14"
trackable_list_wrapper
 "
trackable_list_wrapper
�
Klayer_regularization_losses
trainable_variables

Llayers
	variables
Mnon_trainable_variables
Nmetrics
regularization_losses
�__call__
_default_save_signature
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
-
�serving_default"
signature_map
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Olayer_regularization_losses
trainable_variables

Players
	variables
Qnon_trainable_variables
Rmetrics
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Slayer_regularization_losses
trainable_variables

Tlayers
	variables
Unon_trainable_variables
Vmetrics
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Wlayer_regularization_losses
trainable_variables

Xlayers
	variables
Ynon_trainable_variables
Zmetrics
regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
[trainable_variables
\	variables
]regularization_losses
^	keras_api
�__call__
+�&call_and_return_all_conditional_losses"�
_tf_keras_layer�{"class_name": "LeakyReLU", "name": "leaky_re_lu", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "config": {"name": "leaky_re_lu", "trainable": true, "dtype": "float32", "alpha": 0.009999999776482582}}
:
��2	e1/kernel
:�2e1/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
_layer_regularization_losses
trainable_variables

`layers
 	variables
anon_trainable_variables
bmetrics
!regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
:�2	be1/gamma
 :� (2be1/moving_mean
$:"� (2be1/moving_variance
'
$0"
trackable_list_wrapper
5
$0
%1
&2"
trackable_list_wrapper
 "
trackable_list_wrapper
�
clayer_regularization_losses
'trainable_variables

dlayers
(	variables
enon_trainable_variables
fmetrics
)regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
glayer_regularization_losses
+trainable_variables

hlayers
,	variables
inon_trainable_variables
jmetrics
-regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
:
��2	e2/kernel
:�2e2/bias
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
 "
trackable_list_wrapper
�
klayer_regularization_losses
1trainable_variables

llayers
2	variables
mnon_trainable_variables
nmetrics
3regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
:�2	be2/gamma
 :� (2be2/moving_mean
$:"� (2be2/moving_variance
'
60"
trackable_list_wrapper
5
60
71
82"
trackable_list_wrapper
 "
trackable_list_wrapper
�
olayer_regularization_losses
9trainable_variables

players
:	variables
qnon_trainable_variables
rmetrics
;regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 :	�2z_mean/kernel
:2z_mean/bias
.
=0
>1"
trackable_list_wrapper
.
=0
>1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
slayer_regularization_losses
?trainable_variables

tlayers
@	variables
unon_trainable_variables
vmetrics
Aregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
:2bz/gamma
: (2bz/moving_mean
":  (2bz/moving_variance
'
D0"
trackable_list_wrapper
5
D0
E1
F2"
trackable_list_wrapper
 "
trackable_list_wrapper
�
wlayer_regularization_losses
Gtrainable_variables

xlayers
H	variables
ynon_trainable_variables
zmetrics
Iregularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
f
0
1
2
3
4
5
6
7
	8

9"
trackable_list_wrapper
J
%0
&1
72
83
E4
F5"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
{layer_regularization_losses
[trainable_variables

|layers
\	variables
}non_trainable_variables
~metrics
]regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
70
81"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
E0
F1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�2�
__inference__wrapped_model_3319�
���
FullArgSpec
args� 
varargsjargs
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *I�F
D�A
"�
x_input����������
�
B���������
�2�
?__inference_model_layer_call_and_return_conditional_losses_4250
?__inference_model_layer_call_and_return_conditional_losses_4330
?__inference_model_layer_call_and_return_conditional_losses_3952
?__inference_model_layer_call_and_return_conditional_losses_3983�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
$__inference_model_layer_call_fn_4034
$__inference_model_layer_call_fn_4351
$__inference_model_layer_call_fn_4087
$__inference_model_layer_call_fn_4372�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkwjkwargs
defaults� 

kwonlyargs�

jtraining%
kwonlydefaults�

trainingp 
annotations� *
 
�2�
*__inference_concatenate_layer_call_fn_4385�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
E__inference_concatenate_layer_call_and_return_conditional_losses_4379�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
!__inference_e1_layer_call_fn_4403�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
<__inference_e1_layer_call_and_return_conditional_losses_4396�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
"__inference_be1_layer_call_fn_4475
"__inference_be1_layer_call_fn_4467�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
=__inference_be1_layer_call_and_return_conditional_losses_4438
=__inference_be1_layer_call_and_return_conditional_losses_4459�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
&__inference_dropout_layer_call_fn_4510
&__inference_dropout_layer_call_fn_4505�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
A__inference_dropout_layer_call_and_return_conditional_losses_4495
A__inference_dropout_layer_call_and_return_conditional_losses_4500�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
!__inference_e2_layer_call_fn_4528�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
<__inference_e2_layer_call_and_return_conditional_losses_4521�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
"__inference_be2_layer_call_fn_4592
"__inference_be2_layer_call_fn_4600�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
=__inference_be2_layer_call_and_return_conditional_losses_4563
=__inference_be2_layer_call_and_return_conditional_losses_4584�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
%__inference_z_mean_layer_call_fn_4618�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
@__inference_z_mean_layer_call_and_return_conditional_losses_4611�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2�
!__inference_bz_layer_call_fn_4690
!__inference_bz_layer_call_fn_4682�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
�2�
<__inference_bz_layer_call_and_return_conditional_losses_4674
<__inference_bz_layer_call_and_return_conditional_losses_4653�
���
FullArgSpec)
args!�
jself
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults� 
annotations� *
 
2B0
"__inference_signature_wrapper_4110Bx_input
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�2��
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 |
"__inference_be1_layer_call_fn_4475V&$%4�1
*�'
!�
inputs����������
p 
� "������������
A__inference_dropout_layer_call_and_return_conditional_losses_4495^4�1
*�'
!�
inputs����������
p
� "&�#
�
0����������
� �
?__inference_model_layer_call_and_return_conditional_losses_3952�%&$/0786=>EFD[�X
Q�N
D�A
"�
x_input����������
�
B���������
p

 
� "%�"
�
0���������
� �
<__inference_bz_layer_call_and_return_conditional_losses_4674aFDE3�0
)�&
 �
inputs���������
p 
� "%�"
�
0���������
� {
&__inference_dropout_layer_call_fn_4505Q4�1
*�'
!�
inputs����������
p
� "�����������{
&__inference_dropout_layer_call_fn_4510Q4�1
*�'
!�
inputs����������
p 
� "�����������|
"__inference_be2_layer_call_fn_4600V8674�1
*�'
!�
inputs����������
p 
� "������������
=__inference_be1_layer_call_and_return_conditional_losses_4459c&$%4�1
*�'
!�
inputs����������
p 
� "&�#
�
0����������
� �
=__inference_be2_layer_call_and_return_conditional_losses_4563c7864�1
*�'
!�
inputs����������
p
� "&�#
�
0����������
� �
@__inference_z_mean_layer_call_and_return_conditional_losses_4611]=>0�-
&�#
!�
inputs����������
� "%�"
�
0���������
� �
$__inference_model_layer_call_fn_4034�%&$/0786=>EFD[�X
Q�N
D�A
"�
x_input����������
�
B���������
p

 
� "�����������
?__inference_model_layer_call_and_return_conditional_losses_4250�%&$/0786=>EFDc�`
Y�V
L�I
#� 
inputs/0����������
"�
inputs/1���������
p

 
� "%�"
�
0���������
� �
*__inference_concatenate_layer_call_fn_4385x[�X
Q�N
L�I
#� 
inputs/0����������
"�
inputs/1���������
� "������������
=__inference_be2_layer_call_and_return_conditional_losses_4584c8674�1
*�'
!�
inputs����������
p 
� "&�#
�
0����������
� �
__inference__wrapped_model_3319�&$%/0867=>FDES�P
I�F
D�A
"�
x_input����������
�
B���������
� "'�$
"
bz�
bz����������
?__inference_model_layer_call_and_return_conditional_losses_3983�&$%/0867=>FDE[�X
Q�N
D�A
"�
x_input����������
�
B���������
p 

 
� "%�"
�
0���������
� v
!__inference_e2_layer_call_fn_4528Q/00�-
&�#
!�
inputs����������
� "�����������y
!__inference_bz_layer_call_fn_4682TEFD3�0
)�&
 �
inputs���������
p
� "����������y
%__inference_z_mean_layer_call_fn_4618P=>0�-
&�#
!�
inputs����������
� "����������|
"__inference_be2_layer_call_fn_4592V7864�1
*�'
!�
inputs����������
p
� "������������
"__inference_signature_wrapper_4110�&$%/0867=>FDE^�[
� 
T�Q
-
x_input"�
x_input����������
 
B�
B���������"'�$
"
bz�
bz���������y
!__inference_bz_layer_call_fn_4690TFDE3�0
)�&
 �
inputs���������
p 
� "�����������
?__inference_model_layer_call_and_return_conditional_losses_4330�&$%/0867=>FDEc�`
Y�V
L�I
#� 
inputs/0����������
"�
inputs/1���������
p 

 
� "%�"
�
0���������
� �
$__inference_model_layer_call_fn_4351�%&$/0786=>EFDc�`
Y�V
L�I
#� 
inputs/0����������
"�
inputs/1���������
p

 
� "�����������
E__inference_concatenate_layer_call_and_return_conditional_losses_4379�[�X
Q�N
L�I
#� 
inputs/0����������
"�
inputs/1���������
� "&�#
�
0����������
� �
A__inference_dropout_layer_call_and_return_conditional_losses_4500^4�1
*�'
!�
inputs����������
p 
� "&�#
�
0����������
� �
$__inference_model_layer_call_fn_4087�&$%/0867=>FDE[�X
Q�N
D�A
"�
x_input����������
�
B���������
p 

 
� "�����������
$__inference_model_layer_call_fn_4372�&$%/0867=>FDEc�`
Y�V
L�I
#� 
inputs/0����������
"�
inputs/1���������
p 

 
� "�����������
<__inference_e2_layer_call_and_return_conditional_losses_4521^/00�-
&�#
!�
inputs����������
� "&�#
�
0����������
� v
!__inference_e1_layer_call_fn_4403Q0�-
&�#
!�
inputs����������
� "������������
<__inference_bz_layer_call_and_return_conditional_losses_4653aEFD3�0
)�&
 �
inputs���������
p
� "%�"
�
0���������
� �
=__inference_be1_layer_call_and_return_conditional_losses_4438c%&$4�1
*�'
!�
inputs����������
p
� "&�#
�
0����������
� �
<__inference_e1_layer_call_and_return_conditional_losses_4396^0�-
&�#
!�
inputs����������
� "&�#
�
0����������
� |
"__inference_be1_layer_call_fn_4467V%&$4�1
*�'
!�
inputs����������
p
� "�����������