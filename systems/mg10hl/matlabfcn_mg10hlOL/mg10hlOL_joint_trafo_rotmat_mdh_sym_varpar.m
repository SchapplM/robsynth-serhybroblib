% Calculate homogenous joint transformation matrices for
% mg10hlOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% T_mdh [4x4x15]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:06
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = mg10hlOL_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlOL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlOL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:05:44
% EndTime: 2020-04-11 13:05:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (37->28), mult. (20->12), div. (0->0), fcn. (76->28), ass. (0->273)
unknown=NaN(60,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t5 = cos(pkin(14));
t6 = cos(qJ(3));
t8 = sin(pkin(14));
t9 = sin(qJ(3));
t11 = -t5 * t6 + t8 * t9;
t14 = t5 * t9 + t8 * t6;
t17 = cos(pkin(16));
t18 = -qJ(4) + pkin(15);
t19 = cos(t18);
t21 = sin(pkin(16));
t22 = sin(t18);
t24 = -t17 * t19 - t21 * t22;
t27 = -t17 * t22 + t21 * t19;
t30 = cos(qJ(5));
t31 = sin(qJ(5));
t32 = cos(qJ(6));
t33 = sin(qJ(6));
t34 = sin(qJ(7));
t35 = cos(qJ(7));
t36 = cos(qJ(8));
t37 = sin(qJ(8));
t38 = cos(qJ(9));
t39 = sin(qJ(9));
t40 = cos(qJ(10));
t41 = sin(qJ(10));
t42 = cos(qJ(11));
t43 = sin(qJ(11));
t44 = cos(qJ(12));
t45 = sin(qJ(12));
unknown(1,1) = t1;
unknown(1,2) = -t2;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(2,1) = t2;
unknown(2,2) = t1;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.1e1;
unknown(3,4) = pkin(10);
unknown(4,1) = 0.0e0;
unknown(4,2) = 0.0e0;
unknown(4,3) = 0.0e0;
unknown(4,4) = 0.1e1;
unknown(5,1) = -t3;
unknown(5,2) = -t4;
unknown(5,3) = 0.0e0;
unknown(5,4) = pkin(12);
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.0e0;
unknown(7,1) = t4;
unknown(7,2) = -t3;
unknown(7,3) = 0.0e0;
unknown(7,4) = 0.0e0;
unknown(8,1) = 0.0e0;
unknown(8,2) = 0.0e0;
unknown(8,3) = 0.0e0;
unknown(8,4) = 0.1e1;
unknown(9,1) = t11;
unknown(9,2) = t14;
unknown(9,3) = 0.0e0;
unknown(9,4) = t5 * pkin(1);
unknown(10,1) = -t14;
unknown(10,2) = t11;
unknown(10,3) = 0.0e0;
unknown(10,4) = t8 * pkin(1);
unknown(11,1) = 0.0e0;
unknown(11,2) = 0.0e0;
unknown(11,3) = 0.1e1;
unknown(11,4) = 0.0e0;
unknown(12,1) = 0.0e0;
unknown(12,2) = 0.0e0;
unknown(12,3) = 0.0e0;
unknown(12,4) = 0.1e1;
unknown(13,1) = t24;
unknown(13,2) = t27;
unknown(13,3) = 0.0e0;
unknown(13,4) = t17 * pkin(2);
unknown(14,1) = -t27;
unknown(14,2) = t24;
unknown(14,3) = 0.0e0;
unknown(14,4) = t21 * pkin(2);
unknown(15,1) = 0.0e0;
unknown(15,2) = 0.0e0;
unknown(15,3) = 0.1e1;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = -t30;
unknown(17,2) = t31;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(6);
unknown(18,1) = -t31;
unknown(18,2) = -t30;
unknown(18,3) = 0.0e0;
unknown(18,4) = 0.0e0;
unknown(19,1) = 0.0e0;
unknown(19,2) = 0.0e0;
unknown(19,3) = 0.1e1;
unknown(19,4) = 0.0e0;
unknown(20,1) = 0.0e0;
unknown(20,2) = 0.0e0;
unknown(20,3) = 0.0e0;
unknown(20,4) = 0.1e1;
unknown(21,1) = 0.0e0;
unknown(21,2) = 0.0e0;
unknown(21,3) = 0.1e1;
unknown(21,4) = pkin(11);
unknown(22,1) = t32;
unknown(22,2) = -t33;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = t33;
unknown(23,2) = t32;
unknown(23,3) = 0.0e0;
unknown(23,4) = 0.0e0;
unknown(24,1) = 0.0e0;
unknown(24,2) = 0.0e0;
unknown(24,3) = 0.0e0;
unknown(24,4) = 0.1e1;
unknown(25,1) = t34;
unknown(25,2) = t35;
unknown(25,3) = 0.0e0;
unknown(25,4) = 0.0e0;
unknown(26,1) = 0.0e0;
unknown(26,2) = 0.0e0;
unknown(26,3) = 0.1e1;
unknown(26,4) = 0.0e0;
unknown(27,1) = t35;
unknown(27,2) = -t34;
unknown(27,3) = 0.0e0;
unknown(27,4) = pkin(9);
unknown(28,1) = 0.0e0;
unknown(28,2) = 0.0e0;
unknown(28,3) = 0.0e0;
unknown(28,4) = 0.1e1;
unknown(29,1) = 0.0e0;
unknown(29,2) = 0.0e0;
unknown(29,3) = 0.1e1;
unknown(29,4) = pkin(8);
unknown(30,1) = t36;
unknown(30,2) = -t37;
unknown(30,3) = 0.0e0;
unknown(30,4) = 0.0e0;
unknown(31,1) = t37;
unknown(31,2) = t36;
unknown(31,3) = 0.0e0;
unknown(31,4) = 0.0e0;
unknown(32,1) = 0.0e0;
unknown(32,2) = 0.0e0;
unknown(32,3) = 0.0e0;
unknown(32,4) = 0.1e1;
unknown(33,1) = t38;
unknown(33,2) = -t39;
unknown(33,3) = 0.0e0;
unknown(33,4) = 0.0e0;
unknown(34,1) = t39;
unknown(34,2) = t38;
unknown(34,3) = 0.0e0;
unknown(34,4) = 0.0e0;
unknown(35,1) = 0.0e0;
unknown(35,2) = 0.0e0;
unknown(35,3) = 0.1e1;
unknown(35,4) = 0.0e0;
unknown(36,1) = 0.0e0;
unknown(36,2) = 0.0e0;
unknown(36,3) = 0.0e0;
unknown(36,4) = 0.1e1;
unknown(37,1) = t40;
unknown(37,2) = -t41;
unknown(37,3) = 0.0e0;
unknown(37,4) = -pkin(13);
unknown(38,1) = 0.0e0;
unknown(38,2) = 0.0e0;
unknown(38,3) = -0.1e1;
unknown(38,4) = 0.0e0;
unknown(39,1) = t41;
unknown(39,2) = t40;
unknown(39,3) = 0.0e0;
unknown(39,4) = 0.0e0;
unknown(40,1) = 0.0e0;
unknown(40,2) = 0.0e0;
unknown(40,3) = 0.0e0;
unknown(40,4) = 0.1e1;
unknown(41,1) = -t42;
unknown(41,2) = t43;
unknown(41,3) = 0.0e0;
unknown(41,4) = pkin(5);
unknown(42,1) = -t43;
unknown(42,2) = -t42;
unknown(42,3) = 0.0e0;
unknown(42,4) = 0.0e0;
unknown(43,1) = 0.0e0;
unknown(43,2) = 0.0e0;
unknown(43,3) = 0.1e1;
unknown(43,4) = 0.0e0;
unknown(44,1) = 0.0e0;
unknown(44,2) = 0.0e0;
unknown(44,3) = 0.0e0;
unknown(44,4) = 0.1e1;
unknown(45,1) = -t44;
unknown(45,2) = t45;
unknown(45,3) = 0.0e0;
unknown(45,4) = pkin(4);
unknown(46,1) = -t45;
unknown(46,2) = -t44;
unknown(46,3) = 0.0e0;
unknown(46,4) = 0.0e0;
unknown(47,1) = 0.0e0;
unknown(47,2) = 0.0e0;
unknown(47,3) = 0.1e1;
unknown(47,4) = 0.0e0;
unknown(48,1) = 0.0e0;
unknown(48,2) = 0.0e0;
unknown(48,3) = 0.0e0;
unknown(48,4) = 0.1e1;
unknown(49,1) = 0.0e0;
unknown(49,2) = 0.0e0;
unknown(49,3) = 0.1e1;
unknown(49,4) = qJ(13);
unknown(50,1) = -0.1e1;
unknown(50,2) = 0.0e0;
unknown(50,3) = 0.0e0;
unknown(50,4) = 0.0e0;
unknown(51,1) = 0.0e0;
unknown(51,2) = -0.1e1;
unknown(51,3) = 0.0e0;
unknown(51,4) = 0.0e0;
unknown(52,1) = 0.0e0;
unknown(52,2) = 0.0e0;
unknown(52,3) = 0.0e0;
unknown(52,4) = 0.1e1;
unknown(53,1) = 0.1e1;
unknown(53,2) = 0.0e0;
unknown(53,3) = 0.0e0;
unknown(53,4) = -pkin(3);
unknown(54,1) = 0.0e0;
unknown(54,2) = 0.1e1;
unknown(54,3) = 0.0e0;
unknown(54,4) = 0.0e0;
unknown(55,1) = 0.0e0;
unknown(55,2) = 0.0e0;
unknown(55,3) = 0.1e1;
unknown(55,4) = 0.0e0;
unknown(56,1) = 0.0e0;
unknown(56,2) = 0.0e0;
unknown(56,3) = 0.0e0;
unknown(56,4) = 0.1e1;
unknown(57,1) = 0.0e0;
unknown(57,2) = 0.0e0;
unknown(57,3) = 0.1e1;
unknown(57,4) = -pkin(7);
unknown(58,1) = -0.1e1;
unknown(58,2) = 0.0e0;
unknown(58,3) = 0.0e0;
unknown(58,4) = 0.0e0;
unknown(59,1) = 0.0e0;
unknown(59,2) = -0.1e1;
unknown(59,3) = 0.0e0;
unknown(59,4) = 0.0e0;
unknown(60,1) = 0.0e0;
unknown(60,2) = 0.0e0;
unknown(60,3) = 0.0e0;
unknown(60,4) = 0.1e1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,15);             % numerisch
else,                         T_mdh = sym('xx', [4,4,15]); end % symbolisch

for i = 1:15
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
