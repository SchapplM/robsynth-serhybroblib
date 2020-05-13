% Explicit kinematic constraints of
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% jv [16x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = palh1m2DE1_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_kinconstr_expl_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_kinconstr_expl_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:29
% EndTime: 2020-05-01 20:55:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (90->24), mult. (156->44), div. (0->0), fcn. (261->23), ass. (0->32)
t26 = sin(qJ(3));
t27 = sin(qJ(2));
t30 = cos(qJ(3));
t31 = cos(qJ(2));
t12 = -t27 * t26 + t30 * t31;
t14 = t31 * t26 + t30 * t27;
t28 = sin(pkin(18));
t32 = cos(pkin(18));
t4 = t12 * t32 + t28 * t14;
t33 = cos(pkin(17));
t29 = sin(pkin(17));
t25 = cos(pkin(19));
t24 = cos(pkin(20));
t23 = cos(pkin(22));
t22 = sin(pkin(19));
t21 = sin(pkin(20));
t20 = sin(pkin(22));
t19 = pkin(22) + pkin(21);
t18 = cos(t19);
t17 = sin(t19);
t15 = t29 * t28 + t33 * t32;
t13 = t28 * t33 - t32 * t29;
t11 = t32 * t20 - t28 * t23;
t10 = t28 * t20 + t23 * t32;
t9 = t26 * t22 - t30 * t25;
t8 = t30 * t22 + t26 * t25;
t7 = atan2(t9, -t8);
t5 = -t28 * t12 + t14 * t32;
t3 = t21 * t5 - t24 * t4;
t2 = t21 * t4 + t5 * t24;
t1 = atan2(t17 * t4 + t5 * t18, -t5 * t17 + t4 * t18);
t6 = [qJ(1); qJ(2); qJ(3); atan2(-t3 * t17 + t2 * t18, t2 * t17 + t3 * t18); qJ(4); atan2(t27 * t13 + t15 * t31, t13 * t31 - t15 * t27); atan2(t10 * t31 - t11 * t27, t27 * t10 + t11 * t31); atan2(t9, t8); t7; t1; atan2(t27 * t28 + t31 * t32, -t27 * t32 + t31 * t28); t7; t1; 0; 0; 0;];
jv = t6(:);
